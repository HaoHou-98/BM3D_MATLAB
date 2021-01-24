function [imout] = BM3D_matlab(im,imout,sigma)

[M1,M2]=size(im); 
if sigma < 40
N1    = 8;     %% block size.
N2    = 32;    %% the number of block in a group.
Ns    = 39;    %% the size of block-matching neighbor. 
Nstep = 3;
else
N1    = 11;     %% block size.
N2    = 32;    %% the number of block in a group.
Ns    = 39;    %% the size of block-matching neighbor. 
Nstep = 6;
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Create 2D Block transform matrices.
%%%%
[Tfor, Tinv]   = getTransfMatrix(N1, 'dct', 0);     %% get (normalized) forward and inverse transform matrices
Tfor=single(Tfor);
Tinv=single(Tinv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Create Third Dimension transform matrices.
%%%%
hadper_trans_single_den = cell(1,1000);
inverse_hadper_trans_single_den = cell(1,1000);
for hpow = 0:ceil(log2(N2))
    h = 2^hpow;
    [Tfor3rd, Tinv3rd]   = getTransfMatrix(h, 'haar', 0);
    hadper_trans_single_den{h}         = single(Tfor3rd);
    inverse_hadper_trans_single_den{h} = single(Tinv3rd);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Create Kaiser Window.

Wwin2D    = kaiser(N1, 2) * kaiser(N1, 2)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group1=cell(1,15000);      %% Predefine All Groups.
group2=cell(1,15000);
Coord=zeros(80000,4);     %% Predefine Image Block Coordinate.
weight = zeros(1,10000);  %% Predefine Weight Value.
for i=1:15000
    group1{i}=zeros(N1,N1,N2); %% Predefine Each Group.
    group2{i}=zeros(N1,N1,N2);
end

im_p1=single(padarray(im,[40 40],'symmetric','both')); 
im_p2=single(padarray(imout,[40 40],'symmetric','both'));%% Image Extension.
[M,N]=size(im_p1);                              %% Image Size.                             
 x = zeros(1000,1);                            %% Image Block Horizental Coordinate in Each group.
 y = zeros(1000,1);                            %% Image Block Vertical Coordinate in Each group. 

 l=1;
 c=1;
 
 for a = 41:Nstep:M1+40                              %% Block Matching 3D Transform Image Denoising Main Function
    for b = 41:Nstep:M2+40
     
       k=1;                                    %% 1st Image Block Index.     
       im_n1=im_p1(a:a+N1-1,b:b+N1-1);           %% 1st Image Block,i.e.,Reference Image Block.
       im_n2=im_p2(a:a+N1-1,b:b+N1-1); 

                   group1{l}(:,:,k) = im_n1; 
                   group2{l}(:,:,k) = im_n2;%% lth Group and kth Block, Note that l = 1, k = 1 here.
                   Coord(c,1)=l;               %% 1st Column represents lth Group in Array Coord. 
                   Coord(c,2)=1;               %% Second Column represents 1st block in Array Coord. 
                   Coord(c,3)=a;               %% Third Column represents Horizental Coordinate in Array Coord.
                   Coord(c,4)=b;               %% 4th Column represents Vertical Coordinate in Array Coord.
                   c = c+1; 
                   k = k+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Create Block Matching Neighbor, The size of Neighbor is Ns*Ns.
               
                Q = 1;
                for r = 1:(Ns-1)/2  %% (Ns-1)/2 is Neighbor Radius.  
                      for i = -r:r-1
                        x(Q) = a-r;
                        y(Q) = b+i;
                        Q = Q+1;
                      end   
                      for i = -r:r-1
                        x(Q) = a+i;
                        y(Q) = b+r;
                        Q = Q+1;
                      end   
                      for i = -r:r-1
                        x(Q) = a+r;
                        y(Q) = b-i;
                        Q = Q+1;
                      end   
                      for i = -r:r-1
                        x(Q) = a-i;
                        y(Q) = b-r;
                        Q = Q+1;
                      end 
                end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Block Matching Procedure.
             dis1=zeros(Ns^2-1,3);
             for R=1:length(x)
               ux = im_p2(x(R):x(R)+N1-1,y(R):y(R)+N1-1);   %% Choosing Block in The Neighborhood to be Matched.
                dis1(R,3)=sum(sum((ux-im_n2).*(ux-im_n2)));  %% Computing Euclidean Distance between Reference Block The Blcok to be Matched.
                dis1(R,1)=x(R);                            %% Store Block Horizental Coordinate.
                dis1(R,2)=y(R);                            %% Store Block Vertical Coordinate. 
             end
             
              dis2=sortrows(dis1,3);                       %% Sorting Array Dis1 by the Third Column Index, i.e., The Euclidean Distances.
              
              for R=1:N2-1                                 %% Selecting 15 Blocks and Store them to the group as the second block to 16th block.
                   ux1 = im_p1(dis2(R,1):dis2(R,1)+N1-1,dis2(R,2):dis2(R,2)+N1-1);
                   ux2 = im_p2(dis2(R,1):dis2(R,1)+N1-1,dis2(R,2):dis2(R,2)+N1-1);
                   group1{l}(:,:,k) = ux1;    %% Note that l = 1, k = 2 here.
                   group2{l}(:,:,k) = ux2;
                   Coord(c,1)=l;                            %% Note that l = 2, c = 2 here.
                   Coord(c,2)=k;
                   Coord(c,3)=dis2(R,1);                    %% The Third Column represents matched block Horizental Coordinate in Array Coord.
                   Coord(c,4)=dis2(R,2);                    %% The 4th Column represents matched block Vertical Coordinate in Array Coord.
                   c=c+1;
                   k=k+1;
              end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  %% 2D Block Forward Transform
  for K=1:N2
        group1{l}(:,:,K)=Tfor*group1{l}(:,:,K)*Tfor';
        group2{l}(:,:,K)=Tfor*group2{l}(:,:,K)*Tfor'; 
   end 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 %% The Third Dimension Forward Transform and Hard Thresholding. 
 for j=1:N1
           
           d31=zeros(N2,N1); 
           d32=zeros(N2,N1); 

           for u=1:N2
                    d31(u,:)=group1{l}(:,j,u); %% Each Column in Group gives a temp column vector d3.
                    d32(u,:)=group2{l}(:,j,u);
           end 

                 d31 = hadper_trans_single_den{h}*d31; 
                 d32 = hadper_trans_single_den{h}*d32;%% The Third Dimension Forward Transform.
                 
                 d31 = ((d32.^2)./(d32.^2+(sigma/255)^2)).*d31; %% The Hard Thresholding.
                 
           for u=1:N2
                group1{l}(:,j,u) =d31(u,:);  %% Transformed Column Vector are replaced to Group .  
           end 
           for u=1:N2
                group2{l}(:,j,u) =d32(u,:);  %% Transformed Column Vector are replaced to Group .  
           end   
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
 %% Computing Weight Value in Each Group.

    weight(l)=1/((norm((group2{l}(:).^2)./((group2{l}(:).^2)+(sigma/255)^2))^2)*(sigma/255)^2);
        
     
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 %% The Third Dimension Inverse Transform.    
  for j=1:N1

           d31=zeros(N2,N1); 

           for u=1:N2
               d31(u,:)=group1{l}(:,j,u);
           end 

               d31=inverse_hadper_trans_single_den{h}*d31;
               
           for u=1:N2
                group1{l}(:,j,u) =d31(u,:);
           end   
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %% The 2D Block Inverse Transform.      
   for K=1:N2
        group1{l}(:,:,K)=Tinv*group1{l}(:,:,K)*Tinv';
   end 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
              l=l+1;  %% Implement next Group
                
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Final Aggregation Processing
imr=zeros(M,N); %% Create the same size all zeros array as the extended image.
W=zeros(M,N);   %% Create the same size all zeros array as the extended image as the weight array.

C=length(Coord(:,1));  %% Check the number of group.
for r=1:C
      if Coord(r,1)==0
         break;
      end
      %% Group Aggregation
      imr(Coord(r,3):Coord(r,3)+N1-1,Coord(r,4):Coord(r,4)+N1-1)=imr(Coord(r,3):Coord(r,3)+N1-1,Coord(r,4):Coord(r,4)+N1-1)+group1{Coord(r,1)}(:,:,Coord(r,2))*weight(Coord(r,1)).*Wwin2D; 
      %% Weight Aggregation
      W(Coord(r,3):Coord(r,3)+N1-1,Coord(r,4):Coord(r,4)+N1-1)=W(Coord(r,3):Coord(r,3)+N1-1,Coord(r,4):Coord(r,4)+N1-1)+weight(Coord(r,1))*Wwin2D; 
end


imr=imcrop(imr,[41 41 M2-1 M1-1]); %% Cropping effective parts from the extended image.
W=imcrop(W,[41 41 M2-1 M1-1]);     %% Cropping effective parts from the extended weight array.

imout=imr./W;                      %% Mean of the repetitive pixels.              






 