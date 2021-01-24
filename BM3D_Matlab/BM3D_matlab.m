function [imout] = BM3D_matlab(im,sigma)

[M1,M2]=size(im); %% Get the input image size.
if sigma < 40
N1     = 8;     %% Block size.
N2     = 16;    %% The number of block in a group.
Ns     = 39;    %% The size of block-matching neighborhood. 
Nstep  = 3;
Thr    = 2.7;
else
N1     = 8;     %% Block size.
N2     = 32;    %% The number of block in a group.
Ns     = 39;    %% The size of block-matching neighborhood. 
Nstep  = 4;
Thr    = 2.8;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Create 2D block transform matrices.
%%%%
[Tfor, Tinv]   = getTransfMatrix(N1, 'bior1.5', 0);     %% get (normalized) forward and inverse transform matrices
Tfor=single(Tfor);
Tinv=single(Tinv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Create third dimension transform matrices.
%%%%
hadper_trans_single_den = cell(1,1000);
inverse_hadper_trans_single_den = cell(1,1000);
for hpow = 0:ceil(log2(N2))
    h = 2^hpow;
    [Tfor3rd, Tinv3rd]  = getTransfMatrix(h, 'haar', 0);
    hadper_trans_single_den{h}  = single(Tfor3rd);
    inverse_hadper_trans_single_den{h} = single(Tinv3rd);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Create Kaiser window.

Wwin2D  = kaiser(N1, 2) * kaiser(N1, 2)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group = cell(1,15000);      %% Predefine all groups.
Coord = zeros(80000,4);     %% Predefine image block coordinates.
weight = zeros(1,10000);  %% Predefine weight values.

for i=1:15000
    group{i} = zeros(N1,N1,N2); %% Predefine each group,i.e.,3D array.
end


im_p = single(padarray(im,[40 40],'symmetric','both'));  %% Image extension.
[M,N] = size(im_p);                              %% Get the extended image size.                             
 x = zeros(1000,1);                            %% It is used to store image block horizontal coordinates in each group.
 y = zeros(1000,1);                            %% It is used to store image block vertical coordinates in each group. 

 L=1;
 c=1;
 
 for a = 41:Nstep:M1+40                              %% Block matching 3D transform image denoising main function
    for b = 41:Nstep:M2+40
     
       k=1;                                    %% 1st image block index.     
       im_n=im_p(a:a+N1-1,b:b+N1-1);           %% 1st image block,i.e.,reference image block.

                   group{L}(:,:,k) = im_n;     %% Lth group and kth block, note that L = 1, k = 1 here.
                   Coord(c,1) = L;             %% The 1st column represents Lth group in array Coord. 
                   Coord(c,2) = k;             %% The 2nd column represents kth block in array Coord. 
                   Coord(c,3) = a;             %% The 3rd column represents horizontal coordinate in array Coord.
                   Coord(c,4) = b;             %% The 4th column represents vertical coordinate in array Coord.
                   c = c+1; 
                   k = k+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Create block matching neighborhood, the size of neighborhood is Ns*Ns-1.
               
                Q = 1;
                for r = 1:(Ns-1)/2  %% (Ns-1)/2 is the neighborhood radius.  
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
               ux = im_p(x(R):x(R)+N1-1,y(R):y(R)+N1-1);   %% Choosing blocks in the neighborhood to be matched.
                dis1(R,3)=sum(sum((ux-im_n).*(ux-im_n)));  %% Computing Euclidean distance between reference block and the blcok to be matched.
                dis1(R,1)=x(R);                            %% Store block horizontal coordinate.
                dis1(R,2)=y(R);                            %% Store block vertical coordinate. 
             end
             
              dis2=sortrows(dis1,3);                       %% Sorting array dis1 by the third column index, i.e., The Euclidean distances.
              
              for R=1:N2-1                                 %% Selecting 15 most similar blocks to the reference block and store them to the group as the second block to 16th block.
                   ux = im_p(dis2(R,1):dis2(R,1)+N1-1,dis2(R,2):dis2(R,2)+N1-1);
                   group{L}(:,:,k) = ux;                    %% Note that L = 1, k = 2 here.
                   Coord(c,1)=L;                            %% Note that L = 1, c = 2 here.
                   Coord(c,2)=k;
                   Coord(c,3)=dis2(R,1);                    %% The 3rd column represents matched block horizontal coordinate in array coord.
                   Coord(c,4)=dis2(R,2);                    %% The 4th column represents matched block vertical coordinate in array coord.
                   c=c+1;
                   k=k+1;
              end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  %% 2D block forward transform
  for K=1:N2
        group{L}(:,:,K)=Tfor*group{L}(:,:,K)*Tfor';
   end 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 %% The third dimension forward transform, hard thresholding and inverse transform. 
 
 Non_zero = 0;
 
 for j=1:N1
    for u=1:N1
                 d3=zeros(N2,1); 
               
                 d3(:,1)=group{L}(j,u,:); %% Each column in the group is assigned to a temporary column vector d3.
           

                 d3 = hadper_trans_single_den{h}*d3; %% The 3rd dimension forward transform.
                 
                 d3 = d3.*(abs(d3)>=Thr*(sigma/255)); %% The hard thresholding.
                 Non_zero = Non_zero + nnz(d3);
                 d3 = inverse_hadper_trans_single_den{h}*d3; %% The 3rd dimension inverse transform.

                group{L}(j,u,:) =d3(:,1);  %% Transformed column vector are replaced to group .  
    end   
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %% Computing weight value in each group.

if Non_zero >= 1
    weight(L) = 1/Non_zero;
else
     weight(L)=1;  
end  
    

 %% The 2D block inverse transform.      
   for K=1:N2
        group{L}(:,:,K)=Tinv*group{L}(:,:,K)*Tinv';
   end 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
              L = L+1;  %% Implement next group
                
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
      %% Group aggregation
      imr(Coord(r,3):Coord(r,3)+N1-1,Coord(r,4):Coord(r,4)+N1-1)=imr(Coord(r,3):Coord(r,3)+N1-1,Coord(r,4):Coord(r,4)+N1-1)+group{Coord(r,1)}(:,:,Coord(r,2))*weight(Coord(r,1)).*Wwin2D; 
      %% Weight aggregation
      W(Coord(r,3):Coord(r,3)+N1-1,Coord(r,4):Coord(r,4)+N1-1)=W(Coord(r,3):Coord(r,3)+N1-1,Coord(r,4):Coord(r,4)+N1-1)+weight(Coord(r,1)).*Wwin2D; 
end


imr=imcrop(imr,[41 41 M2-1 M1-1]); %% Cropping effective parts from the extended image.
W=imcrop(W,[41 41 M2-1 M1-1]);     %% Cropping effective parts from the extended weight array.

imout=imr./W;                      %% Mean of the repetitive pixels.              


