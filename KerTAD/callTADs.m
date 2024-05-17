%{

Main function for calling TADs. Below we describe the only required
parameter as well as two hard-coded parameters that more advanced users can
modify depending on their use cases.


Required Parameter:

--Input Map: Input maps should be in the form of a symmetric matrix with nonnegative
values (elements of the matrix do not have to be integer counts). This variable is required.

Optional Parameters: These parameters can be manually edited in the source
code for advanced users. 

--Gamma: This parameter determines how many times each initial segmentation
of the input map is divided into smaller maps. By default, gamma=2
and in most cases this should work well. For very large or heterogenous maps
gamma can be set higher and for "well-behaved" maps gamma should be set to 1.
For the paper we used gamma=2 for all experimental maps and gamma=1 for
synthetic maps. 

--Kappa: This parameter determines how many TAD corners are allowed in each
row (so if kappa=n then the first n TAD corner points starting from the diagonal and
moving right are kept for each row). Kappa=3 by default and this value is used
throughout the paper except for simple maps where kappa=1 (this is done to
restrict KerTAD to working like simple TAD callers to
level the playing field). 

Outputs: callTADs returns 2 variables:

 --output_map: an NxN binary matrix where every 1 represents a TAD corner in the upper triangular matrix
 --TAD_list: a list of TAD corner locations (upper triangular matrix) in column format.

%}

function [output_map, TAD_list] =callTADs(input_map)
 

%gamma is set to 2 by default. For paper gamma=1 for synthetic and =2 for experimental maps. 
%See documentation at top of function. When in doubt keep at default value.
gamma=2;   

%Kappa=3 by default. See documentation at top of function. When in doubt keep at default value.
kappa=3;

%preprocess input_map
map=preprocessing(input_map); 


%Extremely sparse maps can throw errors for the masks thus we use a 
%Gaussian for high levels of sparsity to avoid errors
sp_level=(size(input_map,1).^2-nnz(input_map))/(size(input_map,1).^2);
if((sp_level>0))
    map=imgaussfilt(map,sp_level*1.5);   
end
%}
 
%Find outliers in preprocessed map. If looking at extremely sparse data,
%this can be turned off
n=5;
ratio_zeros=zeros(size(input_map,1),1); %ratio of n elements that are 0 for each row of input map
input_map2=double(input_map==0);
for i=1:round(size(input_map,1)/2)
    ratio_zeros(i,:)=sum(input_map2(i,i:i+n-1))./n;
end
for i=round(size(input_map,1)/2)+1:size(input_map,1)
    ratio_zeros(i,:)=sum(input_map2(i,i-n+1:i))./n;
end
TF2=(ratio_zeros>.8);

%{
TF3 = isoutlier(diag(map),"grubbs");  %alternatively use preprocessing(input_map)
TF4=TF3.*(diag(map)<median(diag(map))); %only include outliers that are due to low counts
TF2=(TF+TF4)>0;
%}

%Segment maps into smaller maps
pos=1; 
min_TAD_size=2*gamma;
indices=[];
for k=1:size(TF2)
    if(pos==1)
        if(TF2(k)==0)
           indices=[indices;k];
           pos=2;
        end
    end
    if(pos==2)
        if(TF2(k)==1)
           indices=[indices;k];
           pos=1;
        end
    end
end
%add last index if odd number
if(mod(size(indices),2)==1)
    indices=[indices; size(map,1)];
end
index_list=[indices(1:2:end) indices(2:2:end)];
index_list=index_list(abs(index_list(:,1)-index_list(:,2))>min_TAD_size,:); %Ignore very small segementations which might cause errors
num_maps=size(index_list,1); 
 


fullmap=zeros(size(input_map));
for k=1:num_maps
    
    map_segment=map(index_list(k,1):index_list(k,2),index_list(k,1):index_list(k,2));  %segment of original map being analyzed
    binary_map=zeros(size(map_segment));
    %further split map_segment into smaller maps based on gamma. 
    for j=1:gamma
        cut=round(size(map_segment,1)./j);  
        hmap=zeros(size(map_segment)); %map of 0s to fill in
        index=1;
        for i=1:(ceil(size(map_segment,1)/cut))
            if(index+cut<size(map_segment,1))
                tmap=map_segment(index:index+cut,index:index+cut);
                hmap(index:index+cut,index:index+cut)=mask_1(tmap).*mask_2(tmap);
                index=index+cut+1;
            else 
                tmap=map_segment(index:end,index:end);
                hmap(index:end,index:end)=mask_1(tmap).*mask_2(tmap);
            end
        end
        binary_map=binary_map+hmap;
    end
    binary_map=double(binary_map>0);
    fullmap(index_list(k,1):index_list(k,2),index_list(k,1):index_list(k,2))=binary_map;
 
end
 

final_map=zeros(size(fullmap));
for k=1:size(fullmap,1)
    count=0;
    for l=k:size(fullmap,2)
        if(fullmap(k,l)==1)
            count=count+1;
            if(count>kappa)
                break
            end
            final_map(k,l)=1;
         end
    end    
end



%Remove TAD calls very close to the diagonal (first superdiagonal) since these are not TADs
output_map=final_map.*(~(diag(diag(ones(size(final_map))))+diag(diag(ones(size(final_map)),1),1))); 


[row,col] = ind2sub(size(input_map),find(output_map));
TAD_list=[row col];
end