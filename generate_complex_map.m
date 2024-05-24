%Generation of complex Hi-C maps for benchmarks. This code is a variation of the technique developed in diffHic (Lun et al)
%which was then implemented by Forcato et al. Here we use a variation of
%these methods. Only run this code if you want to generate new maps, for
%the maps used in the paper see the maps folders in "Complex Maps"

function [final_map,ground_truth,full_TAD_list] = generate_complex_map(noise_level,sparsity_level)
TAD_width_min=5;  %minimum width of TAD for starting layer
TAD_width_max=20; %maximum width of TAD for starting layer
no_TADs=100; %number of TADs
TAD_list=zeros(no_TADs,2); 
index=1;
for i=1:no_TADs
    TAD_list(i,:)=[index,index+randi([TAD_width_min,TAD_width_max])];
    index=TAD_list(i,2)+1;
end
map1=generate_map(TAD_list);
full_TAD_list=TAD_list;
%We now add layers to the map, i.e. nested and overlapping TADs
layers=3;
map2=zeros(size(map1,1),size(map1,1),layers);
for l=1:layers
    pct_removal=0.25;
    cuts=randperm(no_TADs-2)+1;
    cuts=cuts(1:round(pct_removal*no_TADs)).';
    TAD_list2=TAD_list;
    for i=1:size(cuts,1)
        TAD_list2(cuts(i),:)=[0,0];
    end
    TAD_list2=[nonzeros(TAD_list2(:,1)) nonzeros(TAD_list2(:,2))];
    for i=1:size(TAD_list2)-1
       if(TAD_list2(i,2)+1 ~= TAD_list2(i+1,1))
           TAD_list2(i+1,1)=TAD_list2(i,2)+1;
       end
    end
    map2(:,:,l)= generate_map(TAD_list2);
    full_TAD_list=[full_TAD_list;TAD_list2];
end
 
full_TAD_list=unique(full_TAD_list,'rows');
%Create final map and final TAD location mask (for ground truth)
ground_truth=zeros(size(map1));
for i=1:size(full_TAD_list,1)
    ground_truth(full_TAD_list(i,1),full_TAD_list(i,2))=1;
end
final_map=map1+sum(map2,3);



%%%Add noise
%In the original method impulse noise (set at an amplitude of 2) was added
%by sampling across the map (based on a power law decay) with replacement. We use 
%an amplitude of 5. 
if(noise_level>0)  %Set this variable to be between 0 to 1.
    noise_map=1:size(final_map,1)^2;
    x=repmat(1:size(final_map,1),size(final_map,1),1);
    y=x.';
    mu=28.*((abs(x-y)+1).^(-.69));  
    mu(logical(eye(size(mu))))=0;
    mu=mu(:);
    y = randsample(noise_map,round(size(final_map,1).^2*noise_level),true,mu);
    y=y.';
    noise_map2=reshape(noise_map.',[size(final_map,1) size(final_map,1)]).';
    Nx=zeros(size(final_map));
    
    while(size(y)>0)
        Nx=Nx+(double(ismember(noise_map2,y)).* 5);
        [~,ia,~] = unique(y,"first");
        A2=zeros(size(y));
        A2(ia)=1;
        y(logical(A2))=0;
        y=nonzeros(y);
    end
    final_map=final_map+Nx;
end

%%%Add sparsity
%We randomly (uniformly) set values of the Hi-C matrix to 0 depending on
%sparsity_level (the percent of elements of the matrix that are set to 0).
if(sparsity_level>0)  %Set this variable to be between 0 to 1.
    no_deletes=round(size(final_map,1)^2*sparsity_level);
    ix=randperm(size(final_map,1)^2,no_deletes);
    sparse_map=1:size(final_map,1)^2;
    sparse_map2=reshape(sparse_map.',[size(final_map,1) size(final_map,1)]).';
    iy=ismember(sparse_map2,ix);
    final_map(iy)=0;
end

end

 

function final_map=generate_map(TAD_list)
size_map=TAD_list(end,end);

x=repmat(1:size_map,size_map,1);
y=x.';



mu=28.*((abs(x-y)+1).^(-.69));  
 

TAD_map=zeros(size_map,size_map);
for i=1:size(TAD_list,1)
    TAD_map(TAD_list(i,1):TAD_list(i,2),TAD_list(i,1):TAD_list(i,2) )=1;
end

 
mu=mu.*TAD_map;
mu(logical(eye(size(mu))))=35;
dispersion=0.01;
variance=mu+dispersion.*(mu.^2);
n=(mu.^2)./(variance-mu);
p=mu./variance;

final_map=random('Negative Binomial',n,p);
final_map(isnan(final_map))=0;
end