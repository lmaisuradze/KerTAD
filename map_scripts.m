%Script for making and saving the 100 complex maps with no added noise.
%Only use this if you want to create your own maps, the ones used in the
%paper are in the "Normal Maps" folder.

%Uncomment each section to utilize

%{
for i=1:100
    [X,Y,Z]= generate_complex_map(0,0);
    d1=strcat("C",num2str(i));
    d2=strcat("G",num2str(i));
    writematrix(X,d1,'Delimiter','tab');
    writematrix(Z,d2,'Delimiter','tab')
end
%}

%Script for making and saving noisy maps. Only use this if you want to create your own maps, the ones used in the
%paper are in the "Noisy Maps" folder.
%{
for i=1:10
    [X,Y,Z]= generate_complex_map(0,0);
    d2=strcat("G",num2str(i));
    writematrix(Z,d2,'Delimiter','tab')
    index=1;
    for j=0:.05:1
        Nx=add_noise(X,j);
        X2=X+Nx;
        d1=strcat("N_",num2str(i),"_",num2str(index));
        writematrix(X2,d1,'Delimiter','tab');
        index=index+1;
    end
end



function Nx = add_noise(final_map,noise_level)  %Set this variable to be between 0 to 1.
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
end
%}


%Script for making and saving sparse maps. Only use this if you want to create your own maps, the ones used in the
%paper are in the "Sparse Maps" folder.
%{
for i=1:10
    [X,Y,Z]= generate_complex_map(0,0);
    d2=strcat("G",num2str(i));
    writematrix(Z,d2,'Delimiter','tab')
    index=1;
    for j=0:.05:.95
        X2=add_sparsity(X,j);
        d1=strcat("S_",num2str(i),"_",num2str(index));
        writematrix(X2,d1,'Delimiter','tab');
        index=index+1;
    end
end


function sparse_map = add_sparsity(final_map,sparsity_level)  %Set this variable to be between 0 to 1.
    no_deletes=round(size(final_map,1)^2*sparsity_level);
    ix=randperm(size(final_map,1)^2,no_deletes);
    sparse_map=1:size(final_map,1)^2;
    sparse_map2=reshape(sparse_map.',[size(final_map,1) size(final_map,1)]).';
    iy=ismember(sparse_map2,ix);
    final_map(iy)=0;
    sparse_map=final_map;
end
%}