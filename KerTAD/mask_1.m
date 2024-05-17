%Calculates the first mask for a given input map. Returns binary NxN matrix
function [output_map]=mask_1(input_map)

%calculate partial derivative using imfilter
h=[-1 1];
Gy=imfilter(input_map,h','symmetric');
x=repmat(range(Gy,2), 1, size(Gy,2));
x2=repmat(sum(abs(Gy),2), 1, size(abs(Gy),2));
S=x.*x.'.* x2.*x2.';

%Filter local maxima
S2=zeros(size(S));
for i=2:size(S2,1)-1
         pt=[i i];
         idy = max(pt(1)-1,1):min(pt(1)+1,size(S,1));
         idx = max(pt(2)-1,1):min(pt(2)+1,size(S,2));
         nhood = S(idy,idx);
         if(S(i,i)==max(nhood,[],'all'))
            S2(i,i)=1;
         end

end
S2(1,:)=1; S2(size(S2,1),:)=1;
 

[~, threshold] = triangle_method (S);
S=S./max(S,[],'all');
M=double(imbinarize(diag(S./max(S(:))),threshold)).*S2;

M(1,1)=1; 
M(size(M,1),size(M,1))=1;
M2=repmat(diag(M),1, size(input_map,1)) + repmat(diag(M),1, size(input_map,1)).';

empty_map=triu(M2>1);
output_map=empty_map.*(double(~eye(size(input_map))));
first_mask=[output_map(1,:);zeros(1,size(output_map,1));output_map(2:end-1,:)];
output_map= first_mask; 

end