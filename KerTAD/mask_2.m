%Calculates the second mask for a given input map. Returns binary NxN matrix
function [output_map]=mask_2(input_map)

corr_map= double(imgaussfilt(double(corr(input_map)>0))>0); %corr matrix for each row, to use later
%further preprocessing on the input map
input_map=triu(input_map);
input_map(logical(eye(size(input_map))))=0;
input_map(logical(eye(size(input_map))))=max(input_map,[],2);
for i=1:size(input_map,1)
     input_map(i,:)=input_map(i,:)/max(max(input_map(i,:))); 
end
input_map(isnan(input_map))=0;
input_map(isinf(input_map))=0;
input_map = tril(input_map.',-1) + triu(input_map);

%calculate partial derivative using imfilter
h=[-1;1];
Gx = imfilter(input_map,h','circular');
Gx=(triu(Gx<0).*triu(Gx))  + (tril(Gx>0).*tril(Gx));
Gx(isnan(Gx))=0;
Gx=triangle_method(abs(Gx)).*Gx;


P=Gx*Gx.'-Gx.'*Gx;
P(isnan(P))=0;
P2=-1.*((repmat(double(diag(P)),1,size(P,1))+repmat(double(diag(P)),1,size(P,1)).')).*P;
P2=P2./max(P2,[],'all');
P3=triangle_method(P2);                                                                                                                                                                                              
P3=double(imgaussfilt(double(P3))>0);
P3=[P3(1:size(P3,1)); P3(1:size(P3,1)-1,1:size(P3,1))];
output_map=P3.*corr_map;

end