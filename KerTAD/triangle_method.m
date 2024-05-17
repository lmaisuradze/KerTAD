%Triangle method code for thresholding.
function [binary_image, threshold] = triangle_method (input_map)
input_map=input_map./max(input_map(:)); %normalize
hst=imhist(input_map);
[y,x]=max(hst(:));
p=find(hst>y/2e4); 
hst=fliplr(hst');
i=0:(p(end)-x);
j=hst(i+(257-p(end)));
l=(j+i/(y/(p(end)-x)))/((y/(p(end)-x))+1/(y/(p(end)-x)));
m=((((y/(p(end)-x))*l)-j).^2+(l-i).^2).^0.5;
threshold=find(max(m)==m);
threshold=(257-p(end))+mean(threshold);
threshold=(257-threshold)/256;
binary_image=imbinarize(input_map,threshold);
end

