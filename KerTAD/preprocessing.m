function output_map= preprocessing(input_map)

%There are many ways to normalize/represent HiC maps. Unfortunately, these can significantly affect performance
%of TAD callers. We preprocess the input map in several steps to try to ensure consistency. 
 
%First: we check that the map is in NxN format
if(size(input_map,1)~=size(input_map,2))
    error('Error: Input map needs to be in NxN format.')
end

%Next we make sure there are no negative,inf, or NaN values
if (min(input_map,[],'all')<0 || sum(double(isnan(input_map)),'all')>0 || sum(double(isinf(input_map)),'all')>0 )
    error('Error: Invalid elements in input map. Input has either negative, NaN or inf values')
end

 
%We now normalize the map. First, we set every element on the diagonal to be a maximum in its respective
%row. 

input_map(logical(eye(size(input_map))))=max((input_map),[],2);


%To ensure symmetry we take the upper triangular matrix of the 
%input_map, flip it and set it as the lower triangular matrix (or vice versa)

input_map = tril(input_map.',-1) + triu(input_map);
 
%Finally we normalize by row
for i=1:size(input_map,1)
    if(input_map(i,i)~=0)
        input_map(i,:)=(input_map(i,:)-mean(input_map(i,:)))./std(input_map(i,:));
    end
end
output_map=input_map;

end