import numpy as np

#There are many ways to normalize/represent HiC maps. Unfortunately, these can significantly affect performance
#of TAD callers. We preprocess the input map in several steps to try to ensure consistency. 

def preprocessing(input_map):
    
    #First: we check that the map is in NxN format
    num_rows, num_cols = input_map.shape
    if(num_rows!=num_cols):
        raise Exception("Input map needs to be in NxN format") 
    
    #Next we make sure there are no negative,inf, or NaN values
    if (input_map.min()<0 or np.sum(np.isnan(input_map))>0 or np.sum(np.isinf(input_map))>0 ):
        raise Exception("Invalid elements in input map. Input has either negative, NaN or inf values") 
    
    #We now normalize the map. First, we set every element on the diagonal to be a maximum in its respective
    #row. 
    np.fill_diagonal(input_map, np.max(input_map,axis=1), wrap=False)

    #To ensure symmetry we take the upper triangular matrix of the 
    #input_map, flip it and set it as the lower triangular matrix (or vice versa)
    input_map = np.tril(input_map.transpose(),-1) + np.triu(input_map)


    #Finally we normalize by row
    for i in range(0, num_rows):
        if(input_map[i,i]!=0):
            input_map[i]= ((input_map[i,:]-np.mean(input_map[i,:]))/np.std(input_map[i,:], ddof=1))
    
    
    output_map=input_map
    return output_map


    
