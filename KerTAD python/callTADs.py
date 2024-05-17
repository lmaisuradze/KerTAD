import math
import sys
import numpy as np
from preprocessing import *
from mask_1 import *
from mask_2 import *
from scipy.ndimage import gaussian_filter
 
# Main function for calling TADs. Below we describe the only required
# parameter as well as two hard-coded parameters that more advanced users can
# modify depending on their use cases.


# Required Parameter:

# --Input Map: Input maps should be in the form of a symmetric matrix with nonnegative
# values (elements of the matrix do not have to be integer counts). This variable is required.

# Optional Parameters: These parameters can be manually edited in the source
# code for advanced users. 

# --Gamma: This parameter determines how many times each initial segmentation
# of the input map is divided into smaller maps. By default, gamma=2
# and in most cases this should work well. For very large or heterogenous maps
# gamma can be set higher and for "well-behaved" maps gamma should be set to 1.
# For the paper we used gamma=2 for all experimental maps and gamma=1 for
# synthetic maps. 

# --Kappa: This parameter determines how many TAD corners are allowed in each
# row (so if kappa=n then the first n TAD corner points starting from the diagonal and
# moving right are kept for each row). Kappa=3 by default and this value is used
# throughout the paper except for simple maps where kappa=1 (this is done to
# restrict KerTAD to working like simple TAD callers to
# level the playing field). 

# Outputs:

#  --TAD_list: a list of TAD corner locations (upper triangular matrix) in column format.
#  --output_map: an NxN binary matrix where every 1 represents a TAD corner in the upper triangular matrix

 

 

def callTADs(input_map):

    #gamma is set to 2 by default. For paper gamma=1 for synthetic and =2 for experimental maps. 
    #See documentation at top of function. When in doubt keep at default value.
    gamma=2;  

    #Kappa=3 by default. See documentation at top of function. When in doubt keep at default value.
    kappa=3

    original_map= np.copy(input_map) #use copy to avoid overwriting
    #preprocess input_map
    map=preprocessing(input_map)
    #Extremely sparse maps can throw errors for the masks thus we use a 
    #Gaussian for high levels of sparsity to avoid errors
    sp_level=np.count_nonzero(original_map==0)/(original_map.shape[0]**2)
    if(sp_level>0):
        map=gaussian_filter(map.astype(float), sigma=sp_level*1.5, mode = 'nearest', truncate=2.0)
    
    #Find outliers in preprocessed map using grubbs method
    n=5
    ratio_zeros=np.zeros([original_map.shape[0],1]) #ratio of n elements that are 0 for each row of original input
    input_map2=np.zeros(original_map.shape[0])+(original_map==0)*1
    r=int(np.round(original_map.shape[0]/2))
    for i in range(0, r):
        ratio_zeros[i,:]=np.sum(input_map2[i,i:i+n])/n
    for i in range(r, original_map.shape[0]):
        ratio_zeros[i,:]=np.sum(input_map2[i,(i-n)+1:i+1])/n
        
    TF=(ratio_zeros>.8)*1

    #Grubbs outliers not implemented yet, will add shortly
    TF2=TF

    pos=1
    min_TAD_size=2*gamma
    indices=[]
    for k in range(0,TF2.shape[0]):
        if(pos==1):
            if(TF2[k]==0):
                indices=np.append(indices,[k],axis=0)
                pos=2
        if(pos==2):
            if(TF2[k]==1):
                indices=np.append(indices,[k],axis=0)
                pos=1
    indices = np.array(indices)          
        
        
    #add a last index if there are an odd number 
    if(np.mod(indices.shape[0],2)==1):
        indices=np.append(indices,[original_map.shape[0]-1],axis=0)
    #turn indices into two column format
    index_list= np.append(np.expand_dims(indices[0:indices.shape[0]:2],axis=1),np.expand_dims(indices[1:indices.shape[0]:2],axis=1),axis=1)
    #get rid of index pairs less than the minimum TAD size allowed
    index_list=index_list[np.abs(index_list[:,0]-index_list[:,1])>min_TAD_size,:] #Ignore very small segmentations which might cause errors
    #count number of maps in index list
    num_maps=index_list.shape[0] 
 
    fullmap=np.zeros(original_map.shape)
    for k in range(0,num_maps):        
        map_segment=map[int(index_list[k,0]):int(index_list[k,1]+1),int(index_list[k,0]):int(index_list[k,1]+1)];  #segment of the map being analyzed
        binary_map=np.zeros(map_segment.shape)
        #further split map_segment into smaller maps based on gamma. 
        for j in range(1,gamma+1):  
            cut=int(np.floor(map_segment.shape[0]/j + 0.5))
            hmap=np.zeros(map_segment.shape) #map of 0s to fill in
            index=0
            for i in range(1,math.ceil((map_segment.shape[0])/cut)+1):
                if((index+cut)<map_segment.shape[0]):
                    tmap=map_segment[int(index):int(index+cut)+1,int(index):int(index+cut)+1]
                    hmap[int(index):int(index+cut)+1,int(index):int(index+cut)+1]=mask_1(tmap)*mask_2(tmap)
                    index=index+cut+1
                    
                else:
                    tmap=map_segment[int(index):int(map_segment.shape[0]),int(index):int(map_segment.shape[0])]
                     
                    hmap[int(index):int(map_segment.shape[0]),int(index):int(map_segment.shape[0])]=mask_1(tmap)*mask_2(tmap)
                
            
            binary_map=binary_map+hmap
        
        binary_map=(binary_map>0)*1
        fullmap[int(index_list[k,0]):int(index_list[k,1]+1),int(index_list[k,0]):int(index_list[k,1]+1)]=binary_map
        
    
    final_map=np.zeros(fullmap.shape)
    for k in range(0,fullmap.shape[0]):  
        count=0
        for l in range(k,fullmap.shape[1]):
            if(fullmap[k,l]==1):
                count=count+1
                if(count>kappa):
                    break
                
                final_map[k,l]=1

    #Remove TAD calls very close to the diagonal (first superdiagonal) since these are not TADs
    output_map=final_map * ( ~np.array(((np.diag(np.diag(np.ones(final_map.shape)))+np.diag(np.diag(np.ones(final_map.shape),1),1))), dtype=bool))
    
    #convert to list format (remember that Python indices start at 0 not 1 like Matlab. The ordering may be different as well)
    TAD_list=np.argwhere(output_map)
    return TAD_list, output_map
    




#Load in input file and output file names through command line arg
inFile = sys.argv[1]
outFile = sys.argv[2]
#load in map (should be in same folder if path not specified), run KerTAD, and save to output file specified (ideally .txt)
A=np.loadtxt(inFile)
TAD_list,output_map=callTADs(A)
np.savetxt(outFile, TAD_list) #if you want the binary map of TAD corner point locations, change 'TAD_list' to 'output_map'