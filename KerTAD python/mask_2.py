import numpy as np
from scipy.ndimage import gaussian_filter
from triangle_method import *
#calculate second mask for an input matrix and return NxN binary matrix
def mask_2(input_map):
    #calculate correlation matrix of columns in input_map for later
    corr=np.corrcoef(input_map.T)
    corr_map=(corr>0)*1
    corr_map=(gaussian_filter(corr_map.astype(float), sigma=0.5, mode = 'nearest', truncate=2.0) >0)*1
    #further preprocess the input map
    input_map=np.triu(input_map,1)
    input_map=np.diag(np.max(input_map,axis=1))+input_map

    for i in range(0,input_map.shape[0]):
        if(np.max(input_map[i,:]) !=0 ):
            input_map[i,:]=input_map[i,:]/np.max(input_map[i,:])
     
    input_map = np.tril(input_map.T,-1) + np.triu(input_map)

    #calculate partial derivative and further process with triangle method
    Gx=np.zeros(input_map.shape)
    Gx[:,0:input_map.shape[1]-1]=input_map[:,1:input_map.shape[1]]-input_map[:,0:input_map.shape[1]-1]
    Gx[:,-1]=input_map[:,0]-input_map[:,-1]
    Gx=(np.triu(Gx<0)* np.triu(Gx))  + (np.tril(Gx>0)*np.tril(Gx))
    Gx[np.isnan(Gx)]=0
    Gx=triangle_method(np.abs(Gx))[0]*Gx

    #calculate P by matrix multiplying Gx and Gx.T
    P=np.matmul(Gx,Gx.T)-np.matmul(Gx.T,Gx)
    P[np.isnan(P)]=0
    d=np.expand_dims(np.diag(P),axis=1)
    L=np.tile(d, P.shape[0])
    P2=-1*(L+L.T)*P
    P2=P2/np.max(P2)    
    P3=triangle_method(P2)[0]; 
    P3=(gaussian_filter(P3.astype(float), sigma=0.5, mode = 'nearest', truncate=2.0) >0)*1 
    seg1=P3[0:np.shape(P3)[0],0].reshape(1, -1)
    seg2=P3[0:np.shape(P3)[0]-1,0:np.shape(P3)[0]]
    P3=np.append(seg1,seg2,axis=0)
    output_map=P3*corr_map;   
    return output_map
