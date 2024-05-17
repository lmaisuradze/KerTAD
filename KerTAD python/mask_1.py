import numpy as np
from triangle_method import *
#calculate first mask for an input matrix and return NxN binary matrix
def mask_1(input_map):
    #Calculate partial derivative using slicing (alternatively you can use scipy ndimage)
    Gy=np.zeros(input_map.shape)
    Gy[0:input_map.shape[1]-1,:]=input_map[1:input_map.shape[1],:] - input_map[0:-1,:]
    x=np.tile(np.ptp(Gy,1,keepdims=True), (1, Gy.shape[1]))
    x2=np.tile(np.sum(np.absolute(Gy),1,keepdims=True), (1, Gy.shape[1]))
    S=x*x.T* x2*x2.T
    #Go through diagonal of S and keep only local neighborhood maxima
    S2=np.zeros(S.shape)
    for i in range(1, S.shape[1]-1):
         nhood=S[i-1:i+2,i-1:i+2]
         if(S[i,i]==np.max(nhood)):
            S2[i,i]=1
         
    S2[0,:]=1
    S2[S2.shape[1]-1,:]=1

    binary_mask, threshold = triangle_method(S)
    S=S/np.max(S)
    M=((np.diag(S)>threshold)).reshape(-1, 1)*S2   
    M[0,0]=1;   
    M[-1,-1]=1;   
    d=np.expand_dims(np.diag(M),axis=1)
    L=np.tile(d, input_map.shape[0])
    M2=L+L.T
    empty_map=np.triu(M2>1)
    output_map=empty_map*((~np.eye(input_map.shape[0],dtype=bool)))
    first_mask=np.append(np.expand_dims(output_map[0,:],axis=0),np.zeros((1,output_map.shape[0])),axis=0)
    output_map=np.append(first_mask,output_map[1:-1,:],axis=0)
    return output_map