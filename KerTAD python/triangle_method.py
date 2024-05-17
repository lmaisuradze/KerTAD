import numpy as np
#calculates the triangle method threshold
def triangle_method(input_map):
    #normalize input map
    input_map=input_map/np.max(input_map); 
    input_map[input_map<0]=0
    #bin the values of input map into 256 bins (this code ends up being equivalent to imhist in Matlab)
    bins = np.linspace(0, 1, 256)
    bins = bins -bins[1]/2 # offset the bins
    hst=np.zeros([256,1])
    hst[0:256-1,0] =np.histogram(input_map, bins)[0]
    hst[255,0]=(input_map.shape[1] **2) - np.sum(np.histogram(input_map, bins)[0])
    y=np.max(hst)
    x=np.argmax(hst)
    
    
    a=(hst>y/2e4)*1
    p=np.flatnonzero(a)
  
    
    hst=np.fliplr(hst.T)
    i=np.arange((p[-1]-x+1)).reshape(1, -1)
    ind=(i+(257-p[-1])-2)
    j=hst.T[np.asarray(ind)]
    j=np.reshape(j, [1 ,-1])
    b=y/(p[-1]-x)
    l=(j+i/b)/(b+1/b)
    m=((((y/(p[-1]-x))*l)-j)**2+(l-i)**2)**0.5
    t= ((np.max(m)==m)*1)
    threshold=np.nonzero(np.max(m)==m)[:][1]
    threshold=(257-p[-1])+np.mean(threshold)
    threshold=(257-threshold)/256
    binary_image=input_map>threshold *1
    return binary_image,threshold
    