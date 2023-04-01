import numpy as np
import os
import matplotlib.pyplot as plt
X=[]
Y=[]
lx=0.05
ly=0.05
dir=r'D:\\GIT_projects\003_2DHeat_diffusion_transient'
for path in os.listdir(dir):
    if not path.endswith(".out"):
        continue
    T=np.genfromtxt(path, dtype=None)    
    X=np.linspace(0, lx,T.shape[0])
    Y=np.linspace(0, ly,T.shape[1])
    h=plt.contourf(X,Y,T)
    plt.colorbar(h)    
    foo=path.split(".out")[0]
    fname=foo+".png"
    plt.savefig(fname)
    plt.clf()
    
    