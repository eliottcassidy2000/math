# Residue equidistribution of the fixed-kernel statistic in tilted ensembles (an LLT aperiodicity check)
import numpy as np, glob, sys
from collections import Counter
n=int(sys.argv[1])
for f in sorted(glob.glob('s%d_*.txt'%n))[:4]:
    a=np.loadtxt(f,dtype=np.int64,usecols=(0,1))
    for (al,be,name) in [(3,1,'K1'),(5,3,'TT2')]:
        s=al*a[:,0]+be*a[:,1]
        line=[]
        for m in (3,5,7,9,16,25):
            c=np.bincount(s%m,minlength=m); adm=c[c>0]
            line.append('m=%d: %d classes, max/min %.3f'%(m,len(adm),adm.max()/adm.min()))
        print(f,name,len(s),'; '.join(line))
