import numpy as np
M=24; N=1<<M
p=np.zeros(N,dtype=np.int64)
for u in range(1,N):
    w=u>>1; p[u]=p[w] if u&1 else p[w]+w+(w&1)
vals=np.unique(1+p)
for B in (2**10,2**14,2**18,2**21):
    odd=np.arange(1,B,2); miss=odd[~np.isin(odd,vals)]
    print('one-apex transitive tapes (tape length <= %d): odd N < %d not of the form 1 + p(u): %d (%.1f%%); first misses %s'%(M,B,len(miss),100*len(miss)/len(odd),miss[:12].tolist()))
