import numpy as np, glob
V10=np.loadtxt('pairs_10.txt',dtype=np.int64)
S12=np.unique(np.concatenate([np.loadtxt(f,dtype=np.int64,usecols=(0,1)) for f in sorted(glob.glob('s12_*.txt'))]),axis=0)
print('order-12 distinct sampled pairs',len(S12))
for (al,be,name,start) in [(3,1,'K1',36835),(5,3,'TT2',188269)]:
    v10=np.unique(al*V10[:,0]+be*V10[:,1]); v12=np.unique(al*S12[:,0]+be*S12[:,1])
    for r in (1,2,0):
        U=np.unique(np.concatenate([v10,v12]))
        N=np.arange(start, U[-1]+1); N=N[(N%2==1)&(N%3==r)]
        hit=np.isin(N,U)
        miss=N[~hit]
        first=int(miss[0]) if len(miss) else None
        # order-12 alone: its longest run in this residue class
        h12=np.isin(N,v12); idx=np.nonzero(~h12)[0]; e=np.concatenate(([-1],idx,[len(N)])); g=np.diff(e)-1; k=int(np.argmax(g))
        print('%s, N = %d mod 3: orders 10 (exact) + 12 (sampled) cover every odd N from %d up to %s; order-12 longest run [%d, %d]'
              %(name,r,start,(first-6) if first else U[-1],int(N[e[k]+1]),int(N[e[k+1]-1])))
