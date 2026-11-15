import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
import numpy as np, sys
Lmax=int(sys.argv[1])
N=1<<Lmax
P=np.zeros(N,dtype=np.int64)
for u in range(1,N):
    P[u]= P[u>>1] if u&1 else P[u>>1]+(u>>1)+((u>>1)&1)
def rows(c,k,ph,label):
    print(label)
    for L in range(12,Lmax+1,2):
        U=np.arange(1,1<<L,dtype=np.int64)
        G=c*U+k*P[1:1<<L]+ph*(U&1)
        order=np.argsort(G,kind='stable'); Gs=G[order]; Us=U[order]
        # collision pairs: groups of equal values
        eq=np.nonzero(Gs[1:]==Gs[:-1])[0]
        # count ordered off-diagonal pairs via group sizes
        diff=np.diff(np.concatenate(([0],np.nonzero(Gs[1:]!=Gs[:-1])[0]+1,[len(Gs)])))
        offpairs=int((diff*(diff-1)).sum())
        # level profile: top differing bit of adjacent equal pairs (approx: first-neighbour pairs)
        a=Us[eq]; b=Us[eq+1]; lev=np.floor(np.log2((a^b).astype(np.float64))).astype(int)+1
        prof=np.bincount(L-lev,minlength=L+1)[:8]
        print(f'  L={L}: ordered collision pairs {offpairs}, /2^L = {offpairs/2**L:.4f}; adjacent-pair split levels (L-l = 0..7): {list(prof)}')
rows(18,3,2,'row 18U + 3p(U) + 2[U odd]  (3-cycle carrier, Y = 1)')
rows(189,11,0,'row 189U + 11p(U)  (mask-zero order-5 carrier, Y = 1)')
rows(2,1,0,'row 2U + p(U)  (small slope)')
