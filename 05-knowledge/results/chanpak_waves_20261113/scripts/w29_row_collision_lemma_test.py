import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
import numpy as np, sys
from collections import Counter
L=int(sys.argv[1]); c,k,ph=int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])
N=1<<L
P=np.zeros(N,dtype=np.int64)
for u in range(1,N):
    P[u]= P[u>>1] if u&1 else P[u>>1]+(u>>1)+((u>>1)&1)
U=np.arange(0,N,dtype=np.int64)
G=c*U+k*P+ph*(U&1)
order=np.argsort(G,kind='stable'); Gs=G[order]; Us=U[order]
bounds=np.concatenate(([0],np.nonzero(Gs[1:]!=Gs[:-1])[0]+1,[len(Gs)]))
cnt=Counter()
for i in range(len(bounds)-1):
    a,b=bounds[i],bounds[i+1]
    if b-a<2: continue
    grp=Us[a:b]
    for x in range(len(grp)):
        for y in range(len(grp)):
            if x==y: continue
            u,v=int(grp[x]),int(grp[y])
            l=(u^v).bit_length(); w=u>>l
            cnt[(l,w)]+=1
worst=0; arg=None; tot=0
for (l,w),n in cnt.items():
    r=n*(1+w)/2**l; tot+=n
    if r>worst: worst,arg=r,(l,w,n)
print(f'row {c}U+{k}p+{ph}[odd], L={L}: total ordered collisions {tot} ({tot/N:.4f} per U); max over (l,w) of count*(1+w)/2^l = {worst:.4f} at (l,w,count)={arg}')
# profile: for each l, sum over w of count, and sum of count*(1+w)
for l in range(max(1,L-10),L+1):
    s=sum(n for (ll,w),n in cnt.items() if ll==l); s2=max([n*(1+w)/2**l for (ll,w),n in cnt.items() if ll==l] or [0])
    print(f'   l={l}: collisions {s} ({s/2**l:.4f} per 2^l), max_w count(1+w)/2^l = {s2:.4f}')
