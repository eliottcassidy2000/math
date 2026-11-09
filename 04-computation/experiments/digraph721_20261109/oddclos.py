# Odd-graph closure (exploration): O_n = d-values of odd graphs on n vertices; OW_n = ... without universal vertex.
# Bases: exhaustive odd graphs on <= 9 vertices (unions with >= 1 odd connected component).
# Rules: pad (O_{n-1} -> OW_n); cone over OW_{n-1} (odd iff the base is odd); K_2 v G for ANY G on n-2 vertices
# (adjacent twins): 3*2^(n-2) + x, x any odd < 2^(n-2); E_2 v (odd G); K_2 u G: 3 * (any value on n-2 vertices).
import sys, numpy as np
from collections import defaultdict
NB, NM = 9, 22
comps=defaultdict(list)
for line in open(sys.argv[1]):
    p=line.split(); n=int(p[0]); d=int(p[2]); par=p[3]
    if n>NB: continue
    e=p[5] if len(p)>5 else ''
    deg=[0]*n
    for i in range(0,len(e),2): deg[int(e[i],16)]+=1; deg[int(e[i+1],16)]+=1
    U = n==1 or any(x==n-1 for x in deg)
    comps[n].append((d,par=='O',U))
O={n:np.zeros(1<<(n-1),bool) for n in range(1,NM+1)}; OW={n:np.zeros(1<<(n-1),bool) for n in range(1,NM+1)}
clist=[(n,d,odd,U) for n in sorted(comps) for (d,odd,U) in comps[n]]
# unions: enumerate multisets of components (sizes sum <= NB); track (n, d, anyodd, U)
states=set([(0,1,False,False)])
frontier=[(0,1,False,0)]
seen=set()
def addstate(n,d,odd,U):
    if odd:
        O[n][(d-1)//2]=True
        if not U: OW[n][(d-1)//2]=True
# DP over component list index to avoid repeats
import itertools
def rec(start,n,d,odd,ncomp,U):
    if ncomp>0: addstate(n,d,odd, U if ncomp==1 else False)
    for idx in range(start,len(clist)):
        m,dc,oc,Uc=clist[idx]
        if n+m>NB: break
        rec(idx,n+m,d*dc,odd or oc,ncomp+1,Uc)
clist.sort()
rec(0,0,1,False,0,False)
def sh(dst,n,src,m,off):
    h=off//2; D=len(dst); 
    if h>=D: return
    seg=src[:max(0,min(len(src),D-h))]; dst[h:h+len(seg)] |= seg
for n in range(2,NM+1):
    O[n][:len(O[n-1])] |= O[n-1]; OW[n][:len(O[n-1])] |= O[n-1]; OW[n][:len(OW[n-1])] |= OW[n-1]
    if n>NB:
        allv=np.ones(1<<(n-3),bool) if n>=3 else None          # any graph on n-2 vertices: all odd < 2^(n-2)
        sh(O[n],n,allv,n-2,3<<(n-2))                            # K_2 v G  (has universal vertices)
        sh(OW[n],n,OW[n-2],n-2,(3<<(n-2))-2)                    # E_2 v (odd, no universal)
        sh(O[n],n,O[n-2],n-2,(3<<(n-2))-2)
        # K_2 u G: 3*x for x any odd < 2^(n-2) (G on n-2 vertices); no universal vertex
        idx=np.arange(1<<(n-3)); vals=3*(2*idx+1); ok=vals<(1<<n); 
        O[n][(vals[ok]-1)//2]=True; OW[n][(vals[ok]-1)//2]=True
        sh(O[n],n,OW[n-1],n-1,1<<(n-1))                         # cone over odd, no universal
        O[n] |= OW[n]
for n in range(3,NM+1):
    miss=[2*i+1 for i in np.nonzero(~O[n])[0]]
    print(n, 'O_n misses:', miss[:12] if len(miss)<=12 else (miss[:8],'...',len(miss)), ' max-missing-gap-from-top:', (1<<n)-1-max(2*i+1 for i in np.nonzero(O[n])[0]))
