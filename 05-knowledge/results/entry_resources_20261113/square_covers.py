"""Minimum path-cover width and optimal endpoints in all first25 graphs."""
from pathlib import Path
from itertools import combinations
from math import isqrt
import numpy as np
from numba import njit
import json,gzip,sys,time
HERE=Path(__file__).resolve().parent
def need(p,m):
    if not p:raise RuntimeError(m)

@njit(cache=True)
def cover_dp(adj):
    n=len(adj);N=1<<n
    K=np.zeros(N,dtype=np.uint8);E=np.zeros(N,dtype=np.uint32)
    index=np.array([0,1,28,2,29,14,24,3,30,22,20,15,25,17,4,8,31,27,13,23,21,19,16,7,26,12,18,6,11,5,10,9])
    for mask in range(1,N):
        bits=mask;best=255;ends=0
        while bits:
            bit=bits&-bits;bits-=bit
            v=index[((bit*0x077CB531)&0xffffffff)>>27]
            smaller=mask^bit
            cost=int(K[smaller])+int((int(E[smaller])&int(adj[v]))==0)
            if cost<best:best=cost;ends=bit
            elif cost==best:ends|=bit
        K[mask]=best;E[mask]=ends
    return K,E

def reconstruct(mask,adj,K,E,end=None):
    if not mask:return []
    v=(int(E[mask])&-int(E[mask])).bit_length()-1 if end is None else end
    need(int(E[mask])>>v&1,'optimal ending vertex')
    smaller=mask^(1<<v)
    if not smaller:return [[v+1]]
    possible=int(E[smaller])&int(adj[v])
    if int(K[smaller])==int(K[mask]):
        need(possible!=0,'same-width extension has predecessor')
        u=(possible&-possible).bit_length()-1
        paths=reconstruct(smaller,adj,K,E,u);paths[-1].append(v+1)
    else:
        need(int(K[smaller])+1==int(K[mask]),'one-unit reset')
        paths=reconstruct(smaller,adj,K,E);paths.append([v+1])
    return paths

def components(mask,adj):
    number=0
    while mask:
        frontier=mask&-mask;mask^=frontier;number+=1
        while frontier:
            b=frontier&-frontier;frontier-=b;v=b.bit_length()-1
            nxt=int(adj[v])&mask;mask^=nxt;frontier|=nxt
    return number

def run():
    n=25;adj=np.array([sum(1<<j for j in range(n) if i!=j and isqrt(i+j+2)**2==i+j+2) for i in range(n)],dtype=np.uint32)
    start=time.time();K,E=cover_dp(adj);print('DP seconds',round(time.time()-start,2),flush=True)
    out=[]
    for n in range(1,26):
        mask=(1<<n)-1;paths=reconstruct(mask,adj,K,E);need(len(paths)==int(K[mask]),'minimum-width witness')
        need(sorted(v for p in paths for v in p)==list(range(1,n+1)),'each vertex used exactly once')
        need(all(isqrt(a+b)**2==a+b for p in paths for a,b in zip(p,p[1:])),'square-sum path edges')
        cut=None
        for size in range(0,6):
            for S in combinations(range(n),size):
                remaining=mask^sum(1<<v for v in S);c=components(remaining,adj)
                if c-size==len(paths):cut={'vertices':[v+1 for v in S],'components_after_deletion':c};break
            if cut is not None:break
        out.append({'n':n,'minimum_path_cover':len(paths),'optimal_endpoint_mask':int(E[mask]),'paths':paths,'sharp_cut_bound':cut,'new_neighbors':[j+1 for j in range(n-1) if int(adj[n-1])>>j&1]})
    tag='' if __debug__ else '.optimized';result={'scope':'Exact minimum disjoint spanning path-cover width of each induced integer square-sum graph on1..n, n1..25. Full subset DP on2^25 masks; endpoint support, not path multiplicity. Each witness is checked literally; matching vertex cuts give separate optimality certificates when present.','records':out}
    (HERE/('square_covers'+tag+'.json.gz')).write_bytes(gzip.compress((json.dumps(result,separators=(',',':'))+'\n').encode(),mtime=0))
    print([(r['n'],r['minimum_path_cover'],r['sharp_cut_bound']) for r in out],flush=True)
if __name__=='__main__':run()
