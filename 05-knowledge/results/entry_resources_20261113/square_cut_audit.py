"""Independent all-subset component lower bounds, especially the order24 gap."""
from pathlib import Path
import numpy as np
from numba import njit
from math import isqrt
import json,gzip,time,sys
HERE=Path(__file__).resolve().parent

@njit(cache=True)
def cuts(adj):
    n=len(adj);N=1<<n;cc=np.zeros(N,dtype=np.uint8);pop=np.zeros(N,dtype=np.uint8)
    idx=np.array([0,1,28,2,29,14,24,3,30,22,20,15,25,17,4,8,31,27,13,23,21,19,16,7,26,12,18,6,11,5,10,9])
    best=np.zeros(n+1,dtype=np.int64);witness=np.zeros(n+1,dtype=np.int64);maxval=0;maxmask=0;nextn=1
    for mask in range(1,N):
        pop[mask]=pop[mask&(mask-1)]+1
        frontier=mask&-mask;remain=mask^frontier
        while frontier:
            bit=frontier&-frontier;frontier-=bit;v=idx[((bit*0x077CB531)&0xffffffff)>>27]
            new=int(adj[v])&remain;remain^=new;frontier|=new
        cc[mask]=cc[remain]+1
        score=int(cc[mask])+int(pop[mask])
        # Prefer the smallest deleted set among equal maxima.
        if score>maxval or (score==maxval and pop[mask]>pop[maxmask]):maxval=score;maxmask=mask
        if mask==(1<<nextn)-1:
            best[nextn]=maxval-nextn;witness[nextn]=maxmask;nextn+=1
    return best,witness,cc,pop

def run():
    n=24;adj=np.array([sum(1<<j for j in range(n) if i!=j and isqrt(i+j+2)**2==i+j+2) for i in range(n)],dtype=np.uint32)
    start=time.time();best,win,cc,pop=cuts(adj)
    rows=[]
    for n in range(1,25):
        kept=int(win[n]);deleted=[i+1 for i in range(n) if not kept>>i&1]
        rows.append({'n':n,'maximum_components_minus_deleted':int(best[n]),'deleted':deleted,'components':int(cc[kept])})
    out={'scope':'Exhaustive component bound over every vertex deletion set for each Q_n, n1..24; one all-subset recursion. This is a lower bound for path-cover width, not assumed sufficient.','subsets':1<<24,'records':rows}
    tag='' if __debug__ else '.optimized';(HERE/('square_cut_audit'+tag+'.json.gz')).write_bytes(gzip.compress((json.dumps(out,separators=(',',':'))+'\n').encode(),mtime=0))
    print(rows,'seconds',round(time.time()-start,2),flush=True)
if __name__=='__main__':run()
