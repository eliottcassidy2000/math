"""Full single-module visit profiles; no transitive two-run restriction."""
from pathlib import Path
from collections import Counter
import gzip,json
import numpy as np
from numba import njit
from frame_probe import ordinary_bank,hp,clone,strong,need

HERE=Path(__file__).resolve().parent

@njit
def visits(T,r):
    n=len(T);ext=np.array([v for v in range(n) if v!=r],dtype=np.int64);m=n-1
    bound=(n-1)//2+1;size=1<<m
    D=np.zeros((size,m,bound+1),dtype=np.int64)
    R=np.zeros((size,bound+1),dtype=np.int64);R[0,1]=1
    for v in range(m):D[1<<v,v,0]=1
    for S in range(size):
        for v in range(m):
            if S>>v&1 and T[ext[v]]>>r&1:
                for j in range(1,bound+1):R[S,j]+=D[S,v,j-1]
        for w in range(m):
            if S>>w&1:continue
            dest=S|(1<<w)
            if T[r]>>ext[w]&1:
                for j in range(bound+1):D[dest,w,j]+=R[S,j]
            for v in range(m):
                if S>>v&1 and T[ext[v]]>>ext[w]&1:
                    for j in range(bound+1):D[dest,w,j]+=D[S,v,j]
    result=np.zeros(bound,dtype=np.int64)
    for j in range(1,bound+1):
        result[j-1]=R[-1,j]+D[-1,:,j].sum()
    return result

def covers_tt(m,k):
    from math import comb
    return sum((-1)**(k-i)*comb(k,i)*i**m for i in range(k+1))

def padded(v):return tuple(v)+(0,)*(4-len(v))

def run():
    catalog={};occurrences=0;controls=0
    for (n,code),T in sorted(ordinary_bank().items()):
        if n<3 or not strong(T):continue
        for r in range(n):
            v=padded(tuple(map(int,visits(np.array(T,dtype=np.int64),r))))
            occurrences+=1
            if v not in catalog:
                vals=[]
                for m in range(1,5):
                    actual=hp(clone(T,r,m));pred=sum(a*covers_tt(m,k) for k,a in enumerate(v,1))
                    need(actual==pred,'walk DP versus literal module expansion')
                    vals.append(actual);controls+=1
                catalog[v]=dict(n=n,rows=T,root=r,profile=v,checks=vals)
    transitions=[]
    for v,source in sorted(catalog.items()):
        for digit in (2,4):
            w=(3*v[0]+digit,)+tuple(3*x for x in v[1:])
            if w in catalog:transitions.append(dict(digit=digit,source=source,target=catalog[w],max_runs=max(i for i,x in enumerate(v,1) if x)))
    result=dict(status='FINITE-EXACT catalog; each profile is a universal arbitrary-module identity',scope='All strong ordinary frames3..8, every module vertex, full possible run count up to4. One representative per entire response vector. Profile comparison supplies actual core-replacement identities, not arbitrary target coverage.',occurrences=occurrences,profiles=len(catalog),controls=controls,records=list(catalog.values()),transitions=transitions)
    suffix='' if __debug__ else '.optimized'
    (HERE/('visits'+suffix+'.json.gz')).write_bytes(gzip.compress(json.dumps(result,separators=(',',':')).encode(),mtime=0))
    print('PROFILES',len(catalog),'TRANSITIONS',Counter((t['digit'],t['max_runs']) for t in transitions),'CONTROLS',controls,flush=True)
    for h in (3,15,17,25,27,35,41,43,53,55,57,63,67):
        print('HARD_PARENT',h,[(t['digit'],t['source']['profile'],t['source']['rows'],t['source']['root'],t['target']['rows'],t['target']['root']) for t in transitions if t['source']['profile'][0]==h][:4],flush=True)

if __name__=='__main__':run()
