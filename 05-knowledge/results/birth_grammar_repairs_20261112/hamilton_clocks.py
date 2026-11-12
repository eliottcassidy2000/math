"""Exact one-port Hamilton skeletons, without inflating literal graphs."""
from common import *
from numba import njit
from collections import Counter
import csv,time

@njit(cache=True)
def skeleton(rows,p):
    n=len(rows);m=n-1;verts=np.empty(m,dtype=np.int64);z=0
    for v in range(n):
        if v!=p:verts[z]=v;z+=1
    degree=0
    for j in range(n):degree+=(rows[p]>>j)&1
    R=1+min(n-1-degree,degree)
    N=1<<m;dp=np.zeros((N,m+1,R+1),dtype=np.int64)
    dp[0,m,1]=1
    for j in range(m):dp[1<<j,j,0]=1
    for mask in range(N):
        for r in range(R+1):
            for j in range(m+1):
                c=dp[mask,j,r]
                if not c:continue
                v=p if j==m else verts[j]
                if j!=m and r<R and rows[v]>>p&1:dp[mask,m,r+1]+=c
                for k in range(m):
                    if not mask>>k&1 and rows[v]>>verts[k]&1:dp[mask|(1<<k),k,r]+=c
    out=np.zeros(R,dtype=np.int64)
    for r in range(1,R+1):
        for j in range(m+1):out[r-1]+=dp[N-1,j,r]
    return out

def graphs(max_n):
    with gzip.open(HERE.parent/'metagraph_aspect_sequences_20261111/records.jsonl.gz','rt') as f:
        for line in f:
            r=json.loads(line)
            if len(r['rows'])<=max_n:yield r['rows']
    if max_n>=9:
        from two_block_probe import decode
        with gzip.open(HERE.parent/'blue_quotient_20261020b/n9_classes.tsv.gz','rt') as f:
            for r in csv.DictReader(f,delimiter='\t'):yield decode(int(r['code']),9)

def trim(e):
    while e and e[-1]==0:e.pop()
    return tuple(e)

def run(max_n=8):
    found={};stats=Counter();start=time.time();controls=0
    for idx,rows in enumerate(graphs(max_n)):
        n=len(rows);stats[n]+=1
        for p in range(n):
            u=[int(x) for x in skeleton(np.array(rows,dtype=np.int64),p)];e=trim(expcoeff(u))
            if idx<80 or idx%10000==0:
                h,_,_=ib.coefficients(rows,p);need(u==[int(x) for x in h[1:len(u)+1]] and not any(h[len(u)+1:]),'independent full H skeleton')
                controls+=1
            if e not in found:found[e]={'rows':rows,'port':p,'Stirling_coefficients':u,'exponential_coefficients':list(e)}
        if idx%10000==0:print(json.dumps({'sources':idx,'profiles':len(found),'seconds':round(time.time()-start,1)}),flush=True)
    out={'scope':f'All ordinary tournaments through{max_n}, every pointed vertex. Exact complete Hamilton skeleton by repeated-port subset DP; distinct exponential sequences, not marked graph classes.','source_counts':dict(stats),'distinct_profiles':len(found),'independent_skeleton_controls':controls,'records':list(found.values())}
    save(f'hamilton_clocks_{max_n}.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k!='records'}),flush=True)
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 8)
