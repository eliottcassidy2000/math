#!/usr/bin/env python3
"""Independent subset-path / Newton-profile audit. No producer imports or inputs."""
import itertools, math, json, hashlib
from functools import lru_cache

def require(x,s):
    if not x: raise ValueError(s)
def graphs(n):
    edges=list(itertools.combinations(range(n),2))
    for bits in range(1<<len(edges)):
        rows=[0]*n
        for k,(i,j) in enumerate(edges):
            a,b=(i,j) if (bits>>k)&1 else (j,i);rows[a]|=1<<b
        yield tuple(rows)
def transitive(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def substitute(template,parts):
    offsets=[];n=0
    for p in parts:offsets.append(n);n+=len(p)
    rows=[]
    for i,p in enumerate(parts):
        outside=sum(((1<<len(q))-1)<<offsets[j] for j,q in enumerate(parts) if template[i]>>j&1)
        rows.extend((r<<offsets[i])|outside for r in p)
    return tuple(rows)
def path_counts(rows):
    n=len(rows);d=[[0]*n for _ in range(1<<n)]
    for v in range(n):d[1<<v][v]=1
    for m in range(1,1<<n):
        for v in range(n):
            if m>>v&1:
                prev=m^(1<<v)
                if prev:d[m][v]=sum(d[prev][u] for u in range(n) if prev>>u&1 and rows[u]>>v&1)
    return (1,)+tuple(sum(row) for row in d[1:])
def profile(rows):
    h=path_counts(rows)
    @lru_cache(None)
    def partitions(m):
        if m==0:return (1,)
        least=m&-m;out=[0]*(m.bit_count()+1);s=m
        while s:
            if s&least:
                for k,c in enumerate(partitions(m^s)):out[k+1]+=h[s]*c
            s=(s-1)&m
        return tuple(out)
    return partitions(len(h)-1)
def q(pc,t):
    f=1;out=pc[0]
    for k in range(1,len(pc)):
        f*=t-k+1;out+=pc[k]*f
    return out
def from_values(vals):
    out=[];k=0
    while vals:
        require(vals[0]%math.factorial(k)==0,'Newton integral')
        out.append(vals[0]//math.factorial(k));vals=[y-x for x,y in zip(vals,vals[1:])];k+=1
    return tuple(out)
@lru_cache(None)
def ordered(a,b):return from_values([q(a,j)*q(b,j) for j in range(len(a)+len(b)-1)])
@lru_cache(None)
def collar_pair(a,b):
    return from_values([j*(q(a,j)*q(b,j)+(q(a,j+1)-q(a,j))*(q(b,j+1)-q(b,j))+(q(a,j)-q(a,j-1))*(q(b,j)-q(b,j-1))) for j in range(len(a)+len(b))])
def c3_h(a,b,c):
    fs=[tuple(math.factorial(k)*v for k,v in enumerate(p)) for p in (a,b,c)]
    total=0
    for i,j,k in itertools.product(*(range(1,len(f)) for f in fs)):
        spread=max(i,j,k)-min(i,j,k)
        if spread<=1:total+=(3 if spread==0 else 1)*fs[0][i]*fs[1][j]*fs[2][k]
    return total

def compile_old(seed,limit=4095):
    N=1+(limit-1).bit_length();bank=[{} for _ in range(N+1)]
    def put(n,p,d):
        if p[1]<=limit:bank[n][p]=min(d,bank[n].get(p,10**9))
    for p in seed:put(len(p)-1,p,0)
    for n in range(2,N+1):
        for i in range(1,n):
            for a,da in list(bank[i].items()):
                for b,db in list(bank[n-i].items()):
                    if a[1]*b[1]<=limit:put(n,ordered(a,b),max(da,db))
        for i in range(1,n-1):
            for a,da in bank[i].items():
                for b,db in bank[n-i-1].items():
                    h=2*a[1]*b[1]+(q(a,2)-a[1])*(q(b,2)-b[1])
                    if h<=limit:put(n,collar_pair(a,b),1+max(da,db))
    generators={p[1]:0 for p in seed}
    # Every bounded-order expression is valid as an exterior factor; seeds too.
    for row in bank:
        for p,d in row.items():generators[p[1]]=min(d,generators.get(p[1],10**9))
    depth={1:0}
    for h in range(2,limit+1):
        candidates=[generators[h]] if h in generators else []
        for g,d in generators.items():
            if 1<g<h and h%g==0 and h//g in depth:candidates.append(max(d,depth[h//g]))
        if candidates:depth[h]=min(candidates)
    return bank,depth


#!/usr/bin/env python3
"""Independent full C3 compiler via direct quotient-walk convolution."""
import sys
a=sys.modules[__name__]
import itertools,math,json,hashlib
from functools import lru_cache

def kernel(N):
    # Integer coefficient of W^d; no Newton transform used.
    W={}
    for length in range(1,N+1):
        for start in range(3):
            c=[0]*3
            for j in range(length):c[(start+j)%3]+=1
            c=tuple(c);W[c]=W.get(c,0)+1
    rows={};power={(0,0,0):1}
    for d in range(1,N+1):
        nxt={}
        for c,x in power.items():
            for w,y in W.items():
                v=tuple(c[i]+w[i] for i in range(3))
                if sum(v)<=N:nxt[v]=nxt.get(v,0)+x*y
        power=nxt
        for c,x in power.items():
            if min(c)==0:continue
            z=x*math.prod(math.factorial(i) for i in c)
            a.require(z%math.factorial(d)==0,'decorated walk integrality')
            rows.setdefault(c,[0]*(N+1))[d]=z//math.factorial(d)
    return rows
K=kernel(13)
@lru_cache(None)
def full(p,q,r):
    n=len(p)+len(q)+len(r)-3;out=[0]*(n+1)
    for i,j,k in itertools.product(range(1,len(p)),range(1,len(q)),range(1,len(r))):
        factor=p[i]*q[j]*r[k]
        for d,c in enumerate(K[(i,j,k)][:n+1]):out[d]+=factor*c
    return tuple(out)

def compile(seed,X=4095):
    N=1+(X-1).bit_length();bank=[{} for _ in range(N+1)]
    def put(p,d):
        if p[1]<=X:bank[len(p)-1][p]=min(d,bank[len(p)-1].get(p,999))
    for p in seed:put(p,0)
    for n in range(2,N+1):
        for i in range(1,n):
            for p,dp in bank[i].items():
                for q,dq in bank[n-i].items():
                    if p[1]*q[1]<=X:put(a.ordered(p,q),max(dp,dq))
        for i in range(1,n//3+1):
            for j in range(i,(n-i)//2+1):
                k=n-i-j
                for p,dp in bank[i].items():
                    for q,dq in bank[j].items():
                        for r,dr in bank[k].items():
                            if a.c3_h(p,q,r)<=X:put(full(p,q,r),1+max(dp,dq,dr))
    generators={}
    for row in bank:
        for p,d in row.items():generators[p[1]]=min(d,generators.get(p[1],999))
    reachable={1:0}
    for n in range(2,X+1):
        depths=[generators[n]] if n in generators else []
        for g,d in generators.items():
            if 1<g<n and n%g==0 and n//g in reachable:depths.append(max(d,reachable[n//g]))
        if depths:reachable[n]=min(depths)
    return bank,reachable

def main():
    seeds={}
    for n in range(1,6):
        for r in a.graphs(n):seeds.setdefault(a.profile(r),r)
    small=[(p,r) for p,r in seeds.items() if len(r)<=3];controls=0
    for (p,pr),(q,qr),(r,rr) in itertools.product(small,repeat=3):
        a.require(full(p,q,r)==a.profile(a.substitute((2,4,1),(pr,qr,rr))),'literal full-profile control');controls+=1
    extra=a.profile((40,1,3,6,15,30));a.require(extra[1]==19,'H19 seed')
    result={'status':'PASS','method':'direct powers of quotient-walk series, no Newton transform','full_profile_literal_controls':controls,'H19_profile':extra}
    for name,ss in [('seed5',seeds),('seed5_plus_H19',list(seeds)+[extra])]:
        bank,reach=compile(ss)
        result[name]={'scalar_count':len(reach),'profile_counts':[len(b) for b in bank[1:]],'depth_counts':{str(d):sum(v==d for v in reach.values()) for d in sorted(set(reach.values()))},'scalar_depth_sha256':hashlib.sha256(json.dumps(sorted(reach.items()),separators=(',',':')).encode()).hexdigest(),'first_missing':sorted(set(range(1,4096,2))-reach.keys())[:40]}
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
