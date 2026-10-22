#!/usr/bin/env python3
"""Independent OCF orbit certificates; no repository imports or removable asserts."""
import itertools, json, random
from collections import Counter
from pathlib import Path

def require(p, message):
    if not p: raise RuntimeError(message)

def hp_table(out):
    n=len(out); dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,x in enumerate(dp[mask]):
            if x:
                nxt=out[v]&~mask
                while nxt:
                    bit=nxt&-nxt; nxt-=bit
                    dp[mask|bit][bit.bit_length()-1]+=x
    return dp

def cycles(out):
    n=len(out); result=[]
    for k in range(3,n+1,2):
        for vertices in itertools.combinations(range(n),k):
            root=vertices[0]
            for tail in itertools.permutations(vertices[1:]):
                c=(root,)+tail
                if all(out[c[i]]>>c[(i+1)%k]&1 for i in range(k)): result.append(c)
    return result

def canonical(c):
    j=c.index(min(c)); return c[j:]+c[:j]

def packs(masks):
    result=[]
    def visit(start,used,chosen):
        result.append(tuple(chosen))
        for i in range(start,len(masks)):
            if not used&masks[i]: visit(i+1,used|masks[i],chosen+[i])
    visit(0,0,[]); return result

def orientation(n,bits):
    remaining={(i,j) for i in range(n) for j in range(i+1,n)}; orbits=[]
    while remaining:
        p=min(remaining); q=(n-1-p[1],n-1-p[0]); o={p,q};remaining-=o;orbits.append(o)
    out=[0]*n
    for k,o in enumerate(orbits):
        for i,j in o:
            a,b=(i,j) if bits>>k&1 else (j,i);out[a]|=1<<b
    return out,len(orbits)

def inspect(out):
    n=len(out); full=(1<<n)-1; dp=hp_table(out);h=sum(dp[-1]);cs=cycles(out);index={c:i for i,c in enumerate(cs)}
    masks=[sum(1<<v for v in c) for c in cs]
    tau=[index[canonical(tuple(n-1-v for v in reversed(c)))] for c in cs]
    require(all(tau[tau[i]]==i for i in range(len(cs))),'cycle involution')
    ps=packs(masks); alpha=Counter(map(len,ps)); fixed=Counter()
    for p in ps:
        if tuple(sorted(tau[i] for i in p))==p:fixed[len(p)]+=1
    require(h==sum(v*2**k for k,v in alpha.items()),'OCF versus independent HP DP')
    fixed_cycles=[i for i in range(len(cs)) if tau[i]==i]
    pair_orbits=[(i,tau[i]) for i in range(len(cs)) if i<tau[i]]
    admissible=[(i,j) for i,j in pair_orbits if not masks[i]&masks[j]]
    qmasks=[masks[i]|masks[j] for i,j in admissible]
    expected=Counter(2*len(p) for p in packs(qmasks))
    for f in fixed_cycles:
        available=[m for m in qmasks if not m&masks[f]]
        expected.update(1+2*len(p) for p in packs(available))
    require(fixed==expected,'fixed packing quotient polynomial')
    require(all((alpha[k]-fixed[k])%2==0 for k in alpha),'moving packing orbits')
    overlap=sum(bool(masks[i]&masks[j]) for i,j in pair_orbits)
    require((h-1-2*len(fixed_cycles)-4*overlap)%8==0,'mod8 fixed-cycle/overlap law')
    if n%2==0:
        require(not fixed_cycles and h%4==1,'even converse mod4')
        require(all(alpha[k]%2==0 for k in alpha if k%2),'odd packing coefficients even')
    rev=[sum(1<<u for u in range(n) if out[u]>>v&1) for v in range(n)]
    starts=hp_table(rev)[-1]
    for a in range(n):
        R=dp[-1][a]+sum(dp[m][a]*sum(dp[full^m]) for m in range(1,full) if m>>a&1)
        require((R-starts[a])%2==0,'universal port parity')
    return {'H':h,'fixed_odd_cycles':len(fixed_cycles),'overlapping_pair_orbits':overlap,'odd_cycles':len(cs),'alpha':dict(alpha)}

def main():
    rng=random.Random(2908); results={}; examples={};total=0
    for n in range(1,9):
        _,d=orientation(n,0)
        codes=list(range(1<<d)) if n<=6 else sorted(rng.sample(range(1<<d),32))
        hist=Counter()
        for code in codes:
            data=inspect(orientation(n,code)[0]);hist[data['H']%8]+=1
            examples.setdefault(f'{n}:{data["H"]%8}',{'code':code,**data});total+=1
        results[n]={'tested':len(codes),'universe':1<<d,'coverage':'complete fixed reversal-involution universe' if n<=6 else 'deterministic sample','H_mod8':dict(hist)}
    # Hypothesis is essential: C3 followed by a sink is not self-converse and has H=3.
    out=[(1<<1)|(1<<3),(1<<2)|(1<<3),(1<<0)|(1<<3),0]
    require(sum(hp_table(out)[-1])==3,'hostile non-self-converse example')
    report={'status':'PASS','total':total,'orders':results,'examples':examples,'hostile_n4_H3_outmasks':out}
    print(json.dumps(report,sort_keys=True,indent=2))
if __name__=='__main__':main()
