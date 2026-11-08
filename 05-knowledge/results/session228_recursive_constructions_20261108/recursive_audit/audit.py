#!/usr/bin/env python3
"""Independent value/Newton-basis compilation and canonical set-partition covers."""
import json, hashlib, itertools, math
from pathlib import Path
from functools import lru_cache
ROOT=Path(__file__).resolve().parent/'inputs'
X=4095; N=13

def check(b,m):
    if not b: raise RuntimeError(m)

def graphs(n):
    pairs=list(itertools.combinations(range(n),2))
    for z in range(1<<len(pairs)):
        rows=[0]*n
        for k,(a,b) in enumerate(pairs):
            if z&(1<<k): a,b=b,a
            rows[a]|=1<<b
        yield tuple(rows)

def counts(rows):
    n=len(rows)
    check(all(not(r>>i&1) and r<(1<<n) for i,r in enumerate(rows)), 'row validity')
    check(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(n) for j in range(i)), 'tournament validity')
    @lru_cache(None)
    def ending(mask,v):
        rest=mask^(1<<v)
        if not rest:return 1
        return sum(ending(rest,u) for u in range(n) if rest>>u&1 and rows[u]>>v&1)
    return [1]+[sum(ending(s,v) for v in range(n) if s>>v&1) for s in range(1,1<<n)]

def literal_values(rows,k):
    # Unordered path covers: the part containing the least vertex is unique.
    h=counts(rows)
    @lru_cache(None)
    def covers(mask):
        if not mask:return (1,)
        low=mask&-mask; other=mask^low; s=other; ans=[0]*(mask.bit_count()+1)
        while True:
            part=s|low
            for d,c in enumerate(covers(mask^part)):ans[d+1]+=h[part]*c
            if s==0:break
            s=(s-1)&other
        return tuple(ans)
    pc=covers(len(h)-1)
    return tuple(sum(pc[d]*math.factorial(t)//math.factorial(t-d) for d in range(min(t,len(rows))+1)) for t in range(k+1))

def differences(values):
    a=list(values); out=[]
    while a:
        out.append(a[0]); a=[y-x for x,y in zip(a,a[1:])]
    return tuple(out)

@lru_cache(None)
def ev(p,t):
    ans=0; c=1
    for i,v in enumerate(p):
        ans+=c*v; c=c*(t-i)//(i+1)
    return ans

def operation(p,q,n,cyclic):
    if cyclic:
        values=[t*(ev(p,t+1)*ev(q,t+1)-ev(p,t+1)*ev(q,t)-ev(p,t)*ev(q,t+1)+2*ev(p,t)*ev(q,t)+(ev(p,t)-ev(p,t-1))*(ev(q,t)-ev(q,t-1))) for t in range(n+1)]
    else: values=[ev(p,t)*ev(q,t) for t in range(n+1)]
    return differences(values)

def physical(a,b,cyclic):
    n=len(a);m=len(b)
    out=tuple(r+(((1<<m)-1)<<n) for r in a)+tuple((r<<n)+((1<<(n+m)) if cyclic else 0) for r in b)
    return out+(((1<<n)-1,) if cyclic else ())

def compile_bank(seedbank,limit,cutoff):
    bank=[{} for _ in range(limit+1)]
    def insert(n,p,d):
        if ev(p,1)<=cutoff: bank[n][p]=min(d,bank[n].get(p,999))
    for n,p in seedbank:
        if n<=limit:insert(n,p,0)
    for n in range(2,limit+1):
        # Intentionally use BOTH operand orders, unlike producer symmetry reduction.
        for cyclic in (False,True):
            for i in range(1,n-int(cyclic)):
                j=n-int(cyclic)-i
                for p,dp in bank[i].items():
                    for q,dq in bank[j].items():
                        h=(2*ev(p,1)*ev(q,1)+(ev(p,2)-ev(p,1))*(ev(q,2)-ev(q,1))) if cyclic else ev(p,1)*ev(q,1)
                        if h<=cutoff:insert(n,operation(p,q,n,cyclic),int(cyclic)+max(dp,dq))
    return bank

def scalar_closure(bank,seeds,cutoff):
    # All bounded profiles may be generators: joins here are redundant, and
    # independently test outer-join factoring against producer cyclic-only bank.
    d={1:0}
    for n,p in seeds:
        h=ev(p,1)
        if h<=cutoff:d[h]=0
    for level in bank:
        for p,depth in level.items():
            h=ev(p,1); d[h]=min(depth,d.get(h,999))
    # Forward multiplicative relaxation, not producer divisor recursion.
    generators=sorted(d.items())
    for v in range(1,cutoff+1):
        if v not in d:continue
        for g,dg in generators:
            if v*g>cutoff:break
            if g>1:d[v*g]=min(d.get(v*g,999),max(d[v],dg))
    return d

def main():
    receipt=json.loads((ROOT/'CHECK.json').read_text())
    seeds=set(); small=[]; total=0
    for n in range(1,6):
        for r in graphs(n):
            v=literal_values(r,n);p=differences(v);seeds.add((n,p));total+=1
            # Separate permutation Hamilton control for every seed.
            hp=sum(all(r[a]>>b&1 for a,b in zip(order,order[1:])) for order in itertools.permutations(range(n)))
            check(hp==v[1],'seed Hamilton permutation')
            if n<=3:small.append((r,p))
    check(len(seeds)==15 and total==1099,'seed counts')
    literal=0
    for a,p in small:
        for b,q in small:
            for cyclic in (False,True):
                r=physical(a,b,cyclic);v=literal_values(r,len(r));out=operation(p,q,len(r),cyclic)
                check(differences(v)==out,'literal operation')
                literal+=1
    bank=compile_bank(seeds,N,X);depth=scalar_closure(bank,seeds,X)
    digest=hashlib.sha256(json.dumps(sorted(depth.items()),separators=(',',':')).encode()).hexdigest()
    check(digest==receipt['scalar_depth_sha256'],'scalar depth digest')
    check([len(x) for x in bank[1:]]==receipt['profile_counts_by_order'],'all profile counts')
    histogram={str(d):sum(x==d for x in depth.values()) for d in sorted(set(depth.values()))}
    check(histogram==receipt['scalar_depth_counts'],'depth histogram')
    witnesses=[]
    for w in receipt['hostile_witnesses']:
        r=tuple(w['rows']);v=literal_values(r,9)
        check(list(v[1:4])==w['Q1_Q2_Q3'],'collision literal values')
        check(all(sum(c*t**i for i,c in enumerate(w['polynomial']))==v[t] for t in range(10)), 'collision entire polynomials')
        rr=physical(physical(r,(0,),True),(0,),True)
        h=counts(rr)[-1];check(h==w['two_collar_H'],'two collar literal')
        witnesses.append(h)
    r=tuple(sum(1<<j for j in range(i+1,4)) for i in range(4));chain=[]
    for k in range(4):
        chain.append(counts(r)[-1]);r=physical(r,(0,),True)
    check(chain==[1,17,241,4001],'collar chain literal')
    for w in receipt['recursive_witnesses']:check(counts(tuple(w['rows']))[-1]==w['H'],'recursive witness literal')
    # Arbitrarily padded outer joins are mathematically handled by closure.
    # Literal checks through 14 vertices include order exceeding N(X).
    padding=[]
    r=(2,4,1)
    for k in range(12):
        h=counts(r)[-1];check(h==3,'exterior padding');padding.append((len(r),h));r=physical(r,(0,),False)
    # General-bank hostile: seed larger than target N must survive as scalar.
    extra=(14,differences([t**14+2*t**12 for t in range(15)]))
    smallbank=compile_bank({(1,(0,1)),extra},3,3)
    check(scalar_closure(smallbank,{(1,(0,1)),extra},3)=={1:0,3:0},'large seed minimum depth')
    unpruned=compile_bank(seeds,7,10**20)
    pruned=compile_bank(seeds,7,21)
    check(all({p:d for p,d in u.items() if ev(p,1)<=21}==v for u,v in zip(unpruned,pruned)), 'target pruning preserves depths')
    files=['PROOF.md','README.md','verify.py','CHECK.json','CHECK.optimized.json']
    print(json.dumps({'status':'PASS','independent_method':'canonical set partitions + endpoint recursion; Newton finite differences; both operand orders; forward multiplicative closure','input_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},'labelled_seeds':total,'seed_profiles':len(seeds),'literal_binary_join_and_cycle_checks':literal,'profiles_by_order':[len(x) for x in bank[1:]],'scalar_count':len(depth),'depth_histogram':histogram,'scalar_depth_sha256':digest,'smallest_by_depth':{str(d):min(v for v,z in depth.items() if z==d) for d in sorted(set(depth.values()))},'collar_chain':chain,'collision_two_collar_counts':witnesses,'outer_padding':padding,'large_seed_depth_control':'PASS'},sort_keys=True,indent=2))
if __name__=='__main__':main()
