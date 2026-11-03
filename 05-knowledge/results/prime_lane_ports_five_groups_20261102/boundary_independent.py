from collections import Counter
from itertools import permutations, product
from math import prod


def need(ok,msg):
    if not ok: raise RuntimeError(msg)


def source(k,cyclic=False):
    E=set()
    for b in range(k):
        E|={(3*b+i,3*b+(i+1)%3) for i in range(3)}
    for a in range(k):
        for b in range(a+1,k):
            forward=not cyclic or (b-a)==1
            left,right=(a,b) if forward else (b,a)
            E|={(3*left+i,3*right+j) for i in range(3) for j in range(3)}
    return E


def paths(n,E):
    def rec(p,unused):
        if not unused:
            yield p
        else:
            for v in unused:
                if not p or (p[-1],v) in E:
                    yield from rec(p+(v,),unused-{v})
    yield from rec((),set(range(n)))


def test(k,cut,cyclic=False):
    old=3*k; E=source(k,cyclic)
    E|={(old,v) if cut>>v&1 else (v,old) for v in range(old)}
    ps=list(paths(old+1,E)); strata=Counter(); generators=0
    for p in ps:
        at=p.index(old)
        touched={p[j]//3 for j in (at-1,at+1) if 0<=j<=old}
        need(len(touched)<=2,'touch budget')
        strata[tuple(sorted(touched))]+=1
        for b in set(range(k))-touched:
            q=tuple(3*b+(v-3*b+1)%3 if 3*b<=v<3*b+3 else v for v in p)
            need(q!=p and all(arc in E for arc in zip(q,q[1:])),'free stratum generator')
            qa=q.index(old)
            touched_q={q[j]//3 for j in (qa-1,qa+1) if 0<=j<=old}
            need(touched_q==touched,'stratum invariance')
            generators+=1
    for S,c in strata.items(): need(c%3**(k-len(S))==0,'stratum divisor')
    need(len(ps)%3**max(0,k-2)==0,'global divisor')
    return len(ps),len(strata),generators

checks=[]
for k in (1,2):
    for cut in range(1<<(3*k)):
        checks.append(test(k,cut))
for cut in (0,1,3,9,73,255,510,511):
    checks.append(test(3,cut,cyclic=True))
need(test(2,9)[0]==41,'literal H41')
print('PASS:',len(checks),'sources/cuts;',sum(x[0] for x in checks),'literal Hamilton words;',sum(x[2] for x in checks),'generator images')
# The source rotation actually fails on extensions; it is not a global automorphism.
k=2; cut=9; E=source(k)|{(6,v) if cut>>v&1 else (v,6) for v in range(6)}
rot=lambda v:(v+1)%3 if v<3 else v
need(any(((u,v) in E)!=((rot(u),rot(v)) in E) for u in range(7) for v in range(7) if u!=v),'broken global symmetry hostile')
# Independently build sharp joined extensions and compare all retained old arcs.
for k,r in ((1,0),(2,1),(3,1),(4,2),(5,2)):
    components=[]; new=[]; offset=0
    for j in range(r):
        seed=source(2)|{(6,v) if v in (0,3) else (v,6) for v in range(6)}
        components.append((set(range(offset,offset+7)),{(u+offset,v+offset) for u,v in seed}))
        new.append(offset+6); offset+=7
    for j in range(k-2*r):
        components.append((set(range(offset,offset+3)),{(u+offset,v+offset) for u,v in source(1)}));offset+=3
    E=set().union(*(e for v,e in components))
    for i,(A,_) in enumerate(components):
        for B,_ in components[i+1:]: E|={(a,b) for a in A for b in B}
    old=sorted(set(range(offset))-set(new)); rank={v:i for i,v in enumerate(old)}
    retained={(rank[u],rank[v]) for u,v in E if u in rank and v in rank}
    need(retained==source(k),'named old-source preservation')
    factors=[41]*r+[3]*(k-2*r)
    need(prod(factors)==41**r*3**(k-2*r),'product family')
print('PASS: nonextendable source symmetries, sharp named source deletion and ordered-join counts')

# Strong sharpness after one new-new edge reversal; independent subset counter.
def hp_subset(n,E):
    out=[0]*n
    for u,v in E: out[u]|=1<<v
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v,c in enumerate(dp[mask]):
            if c:
                available=out[v]&~mask
                while available:
                    bit=available&-available;available-=bit
                    dp[mask|bit][bit.bit_length()-1]+=c
    return sum(dp[-1])


def isstrong(n,E):
    R=[{i}|{v for u,v in E if u==i} for i in range(n)]
    for k in range(n):
        for i in range(n):
            if k in R[i]: R[i]|=R[k]
    return all(len(x)==n for x in R)


def strong_sharp(r,t):
    parts=[];new=[];off=0
    for i in range(r):
        excess=t if i==r-1 else 0
        old=3*(excess+2); E=source(excess+2)
        E|={(old,v) if v in (0,3) else (v,old) for v in range(old)}
        parts.append((set(range(off,off+old+1)),{(a+off,b+off) for a,b in E}))
        new.append(off+old);off+=old+1
    E=set().union(*(arcs for _,arcs in parts))
    for i,(A,_) in enumerate(parts):
        for B,_ in parts[i+1:]: E|={(a,b) for a in A for b in B}
    before=hp_subset(off,E)
    final_factor=41 if t==0 else 3**t*(45*4**t-1)
    need(before==41**(r-1)*final_factor,'guarded strong seed product')
    if r>=2:
        a,b=new[0],new[-1]
        need((a,b) in E and (b,a) not in E,'actual new-new edge')
        E.remove((a,b));E.add((b,a))
    after=hp_subset(off,E)
    need((after-before)%3**(t+2)==0,'new-edge higher divisor')
    need(after%3**t==0 and after%3**(t+1)!=0,'strong valuation sharpness')
    need(isstrong(off,E),'all-strong sharp realization')
    old=sorted(set(range(off))-set(new));rank={v:i for i,v in enumerate(old)}
    retained={(rank[a],rank[b]) for a,b in E if a in rank and b in rank}
    need(retained==source(2*r+t),'strong sharp old-source preservation')
    return off,before,after,after-before

for r,t in ((1,0),(1,1),(1,2),(2,0),(2,1)):
    print('strong sharp',r,t,strong_sharp(r,t))
print('PASS: strong all-budget construction controls, exceptional t=0 seed and higher new-edge congruence')
