import sys
from itertools import combinations
def pairs(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def run(n):
    P=pairs(n); m=len(P); idx={p:k for k,p in enumerate(P)}
    swm=set()
    for S in range(1<<n):
        msk=0
        for k,(i,j) in enumerate(P):
            if ((S>>i)&1)!=((S>>j)&1): msk|=1<<k
        swm.add(msk)
    swm=sorted(swm)
    def canon(t): return min(t^s for s in swm)
    # generators: adjacent transpositions
    gens=[]
    for a in range(n-1):
        pm=list(range(n)); pm[a],pm[a+1]=pm[a+1],pm[a]
        tab=[]
        for k,(i,j) in enumerate(P):
            x,y=pm[i],pm[j]
            if x<y: tab.append((idx[(x,y)],0))
            else:   tab.append((idx[(y,x)],1))
        gens.append(tab)
    def app(t,tab):
        v=0
        for k,(kk,f) in enumerate(tab):
            v |= (((t>>k)&1)^f)<<kk
        return v
    # enumerate class reps
    reps=set()
    seen=bytearray(1<<m)
    for t in range(1<<m):
        if seen[t]: continue
        c=canon(t)
        for s in swm: seen[t^s]=1
        reps.add(c)
    reps=sorted(reps)
    ridx={r:i for i,r in enumerate(reps)}
    par=list(range(len(reps)))
    def find(x):
        while par[x]!=x:
            par[x]=par[par[x]]; x=par[x]
        return x
    def uni(a,b):
        ra,rb=find(a),find(b)
        if ra!=rb: par[ra]=rb
    for r in reps:
        i=ridx[r]
        for tab in gens:
            uni(i, ridx[canon(app(r,tab))])
    orb=len({find(i) for i in range(len(reps))})
    return len(reps), orb

A=[1,1,1,2,2,6,12,79]
print("n | labelled switching classes | 2^C(n-1,2) | iso classes | A049313(n)")
for n in range(1,9):
    if n>=9: break
    lab,orb=run(n)
    exp=1<<((n-1)*(n-2)//2)
    ok = "OK" if orb==A[n-1] else "MISMATCH"
    print(f"{n} | {lab} | {exp} | {orb} | {A[n-1]}  {ok}")
    sys.stdout.flush()

# ---- appended: invariance tests (H, minFAS, cyclic triangles, triple cocycle) ----
import itertools,sys
from itertools import combinations
def pairs(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def swmasks(n,P):
    s=set()
    for S in range(1<<n):
        m=0
        for k,(i,j) in enumerate(P):
            if ((S>>i)&1)!=((S>>j)&1): m|=1<<k
        s.add(m)
    return sorted(s)
def adj(n,P,t):
    A=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(P):
        if (t>>k)&1: A[i][j]=1
        else: A[j][i]=1
    return A
def ham(n,A):
    full=(1<<n)-1; memo={}
    def go(v,mask):
        if mask==full: return 1
        k=(v,mask)
        if k in memo: return memo[k]
        s=0
        for w in range(n):
            if not (mask>>w)&1 and A[v][w]: s+=go(w,mask|(1<<w))
        memo[k]=s; return s
    return sum(go(v,1<<v) for v in range(n))
def mfas(n,A):
    best=99
    for pm in itertools.permutations(range(n)):
        pos={v:i for i,v in enumerate(pm)}
        b=sum(1 for i in range(n) for j in range(n) if A[i][j] and pos[i]>pos[j])
        if b<best: best=b
        if best==0: break
    return best
def ncyc3(n,A):  # number of cyclic triangles (Kendall-Babington Smith)
    c=0
    for i,j,k in combinations(range(n),3):
        e=A[i][j]+A[j][k]+A[k][i]
        if e==3 or e==0: c+=1
    return c
def triple_inv(n,P,t):
    """cocycle invariant: for i<j<k, product a_ij*a_jk*a_ik (a=+-1). switching-invariant."""
    idx={p:kk for kk,p in enumerate(P)}
    def a(i,j):
        return 1 if (t>>idx[(i,j)])&1 else -1
    return tuple(a(i,j)*a(j,k)*a(i,k) for i,j,k in combinations(range(n),3))

print("=== per-switching-class variation of invariants ===")
print("n | #classes | H varies | minFAS varies | #cyclic-triangles varies | triple-cocycle separates classes")
for n in range(3,7):
    P=pairs(n); m=len(P); N=1<<m; sw=swmasks(n,P)
    seen=bytearray(N); nc=0; vH=0; vF=0; vC=0
    cocyc={}
    exH=None
    for t in range(N):
        if seen[t]: continue
        orb=sorted({t^s for s in sw})
        for x in orb: seen[x]=1
        nc+=1
        As=[adj(n,P,x) for x in orb]
        hs={ham(n,A) for A in As}
        fs={mfas(n,A) for A in As}
        cs={ncyc3(n,A) for A in As}
        if len(hs)>1:
            vH+=1
            if exH is None: exH=sorted(hs)
        if len(fs)>1: vF+=1
        if len(cs)>1: vC+=1
        # cocycle constant on class?
        ivs={triple_inv(n,P,x) for x in orb}
        cocyc[nc]=(len(ivs), triple_inv(n,P,orb[0]))
    sep = (len({v[1] for v in cocyc.values()})==nc) and all(v[0]==1 for v in cocyc.values())
    print(f"{n} | {nc} | {vH}/{nc} | {vF}/{nc} | {vC}/{nc} | {sep}")
    if exH: print(f"      example H-multiset within one class at n={n}: {exH}")
    sys.stdout.flush()

print()
print("=== Redei parity check: is H always ODD, and is H mod 2 switching-invariant? ===")
for n in range(3,7):
    P=pairs(n); m=len(P); N=1<<m; sw=swmasks(n,P)
    allodd=True
    seen=bytearray(N); parvar=0; nc=0
    for t in range(N):
        if seen[t]: continue
        orb=sorted({t^s for s in sw}); 
        for x in orb: seen[x]=1
        nc+=1
        hs=[ham(n,adj(n,P,x)) for x in orb]
        if any(h%2==0 for h in hs): allodd=False
        if len({h%2 for h in hs})>1: parvar+=1
    print(f"n={n}: all H odd (Redei) = {allodd}; classes where H mod 2 varies = {parvar}/{nc}")
    sys.stdout.flush()
