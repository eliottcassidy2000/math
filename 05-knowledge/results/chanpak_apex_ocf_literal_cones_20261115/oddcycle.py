import random, itertools
def ham_and_P(out,n,qs):
    full=(1<<n)-1; dp=[[0]*n for _ in range(1<<n)]; H=[0]*(1<<n)
    for v in range(n): dp[1<<v][v]=1
    for S in range(1,1<<n):
        for v in range(n):
            x=dp[S][v]
            if not x: continue
            H[S]+=x; nb=out[v]&~S
            while nb:
                u=(nb&-nb).bit_length()-1; nb&=nb-1; dp[S|1<<u][u]+=x
    H[0]=1; res={1:H[full]}
    F=H[:]
    for q in range(2,max(qs)+1):
        G=[0]*(1<<n)
        for R in range(1<<n):
            s=0; S=R
            while True:
                s+=H[S]*F[R^S]
                if S==0: break
                S=(S-1)&R
            G[R]=s
        F=G; res[q]=F[full]
    return res
def odd_cycles(out,n):
    cyc=[]
    for L in range(3,n+1,2):
        for S in itertools.combinations(range(n),L):
            s0=S[0]
            for perm in itertools.permutations(S[1:]):
                seq=(s0,)+perm
                if all((out[seq[i]]>>seq[(i+1)%L])&1 for i in range(L)): cyc.append(frozenset(S))
    return cyc   # one entry per directed cycle (vertex set repeated if several cyclic orders)
def packing_poly(cyc,n,q):
    # sum over sets of pairwise disjoint directed odd cycles of 2^|K| q^(n - |V(K)| + |K|)
    m=len(cyc); tot=0
    def rec(i,used,k,nv):
        nonlocal tot
        if i==m: tot+=2**k*q**(n-nv+k); return
        rec(i+1,used,k,nv)
        if not (cyc[i]&used): rec(i+1,used|cyc[i],k+1,nv+len(cyc[i]))
    rec(0,frozenset(),0,0); return tot
rng=random.Random(1); bad=0
for t in range(80):
    n=rng.randint(3,7); out=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if rng.random()<0.5: out[i]|=1<<j
            else: out[j]|=1<<i
    P=ham_and_P(out,n,(1,2,3)); cyc=odd_cycles(out,n)
    for q in (1,2,3):
        if packing_poly(cyc,n,q)!=P[q]: bad+=1
print('odd-cycle packing form P_T(q) = sum_K 2^|K| q^(n-|V(K)|+|K|) on 80 random tournaments (n<=7), q=1,2,3: %d mismatches'%bad)
# apex over a transitive tape: T = TT_m (i -> j for i < j) plus apex v beating exactly X; claim H = 1 + p(u), u = sum_{i not in X} 2^i
def p(u):
    r=0
    while u>1:
        w=u>>1
        if not (u&1): r+=w+(w&1)
        u=w
    return r
bad=0; tested=0
for m in range(1,10):
    for X in range(1<<m):
        n=m+1; out=[0]*n
        for i in range(m):
            for j in range(i+1,m): out[i]|=1<<j
        for i in range(m):
            if (X>>i)&1: out[m]|=1<<i
            else: out[i]|=1<<m
        h=ham_and_P(out,n,(1,))[1]; u=((1<<m)-1)^X
        tested+=1
        if h!=1+p(u): bad+=1
print('apex over a transitive tape: H = 1 + p(u), u = complement of the out-set, all %d cases m <= 9: %d mismatches'%(tested,bad))
