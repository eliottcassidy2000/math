import random, itertools
def ham_and_p2(adj, n):
    # H = Hamiltonian paths; P2 = ordered pairs (Q1,Q2) of disjoint nonempty paths covering all vertices
    full=(1<<n)-1
    cnt=[0]*(1<<n)   # number of paths with vertex set S (any start/end)
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for S in range(1,1<<n):
        for v in range(n):
            c=dp[S][v]
            if not c: continue
            cnt[S]+=c
            nb=adj[v]&~S
            while nb:
                u=(nb&-nb).bit_length()-1; nb&=nb-1; dp[S|1<<u][u]+=c
    H=cnt[full]
    P2=sum(cnt[S]*cnt[full^S] for S in range(1,full))
    return H,P2
def rand_t(n,rng):
    adj=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if rng.random()<0.5: adj[i]|=1<<j
            else: adj[j]|=1<<i
    return adj
def cyclic_join(A,a,B,b):
    n=a+b+1; adj=[0]*n
    for i in range(a): adj[i]=A[i]
    for i in range(b): adj[a+i]=B[i]<<a
    c=a+b
    for i in range(a):
        for j in range(b): adj[i]|=1<<(a+j)      # A => B
    for j in range(b): adj[a+j]|=1<<c            # B => c
    for i in range(a): adj[c]|=1<<i              # c => A
    return adj,n
rng=random.Random(11); ok=0; tot=0
for t in range(200):
    a=rng.randint(1,5); b=rng.randint(1,5)
    A=rand_t(a,rng); B=rand_t(b,rng)
    hA,pA=ham_and_p2(A,a); hB,pB=ham_and_p2(B,b)
    J,n=cyclic_join(A,a,B,b); N,_=ham_and_p2(J,n)
    tot+=1; ok+= (N==(hA+pA)*(hB+pB)+2*hA*hB)
print('cyclic join law N = (hA+PA)(hB+PB) + 2 hA hB:',ok,'of',tot)
# diagonal: N = (h+P)^2 + 2h^2 ; list primes produced by diagonal joins of all tournaments up to order 4..5 (random)
from sympy import isprime
pr=set()
for t in range(3000):
    a=rng.randint(1,6); A=rand_t(a,rng); h,P=ham_and_p2(A,a)
    N=(h+P)**2+2*h*h
    if isprime(N): pr.add(N)
print('diagonal-join primes found:',sorted(pr)[:15],' residues mod 8:',sorted(set(p%8 for p in pr)))
