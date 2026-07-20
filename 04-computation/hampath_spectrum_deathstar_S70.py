import random
def ham(n, arc):
    N=1<<n; dp=[0]*(N*n)
    for v in range(n): dp[(1<<v)*n+v]=1
    for S in range(N):
        b=S*n
        for v in range(n):
            c=dp[b+v]
            if c:
                for w in arc[v]:
                    if not (S>>w)&1: dp[(S|(1<<w))*n+w]+=c
    f=(N-1)*n; return sum(dp[f:f+n])
def tour(n,bits):
    out=[[] for _ in range(n)]; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if (bits>>idx)&1: out[i].append(j)
            else: out[j].append(i)
            idx+=1
    return out
rng=random.Random(3); U=set()
for n in range(2,6):
    for bits in range(1<<(n*(n-1)//2)): U.add(ham(n,tour(n,bits)))
for n,ns in [(6,30000),(7,60000),(8,120000)]:
    pr=n*(n-1)//2
    for _ in range(ns): U.add(ham(n,tour(n,rng.getrandbits(pr))))
    print(f"done n={n}",flush=True)
def isprime(x):
    if x<2: return False
    if x==2: return True
    if x%2==0: return False
    d=3
    while d*d<=x:
        if x%d==0: return False
        d+=2
    return True
forb=[k for k in range(1,150,2) if k not in U]
print("FORBIDDEN odd <=149 (n<=8 sample):", forb, flush=True)
print("  of these, PRIME:", [x for x in forb if isprime(x)], flush=True)
print("  of these, mult of 7:", [x for x in forb if x%7==0], flush=True)
print("achieved mult-of-7 in [35,149]:", sorted(x for x in U if x%7==0 and 35<=x<150), flush=True)
print("every odd prime in [3,149] achieved except 7:",
      all((x in U) for x in range(3,150,2) if isprime(x) and x!=7), flush=True)
print("max sampled value:", max(U), flush=True)
