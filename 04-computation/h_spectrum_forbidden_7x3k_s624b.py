import itertools, random
def make_beats(n,bits):
    beats=[[False]*n for _ in range(n)]; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if (bits>>idx)&1: beats[i][j]=True
            else: beats[j][i]=True
            idx+=1
    return beats
def ham_paths(n,beats):
    # count directed Hamiltonian paths via subset DP
    full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for last in range(n):
            c=dp[mask][last]
            if not c: continue
            for nxt in range(n):
                if not (mask>>nxt)&1 and beats[last][nxt]:
                    dp[mask|(1<<nxt)][nxt]+=c
    return sum(dp[full][v] for v in range(n))
def spectrum(n, exhaustive, samples=200000):
    m=n*(n-1)//2; vals=set()
    if exhaustive:
        for bits in range(1<<m): vals.add(ham_paths(n,make_beats(n,bits)))
    else:
        random.seed(1)
        for _ in range(samples): vals.add(ham_paths(n,make_beats(n,random.getrandbits(m))))
    return vals
print("=== H = #Hamiltonian paths (Redei; = I(Omega,2)); achievable spectrum + ODD GAPS ===")
allv={}
for n in range(3,8):
    ex=(n*(n-1)//2)<=15
    v=spectrum(n,ex, 300000); allv[n]=v; mx=max(v)
    gaps=[h for h in range(1,mx+1,2) if h not in v]
    print(f"  n={n} ({'exhaustive' if ex else 'sampled'}): maxH={mx}; odd GAPS below max: {gaps}")
U=set().union(*allv.values()); mx=max(U)
gaps=[h for h in range(1,min(mx,100)+1,2) if h not in U]
print(f"\n  union(n<=7) odd gaps <= {min(mx,100)}: {gaps}")
print(f"  is 7 a gap? {7 not in U}   is 21 a gap? {21 not in U}   is 63 a gap? {63 not in U}")
print(f"  7*3^k = {[7*3**k for k in range(3)]}")
