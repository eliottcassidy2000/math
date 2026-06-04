import itertools, random
def make_beats(n,bits):
    b=[[False]*n for _ in range(n)]; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if (bits>>idx)&1: b[i][j]=True
            else: b[j][i]=True
            idx+=1
    return b
def ham_paths(n,beats,verts=None):
    if verts is None: verts=list(range(n))
    k=len(verts); idx={v:i for i,v in enumerate(verts)}; full=(1<<k)-1
    dp=[[0]*k for _ in range(1<<k)]
    for i in range(k): dp[1<<i][i]=1
    for mask in range(1<<k):
        for last in range(k):
            c=dp[mask][last]
            if not c: continue
            for nx in range(k):
                if not (mask>>nx)&1 and beats[verts[last]][verts[nx]]:
                    dp[mask|(1<<nx)][nx]+=c
    return sum(dp[full][i] for i in range(k))
def sccs(n,beats):
    # strongly connected components via reachability (small n)
    reach=[[i==j for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            if beats[i][j]: reach[i][j]=True
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if reach[i][k] and reach[k][j]: reach[i][j]=True
    seen=set(); comps=[]
    for i in range(n):
        if i in seen: continue
        comp=[j for j in range(n) if reach[i][j] and reach[j][i]]
        comps.append(comp); seen|=set(comp)
    return comps
def is_strong(n,beats): return len(sccs(n,beats))==1

# 1) verify H(T) = product of H over SCCs
print("=== H(T) = product of H(SCC_i)  (the SCC equidecomposition; H multiplicative) ===")
random.seed(0); ok=True
for n in range(2,8):
    m=n*(n-1)//2
    for _ in range(20000):
        beats=make_beats(n,random.getrandbits(m))
        H=ham_paths(n,beats)
        prod=1
        for comp in sccs(n,beats): prod*=ham_paths(n,beats,comp)
        if H!=prod: ok=False; print("FAIL",n); break
print("  H = prod H(SCC) holds on 20000 random tournaments per n=2..7:", ok)

# 2) ATOM spectrum: H values of STRONGLY-CONNECTED tournaments
print("\n=== ATOM spectrum: H(strongly-connected tournament), by n ===")
atoms=set([1])  # single vertex
for n in range(2,8):
    m=n*(n-1)//2; vals=set()
    if m<=15:
        for bits in range(1<<m):
            beats=make_beats(n,bits)
            if is_strong(n,beats): vals.add(ham_paths(n,beats))
    else:
        random.seed(2)
        for _ in range(800000):
            beats=make_beats(n,random.getrandbits(m))
            if is_strong(n,beats): vals.add(ham_paths(n,beats))
    atoms|=vals
    print(f"  n={n} ({'exh' if m<=15 else 'samp'}): strongly-connected H-values = {sorted(vals)}")

# 3) achievable = multiplicative closure of atoms; forbidden = complement
print("\n=== achievable = multiplicative closure of atoms; FORBIDDEN = odd integers NOT a product of atoms ===")
LIM=200
A=sorted(a for a in atoms if a<=LIM)
reach={1}
changed=True
while changed:
    changed=False
    for v in list(reach):
        for a in A:
            if v*a<=LIM and v*a not in reach: reach.add(v*a); changed=True
forbidden=[h for h in range(1,LIM+1,2) if h not in reach]
print(f"  atoms (<= {LIM}): {A}")
print(f"  FORBIDDEN odd values <= {LIM} (not a product of atoms): {forbidden}")
print(f"  7*3^k = {[7*3**k for k in range(4) if 7*3**k<=LIM]}")
print(f"  is the forbidden set exactly {{7,21,63,...}} = 7*3^k (plus maybe high sparse)? first few forbidden:", forbidden[:12])
