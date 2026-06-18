"""
lrc14_refute_residue-orbit_M4_witness_kps-S2-wf.py

VERIFY the refutation witness for the M4 band-depth orbit-majority forbidden-class claim.

The lean ceiling found the forbidden class (H=9, #3cyc=3, score(1,1,2,3,3)) is realized
when three runners share residue 1 mod 14 (residue multiset (1,1,1,3,5)), because the
three depth-equal runners are oriented purely by the raw-speed tie-break, manufacturing a
transitive 3-chain among them that, combined with the depth-majority arcs to residues 3,5,
yields the forbidden iso class.

Here we:
 (1) Confirm (H,#3cyc,score)=(9,3,(1,1,2,3,3)) is a UNIQUE iso class on 5 vertices
     (so the signature pins down the class -- it IS the claimed-forbidden one).
 (2) Realize it with GENUINE distinct primitive integer speeds (three speeds = 1 mod 14:
     1,15,29; plus 3,5) using the EXACT method4_band_majority map (raw-speed tie-break),
     and verify the resulting tournament's iso class.
 (3) Check whether such a speed set is LRC-LONELY (gap M>=1/14) using the exact M tool,
     so the input is a legitimate LRC configuration (not just an abstract residue pattern).
 (4) Search a broad window of genuine primitive speed sets for ANY realization of each of
     the 4 forbidden classes via exact M4, to see which are reachable with real speeds.

Exact arithmetic throughout.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, product

MOD = 14
UNITS = [a for a in range(1, MOD) if gcd(a, MOD) == 1]
def depth(r):
    r %= MOD
    return min(r, MOD - r)

# ---- EXACT M4 with raw-speed tie-break (verbatim method4_band_majority semantics) ----
def m4_adj(S):
    m = len(S)
    adj = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            wi = wj = 0
            for a in UNITS:
                di = depth(S[i]*a); dj = depth(S[j]*a)
                if di > dj: wi += 1
                elif dj > di: wj += 1
            if wi > wj: adj[i][j] = 1
            elif wj > wi: adj[j][i] = 1
            else:
                if S[i] < S[j]: adj[i][j] = 1
                else: adj[j][i] = 1
    return adj

# ---- invariants ----
def canon(adj, m):
    best = None
    for p in permutations(range(m)):
        b = tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)
        if best is None or b < best: best = b
    return best
def H(adj, m):
    c = 0
    for p in permutations(range(m)):
        if all(adj[p[k]][p[k+1]] for k in range(m-1)): c += 1
    return c
def c3(adj, m):
    c = 0
    for a, b, cc in combinations(range(m), 3):
        if adj[a][b] and adj[b][cc] and adj[cc][a]: c += 1
        if adj[a][cc] and adj[cc][b] and adj[b][a]: c += 1
    return c
def score(adj, m):
    return tuple(sorted(sum(adj[i][j] for j in range(m) if j != i) for i in range(m)))
def sig(adj, m):
    return (H(adj, m), c3(adj, m), score(adj, m))
def valid(adj, m):
    return all(adj[i][j]+adj[j][i]==1 for i in range(m) for j in range(i+1,m))

# ---- EXACT M tool ----
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def gmin(S, t):
    return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def Mgap(S):
    b = F(0); at = None
    for t in cand(S):
        v = gmin(S, t)
        if v > b: b = v; at = t
    return b, at

# ================================================================
def free_iso_classes(m):
    seen = {}
    pairs = list(combinations(range(m), 2))
    for bits in product([0,1], repeat=len(pairs)):
        adj = [[0]*m for _ in range(m)]
        for (i,j), bb in zip(pairs, bits):
            if bb: adj[i][j]=1
            else: adj[j][i]=1
        k = canon(adj, m)
        if k not in seen: seen[k] = sig(adj, m)
    return seen

def main():
    m = 5
    print("="*72)
    print("(1) Is signature (9,3,(1,1,2,3,3)) a UNIQUE iso class on 5 vertices?")
    print("="*72)
    free = free_iso_classes(m)
    by_sig = {}
    for k, s in free.items():
        by_sig.setdefault(s, []).append(k)
    target = (9, 3, (1,1,2,3,3))
    print(f"   #iso classes with this signature: {len(by_sig.get(target, []))}")
    print(f"   (a count of 1 means the forbidden signature pins exactly one iso class)")
    # also confirm all 4 forbidden are unique-signature classes
    FORB = [(9,3,(1,1,2,3,3)),(13,4,(1,2,2,2,3)),(15,4,(1,2,2,2,3)),(15,5,(2,2,2,2,2))]
    for f in FORB:
        print(f"   forbidden {f}: {len(by_sig.get(f, []))} iso class(es) share this signature")
    print()

    print("="*72)
    print("(2)+(3) GENUINE-SPEED witness: residues (1,1,1,3,5) via speeds {1,15,29,3,5}")
    print("="*72)
    S = [1, 15, 29, 3, 5]    # 1,15,29 == 1 mod 14 ; 3==3 ; 5==5
    res = [v % MOD for v in S]
    print(f"   speeds S = {S}, residues mod 14 = {res}, primitive gcd={__import__('math').gcd.__call__ and 1}")
    g = 0
    for v in S: g = gcd(g, v)
    print(f"   gcd(S) = {g} (primitive iff 1)")
    adj = m4_adj(S)
    print(f"   M4 tournament valid={valid(adj,m)}  signature={sig(adj,m)}")
    print(f"   -> matches forbidden {target}? {sig(adj,m)==target}")
    gap, at = Mgap(S)
    print(f"   exact LRC gap M(S) = {gap} = {float(gap):.6f}   (1/14 = {float(F(1,14)):.6f})")
    print(f"   lonely (M>=1/14)? {gap >= F(1,14)}   optimal tau = {at}")
    print()

    print("="*72)
    print("(4) BROAD genuine-speed search: which forbidden classes appear with REAL speeds")
    print("    (primitive 5-sets, speeds 1..VMAX), and are they LRC-lonely?")
    print("="*72)
    def primitive(s):
        gg = 0
        for v in s: gg = gcd(gg, v)
        return gg == 1
    FORBset = set(FORB)
    for VMAX in [30, 45, 60]:
        found = {}          # sig -> first speed witness
        found_lonely = {}   # sig -> first LRC-lonely (M>=1/14) speed witness
        n_used = 0
        for s in combinations(range(1, VMAX+1), m):
            if not primitive(s): continue
            n_used += 1
            adj = m4_adj(list(s))
            if not valid(adj, m): continue
            sg = sig(adj, m)
            if sg in FORBset:
                if sg not in found: found[sg] = s
                if sg not in found_lonely:
                    gp, _ = Mgap(list(s))
                    if gp >= F(1, 14): found_lonely[sg] = (s, gp)
        print(f"   VMAX={VMAX}: used {n_used} primitive 5-sets")
        print(f"     forbidden classes realized (any speed): {sorted(found.keys())}")
        for sg in sorted(found): print(f"        {sg}: witness speeds {found[sg]}")
        print(f"     forbidden classes realized by LRC-LONELY speeds (M>=1/14):")
        if found_lonely:
            for sg in sorted(found_lonely):
                s, gp = found_lonely[sg]
                print(f"        {sg}: speeds {s}, gap={gp}={float(gp):.5f}")
        else:
            print("        (none lonely at this window)")
        print()

if __name__ == "__main__":
    main()
