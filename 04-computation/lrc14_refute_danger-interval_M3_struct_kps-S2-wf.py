# lrc14_refute_danger-interval_M3_struct_kps-S2-wf.py
# Part 2: (a) m=13 induced case on the ACTUAL hard LRC families + random sample.
#         (b) STRUCTURAL explanation: M3 arc(i,j) = [ttd_i<ttd_j] XOR [parity(S_i+S_j) even].
#             Enumerate every achievable (scalar-order, speed-parity-pattern) and show the
#             realizable iso-class set is exactly 4 at n=5, independent of the LRC phase.
#
# Exact Fractions for all decisions.

from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, product
from functools import reduce
import sys, random

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def g(S, t): return min(nrm(v * t) for v in S)
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
def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at
def frac(x):
    r = x - int(x); return r + 1 if r < 0 else r

def time_to_danger_signed(v, tau0):
    p = frac(v * tau0)
    dp = F(13, 14) - p
    if dp < 0: dp += 1
    if dp < 0: dp += 1
    return dp / v
def method3_adj(S, tau0):
    m = len(S)
    ttd = [time_to_danger_signed(S[i], tau0) for i in range(m)]
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j: continue
            base = ttd[i] < ttd[j]
            if ttd[i] == ttd[j]: base = (S[i] < S[j])
            if (S[i] + S[j]) % 2 == 0: base = not base
            adj[i][j] = base
    return adj

def canon_key(adj, m):
    best = None
    for p in permutations(range(m)):
        bits = tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)
        if best is None or bits < best: best = bits
    return best
def h_count(adj, m):
    return sum(1 for p in permutations(range(m)) if all(adj[p[k]][p[k+1]] for k in range(m-1)))
def score_seq(adj, m):
    return tuple(sorted(sum(1 for j in range(m) if j != i and adj[i][j]) for i in range(m)))
def num_3cycles(adj, m):
    c = 0
    for a, b, cc in combinations(range(m), 3):
        if adj[a][b] and adj[b][cc] and adj[cc][a]: c += 1
        if adj[a][cc] and adj[cc][b] and adj[b][a]: c += 1
    return c
def all_iso(m):
    seen = {}
    pairs = list(combinations(range(m), 2))
    for bits in range(2**len(pairs)):
        adj = [[False]*m for _ in range(m)]
        for idx, (i, j) in enumerate(pairs):
            if (bits >> idx) & 1: adj[i][j] = True
            else: adj[j][i] = True
        kk = canon_key(adj, m)
        if kk not in seen: seen[kk] = (h_count(adj, m), num_3cycles(adj, m), score_seq(adj, m))
    return seen
def primitive(S): return reduce(gcd, S) == 1

FREE5 = all_iso(5)
# IMPORTANT: invariant tuple (H,c3,score) is NOT a complete iso invariant at n=5.
# There are TWO distinct classes with tuple (9,3,(1,1,2,3,3)). The claim realizes ONE
# and forbids the OTHER. So we MUST key on canonical iso keys, not tuples.
# First compute the BASELINE realized canonical keys via the actual M3 map at the optimum.
def baseline_realized_keys(vmax=20):
    keys = {}
    for combo in combinations(range(1, vmax+1), 5):
        if not primitive(combo): continue
        _, tau0 = M(list(combo))
        if tau0 is None: continue
        adj = method3_adj(list(combo), tau0)
        ck = canon_key(adj, 5)
        if ck not in keys:
            keys[ck] = (h_count(adj,5), num_3cycles(adj,5), score_seq(adj,5))
    return keys
BASE_REAL_KEYS = baseline_realized_keys(20)
BASELINE_FORBIDDEN_KEYS = set(FREE5) - set(BASE_REAL_KEYS)
print(f"[setup] baseline M3-realized canonical keys: {len(BASE_REAL_KEYS)}/12")
print(f"[setup] baseline-forbidden canonical keys:   {len(BASELINE_FORBIDDEN_KEYS)}/12")
for k in sorted(BASELINE_FORBIDDEN_KEYS, key=lambda k:(FREE5[k][0],FREE5[k][1])):
    print(f"        forbidden-key tuple {FREE5[k]}")
sys.stdout.flush()

print("PART 2: m=13 induced + STRUCTURAL refutation attempt")
print("="*78); sys.stdout.flush()

# ---------------------------------------------------------------------------
# (a) m=13 induced from ACTUAL LRC families + random covering/tight 13-sets.
# ---------------------------------------------------------------------------
print("\n(a) m=13 induced sub-tournaments from real LRC-13 speed sets")
fams = {
    "AP {1..13}": tuple(range(1, 14)),
    "GW {1..11,13,24}": tuple(sorted(set(list(range(1,12))+[13,24]))),
    "cover {1..11,13,84}": tuple(sorted(set(list(range(1,12))+[13,84]))),
    "AP-drop6+28": tuple(sorted(set([1,2,3,4,5,7,8,9,10,11,12,13,28]))),
    "cover {1..11,13,168}": tuple(sorted(set(list(range(1,12))+[13,168]))),
}
realized_induced = {}
witnesses = []
for name, S13 in fams.items():
    if len(set(S13)) != 13: continue
    Sl = list(S13)
    _, tau0 = M(Sl)
    adj = method3_adj(Sl, tau0)
    for sub in combinations(range(13), 5):
        subadj = [[adj[sub[a]][sub[b]] for b in range(5)] for a in range(5)]
        canon = canon_key(subadj,5)
        ck = (h_count(subadj,5), num_3cycles(subadj,5), score_seq(subadj,5))
        if canon not in realized_induced:
            realized_induced[canon] = (ck, tuple(Sl[x] for x in sub), tau0)
        if canon in BASELINE_FORBIDDEN_KEYS:   # TRUE forbidden CLASS (canonical), not tuple
            witnesses.append((name, ck, tuple(Sl[x] for x in sub), tau0))
print(f"   from {len(fams)} named LRC families: realized {len(realized_induced)}/12")
inv_real = set(v[0] for v in realized_induced.values())
print(f"   invariant tuples realized: {sorted(inv_real)}")
print(f"   TRUE forbidden-CLASS hits as induced sub? {len(witnesses)} witnesses")
for w in witnesses[:5]: print("     WITNESS:", w)
sys.stdout.flush()

# random covering 13-sets: C = {1..11,13} U {w}, w in 14*Z (parked runner), plus random
random.seed(424242)
rand_real = dict(realized_induced)
rand_wit = list(witnesses)
core = list(range(1,12)) + [13]
for trial in range(400):
    w = 14 * random.randint(1, 40)
    S13 = tuple(sorted(set(core + [w])))
    if len(S13) != 13: continue
    Sl = list(S13)
    _, tau0 = M(Sl)
    adj = method3_adj(Sl, tau0)
    # sample 60 random 5-subsets per set to bound cost
    subs = list(combinations(range(13), 5))
    for sub in random.sample(subs, min(80, len(subs))):
        subadj = [[adj[sub[a]][sub[b]] for b in range(5)] for a in range(5)]
        canon = canon_key(subadj,5)
        ck = (h_count(subadj,5), num_3cycles(subadj,5), score_seq(subadj,5))
        if canon not in rand_real:
            rand_real[canon] = (ck, tuple(Sl[x] for x in sub), tau0)
        if canon in BASELINE_FORBIDDEN_KEYS:
            rand_wit.append(("randcover", ck, tuple(Sl[x] for x in sub), tau0))
print(f"\n   + 400 random covering 13-sets (parked w in 14Z): realized {len(rand_real)}/12")
print(f"   TRUE forbidden-CLASS hit total? {len(rand_wit)} witnesses")
for w in rand_wit[:5]: print("     WITNESS:", w)
sys.stdout.flush()

# ---------------------------------------------------------------------------
# (b) STRUCTURAL: the M3 tournament depends on tau0 and S ONLY through
#     (i)  the RANK ORDER of ttd values  (a permutation / total preorder), and
#     (ii) the parity pattern  par_{ij} = (S_i+S_j) mod 2  (determined by parities of S_i).
#     Enumerate ALL (linear order on 5 vertices) x (parity vector p in {0,1}^5 for the
#     vertices, par_ij = p_i XOR p_j XOR 1==even... wait (S_i+S_j) even <=> S_i,S_j same parity).
#     => par_ij "even" iff p_i==p_j. Arc reverses when p_i==p_j.
#     This is EXACTLY enumerable: 5! orders x 2^5 parity assignments = 120*32 = 3840 cases,
#     but order matters only up to the induced tournament. Enumerate -> realized classes.
# ---------------------------------------------------------------------------
print("\n(b) STRUCTURAL enumeration of ALL M3-reachable tournaments at n=5")
print("    arc(i,j) = [rank_i < rank_j] XOR [p_i == p_j]   (p = parity of speed)")
struct_real = {}
struct_inv = set()
# rank: a permutation giving strict order of ttd (ties broken by speed, but generic = strict).
# We enumerate all strict total orders (5! ) and all parity assignments (2^5).
for order in permutations(range(5)):
    rank = [0]*5
    for pos, v in enumerate(order): rank[v] = pos
    for pbits in range(32):
        p = [(pbits >> k) & 1 for k in range(5)]
        adj = [[False]*5 for _ in range(5)]
        for i in range(5):
            for j in range(5):
                if i == j: continue
                base = rank[i] < rank[j]
                if p[i] == p[j]:   # (S_i+S_j) even -> reverse
                    base = not base
                adj[i][j] = base
        ck = canon_key(adj, 5)
        if ck not in struct_real:
            struct_real[ck] = (h_count(adj,5), num_3cycles(adj,5), score_seq(adj,5))
        struct_inv.add((h_count(adj,5), num_3cycles(adj,5), score_seq(adj,5)))
print(f"    structurally reachable iso classes: {len(struct_real)}/12")
print(f"    reachable invariant tuples: {sorted(struct_inv)}")
struct_real_keys = set(struct_real)
struct_forbidden_keys = set(FREE5) - struct_real_keys
print(f"    structurally FORBIDDEN CLASSES ({len(struct_forbidden_keys)}) [by canon key]:")
for k in sorted(struct_forbidden_keys, key=lambda k:(FREE5[k][0],FREE5[k][1])):
    print("       ", FREE5[k])
# Does structural reachable set hit any baseline-forbidden CLASS (canonical key)?
hit = struct_real_keys & BASELINE_FORBIDDEN_KEYS
print(f"\n    baseline-forbidden CLASSES reachable structurally: {len(hit)}")
for k in sorted(hit, key=lambda k:(FREE5[k][0],FREE5[k][1])): print("       REACHED:", FREE5[k])
# Cross-check: is the structural reachable key set EXACTLY the baseline realized key set?
print(f"    structural-reachable == baseline-realized key sets? "
      f"{struct_real_keys == set(BASE_REAL_KEYS)}")
sys.stdout.flush()

# ---------------------------------------------------------------------------
# VERDICT
# ---------------------------------------------------------------------------
print("\n" + "="*78)
print("VERDICT (part 2)  [keyed on CANONICAL iso class, not invariant tuple]")
all_wit = rand_wit
struct_breaks = hit
print(f"  m=13 induced TRUE-forbidden-class witnesses: {len(all_wit)}")
print(f"  structural enumeration reaching a forbidden CLASS: {len(struct_breaks)}")
if all_wit or struct_breaks:
    print("  >>> CLAIM REFUTED (a genuinely forbidden iso CLASS was realized) <<<")
else:
    print("  >>> forbidden set UNBROKEN: M3-reachable iso classes = exactly the 4 baseline ones <<<")
    print("  Reachable canonical classes (tuples):")
    for k in sorted(set(BASE_REAL_KEYS), key=lambda k:(FREE5[k][0],FREE5[k][1])):
        print("     ", FREE5[k])
sys.stdout.flush()
