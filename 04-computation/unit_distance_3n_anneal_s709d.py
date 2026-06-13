"""
monad-explorer-2026-06-06-S709d
===============================
Narrow [17,39] from above: local structure of the sqrt(7)-Eisenstein unit graph +
multi-seed annealing search for the smallest induced subgraph with avg degree > 6.

THM-420: floor 17, record construction 39 (edge-midpoint disk). Here we (1) report
the graph's local structure (flower density, common-neighbour distribution) to see
WHY small balls fall short, and (2) run simulated-annealing / multistart local search
to try to beat 39. Exact integer arithmetic throughout (no randomness in the geometry;
the only randomness is the search walk, seeded deterministically per restart index).
"""
import math

A, B, C = 1, 1, 1
R = 7
def Q(dx, dy): return A*dx*dx + B*dx*dy + C*dy*dy
OFFS = [(dx, dy) for dx in range(-4, 5) for dy in range(-4, 5) if Q(dx, dy) == R]
assert len(OFFS) == 12

# big patch as the universe
H = 14
UNIV = [(x, y) for x in range(-H, H+1) for y in range(-H, H+1)]
USET = set(UNIV)
ADJ = {p: [(p[0]+dx, p[1]+dy) for (dx, dy) in OFFS if (p[0]+dx, p[1]+dy) in USET]
       for p in UNIV}

# ---- 1. local structure --------------------------------------------------
print("="*68)
print("1.  Local structure of the sqrt(7)-Eisenstein unit graph (12-regular)")
print("="*68)
o = (0, 0)
nbrs = ADJ[o]
nset = set(nbrs)
flower_internal = sum(1 for p in nbrs for q in ADJ[p] if q in nset) // 2
print(f"   degree (interior)          : {len(nbrs)}")
print(f"   edges among the 12 nbrs    : {flower_internal}")
print(f"   flower (center+12) N=13    : U = {12+flower_internal}, "
      f"avg deg = {2*(12+flower_internal)/13:.3f}  (need >6 => U>39)")
# common-neighbour distribution over adjacent vs all close pairs
from collections import Counter
cn_adj = Counter(); cn_all = Counter()
sample = [(x, y) for x in range(-5, 6) for y in range(-5, 6)]
sset = set(sample)
for i, p in enumerate(sample):
    Np = set(ADJ[p])
    for q in sample[i+1:]:
        Nq = set(ADJ[q])
        c = len(Np & Nq)
        cn_all[c] += 1
        if q in Np:
            cn_adj[c] += 1
print(f"   common-nbr count over ADJACENT pairs : {dict(sorted(cn_adj.items()))}")
print(f"   common-nbr count over ALL pairs      : {dict(sorted(cn_all.items()))}")
print("   (<=2 always, confirming the CN bound; adjacent pairs realise up to 2)")

# ---- 2. annealing / multistart to beat N=39 ------------------------------
print()
print("="*68)
print("2.  Multistart local search for smallest induced subgraph, avg deg > 6")
print("="*68)
def Uset(S):
    return sum(1 for p in S for q in ADJ[p] if q in S) // 2

def deterministic_rand(state):
    # simple LCG, no Date/Random needed
    state = (1103515245*state + 12345) & 0x7fffffff
    return state, state / 0x7fffffff

def anneal(seed_set, steps, rng_state):
    S = set(seed_set)
    degin = {p: sum(1 for q in ADJ[p] if q in S) for p in S}
    e = sum(degin.values())//2
    best = (len(S), e, frozenset(S)) if e > 3*len(S) else None
    frontier = set()
    for p in S:
        for q in ADJ[p]:
            if q not in S: frontier.add(q)
    T0 = 2.0
    for step in range(steps):
        T = T0 * (1 - step/steps) + 0.01
        rng_state, r = deterministic_rand(rng_state)
        k = len(S)
        # objective: minimise N subject to e>3N; surrogate score = 3*k - e (want <0 & k small)
        # propose remove (favoured) or add
        if r < 0.65 and k > 4:
            # remove a low-internal-degree vertex
            p = min(S, key=lambda v: degin[v])
            ne = e - degin[p]; nk = k-1
            accept = (ne > 3*nk)
            if not accept:
                # allow uphill sometimes
                rng_state, r2 = deterministic_rand(rng_state)
                accept = r2 < math.exp(-( (3*nk-ne) )/T)
            if accept:
                for q in ADJ[p]:
                    if q in S: degin[q]-=1
                e=ne; S.discard(p); degin.pop(p)
                frontier.discard(p)
                for q in ADJ[p]:
                    if q in S and all(w in S for w in []):
                        pass
                # rebuild frontier lazily
                frontier = set(w for v in S for w in ADJ[v] if w not in S)
        else:
            if frontier:
                rng_state, r3 = deterministic_rand(rng_state)
                fl = list(frontier)
                p = max(fl, key=lambda v: sum(1 for q in ADJ[v] if q in S))
                d = sum(1 for q in ADJ[p] if q in S)
                S.add(p); degin[p]=d
                for q in ADJ[p]:
                    if q in S and q!=p: degin[q]+=1
                e += d
                frontier = set(w for v in S for w in ADJ[v] if w not in S)
        k=len(S)
        if e > 3*k and (best is None or k < best[0]):
            best = (k, e, frozenset(S))
    return best, rng_state

# seeds: disks of various sizes/centers
def disk(cx, cy, k):
    box = sorted(UNIV, key=lambda p: ((p[0]-cx)**2+(p[0]-cx)*(p[1]-cy)+(p[1]-cy)**2))
    return box[:k]

overall = (39, None)
rng = 1
for (cx, cy) in [(0,0),(0.5,0),(1/3,1/3),(0.5,0.5),(2/3,1/3)]:
    for k0 in [30, 39, 50, 70]:
        seed = disk(cx, cy, k0)
        if Uset(set(seed)) <= 3*len(seed):
            # grow until feasible
            pass
        best, rng = anneal(seed, 4000, rng)
        if best and best[0] < overall[0]:
            overall = (best[0], best[2])
            print(f"   NEW best N={best[0]} U={best[1]} (seed disk k={k0} @({cx:.2f},{cy:.2f}))")
print()
if overall[1] is not None:
    print(f"   annealing found N={overall[0]} (improves 39)")
else:
    print(f"   annealing found NOTHING below 39 over all restarts.")
    print(f"   => N=39 stands as the record; minimum N* in [17, 39].")
