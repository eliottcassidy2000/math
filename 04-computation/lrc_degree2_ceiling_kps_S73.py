#!/usr/bin/env python3
r"""
lrc_degree2_ceiling_kps_S73.py   (kind-pasteur-2026-07-07-S73, HYP-5147 part a)  v2

THE DEGREE-2 CEILING OF THE (A') TAIL: what uniform floor on mu_{1/7}(E_k) can ANY
argument consuming only PAIRWISE data prove?

Degree = number of ELEMENTS involved.  Degree-2 data = the joint law of one phase pair
(frac(e_i x), frac(e_j x)) -- a line measure; by THM-641 its window masses converge to
the PRODUCT-UNIFORM values as the reduced ratio grows, so product data is the legal
generic limit at any fixed window resolution.  A uniform degree-2 lower bound on the
tail must hold at that data.  We exhibit COUPLINGS of k circle variables whose pairwise
DIFFERENCE law is uniform at 14-adic resolution (the resolution of every tool in play)
and whose tail P(maxgap > 1/7) is minimal -- the LP value is the ceiling.

Couplings: mixtures over (unlabeled grid config C on Z_G, uniform rotation, uniform
labeling).  Pair-difference law = sum_C w_C A_C(cell)/(k(k-1)).  Constraint: flat over
the 14 difference cells ((c/14, (c+1)/14], with exact-coincidence d=0 as a 15th
coordinate that must vanish in the continuum-legal run).  Objective: min tail weight.

Grids: G=14 (coincidence pigeonhole appears: safe classes carry at most one doubled
bin) and G=28 (twin pairs at 1/28 supply the near-zero difference mass legally).

Safe classes enumerated via GAP-WORD necklaces: compositions of G into m parts, each
part <= G/7 (bin-center gap G/7 = 1/7 exactly => maxgap NOT > 1/7, safe), m occupied
bins, k-m extra points as multiplicities.  G=28,k=13 uses a curated zoo (any feasible
subset upper-bounds the LP value).
"""
import itertools, math
import numpy as np
from scipy.optimize import linprog

def gapword_necklaces(G, m, maxpart):
    """compositions of G into m parts each in [1, maxpart], up to rotation."""
    seen = set(); out = []
    def rec(parts, rem, slots):
        if slots == 0:
            if rem == 0:
                w = tuple(parts)
                canon = min(w[r:] + w[:r] for r in range(m))
                if canon not in seen:
                    seen.add(canon); out.append(canon)
            return
        lo = max(1, rem - maxpart * (slots - 1))
        hi = min(maxpart, rem - (slots - 1))
        for p in range(lo, hi + 1):
            rec(parts + [p], rem - p, slots - 1)
    rec([], G, m)
    return out

def safe_classes(G, k, cap_mult=None):
    """occupancy tuples (canonical under rotation) of safe configs: k points, occupied
    gaps <= G//7.  cap_mult: max total classes with multiplicities (None = all)."""
    maxpart = G // 7
    classes = []
    m_min = -(-G // maxpart)
    for m in range(m_min, min(G, k) + 1):
        words = gapword_necklaces(G, m, maxpart)
        n_extra = k - m
        for w in words:
            # positions from gap word: bin 0, then cumulative
            pos = [0]
            for g in w[:-1]:
                pos.append(pos[-1] + g)
            # multiplicity placements: n_extra extra points on m bins (combinations with repetition)
            if n_extra == 0:
                placements = [()]
            else:
                placements = itertools.combinations_with_replacement(range(m), n_extra)
            cnt = 0
            for extra in placements:
                occ = [0] * G
                for b in pos:
                    occ[b] = 1
                for e in extra:
                    occ[pos[e]] += 1
                canon = min(tuple(occ[(i + r) % G] for i in range(G)) for r in range(G))
                classes.append(canon)
                cnt += 1
                if cap_mult and cnt >= cap_mult:
                    break
    return list(dict.fromkeys(classes))

def tail_zoo(G, k):
    import random as _rnd
    rng = _rnd.Random(73)
    zoo = []
    occ = [0] * G; occ[0] = k; zoo.append(tuple(occ))                      # all coincident
    occ = [0] * G; occ[0] = k // 2; occ[1] = k - k // 2; zoo.append(tuple(occ))
    for spread in range(2, G - 2):                                          # k points in a block
        occ = [0] * G
        for j in range(k):
            occ[round(j * spread / (k - 1))] += 1
        zoo.append(tuple(occ))
    for hole in range(G // 7 + 1, G // 2 + 1):                              # even spread + one big hole
        rest = G - hole
        occ = [0] * G
        for j in range(k):
            occ[(hole + round(j * (rest - 1) / max(1, k - 1))) % G] += 1
        zoo.append(tuple(occ))
    # rich random multisets (profiles span the tail polytope well)
    for _ in range(4000):
        occ = [0] * G
        for j in range(k):
            occ[rng.randrange(G)] += 1
        # ensure it IS a tail config (some occupied-gap > G//7)
        pos = [b for b in range(G) if occ[b]]
        gmax = max(((pos[(i + 1) % len(pos)] - pos[i]) % G) or G for i in range(len(pos))) if len(pos) > 1 else G
        if gmax > G // 7:
            zoo.append(tuple(occ))
    # twin-rich tails: t twin pairs evenly spread + remaining points clustered
    for t in range(1, k // 2 + 1):
        occ = [0] * G
        for j in range(t):
            b = (j * G) // t
            occ[b] += 1; occ[(b + 1) % G] += 1
        left = k - 2 * t
        for j in range(left):
            occ[(2 + j) % G] += 1
        zoo.append(tuple(occ))
    return list(dict.fromkeys(zoo))

def diff_profile_cells(occ, G):
    """15-vector: 14-adic difference cells (c/14,(c+1)/14] for d!=0, [14]=coincidence."""
    A = [0] * 15
    for d in range(G):
        cnt = 0
        for b in range(G):
            if occ[b]:
                cnt += occ[b] * occ[(b + d) % G] if d != 0 else occ[b] * (occ[b] - 1)
        if cnt == 0:
            continue
        if d == 0:
            A[14] += cnt
        else:
            c = math.ceil(d * 14 / G) - 1
            A[c] += cnt
    return A

def lp_ceiling(G, k, cap_mult=None, curated_safe=None, label=""):
    safe = curated_safe if curated_safe is not None else safe_classes(G, k, cap_mult)
    tails = tail_zoo(G, k)
    classes = [(c, True) for c in safe] + [(c, False) for c in tails]
    profs = [(diff_profile_cells(c, G), is_safe) for c, is_safe in classes]
    kk = k * (k - 1)
    for fold in ([False, True] if G == 14 else [False]):
        n = len(profs)
        Aeq = np.zeros((14, n)); c_obj = np.zeros(n)
        bounds = [None] * n
        for j, (A, is_safe) in enumerate(profs):
            cells = list(A[:14])
            if fold:
                cells[0] += A[14]
            for i in range(14):
                Aeq[i, j] = cells[i]
            c_obj[j] = 0.0 if is_safe else 1.0
            bounds[j] = (0, 0) if (A[14] > 0 and not fold) else (0, None)
        beq = np.full(14, kk / 14.0)
        res = linprog(c_obj, A_eq=Aeq, b_eq=beq, bounds=bounds, method="highs")
        tag = "grid-legal (coincidence folded)" if fold else "continuum-legal (0 coincidence)"
        if res.status == 0:
            print(f"    G={G} k={k} {label}[{tag}]: ceiling <= {res.fun:.6f}   "
                  f"({len(safe)} safe / {len(tails)} tail classes)")
        else:
            print(f"    G={G} k={k} {label}[{tag}]: LP infeasible over this class set "
                  f"({len(safe)} safe) -- enlarge zoo")

print("=" * 96)
print("THE DEGREE-2 CEILING: min tail weight of pairwise-uniform couplings (14-adic data)")
print("=" * 96)
HONEST_T = {8: 1702763 / 2522520, 13: 14249 / 252252}
print(f"\n  k=8   (honest bar T_8 = {HONEST_T[8]:.4f}, MISTAKE-123)")
lp_ceiling(14, 8)
lp_ceiling(28, 8)
print(f"\n  k=13  (honest bar T_13 = {HONEST_T[13]:.4f})")
lp_ceiling(14, 13, cap_mult=200)
# curated G=28 k=13 safe zoo: gap words of 1..4 summing 28 with m in {13,14... m>=7};
# use m=13,14 twin-rich and spread patterns only (feasible subset suffices)
cur = []
G = 28
for w in gapword_necklaces(G, 13, 4)[:4000]:
    pos = [0]
    for g in w[:-1]:
        pos.append(pos[-1] + g)
    occ = [0] * G
    for b in pos:
        occ[b] = 1
    canon = min(tuple(occ[(i + r) % G] for i in range(G)) for r in range(G))
    cur.append(canon)
cur = list(dict.fromkeys(cur))
print(f"    (curated m=13 gap-word zoo: {len(cur)} classes)")
lp_ceiling(28, 13, curated_safe=cur, label="curated ")
print()
print("READING: the ceiling is the max uniform floor ANY pairwise-only argument can give.")
print("Compare: honest bars T_8 = 0.675, T_13 = 0.0565; degree-3 endpoint-Hunter ceiling")
print("3/7 (lrc_adaptive_hunter_kps_S73).  Anchor-event routes (2-anchor PA_2, CE-on-")
print("anchors) are k-body and NOT graded by this ceiling.")
