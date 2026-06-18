"""
lrc14_tourmap_diophantine_control_kps-S2-wf.py  (FAST version)

CONTROL EXPERIMENT for the diophantine tournament maps.

Census found M2 (forbids 1 class at n=5) and M5 (forbids 9 at n=5).
KEY QUESTION: is the forbidding caused by the LRC/loneliness constraint on tau,
or is it INTRINSIC to the arc rule (the rule can't express those classes for ANY
tau)?

We compare two input families over the SAME speed sets:
  (A) LRC-constrained: tau in cand(S) (binding crossings) + the argmax lonely tau.
  (B) UNCONSTRAINED control: arbitrary tau = a/Q for a wide grid (Q=2..Qmax, all a).

If a class is forbidden under (A) AND (B): RULE ARTIFACT (not a loneliness signal).
If realized under (B) but forbidden under (A): LONELINESS truly forbids it.

SPEED: precompute, for each n, a lookup  adjacency-bitmask -> canonical class id,
so the inner loop is O(1) per tournament (no n! canonicalization in the hot loop).
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, product

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def signed_frac(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else r - 1
def nrm_frac01(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r
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
def g(S, t): return min(nrm(v*t) for v in S)
def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at
def gcd_list(xs):
    g0 = 0
    for x in xs: g0 = gcd(g0, x)
    return g0
def speed_sets(n, maxspeed):
    return [S for S in combinations(range(1, maxspeed+1), n) if gcd_list(S) == 1]

# ---- Fast tournament encoding ----
# Encode a tournament on 0..n-1 as a bitmask over ordered pairs (i<j):
#   bit p (p indexes the pair list) = 1 if i->j, 0 if j->i.
PAIRS = {n: list(combinations(range(n), 2)) for n in (3,4,5)}

def arcs_to_mask(arcs, n):
    pairs = PAIRS[n]
    idx = {pr: k for k, pr in enumerate(pairs)}
    mask = 0
    for (i, j) in arcs:
        a, b = (i, j) if i < j else (j, i)
        k = idx[(a, b)]
        if (i, j) == (a, b):   # i->j with i<j  => bit 1
            mask |= (1 << k)
    return mask

def mask_to_adj(mask, n):
    pairs = PAIRS[n]
    A = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(pairs):
        if (mask >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    return A

def adj_num_3cyc(A, n):
    c = 0
    for a, b, cc in combinations(range(n), 3):
        outs = {a:0, b:0, cc:0}
        for (x, y) in combinations([a,b,cc], 2):
            if A[x][y] == 1: outs[x]+=1
            else: outs[y]+=1
        if sorted(outs.values()) == [1,1,1]: c += 1
    return c

def adj_score(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def relabel_mask(mask, perm, n):
    A = mask_to_adj(mask, n)
    B = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            B[perm[i]][perm[j]] = A[i][j]
    # encode B back
    pairs = PAIRS[n]; m2 = 0
    for k, (i, j) in enumerate(pairs):
        if B[i][j] == 1:
            m2 |= (1 << k)
    return m2

def build_canon_table(n):
    """Return (mask->classid, classid->(c3,score), num_classes)."""
    perms = list(permutations(range(n)))
    npairs = len(PAIRS[n])
    canon_of = {}
    class_info = {}
    next_id = 0
    for mask in range(1 << npairs):
        if mask in canon_of:
            continue
        # find canonical (minimum) mask in orbit
        orbit = set()
        for perm in perms:
            orbit.add(relabel_mask(mask, perm, n))
        rep = min(orbit)
        if rep not in class_info:
            A = mask_to_adj(rep, n)
            class_info[rep] = (adj_num_3cyc(A, n), adj_score(A, n))
        cid = class_info[rep]
        for mm in orbit:
            canon_of[mm] = rep
    # assign ids
    reps = sorted(set(canon_of.values()))
    repid = {r:i for i,r in enumerate(reps)}
    mask_to_id = {m: repid[canon_of[m]] for m in canon_of}
    id_info = {repid[r]: class_info[r] for r in reps}
    return mask_to_id, id_info, len(reps)

# ---- the two forbidding maps ----
def method2(S, tau):
    n = len(S)
    dist = [nrm(v*tau) for v in S]
    side = [1 if signed_frac(v*tau) > 0 else (-1 if signed_frac(v*tau) < 0 else 0) for v in S]
    q = [v * nrm(v*tau) for v in S]
    arcs = []
    for i in range(n):
        for j in range(i+1, n):
            if side[i] == side[j]:
                if dist[i] == dist[j]: return None
                arcs.append((i, j) if dist[i] < dist[j] else (j, i))
            else:
                if q[i] == q[j]: return None
                arcs.append((i, j) if q[i] < q[j] else (j, i))
    return arcs
def method5(S, tau):
    n = len(S)
    third = []
    for v in S:
        f = nrm_frac01(v*tau)
        third.append(0 if f < F(1,3) else (1 if f < F(2,3) else 2))
    dist = [nrm(v*tau) for v in S]
    arcs = []
    for i in range(n):
        for j in range(i+1, n):
            ti, tj = third[i], third[j]
            if ti == tj:
                if dist[i] == dist[j]: return None
                arcs.append((i, j) if dist[i] < dist[j] else (j, i))
            else:
                arcs.append((i, j) if (ti - tj) % 3 == 2 else (j, i))
    return arcs

def lrc_ids(method, n, maxspeed, mask_to_id):
    seen = set()
    for S in speed_sets(n, maxspeed):
        taus = set(cand(S)); taus.add(M(S)[1])
        for tau in taus:
            arcs = method(S, tau)
            if arcs is None: continue
            mask = arcs_to_mask(arcs, n)
            seen.add(mask_to_id[mask])
    return seen

def unc_ids(method, n, maxspeed, Qmax, mask_to_id):
    seen = set()
    for S in speed_sets(n, maxspeed):
        for Q in range(2, Qmax+1):
            for a in range(1, Q):
                tau = F(a, Q)
                arcs = method(S, tau)
                if arcs is None: continue
                mask = arcs_to_mask(arcs, n)
                seen.add(mask_to_id[mask])
    return seen

def main():
    print("CONTROL: loneliness-forbidding vs rule-artifact forbidding")
    print("Building canonical tables...")
    tables = {}
    for n in (4, 5):
        mt, info, nc = build_canon_table(n)
        tables[n] = (mt, info, nc)
        print(f"  n={n}: {nc} iso classes")
    print()

    for name, method in (("M2_side_then_quality", method2),
                         ("M5_thirds_rockpaper", method5)):
        print("="*66); print(name); print("="*66)
        for n in (4, 5):
            ms = 12 if n == 5 else 14
            mt, info, nc = tables[n]
            lrc = lrc_ids(method, n, ms, mt)
            unc = unc_ids(method, n, ms, 30, mt)
            allids = set(range(nc))
            forb_lrc = allids - lrc
            forb_unc = allids - unc
            lonely_forbidden = forb_lrc - forb_unc
            rule_forbidden = forb_lrc & forb_unc
            print(f" n={n}: free={nc}  LRC-realized={len(lrc)}  "
                  f"unconstrained-realized={len(unc)}")
            print(f"      forbidden under LRC: {len(forb_lrc)} | "
                  f"RULE-artifact (forbidden even unconstrained): {len(rule_forbidden)} | "
                  f"LONELINESS-forbidden (rule CAN, loneliness CANNOT): {len(lonely_forbidden)}")
            for cid in sorted(lonely_forbidden):
                c3, ss = info[cid]
                print(f"        >>> LONELINESS-FORBIDDEN: #3cyc={c3}, scores={ss}")
            for cid in sorted(rule_forbidden):
                c3, ss = info[cid]
                print(f"        (rule-artifact: #3cyc={c3}, scores={ss})")
        print()

if __name__ == "__main__":
    main()
