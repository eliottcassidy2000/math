#!/usr/bin/env python3
"""
ADVERSARIAL REFUTATION of the character-spectral M3 (Walsh-weighted lead) forbidden-class claim.

CONSTRUCTION (M3, reconstructed verbatim from the task spec):
  Vertices = n runners (speeds). Reference v0 = slowest speed.
  Rademacher weight w(a) = +1 if v0 in near-band sections {1,13} at time a/14, else -1,
      where "section of speed v at a/14" = (v*a) mod 14, and a ranges over units (Z/14)*.
  Arc x->y iff D_xy = sum_{a in units} w(a) * sign( ||x*a/14|| - ||y*a/14|| ) > 0 ; tie-break x>y.
      (D antisymmetric: D_xy = -D_yx.)

CLAIMED FORBIDDEN iso classes (under M3, over the stated exhaustive inputs):
  n=4: (H=5, #3cyc=2, score(1,1,2,2))  -- the regular 4-vertex near-3-cycle.
  n=5: 6 forbidden classes: both (9,3,(1,1,2,3,3)); all three (11/13/15,4,(1,2,2,2,3)); the regular (15,5,(2,2,2,2,2)).

GOAL: realize any "forbidden" iso class with this EXACT map over a broad LRC-constrained search.
  If realized even once -> claim REFUTED. Else CONFIRMED (with honest search-bound statement).
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
import sys, random

UNITS = [a for a in range(1,14) if gcd(a,14)==1]  # (Z/14)* = {1,3,5,9,11,13}

def nrm(x):
    # x is a Fraction; ||x|| = dist to nearest integer
    r = x - int(x)
    if r < 0: r += 1
    return r if r <= F(1,2) else 1-r

def section(v, a):
    return (v*a) % 14

def rademacher_weights(v0):
    # w(a) = +1 if v0 in near-band sections {1,13} at a/14 else -1
    return {a: (1 if section(v0, a) in (1,13) else -1) for a in UNITS}

def frac_at(v, a):
    # ||v * a/14||
    return nrm(F(v*a, 14))

def build_tournament_M3(speeds):
    """Return adjacency: adj[i][j]=1 if arc i->j. Vertices indexed by position in 'speeds' (sorted)."""
    S = sorted(set(speeds))
    n = len(S)
    v0 = S[0]  # slowest = reference
    w = rademacher_weights(v0)
    # Precompute frac for each speed at each unit
    fr = {v: {a: frac_at(v, a) for a in UNITS} for v in S}
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            x = S[i]; y = S[j]
            D = 0
            for a in UNITS:
                dx = fr[x][a]; dy = fr[y][a]
                if dx > dy: s = 1
                elif dx < dy: s = -1
                else: s = 0
                D += w[a]*s
            # arc x->y iff D_xy > 0 ; tie-break x>y  (D_xy with x as 'first', y as 'second')
            if D > 0:
                adj[i][j] = 1  # x->y
            elif D < 0:
                adj[j][i] = 1  # y->x
            else:
                # tie-break x>y : x->y iff x>y. Here x=S[i]<S[j]=y always (sorted), so y->x.
                if x > y:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1
    return S, adj

# --- tournament invariants ---
def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))

def num_3cycles(adj, n):
    c = 0
    for a in range(n):
        for b in range(n):
            for cc in range(n):
                if a<b<cc:
                    # count cyclic triangles among {a,b,cc}
                    pass
    # standard: number of cyclic triangles = C(n,3) - sum C(score_i,2)
    total = n*(n-1)*(n-2)//6
    s = sum((sum(adj[i][j] for j in range(n)))*(sum(adj[i][j] for j in range(n))-1)//2 for i in range(n))
    return total - s

def count_hampaths(adj, n):
    # count directed Hamiltonian paths
    cnt = 0
    for perm in permutations(range(n)):
        ok = True
        for k in range(n-1):
            if not adj[perm[k]][perm[k+1]]:
                ok = False; break
        if ok: cnt += 1
    return cnt

# canonical iso-class signature: brute force over all relabelings -> canonical adjacency bitstring
def canon_sig(adj, n):
    best = None
    for perm in permutations(range(n)):
        bits = 0
        idx = 0
        for i in range(n):
            for j in range(n):
                if i!=j:
                    pi, pj = perm[i], perm[j]
                    bit = adj[pi][pj]
                    bits = (bits<<1) | bit
        if best is None or bits < best:
            best = bits
    return best

def class_descriptor(adj, n):
    H = count_hampaths(adj, n)
    c3 = num_3cycles(adj, n)
    sc = score_seq(adj, n)
    sig = canon_sig(adj, n)
    return (H, c3, sc, sig)

# --- forbidden targets ---
# We identify "forbidden iso class" by (H, c3, score) since the claim is phrased that way,
# but we ALSO track the canonical sig to distinguish the two distinct (9,3,(1,1,2,3,3)) classes etc.

FORBIDDEN_N4 = {(5, 2, (1,1,2,2))}  # the regular near-3-cycle
FORBIDDEN_N5_BY_HCS = {
    (9, 3, (1,1,2,3,3)),
    (11,4,(1,2,2,2,3)),
    (13,4,(1,2,2,2,3)),
    (15,4,(1,2,2,2,3)),
    (15,5,(2,2,2,2,2)),
}

def is_forbidden(desc, n):
    H,c3,sc,sig = desc
    if n==4:
        return (H,c3,sc) in FORBIDDEN_N4
    if n==5:
        return (H,c3,sc) in FORBIDDEN_N5_BY_HCS
    return False

def primitive(S):
    from math import gcd as _g
    g = 0
    for v in S: g = _g(g, v)
    return g == 1

def search(n, speed_max, forbidden_check, label, sample=None, seed=0):
    """Exhaustive over all n-subsets of {1..speed_max} (or sampled). Report any realization of forbidden class."""
    found = []
    classes_seen = {}
    rng = random.Random(seed)
    universe = list(range(1, speed_max+1))
    if sample is None:
        iterator = combinations(universe, n)
        total = None
    else:
        def gen():
            seen=set()
            for _ in range(sample):
                S = tuple(sorted(rng.sample(universe, n)))
                if S in seen: continue
                seen.add(S)
                yield S
        iterator = gen()
    count = 0
    for S in iterator:
        if not primitive(S):
            continue
        count += 1
        Ssrt, adj = build_tournament_M3(S)
        desc = class_descriptor(adj, n)
        key = (desc[0], desc[1], desc[2])
        classes_seen.setdefault(key, 0)
        classes_seen[key]+=1
        if forbidden_check(desc, n):
            found.append((S, desc))
    return found, classes_seen, count

def main():
    print("UNITS (Z/14)* =", UNITS)
    print("="*70)

    # ---- Phase 0: reproduce the claim's baseline ----
    print("\n[Phase 0] Reproduce baseline exhaustive search as in claim")
    print("n=4, speeds in {1..15}:")
    f4, seen4, cnt4 = search(4, 15, is_forbidden, "n4")
    print(f"  primitive 4-subsets tested: {cnt4}")
    print(f"  distinct (H,c3,score) classes realized: {len(seen4)}")
    for k in sorted(seen4): print(f"    {k}: {seen4[k]}")
    print(f"  FORBIDDEN n=4 target {FORBIDDEN_N4} realized? {len(f4)>0}")
    for S,d in f4[:5]: print("    WITNESS:", S, d[:3])

    print("\nn=5, speeds in {1..13}:")
    f5, seen5, cnt5 = search(5, 13, is_forbidden, "n5")
    print(f"  primitive 5-subsets tested: {cnt5}")
    print(f"  distinct (H,c3,score) classes realized: {len(seen5)}")
    for k in sorted(seen5): print(f"    {k}: {seen5[k]}")
    realized_forb5 = set(k for k in seen5 if k in FORBIDDEN_N5_BY_HCS)
    print(f"  FORBIDDEN n=5 (H,c3,sc) realized in baseline: {realized_forb5}")
    for S,d in f5[:10]: print("    WITNESS:", S, d[:3])

    # ---- Phase 1: broaden n=4 search ----
    print("\n" + "="*70)
    print("[Phase 1] Broaden n=4: speeds {1..40} exhaustive")
    f4b, seen4b, cnt4b = search(4, 40, is_forbidden, "n4-wide")
    print(f"  primitive 4-subsets tested: {cnt4b}")
    print(f"  distinct (H,c3,score) classes: {len(seen4b)}")
    for k in sorted(seen4b): print(f"    {k}: {seen4b[k]}")
    print(f"  FORBIDDEN realized? {len(f4b)>0}")
    for S,d in f4b[:10]: print("    WITNESS:", S, d[:3])

    # ---- Phase 2: broaden n=5 search ----
    print("\n" + "="*70)
    print("[Phase 2] Broaden n=5: speeds {1..22} exhaustive")
    f5b, seen5b, cnt5b = search(5, 22, is_forbidden, "n5-wide")
    print(f"  primitive 5-subsets tested: {cnt5b}")
    print(f"  distinct (H,c3,score) classes: {len(seen5b)}")
    for k in sorted(seen5b): print(f"    {k}: {seen5b[k]}")
    realized_forb5b = set(k for k in seen5b if k in FORBIDDEN_N5_BY_HCS)
    print(f"  FORBIDDEN n=5 (H,c3,sc) realized: {realized_forb5b}")
    for S,d in f5b[:10]: print("    WITNESS:", S, d[:3])

    print("\nDONE phase 0-2.")

    # ---- Phase 3: CANONICAL EXHAUSTIVE over residue configs ----
    # M3 sign-votes depend only on speeds mod 14. The ONLY dependence on actual speed values is:
    #   (1) reference v0 = the slowest ACTUAL speed (so its residue is whatever the slowest speed's residue is),
    #   (2) tie-break "x>y" when D_xy==0 (which can only happen via frac ties; uses actual speed order).
    # Therefore EVERY achievable M3 tournament for ANY primitive LRC input arises from some
    #   (residue-tuple r_0<-(actual order), reference residue, total speed order) combination.
    # We enumerate ALL of these abstractly:
    #   - choose a multiset of n residues from {0..13} (with repetition) -> the n runners' residues
    #   - the runners are ordered by actual speed (we don't know the speeds, so we enumerate ALL
    #     consistent (residue, position) assignments). Reference = position 0 (slowest).
    #   - tie-break uses position order (actual speed). We enumerate all orderings by permuting which
    #     residue sits at which position.
    # This is a SUPERSET of all real inputs (some residue/order combos may be unrealizable by primitivity,
    # but if a forbidden class can't even appear here it certainly can't appear for real inputs).
    print("\n" + "="*70)
    print("[Phase 3] CANONICAL EXHAUSTIVE over residue configurations (covers ALL LRC inputs, any size)")
    run_residue_exhaustive(4, FORBIDDEN_N4, lambda d: (d[0],d[1],d[2]) in FORBIDDEN_N4)
    run_residue_exhaustive(5, FORBIDDEN_N5_BY_HCS, lambda d: (d[0],d[1],d[2]) in FORBIDDEN_N5_BY_HCS)

    print("\nDONE phase 3.")

def build_tournament_from_residues(res_by_pos):
    """res_by_pos: list of (residue) indexed by speed-position 0..n-1 (0=slowest).
       Reference = position 0. Tie-break: position with larger index = larger speed -> 'x>y' means larger position.
       Returns adj indexed by position."""
    n = len(res_by_pos)
    v0res = res_by_pos[0]
    w = {a: (1 if (v0res*a)%14 in (1,13) else -1) for a in UNITS}
    fr = {}
    for p in range(n):
        r = res_by_pos[p]
        fr[p] = {a: nrm(F((r%14)*a, 14)) for a in UNITS}
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            # 'x' = position i, 'y' = position j. But the spec's D_xy uses speeds sorted ascending;
            # position i<j means speed_i < speed_j (since positions ordered by actual speed).
            D = 0
            for a in UNITS:
                dx = fr[i][a]; dy = fr[j][a]
                s = 1 if dx>dy else (-1 if dx<dy else 0)
                D += w[a]*s
            if D > 0:
                adj[i][j] = 1
            elif D < 0:
                adj[j][i] = 1
            else:
                # tie-break x>y: x is speed_i, y is speed_j; speed_i<speed_j so x>y is FALSE -> y->x
                adj[j][i] = 1
    return adj

def run_residue_exhaustive(n, forbidden_set, check):
    from itertools import product as _product
    found = []
    classes = {}
    total = 0
    # All ordered residue tuples of length n over 0..13 (14^n). For n=4: 38416, n=5: 537824. Feasible.
    # Position 0 = slowest (reference). Order matters for tie-break, so use full product (ordered).
    for res_tuple in _product(range(14), repeat=n):
        total += 1
        adj = build_tournament_from_residues(list(res_tuple))
        desc = class_descriptor(adj, n)
        key = (desc[0], desc[1], desc[2])
        classes[key] = classes.get(key, 0) + 1
        if check(desc):
            found.append((res_tuple, desc))
    print(f"  n={n}: enumerated {total} ordered residue-tuples (14^{n})")
    print(f"  distinct (H,c3,score) classes realized: {len(classes)}")
    for k in sorted(classes): print(f"    {k}: {classes[k]}")
    print(f"  FORBIDDEN set targeted: {sorted(forbidden_set)}")
    print(f"  FORBIDDEN realized? {len(found)>0}")
    for res, d in found[:20]:
        print(f"    WITNESS residues={res} -> (H,c3,score)={(d[0],d[1],d[2])}")
    return found, classes

if __name__ == "__main__":
    main()
