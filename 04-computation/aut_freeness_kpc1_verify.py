#!/usr/bin/env python3
"""
ADVERSARIAL VERIFICATION (kind-pasteur-2026-06-10-S1, thread B-universal-freeness).

Independent re-verification of claims B1/B2/B3 about LEM-003:
  Aut(D) acts freely on directed Hamiltonian paths, hence |Aut(T)| | H(T).

FRESH METHODS (deliberately different from worker's Held-Karp + per-mask S_n loop):
  * H(T) counted by brute force over all n! vertex orderings, vectorized over
    ALL 2^C(n,2) labeled tournaments simultaneously with numpy:
    ordering p is a Ham path of mask m  iff  (m & B_p) == V_p,
    where B_p = bits of the 5 consecutive pairs, V_p = required orientations.
  * |Aut| via the relabeling action p.m computed by chunked lookup tables
    (low byte / high bits), vectorized over all masks; p in Aut(T_m) iff p.m == m.
  * Freeness: for EVERY mask with |Aut|>1 (n=5 and n=6), the action of Aut on
    Ham-path vertex sequences is constructed explicitly; assert every
    non-identity automorphism fixes NO path, every orbit has size |Aut|,
    and #orbits * |Aut| == H.
  * Cycle caveat: enumerate Ham CYCLE arc sets directly for C3 and the
    circulant RQ5 (i -> i+1, i+2 mod 5); exhibit rotation-fixed cycles.

Bit convention: pairs (i,j), i<j in lexicographic order, bit set <=> arc i->j.
All arithmetic exact (numpy int64 well within range; final tallies as Python ints).
"""
import itertools, sys
import numpy as np

def pair_index(n):
    idx = {}
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            idx[(i, j)] = k
            k += 1
    return idx, k

def arc_present(m, i, j, idx):
    """True iff arc i->j present in tournament mask m."""
    if i < j:
        return (m >> idx[(i, j)]) & 1 == 1
    return (m >> idx[(j, i)]) & 1 == 0

def act_mask(p, m, n, idx):
    """Relabeled tournament p.m : arc p(i)->p(j) iff arc i->j in m."""
    out = 0
    for (i, j), k in idx.items():
        # orientation of pair (i,j) in p.m comes from preimage pair
        a, b = p.index(i), p.index(j)   # p(a)=i, p(b)=j
        if arc_present(m, a, b, idx):
            out |= 1 << k               # i->j in p.m
    return out

def analyze_n(n, expected_classes):
    idx, C = pair_index(n)
    NM = 1 << C
    masks = np.arange(NM, dtype=np.int64)
    perms = list(itertools.permutations(range(n)))
    fact = len(perms)

    # ---- H by brute force over orderings, vectorized over all masks ----
    H = np.zeros(NM, dtype=np.int64)
    for p in perms:
        B = 0; V = 0
        for t in range(n - 1):
            a, b = p[t], p[t + 1]
            if a < b:
                k = idx[(a, b)]; B |= 1 << k; V |= 1 << k
            else:
                k = idx[(b, a)]; B |= 1 << k
        H += ((masks & B) == V)

    # ---- |Aut| via chunked lookup tables, vectorized over all masks ----
    LOW = 8 if C > 8 else C
    nlow = 1 << LOW
    nhigh = 1 << (C - LOW)
    AUT = np.zeros(NM, dtype=np.int64)
    aut_elements = {}   # filled later per-mask only where needed
    for p in perms:
        # output bit k (pair (i,j)) takes value of arc p^{-1}(i)->p^{-1}(j) in m
        pinv = [0] * n
        for a, v in enumerate(p):
            pinv[v] = a
        # source spec per output bit: (source_bit, negate?)
        src = []
        for (i, j), k in idx.items():
            a, b = pinv[i], pinv[j]
            if a < b:
                src.append((k, idx[(a, b)], False))
            else:
                src.append((k, idx[(b, a)], True))
        T1 = np.zeros(nlow, dtype=np.int64)
        T2 = np.zeros(nhigh, dtype=np.int64)
        for c in range(nlow):
            v = 0
            for k, s, neg in src:
                if s < LOW:
                    bit = (c >> s) & 1
                    if neg: bit ^= 1
                    v |= bit << k
            T1[c] = v
        for c in range(nhigh):
            v = 0
            for k, s, neg in src:
                if s >= LOW:
                    bit = (c >> (s - LOW)) & 1
                    if neg: bit ^= 1
                    v |= bit << k
            T2[c] = v
        pm = T1[masks & (nlow - 1)] | T2[masks >> LOW]
        AUT += (pm == masks)

    sumaut = int(AUT.sum())
    print(f"=== n={n}: {NM} labeled tournaments ===")
    print(f"  H=0 count: {int((H == 0).sum())}   (Redei: expect 0)")
    print(f"  all H odd: {bool(((H & 1) == 1).all())}   (Redei: expect True)")
    print(f"  max H: {int(H.max())}")
    div_fail = int((H % AUT != 0).sum())
    print(f"  divisibility |Aut| | H failures: {div_fail}")
    print(f"  sum |Aut| over all masks: {sumaut}  "
          f"(Burnside check: {fact}*{expected_classes} = {fact*expected_classes}, "
          f"match={sumaut == fact*expected_classes})")
    aut_vals = sorted(set(int(a) for a in np.unique(AUT)))
    print(f"  |Aut| value set: {aut_vals}")
    nontriv = [int(m) for m in np.nonzero(AUT > 1)[0]]
    print(f"  masks with |Aut|>1: {len(nontriv)}")

    # ---- FREENESS: explicit action on path sequences, every nontrivial mask ----
    free_fail = 0
    profiles = {}
    for m in nontriv:
        auts = [p for p in perms if act_mask(p, m, n, idx) == m]
        paths = []
        for p in perms:
            if all(arc_present(m, p[t], p[t + 1], idx) for t in range(n - 1)):
                paths.append(p)
        a = len(auts)
        h = len(paths)
        assert a == int(AUT[m]) and h == int(H[m]), "internal cross-check failed"
        pset = set(paths)
        # stabilizers must be trivial: no non-identity g fixes any path
        for g in auts:
            if g == tuple(range(n)):
                continue
            for P in paths:
                gp = tuple(g[v] for v in P)
                assert gp in pset, "automorphism image not a path?!"
                if gp == P:
                    free_fail += 1
        # orbit sizes
        seen = set()
        orbits = 0
        for P in paths:
            if P in seen:
                continue
            orb = set(tuple(g[v] for v in P) for g in auts)
            if len(orb) != a:
                free_fail += 1
            seen |= orb
            orbits += 1
        if orbits * a != h:
            free_fail += 1
        profiles.setdefault((a, h, orbits), 0)
        profiles[(a, h, orbits)] += 1
    print(f"  freeness failures over all {len(nontriv)} nontrivial-Aut masks: {free_fail}")
    print(f"  (|Aut|, H, #orbits) profile (nontrivial masks): "
          f"{sorted(profiles.items())}")
    return div_fail == 0 and free_fail == 0

def ham_cycles_arcsets(n, mask, idx):
    """All directed Ham cycles as frozensets of arcs (i,j)."""
    out = set()
    for p in itertools.permutations(range(1, n)):
        seq = (0,) + p
        arcs = [(seq[t], seq[(t + 1) % n]) for t in range(n)]
        if all(arc_present(mask, a, b, idx) for a, b in arcs):
            out.add(frozenset(arcs))
    return sorted(out, key=sorted)

def cycle_caveat():
    print("=== CYCLE CAVEAT ===")
    # C3: 0->1, 1->2, 2->0
    n = 3
    idx, C = pair_index(n)
    m = 0
    for (i, j) in [(0, 1), (1, 2), (2, 0)]:
        if i < j:
            m |= 1 << idx[(i, j)]
    perms = list(itertools.permutations(range(n)))
    auts = [p for p in perms if act_mask(p, m, n, idx) == m]
    cycles = ham_cycles_arcsets(n, m, idx)
    print(f"  C3: |Aut|={len(auts)}, #HamCycles={len(cycles)}, "
          f"3 divides 1: {len(cycles) % len(auts) == 0}")
    rot = (1, 2, 0)
    cyc = cycles[0]
    rcyc = frozenset((rot[a], rot[b]) for a, b in cyc)
    print(f"  C3: rotation (0 1 2) fixes the unique Ham cycle: {rcyc == cyc}")

    # RQ5 circulant: i -> i+1, i+2 (mod 5)
    n = 5
    idx, C = pair_index(n)
    m = 0
    for i in range(5):
        for d in (1, 2):
            j = (i + d) % 5
            if i < j:
                m |= 1 << idx[(i, j)]
    perms = list(itertools.permutations(range(n)))
    auts = [p for p in perms if act_mask(p, m, n, idx) == m]
    cycles = ham_cycles_arcsets(n, m, idx)
    print(f"  RQ5: |Aut|={len(auts)}, aut group = rotations only: "
          f"{sorted(auts) == sorted(tuple((v + s) % 5 for v in range(5)) for s in range(5))}")
    print(f"  RQ5: #HamCycles={len(cycles)} (claim: 2), "
          f"5 divides {len(cycles)}: {len(cycles) % 5 == 0}")
    rot = tuple((v + 1) % 5 for v in range(5))
    orbit_sizes = []
    seen = set()
    for c in cycles:
        if c in seen:
            continue
        orb = set()
        cc = c
        for g in auts:
            orb.add(frozenset((g[a], g[b]) for a, b in c))
        seen |= orb
        orbit_sizes.append(len(orb))
    print(f"  RQ5: cycle-orbit sizes under Aut: {sorted(orbit_sizes)} (claim: [1, 1])")
    fixed_by_rot = sum(1 for c in cycles
                       if frozenset((rot[a], rot[b]) for a, b in c) == c)
    print(f"  RQ5: cycles fixed by rotation +1: {fixed_by_rot} of {len(cycles)}")
    # identify the two cycles by step pattern
    for c in cycles:
        steps = sorted(((b - a) % 5) for a, b in c)
        print(f"    cycle arc-steps (mod 5): {steps}")
    # path action on the same tournament IS free
    paths = [p for p in perms
             if all(arc_present(m, p[t], p[t + 1], idx) for t in range(4))]
    h = len(paths)
    pset = set(paths)
    stab_viol = 0
    for g in auts:
        if g == tuple(range(5)):
            continue
        for P in paths:
            if tuple(g[v] for v in P) == P:
                stab_viol += 1
    seen = set(); orbs = []
    for P in paths:
        if P in seen:
            continue
        orb = set(tuple(g[v] for v in P) for g in auts)
        seen |= orb
        orbs.append(len(orb))
    print(f"  RQ5 paths: H={h} (claim 15), stabilizer violations={stab_viol}, "
          f"path-orbit sizes={sorted(orbs)} (claim [5,5,5])")

def main():
    ok5 = analyze_n(5, 12)   # A000568(5) = 12
    ok6 = analyze_n(6, 56)   # A000568(6) = 56
    cycle_caveat()
    print(f"\nVERDICT: n=5 ok={ok5}, n=6 ok={ok6}")

if __name__ == "__main__":
    main()
