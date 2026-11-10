# a01: (1) validate a_paths.exe against definition-level brute force (permutations, explicit path DFS);
#      (2) check PROOF.md section 2 table (Theorem B) AND THM-4879 (B) read through the claimed convention map
#          (seam index s = m - i for direct position i; SEAM's S = direct e; SEAM's (H,S,R) = direct (h,e,R))
#          against brute-force one-arc joins, all p, q (p = 1 included), random bases of order 1..7, m <= 8.
import random, sys
from a_lib import *

rng = random.Random(509)
# ---------------- (1) validation ----------------
tours = []
for n in range(0, 8):
    for _ in range(6 if n >= 2 else 1):
        tours.append(rand_tour(n, rng))
tours.append(TT(5)); tours.append(TT(1))
res = run(tours)
bad = 0
for A, r in zip(tours, res):
    nv = naive_stats(A)
    if (nv['H'], nv['Q'], nv['e'], nv['R'], nv['s']) != (r['H'], r['Q'], r['e'], r['R'], r['s']):
        bad += 1; print('VALIDATION MISMATCH', A, nv, r)
print('(1) a_paths.exe vs naive definitions: %d tournaments (orders 0-7), mismatches %d' % (len(tours), bad))

# ---------------- (2) Theorem B ----------------
def table_B(h, e, R, m, p, q):
    """PROOF.md section 2 (direct convention), positions 1-indexed"""
    hp = h + 2 ** (p - 1) * R - (1 if p == 1 else 0) * e
    if q == m and m > p: Euz = 2 ** (p - 1)
    elif q < p and p == m: Euz = 2 ** (q - 1)
    else: Euz = 0
    ep = (1 if q == m else 0) * h - (1 if (p == 1 and q == m) else 0) * e + Euz * R
    if p == 1 and q == 1: Rp = 2 ** (m - 1) * (h + R)
    elif p == 1: Rp = 2 ** (m - q) * hp
    elif q >= p: Rp = 2 ** (m - q) * hp
    elif q >= 2: Rp = 2 ** (m - q) * h + 3 ** (p - q - 1) * 2 ** (m - p + q) * R
    else: Rp = 2 ** (m - 1) * h - 2 ** (m - p) * e + 3 ** (p - 2) * 2 ** (m - p + 1) * R
    return hp, ep, Rp

def thm4879_B(H, S, R, k, u, z):
    """THM-4879 (B) verbatim, SEAM convention, vertices 0..k-1 of TT_k (0 the source), entry u, exit z."""
    Hrow = (1, -(1 if u == k - 1 else 0), 2 ** (k - 1 - u))
    if k == 1: Srow = (1, -1, 0)
    elif z == 0: Srow = (1, 0, 0) if u == 0 else Hrow
    else: Srow = (0, 0, 2 ** (k - 1 - z)) if u == 0 else (0, 0, 0)
    if z == k - 1 and k - 1 > u: mu = 2 ** u
    elif u == k - 1 and k - 1 > z: mu = 2 ** z
    else: mu = 0
    c = 2 ** (k - 1 - u + z) if z <= u else 2 ** (k - z + u) * 3 ** (z - u - 1)
    Rrow = (2 ** z, -mu, c)
    v = (H, S, R)
    dot = lambda row: sum(a * b for a, b in zip(row, v))
    return dot(Hrow), dot(Srow), dot(Rrow)

cases = 0; badB = 0; bad4879 = 0; cases_p1 = 0; cases_pge2 = 0
bases = []
for n in range(1, 8):
    for _ in range(5):
        A = rand_tour(n, rng); bases.append((A, rng.randrange(n)))
# add structured bases: transitive at sink / at source, 3-cycle
bases += [(TT(1), 0), (TT(2), 0), (TT(2), 1), (TT(4), 1)]
for (T, v) in bases:
    r0 = run([T])[0]; h, e, R = r0['H'], r0['e'][v], r0['R'][v]
    jobs = []; keys = []
    for m in range(1, 9):
        if len(T) + m > 15: continue
        for p in range(1, m + 1):
            N = one_arc_join(T, v, TT(m), p - 1)
            jobs.append(N); keys.append((m, p))
    outs = run(jobs)
    for (m, p), rr in zip(keys, outs):
        for q in range(1, m + 1):
            z = len(T) + q - 1
            got = (rr['H'], rr['e'][z], rr['R'][z])
            cases += 1
            if p == 1: cases_p1 += 1
            else: cases_pge2 += 1
            if table_B(h, e, R, m, p, q) != got:
                badB += 1
                if badB < 6: print('TABLE-B MISMATCH', (h, e, R), (m, p, q), got, table_B(h, e, R, m, p, q))
            # convention map: SEAM gadget = G^op; seam index of direct position i is m - i
            if thm4879_B(h, e, R, m, m - p, m - q) != got:
                bad4879 += 1
                if bad4879 < 6: print('THM-4879 MAP MISMATCH', (h, e, R), (m, p, q), got, thm4879_B(h, e, R, m, m - p, m - q))
print('(2) Theorem B: %d (base, m, p, q) cases from %d bases (orders 1-7, m <= 8, N <= 15): p = 1: %d, p >= 2: %d'
      % (cases, len(bases), cases_p1, cases_pge2))
print('    PROOF.md table mismatches: %d;  THM-4879 (B) through the map (s = m - i, S = e) mismatches: %d' % (badB, bad4879))
