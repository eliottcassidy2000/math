#!/usr/bin/env python3
"""Recon for HYP-7260 (death-star-S51): relation lock by coefficient weight.

(a) RELATION LOCK: sum alpha_i v_i = 0, sum |alpha_i| <= 14, all fail
    => sum alpha_i w_i = 0.  Boundary: weight 14 locks (strict); breaks at 15.
(b) WEIGHT-14 PAIRS of {1..13} -- (1,13),(3,11),(5,9): branches EMPTY
    (correcting the S48 locked/sparse boundary: lock <=> weight <= 14).
(c) SUM-TRIPLE LOCK: {a,b,a+b}: w_c = w_a + w_b always (weight 3) --
    including triples whose PAIRS are sparse, e.g. (5,6,11).
(d) MEDIANT TRIPLE COUNT: c = a+b <= 13, g = gcd structure: N_triple =
    2*floor((q-1)/(14*(a+b)/gcd(a,b,c)...)) -- pin the exact form.
"""
from math import gcd
import random
random.seed(51)

def fw(v, q, p):
    r = (v * p) % q
    if 14 * r < q: return (v * p) // q
    if 14 * r > 13 * q: return (v * p) // q + 1
    return None

# (a) random relations: pick speeds, coefficients with sum alpha v = 0
lock_ok = lock_bad = brk15 = 0
for _ in range(600000):
    # construct a vanishing relation: alpha*va + beta*vb = gamma*vc
    al = random.randint(1, 6); be = random.randint(1, 6); 
    va = random.randint(1, 40); vb = random.randint(1, 40)
    # choose gamma, vc with gamma*vc = al*va + be*vb
    tot = al * va + be * vb
    ga = random.choice([d for d in range(1, 8) if tot % d == 0])
    vc = tot // ga
    w = al + be + ga
    q = random.randint(2, 600); p = random.randint(1, q - 1)
    wa, wb, wc = fw(va, q, p), fw(vb, q, p), fw(vc, q, p)
    if wa is None or wb is None or wc is None: continue
    k = al * wa + be * wb - ga * wc
    if w <= 14:
        if k == 0: lock_ok += 1
        else:
            lock_bad += 1
            if lock_bad <= 3: print("RELLOCK FAIL", al, be, ga, va, vb, vc, q, p, k)
    elif w == 15 and k != 0:
        brk15 += 1
print(f"(a) relation lock (weight<=14): {lock_ok} ok, {lock_bad} bad",
      "PASS" if lock_bad == 0 else "FAIL", f"| weight-15 breaks seen: {brk15}")

# (b) weight-14 pairs: branches empty
for (a, b) in [(1, 13), (3, 11), (5, 9)]:
    tot_branch = 0
    for q in range(300, 1000, 97):
        for p in range(1, q):
            wa, wb = fw(a, q, p), fw(b, q, p)
            if wa is None or wb is None: continue
            if b * wa - a * wb != 0: tot_branch += 1
    print(f"(b) pair ({a},{b}) weight 14: nonzero-branch occurrences = {tot_branch}",
          "PASS(empty)" if tot_branch == 0 else "FAIL")

# (c) sum-triple lock incl. sparse-pair triples
tri_ok = tri_bad = 0
triples = [(a, b, a + b) for a in range(1, 13) for b in range(a + 1, 14 - a)]
for (a, b, c) in triples:
    for q in random.sample(range(200, 2000), 5):
        for p in range(1, q):
            wa, wb, wc = fw(a, q, p), fw(b, q, p), fw(c, q, p)
            if wa is None or wb is None or wc is None: continue
            if wa + wb == wc: tri_ok += 1
            else:
                tri_bad += 1
                if tri_bad <= 3: print("TRIPLE FAIL", a, b, c, q, p)
print(f"(c) sum-triple lock over {len(triples)} triples: {tri_ok} ok, {tri_bad} bad",
      "PASS" if tri_bad == 0 else "FAIL")

# (d) mediant triple count: N = 2*floor((q-1)/(14*L)) -- pin L
print("(d) mediant triple counts (find L such that N = 2*floor((q-1)/(14*L))):")
for (a, b, c) in [(1, 2, 3), (2, 3, 5), (3, 4, 7), (2, 4, 6), (4, 6, 10), (5, 6, 11), (3, 6, 9), (2, 5, 7), (6, 7, 13)]:
    g = gcd(a, gcd(b, c))
    hits = []
    for q in [911, 1301, 2003]:
        N = sum(1 for p in range(1, q)
                if fw(a, q, p) is not None and fw(b, q, p) is not None
                and fw(c, q, p) is not None)
        # candidate L = c/g (top member over global gcd)
        pred = 2 * ((q - 1) // (14 * (c // g)))
        hits.append((N, pred))
    ok = all(n == pr for n, pr in hits)
    print(f"    ({a},{b},{c}) g={g}: L = c/g = {c//g}: {hits}",
          "PASS" if ok else "CHECK")
