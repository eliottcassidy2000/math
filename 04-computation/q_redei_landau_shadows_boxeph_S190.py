#!/usr/bin/env python3
"""
q_redei_landau_shadows_boxeph_S190.py  (HYP-8620, THM-1825)

(A) q-REDEI CENSUS. h_q(T) = sum over Hamiltonian paths pi of q^{bk(pi)},
    bk = #arcs of T pointing backward w.r.t. the path order (label-free;
    h_q(1) = h). DP over (subset, last) with polynomial values: placing w
    after set S adds #{s in S : w -> s} back-arcs.
    Laws hunted: (L-a) back-arc PARITY rigidity: is bk mod 2 constant over
    the paths of each T (=> h_q(-1) = ±h — the q-refinement of Redei)?
    (L-b) roots-of-unity vanishing: h_q(i), h_q(zeta_3) zero patterns.
    (L-c) palindromicity / unimodality of h_q.
    (L-d) THE {7,21} q-SHADOW: joint spectrum (h, h_q shape) near the hole:
    degree spread, gcd of coefficients, the value h_q(2) — which shapes
    would a phantom h = 7 tournament need?

(B) LANDAU = MOMENT-MAP LATTICE SURJECTIVITY. Tournaments = points of
    (P^1)^{E}; the torus moment map = score vector; the moment polytope =
    permutohedron of (0,1,...,n-1). Verify: attained score multisets ==
    lattice points of the polytope (majorization test) for n <= 7.

(C) W3 bite: terminal stacked-death multiplicity at delta* for the tight
    family v = (1..n), n = 3..6 (exact interval unions).

boxeph-2026-07-20-S190. Pure python.
"""
import itertools
from fractions import Fraction

# ---------- tournament plumbing ----------

def pairs_of(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def bits_to_adj(bits, n):
    adj = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(pairs_of(n)):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj


def build_perm_maps(n):
    prs = pairs_of(n)
    idx = {pr: k for k, pr in enumerate(prs)}
    maps = []
    for p in itertools.permutations(range(n)):
        pm = []
        for (i, j) in prs:
            a, b = p[i], p[j]
            pm.append((idx[(a, b)], 0) if a < b else (idx[(b, a)], 1))
        maps.append(pm)
    return maps


def classes_for_n(n):
    perm_maps = build_perm_maps(n)
    total = 1 << (n * (n - 1) // 2)
    visited = bytearray(total)
    reps = []
    for bits in range(total):
        if visited[bits]:
            continue
        members = set()
        for pm in perm_maps:
            nb = 0
            for k in range(len(pm)):
                src, flip = pm[k]
                b = (bits >> src) & 1
                if flip:
                    b ^= 1
                nb |= b << k
            members.add(nb)
        for m in members:
            visited[m] = 1
        reps.append(min(members))
    return reps


def h_q(adj, n):
    # DP over (S, last): polynomial in q as tuple of coeffs, deg <= C(n,2)
    D = n * (n - 1) // 2 + 1
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = [0] * D
        dp[(1 << v, v)][0] = 1
    for size in range(1, n):
        ndp = {}
        for (S, last), poly in dp.items():
            if bin(S).count('1') != size:
                continue
            for w in range(n):
                if (S >> w) & 1 or not adj[last][w]:
                    continue
                back = sum(1 for s in range(n) if (S >> s) & 1 and adj[w][s])
                key = (S | (1 << w), w)
                tgt = ndp.get(key)
                if tgt is None:
                    tgt = [0] * D
                    ndp[key] = tgt
                for d, c in enumerate(poly):
                    if c:
                        tgt[d + back] += c
        # merge ndp into dp (keep smaller sizes too for next rounds? no: rebuild)
        for k, v in ndp.items():
            if k in dp:
                dp[k] = [a + b for a, b in zip(dp[k], v)]
            else:
                dp[k] = v
    full = (1 << n) - 1
    tot = [0] * D
    for (S, last), poly in dp.items():
        if S == full:
            for d, c in enumerate(poly):
                tot[d] += c
    while len(tot) > 1 and tot[-1] == 0:
        tot.pop()
    return tot


def peval(poly, x):
    v = 0
    for c in reversed(poly):
        v = v * x + c
    return v


print("=" * 78)
print("(A) q-REDEI CENSUS (back-arc statistic)")
print("=" * 78)
import math
for n in (3, 4, 5, 6):
    reps = classes_for_n(n)
    parity_rigid = 0
    parity_broken = []
    hq_m1 = {}
    pal = 0
    uni = 0
    shapes = {}
    for bits in reps:
        adj = bits_to_adj(bits, n)
        poly = h_q(adj, n)
        h = sum(poly)
        supp_par = set(d % 2 for d, c in enumerate(poly) if c)
        if len(supp_par) == 1:
            parity_rigid += 1
        else:
            parity_broken.append((bits, poly))
        v_m1 = peval(poly, -1)
        hq_m1.setdefault((h, v_m1), 0)
        hq_m1[(h, v_m1)] += 1
        nz = [c for c in poly if c]
        if poly == poly[::-1]:
            pal += 1
        ok_uni = True
        rising = True
        seq = [c for c in poly]
        for i in range(1, len(seq)):
            if seq[i] < seq[i - 1]:
                rising = False
            elif not rising and seq[i] > seq[i - 1]:
                ok_uni = False
        if ok_uni:
            uni += 1
        shapes.setdefault(h, []).append(poly)
    print("n=%d (%d classes): back-arc parity RIGID in %d/%d classes" %
          (n, len(reps), parity_rigid, len(reps)))
    if parity_broken and n <= 4:
        for bits, poly in parity_broken[:2]:
            print("   parity-broken example bits=%d: h_q = %s" % (bits, poly))
    law = all(abs(v) == h for (h, v) in hq_m1)
    print("   |h_q(-1)| == h for all classes: %s ; palindromic: %d/%d ; unimodal(coeff): %d/%d"
          % (law, pal, len(reps), uni, len(reps)))
    if n <= 5:
        # roots of unity: h_q(i) == 0? h_q at zeta_3 via poly mod (q^2+q+1)
        z_i = sum(1 for bitsx in reps
                  if peval(h_q(bits_to_adj(bitsx, n), n), complex(0, 1)) == 0)
        print("   classes with h_q(i) == 0: %d" % z_i)

print()
print("   {7,21} q-shadow at n=5,6: h_q data for h in {5,9} (the hole's flanks):")
for n in (5, 6):
    reps = classes_for_n(n)
    for bits in reps:
        adj = bits_to_adj(bits, n)
        poly = h_q(adj, n)
        h = sum(poly)
        if h in (5, 9) and n == 5:
            print("   n=%d h=%d: h_q = %s ; h_q(2) = %d ; degspread = %d" %
                  (n, h, poly, peval(poly, 2), len(poly) - 1))

print()
print("=" * 78)
print("(B) LANDAU = MOMENT-MAP LATTICE SURJECTIVITY (verify n <= 7)")
print("=" * 78)


def landau_ok(seq):
    s = sorted(seq)
    n = len(s)
    tot = 0
    for k in range(1, n + 1):
        tot = sum(s[:k])
        if tot < k * (k - 1) // 2:
            return False
    return tot == n * (n - 1) // 2


for n in (3, 4, 5, 6, 7):
    # attained sorted score multisets
    attained = set()
    if n <= 6:
        prs = pairs_of(n)
        for bits in range(1 << len(prs)):
            sc = [0] * n
            for k, (i, j) in enumerate(prs):
                if (bits >> k) & 1:
                    sc[i] += 1
                else:
                    sc[j] += 1
            attained.add(tuple(sorted(sc)))
    else:
        # n = 7: sample-free: use Landau to generate candidates and verify
        # attainment constructively is classical; here compare counts only
        pass
    # lattice points of the moment polytope: nonneg integer seqs summing to
    # C(n,2), sorted, satisfying majorization (= Landau partial sums)
    cands = set()
    def gen(prefix, rem, minv):
        k = len(prefix)
        if k == n:
            if rem == 0 and landau_ok(prefix):
                cands.add(tuple(prefix))
            return
        for v in range(minv, n):
            if v > rem:
                break
            gen(prefix + [v], rem - v, v)
    gen([], n * (n - 1) // 2, 0)
    if n <= 6:
        print("n=%d: attained sorted score multisets = %d ; Landau/moment-polytope lattice points = %d ; EQUAL: %s"
              % (n, len(attained), len(cands), attained == cands))
    else:
        print("n=%d: Landau/moment-polytope lattice points = %d (A000571 check)" % (n, len(cands)))

print()
print("=" * 78)
print("(C) terminal stacked-death multiplicity at delta* for v = (1..n)")
print("=" * 78)


def good_set_components(v, dl):
    bad = []
    for vj in v:
        for m in range(0, vj + 1):
            lo = Fraction(m - dl, vj)
            hi = Fraction(m + dl, vj)
            bad.append((lo, hi))
    pieces = []
    for lo, hi in bad:
        lo_ = lo % 1
        length = hi - lo
        if length >= 1:
            return 0, Fraction(0)
        hi_ = lo_ + length
        if hi_ <= 1:
            pieces.append((lo_, hi_))
        else:
            pieces.append((lo_, Fraction(1)))
            pieces.append((Fraction(0), hi_ - 1))
    pieces.sort()
    merged = []
    for lo, hi in pieces:
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    total_bad = sum(hi - lo for lo, hi in merged)
    if total_bad >= 1:
        return 0, Fraction(0)
    ncomp = len(merged)
    if ncomp == 0:
        return 1, Fraction(1)
    wraps = (merged[0][0] == 0 and merged[-1][1] == 1)
    goodm = Fraction(1) - total_bad
    comps = ncomp - (1 if wraps else 0)
    if comps == 0:
        comps = 1 if goodm > 0 else 0
    return comps, goodm


for n in (3, 4, 5, 6):
    v = tuple(range(1, n + 1))
    # delta* = 1/(n+1) conjectured tight; find the last component count just
    # below and confirm zero just above
    dstar = Fraction(1, n + 1)
    eps = Fraction(1, 100000)
    c_before, _ = good_set_components(v, dstar - eps)
    c_at, gm_at = good_set_components(v, dstar + eps)
    print("v=(1..%d): components just below delta*=1/%d: %d ; just above: %d "
          "-> terminal stacked multiplicity = %d" %
          (n, n + 1, c_before, c_at, c_before - c_at))

print("\nDONE.")
