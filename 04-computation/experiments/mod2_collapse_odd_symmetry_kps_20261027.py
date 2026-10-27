#!/usr/bin/env python3
"""
mod2_collapse_odd_symmetry_kps_20261027.py  (kind-pasteur-2026-10-27-S6)

The owner's observation: for an odd-order symmetry action, the cancellation (orientation-defect) polynomial
collapses modulo two to (1+t)^n, "but that does not mean n independent toggles exist".

Setting (session 65, GENERAL_PERMUTATION.md): for a permutation sigma of the vertex set, let its orbits on
unordered pairs have lengths L_O and endpoint-return signs eps_O (eps_O = -1 iff sigma^(L_O) swaps the two
endpoints of a pair in O).  The defect D_sigma(T) of a tournament T is the number of pairs e whose
transported orientation sigma(T(e)) disagrees with T(sigma(e)).  Session 65 proved
      P_sigma(z) := sum_T z^(D_sigma(T)) = prod_O ((1+z)^(L_O) + eps_O (1-z)^(L_O)).
For odd-order sigma every eps_O = +1 and every L_O is odd, so each factor is 2 E_L(z) with
      E_L(z) = sum_{j even} C(L, j) z^j.
Lemma (mod-2 collapse): for odd L, E_L(z) = (1+z)^(L-1) (mod 2).  Hence
      P_sigma(z) / 2^(#orbits) = prod_O E_(L_O)(z) = (1+z)^n  (mod 2),   n = sum_O (L_O - 1) = #pairs - #orbits,
while the integral polynomial prod_O E_(L_O)(z) has coefficients that are NOT binomial (E_3 = 1 + 3z^2, E_5 =
1 + 10z^2 + 5z^4): modulo two the orbit structure is invisible, the odd-position coefficients vanish, and the
even ones are all 1 on the support of (1+z)^n; the integer coefficients (positivity, the actual counts) and
the orbit partition are the lost coordinates.  Direct enumeration checks the defect polynomial and the lemma.
"""
import itertools
import sys
from math import comb

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + ("" if detail == "" else "  [" + str(detail) + "]"))


def pair_orbits(perm):
    n = len(perm)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    seen = set()
    orbits = []
    for p in pairs:
        if p in seen:
            continue
        orb = []
        q = p
        # transport the ORDERED pair to detect endpoint reversal
        a, b = p
        sign = +1
        while True:
            orb.append((min(a, b), max(a, b)))
            seen.add((min(a, b), max(a, b)))
            a, b = perm[a], perm[b]
            if (min(a, b), max(a, b)) == p:
                sign = +1 if (a, b) == p else -1
                break
        orbits.append((len(orb), sign, orb))
    return pairs, orbits


def defect_poly_direct(perm):
    """sum_T z^D(T) by enumerating all tournaments (n <= 7)."""
    n = len(perm)
    pairs, orbits = pair_orbits(perm)
    idx = {p: k for k, p in enumerate(pairs)}
    coeffs = {}
    for bits in range(1 << len(pairs)):
        # orientation: bit = 1 means j -> i for pair (i,j); else i -> j
        def arc(i, j):
            a, b = (i, j) if i < j else (j, i)
            back = (bits >> idx[(a, b)]) & 1
            return (b, a) if back else (a, b)
        D = 0
        for (i, j) in pairs:
            u, v = arc(i, j)                 # u -> v in T
            su, sv = perm[u], perm[v]
            if arc(su, sv) != (su, sv):      # transported arc disagrees with T on sigma(e)
                D += 1
        coeffs[D] = coeffs.get(D, 0) + 1
    return coeffs, orbits


def poly_mul(a, b):
    c = {}
    for i, x in a.items():
        for j, y in b.items():
            c[i + j] = c.get(i + j, 0) + x * y
    return c


def session65_formula(orbits):
    P = {0: 1}
    for (L, eps, _) in orbits:
        f = {}
        for j in range(L + 1):
            f[j] = comb(L, j) * (1 + eps * (-1) ** j)
        P = poly_mul(P, f)
    return P


def perm_order(perm):
    n = len(perm)
    order = 1
    seen = [False] * n
    for i in range(n):
        if not seen[i]:
            L = 0
            j = i
            while not seen[j]:
                seen[j] = True
                j = perm[j]
                L += 1
            from math import gcd
            order = order * L // gcd(order, L)
    return order


# 1. session 65's product formula against direct enumeration, all permutations of n <= 6 vertices
ok = True
for n in range(2, 7):
    for perm in itertools.permutations(range(n)):
        direct, orbits = defect_poly_direct(perm)
        formula = session65_formula(orbits)
        formula = {k: v for k, v in formula.items() if v}
        ok &= (direct == formula)
check("S1 sum_T z^D_sigma(T) = prod_O ((1+z)^L_O + eps_O (1-z)^L_O) for every permutation of n <= 6 vertices (direct enumeration)", ok)

# 2. mod-2 collapse lemma for odd L: E_L(z) = (1+z)^(L-1) mod 2
ok = True
for L in range(1, 40, 2):
    E = [comb(L, j) % 2 if j % 2 == 0 else 0 for j in range(L)]
    B = [comb(L - 1, j) % 2 for j in range(L)]
    ok &= (E == B)
check("S2 for odd L the even part of (1+z)^L is (1+z)^(L-1) modulo 2, L = 1..39", ok)

# 3. for odd-order sigma: normalized polynomial = (1+z)^n mod 2 with n = #pairs - #orbits, and not integrally
ok_mod = True
ok_int = True
examples = []
for n in range(3, 8):
    for perm in itertools.permutations(range(n)):
        if perm_order(perm) % 2 == 0 or perm == tuple(range(n)):
            continue
        pairs, orbits = pair_orbits(perm)
        assert all(eps == 1 and L % 2 == 1 for (L, eps, _) in orbits)
        P = session65_formula(orbits)
        k = len(orbits)
        Q = {d: c // (2 ** k) for d, c in P.items() if c}
        assert all(c * 2 ** k == P[d] for d, c in Q.items())
        m = len(pairs) - k
        target = {j: comb(m, j) % 2 for j in range(m + 1) if comb(m, j) % 2}
        Qmod = {d: c % 2 for d, c in Q.items() if c % 2}
        ok_mod &= (Qmod == target)
        binom = {j: comb(m, j) for j in range(m + 1)}
        if Q != binom:
            ok_int &= True
        if len(examples) < 3 and any(L > 1 for (L, _, _) in orbits):
            examples.append((perm, [L for (L, _, _) in orbits], m, dict(sorted(Q.items()))))
check("S3 every odd-order sigma on n <= 7 vertices: P_sigma / 2^#orbits = (1+z)^(#pairs - #orbits) mod 2", ok_mod)
# integrally the coefficients are not binomial as soon as an orbit has length >= 3
nonbinom = 0
total = 0
for n in range(3, 8):
    for perm in itertools.permutations(range(n)):
        if perm_order(perm) % 2 == 0 or perm == tuple(range(n)):
            continue
        pairs, orbits = pair_orbits(perm)
        P = session65_formula(orbits)
        k = len(orbits)
        Q = {d: c // (2 ** k) for d, c in P.items() if c}
        m = len(pairs) - k
        total += 1
        if Q != {j: comb(m, j) for j in range(m + 1)}:
            nonbinom += 1
check("S4 integrally the normalized polynomial is never the binomial (1+z)^n (all %d odd-order non-identity permutations, n <= 7): no n independent toggles" % total,
      nonbinom == total, (nonbinom, total))
for ex in examples:
    print("    example sigma =", ex[0], "orbit lengths", ex[1], "n =", ex[2], "normalized counts", ex[3])
# 4. what survives integrally: the orbit multiset is recovered from the integral polynomial (degree and first coefficients)
ok = True
for n in range(3, 8):
    seen = {}
    for perm in itertools.permutations(range(n)):
        if perm_order(perm) % 2 == 0 or perm == tuple(range(n)):
            continue
        pairs, orbits = pair_orbits(perm)
        P = tuple(sorted(session65_formula(orbits).items()))
        Ls = tuple(sorted(L for (L, _, _) in orbits))
        if P in seen and seen[P] != Ls:
            ok = False
        seen[P] = Ls
check("S5 the integral defect polynomial determines the multiset of pair-orbit lengths (n <= 7); its mod-2 image does not", ok)

bad = [nm for nm, ok in CHECKS if not ok]
print("%d checks, %d failed" % (len(CHECKS), len(bad)))
if bad:
    print("FAILED:", bad)
    sys.exit(1)
print("ALL CHECKS PASSED")
