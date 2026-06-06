"""SIGNED LRC sign-orbit collisions = CYCLIC HOMOMETRY of half-system signed sets.
monad-explorer-2026-06-06-S705.

REFRAME (the new lens).  AP_n runners = {1,...,n-1} = the half-system {1,...,(C-1)/2} of Z/C,
C=2n-1.  A sign vector eps maps runner i -> point u_i = eps_i*i in Z/C.  Since
{0} U {+-1,...,+-(n-1)} = Z/C, the point set S_eps = {eps_i*i} is an (n-1)-subset of Z/C\{0}
picking exactly one of {i, C-i} per magnitude (a "half-system selection").  The folded pair
spectrum is EXACTLY the multiset of circular distances rho(u_i-u_j).  Hence:

   two cuts collide  <=>  S_eps and S_eps' are HOMOMETRIC (same difference multiset in Z/C)
                     <=>  same |DFT|^2 (same Patterson / autocorrelation).

This script:
  (1) VERIFIES the three notions coincide (folded-spectrum = diff-multiset = |DFT|^2) exhaustively.
  (2) Tabulates deficiency = 2^{n-2} - folded_orbit and the class-size distribution.
  (3) Tests whether collision classes are COSETS of a group H of universally-silent flips.
  (4) Prints |DFT|^2 null structure (which frequencies t can ever be killed) vs factorization of C.
"""
from itertools import product
import cmath, math, sys
from collections import Counter


def folded_spectrum(V, eps, C):
    out = []
    n1 = len(V)
    for i in range(n1):
        for j in range(i + 1, n1):
            f = (eps[i] * V[i] - eps[j] * V[j]) % C
            out.append(min(f, C - f))
    return tuple(sorted(out))


def diff_multiset(pts, C):
    """full ordered difference multiset in Z/C (autocorrelation support)."""
    out = []
    for a in pts:
        for b in pts:
            if a != b:
                out.append((a - b) % C)
    return tuple(sorted(out))


def dft_power(pts, C):
    """|chi_hat(t)|^2 for t=0..C-1, rounded -- the Patterson spectrum."""
    out = []
    for t in range(C):
        s = sum(cmath.exp(2j * math.pi * t * p / C) for p in pts)
        out.append(round(abs(s) ** 2, 6))
    return tuple(out)


def analyze(V, triple=False):
    n = len(V) + 1
    C = 2 * n - 1
    cuts = []
    for bits in product([0, 1], repeat=len(V) - 1):
        eps = [1] + [1 if b else -1 for b in bits]  # fix runner 0 sign=+ (global swap quotient)
        cuts.append(eps)
    fold_map, diff_map, dft_map = {}, {}, {}
    for eps in cuts:
        fs = folded_spectrum(V, eps, C)
        fold_map.setdefault(fs, []).append(tuple(eps))
        if triple:
            pts = [(eps[i] * V[i]) % C for i in range(len(V))]
            diff_map.setdefault(diff_multiset(pts, C), []).append(tuple(eps))
            dft_map.setdefault(dft_power(pts, C), []).append(tuple(eps))
    triple_consistent = None
    if triple:
        def partition_of(m):
            return frozenset(frozenset(g) for g in m.values())
        triple_consistent = (partition_of(fold_map) == partition_of(diff_map) == partition_of(dft_map))
    sizes = sorted((len(g) for g in fold_map.values()), reverse=True)
    deficiency = len(cuts) - len(fold_map)
    return dict(n=n, C=C, total=len(cuts), orbit=len(fold_map), deficiency=deficiency,
                sizes=Counter(sizes), consistent=triple_consistent, fold_map=fold_map, cuts=cuts)


def is_prime(m):
    return m > 1 and all(m % d for d in range(2, int(m**0.5) + 1))


def group_test(V, r):
    """Is the collision-equivalence a union of cosets of the universally-silent subgroup H?
    H = {T subset of runners : flipping T preserves the AP folded spectrum for ALL cuts}.
    A cleaner direct test: take the equivalence class of the all-+ cut; is it a subgroup of
    (Z/2)^{n-1} (closed under symmetric difference of the sign-flip patterns vs all-+)?"""
    C = r['C']
    n1 = len(V)
    # baseline cut = all +
    base = tuple([1] * n1)
    base_sig = folded_spectrum(V, list(base), C)
    # silent flips of the base = sign patterns giving same spectrum, as flip-sets vs base
    flipsets = []
    for eps, in [(e,) for e in r['cuts']]:
        if folded_spectrum(V, list(eps), C) == base_sig:
            T = frozenset(i for i in range(n1) if eps[i] != base[i])
            # canonicalize under global swap (T or complement)
            comp = frozenset(range(n1)) - T
            flipsets.append(min(T, comp, key=lambda s: (len(s), sorted(s))))
    flipsets = set(flipsets)
    # closure under symmetric difference?
    is_group = True
    fl = list(flipsets)
    for a in fl:
        for b in fl:
            sd = a ^ b
            comp = frozenset(range(n1)) - sd
            canon = min(sd, comp, key=lambda s: (len(s), sorted(s)))
            if canon not in flipsets:
                is_group = False
                break
        if not is_group:
            break
    return len(flipsets), is_group


print("=" * 96)
print("PART 1: collision = homometry = |DFT|^2.  Verify three notions coincide; tabulate deficiency.")
print("=" * 96)
print(f"{'n':>3} {'C':>4} {'C-fact':>10} {'2^(n-2)':>8} {'orbit':>6} {'defic':>6} "
      f"{'consistent':>10} {'class-sizes':>22} {'silentH':>8} {'grp':>4}", flush=True)
for n in range(3, 21):
    V = list(range(1, n))
    do_triple = (n <= 13)   # DFT/diff verification only where cheap (C<=25)
    r = analyze(V, triple=do_triple)
    C = r['C']
    fact = []
    m = C
    d = 2
    while d * d <= m:
        while m % d == 0:
            fact.append(d)
            m //= d
        d += 1
    if m > 1:
        fact.append(m)
    factstr = "x".join(map(str, fact)) if not is_prime(C) else f"{C}(P)"
    hsize, isgrp = group_test(V, r)
    szs = dict(r['sizes'])
    print(f"{n:>3} {C:>4} {factstr:>10} {r['total']:>8} {r['orbit']:>6} {r['deficiency']:>6} "
          f"{str(r['consistent']):>10} {str(szs):>22} {hsize:>8} {str(isgrp):>4}", flush=True)
