"""SIGNED-LRC sign-orbit: CLOSING the converse of HYP-2270 + uniform count law.
monad-explorer-2026-06-06-S707.

Setup (S705 homometry lens). C=2n-1, runners V={1,..,n-1}=half-system of Z/C.
A cut eps in {+-1}^{n-1} maps runner i -> eps_i*i.  Folded clock multiset = multiset of
circular distances rho(eps_i*i - eps_j*j).  Decompose |f_eps(zeta^t)|^2 = A(t)^2 + Phi_eps(t)^2,
A(t)=sum_i cos(2pi t i/C) cut-INDEP, Phi_eps(t)=sum_i eps_i sin(2pi t i/C) cut-DEP (DST).
COLLISION  eps ~ eps'  <=>  Phi_eps(t)^2 = Phi_eps'(t)^2  for all t  (per-freq sign).

Goal 1 (CONVERSE of HYP-2270): C composite => orbit < 2^{n-2} (a collision EXISTS).
  We exhibit an explicit pair and identify the exact "silent condition" on eps.
Goal 2 (COUNT LAW): tabulate deficiency vs factorization, find uniform formula beyond 3p.

This script: brute exact for small n; extracts the silent-condition; tests a constructive
collision pair for every composite C (no brute force) to support the converse proof.
"""
import math
from itertools import product
from collections import Counter, defaultdict


def folded_spectrum(V, eps, C):
    out = []
    n1 = len(V)
    for i in range(n1):
        for j in range(i + 1, n1):
            f = (eps[i] * V[i] - eps[j] * V[j]) % C
            out.append(min(f, C - f))
    return tuple(sorted(out))


def factor(C):
    f = []
    m, d = C, 2
    while d * d <= m:
        while m % d == 0:
            f.append(d); m //= d
        d += 1
    if m > 1:
        f.append(m)
    return f


def is_prime(C):
    return C > 1 and all(C % d for d in range(2, int(C**0.5) + 1))


def brute_collisions(n):
    """Return (deficiency, class-size Counter, list of (flip-set, base-eps) for nontrivial classes)."""
    C = 2 * n - 1
    V = list(range(1, n))
    n1 = len(V)
    fold_map = defaultdict(list)
    for bits in product([0, 1], repeat=n1 - 1):
        eps = (1,) + tuple(1 if b else -1 for b in bits)  # fix runner-0 sign +
        fold_map[folded_spectrum(V, eps, C)].append(eps)
    sizes = Counter(len(g) for g in fold_map.values())
    deficiency = sum(len(g) for g in fold_map.values()) - len(fold_map)
    # extract nontrivial classes
    nontrivial = [g for g in fold_map.values() if len(g) > 1]
    return C, deficiency, sizes, nontrivial


def silent_flipsets(nontrivial, n):
    """For each nontrivial class, record (flip-set vs first member, deleted/flipped runners)."""
    out = []
    for g in nontrivial:
        base = g[0]
        for other in g[1:]:
            T = tuple(i + 1 for i in range(len(base)) if base[i] != other[i])  # runner labels
            out.append((T, base, other))
    return out


print("=" * 100)
print("PART A: exact collision pairs + flip-sets for small composite C (extract silent condition)")
print("=" * 100)
for n in range(5, 18):
    C = 2 * n - 1
    if is_prime(C):
        continue
    if (n - 2) > 16:   # keep brute force <= 2^16
        continue
    C, defic, sizes, nontrivial = brute_collisions(n)
    flips = silent_flipsets(nontrivial, n)
    # summarize distinct flip-sets and how many cuts each is silent for
    flip_counter = Counter(T for (T, b, o) in flips)
    print(f"\nC={C} (n={n}) fact={'x'.join(map(str,factor(C)))}  deficiency={defic}  sizes={dict(sizes)}")
    print(f"   distinct flip-sets (runner labels) -> #pairs:")
    for T, cnt in sorted(flip_counter.items(), key=lambda kv: (len(kv[0]), kv[0])):
        # is T a half-system of a subgroup of order q=C/gcd? check gcd content
        g = math.gcd(C, T[0]) if T else 0
        allg = set(math.gcd(C, t) for t in T)
        q = C // g if g else None
        print(f"      T={T}  gcd-content={sorted(allg)}  (subgroup order q=C/gcd={q})  pairs={cnt}")
