# Divisor lattice of a squarefree n = p_1...p_k as cube, simplex, sphere and hyperplane arrangement, with prime weights.
#   all divisors        <-> vertices of the k-cube {0,1}^k                 <-> 2^k regions of k great spheres on S^(k-1)
#   proper divisors     <-> 2^k - 1 regions of k generic affine hyperplanes in R^(k-1) (the facet hyperplanes of a
#                           (k-1)-simplex: region "lambda_i < 0 exactly for i in S" is nonempty iff S is proper)
#   divisors > 1, < n   <-> faces of the boundary of the (k-1)-simplex, a (k-2)-sphere (a circle for k = 3)
#   pairs {d, n/d}      <-> 2^(k-1) regions of k hyperplanes in RP^(k-1)
# Weights: unit -> counts; Moebius -> (reduced) Euler characteristics; d itself -> sigma(n), s(n) = sigma(n) - n.
# mac-mini-2026-10-25-pyth (2026-10-26).
import itertools, random
from math import prod
from sympy import primerange, mobius, divisor_sigma, factorint

def realized(eps):
    """exact: a point with barycentric sign vector eps (sum of coordinates 1) exists iff some entry is +.
    Construction: lambda_i = eps_i * t for i != j (j an index with eps_j = +), lambda_j = 1 - sum; t = 1/(2k)."""
    k = len(eps)
    if all(e < 0 for e in eps): return None
    j = eps.index(1); t = 1.0 / (2 * k)
    lam = [e * t for e in eps]; lam[j] = 1 - sum(lam[i] for i in range(k) if i != j)
    assert abs(sum(lam) - 1) < 1e-12 and all((x > 0) == (e > 0) for x, e in zip(lam, eps))
    return lam

print("affine arrangement of the k facet hyperplanes of a (k-1)-simplex: realized sign vectors (exact construction)")
for k in range(2, 9):
    vecs = list(itertools.product((1, -1), repeat=k))
    ok = [e for e in vecs if realized(e) is not None]
    print(f"  k = {k}: {len(ok)} of {2 ** k} realized = 2^k - 1, the missing one all-negative: {len(ok) == 2 ** k - 1}")

def report(ps):
    n = prod(ps); k = len(ps)
    divs = sorted(prod(S) for r in range(k + 1) for S in itertools.combinations(ps, r))
    proper = [d for d in divs if d < n]
    bdry = [d for d in proper if d > 1]
    s = sum(proper)
    e = [sum(prod(S) for S in itertools.combinations(ps, j)) for j in range(k + 1)]
    print(f"n = {n} = {'*'.join(map(str, ps))}: proper divisors {proper}")
    print(f"   s(n) = {s} = 1 + " + " + ".join(f"e{j}" for j in range(1, k)) + f" = {' + '.join(str(e[j]) for j in range(k))}  (e_j elementary symmetric in the primes)")
    print(f"   counts: divisors {len(divs)} = 2^k, proper {len(proper)} = 2^k - 1, boundary faces of the simplex {len(bdry)} = 2^k - 2")
    mu = sum(mobius(d) for d in proper)
    print(f"   Moebius-weighted proper sum {mu} = -mu(n) = {-mobius(n)}; reduced Euler characteristic of S^(k-2) = {(-1) ** k}")
    print(f"   pairs {{d, n/d}}: {sorted({tuple(sorted((d, n // d))) for d in divs})}  -> {2 ** (k - 1)} projective regions")
    return s

print()
for ps in ([2, 3, 5], [2, 3, 7], [3, 5, 7], [2, 3, 5, 7]):
    report(ps)
# What the geometric picture loses: all n = pqr share the same complex; s(n) = 1 + e1 + e2.  Which n = pqr have the same
# s(n)?  And which other n (any shape) share s(30) = 42?
N = 200000
from collections import defaultdict
sig = [0] * (N + 1)
for d in range(1, N + 1):
    for m in range(d, N + 1, d): sig[m] += d
s_of = [sig[m] - m for m in range(N + 1)]
pre42 = [m for m in range(2, N + 1) if s_of[m] == 42]
print(f"\ns^-1(42) (complete, since composite preimages are <= 41^2): {pre42}")
for m in pre42: print(f"   {m} = {factorint(m)}")
# three-prime squarefree n sharing s(n): collisions
byS = defaultdict(list)
ps = list(primerange(2, 60))
for a, b, c in itertools.combinations(ps, 3):
    byS[1 + (a + b + c) + (a * b + a * c + b * c)].append((a, b, c))
coll = {v: L for v, L in byS.items() if len(L) > 1}
print(f"three-prime squarefree n (primes < 60) with equal s(n): {len(coll)} values; smallest: " +
      ", ".join(f"s = {v}: {L}" for v, L in sorted(coll.items())[:4]))
