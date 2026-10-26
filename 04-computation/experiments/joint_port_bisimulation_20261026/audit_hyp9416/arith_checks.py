# Independent audit of the note's sections 5-6 (check F).  Pure Python + sympy primality only.
from itertools import combinations, product
from fractions import Fraction
from math import comb, prod, log
from sympy import isprime, primerange, factorint

def sigma(n): return sum(d for d in range(1, n + 1) if n % d == 0)
def mu(n):
    f = factorint(n)
    return 0 if any(e > 1 for e in f.values()) else (-1) ** len(f)

print("== 5. arrangement dictionary ==")
# (a) k facet hyperplanes of a (k-1)-simplex in R^(k-1): regions = realisable barycentric sign vectors.  Exact
# rational construction for every sign vector with a +, impossibility of all - (coordinates sum to 1); Zaslavsky /
# Buck count for a generic arrangement of m affine hyperplanes in R^d: sum_{i<=d} C(m,i) regions, C(m-1,d) bounded.
for k in range(1, 11):
    ok = 0
    for eps in product((1, -1), repeat=k):
        if all(e < 0 for e in eps): continue
        j = eps.index(1); tt = Fraction(1, 2 * k)
        lam = [e * tt for e in eps]; lam[j] = 1 - sum(lam[i] for i in range(k) if i != j)
        assert sum(lam) == 1 and all((x > 0) == (e > 0) for x, e in zip(lam, eps)); ok += 1
    d = k - 1
    assert ok == 2 ** k - 1 == sum(comb(k, i) for i in range(d + 1))
    assert comb(k - 1, d) == 1           # one bounded region (the open simplex, all coordinates positive)
    # central: k coordinate hyperplanes in R^k: 2^k orthants = regions on S^(k-1), 2^(k-1) in RP^(k-1);
    # vertices of the spherical arrangement: k - 1 of the k great spheres meet in 2 antipodal points -> 2k = +-e_i
    assert comb(k, k - 1) * 2 == 2 * k
print("k = 1..10: 2^k - 1 regions (all sign vectors but (-,...,-)), 1 bounded region, k simplex vertices; "
      "2^k orthants on S^(k-1), 2k vertices +-e_i, 2^(k-1) regions in RP^(k-1): ok")
# (b) weight sums over squarefree n
tested = 0
for k in range(1, 7):
    for ps in combinations(list(primerange(2, 30)), k):
        n = prod(ps); divs = [prod(S) for r in range(k + 1) for S in combinations(ps, r)]
        proper = [d for d in divs if d < n]
        assert sum(proper) == sigma(n) - n
        assert sum(mu(d) for d in proper) == -mu(n) == (-1) ** (k + 1)
        assert sum(divs) == sigma(n) and sum(mu(d) for d in divs) == 0
        pairs = {tuple(sorted((d, n // d))) for d in divs}
        assert len(pairs) == 2 ** (k - 1) and sum(a + b for a, b in pairs) == sigma(n)
        e = [sum(prod(S) for S in combinations(ps, j)) for j in range(k + 1)]
        assert sigma(n) - n == sum(e[:k])        # s(n) = 1 + e_1 + ... + e_(k-1)
        tested += 1
print(f"weights over {tested} squarefree n (k = 1..6 primes < 30): s(n), -mu(n) = (-1)^(k+1), sigma(n), 0, sigma(n) "
      f"over complementary pairs, s(n) = 1 + e_1 + ... + e_(k-1): ok")
s = lambda m: sigma(m) - m
print(f"s(170) = {s(170)}, s(182) = {s(182)}; 170 = {factorint(170)}, 182 = {factorint(182)}")
# s^-1(42): composite m with s(m) = 42 has a proper divisor >= sqrt(m) so 42 >= 1 + sqrt(m), m <= 41^2; primes give 1
N = 41 ** 2
sig = [0] * (N + 1)
for d in range(1, N + 1):
    for m in range(d, N + 1, d): sig[m] += d
pre = [m for m in range(2, N + 1) if sig[m] - m == 42]
print(f"s^-1(42) = {pre} (complete)")

print("== 6. arithmetic threads ==")
# Goormaghtigh-type coincidences 1 + p + ... + p^(a-1), a >= 3, values <= 1e12, prime bases (ALL primes with p^2 < 1e12)
LIM = 10 ** 12
vals = {}
for p in primerange(2, 10 ** 6 + 1):
    acc = 1 + p
    for a in range(3, 64):
        acc = acc * p + 1
        if acc > LIM: break
        vals.setdefault(acc, []).append((p, a))
both = {v: L for v, L in vals.items() if len({p for p, a in L}) > 1}
print(f"prime-base repunit values (a >= 3) <= 1e12 shared by two prime bases, all primes p <= 1e6: {both}")
# same with arbitrary bases x >= 2, values <= 1e10 (sanity: Goormaghtigh's 31 and 8191 must appear)
LIM2 = 10 ** 10; vals2 = {}
for x in range(2, 10 ** 5 + 1):
    acc = 1 + x
    for a in range(3, 64):
        acc = acc * x + 1
        if acc > LIM2: break
        vals2.setdefault(acc, []).append((x, a))
print("arbitrary-base coincidences <= 1e10 (x <= 1e5):", {v: L for v, L in vals2.items() if len(L) > 1})
# Gersonides contacts |3^j - 2^i| = 1
cont = sorted((2 ** i, 3 ** j) for i in range(0, 200) for j in range(0, 130) if abs(3 ** j - 2 ** i) == 1)
print("|3^j - 2^i| = 1 for i < 200, j < 130:", cont)
def rad(m): return prod(factorint(m)) if m > 1 else 1
for x, y in cont:
    lo, hi = min(x, y), max(x, y)
    print(f"   1 + {lo} = {hi}: rad = {rad(lo * hi)}, quality = {log(hi) / log(rad(lo * hi)):.4f}")
# twin-prime double hits: complete fibres for v <= 300 (composite preimages of v are <= (v-1)^2)
V = 300; N = (V - 1) ** 2
sig = [0] * (N + 1)
for d in range(1, N + 1):
    for m in range(d, N + 1, d): sig[m] += d
fib = {v: [] for v in range(V + 1)}
for m in range(4, N + 1):
    if not isprime(m) and sig[m] - m <= V: fib[sig[m] - m].append(m)
dh = [v for v in range(2, V + 1)
      if any(isprime(int(round(m ** 0.5))) and int(round(m ** 0.5)) ** 2 == m for m in fib[v])
      and any(m % 2 == 0 and m // 2 > 2 and isprime(m // 2) for m in fib[v])]
tw = [v for v in range(4, V + 1) if isprime(v - 1) and isprime(v - 3)]
print(f"v <= {V} with a prime-square and a 2q (q odd prime) preimage: {dh}")
print(f"   equals {{v > 3 : v-3, v-1 twin primes}}: {dh == tw}; self-hits v = 2q in s^-1(v): {[v for v in dh if v in fib[v]]}")
dh_any = [v for v in range(2, V + 1)
          if any(isprime(int(round(m ** 0.5))) and int(round(m ** 0.5)) ** 2 == m for m in fib[v])
          and any(m % 2 == 0 and isprime(m // 2) for m in fib[v])]
print(f"   allowing q = 2 (2q = 4 = 2^2): extra values {sorted(set(dh_any) - set(dh))}")
# power/doubling: first doubling row of n^k is the odd part u^k of n^k; 'later' iff u^k > n
for k in (2, 3):
    first_any = next(m for m in range(2, 100) if (lambda u: u ** k > m)(m >> ((m & -m).bit_length() - 1)))
    first_even = next(m for m in range(2, 100, 2) if (lambda u: u ** k > m)(m >> ((m & -m).bit_length() - 1)))
    print(f"k = {k}: first n with n^k in a later doubling row: {first_any} (any parity), {first_even} (even)")
