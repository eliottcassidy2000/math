#!/usr/bin/env python3
"""Independent audit of section 3 (F, U, S; level sets; recursion; Klein quadrant; cusps and
Atkin-Lehner action) of supersingular_divisor_shapes_qseries_20261022.md.  (audit_ssj)

Methods deliberately different from fus_shapes.py / dr_genus.gp:
  * F, U, S from sympy.divisors / sympy.isprime / sympy.factorint (no closed forms used);
  * level sets: (i) exhaustive with sympy.divisor_count-free direct divisor lists for N <= 3*10^5,
    (ii) a PROOF-style enumeration over exponent profiles: delta = tau - 2^r - r - 1 + [sqfree] + [prime]
        and the bound tau >= 3*2^(r-1) for non-squarefree N, which closes the classification;
  * recursion: symbolic proof check over all exponent profiles with r <= 6, exps <= 6, plus random N;
  * cusps of X_0(N): enumerated directly as Gamma_0(N)-classes of P^1(Q) (Shimura's criterion
    a/c ~ a'/c' iff c' = u c mod N-ish is NOT used: we brute-force reduce by explicit matrices),
    then the Atkin-Lehner matrix W_Q is applied to each cusp and the denominator class recorded.
Run: python3 -B audit_fus_cusps.py > audit_fus_cusps.out
"""
from __future__ import annotations

import itertools
import random
from math import gcd

import sympy
from sympy import divisors, factorint, isprime, primerange

bad = 0


def fail(msg):
    global bad
    bad += 1
    print("FAIL:", msg)


def FUS(N):
    D = [d for d in divisors(N) if 1 < d < N]
    F = len(D)
    S = sum(1 for d in D if all(e == 1 for e in factorint(d).values()))
    U = sum(1 for d in D if isprime(d))
    return F, U, S


# ---------- 1. values at the three shapes ----------
P = list(primerange(2, 400))
cnt = {"p^2": 0, "pq": 0, "p^2q": 0}
for p in P:
    if FUS(p * p) != (1, 1, 1):
        fail(f"p^2 {p}")
    cnt["p^2"] += 1
for p, q in itertools.permutations(P[:40], 2):      # 40*39 ordered pairs (sympy is slow enough)
    if FUS(p * q) != (2, 2, 2):
        fail(f"pq {p} {q}")
    if FUS(p * p * q) != (4, 2, 3):
        fail(f"p^2q {p} {q}")
    cnt["pq"] += 1
    cnt["p^2q"] += 1
print("1. (F,U,S) = (1,1,1), (2,2,2), (4,2,3) at p^2, pq, p^2 q via sympy.divisors:", cnt)

# ---------- 2. level sets ----------
# (i) exhaustive, direct divisor lists
NMAX = 300_000
levels = {}
for N in range(2, NMAX + 1):
    F, U, S = FUS(N)
    d = F - U - S
    if d <= 0:
        prof = tuple(sorted(factorint(N).values(), reverse=True))
        levels.setdefault(d, set()).add(prof)
print(f"2(i). defect level sets d <= 0 for 2 <= N <= {NMAX} (exponent profiles):")
for d in sorted(levels):
    print("     d =", d, sorted(levels[d]))
# (ii) proof-style enumeration over profiles.  delta(profile) = tau - 2^r - r - 1 + [sqfree] + [prime]
def delta_prof(exps):
    r = len(exps)
    tau = 1
    for e in exps:
        tau *= e + 1
    sq = all(e == 1 for e in exps)
    pr = (r == 1 and exps[0] == 1)
    return tau - 2 ** r - r - 1 + sq + pr
# for non-squarefree profiles tau >= 3*2^(r-1), so delta >= 2^(r-1) - r - 1 > 0 once r >= 4:
assert all(2 ** (r - 1) - r - 1 > 0 for r in range(4, 200))
found = {}
for r in range(1, 4):
    for exps in itertools.combinations_with_replacement(range(1, 8), r):
        exps = tuple(sorted(exps, reverse=True))
        d = delta_prof(exps)
        if d <= 0 and not all(e == 1 for e in exps):
            found.setdefault(d, []).append(exps)
# exponents >= 8 only increase tau, so the list above is complete for r <= 3
print("2(ii). non-squarefree profiles with delta <= 0 (complete: tau monotone in exponents, r >= 4 impossible):", found)
print("       squarefree: r = 1 -> delta 0 (primes); r >= 2 -> delta = -r  (tau = 2^r)")
want = {-1: [(2,), (2, 1)], 0: [(3,), (2, 1, 1)]}
if {k: sorted(v) for k, v in found.items()} != {k: sorted(v) for k, v in want.items()}:
    fail("profile classification")

# ---------- 3. recursion delta(N r) = 2 delta(N) + U(N) - [N sqfree] - [N prime] ----------
def delta(N):
    F, U, S = FUS(N)
    return F - U - S
nrec = 0
random.seed(20261022)
for _ in range(3000):
    N = random.randint(2, 50_000)
    r = sympy.randprime(2, 200)
    if N % r == 0:
        continue
    sq = all(e == 1 for e in factorint(N).values())
    lhs = delta(N * r)
    rhs = 2 * delta(N) + FUS(N)[1] - sq - isprime(N)
    nrec += 1
    if lhs != rhs:
        fail(f"recursion N={N} r={r}")
# symbolic check over profiles (adjoining a new exponent-1 prime)
for rr in range(1, 7):
    for exps in itertools.combinations_with_replacement(range(1, 7), rr):
        e2 = tuple(sorted(exps + (1,), reverse=True))
        r = len(exps)
        tau = 1
        for e in exps:
            tau *= e + 1
        sq = all(e == 1 for e in exps)
        pr = (r == 1 and exps[0] == 1)
        U = r - pr
        if delta_prof(e2) != 2 * delta_prof(exps) + U - sq - pr:
            fail(f"symbolic recursion {exps}")
print(f"3. recursion: {nrec} random (N, r) pairs and all profiles r <= 6, e <= 6: ok" if bad == 0 else "3. recursion FAILED")

# ---------- 4. Klein orbits on D(p^2 q r), THM-4733 conventions ----------
pts = [(a, b, c) for a in range(3) for b in range(2) for c in range(2) if (a, b, c) not in ((0, 0, 0), (2, 1, 1))]
C = lambda v: (2 - v[0], 1 - v[1], 1 - v[2])
Q = lambda v: (v[0], v[2], v[1])
orbits = []
seen = set()
for v in pts:
    if v in seen:
        continue
    o = {v, C(v), Q(v), C(Q(v))}
    seen |= o
    orbits.append(o)
reps = []
for o in orbits:
    quad = [w for w in o if w[0] - w[1] - w[2] >= 0 and w[1] - w[2] >= 0]
    reps.append(quad)
print("4. Klein orbits (a,b,c):", [sorted(o) for o in orbits])
print("   quadrant representatives:", reps)
if sorted(r[0] for r in reps if len(r) == 1) != sorted([(1, 0, 0), (2, 0, 0), (1, 1, 0), (2, 1, 0)]) or any(len(r) != 1 for r in reps):
    fail("quadrant reps")

# ---------- 5. cusps of X_0(N) and the Atkin-Lehner action ----------
# Ground truth equivalence: a/c ~ a'/c' under Gamma_0(N) iff there are g1 = [[a,b],[c,e]], g2 = [[a',b'],[c',e']]
# in SL2(Z) and k with g2 T^k g1^{-1} (or its negative) in Gamma_0(N); its lower-left entry is
# c' e - c e' - c c' k, so the test is the linear congruence  c c' k == c' e - c e'  (mod N).
def egcd(a, b):
    """(g, s, t) with s a + t b = g = gcd(a, b) >= 0 (own implementation)"""
    s0, s1, t0, t1 = 1, 0, 0, 1
    while b != 0:
        q, r = divmod(a, b)
        a, b = b, r
        s0, s1 = s1, s0 - q * s1
        t0, t1 = t1, t0 - q * t1
    if a < 0:
        a, s0, t0 = -a, -s0, -t0
    return a, s0, t0


def complete(a, c):
    """b, e with a e - b c = 1"""
    if c == 0:
        return 0, 1 if a == 1 else -1
    g, s, t = egcd(a, -c)             # s a + t (-c) = 1  ->  e = s, b = t
    assert g == 1 and a * s - t * c == 1
    return t, s


def equivalent(a, c, a2, c2, N):
    b, e = complete(a, c)
    b2, e2 = complete(a2, c2)
    A = (c * c2) % N
    B = (c2 * e - c * e2) % N
    g = gcd(A, N)
    return B % g == 0


def key(a, c, N):
    """invariant (d, a*(c/d) mod g), g = gcd(d, N/d); invariance: for gamma = [[x,y],[z,w]] in Gamma_0(N),
    d' = d and a'(c'/d) == x w a (c/d) == a (c/d) mod g."""
    if c == 0:
        a, c = 1, N
    d = gcd(c, N)
    g = gcd(d, N // d)
    return (d, (a * (c // d)) % g if g > 1 else 0)


def apply(M, a, c):
    (x, y), (z, w) = M
    na, nc = x * a + y * c, z * a + w * c
    g = gcd(na, nc)                      # det may be Q > 1: reduce to lowest terms
    na, nc = na // g, nc // g
    if nc < 0:
        na, nc = -na, -nc
    if nc == 0:
        return 1, 0
    return na, nc


def W_matrix(Q, N):
    """an Atkin-Lehner matrix [[Q x, y], [N z, Q w]] of determinant Q (Q || N)"""
    M = N // Q
    for x in range(1, 500):
        for z in range(1, 500):
            if gcd(Q * x, M * z) != 1:
                continue
            _, u, v = egcd(Q * x, -M * z)            # u Q x - v M z = 1
            w, y = u, v
            assert Q * x * w - y * M * z == 1
            return ((Q * x, y), (N * z, Q * w))
    raise ValueError


print("5. cusps of X_0(N) (explicit Gamma_0(N)-equivalence) and Atkin-Lehner action on denominator classes")
rnd = random.Random(4733)
for N in (6, 9, 12, 18, 25, 30, 36, 49, 50, 60, 75, 98, 100, 147, 150, 180, 245, 275, 300, 637):
    # candidates: every a/c with c <= 2N (a mod c, coprime), plus infinity
    buckets = {}
    for c in range(1, 2 * N + 1):
        for a in range(0, c):
            if gcd(a, c) == 1:
                buckets.setdefault(key(a, c, N), []).append((a, c))
    buckets.setdefault(key(1, 0, N), []).append((1, N))
    # every candidate equivalent to its bucket's first element (explicit congruence test)
    okb = all(equivalent(*b[0], *x, N) for b in buckets.values() for x in (b if len(b) < 60 else rnd.sample(b, 60)))
    # distinct buckets are inequivalent (sampled pairs), confirming the invariant
    ks = list(buckets)
    oks = all(not equivalent(*buckets[k1][0], *buckets[k2][0], N) for k1 in ks for k2 in ks if k1 != k2)
    ncusps = len(buckets)
    formula = sum(int(sympy.totient(gcd(d, N // d))) for d in divisors(N))
    if not (okb and oks) or ncusps != formula:
        fail(f"cusps N={N}: {ncusps} classes vs formula {formula}, okb={okb}, oks={oks}")
    # Atkin-Lehner W_Q (Q || N): denominator class d -> (d/gcd(d,Q)) * (Q/gcd(d,Q))
    ok = True
    for pr, e in factorint(N).items():
        Qp = pr ** e
        W = W_matrix(Qp, N)
        for (dd, cls), b in buckets.items():
            a, c = b[0]
            a1, c1 = apply(W, a, c)
            d1 = N if c1 == 0 else gcd(c1, N)
            g = gcd(dd, Qp)
            if d1 != (dd // g) * (Qp // g):
                ok = False
    Wn = W_matrix(N, N)
    okF = True
    for (dd, cls), b in buckets.items():
        a1, c1 = apply(Wn, *b[0])
        if (N if c1 == 0 else gcd(c1, N)) != N // dd:
            okF = False
    if not (ok and okF):
        fail(f"Atkin-Lehner rule N={N}")
    tau = int(sympy.divisor_count(N))
    sqf = all(e == 1 for e in factorint(N).values())
    print(f"   N={N:4d}: #cusp classes = {ncusps:3d} = sum_d phi(gcd(d,N/d)) = {formula:3d}; tau(N) = {tau:2d}; "
          f"squarefree {sqf}; W_Q toggles Q-part of d: {ok}; Fricke d -> N/d: {okF}")
# the prime exchange q <-> r and the Atkin-Lehner group as permutations of the exponent box of p^2 q r
box = [(a, b, c) for a in range(3) for b in range(2) for c in range(2)]
refl = [lambda v: (2 - v[0], v[1], v[2]), lambda v: (v[0], 1 - v[1], v[2]), lambda v: (v[0], v[1], 1 - v[2])]
ALgroup = set()
for bits in itertools.product((0, 1), repeat=3):
    def f(v, bits=bits):
        for i, bt in enumerate(bits):
            if bt:
                v = refl[i](v)
        return v
    ALgroup.add(tuple(f(v) for v in box))
Cperm = tuple((2 - v[0], 1 - v[1], 1 - v[2]) for v in box)
Qperm = tuple((v[0], v[2], v[1]) for v in box)
CQperm = tuple((2 - v[0], 1 - v[2], 1 - v[1]) for v in box)
print("   Klein group <C,Q> meets the Atkin-Lehner reflection group in:",
      [nm for nm, pm in (("C", Cperm), ("Q", Qperm), ("CQ", CQperm)) if pm in ALgroup], "(+ identity)")
if [nm for nm, pm in (("C", Cperm), ("Q", Qperm), ("CQ", CQperm)) if pm in ALgroup] != ["C"]:
    fail("Klein / AL intersection")

print("\nfailures:", bad)
print("ALL AUDIT CHECKS PASSED" if bad == 0 else "AUDIT FOUND FAILURES")
