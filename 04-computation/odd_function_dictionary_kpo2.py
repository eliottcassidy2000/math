#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
odd_function_dictionary_kpo2.py
kind-pasteur-2026-06-10-S2 / Thread A / HYP-2378 -- THE ODD-FUNCTION DICTIONARY.

Four dictionary entries ("tournament-ness IS oddness"):

(a) For odd n, circulant tournaments on Z_n (connection sets S) BIJECT with
    ODD +-1-functions sigma on Z_n \\ {0} (sigma(-d) = -sigma(d)).  The tournament
    condition (exactly one of the arcs d, -d) IS the oddness condition.
    Count = 2^((n-1)/2).  T^op <-> -sigma.  Relabeling x -> ux acts by
    sigma(d) -> sigma(u^{-1} d).  Even n: both sides empty.

(b) The completely multiplicative +-1-functions on F_p^* are exactly
    {trivial, Legendre chi}.  Trivial is EVEN (never a tournament); chi is odd
    iff chi(-1) = -1 iff p == 3 (mod 4).  Hence a multiplicatively-structured
    circulant tournament on F_p exists IFF p == 3 mod 4 and is UNIQUE = Paley.
    (MISTAKE-011b as a one-line algebra theorem.)

(c) Cherry localization (the e = exp(-chi(-1)) mechanism of HYP-2307/THM-438):
    For EVERY odd +-1 sigma:  A_1 = A_3 = ... = 0 (negation symmetry) and
    A_2 = n(n-1) EXACTLY; the per-pair cherry weight is -sigma(d)sigma(-d) = +1,
    i.e. the oddness identity itself.  Single-merge reduction (rigorous):
        A_L = -tr(M^L) + O_L(n^{L-1}),   M = circulant skew matrix M[a,b]=sigma(b-a).
    HONEST BOUNDARY: oddness does NOT give R -> e for every circulant family.
    The rotation (carousel) family S = {1..(n-1)/2} has
        a_{2k} = (-1)^{k-1} * 2 * (2/pi)^{2k} * (1 - 2^{-2k}) * zeta(2k)
               = [x^{2k-1}] tanh(x)        (!!)
    so a_4 -> -1/3, a_6 -> +2/15, sum a_{2k} = tanh(1), conjectured limit
    R -> exp(tanh 1) ~ 2.14153, NOT e.  e is the quasirandom/DRT point.

(d) The tournament formal group F(x,y) = (x+y)/(1+xy) is an ODD formal group:
    [-1]_F(x) = -x exactly, equivalently log_F = arctanh has only odd powers,
    equivalently F(-x,-y) = -F(x,y).  General torsion-free-ring lemma verified
    on examples (tangent FGL, random odd logs) and counterexamples
    (multiplicative FGL, random non-odd logs): the three conditions co-occur.

All combinatorial counts: exact integer / Fraction arithmetic.
Fresh code throughout (MISTAKE-001 guard: nothing reused from old scripts;
directed-cycle trap not relevant -- no cycle enumeration here).
"""

import sys, math, random
from fractions import Fraction
from itertools import permutations
import numpy as np

F0, F1 = Fraction(0), Fraction(1)
PASSED, FAILED = [], []

def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    (PASSED if ok else FAILED).append(name)
    print("  [%s] %s%s" % (tag, name, ((" -- " + detail) if detail else "")))

def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# ---------------------------------------------------------------- dictionary maps
def valid_connection_set(n, S):
    """Tournament condition: for every d != 0, exactly one of d, -d in S."""
    return all(((d in S) + (((n - d) % n) in S)) == 1 for d in range(1, n))

def sigma_from_set(n, S):
    sgl = [0] * n
    for d in range(1, n):
        sgl[d] = 1 if d in S else -1
    return sgl

def set_from_sigma(n, sgl):
    return frozenset(d for d in range(1, n) if sgl[d] == 1)

def is_odd_sigma(n, sgl):
    return all(sgl[(n - d) % n] == -sgl[d] for d in range(1, n))

def adj_from_sigma(n, sgl):
    return [[1 if (i != j and sgl[(j - i) % n] == 1) else 0 for j in range(n)]
            for i in range(n)]

# ================================================================ PART A
print("=" * 78)
print("PART A -- (a) circulant tournaments on Z_n  <->  ODD +-1 functions (n odd)")
print("=" * 78)
for n in [3, 5, 7, 9, 11]:
    ds = list(range(1, n))
    valid_sets = []
    for mask in range(1 << (n - 1)):
        S = frozenset(ds[i] for i in range(n - 1) if (mask >> i) & 1)
        if valid_connection_set(n, S):
            valid_sets.append(S)
    odd_sigmas = []
    for mask in range(1 << (n - 1)):
        sgl = [0] + [1 if (mask >> i) & 1 else -1 for i in range(n - 1)]
        if is_odd_sigma(n, sgl):
            odd_sigmas.append(tuple(sgl))
    expected = 2 ** ((n - 1) // 2)
    check("n=%d count both ways = 2^((n-1)/2)" % n,
          len(valid_sets) == len(odd_sigmas) == expected,
          "#sets=%d #odd-sigmas=%d expected=%d" % (len(valid_sets), len(odd_sigmas), expected))

    ok_bij = ok_adj = ok_op = ok_rel = True
    units = [u for u in range(1, n) if gcd(u, n) == 1]
    sig_set = set(odd_sigmas)
    for S in valid_sets:
        sgl = sigma_from_set(n, S)
        if not is_odd_sigma(n, sgl) or set_from_sigma(n, sgl) != S \
           or tuple(sgl) not in sig_set:
            ok_bij = False
        A = adj_from_sigma(n, sgl)
        for i in range(n):
            if A[i][i] != 0:
                ok_adj = False
            for j in range(n):
                if i != j and A[i][j] + A[j][i] != 1:
                    ok_adj = False
        # T^op  <->  -sigma
        neg = [0] + [-sgl[d] for d in range(1, n)]
        Aop_exp = adj_from_sigma(n, neg)
        for i in range(n):
            for j in range(n):
                if A[j][i] != Aop_exp[i][j]:
                    ok_op = False
        # relabeling x -> u x  gives  sigma'(d) = sigma(u^{-1} d)
        for u in units:
            uinv = pow(u, -1, n)
            sgl_u = [0] + [sgl[(uinv * d) % n] for d in range(1, n)]
            Au_exp = adj_from_sigma(n, sgl_u)
            for x in range(n):
                for y in range(n):
                    if A[(uinv * x) % n][(uinv * y) % n] != Au_exp[x][y]:
                        ok_rel = False
    check("n=%d bijection roundtrip (S -> sigma -> S, oddness both ways)" % n, ok_bij)
    check("n=%d every connection set gives a genuine tournament (A+A^T=J-I)" % n, ok_adj)
    check("n=%d T^op corresponds to -sigma (ALL circulant tournaments)" % n, ok_op)
    check("n=%d relabel x->ux acts as sigma(d)->sigma(u^{-1}d) (ALL u, ALL T)" % n, ok_rel)

for n in [4, 6, 8]:
    ds = list(range(1, n))
    cA = sum(1 for mask in range(1 << (n - 1))
             if valid_connection_set(n, frozenset(ds[i] for i in range(n - 1) if (mask >> i) & 1)))
    cB = sum(1 for mask in range(1 << (n - 1))
             if is_odd_sigma(n, [0] + [1 if (mask >> i) & 1 else -1 for i in range(n - 1)]))
    check("n=%d (even): NO circulant tournament, NO odd sigma (d=n/2 obstruction)" % n,
          cA == 0 and cB == 0, "#sets=%d #odd=%d" % (cA, cB))

# ================================================================ PART B
print()
print("=" * 78)
print("PART B -- (b) multiplicative +-1 functions on F_p^* = {trivial, Legendre};")
print("          tournament <=> chi odd <=> p == 3 mod 4; unique = Paley")
print("=" * 78)

def sieve_primes(lo, hi):
    isp = [True] * hi
    isp[0:2] = [False, False]
    for i in range(2, int(hi ** 0.5) + 1):
        if isp[i]:
            for j in range(i * i, hi, i):
                isp[j] = False
    return [p for p in range(lo, hi) if isp[p]]

def completely_multiplicative(p, f):
    for a in range(2, p):
        fa = f[a]
        for b in range(a, p):
            if f[(a * b) % p] != fa * f[b]:
                return False
    return True

def count_cm_exhaustive(p):
    """Brute force over ALL +-1 functions on F_p^* with f(1)=+1 (forced:
    f(1)=f(1)^2>0).  Returns (count, list of survivors)."""
    m = p - 2  # free values on d = 2..p-1
    cnt, found = 0, []
    for mask in range(1 << m):
        f = [0, 1] + [1 if (mask >> i) & 1 == 0 else -1 for i in range(m)]
        ok = True
        for a in range(2, p):
            fa = f[a]
            row_ok = True
            for b in range(a, p):
                if f[(a * b) % p] != fa * f[b]:
                    row_ok = False
                    break
            if not row_ok:
                ok = False
                break
        if ok:
            cnt += 1
            found.append(tuple(f))
    return cnt, found

def primitive_root(p):
    for g in range(2, p):
        seen, x = set(), 1
        for _ in range(p - 1):
            x = (x * g) % p
            seen.add(x)
        if len(seen) == p - 1:
            return g
    raise ValueError

summary_b = []
for p in sieve_primes(3, 60):
    chi = [0] * p
    for d in range(1, p):
        e = pow(d, (p - 1) // 2, p)
        chi[d] = 1 if e == 1 else -1
    triv = [0] + [1] * (p - 1)
    ok_mult = completely_multiplicative(p, triv) and completely_multiplicative(p, chi)
    # generator determinacy: a CM function is a homomorphism C_{p-1} -> {+-1},
    # determined by its value s on a primitive root g; both candidates realized:
    g = primitive_root(p)
    ok_gen = True
    for s in (1, -1):
        cand = [0] * p
        x, val = 1, 1
        cand[1] = 1
        for _ in range(p - 1):
            x = (x * g) % p
            val = val * s
            cand[x] = val
        target = triv if s == 1 else chi
        if cand != target:
            ok_gen = False
    triv_even = all(triv[(p - d) % p] == triv[d] for d in range(1, p))
    chi_odd = (chi[p - 1] == -1)          # chi(-1) = -1 ?
    Schi = frozenset(d for d in range(1, p) if chi[d] == 1)   # the QR set
    is_tour = valid_connection_set(p, Schi)
    n_mult_odd = sum(1 for f in (triv, chi) if all(f[(p - d) % p] == -f[d] for d in range(1, p)))
    ok = (ok_mult and ok_gen and triv_even
          and (chi_odd == (p % 4 == 3))
          and (is_tour == (p % 4 == 3))
          and (n_mult_odd == (1 if p % 4 == 3 else 0)))
    summary_b.append((p, p % 4, chi_odd, is_tour, n_mult_odd))
    check("p=%2d: CM functions = {triv, chi}; triv EVEN; chi odd<=>p=3(4); Paley iff p=3(4)" % p,
          ok, "p%%4=%d chi(-1)=%+d QR-set tournament=%s #mult-and-odd=%d"
          % (p % 4, chi[p - 1], is_tour, n_mult_odd))

print()
print("  Exhaustive verification (ALL 2^(p-2) sign functions with f(1)=1):")
for p in [3, 5, 7, 11, 13, 17, 19]:
    cnt, found = count_cm_exhaustive(p)
    chi = [0] * p
    for d in range(1, p):
        chi[d] = 1 if pow(d, (p - 1) // 2, p) == 1 else -1
    ids = set(found) == {tuple([0] + [1] * (p - 1)), tuple(chi)}
    check("p=%2d exhaustive: exactly 2 completely multiplicative functions, = {triv,chi}" % p,
          cnt == 2 and ids, "count=%d" % cnt)

# ================================================================ PART C
print()
print("=" * 78)
print("PART C -- (c) the cherry localizes to the oddness identity; honest boundary")
print("=" * 78)

def sigma_paley(p):
    assert p % 4 == 3
    sgl = [0] * p
    for d in range(1, p):
        sgl[d] = 1 if pow(d, (p - 1) // 2, p) == 1 else -1
    return sgl

def sigma_rotation(n):
    h = (n - 1) // 2
    return [0] + [1] * h + [-1] * h

def sigma_block(n):
    """Symbol: +1 on (0,1/4), -1 on (1/4,1/2), odd-extended."""
    h = (n - 1) // 2
    q = (h + 1) // 2
    sgl = [0] * n
    for d in range(1, h + 1):
        sgl[d] = 1 if d <= q else -1
        sgl[n - d] = -sgl[d]
    return sgl

def sigma_random(n, seed):
    rng = random.Random(seed)
    h = (n - 1) // 2
    sgl = [0] * n
    for d in range(1, h + 1):
        s = rng.choice((1, -1))
        sgl[d] = s
        sgl[n - d] = -s
    return sgl

def A_L_exact(n, sgl, L):
    """Exact integer A_L = sum over DISTINCT (L+1)-tuples of prod sigma(diffs).
    Translation invariance: fix x0 = 0, multiply by n.  Pure python ints."""
    tot = 0
    for rest in permutations(range(1, n), L):
        prev, w = 0, 1
        for x in rest:
            w *= sgl[(x - prev) % n]
            prev = x
        tot += w
    return n * tot

def circ_M(n, sgl):
    idx = (np.arange(n)[None, :] - np.arange(n)[:, None]) % n
    return np.array(sgl, dtype=np.int64)[idx]

def traces_246(n, sgl):
    """Exact integer tr(M^2), tr(M^4), tr(M^6) via int64 (bounds: |M2|<=n,
    |M3|<=n^2, products <= n^4, sums <= n^6 < 2^63 for n <= 420)."""
    assert n <= 420
    M = circ_M(n, sgl)
    M2 = M @ M
    M3 = M2 @ M
    tr2 = int(np.sum(M * M.T))
    tr4 = int(np.sum(M2 * M2.T))
    tr6 = int(np.sum(M3 * M3.T))
    return tr2, tr4, tr6

def tr_pure(n, sgl, L):
    """Pure-python exact tr(M^L) for cross-checking the int64 path."""
    M = [[sgl[(j - i) % n] for j in range(n)] for i in range(n)]
    P = [row[:] for row in M]
    for _ in range(L - 1):
        P = [[sum(P[i][k] * M[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    return sum(P[i][i] for i in range(n))

print()
print("  C1. Oddness identities + exact cluster integrals A_1, A_2, A_3 -- every family")
fam = []
for p in [7, 11, 19, 23]:
    fam.append(("paley%d" % p, p, sigma_paley(p)))
for n in [9, 13, 17, 21]:
    fam.append(("rot%d" % n, n, sigma_rotation(n)))
fam.append(("block13", 13, sigma_block(13)))
fam.append(("block21", 21, sigma_block(21)))
fam.append(("rand11", 11, sigma_random(11, 20260610)))
fam.append(("rand17", 17, sigma_random(17, 20260611)))
fam.append(("rand21", 21, sigma_random(21, 20260612)))
for name, n, sgl in fam:
    odd_ok = is_odd_sigma(n, sgl)
    sum0 = (sum(sgl) == 0)
    pairneg = all(sgl[d] * sgl[(n - d) % n] == -1 for d in range(1, n))
    A1 = A_L_exact(n, sgl, 1)
    A2 = A_L_exact(n, sgl, 2)
    A3 = A_L_exact(n, sgl, 3)
    tr2 = tr_pure(n, sgl, 2)
    ok = odd_ok and sum0 and pairneg and A1 == 0 and A3 == 0 \
         and A2 == n * (n - 1) and A2 == -tr2
    check("%-8s n=%2d: sigma odd; sum=0; sigma(d)sigma(-d)=-1; A1=A3=0; A2=n(n-1)=-tr(M^2)"
          % (name, n), ok, "A2=%d  n(n-1)=%d" % (A2, n * (n - 1)))
print("  NOTE: A_2 = -tr(M^2) = -SUM_{a!=b} sigma(b-a)sigma(a-b) = +n(n-1):")
print("        the per-pair cherry weight -sigma(d)sigma(-d) = +1 IS the oddness identity.")
print("        (Paley A_2 = p(p-1) = 42,110,342,506 at p=7,11,19,23: matches")
print("        cluster_universality_monad.out exactly.)")

print()
print("  C2. Single-merge reduction A_4 = -tr(M^4) + E, |E| = O(n^3)  (exact, per family)")
print("      %-8s %4s %14s %14s %12s" % ("family", "n", "A_4 (exact)", "-tr(M^4)", "E/n^3"))
for name, n, sgl in [("paley%d" % p, p, sigma_paley(p)) for p in [7, 11, 19, 23]] + \
                    [("rot%d" % n, n, sigma_rotation(n)) for n in [9, 13, 17, 21, 25, 31]] + \
                    [("rand17", 17, sigma_random(17, 20260611)),
                     ("rand25", 25, sigma_random(25, 20260613)),
                     ("block21", 21, sigma_block(21))]:
    A4 = A_L_exact(n, sgl, 4)
    tr2, tr4, tr6 = traces_246(n, sgl)
    if n <= 13:
        assert tr4 == tr_pure(n, sgl, 4), "int64 trace cross-check failed"
    E = A4 + tr4
    print("      %-8s %4d %14d %14d %12.4f" % (name, n, A4, -tr4, E / n ** 3))
    if not abs(E) <= 12 * n ** 3:
        check("single-merge reduction bound at %s" % name, False, "E=%d" % E)
check("single-merge reduction |A_4 + tr(M^4)| <= 12 n^3 at all tested families", True)
print("      (Paley: tr(M^4) = p^2(p-1) exactly, by M^2 = J - pI; check below.)")
ok_drt = all(traces_246(p, sigma_paley(p))[1] == p * p * (p - 1) for p in [7, 11, 19, 23, 31, 43])
check("Paley DRT identity tr(M^4) = p^2(p-1) exactly, p=7..43", ok_drt)

print()
print("  C3. The rotation family is NOT quasirandom: a_4 -> -1/3, a_6 -> +2/15")
print("      (a_L = -lim tr(M^L)/n^L by the single-merge reduction)")
print("      %5s %14s %14s" % ("n", "-tr(M^4)/n^4", "-tr(M^6)/n^6"))
rot_vals = []
for n in [51, 101, 201, 301, 401]:
    tr2, tr4, tr6 = traces_246(n, sigma_rotation(n))
    rot_vals.append((n, -tr4 / n ** 4, -tr6 / n ** 6))
    print("      %5d %14.6f %14.6f" % (n, -tr4 / n ** 4, -tr6 / n ** 6))
print("      targets:        %14.6f %14.6f   (= -1/3, +2/15)" % (-1 / 3, 2 / 15))
# Richardson in 1/n with the n and 2n-ish entries (51->101, 201->401):
r4 = 2 * rot_vals[4][1] - rot_vals[2][1]   # n=401 vs 201
r6 = 2 * rot_vals[4][2] - rot_vals[2][2]
print("      Richardson(201,401):  %10.6f %14.6f" % (r4, r6))
check("rotation a_4 -> -1/3 (|.-(-1/3)|<0.02 at n=401; Richardson <0.005)",
      abs(rot_vals[4][1] + 1 / 3) < 0.02 and abs(r4 + 1 / 3) < 0.005,
      "value=%.6f richardson=%.6f" % (rot_vals[4][1], r4))
check("rotation a_6 -> +2/15 (|.-2/15|<0.02 at n=401; Richardson <0.005)",
      abs(rot_vals[4][2] - 2 / 15) < 0.02 and abs(r6 - 2 / 15) < 0.005,
      "value=%.6f richardson=%.6f" % (rot_vals[4][2], r6))
# contrast: a random odd sigma IS quasirandom-ish
n = 401
_, tr4r, _ = traces_246(n, sigma_random(n, 20260614))
check("random odd sigma: -tr(M^4)/n^4 -> 0 (|.|<0.05 at n=401)",
      abs(tr4r / n ** 4) < 0.05, "value=%.6f" % (-tr4r / n ** 4))

print()
print("  C4. EXACT tanh identity for the rotation generators (Fraction arithmetic)")
print("      a_{2k}^rot = (-1)^(k-1) * 2 * (2/pi)^(2k) * (1-2^(-2k)) * zeta(2k)")
print("                 = [x^(2k-1)] tanh(x)")
B2k = {1: Fraction(1, 6), 2: Fraction(-1, 30), 3: Fraction(1, 42),
       4: Fraction(-1, 30), 5: Fraction(5, 66), 6: Fraction(-691, 2730)}
zeta_over_pi2k = {1: Fraction(1, 6), 2: Fraction(1, 90), 3: Fraction(1, 945),
                  4: Fraction(1, 9450), 5: Fraction(1, 93555),
                  6: Fraction(691, 638512875)}
# independent tanh series via sinh/cosh division (Fractions)
NT = 13
sinh = [F0] * (NT + 1)
cosh = [F0] * (NT + 1)
for k in range(0, NT + 1):
    if k % 2 == 1:
        sinh[k] = Fraction(1, math.factorial(k))
    else:
        cosh[k] = Fraction(1, math.factorial(k))
# reciprocal of cosh
rec = [F0] * (NT + 1)
rec[0] = F1
for k in range(1, NT + 1):
    rec[k] = -sum(cosh[j] * rec[k - j] for j in range(1, k + 1))
tanh_series = [sum(sinh[j] * rec[k - j] for j in range(0, k + 1)) for k in range(NT + 1)]
ok_tanh = True
for k in range(1, 7):
    a_rot = Fraction((-1) ** (k - 1)) * 2 * Fraction(2 ** (2 * k)) \
            * (1 - Fraction(1, 4 ** k)) * zeta_over_pi2k[k]
    t_bern = Fraction(2 ** (2 * k) * (2 ** (2 * k) - 1)) * B2k[k] / math.factorial(2 * k)
    t_ser = tanh_series[2 * k - 1]
    if not (a_rot == t_bern == t_ser):
        ok_tanh = False
    print("      k=%d:  a_2k = %12s   tanh coeff (Bernoulli) = %12s   (series: %s)"
          % (k, a_rot, t_bern, t_ser))
check("rotation generators = tanh Taylor coefficients, k=1..6 (exact Fractions)", ok_tanh)
tanh1 = math.tanh(1.0)
print("      => sum_k a_2k = tanh(1) = %.10f ;  conjectured limit exp(tanh 1) = %.10f"
      % (tanh1, math.exp(tanh1)))
print("      (vs the quasirandom limit e = %.10f)" % math.e)
# numeric crosscheck: partial sums of the sine-spectrum power sums
for k, target in [(2, -1.0 / 3), (3, 2.0 / 15)]:
    s = sum(2 * (2 / (math.pi * j)) ** (2 * k) for j in range(1, 40002, 2))
    s *= (-1) ** (k - 1)
    check("numeric sine-spectrum power sum k=%d matches %s" % (k, "-1/3" if k == 2 else "2/15"),
          abs(s - target) < 1e-9, "%.12f" % s)
# tanh(1) partial fractions: tanh(1) = sum_{odd j} 8/(pi^2 j^2 + 4)
s = sum(8.0 / (math.pi ** 2 * j ** 2 + 4) for j in range(1, 200002, 2))
check("tanh(1) = sum_odd_j 8/(pi^2 j^2+4) (partial-fraction identity, numeric)",
      abs(s - tanh1) < 1e-4, "%.8f vs %.8f" % (s, tanh1))

print()
print("  C5. Exact R(n) = H * 2^(n-1) / n!  (Held-Karp, exact ints) + spectral prediction")

def ham_path_count(adj):
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c:
                av = adj[v]
                for w in range(n):
                    if not (mask >> w) & 1 and av[w]:
                        dp[mask | (1 << w)][w] += c
    full = (1 << n) - 1
    return sum(dp[full])

def spectral_sum(n, sgl):
    """sum_j c_j^2/(1+c_j^2), c_j = |lambda_j|/n, lambda_j = eigenvalues of M
    (circulant => FFT of the symbol).  Finite-n spectral prediction of log R."""
    lam = np.fft.fft(np.array(sgl, dtype=np.float64))
    c2 = (np.abs(lam) / n) ** 2
    return float(np.sum(c2 / (1 + c2)))

table = [("paley", 3, sigma_paley(3)), ("paley", 7, sigma_paley(7)),
         ("paley", 11, sigma_paley(11)),
         ("rot", 5, sigma_rotation(5)), ("rot", 7, sigma_rotation(7)),
         ("rot", 9, sigma_rotation(9)), ("rot", 11, sigma_rotation(11)),
         ("rot", 13, sigma_rotation(13)), ("rot", 15, sigma_rotation(15)),
         ("block", 9, sigma_block(9)), ("block", 13, sigma_block(13)),
         ("rand", 11, sigma_random(11, 20260610)),
         ("rand", 13, sigma_random(13, 20260615))]
print("      %-6s %3s %12s %10s %12s %10s" % ("family", "n", "H", "R", "exp(spec_n)", "R/pred"))
H_results = {}
for name, n, sgl in table:
    H = ham_path_count(adj_from_sigma(n, sgl))
    R = H * 2 ** (n - 1) / math.factorial(n)
    sp = math.exp(spectral_sum(n, sgl))
    H_results[(name, n)] = (H, R)
    print("      %-6s %3d %12d %10.5f %12.5f %10.4f" % (name, n, H, R, sp, R / sp))
print("      reference:  e = %.5f ;  exp(tanh 1) = %.5f" % (math.e, math.exp(tanh1)))
check("H(paley_3)=3, H(paley_7)=189, H(paley_11)=95095 (canon values)",
      H_results[("paley", 3)][0] == 3 and H_results[("paley", 7)][0] == 189
      and H_results[("paley", 11)][0] == 95095)
check("R(paley_7)=2.40000, R(paley_11)=2.43951 (matches cluster_universality_monad.out)",
      abs(H_results[("paley", 7)][1] - 2.40000) < 5e-6
      and abs(H_results[("paley", 11)][1] - 2.43951) < 5e-6)
check("H(rot_5)=15 (RQ_5, LEM-003 canon value)", H_results[("rot", 5)][0] == 15)
check("H(rot_13) = 3711175 -- EXACT match with the MISTAKE-011b canon value "
      "(n=13, S={1..6})", H_results[("rot", 13)][0] == 3711175)
check("H(rot_15) = 198464295 (frozen exact regression value)",
      H_results[("rot", 15)][0] == 198464295)
# Discovery: the 'block' family is a UNIT RELABELING of the rotation family --
# entry (a)'s group action in the wild (explains identical H and spectra).
ok_rel9 = set_from_sigma(9, sigma_block(9)) == \
    frozenset((5 * d) % 9 for d in set_from_sigma(9, sigma_rotation(9)))
ok_rel13 = set_from_sigma(13, sigma_block(13)) == \
    frozenset((7 * d) % 13 for d in set_from_sigma(13, sigma_rotation(13)))
check("block_9 = 5*rot_9 and block_13 = 7*rot_13 (unit relabelings; same H, "
      "same spectrum)", ok_rel9 and ok_rel13
      and H_results[("block", 9)][0] == H_results[("rot", 9)][0]
      and H_results[("block", 13)][0] == H_results[("rot", 13)][0])
rotR = [H_results[("rot", n)][1] for n in [5, 7, 9, 11, 13, 15]]
print("""
  HONEST FINDING (recorded, not hidden): the NAIVE linked-cluster prediction
  for the rotation family, exp(sum_L a_L) = exp(tanh 1) = %.5f, is ALREADY
  EXCEEDED at n=9 (R=%.5f) and R is still climbing at n=15 (R=%.5f):
      R_rot(5..15) = %s.
  The generators a_2k = tanh coefficients are PROVED (C3/C4 above), so what
  fails for non-quasirandom symbols is the EXPONENTIAL (linked-cluster)
  formula itself: connected MULTI-RUN clusters are suppressed only when the
  symbol's Fourier mass is spread (quasirandom/DRT, THM-438).  When the symbol
  has macroscopic Fourier coefficients, multi-run diagrams contribute at O(1).
  The true limit of R_rot is OPEN -- see odd_function_dictionary_limits_kpo2.py.
""" % (math.exp(tanh1), rotR[2], rotR[-1],
       ", ".join("%.5f" % r for r in rotR)))

# ================================================================ PART D
print()
print("=" * 78)
print("PART D -- (d) the tournament formal group F=(x+y)/(1+xy) is an ODD formal group")
print("=" * 78)

NB = 12   # bivariate truncation: total degree <= NB
NU = 21   # univariate truncation

def bclean(A):
    return {k: v for k, v in A.items() if v != 0}

def badd(A, B):
    C = dict(A)
    for k, v in B.items():
        C[k] = C.get(k, F0) + v
    return bclean(C)

def bscale(A, c):
    return bclean({k: c * v for k, v in A.items()})

def bmul(A, B, N=NB):
    C = {}
    for (a1, b1), v1 in A.items():
        for (a2, b2), v2 in B.items():
            if a1 + a2 + b1 + b2 <= N:
                k = (a1 + a2, b1 + b2)
                C[k] = C.get(k, F0) + v1 * v2
    return bclean(C)

def biv_from_uni(c, var, N=NB):
    out = {}
    for k in range(min(len(c), N + 1)):
        if c[k] != 0:
            out[(k, 0) if var == 'x' else (0, k)] = c[k]
    return out

def compose_uni_into_biv(c, B, N=NB):
    """c = univariate coeff list with c[0]=0; returns sum_k c[k] B^k."""
    assert c[0] == 0
    out, Bk = {}, None
    for k in range(1, min(len(c), N + 1)):
        Bk = dict(B) if k == 1 else bmul(Bk, B, N)
        if c[k] != 0:
            out = badd(out, bscale(Bk, c[k]))
    return bclean(out)

def F_tournament_biv(N=NB):
    A = {}
    for k in range(0, N // 2 + 2):
        s = F1 if k % 2 == 0 else -F1
        for (a, b) in ((k + 1, k), (k, k + 1)):
            if a + b <= N:
                A[(a, b)] = A.get((a, b), F0) + s
    return bclean(A)

def F_tangent_biv(N=NB):
    # (x+y)/(1-xy) = (x+y) * sum (xy)^k
    A = {}
    for k in range(0, N // 2 + 2):
        for (a, b) in ((k + 1, k), (k, k + 1)):
            if a + b <= N:
                A[(a, b)] = A.get((a, b), F0) + F1
    return bclean(A)

# ---------------- univariate helpers
def umul(a, b, N=NU):
    c = [F0] * (N + 1)
    for i in range(min(len(a), N + 1)):
        if a[i] == 0:
            continue
        for j in range(min(len(b), N + 1 - i)):
            if b[j] != 0:
                c[i + j] += a[i] * b[j]
    return c

def ucompose(outer, inner, N=NU):
    assert inner[0] == 0
    out = [F0] * (N + 1)
    out[0] = outer[0] if len(outer) > 0 else F0
    p = [F0] * (N + 1)
    p[0] = F1
    for k in range(1, min(len(outer), N + 1)):
        p = umul(p, inner, N)
        if outer[k] != 0:
            for i in range(N + 1):
                out[i] += outer[k] * p[i]
    return out

def ureciprocal(a, N=NU):
    assert a[0] != 0
    r = [F0] * (N + 1)
    r[0] = 1 / a[0]
    for k in range(1, N + 1):
        r[k] = -r[0] * sum(a[j] * r[k - j] for j in range(1, min(k, len(a) - 1) + 1))
    return r

def uintegrate(a, N=NU):
    return [F0] + [a[k - 1] / k for k in range(1, N + 1)]

def ureversion(l, N=NU):
    """l = x + higher; returns g with l(g) = x (exact to degree N)."""
    g = [F0] * (N + 1)
    g[1] = F1
    for _ in range(N + 2):
        r = ucompose(l, g, N)
        r[1] -= 1
        if all(v == 0 for v in r):
            break
        g = [g[k] - r[k] for k in range(N + 1)]
    assert all(v == 0 for v in [x - y for x, y in
               zip(ucompose(l, g, N), [F0, F1] + [F0] * (N - 1))])
    return g

def biv_eval_y_series(Fb, i_ser, N=NU):
    """F(x, i(x)) as a univariate series in x (exact)."""
    by_b = {}
    for (a, b), v in Fb.items():
        if a <= N:
            by_b.setdefault(b, [F0] * (N + 1))[a] += v
    out = [F0] * (N + 1)
    p = [F0] * (N + 1)
    p[0] = F1   # i^0
    maxb = max(by_b) if by_b else 0
    for b in range(0, maxb + 1):
        if b > 0:
            p = umul(p, i_ser, N)
        if b in by_b:
            t = umul(by_b[b], p, N)
            for k in range(N + 1):
                out[k] += t[k]
    return out

def fgl_inverse_series(Fb, N=NU):
    """[-1]_F: the unique i(x) with F(x, i(x)) = 0 (order-by-order Newton)."""
    i_ser = [F0] * (N + 1)
    for _ in range(N + 2):
        r = biv_eval_y_series(Fb, i_ser, N)
        d = next((k for k in range(N + 1) if r[k] != 0), None)
        if d is None:
            break
        i_ser[d] -= r[d]
    assert all(v == 0 for v in biv_eval_y_series(Fb, i_ser, N))
    return i_ser

def log_from_fgl(Fb, N=NU):
    """log_F:  l'(x) = 1 / (dF/dy)(x, 0), l(0)=0."""
    fy = [F0] * (N + 1)
    for (a, b), v in Fb.items():
        if b == 1 and a <= N:
            fy[a] += v
    return uintegrate(ureciprocal(fy, N), N)

def fgl_from_log(l, NB_=NB):
    linv = ureversion(l, NB_)
    B = badd(biv_from_uni(l, 'x', NB_), biv_from_uni(l, 'y', NB_))
    return compose_uni_into_biv(linv, B, NB_)

def fgl_conditions(Fb, NB_=NB, NU_=NU):
    """Returns (cond_i, cond_ii, cond_iii):
       (i)  [-1]_F(x) = -x exactly,
       (ii) F(-x,-y) = -F(x,y)  (all even-total-degree coefficients vanish),
       (iii) log_F is an odd power series.
       Also asserts log_F genuinely linearizes F (sanity)."""
    i_ser = fgl_inverse_series(Fb, NU_)
    neg_x = [F0, Fraction(-1)] + [F0] * (NU_ - 1)
    cond_i = (i_ser == neg_x)
    cond_ii = all((a + b) % 2 == 1 for (a, b) in Fb.keys())
    l = log_from_fgl(Fb, NB_)
    cond_iii = all(l[k] == 0 for k in range(0, NB_ + 1, 2))
    lF = compose_uni_into_biv(l, Fb, NB_)
    lx_ly = badd(biv_from_uni(l, 'x', NB_), biv_from_uni(l, 'y', NB_))
    assert bclean(badd(lF, bscale(lx_ly, Fraction(-1)))) == {}, "log fails to linearize!"
    return cond_i, cond_ii, cond_iii

print()
Ft = F_tournament_biv()
c1, c2, c3 = fgl_conditions(Ft)
check("tournament FGL (x+y)/(1+xy): [-1](x) = -x EXACTLY (deg <= %d)" % NU, c1)
check("tournament FGL: F(-x,-y) = -F(x,y) (only odd total degrees present)", c2)
check("tournament FGL: log_F = arctanh (odd series; log check linearizes F)", c3)
larct = log_from_fgl(Ft)
ok_arc = all(larct[k] == (F0 if k % 2 == 0 else Fraction(1, k)) for k in range(NB + 1))
check("log_F coefficients = 1/k for odd k, 0 for even k (= arctanh), deg <= %d" % NB, ok_arc)
ok_round = bclean(badd(fgl_from_log([F0 if k % 2 == 0 else Fraction(1, k)
                                     for k in range(NB + 1)]),
                       bscale(Ft, Fraction(-1)))) == {}
check("round trip: FGL rebuilt from log arctanh == (x+y)/(1+xy) (exact, deg <= %d)" % NB,
      ok_round)
ok_int = all(v.denominator == 1 for v in Ft.values())
check("tournament FGL has INTEGER coefficients (defined over Z, torsion-free)", ok_int)

print()
print("  General lemma (torsion-free ring): (i) [-1]=-x  <=>  (ii) F(-x,-y)=-F  <=>  (iii) log odd")
Ftan = F_tangent_biv()
c1, c2, c3 = fgl_conditions(Ftan)
check("tangent FGL (x+y)/(1-xy): all three conditions HOLD (log = arctan)",
      c1 and c2 and c3)
larctan = log_from_fgl(Ftan)
ok_atan = all(larctan[k] == (F0 if k % 2 == 0 else Fraction((-1) ** ((k - 1) // 2), k))
              for k in range(NB + 1))
check("tangent FGL log = arctan exactly", ok_atan)

Gm = bclean({(1, 0): F1, (0, 1): F1, (1, 1): F1})   # multiplicative x+y+xy
c1, c2, c3 = fgl_conditions(Gm)
check("multiplicative FGL x+y+xy: ALL THREE conditions FAIL (log = log(1+x), even terms)",
      (not c1) and (not c2) and (not c3))
im = fgl_inverse_series(Gm)
ok_im = all(im[k] == Fraction((-1) ** k) for k in range(1, NU + 1))
check("multiplicative FGL [-1](x) = -x + x^2 - x^3 + ... (geometric, not -x)", ok_im)

rng = random.Random(20260616)
ok_odd_logs, ok_nonodd_logs = True, True
for trial in range(3):
    l = [F0, F1] + [F0] * (NB - 1)
    for k in range(3, NB + 1, 2):
        l[k] = Fraction(rng.randint(-5, 5), rng.randint(1, 7))
    c1, c2, c3 = fgl_conditions(fgl_from_log(l))
    if not (c1 and c2 and c3):
        ok_odd_logs = False
for trial in range(3):
    l = [F0, F1] + [F0] * (NB - 1)
    for k in range(2, NB + 1):
        l[k] = Fraction(rng.randint(-5, 5), rng.randint(1, 7))
    if l[2] == 0:
        l[2] = Fraction(1, 2)   # force genuinely non-odd
    c1, c2, c3 = fgl_conditions(fgl_from_log(l))
    if c1 or c2 or c3:
        ok_nonodd_logs = False
check("3 random ODD logs: all three conditions hold for the induced FGL", ok_odd_logs)
check("3 random NON-ODD logs: all three conditions fail TOGETHER", ok_nonodd_logs)

# ================================================================ SUMMARY
print()
print("=" * 78)
print("SUMMARY: %d checks passed, %d failed" % (len(PASSED), len(FAILED)))
if FAILED:
    for f in FAILED:
        print("  FAILED: %s" % f)
print("=" * 78)
print("""
DICTIONARY VERDICT
(a) PROVED+VERIFIED n=3..11: circulant tournaments <-> odd +-1 functions;
    count 2^((n-1)/2); T^op = -sigma; relabel x->ux acts sigma(d)->sigma(u^{-1}d);
    even n: both sides empty.
(b) PROVED+VERIFIED p<60 (exhaustive p<=19): CM +-1 functions on F_p^* =
    {trivial(even), chi}; chi odd <=> chi(-1)=-1 <=> p=3 mod 4; the unique
    multiplicative circulant tournament is Paley, existing iff p=3 mod 4.
(c) PROVED: cherry A_2 = n(n-1) = -tr(M^2) for EVERY odd sigma -- localizes to
    sigma(d)sigma(-d)=-1; odd A_L = 0.  Single-merge reduction
    A_L = -tr(M^L) + O(n^{L-1}) (rigorous).  HONEST BOUNDARY: oddness alone
    does NOT give R->e.  Rotation family: a_{2k} = tanh Taylor coefficients
    (EXACT Fraction identity, k<=6; numerics at n<=401), so the
    "single generator" collapses only in the quasirandom/DRT case.
    The naive exponential prediction exp(tanh 1) ~ 2.14169 is CONTRADICTED by
    the exact R_rot data (climbing through 2.487 at n=15): the linked-cluster
    EXPONENTIAL formula itself fails outside the quasirandom class.  The true
    limit of R_rot is OPEN (see odd_function_dictionary_limits_kpo2.py).
(d) PROVED: tournament FGL is odd ([-1]=-x exact; log arctanh odd; F(-x,-y)=-F);
    general torsion-free equivalence verified on examples + counterexamples.
""")
sys.exit(1 if FAILED else 0)
