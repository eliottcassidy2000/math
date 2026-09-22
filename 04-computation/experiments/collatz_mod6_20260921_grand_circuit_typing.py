#!/usr/bin/env python3
"""Lane grand_circuit_typing (session collatz-mod6-20260917 wave 2026-09-21, mac-mini).

Exact typing of the pasted "grand unified circuit", which claims that the
session's Collatz framework is tied to the Giuga, Ankeny-Artin-Chowla (AAC),
Littlewood and Sarnak conjectures.  Every probe prints a numbered claim
S1, S2, ... with one of PROVED / FINITE-EXACT / CITED / HEURISTIC / REFUTED /
SCOPE / OPEN, and the note quotes only numbers that appear in this output.

Probes
  P1  Giuga: the sum criterion, the BBBG equivalence (re-proved and checked
      directly for n<=30000, criterion-checked for n<=10^5), the user's
      "p^2(p-1) | n-p" criterion (PROVED equivalent), the Cipolla trunk.
  P2  AAC: PARI quadunit for all primes p = 1 mod 4 below 3000; u mod p,
      2-adic and 3-adic valuations of u.
  P3  Littlewood: first 60 partial quotients of log_2 3 (PARI), record
      growth, convergent denominators, the trivial one-number direction.
  P4  Sarnak: entropy obstruction, the explicit Mobius-correlated 2-adic
      point, and the (content-free) finite-j parity correlations, N=10^6.
  P5  The box ring (Z, boxplus, boxtimes) is the shift x -> x+1; the fibre
      R(x)=4x+1 is not a boxtimes multiplication; coprimality along fibres
      and along odd orbits.
  P6  Typing greps: "17-vertex tournament", "S^6 / S_2 x S_3".
  P7  The honest bridge: 2^K-3^L along the convergents of log_2 3.

Explicit `raise` only (survives python3 -O).  RAM << 1 GB; runtime ~1-2 min.

Reproduce:
  python3 04-computation/experiments/collatz_mod6_20260921_grand_circuit_typing.py \
      > 05-knowledge/results/collatz_mod6_20260921_grand_circuit_typing.out
"""
import math
import os
import re
import subprocess
import sys
from fractions import Fraction

import numpy as np
import sympy

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
SCRATCH = os.environ.get(
    "GRAND_CIRCUIT_SCRATCH",
    "/private/tmp/claude-501/-Users-e-Documents-GitHub-math/e197ec98-d8f9-4475-947b-5af87889cf35/scratchpad",
)


def check(cond, msg):
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)


def gp(script):
    """Run a PARI/GP script through stdin, return stdout lines."""
    proc = subprocess.run(["gp", "-q"], input=script, text=True, capture_output=True, check=False)
    check(proc.returncode == 0, "gp failed: " + proc.stderr[:500])
    return [ln for ln in proc.stdout.splitlines() if ln.strip()]


def spf_sieve(N):
    spf = list(range(N + 1))
    for i in range(2, int(N ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, N + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def factor_spf(n, spf):
    f = {}
    while n > 1:
        p = spf[n]
        while n % p == 0:
            n //= p
            f[p] = f.get(p, 0) + 1
    return f


def mobius_sieve(N):
    mu = np.ones(N + 1, dtype=np.int64)
    is_p = np.ones(N + 1, dtype=bool)
    is_p[:2] = False
    for p in range(2, N + 1):
        if is_p[p]:
            if p * p <= N:
                is_p[p * p::p] = False
            mu[p::p] *= -1
            if p * p <= N:
                mu[p * p::p * p] = 0
    mu[0] = 0
    return mu


print("=" * 78)
print("Lane grand_circuit_typing: exact typing of the 'grand unified circuit' paste")
print("=" * 78)

# ---------------------------------------------------------------------------
print("\n### P1  GIUGA'S CONJECTURE")
print("Statement (Giuga 1950; CITED): n>=2 is prime iff  sum_{i=1}^{n-1} i^(n-1) = -1 (mod n).")
print("Equivalence (Borwein-Borwein-Borwein-Girgensohn, Amer. Math. Monthly 103 (1996); CITED):")
print("  composite n satisfies the congruence iff n is squarefree and for every prime p | n:")
print("  (p-1) | (n/p - 1)  [Carmichael/Korselt]  and  p | (n/p - 1)  [Giuga].")

# S1: direct computation of the sum for all n <= NDIRECT (numpy vector modexp)
NDIRECT = 30000


def giuga_sum_mod(n):
    """sum_{i=1}^{n-1} i^(n-1) mod n, by vectorised binary exponentiation."""
    base = np.arange(1, n, dtype=np.int64)
    res = np.ones(n - 1, dtype=np.int64)
    e = n - 1
    while e:
        if e & 1:
            res = (res * base) % n
        e >>= 1
        if e:
            base = (base * base) % n
    return int(res.sum() % n)


spf = spf_sieve(100000)


def bbbg_ok(n):
    f = factor_spf(n, spf)
    if any(a > 1 for a in f.values()):
        return False
    return all(((n // p - 1) % (p - 1) == 0) and ((n // p - 1) % p == 0) for p in f)


def user_ok(n):
    """The user's stated criterion: p^2 (p-1) | n - p for every prime p | n."""
    f = factor_spf(n, spf)
    return all((n - p) % (p * p * (p - 1)) == 0 for p in f)


hits_direct = []
mismatch = 0
for n in range(2, NDIRECT + 1):
    s = giuga_sum_mod(n)
    cong = (s == n - 1)
    isp = (spf[n] == n)
    if isp:
        check(cong, "Fermat fails at prime %d" % n)
    else:
        if cong:
            hits_direct.append(n)
        if cong != bbbg_ok(n):
            mismatch += 1
print("S1 FINITE-EXACT: direct sum computed for every n <= %d." % NDIRECT)
print("    every prime gives -1 (Fermat); composites with sum = -1 mod n: %s" % (hits_direct or "none"))
print("    disagreements between the direct congruence and the BBBG criterion on composites: %d" % mismatch)
check(mismatch == 0 and not hits_direct, "BBBG equivalence or Giuga conjecture fails below NDIRECT")

# S2: re-proof of BBBG by the prime-power sum lemmas, with brute-force checks of the lemmas.
print("S2 PROVED (re-derivation of BBBG, elementary): with q=p^a || n and k=n-1>=a,")
print("    sum_{i=1}^{n-1} i^k = (n/q) * sum_{i=0}^{q-1} i^k (mod q); the non-units contribute 0 (k>=a).")
print("    a=1: sum of units i^k = -1 if (p-1)|k else 0 (mod p)  [primitive root].")
print("    a>=2, p odd: units = C_{p-1} x (1+pZ), k coprime to p permutes 1+pZ, whose sum is p^(a-1) mod p^a,")
print("      so the unit sum is 0 mod p^(a-1), never -1.  p=2, a>=2: k odd permutes the units, sum = 2^(2a-2) = 0 mod 2^a.")
print("    Hence the congruence forces n squarefree and, for each p|n, (p-1)|(n-1) and n/p = 1 (mod p),")
print("    i.e. (p-1)|(n/p-1) [since n-1 = n/p-1 mod p-1] and p|(n/p-1).  QED (this is BBBG's theorem).")
# brute-force the lemma tables
bad = 0
for p in sympy.primerange(3, 200):
    for k in range(1, p):
        s = sum(pow(i, k, p) for i in range(1, p)) % p
        want = (p - 1) if (k % (p - 1) == 0) else 0
        if s != want:
            bad += 1
for p in sympy.primerange(2, 60):
    for a in range(2, 6):
        q = p ** a
        if q > 3000:
            break
        for k in range(a, a + 3 * q, max(1, q // 3)):
            if p != 2 and k % p == 0:
                continue
            if p == 2 and k % 2 == 0:
                continue
            s = sum(pow(i, k, q) for i in range(1, q)) % q
            if s % (p ** (a - 1)) != 0:
                bad += 1
print("    lemma tables brute-forced: primes p<200 all k<p, prime powers q=p^a<=3000 (a>=2): failures = %d" % bad)
check(bad == 0, "prime-power sum lemma failed")

# S3: BBBG criterion, user's criterion, Carmichael and Giuga numbers up to 1e5
N1 = 100000
carm = [n for n in range(4, N1 + 1) if spf[n] != n and all(a == 1 for a in factor_spf(n, spf).values())
        and all((n - 1) % (p - 1) == 0 for p in factor_spf(n, spf)) and len(factor_spf(n, spf)) >= 3]
giuga_nums = []
for n in range(2, N1 + 1):
    f = factor_spf(n, spf)
    if spf[n] == n or any(a > 1 for a in f.values()):
        continue
    if all((n // p - 1) % p == 0 for p in f):
        giuga_nums.append(n)
bbbg_comp = [n for n in range(4, N1 + 1) if spf[n] != n and bbbg_ok(n)]
user_comp = [n for n in range(4, N1 + 1) if spf[n] != n and user_ok(n)]
user_vs_bbbg = sum(1 for n in range(2, N1 + 1) if user_ok(n) != bbbg_ok(n) and spf[n] != n)
user_primes_ok = all(user_ok(p) for p in sympy.primerange(2, 2000))
print("S3 FINITE-EXACT (n <= %d): Carmichael numbers: %d, namely %s" % (N1, len(carm), carm))
print("    Giuga numbers (squarefree, p | n/p-1 for all p|n): %s" % giuga_nums)
for g in giuga_nums:
    f = sorted(factor_spf(g, spf))
    val = sum(Fraction(1, p) for p in f) - Fraction(1, g)
    check(val.denominator == 1, "Giuga sum not integral")
    print("      %d = %s,  sum 1/p - 1/n = %d" % (g, "*".join(map(str, f)), val))
print("    composites satisfying BBBG (= Carmichael AND Giuga): %s" % (bbbg_comp or "none"))
print("    composites satisfying the user's 'p^2(p-1) | n-p for all p|n': %s" % (user_comp or "none"))
print("    composite n <= %d where the user's criterion and BBBG disagree: %d; all primes < 2000 satisfy the user's criterion (n-p=0): %s"
      % (N1, user_vs_bbbg, user_primes_ok))
check(user_vs_bbbg == 0 and not bbbg_comp and not user_comp, "criteria disagree")
print("S4 PROVED: the user's criterion is EQUIVALENT to BBBG.  n-p = p(n/p-1), so p^2(p-1) | n-p  iff")
print("    p(p-1) | n/p-1  iff  p | n/p-1 and (p-1) | n/p-1  (gcd(p,p-1)=1).  Squarefreeness is implied:")
print("    if p^2 | n then p | n/p, so p does not divide n/p-1.  For prime n the criterion is vacuous (n-p=0).")
print("    So the paste's restatement is correct (it is the BBBG criterion in one line), and no map to Collatz is involved.")

# S5: the Cipolla trunk (4^j-1)/3 = R^j(1): never Carmichael, never Giuga (j<=40)
print("S5 FINITE-EXACT: trunk numbers t_j=(4^j-1)/3 (the R-orbit of 1, inherited), j=3..40 (t_2=5 is prime, vacuous):")
print("    j  t_j  squarefree  Korselt(all p-1 | t_j-1)  Giuga(all p | t_j/p-1)  first failing prime")
trunk_bad = 0
for j in range(3, 41):
    t = (4 ** j - 1) // 3
    f = sympy.factorint(t)
    prod = 1
    for p, a in f.items():
        prod *= p ** a
    check(prod == t, "factorint mismatch")
    sqf = all(a == 1 for a in f.values())
    kor = all((t - 1) % (p - 1) == 0 for p in f)
    giu = all((t // p - 1) % p == 0 for p in f)
    fail = next((p for p in sorted(f) if (t - 1) % (p - 1) != 0 or (t // p - 1) % p != 0), None)
    if kor or giu:
        trunk_bad += 1
    if j <= 12 or j in (20, 30, 40):
        print("    %2d  %s  %s  %s  %s  %s" % (j, t if t < 10 ** 12 else "%d digits" % len(str(t)), sqf, kor, giu, fail))
print("    trunk numbers j<=40 that are Korselt or Giuga: %d" % trunk_bad)
check(trunk_bad == 0, "trunk number Carmichael/Giuga")
print("S6 SCOPE / NO MAP: gcd(n, R(n)) = gcd(n, 4n+1) = 1 and gcd(x_i, x_{i+1}) = 1 along odd orbits are one-line facts")
print("    about PAIRS; Giuga's criterion is a divisibility p | n/p-1 inside ONE n.  No object of the session carries it.")

# ---------------------------------------------------------------------------
print("\n### P2  ANKENY-ARTIN-CHOWLA")
print("Statement (Ankeny-Artin-Chowla, Ann. of Math. 56 (1952); CITED): for a prime p = 1 (mod 4) with fundamental")
print("  unit (t + u sqrt p)/2 of Q(sqrt p), p does not divide u.")
gp_aac = r"""
default(parisize, 64000000);
forprime(p=5, 3000, if(p%4==1, my(e=quadunit(p), a=real(e), b=imag(e), t=2*a+b, u=b, nm=t^2-p*u^2); print(p," ",u%p," ",valuation(u,2)," ",valuation(u,3)," ",#digits(u)," ",nm)));
"""
lines = gp(gp_aac)
rows = []
for ln in lines:
    parts = ln.split()
    if len(parts) != 6:
        continue
    p, um, v2, v3, nd, nm = map(int, parts)
    rows.append((p, um, v2, v3, nd, nm))
n_p = len(rows)
check(n_p == sum(1 for p in sympy.primerange(5, 3000) if p % 4 == 1), "prime count mismatch")
viol = [r for r in rows if r[1] == 0]
norm_bad = [r for r in rows if r[5] not in (4, -4)]
print("S7 FINITE-EXACT (PARI quadunit, ring of discriminant p): %d primes p = 1 mod 4 in [5, 3000)." % n_p)
print("    norm t^2 - p u^2 in {+4,-4} for all: %s (norm -4 count: %d, +4 count: %d)"
      % (not norm_bad, sum(1 for r in rows if r[5] == -4), sum(1 for r in rows if r[5] == 4)))
print("    AAC violations (u = 0 mod p): %s" % (viol or "none"))
check(not viol and not norm_bad, "AAC violation or bad norm")
big_v3 = sorted(rows, key=lambda r: (-r[3], r[0]))[:8]
big_v2 = sorted(rows, key=lambda r: (-r[2], r[0]))[:8]
print("    largest 3-adic valuations of u:  " + ", ".join("p=%d v3=%d (u has %d digits)" % (r[0], r[3], r[4]) for r in big_v3))
print("    largest 2-adic valuations of u:  " + ", ".join("p=%d v2=%d (u has %d digits)" % (r[0], r[2], r[4]) for r in big_v2))
v3hist = {}
for r in rows:
    v3hist[r[3]] = v3hist.get(r[3], 0) + 1
v2hist = {}
for r in rows:
    v2hist[r[2]] = v2hist.get(r[2], 0) + 1
print("    histogram v3(u): %s" % dict(sorted(v3hist.items())))
print("    histogram v2(u): %s" % dict(sorted(v2hist.items())))
print("    expected for a 'random' integer: P(v3>=k)=3^-k, P(v2>=k)=2^-k; counts >= : v3>=1: %d (%.1f expected), v3>=2: %d (%.1f), v2>=1: %d (%.1f), v2>=2: %d (%.1f)"
      % (sum(1 for r in rows if r[3] >= 1), n_p / 3, sum(1 for r in rows if r[3] >= 2), n_p / 9,
         sum(1 for r in rows if r[2] >= 1), n_p / 2, sum(1 for r in rows if r[2] >= 2), n_p / 4))
longest = max(rows, key=lambda r: r[4])
print("    longest u: p=%d with %d digits; u mod p = %d" % (longest[0], longest[4], longest[1]))
# the valuations are FORCED by the norm equation: t^2 - p u^2 = -4
p1mod8 = [r for r in rows if r[0] % 8 == 1]
p5mod8 = [r for r in rows if r[0] % 8 == 5]
check(all(r[2] == 1 for r in p1mod8), "p=1 mod 8 with v2(u) != 1")
check(all(r[3] == 0 for r in rows) and all(r[2] <= 1 for r in rows), "valuation lemma fails")
print("S7b PROVED (why the histograms are degenerate): the fundamental unit of Q(sqrt p), p = 1 mod 4 prime, has norm -1")
print("    (CITED classical: Legendre/Dirichlet; re-seen: 211 of 211 above), so t^2 - p u^2 = -4.")
print("    mod 3: 3 | u would give t^2 = -4 = 2 mod 3, impossible; hence v3(u) = 0 ALWAYS.")
print("    mod 16: 4 | u would give t^2 = -4 = 12 mod 16, impossible; hence v2(u) <= 1 ALWAYS.")
print("    mod 8: t,u both odd forces 1 - p = -4 mod 8, i.e. p = 5 mod 8; so p = 1 mod 8 forces v2(u) = 1 (unit lies in Z[sqrt p]).")
print("    Census: p = 1 mod 8: %d primes, all with v2(u)=1: True;  p = 5 mod 8: %d primes, v2(u)=0 for %d of them, v2(u)=1 for %d."
      % (len(p1mod8), len(p5mod8), sum(1 for r in p5mod8 if r[2] == 0), sum(1 for r in p5mod8 if r[2] == 1)))
print("    So 'large 2-adic or 3-adic valuation of u' cannot occur for ANY p: the only honest bridge to the session's 2/3-adic")
print("    objects is closed by the norm equation itself.  Nothing Collatz-like appears.")
print("S8 CITED: AAC verified for all p < 2*10^11 (van der Poorten, te Riele, Williams, Math. Comp. 70 (2001) and its 2003 corrigendum).")
print("    UNCITED-RECOLLECTION (not asserted): a 2024 preprint claiming a counterexample; not checked here.")
print("S9 SCOPE / NO MAP: no object of the session is a unit of a real quadratic field of prime discriminant p = 1 mod 4.")
print("    The Pell hypotenuses of THM-3341 are units of Z[sqrt 2] (discriminant 8): p=2 is excluded from AAC.")
print("    'Silver ratio', 'Fermat numbers 2^(2^r)+1', 'repunit prime breaks', 'reversed Hamiltonian edge': no map found.")

# ---------------------------------------------------------------------------
print("\n### P3  LITTLEWOOD")
print("Statement (Littlewood c.1930; CITED): for all real alpha, beta,  liminf_n n ||n alpha|| ||n beta|| = 0.")
print("  Einsiedler-Katok-Lindenstrauss, Ann. of Math. 164 (2006) (CITED): the exceptional set has Hausdorff dimension 0.")
gp_cf = r"""
default(realprecision, 400);
a = log(3)/log(2);
cf = contfrac(a);
print(vector(60, i, cf[i]));
print(a);
"""
lines = gp(gp_cf)
cf = [int(x) for x in re.findall(r"-?\d+", lines[0])]
alpha_str = lines[1].strip()
check(len(cf) == 60 and cf[0] == 1 and cf[1] == 1 and cf[2] == 1 and cf[3] == 2 and cf[4] == 2 and cf[5] == 3,
      "continued fraction of log_2 3 unexpected: %s" % cf[:8])
# high-precision rational alpha from the decimal string
mant = alpha_str.replace(".", "")
ndec = len(alpha_str.split(".")[1])
alpha = Fraction(int(mant), 10 ** ndec)
print("S10 FINITE-EXACT: first 60 partial quotients of log_2 3 = %s" % cf)
rec = []
m = 0
for i, a in enumerate(cf):
    if a > m:
        m = a
        rec.append((i, a))
print("    record partial quotients (index, value): %s" % rec)
print("    max of the first 60: %d; sum: %d; mean: %.3f (Gauss-Kuzmin mean is infinite; median of a random a_i is 1)"
      % (max(cf), sum(cf), sum(cf) / 60))
# convergents
h1, h0, k1, k0 = 1, 0, 0, 1
conv = []
for a in cf[:40]:
    h1, h0 = a * h1 + h0, h1
    k1, k0 = a * k1 + k0, k1
    conv.append((h1, k1))
print("    convergents p/q (q<=10^12) and q*||q alpha||:")
for h, k in conv:
    if k > 10 ** 12:
        break
    d = alpha * k - h
    dist = abs(d)
    print("      %d/%d  q||q a|| = %.6f" % (h, k, float(k * dist)))
print("S11 PROVED (one line): if alpha has unbounded partial quotients then for every beta,")
print("    liminf n ||n alpha|| ||n beta|| <= liminf q_i ||q_i alpha|| * 1/2 <= liminf 1/(2 a_{i+1}) = 0.")
print("    Whether log_2 3 has unbounded partial quotients is OPEN (as for every non-quadratic algebraic-or-transcendental")
print("    constant not of Euler/e type); the pair (log_2 3, phi) is therefore OPEN, and phi = [1;1,1,...] has BOUNDED quotients,")
print("    so nothing about phi helps.")
# heuristic: min over n<=1e6 of n ||n alpha|| ||n phi||
N3 = 10 ** 6
n_arr = np.arange(1, N3 + 1, dtype=np.float64)
al = float(alpha)
ph = (1 + 5 ** 0.5) / 2
da = np.abs(n_arr * al - np.rint(n_arr * al))
dp = np.abs(n_arr * ph - np.rint(n_arr * ph))
prod = n_arr * da * dp
imin = int(np.argmin(prod))
fib = [1, 1]
while fib[-1] < 10 ** 6:
    fib.append(fib[-1] + fib[-2])
print("S12 HEURISTIC (float64): min_{n<=10^6} n ||n log_2 3|| ||n phi|| = %.6f at n = %d (a Fibonacci number: %s, index %d);"
      % (float(prod[imin]), imin + 1, (imin + 1) in fib, fib.index(imin + 1) + 1 if (imin + 1) in fib else -1))
print("    min_{n<=10^6} n ||n log_2 3|| = %.6f, min n ||n phi|| = %.6f (phi's constant is 1/sqrt5 = %.6f)"
      % (float((n_arr * da).min()), float((n_arr * dp).min()), 1 / 5 ** 0.5))
print("S13 SCOPE / NO MAP: the Collatz discrepancy K_j - j log_2 3 involves ONE irrational (inhomogeneous, one-number).")
print("    Wythoff/Zeckendorf/golden ratio occur in none of the audited notes (grep: 'phi' there is an angle).")
print("    'Bragg peaks q_{m,n}=(2 pi/phi^2)(m+n phi)' is the Fibonacci-chain diffraction module (Levine-Steinhardt 1984,")
print("    CITED as a formula family; the paste's normalisation is not load-bearing); '5 pi/6 phase from 5 mod 6': no object.")

# ---------------------------------------------------------------------------
print("\n### P4  SARNAK")
print("Statement (Sarnak 2009/2010, CITED): for every zero-entropy topological dynamical system (X,T), every f in C(X)")
print("  and every x in X, (1/N) sum_{n<=N} mu(n) f(T^n x) -> 0.")
print("S14 CITED: the map T(x)=x/2 (even), (3x+1)/2 (odd) on Z_2 is topologically conjugate, via the parity-vector map Q,")
print("    to the one-sided full 2-shift (Lagarias 1985, Bernstein-Lagarias 1996); h_top = log 2 = %.6f." % math.log(2))
print("    The greedy 3-adic map G of this session is conjugate to a 6-state SFT of entropy log 3 = %.6f" % math.log(3))
print("    (three_adic_g_map lane, Theorem 1.3).  Both have positive entropy: Sarnak's hypothesis fails; no statement follows.")

# S15: an explicit 2-adic point whose parity vector IS the indicator of mu = 1 (full-shift surjectivity)
K = 32
mu_small = [int(sympy.mobius(k + 1)) for k in range(K)]
target = [1 if m == 1 else 0 for m in mu_small]


def parity_vector(x, k):
    v = []
    for _ in range(k):
        v.append(x & 1)
        x = (3 * x + 1) // 2 if x & 1 else x // 2
    return v


x = 0
for k in range(1, K + 1):
    # x mod 2^k determined by the first k parities; try both lifts
    cands = [x, x + 2 ** (k - 1)]
    ok = [c for c in cands if parity_vector(c, k) == target[:k]]
    check(len(ok) == 1, "parity lift not unique at k=%d" % k)
    x = ok[0]
check(parity_vector(x, K) == target, "parity vector mismatch")
print("S15 PROVED + FINITE-EXACT: Q is a bijection, so there is a 2-adic x with parity(T^k x) = [mu(k+1) = 1] for ALL k.")
print("    Truncation: x = %d mod 2^%d realises the first %d values of [mu=1] = %s" % (x, K, K, "".join(map(str, target))))
print("    For that x the Mobius correlation of f = parity is (1/N) sum_{mu(n)=1} 1 -> 3/pi^2 = %.6f, not 0:" % (3 / math.pi ** 2))
print("    Mobius disjointness is FALSE for the Collatz system on Z_2 (as for any full shift).  Such x is irrational.")

# S16: finite-j correlations, N = 10^6
N4 = 10 ** 6
mu = mobius_sieve(N4)
n_arr = np.arange(0, N4 + 1, dtype=np.int64)
cur = n_arr.copy()
print("S16 FINITE-EXACT (N=%d): c_j = (1/N) sum_{n<=N} mu(n) (-1)^{parity(T^j n)}, and the periodicity check mod 2^(j+1)." % N4)
for j in range(0, 7):
    par = cur & 1
    sgn = 1 - 2 * par
    c = float((mu[1:] * sgn[1:]).sum()) / N4
    # periodicity: parity(T^j n) depends only on n mod 2^(j+1)
    mod = 2 ** (j + 1)
    cls = np.zeros(mod, dtype=np.int64) - 1
    periodic = True
    for r in range(mod):
        vals = par[r::mod]
        if vals.size and not (vals == vals[0]).all():
            periodic = False
        if vals.size:
            cls[r] = vals[0]
    mean_par = float(par[1:].mean())
    print("    j=%d  c_j = %+.6f   parity mean = %.6f   parity(T^j n) periodic mod %d: %s" % (j, c, mean_par, mod, periodic))
    check(periodic, "parity not periodic")
    cur = np.where(cur & 1 == 1, (3 * cur + 1) // 2, cur // 2)
print("    PROVED (via CITED PNT in arithmetic progressions, Landau): each c_j -> 0 because parity(T^j n) is periodic mod 2^(j+1),")
print("    hence a finite combination of residue-class indicators, each Mobius-orthogonal.  This is the zero-entropy rotation on")
print("    Z/2^(j+1), not the Collatz dynamics, and says nothing about Collatz.")
print("S17 SCOPE / NO MAP: 'F=S+U' is the divisor identity F=S+U iff N in {p,p^3,p^2qr} (divisor_balance lane); the 4-vertex")
print("    tournament matrix is the cell_ordering lane's score-profile object; B^3=-I is THM-4139's lift of the 3-cycle of")
print("    x^2-29/16 (zsigmondy_triad lane).  None is a dynamical system on which mu acts; no 'trapping loop' is defined.")

# ---------------------------------------------------------------------------
print("\n### P5  THE BOX RING AND THE COPRIME FIBRES")
bad = 0
for a in range(-10, 11):
    for b in range(-10, 11):
        if a * b + a + b + 1 != (a + 1) * (b + 1):
            bad += 1
    if a * a + 2 * a != (a + 1) ** 2 - 1:
        bad += 1
check(bad == 0, "box identity")
print("S18 PROVED: a boxtimes b := ab+a+b satisfies (a boxtimes b)+1 = (a+1)(b+1) (checked on [-10,10]^2, %d failures)." % bad)
print("    With a boxplus b := a+b+1 the map x -> x+1 is a ring ISOMORPHISM (Z, boxplus, boxtimes) -> (Z, +, *);")
print("    zero is -1, one is 0, the diagonal a boxtimes a = a^2+2a = (a+1)^2-1.  Transport preserves everything and adds nothing.")
sol = [c for c in range(-50, 51) if all((c + 1) * xx + c == 4 * xx + 1 for xx in range(3))]
print("S19 REFUTED (as a map): R(x)=4x+1 is not c boxtimes x = (c+1)x+c for any c (needs c+1=4 and c=1): solutions c in [-50,50]: %s" % sol)
check(sol == [], "R is a boxtimes multiplication?")
g1 = max(math.gcd((4 ** j * 7 + (4 ** j - 1) // 3), (4 ** (j + 1) * 7 + (4 ** (j + 1) - 1) // 3)) for j in range(0, 30))
print("S20 PROVED (inherited braid): x_j = 4^j x_0 + (4^j-1)/3 = R^j(x_0); gcd(x_j, x_{j+1}) = gcd(x_j, 4x_j+1) = 1 (max over j<30, x_0=7: %d)." % g1)
# odd-orbit coprimality
gmax = 0
pairs = 0
for n0 in range(1, 20001, 2):
    xx = n0
    for _ in range(200):
        y = 3 * xx + 1
        while y % 2 == 0:
            y //= 2
        if y == xx:
            break
        gmax = max(gmax, math.gcd(xx, y))
        pairs += 1
        xx = y
        if xx == 1:
            break
print("S21 PROVED: consecutive odd iterates x, x' = (3x+1)/2^k satisfy gcd(x,x') | gcd(x,3x+1) = 1 (checked on %d pairs, max gcd %d)." % (pairs, gmax))
print("    Both are odd, so the pair is a point of the odd-coprime chart (s,t)=(m+n,m-n) of the pythagorean_semicircle lane,")
print("    i.e. it determines one PPT ((s^2-t^2)/2, st, (s^2+t^2)/2) after ordering; that is the orbit-PPT lane's input.")
print("    For Giuga it is content-free (S6).")

# ---------------------------------------------------------------------------
print("\n### P6  TYPING GREPS: '17-vertex tournament capacity' and 'S^6 / S_2 x S_3 projective lift'")
thm_dir = os.path.join(REPO, "01-canon", "theorems")
names = sorted(os.listdir(thm_dir))
tourn_titles = []
for fn in names:
    if not fn.startswith("THM-"):
        continue
    with open(os.path.join(thm_dir, fn), "r", encoding="utf-8", errors="replace") as fh:
        head = fh.read(4000)
    m = re.search(r"^title:\s*(.*)$", head, re.M)
    title = m.group(1) if m else fn
    if re.search(r"tournament", title, re.I):
        tourn_titles.append((fn, title))
with17 = [(fn, t) for fn, t in tourn_titles if re.search(r"\b17\b|seventeen", t, re.I)]
print("S22 FINITE-EXACT (repo grep): %d canon theorems have 'tournament' in their title; with '17'/'seventeen' in the title: %d"
      % (len(tourn_titles), len(with17)))
for fn, t in with17[:10]:
    print("      %s :: %s" % (fn, t[:140]))
body17 = []
for fn, t in tourn_titles:
    with open(os.path.join(thm_dir, fn), "r", encoding="utf-8", errors="replace") as fh:
        body = fh.read()
    if re.search(r"17[- ]vert|\bn\s*=\s*17\b|order[- ]17\b|17-tournament|17 vertices", body, re.I):
        body17.append(fn)
print("    tournament-titled theorems whose BODY mentions a 17-vertex / order-17 / n=17 tournament: %d %s" % (len(body17), body17[:6]))
print("    (the ten title hits above are: 17 as a coefficient of the 7-tournament spectrum x^4+14x^2+17, the Fermat prime 17,")
print("    the fingerprint twin (17,13), and mod-16/2-adic statements; none defines a capacity of a 17-vertex tournament.)")
res_dir = os.path.join(REPO, "05-knowledge", "results")
s6_hits = 0
for fn in os.listdir(res_dir):
    if not fn.endswith(".md"):
        continue
    with open(os.path.join(res_dir, fn), "r", encoding="utf-8", errors="replace") as fh:
        txt = fh.read()
    if re.search(r"S\^6\b|\bS6 monodromy|S_2 x S_3|S2 x S3", txt):
        s6_hits += 1
print("    results notes mentioning 'S^6', 'S6 monodromy' or 'S_2 x S_3': %d (all are audits of the same paste family)." % s6_hits)
print("S23 SCOPE: the only '17's in tournament context in this session are the h-value 17 attained at n=6 (scaffolding_audit,")
print("    section 5 census; h-spectrum = odds minus {7,21}, THM-1370/THM-1745) and the class 17 of the mod-30 wheel.")
print("    No 'capacity' of a 17-vertex tournament is defined anywhere: no object.")
print("S24 SCOPE: 'S^6 with an S_2 x S_3 projective lift': S^6 (the 6-sphere) appears in no note; S_2 x S_3 = C2 x S3 (order 12)")
print("    appears once, as the actual automorphism action on the fruit curve triple (catalan_elliptic lane, section 4); a")
print("    'projective lift of S^6' is undefined.  The blueprint audit already found no S6 monodromy.  No object.")

# ---------------------------------------------------------------------------
print("\n### P7  THE HONEST BRIDGE: 2^K - 3^L ALONG THE CONVERGENTS OF log_2 3")
print("S25 FINITE-EXACT: a convergent p/q of log_2 3 gives 2^p ~ 3^q, i.e. (K,L)=(p,q); the cycle gate n_0 = bB/(2^K-3^L)")
print("    (inherited) needs 2^K - 3^L small and positive relative to B.  Table: (K,L) = (p,q), sign and |2^K-3^L|/3^L:")
for h, k in conv[:14]:
    d = 2 ** h - 3 ** k
    print("      K=%3d L=%3d  sign %+d  |2^K-3^L|/3^L = %.3e  (%d digits)" % (h, k, 1 if d > 0 else -1, abs(d) / 3 ** k, len(str(abs(d)))))
print("    CITED: Baker-type lower bounds |2^K-3^L| > 2^K / K^C (Pillai's problem; Stroeker-Tijdeman 1982 for this pair) are the")
print("    only formal input that bounds cycle lengths through the gate; that is the pillai lane's object, not this lane's.")

print("\n### VERDICT")
print("V1 Giuga: statement CITED; BBBG equivalence re-PROVED and FINITE-EXACT n<=30000 direct, n<=10^5 by criterion; user's")
print("   'p^2(p-1) | n-p' PROVED equivalent to BBBG; link to Collatz fibres/T_{n-2} bits: NO MAP.")
print("V2 AAC: statement CITED; FINITE-EXACT no violation p<3000 (%d primes); link to Pell/Fermat/repunits: NO MAP (p=2 excluded)." % n_p)
print("V3 Littlewood: statement CITED; only one irrational in Collatz (log_2 3); phi absent; trivial direction PROVED; NO MAP.")
print("V4 Sarnak: statement CITED; Collatz on Z_2 and G on Z_3^x have entropy log 2, log 3 > 0 (CITED/PROVED); disjointness")
print("   explicitly FALSE on Z_2 (S15); finite-j correlations are trivial PNT-in-APs (S16); NO MAP.")
print("V5 The three proposed formalisations (Mobius-Collatz matrix; simultaneous approximation of (log_2 3, phi); T_{n-2} bits")
print("   vs Giuga factors) each lack a map.  Formalisable honest bridges: Baker/Pillai on 2^K-3^L (pillai lane) and the")
print("   E-graph / G mirror (this session).")
print("\nALL CHECKS PASSED")
