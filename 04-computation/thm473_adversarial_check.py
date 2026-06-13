#!/usr/bin/env python3
"""Adversarial verification of THM-473 (E[det(xI+S)] = Hermite / involutions).

Checks:
  (1) E[Pf(S[K])^2] = (|K|-1)!! by BRUTE FORCE over all sign assignments
      (|K| = 4 and 6), using the full Pfaffian definition (not assuming the
      cancellation argument). This tests whether the s_ji = -s_ij dependence
      breaks the monomial-orthogonality argument.
  (2) He_n(ix) = i^n * sum_k C(n,2k)(2k-1)!! x^(n-2k), symbolic n=3,4,5 (sympy).
  (3) sum_k C(n,2k)(2k-1)!! = I(n) (involutions, A000085), n=0..14.
  (4) E[det(I+S)] = I(n): exhaustive n=5,6 (exact integers), Monte Carlo n=7
      (500k samples).  Plus full expected char poly vs signless matching poly
      of K_n exhaustively at n=5,6, and E[det S] = (n-1)!! at n=4,6.
"""
import itertools, math
import numpy as np
import sympy as sp

ok = True
def check(label, cond, detail=""):
    global ok
    status = "PASS" if cond else "FAIL"
    if not cond: ok = False
    print(f"[{status}] {label} {detail}")

# ---------- helpers ----------
def dfact(m):  # double factorial, (-1)!! = 1
    return 1 if m <= 0 else m * dfact(m - 2)

def involutions(n):
    a, b = 1, 1  # I(0), I(1)
    if n == 0: return 1
    for k in range(2, n + 1):
        a, b = b, b + (k - 1) * a
    return b if n >= 1 else 1

def pfaffian_bruteforce(A):
    """Pfaffian via the full definition: sum over perfect matchings with sign,
    sign computed from the permutation (i1 j1 i2 j2 ...), pairs i<j sorted by i."""
    n = A.shape[0]
    assert n % 2 == 0
    idx = list(range(n))
    def rec(rem):
        if not rem: return [[]]
        i = rem[0]
        out = []
        for j in rem[1:]:
            rest = [x for x in rem if x != i and x != j]
            for m in rec(rest):
                out.append([(i, j)] + m)
        return out
    total = 0
    for m in rec(idx):
        perm = [x for pair in m for x in pair]
        # sign of the permutation (i1 j1 i2 j2 ...) via inversion parity
        inv = sum(1 for a in range(n) for b in range(a+1, n) if perm[a] > perm[b])
        sgn = -1 if inv % 2 else 1
        total += sgn * math.prod(A[i, j] for (i, j) in m)
    return total

def skew_from_bits(n, bits):
    S = np.zeros((n, n), dtype=np.int64)
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            s = 1 if (bits >> k) & 1 else -1
            S[i, j] = s; S[j, i] = -s
            k += 1
    return S

# ---------- (1) E[Pf^2] brute force ----------
print("== (1) E[Pf^2] = (|K|-1)!! brute force ==")
for m in (4, 6):
    e = m * (m - 1) // 2
    tot = 0
    for bits in range(1 << e):
        S = skew_from_bits(m, bits)
        p = pfaffian_bruteforce(S)
        # cross-check Pf^2 = det
        d = round(np.linalg.det(S.astype(float)))
        assert p * p == d, (m, bits, p, d)
        tot += p * p
    avg = tot / (1 << e)
    check(f"|K|={m}: E[Pf^2]", avg == dfact(m - 1), f"avg={avg} vs {dfact(m-1)}  (Pf^2==det verified all {1<<e} cases)")

# ---------- (2) Hermite identity, symbolic ----------
print("== (2) He_n(ix) = i^n sum C(n,2k)(2k-1)!! x^(n-2k) ==")
x = sp.symbols('x')
for n in (3, 4, 5):
    lhs = sp.expand(sp.functions.special.polynomials.hermite_prob(n, sp.I * x))
    rhs = sp.expand(sp.I**n * sum(sp.binomial(n, 2*k) * dfact(2*k - 1) * x**(n - 2*k)
                                  for k in range(n // 2 + 1)))
    check(f"n={n}", sp.simplify(lhs - rhs) == 0, f"He_{n}(ix)={lhs}")

# ---------- (3) involution identity ----------
print("== (3) sum_k C(n,2k)(2k-1)!! = I(n) ==")
for n in range(0, 15):
    s = sum(math.comb(n, 2*k) * dfact(2*k - 1) for k in range(n // 2 + 1))
    # independent I(n): count involutions by brute permanence for n<=7, recurrence otherwise
    if n <= 7:
        cnt = sum(1 for p in itertools.permutations(range(n))
                  if all(p[p[i]] == i for i in range(n)))
    else:
        cnt = involutions(n)
    check(f"n={n}", s == cnt, f"sum={s} I(n)={cnt}")

# ---------- (4) E[det(I+S)] etc ----------
def all_skews(n):
    """All 2^C(n,2) skew sign matrices, shape (N, n, n), float64."""
    e = n * (n - 1) // 2
    N = 1 << e
    bits = np.arange(N, dtype=np.uint32)
    iu = np.triu_indices(n, 1)
    S = np.zeros((N, n, n))
    for k in range(e):
        s = ((bits >> k) & 1) * 2.0 - 1.0
        S[:, iu[0][k], iu[1][k]] = s
        S[:, iu[1][k], iu[0][k]] = -s
    return S

def int_dets(M):
    """Exact integer dets of small-entry integer matrices via float det + round."""
    d = np.linalg.det(M)
    r = np.rint(d)
    assert np.max(np.abs(d - r)) < 1e-6, "float det not safely integral"
    return r.astype(np.int64)

print("== (4a) exhaustive n=5,6: E[det(I+S)], full char poly, E[det S] ==")
for n in (5, 6):
    e = n * (n - 1) // 2
    N = 1 << e
    S = all_skews(n)
    tot_det = int(int_dets(np.eye(n)[None] + S).sum())
    Iexp = involutions(n)
    check(f"n={n}: E[det(I+S)]", tot_det % N == 0 and tot_det // N == Iexp,
          f"{tot_det}/{N} = {sp.Rational(tot_det, N)} vs I({n})={Iexp}")
    # char poly coeff of x^(n-k) = sum over |K|=k principal minors of S
    coeff_sums = [N if k == 0 else 0 for k in range(n + 1)]
    for k in range(1, n + 1):
        for K in itertools.combinations(range(n), k):
            sub = S[np.ix_(range(N), K, K)]
            coeff_sums[k] += int(int_dets(sub).sum())
    expect = [0] * (n + 1)
    for k in range(n // 2 + 1):
        expect[2*k] = math.comb(n, 2*k) * dfact(2*k - 1)  # coeff of x^(n-2k)
    got = [sp.Rational(c, N) for c in coeff_sums]
    check(f"n={n}: E[char poly] = signless matching poly K_{n}",
          got == [sp.Rational(v) for v in expect], f"got={got}")
    if n % 2 == 0:
        tot_detS = int(int_dets(S).sum())
        check(f"n={n}: E[det S] = (n-1)!!", tot_detS % N == 0 and tot_detS // N == dfact(n-1),
              f"{tot_detS}/{N} vs {dfact(n-1)}")
    # exact-arithmetic spot check (guard against float det issues)
    rng0 = np.random.default_rng(7)
    lam = sp.symbols('lam')
    for bits in rng0.integers(0, N, size=20):
        Sx = skew_from_bits(n, int(bits))
        cp = sp.Matrix(Sx.tolist()).charpoly(lam).all_coeffs()
        # compare against minor-based coeffs for this single matrix
        for k in range(1, n + 1):
            mk = sum(int(round(np.linalg.det(Sx[np.ix_(K, K)].astype(float))))
                     for K in itertools.combinations(range(n), k))
            assert mk == int(cp[k]), (n, int(bits), k, mk, int(cp[k]))
    print(f"      n={n}: sympy exact char-poly spot check on 20 random matrices OK")

print("== (4a') exhaustive n=7: E[det(I+S)] (bonus, confirms exact claim 486539264/2097152) ==")
n = 7
e = n * (n - 1) // 2
N = 1 << e
iu = np.triu_indices(n, 1)
tot7 = 0
chunk = 1 << 17
for start in range(0, N, chunk):
    bits = np.arange(start, start + chunk, dtype=np.uint32)
    S = np.zeros((chunk, n, n))
    for k in range(e):
        s = ((bits >> k) & 1) * 2.0 - 1.0
        S[:, iu[0][k], iu[1][k]] = s
        S[:, iu[1][k], iu[0][k]] = -s
    tot7 += int(int_dets(np.eye(n)[None] + S).sum())
check("n=7: E[det(I+S)] exhaustive", tot7 == 486539264 and tot7 // N == involutions(7) and tot7 % N == 0,
      f"{tot7}/{N} = {sp.Rational(tot7, N)} vs I(7)={involutions(7)}")

print("== (4b) Monte Carlo n=7, 500k samples ==")
rng = np.random.default_rng(20260611)
n, NS = 7, 500_000
iu = np.triu_indices(n, 1)
signs = rng.integers(0, 2, size=(NS, len(iu[0]))) * 2 - 1
S = np.zeros((NS, n, n))
S[:, iu[0], iu[1]] = signs
S -= S.transpose(0, 2, 1)
M = np.eye(n)[None, :, :] + S
dets = np.linalg.det(M)
mean = dets.mean()
sem = dets.std(ddof=1) / math.sqrt(NS)
Iexp = involutions(7)
zin = abs(mean - Iexp) / sem
check(f"n=7 MC: E[det(I+S)] ~ I(7)={Iexp}", zin < 4,
      f"mean={mean:.3f} sem={sem:.3f} z={zin:.2f}")

print("ALL PASS" if ok else "SOME CHECKS FAILED")
