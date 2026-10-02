#!/usr/bin/env python3
"""drop_multiplicities_supplement_20261002_run.py -- the computations behind
05-knowledge/results/drop_multiplicities_supplement_20261002.md (Collatz drops A - Syr(A)), a supplement to
THM-4530. Written without reading the THM-4530 code; it also re-derives THM-4527 items 1-2 and much of THM-4530.

Results go to stdout (deterministic); elapsed times go to stderr.  Exits with status 1 on any
failure.  Last stdout line: ALL CHECKS PASSED.  Single process, single core, about 45 s and 1.1 GB;
numpy + sympy + mpmath + fractions.

Reproduce:  nice -n 10 python3 04-computation/experiments/drop_multiplicities_supplement_20261002_run.py \
                > 05-knowledge/results/drop_multiplicities_supplement_20261002.out

Notation: Syr(A) = (3A+1)/2^v, v = v_2(3A+1), A odd;  d(A) = (A - Syr(A))/2;
q_v = 2^v - 3;  m(d) = #{v >= 2 : q_v | 6d+1}  (d >= 1).
"""
import sys, time, itertools, math
from math import gcd, comb
from fractions import Fraction
import numpy as np
import sympy
import mpmath

T_START = time.time()
NCHECK = 0


def check(cond, msg):
    global NCHECK
    NCHECK += 1
    if not cond:
        print("FAIL", msg)
        print("CHECKS FAILED")
        sys.exit(1)
    print("PASS", msg)


def info(msg):
    print("   ", msg)


def section(title):
    print("=== %s" % title)
    sys.stderr.write("[t=%.1fs] %s\n" % (time.time() - T_START, title))


def r5(x):
    """round to 5 significant digits (for comparing with the note's tables)"""
    return float("%.5g" % float(x))


def lcm(a, b):
    return a // gcd(a, b) * b


def v2(x):
    return (x & -x).bit_length() - 1


def syr(A):
    """Syracuse step for an odd integer A (any sign, A != -1/3): returns (Syr(A), v)."""
    t = 3 * A + 1
    v = v2(t)
    return t >> v, v


def drop(A):
    s, v = syr(A)
    return (A - s) // 2


def q(v):
    return 2 ** v - 3


def m_formula(d):
    """multiplicity predicted by the theorem, for any integer d"""
    if d < 0:
        return 1
    n = 6 * d + 1
    cnt = 0
    v = 2
    while q(v) <= n:
        if n % q(v) == 0:
            cnt += 1
        v += 1
    return cnt


def m_sieve(X):
    """numpy sieve: m(d) for 0 <= d <= X (index d); uses pow(6,-1,q) for the residue"""
    m = np.ones(X + 1, dtype=np.uint8)  # v = 2, q_2 = 1 divides everything
    v = 3
    while q(v) <= 6 * X + 1:
        qq = q(v)
        d0 = (-pow(6, -1, qq)) % qq
        if d0 == 0:
            d0 = qq
        m[d0::qq] += 1
        v += 1
    return m


# ---------------------------------------------------------------------------
section("A. The drop identity and the multiplicity theorem")

# A1: identity 6d+1 = (2^v-3) Syr(A) for all odd A < 4*10^6 (vectorised, int64)
A = np.arange(1, 4 * 10 ** 6, 2, dtype=np.int64)
t = 3 * A + 1
low = t & (-t)
v = np.log2(low.astype(np.float64)).round().astype(np.int64)
assert np.all((np.int64(1) << v) == low)
S = t >> v
d = (A - S) // 2
check(np.all((A - S) % 2 == 0) and np.all(6 * d + 1 == ((np.int64(1) << v) - 3) * S),
      "6d(A)+1 = (2^v-3) Syr(A) for every odd A < 4*10^6")
check(np.all((d < 0) == (A % 4 == 3)) and np.all(d[A % 4 == 3] == -(A[A % 4 == 3] + 1) // 4),
      "d(A) < 0 iff A = 3 mod 4, and then d(A) = -(A+1)/4 (NOT -(A+1)/2)")
del t, low

# A2: multiplicity theorem, direct enumeration of A versus the divisibility formula
D0 = 200000
Amax = 8 * D0 + 1          # proved: d(A) = d >= 0 forces A <= 8d+1; d < 0 forces A = -4d-1
sel = A <= Amax
dd = d[sel]
inr = (dd >= -D0) & (dd <= D0)
cnt = np.bincount((dd[inr] + D0).astype(np.int64), minlength=2 * D0 + 1)
ok = True
for dv in range(-D0, D0 + 1):
    if cnt[dv + D0] != m_formula(dv):
        ok = False
        print("   mismatch at d =", dv, cnt[dv + D0], m_formula(dv))
        break
check(ok, "multiplicity theorem: #{odd A>0 : d(A)=d} = m(d) for every |d| <= %d "
          "(direct orbit count over A <= %d versus divisibility formula)" % (D0, Amax))
# the preimages themselves
ok = True
As_, ds_ = A[sel], d[sel]
order_ = np.argsort(ds_, kind="stable")
ds_sorted, As_sorted = ds_[order_], As_[order_]
for dv in list(range(1, 3000)) + [24, 314, 7854]:
    i0, i1 = np.searchsorted(ds_sorted, dv, "left"), np.searchsorted(ds_sorted, dv, "right")
    pre = sorted(int(a) for a in As_sorted[i0:i1])
    pred = sorted((2 ** (vv + 1) * dv + 1) // q(vv) for vv in range(2, 40)
                  if q(vv) <= 6 * dv + 1 and (6 * dv + 1) % q(vv) == 0)
    if pre != pred:
        ok = False
        break
check(ok, "the preimages are exactly A = (2^(v+1) d + 1)/(2^v - 3) = 2d + (6d+1)/(2^v-3), "
          "one per admissible v (d < 3000 and d = 24, 314, 7854)")
check(all(np.searchsorted(ds_sorted, -e, "right") - np.searchsorted(ds_sorted, -e, "left") == 1 and
          int(As_sorted[np.searchsorted(ds_sorted, -e, "left")]) == 4 * e - 1 for e in range(1, 2000))
      and int(np.sum(ds_ == 0)) == 1,
      "d < 0 occurs once (A = -4d-1), d = 0 once (A = 1)")
del dd, inr, cnt

# A3: owner's data
K_owner = [0, 2, 1, 4, 2, 10, 3, 9, 4, 15, 5]
check([drop(4 * N - 3) for N in range(1, 12)] == K_owner,
      "owner's K = 0,2,1,4,2,10,3,9,4,15,5 is K_N = d(4N-3)")
check([syr(a)[0] for a in (1, 5, 9, 13)] == [1, 1, 7, 5],
      "owner's examples: 1 -> 1, 5 -> 1, 9 -> 7, 13 -> 5")
check(all(2 * N - 1 - K_owner[N - 1] == (syr(4 * N - 3)[0] + 1) // 2 for N in range(1, 12)),
      "owner's rule: label 2N-1 is sent to 2N-1-K_N, the label of Syr(4N-3)")

# ---------------------------------------------------------------------------
section("B. Independent audit of THM-4527, Statement items 1 and 2")


def F_closed(M):            # (oddpart(3M-1)+1)/2
    x = 3 * M - 1
    while x % 2 == 0:
        x //= 2
    return (x + 1) // 2


def F_rules(M):             # three rules: F(2N)=3N, F(4j+1)=3j+1, F(4n-1)=F(n)
    while M % 4 == 3:
        M = (M + 1) // 4
    if M % 2 == 0:
        return 3 * M // 2
    return 3 * ((M - 1) // 4) + 1


ok1 = ok2 = True
for M in range(1, 10 ** 6 + 1):
    s_ = (syr(2 * M - 1)[0] + 1) // 2
    if F_closed(M) != s_:
        ok1 = False
        break
    if F_rules(M) != s_:
        ok2 = False
        break
check(ok1, "item 1: F(M) := (Syr(2M-1)+1)/2 equals (oddpart(3M-1)+1)/2 for M <= 10^6")
check(ok2, "item 1: the three rules F(2N)=3N, F(4j+1)=3j+1, F(4n-1)=F(n) reproduce F for M <= 10^6")
Kt = [None] + [drop(4 * N - 3) for N in range(1, 10 ** 6 + 1)]
ok = all(Kt[2 * j + 1] == j for j in range(0, 499999)) and \
    all(Kt[4 * m] == 5 * m - 1 for m in range(1, 250001)) and \
    all(Kt[4 * m - 2] == Kt[m] + 6 * m - 4 for m in range(1, 250001))
check(ok, "item 1: K_(2j+1) = j, K_(4m) = 5m-1, K_(4m-2) = K_m + 6m - 4 for N <= 10^6")
check(all(F_closed(4 * n - 1) == F_closed(n) for n in range(1, 200001)) and
      all(syr(4 * a + 1)[0] == syr(a)[0] for a in range(1, 400001, 2)),
      "item 1: microcosm F(4n-1) = F(n), i.e. Syr(4A+1) = Syr(A)")
# item 2: histogram on [0, 10^5]: direct orbit count (only A = 1 mod 4 give d >= 0)
Xh = 10 ** 5
A1 = np.arange(1, 8 * Xh + 2, 4, dtype=np.int64)
t1 = 3 * A1 + 1
lw = t1 & (-t1)
vv1 = np.log2(lw.astype(np.float64)).round().astype(np.int64)
d1 = (A1 - (t1 >> vv1)) // 2
c1 = np.bincount(d1[d1 <= Xh], minlength=Xh + 1)
hist_direct = {}
for x in c1:
    hist_direct[int(x)] = hist_direct.get(int(x), 0) + 1
hist_formula = {}
for dv in range(0, Xh + 1):
    k_ = m_formula(dv)
    hist_formula[k_] = hist_formula.get(k_, 0) + 1
info("histogram of m(d), d in [0,10^5]: direct %s, formula %s" % (dict(sorted(hist_direct.items())),
                                                                   dict(sorted(hist_formula.items()))))
check(hist_direct == hist_formula == {1: 69619, 2: 26617, 3: 3549, 4: 210, 5: 6},
      "item 2 data: histogram on [0,10^5] is {1: 69619, 2: 26617, 3: 3549, 4: 210, 5: 6} (THM-4527 confirmed)")
check(int(c1[0]) == 1 and int(c1[1]) == 1 and int(c1[2]) == 2 and int(c1[24]) == 3,
      "item 2: 0 and 1 occur once, 2 twice (13 = 2^4-3), 24 three times (145 = 5*29)")
# |F(M) - M| over labels M != 3 mod 4: 0 once, every e >= 1 exactly twice
E = 50000
cntabs = [0] * (E + 1)
for M in range(1, 4 * E + 2):
    if M % 4 != 3:
        e_ = abs(F_closed(M) - M)
        if e_ <= E:
            cntabs[e_] += 1
check(cntabs[0] == 1 and all(c == 2 for c in cntabs[1:]),
      "item 2(c): over labels M != 3 mod 4, |F(M)-M| takes 0 once and each 1 <= e <= 50000 exactly twice")
check(all(F_closed(4 * n - 1) - (4 * n - 1) == F_closed(n) - n - 3 * n + 1 for n in range(1, 250001)),
      "note Thm 2(c): Delta(4n-1) = Delta(n) - 3n + 1, Delta(M) = F(M) - M")
check(all((A_ % 8 == 5) == (v_ >= 3) for A_, v_ in ((a, syr(a)[1]) for a in range(1, 200001, 4))),
      "deep descents (v >= 3) are exactly A = 5 mod 8, i.e. labels M = 3 mod 4 (the microcosm quarter)")
# mean value: sum_{v>=2} 1/(2^v-3) with a rigorous enclosure
VQ = 200
P = sum(Fraction(1, q(vq)) for vq in range(2, VQ + 1))
lo = P + Fraction(1, 2 ** VQ)
hi = P + Fraction(1, 2 ** VQ) / (1 - Fraction(3, 2 ** (VQ + 1)))
mpmath.mp.dps = 60
info("sum_(v>=2) 1/(2^v-3) in [%s, %s]" % (mpmath.nstr(mpmath.mpf(lo.numerator) / lo.denominator, 45),
                                          mpmath.nstr(mpmath.mpf(hi.numerator) / hi.denominator, 45)))
QONE_LO, QONE_HI = lo, hi
check(Fraction(134367, 100000) <= lo and hi < Fraction(134368, 100000) and hi - lo < Fraction(1, 10 ** 61),
      "item 2: mean multiplicity sum_(v>=2) 1/(2^v-3) = 1.34367... (enclosure width < 1e-61); not 2")
check(abs(Fraction(sum(hist_formula[k] * k for k in hist_formula), Xh + 1) - Fraction(134369, 100000)) < Fraction(1, 10 ** 5),
      "item 2: empirical mean on [0,10^5] is 1.34369 (THM-4527 value confirmed)")
del Kt, A1, t1, lw, vv1, d1, c1

# ---------------------------------------------------------------------------
section("C. Natural densities of {d >= 1 : m(d) = k}: exact truncations + proved tail")


def shared_structure(T, a=3, vmin=3):
    """moduli Q_v = 2^v - a (vmin <= v <= T); primes shared by two moduli; exponents; private parts"""
    Q = {v_: 2 ** v_ - a for v_ in range(vmin, T + 1)}
    vs = list(range(vmin, T + 1))
    shared = {}
    for i, v_ in enumerate(vs):
        for w_ in vs[i + 1:]:
            g = gcd(Q[v_], Q[w_])
            if g > 1:
                for p_ in sympy.factorint(g):
                    shared.setdefault(p_, set()).update([v_, w_])
    e = {v_: {} for v_ in vs}
    u = {}
    for v_ in vs:
        r = Q[v_]
        for p_ in shared:
            k_ = 0
            while r % p_ == 0:
                r //= p_
                k_ += 1
            if k_:
                e[v_][p_] = k_
        u[v_] = r
    amax = {p_: max(e[v_].get(p_, 0) for v_ in vs) for p_ in shared}
    return Q, vs, shared, e, u, amax


def pmul_lin(Pl, c):
    R = [0] * (len(Pl) + 1)
    for i, x in enumerate(Pl):
        R[i] += x * c
        R[i + 1] += x
    return R


def padd(Pl, Ql):
    if len(Pl) < len(Ql):
        Pl, Ql = Ql, Pl
    R = Pl[:]
    for i, x in enumerate(Ql):
        R[i] += x
    return R


def pmul(Pl, Ql):
    R = [0] * (len(Pl) + len(Ql) - 1)
    for i, x in enumerate(Pl):
        if x:
            for j, y in enumerate(Ql):
                R[i + j] += x * y
    return R


def exact_pgf(T, a=3, vmin=3):
    """E[t^N_T] = numer(t)/Dn exactly, N_T(n) = #{vmin <= v <= T : (2^v - a) | n}, n uniform mod lcm.
    Independence: by CRT the levels min(v_p(n), a_p) of the shared primes and the events
    [u_v | n] for the pairwise coprime private parts u_v are independent."""
    Q, vs, shared, e, u, amax = shared_structure(T, a, vmin)
    touched = [v_ for v_ in vs if e[v_]]
    untouched = [v_ for v_ in vs if not e[v_]]
    primes = sorted(shared)
    parent = {p_: p_ for p_ in primes}

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for v_ in touched:
        ps = list(e[v_])
        for p_ in ps[1:]:
            parent[find(p_)] = find(ps[0])
    comps = {}
    for p_ in primes:
        comps.setdefault(find(p_), []).append(p_)
    Dn = 1
    for p_ in primes:
        Dn *= p_ ** amax[p_]
    for v_ in vs:
        Dn *= u[v_]
    total = [1]
    for v_ in untouched:
        total = pmul_lin(total, u[v_] - 1)
    ncfg = 0
    for root, ps in comps.items():
        cv = [v_ for v_ in touched if find(next(iter(e[v_]))) == root]
        comp_poly = [0]
        for cfg in itertools.product(*[range(amax[p_] + 1) for p_ in ps]):
            ncfg += 1
            wgt = 1
            for p_, l_ in zip(ps, cfg):
                wgt *= (p_ - 1) * p_ ** (amax[p_] - l_ - 1) if l_ < amax[p_] else 1
            lev = dict(zip(ps, cfg))
            Pl = [wgt]
            for v_ in cv:
                if all(lev[p_] >= k_ for p_, k_ in e[v_].items()):
                    Pl = pmul_lin(Pl, u[v_] - 1)
                else:
                    Pl = [x * u[v_] for x in Pl]
            comp_poly = padd(comp_poly, Pl)
        total = pmul(total, comp_poly)
    L = 1
    for v_ in vs:
        L = lcm(L, Q[v_])
    return total, Dn, ncfg, shared, L


TD = 100
numer, Dn, ncfg, shared100, L100 = exact_pgf(TD)
check(Dn == L100 and sum(numer) == Dn, "exact PGF at T = %d: common denominator = lcm(q_3..q_100) "
      "(%d digits), coefficients sum to 1; %d shared primes, %d level configurations"
      % (TD, len(str(Dn)), len(shared100), ncfg))
info("shared primes up to T=100 (prime: v's): " + ", ".join(
    "%d:%s" % (p_, sorted(shared100[p_])) for p_ in sorted(shared100)))
distT = [Fraction(x, Dn) for x in numer]
check(sum(i * x for i, x in enumerate(distT)) == sum(Fraction(1, q(v_)) for v_ in range(3, TD + 1)),
      "exact PGF at T = 100: mean of N_T equals sum_(v=3..100) 1/q_v exactly")
epsT = Fraction(1, 2 ** TD) / (1 - Fraction(3, 2 ** (TD + 1)))   # >= sum_{v>T} 1/q_v
check(sum(Fraction(1, q(v_)) for v_ in range(TD + 1, TD + 400)) < epsT,
      "tail bound eps_T = 2^-T/(1-3*2^-(T+1)) dominates sum_(v>T) 1/q_v (T = 100, ~7.9e-31)")
mpmath.mp.dps = 50
DELTA = {}
for k_ in range(1, 12):
    x = distT[k_ - 1] if k_ - 1 < len(distT) else Fraction(0)
    DELTA[k_] = x
    info("delta_%d in [%s, %s]" % (k_, mpmath.nstr(mpmath.mpf(x.numerator) / x.denominator - mpmath.mpf(epsT.numerator) / epsT.denominator, 28),
                                   mpmath.nstr(mpmath.mpf(x.numerator) / x.denominator + mpmath.mpf(epsT.numerator) / epsT.denominator, 28)))
check(abs(DELTA[1] - Fraction(6961749034401326, 10 ** 16)) < Fraction(1, 10 ** 15) and
      abs(DELTA[2] - Fraction(2662380366777384, 10 ** 16)) < Fraction(1, 10 ** 15) and
      abs(DELTA[3] - Fraction(3538955883931588, 10 ** 17)) < Fraction(1, 10 ** 16),
      "delta_1 = 0.6961749034401326..., delta_2 = 0.2662380366777384..., delta_3 = 0.0353895588393158...")
# truncation stability (the limit exists; difference between T=64 and T=100 below eps_64)
numer64, Dn64, _, _, _ = exact_pgf(64)
eps64 = Fraction(1, 2 ** 64) / (1 - Fraction(3, 2 ** 65))
check(all(abs(Fraction(numer64[i] if i < len(numer64) else 0, Dn64) - distT[i]) <= eps64 for i in range(len(distT))),
      "|delta_k^(64) - delta_k^(100)| <= eps_64 for every k (as the proof predicts)")

# C2: brute force over a full period for T = 7 (independent path: residues of d)
L7 = 1
for v_ in range(3, 8):
    L7 = lcm(L7, q(v_))
cnt7 = [0] * 6
for dv in range(L7):
    n_ = 6 * dv + 1
    cnt7[sum(1 for v_ in range(3, 8) if n_ % q(v_) == 0)] += 1
num7, D7, _, _, _ = exact_pgf(7)
check(L7 == 2874625 and all(Fraction(cnt7[i], L7) == Fraction(num7[i] if i < len(num7) else 0, D7) for i in range(6)),
      "T = 7: exact PGF equals the brute-force count over a full period d mod 2874625 = lcm(5,13,29,61,125) = 6*479104+1")
info("T=7 period counts: %s (note 125 = 5^3: X_7 = 1 forces X_3 = 1)" % cnt7)

# C3: inclusion-exclusion over all subsets for T = 16 (independent path: lcm of subsets)
T16 = 16
vs16 = list(range(3, T16 + 1))
sig = [Fraction(0)] * (len(vs16) + 1)


def rec(i, Lc, size):
    if i == len(vs16):
        sig[size] += Fraction(1, Lc)
        return
    rec(i + 1, Lc, size)
    rec(i + 1, lcm(Lc, q(vs16[i])), size + 1)


rec(0, 1, 0)
num16, D16, _, _, _ = exact_pgf(T16)
ok = True
for j in range(len(vs16) + 1):
    ie = sum((-1) ** (i - j) * comb(i, j) * sig[i] for i in range(j, len(vs16) + 1))
    if ie != Fraction(num16[j] if j < len(num16) else 0, D16):
        ok = False
check(ok, "T = 16: delta_k^(T) = sum_i (-1)^(i-k+1) C(i,k-1) sigma_i^(T), sigma_i = sum_{|S|=i} 1/lcm(q_S) "
          "(all 2^14 subsets) equals the exact PGF")

# C4: empirical frequencies up to 10^8 (sieve)
XS = 10 ** 8
msv = m_sieve(XS)
hist8 = np.bincount(msv[1:])
info("histogram of m(d), 1 <= d <= 10^8: %s" % {k_: int(c) for k_, c in enumerate(hist8) if c})
check(all(abs(int(hist8[k_]) / XS - float(DELTA[k_])) < 2e-6 for k_ in range(1, len(hist8))),
      "EMPIRICAL: frequencies of m(d) = k on [1,10^8] agree with delta_k to within 2e-6")
check(np.array_equal(msv[:D0 + 1].astype(np.int64), np.array([m_formula(x) for x in range(0, D0 + 1)])),
      "sieve m(d) agrees with the divisibility formula for 0 <= d <= 200000")
# C5: proved tail bound P(m >= k+1) <= (k + 5/4)/(2*3^k), k >= 2
tail_ok = True
for k_ in range(2, 10):
    Pge = sum(distT[k_:])            # P(N_T >= k) = P(m_T >= k+1)
    if not (Pge <= Fraction(4 * k_ + 5, 8 * 3 ** k_)):
        tail_ok = False
    # closed form of the bound series sum_{w>=k+1} C(w-3,k-2)(w+1)2^(1-2w)
    ser = sum(Fraction(comb(w_ - 3, k_ - 2) * (w_ + 1), 2 ** (2 * w_ - 1)) for w_ in range(k_ + 1, 400))
    if not (ser <= Fraction(4 * k_ + 5, 8 * 3 ** k_) and Fraction(4 * k_ + 5, 8 * 3 ** k_) - ser < Fraction(1, 10 ** 50)):
        tail_ok = False
check(tail_ok, "tail lemma: sum_{w>=k+1} C(w-3,k-2)(w+1)2^(1-2w) = (k+5/4)/(2*3^k), and the exact "
               "P(m >= k+1) at T=100 respects it (k = 2..9)")
info("P(m >= k), k=2..9: " + ", ".join("%.3e" % float(sum(distT[k_ - 1:])) for k_ in range(2, 10)))
# C6: the note's table: dens{m >= k} to 5 significant digits (both ends of the proved enclosure)
TAILSTR = {2: "0.30383", 3: "0.037587", 4: "2.1975e-3", 5: "6.2916e-5", 6: "8.5447e-7", 7: "5.3196e-9",
           8: "1.6829e-11", 9: "2.8693e-14", 10: "2.6987e-17"}
okt = True
for k_, sval in TAILSTR.items():
    x = sum(distT[k_ - 1:])          # P_T(1 + N_T >= k); |dens{m >= k} - x| <= eps_T as in Theorem 3.4(1)
    okt = okt and r5(x - epsT) == float(sval) == r5(x + epsT)
    info("dens{m >= %d} = %s" % (k_, mpmath.nstr(mpmath.mpf(x.numerator) / x.denominator, 12)))
check(okt, "note table: dens{m >= k} = 0.30383, 0.037587, 2.1975e-3, 6.2916e-5, 8.5447e-7, 5.3196e-9, 1.6829e-11, "
           "2.8693e-14, 2.6987e-17 (k = 2..10; both ends of the enclosure round to these)")
DSTR = {1: "0.6961749034401325952680065255", 2: "0.2662380366777383738915135713",
        3: "0.03538955883931588244282135104", 4: "0.00213458515130198837708089211",
        5: "6.20614225702492134536025e-5", 6: "8.4914929170526775682603e-7", 7: "5.3028197674900875448e-9",
        8: "1.680074489300358645e-11", 9: "2.866616892777891e-14", 10: "2.6972975032442e-17",
        11: "1.436873267e-20"}
okd = True
for k_, sval in DSTR.items():
    mant, _, ex = sval.partition("e")
    last_place = Fraction(10) ** (int(ex or 0) - len(mant.split(".")[1]))
    okd = okd and abs(Fraction(sval) - DELTA[k_]) + epsT <= last_place
check(okd, "note table: every printed digit string of delta_1..delta_11 is within one unit of its last place "
           "of the true delta_k (|string - delta_k^(100)| + eps_T <= one unit)")
# C7: the independent-events model (X_v independent Bernoulli(1/q_v), v = 3..T; coupling error <= eps_T)
indep = [Fraction(1)]
for v_ in range(3, TD + 1):
    xq = Fraction(1, q(v_))
    nw = [Fraction(0)] * (len(indep) + 1)
    for i, c in enumerate(indep):
        nw[i] += c * (1 - xq)
        nw[i + 1] += c * xq
    indep = nw
INDSTR = {1: "0.69023", 2: "0.27725", 3: "0.031150", 4: "1.3412e-3", 5: "2.5082e-5", 6: "2.1688e-7",
          7: "8.9705e-10", 8: "1.8091e-12"}
check(all(r5(indep[k_ - 1] - epsT) == float(sv) == r5(indep[k_ - 1] + epsT) for k_, sv in INDSTR.items()),
      "independent-events model P(1 + sum_v X_v = k) = 0.69023, 0.27725, 0.031150, 1.3412e-3, 2.5082e-5, "
      "2.1688e-7, 8.9705e-10, 1.8091e-12 (k = 1..8); true/model ratio at k = 5 is %.2f" %
      float(DELTA[5] / indep[4]))
del indep

# ---------------------------------------------------------------------------
section("D. Generating function and Dirichlet series")
# D1: sum_{d>=1} m(d) x^d = x/(1-x) + sum_{v even >= 4} x^((2^v-4)/6)/(1-x^(2^v-3))
#                                   + sum_{v odd >= 3} x^((5*2^v-16)/6)/(1-x^(2^v-3))
XG = 200000
coef = np.zeros(XG + 1, dtype=np.int64)
coef[1:] += 1
v_ = 3
while q(v_) <= 6 * XG + 1:        # note: odd v has a larger first exponent than v+1, so no early break
    first = (2 ** v_ - 4) // 6 if v_ % 2 == 0 else (5 * 2 ** v_ - 16) // 6
    assert (6 * first + 1) % q(v_) == 0 and (6 * first + 1) // q(v_) in (1, 5)
    coef[first::q(v_)] += 1
    v_ += 1
check(np.array_equal(coef[1:], msv[1:XG + 1].astype(np.int64)),
      "generating function: the geometric-series formula reproduces m(d) for 1 <= d <= 200000 "
      "(first exponent (2^v-4)/6 for even v, (5*2^v-16)/6 for odd v)")
check(all(((6 * ((2 ** w - 4) // 6) + 1) == q(w)) for w in range(4, 200, 2)) and
      all(((6 * ((5 * 2 ** w - 16) // 6) + 1) == 5 * q(w)) for w in range(3, 200, 2)),
      "2^v - 3 = 1 (mod 6) for even v and = 5 (mod 6) for odd v, so 6d+1 = q_v * s with s = 1 or 5 (mod 6)")
# D2: Dirichlet series F(s) = sum_{n = 1 mod 6} m'(n) n^-s, m'(n) = #{v>=2 : q_v | n}
mpmath.mp.dps = 30


def Qs(s, sign=1):
    return mpmath.nsum(lambda vv: (sign ** int(vv)) * (2 ** vv - 3) ** (-s), [2, mpmath.inf])


def Fclosed(s):
    L0 = mpmath.zeta(s) * (1 - mpmath.mpf(2) ** (-s)) * (1 - mpmath.mpf(3) ** (-s))
    Lchi3 = mpmath.mpf(3) ** (-s) * (mpmath.zeta(s, mpmath.mpf(1) / 3) - mpmath.zeta(s, mpmath.mpf(2) / 3))
    Lchi = (1 + mpmath.mpf(2) ** (-s)) * Lchi3
    return (L0 * Qs(s) + Lchi * Qs(s, -1)) / 2


CH = 10 ** 7
for s_ in (2, 3):
    part = 0.0
    for c0 in range(0, XS + 1, CH):
        c1_ = min(XS + 1, c0 + CH)
        dg = np.arange(c0, c1_, dtype=np.float64)
        part += float(np.sum(msv[c0:c1_].astype(np.float64) / (6 * dg + 1) ** s_))
    del dg
    # rigorous tail: m(d) <= log2(6d+4) - 1 < log2(7d) for d > 10^8, and sum_{d>X} log2(7d)/(6d)^s <= ...
    tail = float(mpmath.quad(lambda x: mpmath.log(7 * x, 2) / (6 * x) ** s_, [XS, mpmath.inf]) +
                 mpmath.log(7 * XS, 2) / (6 * XS) ** s_)
    fc = Fclosed(s_)
    info("s=%d: partial sum (d<=1e8) %.15f, tail <= %.2e, closed form %s" % (s_, part, tail, mpmath.nstr(fc, 18)))
    check(-1e-12 <= float(fc) - part <= tail + 1e-12,
          "Dirichlet series identity F(s) = (L(s,chi0)Q(s) + L(s,chi)Q^-(s))/2 holds numerically at s = %d" % s_)


def Q_cont(s, sign=1, J=200):
    s = mpmath.mpmathify(s)
    tot = mpmath.mpf(1)
    for j in range(J):
        coefj = mpmath.rf(s, j) / mpmath.factorial(j) * mpmath.mpf(3) ** j
        g = mpmath.mpf(2) ** (-(s + j))
        if sign == 1:
            tot += coefj * mpmath.mpf(2) ** (-3 * (s + j)) / (1 - g)
        else:
            tot += coefj * (-mpmath.mpf(2) ** (-3 * (s + j))) / (1 + g)
    return tot


check(abs(Q_cont(2) - Qs(2)) < mpmath.mpf(10) ** -25 and abs(Q_cont(1) - Qs(1)) < mpmath.mpf(10) ** -25
      and abs(Q_cont(2, -1) - Qs(2, -1)) < mpmath.mpf(10) ** -25,
      "continuation formula Q(s) = 1 + sum_j (s)_j/j! 3^j 2^(-3(s+j))/(1-2^(-(s+j))) matches the series at s = 1, 2 (and Q^- at s = 2)")
check(abs(Q_cont(0, -1) - mpmath.mpf(1) / 2) < mpmath.mpf(10) ** -25 and
      abs(mpmath.mpf(10) ** -12 * Q_cont(mpmath.mpf(10) ** -12) - 1 / mpmath.log(2)) < mpmath.mpf(10) ** -9 and
      abs((-1 + mpmath.mpf(10) ** -12 - (-1)) * Q_cont(-1 + mpmath.mpf(10) ** -12) - (-3) / mpmath.log(2)) < mpmath.mpf(10) ** -8,
      "Q^-(0) = 1/2; Q has residue 1/ln2 at s = 0 and (-3)^1/ln2 at s = -1 (numerically)")


def Fcont(s):
    s = mpmath.mpmathify(s)
    L0 = mpmath.zeta(s) * (1 - mpmath.mpf(2) ** (-s)) * (1 - mpmath.mpf(3) ** (-s))
    Lchi3 = mpmath.mpf(3) ** (-s) * (mpmath.zeta(s, mpmath.mpf(1) / 3) - mpmath.zeta(s, mpmath.mpf(2) / 3))
    return (L0 * Q_cont(s) + (1 + mpmath.mpf(2) ** (-s)) * Lchi3 * Q_cont(s, -1)) / 2


eps_s = mpmath.mpf(10) ** -15
check(abs(Fcont(2) - Fclosed(2)) < mpmath.mpf(10) ** -25 and abs(Fcont(eps_s) - mpmath.mpf(1) / 6) < mpmath.mpf(10) ** -12
      and abs(Fcont(-eps_s) - mpmath.mpf(1) / 6) < mpmath.mpf(10) ** -12
      and abs(3 ** 0 * (mpmath.zeta(0, mpmath.mpf(1) / 3) - mpmath.zeta(0, mpmath.mpf(2) / 3)) - mpmath.mpf(1) / 3) < mpmath.mpf(10) ** -25
      and abs(eps_s * Fcont(1 + eps_s) - Qs(1) / 6) < mpmath.mpf(10) ** -12,
      "F(0) = 1/6 (F(+-1e-15) through the continuation; L(0,chi_-3) = 1/3), and (s-1)F(s) -> Q(1)/6 at s = 1")
info("Q(1) = %s ; Q^-(1) = %s ; residue of F at s=1 is Q(1)/6 = %s" % (
    mpmath.nstr(Qs(1), 20), mpmath.nstr(Qs(1, -1), 20), mpmath.nstr(Qs(1) / 6, 20)))
# summatory function: sum_{d<=X} m(d) = sum_v floor-counts; check |sum - X Q(1)| <= log2(6X+4)
okS = True
Q1f = float(Qs(1))
for X_ in [10, 100, 1000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, 10 ** 8]:
    exact_sum = 0
    w = 2
    while q(w) <= 6 * X_ + 1:
        first = 1 if w == 2 else ((2 ** w - 4) // 6 if w % 2 == 0 else (5 * 2 ** w - 16) // 6)
        if first <= X_:
            exact_sum += (X_ - first) // q(w) + 1
        w += 1
    if exact_sum != int(np.sum(msv[1:X_ + 1], dtype=np.int64)) or abs(exact_sum - X_ * Q1f) > math.log2(6 * X_ + 4) + 1:
        okS = False
    info("X=%d: sum_(d<=X) m(d) = %d, X*Q(1) = %.3f" % (X_, exact_sum, X_ * Q1f))
check(okS, "summatory function: sum_(d<=X) m(d) = sum_v (floor((X - e_v)/q_v) + 1) and |sum - X Q(1)| <= log2(6X+4)+1")

# ---------------------------------------------------------------------------
section("E. Maximal order of m(d)")
ok = all(gcd(q(a_), q(b_)) == gcd(q(a_), 2 ** (b_ - a_) - 1) for a_ in range(2, 151) for b_ in range(a_ + 1, 151))
check(ok, "gcd(q_v, q_w) = gcd(q_v, 2^(w-v) - 1) for 2 <= v < w <= 150")
ok = all(lcm(q(a_), q(b_)) >= (2 ** a_ + 1) * q(a_) for a_ in range(2, 151) for b_ in range(a_ + 1, 151))
check(ok, "gap principle: lcm(q_v, q_w) >= (2^v + 1) q_v for 2 <= v < w <= 150")
okb = True
for c0 in range(1, XS + 1, CH):
    c1_ = min(XS + 1, c0 + CH)
    mm = msv[c0:c1_].astype(np.int64)
    dd8 = np.arange(c0, c1_, dtype=np.int64)
    okb = okb and bool(np.all(((np.int64(1) << mm) - 1) ** 2 <= 6 * dd8 + 5))
del mm, dd8
check(okb, "upper bound m(d) <= log2(1 + sqrt(6d+5)) holds for every 1 <= d <= 10^8")
# (1 + sqrt(6d+5))^2 / d = 6 + 6/d + 2 sqrt(6d+5)/d is decreasing, = 12 + 2 sqrt(11) < 18.7 at d = 1
gvals = [(1 + math.sqrt(6 * dv + 5)) ** 2 / dv for dv in range(1, 200001)]
check(all(gvals[i] > gvals[i + 1] for i in range(len(gvals) - 1)) and gvals[0] < 18.7 and
      0.5 * math.log2(18.7) < 2.2 and int(msv[1:].max()) <= 0.5 * math.log2(XS) + 2.2,
      "m(d) <= (1/2) log2 d + 2.2: (1+sqrt(6d+5))^2/d decreases from 12 + 2 sqrt(11) = %.3f < 18.7 "
      "(d <= 2*10^5 checked; monotone in general), and log2(18.7)/2 = %.4f" % (gvals[0], 0.5 * math.log2(18.7)))
del gvals
first_sieve = {}
for k_ in range(1, int(msv[1:].max()) + 1):
    first_sieve[k_] = int(np.argmax(msv[1:] >= k_)) + 1


def dmin_of_L(Lc):
    for t_ in range(1, 13):
        if (t_ * Lc - 1) % 6 == 0 and t_ * Lc > 1:
            return (t_ * Lc - 1) // 6


def first_d_with_m_at_least(k_):
    j = k_ - 1
    if j == 0:
        return 1, []
    L0 = 1
    for w in range(3, 3 + j):
        L0 = lcm(L0, q(w))
    best = [dmin_of_L(L0), list(range(3, 3 + j))]
    vmax = (6 * best[0] + 4).bit_length() + 1

    def dfs(start, Lc, chosen):
        if len(chosen) == j:
            dv = dmin_of_L(Lc)
            if dv < best[0] or (dv == best[0] and chosen < best[1]):
                best[0], best[1] = dv, chosen[:]
            return
        for w in range(start, vmax + 1):
            if q(w) > 6 * best[0] + 1:
                break
            L2 = lcm(Lc, q(w))
            if L2 > 6 * best[0] + 1:
                continue
            chosen.append(w)
            dfs(w + 1, L2, chosen)
            chosen.pop()
    dfs(3, 1, [])
    return best[0], best[1]


RECORDS = {}
for k_ in range(1, 12):
    RECORDS[k_] = first_d_with_m_at_least(k_)
    LS = 1
    for w in RECORDS[k_][1]:
        LS = lcm(LS, q(w))
    info("smallest d with m(d) >= %d: d = %d; S = {v >= 3 : q_v | 6d+1} = %s; 6d+1 = %d * lcm(q_S)" % (
        k_, RECORDS[k_][0], RECORDS[k_][1], (6 * RECORDS[k_][0] + 1) // LS))
check(all(RECORDS[k_][0] == first_sieve[k_] for k_ in first_sieve),
      "records k <= 6: exhaustive lcm search agrees with the 10^8 sieve: d = 1, 2, 24, 314, 7854, 479104")
check([RECORDS[k_][0] for k_ in range(7, 12)] == [121213354, 49576261854, 50617363353104,
                                                  115043252496711229, 117459160799142164979],
      "records k = 7..11 (exhaustive lcm search): 121213354, 49576261854, 50617363353104, "
      "115043252496711229, 117459160799142164979")
check(all(m_formula(RECORDS[k_][0]) == k_ for k_ in range(1, 12)),
      "each record d has m(d) exactly k (direct divisibility test over all v)")
okc = True
for k_ in range(2, 41):
    Lc = 1
    for w in range(3, k_ + 2):
        Lc = lcm(Lc, q(w))
    dk = dmin_of_L(Lc)
    if not (m_formula(dk) >= k_ and dk < 2 ** ((k_ + 1) * (k_ + 2) // 2 - 3)):
        okc = False
check(okc, "CRT construction: for 2 <= k <= 40, d_k = least d with lcm(q_3..q_(k+1)) | 6d+1 has m(d_k) >= k "
           "and d_k < 2^((k+1)(k+2)/2 - 3), hence max_(d<=X) m(d) >= sqrt(2 log2 X) - 3/2 - o(1)")
for k_ in range(1, 12):
    info("k=%2d: log2(record) = %7.3f, sqrt(2 log2 record) = %6.3f, (log2 record)/2 = %6.3f" % (
        k_, math.log2(RECORDS[k_][0]) if RECORDS[k_][0] > 1 else 0.0,
        math.sqrt(2 * math.log2(RECORDS[k_][0])) if RECORDS[k_][0] > 1 else 0.0,
        math.log2(RECORDS[k_][0]) / 2 if RECORDS[k_][0] > 1 else 0.0))
# AP lemma for prime powers Q | q_v, v <= 64
okap = True
nQ = 0
distinctQ = set()
for w in range(3, 65):
    for p_, ex in sympy.factorint(q(w)).items():
        for e_ in range(1, ex + 1):
            Qp = p_ ** e_
            o_ = sympy.n_order(2, Qp)
            hits = [u_ for u_ in range(2, 400) if q(u_) % Qp == 0]
            r_ = hits[0]
            nQ += 1
            distinctQ.add(Qp)
            if not (r_ < o_ and hits == list(range(r_, 400, o_)) and (3 * pow(2, o_ - r_, Qp) - 1) % Qp == 0):
                okap = False
check(okap, "AP lemma: for each of the %d distinct prime powers Q dividing some q_v, v <= 64 (%d pairs (v, Q)), "
            "{u : Q | q_u} = r + o*N with o = ord_Q(2) > r >= 3, and Q | 3*2^(o-r) - 1" % (len(distinctQ), nQ))
# data for the OPEN problem: primitive parts
Lrun = 1
worst = (10.0, None)
for w in range(3, 201):
    g = gcd(q(w), Lrun)
    Pw = q(w) // g
    if w >= 10:
        worst = min(worst, (math.log2(Pw) / w, w))
    Lrun = lcm(Lrun, q(w))
ratios = []
for V_ in (20, 50, 100, 200):
    Lc = 1
    sl = 0
    prodq = 1
    for w in range(3, V_ + 1):
        Lc = lcm(Lc, q(w))
        prodq *= q(w)
    loss = math.log2(prodq) - math.log2(Lc)
    ratios.append((V_, round(math.log2(Lc) / math.log2(prodq), 4), round(loss, 1), round(loss / (V_ * math.log2(V_)), 3)))
info("EMPIRICAL: min_(10<=v<=200) log2(P_v)/v = %.4f at v = %d (P_v = primitive part of q_v)" % worst)
info("EMPIRICAL: (V, log2 lcm(q_3..q_V)/log2 prod q_v, loss = log2(prod/lcm), loss/(V log2 V)) = %s" % ratios)
check(worst[0] > 0.5 and all(r_[3] < 0.5 for r_ in ratios),
      "EMPIRICAL (data for the OPEN sqrt(log) bound): primitive parts P_v > 2^(v/2) for 10 <= v <= 200, "
      "and the lcm loss for initial segments is < 0.5 V log2 V for V = 20, 50, 100, 200")

# ---------------------------------------------------------------------------
section("F. The label map ('scaled down to N'): N = (A+1)/2, F(N) = (Syr(2N-1)+1)/2")
NF = 10 ** 6
Fv = [0] * (NF + 1)
okf = True
for N_ in range(1, NF + 1):
    x = 3 * N_ - 1
    w = v2(x)
    f_ = (x + 2 ** w) // 2 ** (w + 1)
    Fv[N_] = f_
    if N_ % 2 == 0:
        okf = okf and f_ == 3 * N_ // 2
    elif N_ % 4 == 1:
        okf = okf and f_ == (3 * N_ + 1) // 4
    else:
        okf = okf and f_ == Fv[(N_ + 1) // 4]
    okf = okf and (N_ - f_) * 2 ** (w + 1) == (2 ** (w + 1) - 3) * N_ - (2 ** w - 1)
check(okf, "F(N) = (3N-1+2^w)/2^(w+1), w = v2(3N-1): even N -> 3N/2, N = 1 mod 4 -> (3N+1)/4, "
           "N = 3 mod 4 -> F((N+1)/4); N - F(N) = ((2^(w+1)-3)N - (2^w-1))/2^(w+1)  (N <= 10^6)")
pre = {}
for N_ in range(1, NF + 1):
    if Fv[N_] <= 3000:
        pre.setdefault(Fv[N_], []).append(N_)
okp = True
for n_ in range(1, 3001):
    if n_ % 3 == 2:
        pred = []
    else:
        N0 = 2 * n_ // 3 if n_ % 3 == 0 else (4 * n_ - 1) // 3
        pred = []
        j = 0
        while 4 ** j * N0 - (4 ** j - 1) // 3 <= NF:
            pred.append(4 ** j * N0 - (4 ** j - 1) // 3)
            j += 1
    okp = okp and sorted(pre.get(n_, [])) == pred
check(okp, "preimages: label n = 2 (mod 3) has none; otherwise exactly N_j = 4^j N_0 - (4^j-1)/3 (j >= 0) with "
           "N_0 = 2n/3 (n = 0 mod 3, an even label) or (4n-1)/3 (n = 1 mod 3, a label = 1 mod 4)  (n <= 3000, N <= 10^6)")
EE = 10000
cpos = [0] * (EE + 1)
cneg = [0] * (EE + 1)
for N_ in range(1, 4 * EE + 2):
    df = N_ - Fv[N_]
    if 0 <= df <= EE:
        cpos[df] += 1
    elif -EE <= df < 0:
        cneg[-df] += 1
check(all(cpos[e_] == m_formula(e_) for e_ in range(0, EE + 1)) and all(cneg[e_] == 1 for e_ in range(1, EE + 1)),
      "difference multiset of the label map: N - F(N) = e occurs m(e) times (e >= 0; m(0) = 1), = -e exactly once")
tot = [cpos[e_] + (cneg[e_] if e_ else 0) for e_ in range(EE + 1)]
first_bad = next(e_ for e_ in range(1, EE + 1) if tot[e_] != 2)
check(first_bad == 2 and sorted(N_ for N_ in range(1, 20) if abs(N_ - Fv[N_]) == 2) == [3, 4, 9],
      "'exactly two copies of each |difference|' first fails at e = 2: |N - F(N)| = 2 for N = 3, 4, 9 (three copies)")
frac_two = sum(1 for e_ in range(1, EE + 1) if tot[e_] == 2) / EE
info("fraction of 1 <= e <= 10^4 with exactly two copies of |difference|: %.4f (density delta_1 = %.6f)" % (
    frac_two, float(DELTA[1])))
kcounts = [cpos[e_] for e_ in range(0, EE + 1)]
check(kcounts[0] == 1 and kcounts[1] == 1 and kcounts[2] == 2 and min(e_ for e_ in range(EE + 1) if kcounts[e_] >= 3) == 24,
      "in the owner's K (odd labels only): 0 and 1 occur once, 2 is the first value occurring twice, 24 the first occurring three times")
# summatory function of K: S(X) = sum_{N<=X} K_N
Ssum = [0] * (10 ** 5 + 1)
for N_ in range(1, 10 ** 5 + 1):
    Ssum[N_] = Ssum[N_ - 1] + drop(4 * N_ - 3)


def S_rec(X_):
    if X_ <= 0:
        return 0
    J = (X_ - 1) // 2
    M1 = X_ // 4
    M2 = (X_ + 2) // 4
    return J * (J + 1) // 2 + 5 * M1 * (M1 + 1) // 2 - M1 + 3 * M2 * (M2 + 1) - 4 * M2 + S_rec(M2)


check(all(S_rec(X_) == Ssum[X_] for X_ in range(0, 10 ** 5 + 1)),
      "summatory function of K: S(X) = J(J+1)/2 + 5M1(M1+1)/2 - M1 + 3M2(M2+1) - 4M2 + S(M2), "
      "J = floor((X-1)/2), M1 = floor(X/4), M2 = floor((X+2)/4)  (X <= 10^5)")
rr = [(S_rec(X_) - X_ * X_ / 2) / X_ for X_ in range(1, 10 ** 6 + 1, 997)]
info("EMPIRICAL: (S(X) - X^2/2)/X ranges in [%.4f, %.4f] for X <= 10^6 (sampled)" % (min(rr), max(rr)))
check(all(abs(x) <= 2 for x in rr), "S(X) = X^2/2 + O(X): |S(X) - X^2/2| <= 2X on the sample (the proof gives O(X))")
# 2-regularity: a(n) = K_(n+1), b(n) = K_(2n+2)
Kl = [None] + [drop(4 * N_ - 3) for N_ in range(1, 4 * 10 ** 5 + 8)]
aK = lambda n_: Kl[n_ + 1]
bK = lambda n_: Kl[2 * n_ + 2]
check(all(aK(2 * n_) == n_ and aK(2 * n_ + 1) == bK(n_) and bK(2 * n_) == aK(n_) + 6 * n_ + 2 and bK(2 * n_ + 1) == 5 * n_ + 4
          for n_ in range(0, 10 ** 5)),
      "K is 2-regular: a(n) = K_(n+1), b(n) = K_(2n+2) satisfy a(2n) = n, a(2n+1) = b(n), b(2n) = a(n) + 6n + 2, "
      "b(2n+1) = 5n + 4  (n < 10^5)")
del Kl
del Fv, pre, Ssum

# ---------------------------------------------------------------------------
section("G. k-step drops: Syr^k(A) = A - 2D,  A(2^p - 3^k) = c(w) + 2^(p+1) D")


def orbit_word(A_, k_):
    w = []
    x = A_
    for _ in range(k_):
        x, vv_ = syr(x)
        w.append(vv_)
    return x, tuple(w)


def cnum(w):
    """c(w) = sum_{i=1..k} 3^(k-i) 2^(v_1+...+v_(i-1))"""
    k_ = len(w)
    c_ = 0
    pp = 0
    for i, vv_ in enumerate(w):
        c_ += 3 ** (k_ - 1 - i) * 2 ** pp
        pp += vv_
    return c_


okg = True
for A_ in range(1, 200001, 2):
    for k_ in (1, 2, 3, 5, 8):
        sk, w = orbit_word(A_, k_)
        p_ = sum(w)
        D_ = (A_ - sk) // 2
        c_ = cnum(w)
        if not ((A_ - sk) % 2 == 0 and A_ * (2 ** p_ - 3 ** k_) == c_ + 2 ** (p_ + 1) * D_ and
                2 * 3 ** k_ * D_ + c_ == (2 ** p_ - 3 ** k_) * sk and c_ <= 2 ** (p_ - k_) * (3 ** k_ - 2 ** k_)):
            okg = False
check(okg, "k-step identities A(2^p-3^k) = c(w) + 2^(p+1) D and 2*3^k D + c(w) = (2^p-3^k) Syr^k(A), "
           "and c(w) <= 2^(p-k)(3^k - 2^k)  (odd A < 2*10^5, k = 1,2,3,5,8)")


def compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first_ in range(1, total - parts + 2):
        for rest in compositions(total - first_, parts - 1):
            yield (first_,) + rest


def p0(k_):
    return next(pp for pp in range(1, 10 * k_ + 10) if 2 ** pp > 3 ** k_)


def Amax_bound(k_, D_):
    P0 = p0(k_)
    r0 = Fraction(2 ** P0, 2 ** P0 - 3 ** k_)
    if D_ >= 0:
        return (Fraction(3 ** k_, 2 ** k_) - 1 + 2 * D_) * r0
    return max(Fraction(2 * 3 ** k_ * (-D_)), (Fraction(3 ** k_, 2 ** k_) - 1) * r0)


DLO, DHI = -30, 300
DIRECT = {}
okm = True
for k_ in (2, 3, 4):
    Abig = int(max(Amax_bound(k_, DLO), Amax_bound(k_, DHI))) + 2
    direct = {}
    pmax_seen = 0
    for A_ in range(1, Abig + 1, 2):
        sk, w = orbit_word(A_, k_)
        D_ = (A_ - sk) // 2
        pmax_seen = max(pmax_seen, sum(w)) if DLO <= D_ <= DHI else pmax_seen
        if DLO <= D_ <= DHI:
            direct.setdefault(D_, set()).add((A_, w))
    PM = 40
    assert pmax_seen < PM
    byword = {}
    for p_ in range(k_, PM + 1):
        Dl = 2 ** p_ - 3 ** k_
        for w in compositions(p_, k_):
            c_ = cnum(w)
            for D_ in range(DLO, DHI + 1):
                num_ = 2 * 3 ** k_ * D_ + c_
                if num_ % Dl == 0:
                    s_ = num_ // Dl
                    if s_ > 0 and s_ + 2 * D_ > 0:
                        byword.setdefault(D_, set()).add((s_ + 2 * D_, w))
    for D_ in range(DLO, DHI + 1):
        if direct.get(D_, set()) != byword.get(D_, set()):
            okm = False
            print("   k-step mismatch", k_, D_)
            break
        if D_ >= 1:
            # for D >= 1 the sign conditions are automatic: only words with 2^p > 3^k occur
            cntD = sum(1 for p_ in range(k_, PM + 1) if 2 ** p_ > 3 ** k_
                       for w in compositions(p_, k_) if (2 * 3 ** k_ * D_ + cnum(w)) % (2 ** p_ - 3 ** k_) == 0) \
                if D_ <= 40 else len(direct.get(D_, ()))
            if cntD != len(direct.get(D_, ())):
                okm = False
    DIRECT[k_] = direct
    info("k=%d: A <= %d enumerated; multiplicities M_k(D) for D = 1..12: %s; M_k(0) = %d" % (
        k_, Abig, [len(direct.get(D_, ())) for D_ in range(1, 13)], len(direct.get(0, ()))))
check(okm, "k-step multiplicity theorem (k = 2,3,4; -30 <= D <= 300): odd A > 0 with Syr^k(A) = A - 2D "
           "<-> words w with (2^p-3^k) | 2*3^k D + c(w), s = quotient > 0, s + 2D > 0 (A = s + 2D); "
           "for D >= 1 exactly the words with 2^p > 3^k and the divisibility (D <= 40 recounted)")
# D = 0: periodic points; the size bound A <= ((3/2)^k - 1) 2^p0/(2^p0 - 3^k)
okc = True
for k_ in range(1, 13):
    bnd = Amax_bound(k_, 0)
    per = [A_ for A_ in range(1, int(bnd) + 2, 2) if orbit_word(A_, k_)[0] == A_]
    okc = okc and per == [1]
    info("k=%2d: every positive A with Syr^k(A) = A satisfies A <= %.2f; found %s" % (k_, float(bnd), per))
check(okc, "D = 0 (the cycle equation): for k <= 12 the only positive periodic point is A = 1 (classical, tiny range)")
# negative side: words with 2^p < 3^k and Delta | c give the negative cycles
negc = set()
for k_ in range(1, 8):
    for p_ in range(k_, 2 * k_ + 1):
        Dl = 2 ** p_ - 3 ** k_
        for w in compositions(p_, k_):
            c_ = cnum(w)
            if c_ % Dl == 0:
                negc.add(c_ // Dl if Dl < 0 else None) if Dl < 0 else None
check(sorted(x for x in negc if x is not None) == [-91, -61, -55, -41, -37, -25, -17, -7, -5, -1],
      "control: D = 0 words with 2^p < 3^k and (2^p-3^k) | c(w) (k <= 7) give exactly the negative cycles "
      "{-1}, {-5,-7}, {-17,...,-91} of 3A+1")
# mean multiplicity mu_k = sum_{w: 2^p > 3^k} 1/(2^p - 3^k) = sum_p C(p-1,k-1)/(2^p-3^k): rigorous enclosure
MB = 200   # fixed-point bits


def mu_enclosure(k_):
    P0 = p0(k_)
    PE = P0 + 40 * k_ + 400
    lo_ = 0
    for p_ in range(P0, PE + 1):
        lo_ += (comb(p_ - 1, k_ - 1) << MB) // (2 ** p_ - 3 ** k_)
    hi_ = lo_ + (PE - P0 + 1)
    # tail p > PE: 1/(2^p - 3^k) <= 2^(1-p) since 2^(p-1) >= 3^k; t_p = C(p-1,k-1) 2^(1-p) has ratio p/(2(p-k+1)) <= r
    r_ = Fraction(PE + 1, 2 * (PE + 2 - k_))
    assert r_ < 1 and 2 ** PE >= 3 ** k_
    hi_ += math.ceil(Fraction(comb(PE, k_ - 1), 2 ** PE) / (1 - r_) * 2 ** MB)
    return Fraction(lo_, 2 ** MB), Fraction(hi_, 2 ** MB)


def nu_exact(k_):
    """mean multiplicity of negative k-step drops (beyond the exceptional range): words with 2^p < 3^k"""
    return sum(Fraction(comb(p_ - 1, k_ - 1), 3 ** k_ - 2 ** p_) for p_ in range(k_, p0(k_)))


MU = {}
NU = {}
KLIST = list(range(1, 21)) + [29, 41, 53, 65, 100, 306]
for k_ in KLIST:
    MU[k_] = mu_enclosure(k_)
    NU[k_] = nu_exact(k_)
    P0 = p0(k_)
    spike = Fraction(comb(P0 - 1, k_ - 1), 2 ** P0 - 3 ** k_)
    info("k=%3d: p0 = %3d, 2^p0/3^k = %.6f, mu_k = %s (enclosure width %.0e), near-gate term "
         "g_k = C(p0-1,k-1)/(2^p0-3^k) = %.3e, P(p_k = p0) = %.3e, nu_k = %.6f" % (
             k_, P0, 2 ** P0 / 3 ** k_, mpmath.nstr(mpmath.mpf(MU[k_][0].numerator) / MU[k_][0].denominator, 12),
             float(MU[k_][1] - MU[k_][0]), float(spike), comb(P0 - 1, k_ - 1) / 2 ** P0, float(NU[k_])))
check(MU[1][0] <= QONE_HI and QONE_LO <= MU[1][1] and all(MU[k_][1] - MU[k_][0] < Fraction(1, 10 ** 50) for k_ in KLIST),
      "mu_1 = sum_(v>=2) 1/(2^v-3) (two independent enclosures overlap); all mu_k enclosures have width < 1e-50")
MUSTR = {2: "0.807697", 3: "1.859436", 4: "1.020035", 5: "3.512549", 6: "1.157888", 7: "0.922874", 12: "0.944261",
         17: "2.002368", 29: "1.618741", 41: "1.589768", 53: "0.997363", 100: "1.000274", 306: "1.000001"}
NUSTR = {2: "2.200000", 3: "0.325359", 7: "1.603846", 12: "4.509745", 53: "1.474409", 100: "0.000441"}
check(all("%.6f" % float(MU[k_][0]) == sv == "%.6f" % float(MU[k_][1]) for k_, sv in MUSTR.items()) and
      all("%.6f" % float(NU[k_]) == sv for k_, sv in NUSTR.items()) and NU[2] == Fraction(11, 5) and NU[1] == 1,
      "mu_k (k = 2..7, 12, 17, 29, 41, 53, 100, 306) and nu_k (k = 2, 3, 7, 12, 53, 100) as tabulated in the note; "
      "spikes of mu_k at upper approximations p0/k of log2 3 (5, 17, 29, 41), of nu_k at lower ones (2, 7, 12, 53)")
check(Fraction(1, 10 ** 6) < MU[306][0] - 1 and MU[306][1] - 1 < Fraction(2, 10 ** 6) and
      abs(MU[100][0] - 1) < Fraction(3, 10 ** 4),
      "mu_100 = 1.000274, mu_306 = 1 + 1.14e-6: consistent with mu_k -> 1 (the limit itself is CONDITIONAL on a CITED Baker-type bound)")
okmu = True
ZEROFRAC = {}
for k_ in (2, 3, 4):
    XD = 3000
    Abig = int(Amax_bound(k_, XD)) + 2
    tot_ = 0
    hitD = set()
    for A_ in range(1, Abig + 1, 2):
        sk, _ = orbit_word(A_, k_)
        D_ = (A_ - sk) // 2
        if 1 <= D_ <= XD:
            tot_ += 1
            hitD.add(D_)
    ratio = tot_ / XD
    ZEROFRAC[k_] = 1 - len(hitD) / XD
    info("EMPIRICAL k=%d: (1/X) sum_(D<=X) M_k(D) = %.4f at X = %d; mu_%d = %.4f; share of 1 <= D <= X that are "
         "not k-step drops = %.4f" % (k_, ratio, XD, k_, float(MU[k_][0]), ZEROFRAC[k_]))
    okmu = okmu and abs(ratio / float(MU[k_][0]) - 1) < 0.05
check(okmu, "EMPIRICAL: average k-step multiplicity over 1 <= D <= 3000 within 5% of mu_k (k = 2,3,4)")
# unit gates: |2^p - 3^k| = 1 only for (p,k) = (1,1), (2,1), (3,2)
units = [(p_, k_) for k_ in range(1, 60) for p_ in range(1, 100) if abs(2 ** p_ - 3 ** k_) == 1]
check(units == [(1, 1), (2, 1), (3, 2)], "unit gates |2^p - 3^k| = 1 (p < 100, k < 60): only (1,1), (2,1), (3,2)")
okneg2 = True
for D_ in range(-200, 201):
    for w, c_ in (((1, 2), 5), ((2, 1), 7)):
        A_ = -16 * D_ - c_
        sk, ww = orbit_word(A_, 2)
        okneg2 = okneg2 and ww == w and (A_ - sk) // 2 == D_ and cnum(w) == c_
okunitcyc = [Fraction(cnum(w), 2 ** sum(w) - 3 ** len(w)) for w in ((1,), (2,), (1, 2), (2, 1))] == [-1, 1, -5, -7]
check(okneg2 and okunitcyc, "k = 2: the unit gate 2^3 - 3^2 = -1 makes every D a 2-step drop of A = -16D-5 (word (1,2)) and "
              "A = -16D-7 (word (2,1)) over the odd integers, positive exactly when D <= -1  (|D| <= 200); the unit words "
              "(1), (2), (1,2), (2,1) have x_w = -1, 1, -5, -7: the unit-gate cycles")
Abig = int(Amax_bound(2, -60)) + 2
cnt2 = {}
cnt3 = {}
for A_ in range(1, max(Abig, int(Amax_bound(3, -60)) + 2) + 1, 2):
    for k_, cdict in ((2, cnt2), (3, cnt3)):
        sk, _ = orbit_word(A_, k_)
        D_ = (A_ - sk) // 2
        if -60 <= D_ <= -1 and A_ <= int(Amax_bound(k_, -60)) + 2:
            cdict[D_] = cdict.get(D_, 0) + 1
info("M_2(D) for D = -1..-12: %s" % [cnt2.get(D_, 0) for D_ in range(-1, -13, -1)])
info("M_3(D) for D = -1..-12: %s" % [cnt3.get(D_, 0) for D_ in range(-1, -13, -1)])
check(min(cnt2.get(D_, 0) for D_ in range(-60, 0)) >= 2 and 0 in [cnt3.get(D_, 0) for D_ in range(-60, 0)] and
      len(DIRECT[2].get(3, ())) == 0 and len(DIRECT[3].get(5, ())) == 0 and len(DIRECT[4].get(7, ())) == 0 and
      ZEROFRAC[2] >= 1 - float(MU[2][1]),
      "k = 2 net ascents always have >= 2 copies; for k = 3 (no unit gate) some D < 0 have none; some D >= 1 are "
      "not k-step drops: M_2(3) = M_3(5) = M_4(7) = 0; share of non-2-step-drops in [1,3000] >= 1 - mu_2")
# negative k-step drops beyond the exceptional range |D| <= ((3/2)^k - 1)/2: only words with 2^p < 3^k, periodic in D
okneg = True
for k_ in (2, 3, 4):
    DN = 3000
    Abig = 2 * 3 ** k_ * DN + 2
    cntN = {}
    for A_ in range(1, Abig + 1, 2):
        sk, _ = orbit_word(A_, k_)
        D_ = (A_ - sk) // 2
        if -DN <= D_ <= -1:
            cntN[D_] = cntN.get(D_, 0) + 1
    words_lo = [w for p_ in range(k_, p0(k_)) for w in compositions(p_, k_)]
    thr = Fraction(3 ** k_ - 2 ** k_, 2 ** (k_ + 1))
    for D_ in range(-DN, 0):
        if -D_ > thr:
            pred = sum(1 for w in words_lo if (2 * 3 ** k_ * D_ + cnum(w)) % (3 ** k_ - 2 ** sum(w)) == 0)
            okneg = okneg and cntN.get(D_, 0) == pred
    if k_ == 2:
        okneg = okneg and all(cntN.get(D_, 0) == 2 + (D_ % 5 == 0) for D_ in range(-DN, 0))
    if k_ == 3:
        okneg = okneg and sum(cntN.get(D_, 0) for D_ in range(-210, -1)) == 68 == nu_exact(3) * 209
    info("k=%d: negative k-step drops, -%d <= D <= -1: mean multiplicity %.4f (nu_k = %.4f), share with none %.4f" % (
        k_, DN, sum(cntN.values()) / DN, float(nu_exact(k_)), sum(1 for D_ in range(-DN, 0) if D_ not in cntN) / DN))
check(okneg, "negative k-step drops (k = 2,3,4; -3000 <= D <= -1): for |D| > ((3/2)^k-1)/2 the multiplicity is "
             "#{w : 2^p < 3^k, (3^k-2^p) | 2*3^k D + c(w)} (periodic in D); M_2(D) = 2 + [5 | D] for all D <= -1; "
             "k = 3: period 209 = 19*11 carries exactly 68 = 209*nu_3 drops")
okcw = all(((cnum(w) == 2 ** (p_ - k_) * (3 ** k_ - 2 ** k_)) == (p_ == k_)) and cnum(w) <= 2 ** (p_ - k_) * (3 ** k_ - 2 ** k_)
           for k_ in range(1, 7) for p_ in range(k_, 2 * k_ + 5) for w in compositions(p_, k_))
okdist = True
for A_ in range(1, 4001, 2):
    for k_ in range(1, 6):
        sk, w = orbit_word(A_, k_)
        p_ = sum(w)
        xw = Fraction(cnum(w), 2 ** p_ - 3 ** k_)
        okdist = okdist and Fraction(A_ - sk) == Fraction(2 ** p_ - 3 ** k_, 2 ** p_) * (A_ - xw)
check(okcw and okdist, "c(w) <= 2^(p-k)(3^k-2^k) with equality iff w = (1,...,1) (all compositions, k <= 6, p <= 2k+4); "
                       "A - Syr^k(A) = (2^p-3^k)(A - x_w)/2^p with x_w = c(w)/(2^p-3^k) the rational cycle of the word "
                       "(odd A < 4000, k <= 5)")

# ---------------------------------------------------------------------------
section("H. Controls: 3A-1, 5A+1, and what is specific to the multiplier 3")


def syr_ab(A_, a_, b_):
    t_ = a_ * A_ + b_
    vv_ = v2(t_)
    return t_ >> vv_, vv_


# H1: 3A - 1
okh = True
for A_ in range(1, 400001, 2):
    s_, vv_ = syr_ab(A_, 3, -1)
    dm = (A_ - s_) // 2
    okh = okh and 6 * dm - 1 == (2 ** vv_ - 3) * s_ and dm == -drop(-A_)
check(okh, "3A-1: 6d - 1 = (2^v - 3) Syr^-(A), and d^-(A) = -d(-A) (conjugacy A -> -A), odd A < 4*10^5")
DM = 100000
cntm = {}
for A_ in range(1, 8 * DM + 2, 2):
    s_, vv_ = syr_ab(A_, 3, -1)
    dm = (A_ - s_) // 2
    if -DM <= dm <= DM:
        cntm[dm] = cntm.get(dm, 0) + 1


def m_minus(dv):
    if dv <= 0:
        return 1
    n_ = 6 * dv - 1
    return sum(1 for w in range(2, n_.bit_length() + 2) if n_ % q(w) == 0)


check(all(cntm.get(dv, 0) == m_minus(dv) for dv in range(-DM, DM + 1)) and
      all(syr_ab(1 - 4 * dv, 3, -1)[1] == 1 for dv in range(-1000, 1)),
      "3A-1 multiplicity: d <= 0 once (A = 1 - 4d, A = 1 mod 4), d >= 1 exactly #{v >= 2 : (2^v-3) | 6d-1} times (|d| <= 10^5)")
XM = 10 ** 7
mneg = np.ones(XM + 1, dtype=np.uint8)
w = 3
while q(w) <= 6 * XM:
    d0 = pow(6, -1, q(w)) % q(w)          # 6d - 1 = 0 mod q
    if d0 == 0:
        d0 = q(w)
    mneg[d0::q(w)] += 1
    w += 1
check(all(int(mneg[dv]) == m_minus(dv) for dv in range(1, 50001)), "3A-1 sieve agrees with the formula for d <= 50000")
hneg = np.bincount(mneg[1:])
info("3A-1: histogram of m^-(d), 1 <= d <= 10^7: %s" % {k_: int(c) for k_, c in enumerate(hneg) if c})
check(all(abs(int(hneg[k_]) / XM - float(DELTA[k_])) < 5e-5 for k_ in range(1, len(hneg))),
      "EMPIRICAL: 3A-1 multiplicity frequencies on [1,10^7] match the SAME densities delta_k (PROVED equal: same moduli)")
del mneg
per2 = [A_ for A_ in range(1, 200, 2) if syr_ab(syr_ab(A_, 3, -1)[0], 3, -1)[0] == A_]
# sign symmetry: over ALL odd integers every d is a drop exactly 2 + #{v >= 3 : q_v | 6d+1} times
DZ = 10000
cntZ = {}
for A_ in range(-8 * DZ - 1, 8 * DZ + 2, 2):
    dz = drop(A_)
    if -DZ <= dz <= DZ:
        cntZ[dz] = cntZ.get(dz, 0) + 1
check(all(cntZ.get(dz, 0) == 2 + sum(1 for w in range(3, 70) if (6 * dz + 1) % q(w) == 0) for dz in range(-DZ, DZ + 1))
      and all(drop(8 * dz + 1) == dz and drop(-4 * dz - 1) == dz for dz in range(-DZ, DZ + 1)),
      "over all odd integers A (both signs) every d in Z is a drop exactly 2 + #{v >= 3 : (2^v-3) | 6d+1} times; "
      "the two unit copies are A = 8d+1 (v = 2) and A = -4d-1 (v = 1)  (|d| <= 10^4)")
del cntZ
check(per2 == [1, 5, 7], "3A-1 has 2-step periodic points 1, 5, 7 (its 2-cycle {5,7}): identical drop statistics, "
                         "different cycles")
# H2: 5A + 1
DF = 30000
cnt5 = {}
for A_ in range(1, 8 * DF + 2, 2):
    s_, vv_ = syr_ab(A_, 5, 1)
    d5 = (A_ - s_) // 2
    assert 10 * d5 + 1 == (2 ** vv_ - 5) * s_
    if -DF <= d5 <= DF:
        cnt5[d5] = cnt5.get(d5, 0) + 1


def m5(dv):
    n_ = 10 * dv + 1
    return sum(1 for w in range(1, abs(n_).bit_length() + 4)
               if n_ % (2 ** w - 5) == 0 and n_ // (2 ** w - 5) > 0)


check(all(cnt5.get(dv, 0) == m5(dv) for dv in range(-DF, DF + 1)),
      "5A+1: 10d + 1 = (2^v - 5) Syr_5(A); d occurs #{v >= 1 : (2^v-5) | 10d+1, quotient > 0} times (|d| <= 30000)")
check(cnt5.get(0, 0) == 0 and all(cnt5.get(dv, 0) == 1 + (dv % 3 == 2) for dv in range(-DF, 0)),
      "5A+1: d = 0 never occurs; each d < 0 occurs 1 + [d = 2 mod 3] times (units: only 2^2 - 5 = -1)")
num5, D5, ncfg5, sh5, L5 = exact_pgf(80, a=5, vmin=3)
dist5 = [Fraction(x, D5) for x in num5]
eps5 = Fraction(1, 2 ** 80) / (1 - Fraction(5, 2 ** 81))
check(sum(dist5) == 1 and D5 == L5 and sum(i * x for i, x in enumerate(dist5)) == sum(Fraction(1, 2 ** w - 5) for w in range(3, 81)),
      "5A+1 exact PGF (T = 80): shared primes %s; mean = sum_(v=3..80) 1/(2^v-5)" % sorted(sh5))
info("5A+1, d >= 1: density of m_5(d) = 0, 1, 2, 3, 4: %s (each +- %.1e)" % (
    ", ".join("%.12f" % float(dist5[k_]) for k_ in range(5)), float(eps5)))
XF = 10 ** 6
m5s = np.zeros(XF + 1, dtype=np.uint8)
w = 3
while 2 ** w - 5 <= 10 * XF + 1:
    Qw = 2 ** w - 5
    d0 = (-pow(10, -1, Qw)) % Qw
    if d0 == 0:
        d0 = Qw
    m5s[d0::Qw] += 1
    w += 1
check(all(int(m5s[dv]) == m5(dv) for dv in range(1, 30001)), "5A+1 sieve agrees with the formula on 1 <= d <= 30000")
h5 = np.bincount(m5s[1:])
check(all(abs(int(h5[k_]) / XF - float(dist5[k_])) < 2e-3 for k_ in range(len(h5))),
      "EMPIRICAL: 5A+1 frequencies on [1,10^6] match the exact densities (P(m_5 = 0) = %.6f: most positive d never occur)"
      % float(dist5[0]))
Q5 = sum(Fraction(1, 2 ** w - 5) for w in range(3, 200))
info("5A+1 mean multiplicities: positive side sum_(v>=3) 1/(2^v-5) = %.12f ; negative side 1/3 + 1 = 4/3" % float(Q5))
# H3: only a = 3 has two units
two_units = [a_ for a_ in range(3, 100001, 2) if sum(1 for w in range(1, 20) if abs(2 ** w - a_) == 1) == 2]
check(two_units == [3], "among odd a in [3, 10^5], only a = 3 has two v with |2^v - a| = 1 (proof: 2^v' - 2^v = 2 forces v=1, v'=2)")
# H3b: for odd a >= 5 (b = 1) the side without a unit misses a positive density of drops
okmiss = True
for a_ in range(5, 4097, 2):
    j_ = a_.bit_length() - 1                       # 2^j < a < 2^(j+1)
    neg_unit = (a_ == 2 ** j_ + 1)
    pos_unit = (a_ == 2 ** (j_ + 1) - 1)
    okmiss = okmiss and not (neg_unit and pos_unit)
    if not neg_unit:
        sn = sum(Fraction(1, a_ - 2 ** w) for w in range(1, j_ + 1))
        okmiss = okmiss and sn <= Fraction(11, 15)
    if not pos_unit:
        sp = sum(1.0 / (2 ** w - a_) for w in range(j_ + 1, j_ + 200))
        okmiss = okmiss and sp < 0.54
missed = {}
for a_ in (5, 7, 9, 11, 13, 15, 17, 31, 33):
    hit = set()
    for A_ in range(1, 8 * 2000 * a_, 2):
        s_, vv_ = syr_ab(A_, a_, 1)
        dv = (A_ - s_) // 2
        if -2000 <= dv <= 2000:
            hit.add(dv)
    missed[a_] = (sum(1 for dv in range(-2000, 0) if dv not in hit), sum(1 for dv in range(1, 2001) if dv not in hit))
info("odd a: (# of -2000 <= d <= -1, # of 1 <= d <= 2000) that are NOT drops of (aA+1)/2^v: %s" % missed)
check(okmiss and all(max(x) > 0 for x in missed.values()),
      "for odd 5 <= a < 4097 (b = 1): a = 3 is the only case with both units; the side without a unit has covering "
      "density <= 11/15 (negative side) or < 0.54 (positive side), and drops are missed in practice (a <= 33 sample)")
# H4: general aA + b identity
okab = True
for (a_, b_) in ((3, 1), (3, -1), (5, 1), (5, 3), (7, 1), (3, 5), (9, 1)):
    for A_ in range(1, 20001, 2):
        if a_ * A_ + b_ == 0:
            continue
        s_, vv_ = syr_ab(A_, a_, b_)
        if (A_ - s_) % 2:
            okab = False
            break
        dv = (A_ - s_) // 2
        okab = okab and 2 * a_ * dv + b_ == (2 ** vv_ - a_) * s_
check(okab, "general map (aA+b)/2^v: 2a d + b = (2^v - a) Syr(A) for (a,b) in {(3,+-1),(5,1),(5,3),(7,1),(3,5),(9,1)}")

# ---------------------------------------------------------------------------
section("I. 2-adic structure; consecutive drops; relation to the gates")
KB = 12
okb2 = True
for vv_ in range(1, 9):
    mod = 2 ** (KB + vv_ + 1)
    seen = set()
    for A_ in range(1, mod, 2):
        if v2(3 * A_ + 1) == vv_:
            dv = (((2 ** vv_ - 3) * A_ - 1) >> (vv_ + 1)) % (2 ** KB)
            seen.add(dv)
    okb2 = okb2 and len(seen) == 2 ** KB
check(okb2, "2-adically, d restricted to the branch {v_2(3A+1) = v} is an affine bijection onto Z_2: "
            "A mod 2^(K+v+1) -> d mod 2^K is onto Z/2^K (K = 12, v = 1..8)")
MJ = 16
cnts = {}
for A_ in range(1, 2 ** MJ, 2):
    s1, a1 = syr(A_)
    s2, b1 = syr(s1)
    cnts[(a1, b1)] = cnts.get((a1, b1), 0) + 1
check(all(cnts.get((a_, b_), 0) == 2 ** (MJ - 1 - a_ - b_) for a_ in range(1, 8) for b_ in range(1, 8) if a_ + b_ <= MJ - 2),
      "consecutive valuations (v_1, v_2) of odd A mod 2^16 are exactly i.i.d. geometric: #{(a,b)} = 2^(15-a-b)")
signs = {(x, y): sum(c for (a_, b_), c in cnts.items() if (a_ == 1) == x and (b_ == 1) == y) for x in (0, 1) for y in (0, 1)}
check(all(abs(signs[k_] / 2 ** (MJ - 1) - 0.25) < 2 ** -10 for k_ in signs),
      "signs of consecutive drops (ascent iff v = 1) are independent fair coins (exact 2-adically from P(a,b) = 2^(-a-b); "
      "frequencies over odd A < 2^16 within 2^-10 of 1/4)")
okc2 = True
for A_ in range(1, 200001, 2):
    s1, a1 = syr(A_)
    s2, b1 = syr(s1)
    n1 = 6 * ((A_ - s1) // 2) + 1
    n2 = 6 * ((s1 - s2) // 2) + 1
    okc2 = okc2 and n2 * 2 ** b1 == q(b1) * (3 * (n1 // q(a1)) + 1) and n1 % q(a1) == 0
check(okc2, "consecutive drops: with n_i = 6d_i + 1, n_2 = q_(v_2) (3 n_1/q_(v_1) + 1)/2^(v_2) (odd A < 2*10^5)")
# gates: the k = 1 moduli are the gates 2^p - 3^1; the k-step class of drops of a word is D = -c(w)/(2*3^k) mod |2^p - 3^k|
okgate = True
for k_ in range(1, 6):
    for p_ in range(k_, 2 * k_ + 3):
        Dl = 2 ** p_ - 3 ** k_
        for w in compositions(p_, k_):
            c_ = cnum(w)
            Dw = (-c_ * pow(2 * 3 ** k_, -1, abs(Dl))) % abs(Dl)
            okgate = okgate and ((Dw == 0) == (c_ % Dl == 0))
check(okgate, "each word w carries one class of k-step drops D = -c(w)/(2*3^k) mod |2^p-3^k|; the class contains 0 "
              "iff (2^p - 3^k) | c(w) (the repo's cycle gate q = 1)  (k <= 5)")
okord = True
for k_ in range(1, 7):
    for p_ in range(k_, 2 * k_ + 3):
        Dl = abs(2 ** p_ - 3 ** k_)
        for w in compositions(p_, k_):
            c_ = cnum(w)
            Dw = (-c_ * pow(2 * 3 ** k_, -1, Dl)) % Dl
            okord = okord and Dl // gcd(Dw, Dl) == Dl // gcd(c_, Dl)
check(okord, "the repo's minimal gate parameter q(w) = |Delta|/gcd(c(w), |Delta|) is the additive order of the "
             "word's drop class D_w in Z/|Delta|  (k <= 6, p <= 2k+2)")
ok3 = True
for j in (1, 2, 3):
    tab = {}
    for A_ in range(1, 200001, 2):
        s_, vv_ = syr(A_)
        key = (A_ % 3 ** j, vv_ % (2 * 3 ** (j - 1)))
        val = ((A_ - s_) // 2) % 3 ** j
        if tab.setdefault(key, val) != val:
            ok3 = False
    if j <= 2:
        ok3 = ok3 and len(tab) == 3 ** j * 2 * 3 ** (j - 1)
ok3 = ok3 and all(syr(A_)[0] % 3 == (1 if syr(A_)[1] % 2 == 0 else 2) for A_ in range(1, 200001, 2))
check(ok3, "3-adic: d mod 3^j is a function of (A mod 3^j, v mod 2*3^(j-1)) (j = 1,2,3; for j <= 2 all pairs occur); "
           "Syr(A) = (-1)^v mod 3")

print("number of checks: %d" % NCHECK)
sys.stderr.write("total time %.1fs\n" % (time.time() - T_START))
print("ALL CHECKS PASSED")
