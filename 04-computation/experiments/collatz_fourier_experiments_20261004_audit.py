#!/usr/bin/env python
"""
collatz_fourier_experiments_20261004_audit.py  --  independent audit (blind re-derivation)

Audits the numbers of 05-knowledge/results/collatz_fourier_experiments_20261004.md.
Written WITHOUT reading the session's scripts (collatz_fixed_frequency_*, collatz_phase_tower_*,
collatz_syracuse_law_*, collatz_fourier_mass_*, collatz_digit_phase_*).  Python 3 + numpy only.

Objects, re-derived from the definitions:
  q-adic Syracuse law   Y_0 = 0,  Y_n = 2^-a (q Y_{n-1} + 1) mod q^n,  P(a = k) = 2^-k (k >= 1).
  Fourier coefficient   mu_hat_n(t) = E e(t Y_n / q^n),  e(x) = exp(2 pi i x).
  Recursion             mu_hat_n(t) = sum_a 2^-a e((t 2^-a mod q^n)/q^n) mu_hat_{n-1}(t 2^-a mod q^{n-1})
                        (derived: t Y_n/q^n = s(qY_{n-1}+1)/q^n with s = t 2^-a mod q^n, and
                         s Y_{n-1}/q^{n-1} mod 1 depends on s mod q^{n-1} only).
  Window family         f_n(k) = mu_hat_n(u 2^k), k <= 0:  f_n(k) = sum_a 2^-a omega_n(k-a) f_{n-1}(k-a),
                        omega_n(j) = e((u 2^j mod q^n)/q^n), f_0 = 1;  window k in [-(N-n)A, 0], a <= A.
  Phase device          theta_{n,d} := (u 2^-d mod q^n)/q^n = (R mod 2^d)/2^d + u 2^-d q^-n  (mod 1),
                        R = -u q^-n mod 2^(d+1).  Derivation: with x = u 2^-d mod q^n, 2^d x = u + m q^n,
                        so x/q^n = m/2^d + u/(2^d q^n) and m = -u q^-n mod 2^d.  Verified exactly below
                        BEFORE it is used; for large n the bits of R are read by a 62-tap convolution.
Usage:  python collatz_fourier_experiments_20261004_audit.py [--quick]
"""
import sys, time, math, itertools
from fractions import Fraction
from math import comb
import numpy as np

QUICK = "--quick" in sys.argv
OUTPATH = __file__[:-3] + ".out"
_outf = open(OUTPATH, "w", encoding="ascii", errors="replace")
T0 = time.time()


def log(*args):
    s = " ".join(str(a) for a in args)
    print(s)
    _outf.write(s + "\n")
    _outf.flush()


def e(x):
    return np.exp(2j * np.pi * np.asarray(x, dtype=float))


LN3 = math.log(3.0)

# =====================================================================================
# A. The phase device, verified exactly (Fractions) before use
# =====================================================================================
log("=" * 100)
log("A. Exact check of the phase device theta_{n,d} = (R mod 2^d)/2^d + u 2^-d q^-n, R = -u q^-n mod 2^(d+1)")
log("=" * 100)


def theta_exact(u, q, n, d):
    qn = q ** n
    return Fraction((u * pow(2, -d, qn)) % qn, qn)


def theta_device(u, q, n, d):
    qn = q ** n
    R = (-u * pow(qn, -1, 1 << (d + 1))) % (1 << (d + 1))
    return Fraction(R % (1 << d), 1 << d) + Fraction(u, (1 << d) * qn)


bad = 0
cnt = 0
for q in (3, 5, 7, 11):
    for u in (1, 2, 5, 7, 11, 13):
        for n in range(1, 9):
            for d in range(1, 31):
                cnt += 1
                if (theta_exact(u, q, n, d) - theta_device(u, q, n, d)) % 1 != 0:
                    bad += 1
log(f"  checked {cnt} cases (q in 3,5,7,11; u in 1,2,5,7,11,13; n<=8; d<=30): violations = {bad}")
# the bits: b_i (i-th binary digit of -u q^-n) is bit i-1 of R; (R mod 2^d)/2^d = 0.b_d b_{d-1} ... b_1
u, q, n, d = 1, 3, 4, 12
R = (-u * pow(q ** n, -1, 1 << 40)) % (1 << 40)
bits = [(R >> i) & 1 for i in range(40)]
val = sum(bits[i] * Fraction(1, 2 ** (d - i)) for i in range(d))
log(f"  digit reading at (u,q,n,d)=({u},{q},{n},{d}): 0.b_d...b_1 = {val} ; (R mod 2^d)/2^d = {Fraction(R % (1 << d), 1 << d)} ; equal = {val == Fraction(R % (1 << d), 1 << d)}")
log(f"  elapsed {time.time() - T0:.1f}s")

# =====================================================================================
# B. Direct enumeration of the law (dynamic programming on Z/q^m), and its Fourier coefficients
# =====================================================================================
log("=" * 100)
log("B. Direct enumeration of the law on Z/q^n (DP, valuations <= A, or exact via the period of 2 mod q^n)")
log("=" * 100)


def mult_order_2(m):
    k, x = 1, 2 % m
    while x != 1:
        x = (2 * x) % m
        k += 1
    return k


def dense_laws(q, n, A=30, exact=False, keep_all=True):
    """Returns the list of laws [level 0 .. n]; law[m] is a float array on Z/q^m (y -> P(Y_m = y))."""
    laws = [np.array([1.0])]
    for m in range(1, n + 1):
        qm = q ** m
        prev = laws[-1]
        x = np.arange(len(prev), dtype=np.int64)
        y = (q * x + 1) % qm
        inv2 = pow(2, -1, qm)
        law = np.zeros(qm)
        if exact:
            L = mult_order_2(qm)
            norm = 1.0 / (1.0 - 2.0 ** (-L))
            rng_a = range(1, L + 1)
        else:
            rng_a = range(1, A + 1)
            norm = 1.0
        for a in rng_a:
            y = (y * inv2) % qm
            law += (2.0 ** (-a) * norm) * np.bincount(y, weights=prev, minlength=qm)
        if keep_all:
            laws.append(law)
        else:
            laws = [law]
    return laws


def fourier_direct(law, t, q, m):
    qm = q ** m
    return np.sum(law * e(t * np.arange(qm) / qm))


CLAIMED_A = {(3, 1): 0.1561397846, (5, 1): 0.0860179286, (7, 1): 0.0212430711, (3, 5): 0.1376897005, (7, 7): 0.0129696631}
laws3_A30 = dense_laws(3, 7, A=30)
laws3_exact = dense_laws(3, 7, exact=True)
laws3_A40 = dense_laws(3, 7, A=40)
log("  q = 3: |mu_hat_n(u)| by direct enumeration; A = 30 (as in the note's control), A = 40, and exact (period of 2 mod 3^n)")
log("   n  u   A=30            A=40            exact           claimed")
direct_vals = {}
for n in (3, 5, 7):
    for u in (1, 5, 7):
        v30 = abs(fourier_direct(laws3_A30[n], u, 3, n))
        v40 = abs(fourier_direct(laws3_A40[n], u, 3, n))
        vx = abs(fourier_direct(laws3_exact[n], u, 3, n))
        direct_vals[(n, u)] = (v30, v40, vx)
        cl = CLAIMED_A.get((n, u))
        log(f"   {n}  {u}   {v30:.12f}  {v40:.12f}  {vx:.12f}  {cl if cl is not None else '-'}")
log(f"  elapsed {time.time() - T0:.1f}s")

# =====================================================================================
# C. The window recursion with my phase device (large n), cross-checked against the direct law
# =====================================================================================
log("=" * 100)
log("C. Window recursion f_n(k), k in [-(N-n)A, 0], phases from the binary digits of -u q^-n")
log("=" * 100)

KTAPS = 62
KERNEL = 2.0 ** (-np.arange(1, KTAPS + 1))


def phases_from_bits(bits_float):
    """theta_m = 0.b_m b_{m-1} ... b_1 for m = 1..M from bits[0..M-1] = b_1..b_M (top 62 bits)."""
    c = np.convolve(bits_float, KERNEL)
    return c[: len(bits_float)]


def level_phases(u, q, n, M, qn):
    """array P[m-1] = e(theta_{n,m}), m = 1..M, theta_{n,m} = (u 2^-m mod q^n)/q^n."""
    mod = 1 << (M + 1)
    R = (-u * pow(q, -n, mod)) % mod
    nbytes = (M + 1 + 7) // 8
    bits = np.unpackbits(np.frombuffer(R.to_bytes(nbytes, "little"), dtype=np.uint8), bitorder="little")[:M].astype(np.float64)
    theta = phases_from_bits(bits)
    corr0 = u / qn if qn.bit_length() < 1000 else 0.0   # u q^-n as float (0 when below 2^-1000)
    with np.errstate(under="ignore"):
        corr = np.ldexp(corr0, -np.arange(1, M + 1))
    return e(theta + corr)


# exact cross-check of the floating phases against modular inverses at a moderate level
u, q, n, M = 7, 3, 40, 3000
P = level_phases(u, q, n, M, q ** n)
qn = q ** n
err = max(abs(P[m - 1] - complex(e(float(Fraction((u * pow(2, -m, qn)) % qn, qn))))) for m in range(1, M + 1))
log(f"  float phases vs modular-inverse phases at (u,q,n)=({u},{q},{n}), m<=3000: max |diff| = {err:.2e}")
u, q, n, M = 5, 5, 30, 2000
P = level_phases(u, q, n, M, q ** n)
qn = q ** n
err = max(abs(P[m - 1] - complex(e(float(Fraction((u * pow(2, -m, qn)) % qn, qn))))) for m in range(1, M + 1))
log(f"  float phases vs modular-inverse phases at (u,q,n)=({u},{q},{n}), m<=2000: max |diff| = {err:.2e}")


def window_run(u, q, N, A, rescale=True, phase_source="real", rng=None, record=None):
    """Runs the window recursion to level N.  Returns logabs[n] = ln|f_n(0)| = ln|mu_hat_n(u)|, n = 0..N.
    phase_source: 'real' (digits of -u q^-n) or 'iid' (fresh uniform bits at every level, no correction term).
    Internally v_n = 3^{n/2} f_n (rescale) to keep numbers O(1); the returned logabs REMOVES the rescaling."""
    w = 2.0 ** (-np.arange(1, A + 1))
    g = np.ones(N * A + 1, dtype=complex)      # level 0: f_0(k) = 1 for every k
    logabs = np.zeros(N + 1)
    scale = math.sqrt(3.0) if rescale else 1.0
    qn = 1
    for n in range(1, N + 1):
        qn *= q
        M = (N - n + 1) * A
        if phase_source == "real":
            Pn = level_phases(u, q, n, M, qn)
        else:
            bits = rng.integers(0, 2, size=M).astype(np.float64)
            Pn = e(phases_from_bits(bits))
        H = Pn * g[1 : M + 1]
        Wn = (N - n) * A
        new = np.zeros(Wn + 1, dtype=complex)
        for a in range(1, A + 1):
            new += w[a - 1] * H[a - 1 : a + Wn]
        new *= scale
        g = new
        if record is not None and n in record:
            record[n] = g.copy()
        logabs[n] = math.log(abs(g[0])) - (0.5 * n * LN3 if rescale else 0.0)
    return logabs


# (a) agreement of the window recursion with the direct enumeration at n = 3, 5, 7
log("  (a) window recursion vs direct enumeration (q = 3):")
log("   n  u   window A=30     window A=40     |win30-DP30|  |win40-exact|  claimed         |win40-claimed|")
for u in (1, 5, 7):
    la30 = window_run(u, 3, 7, 30)
    la40 = window_run(u, 3, 7, 40)
    la40_nr = window_run(u, 3, 7, 40, rescale=False)
    for n in (3, 5, 7):
        v30, v40, vx = direct_vals[(n, u)]
        w30, w40 = math.exp(la30[n]), math.exp(la40[n])
        cl = CLAIMED_A.get((n, u))
        extra = f"{cl:.10f}      {abs(w40 - cl):.1e}" if cl is not None else "-"
        log(f"   {n}  {u}   {w30:.12f}  {w40:.12f}  {abs(w30 - v30):.1e}       {abs(w40 - vx):.1e}        {extra}")
    log(f"      rescaled vs unrescaled recursion (u={u}, n<=7): max |diff of ln| = {np.max(np.abs(la40 - la40_nr)):.1e}")
# q = 5, 7 small-level controls
for q in (5, 7):
    laws = dense_laws(q, 5, A=40)
    for u in (1, 2, 5):
        la = window_run(u, q, 5, 40)
        diffs = [abs(math.exp(la[n]) - abs(fourier_direct(laws[n], u, q, n))) for n in (3, 4, 5)]
        log(f"  q={q} u={u}: |mu_hat_n(u)| n=3,4,5 = {[f'{math.exp(la[n]):.10f}' for n in (3, 4, 5)]}; window-vs-DP max diff {max(diffs):.1e}")
log(f"  elapsed {time.time() - T0:.1f}s")


def ls_rate(logabs, n0, n1):
    ns = np.arange(n0, n1 + 1)
    slope = np.polyfit(ns, logabs[ns], 1)[0]
    return math.exp(slope)


def local_maxima(v, thresh=1.0):
    out = []
    for n in range(1, len(v) - 1):
        if v[n] >= thresh and v[n] > v[n - 1] and v[n] > v[n + 1]:
            out.append((n, round(float(v[n]), 2)))
    return out


# (b) q = 3, u = 1 to N = 800 (and the other five units of E1)
log("=" * 100)
log("C(b). Rates. q = 3, N = 800, A = 40; LS rate = exp(slope of ln|mu_hat_n(u)| vs n), endpoints inclusive")
log("=" * 100)
N800 = 800 if not QUICK else 300
E1_CLAIMS = {1: (0.5680, "(127, 3.88), (129, 7.22), (131, 10.38), (135, 2.62)"), 5: (0.5695, "(68, 2.69)"), 7: (0.5687, "(11,1.80),(49,1.27)"),
             11: (0.5709, "(24,2.57),(32,2.06),(180,1.57)"), 13: (0.5703, "(34,6.89),(36,6.98),(72,1.62)"), 17: (0.5704, "(27,2.08),(92,1.28)")}
logabs_store = {}
for u in (1, 5, 7, 11, 13, 17):
    t1 = time.time()
    la = window_run(u, 3, N800, 40)
    logabs_store[(3, u)] = la
    v = np.exp(la + 0.5 * np.arange(N800 + 1) * LN3)   # 3^{n/2}|mu_hat_n(u)|
    if N800 >= 800:
        meds = [np.median(v[1:201]), np.median(v[201:401]), np.median(v[401:601]), np.median(v[601:801])]
        r = (ls_rate(la, 200, 400), ls_rate(la, 400, 600), ls_rate(la, 600, 800), ls_rate(la, 200, 800))
        log(f"  u={u:2d}: block medians of 3^(n/2)|mu_hat| = {['%.1e' % m for m in meds]}; rates 200-400 {r[0]:.4f}, 400-600 {r[1]:.4f}, 600-800 {r[2]:.4f}, 200..800 {r[3]:.4f} (claimed {E1_CLAIMS[u][0]})")
        log(f"         also 200..600 {ls_rate(la, 200, 600):.4f}, 100..300 {ls_rate(la, 100, 300):.4f}, 300..600 {ls_rate(la, 300, 600):.4f}, 201..800 {ls_rate(la, 201, 800):.4f}")
        log(f"         3^(300)|mu_hat_600| = {v[600]:.2e};  local maxima >= 1 of 3^(n/2)|mu_hat_n|: {local_maxima(v)}  (claimed {E1_CLAIMS[u][1]})")
    else:
        log(f"  u={u:2d}: rate 100..{N800} {ls_rate(la, 100, N800):.4f}; local maxima >= 1: {local_maxima(v)}")
    log(f"         ({time.time() - t1:.1f}s)")
if N800 >= 800:
    la = logabs_store[(3, 1)]
    v = np.exp(la + 0.5 * np.arange(N800 + 1) * LN3)
    log(f"  u=1 wave: 3^(n/2)|mu_hat_n(1)| at n = 125..136: {[f'{x:.3f}' for x in v[125:137]]}")
    # truncation control
    la60 = window_run(1, 3, 400, 60)
    la40 = window_run(1, 3, 400, 40)
    rel = np.max(np.abs(np.exp(la60[1:401] - la40[1:401]) - 1.0))
    log(f"  truncation control u=1: max relative difference |mu_hat_n| A=60 vs A=40, n<=400: {rel:.1e} (claimed < 3e-9)")
    la_nr = window_run(1, 3, 400, 40, rescale=False)
    log(f"  rescaled vs unrescaled recursion to n=400: max |diff of ln| = {np.max(np.abs(la_nr - la40)):.1e}")
log(f"  elapsed {time.time() - T0:.1f}s")

# (b) q = 5, 7 (and 11) with u = 1 (and the second unit of E4) to N = 600
log("=" * 100)
log("C(b). q = 5, 7, 11: N = 600, A = 40; rate over 200..600 of |mu_hat_n(u)| itself; q^(n/2)|mu_hat_600(u)| via logs")
log("=" * 100)
N600 = 600 if not QUICK else 250
E4 = [(3, 1, 0.5721), (3, 5, 0.5698), (5, 1, 0.5735), (5, 2, 0.5738), (7, 1, 0.5716), (7, 5, 0.5718), (11, 1, 0.5720), (11, 5, 0.5734)]
for q, u, claim in E4:
    t1 = time.time()
    la = window_run(u, q, N600, 40)
    logabs_store[(q, u, 600)] = la
    if N600 >= 600:
        r1, r2, r3 = ls_rate(la, 100, 300), ls_rate(la, 300, 600), ls_rate(la, 200, 600)
        l10 = (la[600] + 300 * math.log(q)) / math.log(10)
        l10_3 = (la[600] + 300 * LN3) / math.log(10)
        log(f"  q={q:2d} u={u}: rates 100..300 {r1:.4f}, 300..600 {r2:.4f}, 200..600 {r3:.4f} (claimed {claim}); Parseval q^-1/2 = {q ** -0.5:.4f}; "
            f"ratio {r3 / q ** -0.5:.3f}; q^(300)|mu_hat_600| = 10^{l10:.2f}; 3^(300)|mu_hat_600| = 10^{l10_3:.2f}; |mu_hat_600| = 10^{la[600] / math.log(10):.2f}  ({time.time() - t1:.1f}s)")
    else:
        log(f"  q={q:2d} u={u}: rate 100..{N600} {ls_rate(la, 100, N600):.4f} (claimed {claim} on 200..600)")
log(f"  elapsed {time.time() - T0:.1f}s")

# =====================================================================================
# D. Exact identities of section 4j (Fractions)
# =====================================================================================
log("=" * 100)
log("D. Section 4j identities in exact rational arithmetic: theta_{n,d} = (u 2^-d mod q^n)/q^n in Q/Z")
log("=" * 100)


def th(u, q, n, d):
    qn = q ** n
    return Fraction((u * pow(2, -d, qn)) % qn, qn)


def is_int(x):
    return x.denominator == 1


for q in (3, 5):
    bad_i = bad_ii = bad_iii = 0
    n_i = n_ii = n_iii = 0
    for u in (1, 5, 7, 11):
        for n in range(2, 12):
            for d in range(1, 25):
                n_i += 1
                if not is_int(q * th(u, q, n + 1, d) - th(u, q, n, d)):
                    bad_i += 1
                if d >= 2:
                    n_ii += 1
                    if not is_int(th(u, q, n, d - 1) - 2 * th(u, q, n, d)):
                        bad_ii += 1
        N = 14
        for m in range(0, 7):
            for d in range(m + 1, 25):
                if q == 3:
                    s = sum(comb(m, i) * th(u, q, N, d - i) for i in range(m + 1))
                else:  # q = 5 = 1 + 2^2: step 2  (needs d - 2i >= 0; th(.,0) = (u mod q^N)/q^N is the natural extension)
                    s = sum(comb(m, i) * th(u, q, N, d - 2 * i) for i in range(m + 1) if d - 2 * i >= 0)
                    if any(d - 2 * i < 0 for i in range(m + 1)):
                        continue
                n_iii += 1
                if not is_int(th(u, q, N - m, d) - s):
                    bad_iii += 1
    log(f"  q={q}: (i) q theta_(n+1,d) = theta_(n,d): {bad_i} violations / {n_i};  (ii) theta_(n,d-1) = 2 theta_(n,d): {bad_ii} / {n_ii};  "
        f"(iii) Pascal (N=14, m<=6, m<d<=24): {bad_iii} / {n_iii}")
# q = 7 = 2^3 - 1: alternating signs  theta_{N-m,d} = sum_i C(m,i) (-1)^(m-i) theta_{N,d-3i}
bad = tot = 0
for u in (1, 5, 11):
    for m in range(0, 5):
        for d in range(3 * m + 1, 25):
            s = sum(comb(m, i) * (-1) ** (m - i) * th(u, 7, 14, d - 3 * i) for i in range(m + 1))
            tot += 1
            if not is_int(th(u, 7, 14 - m, d) - s):
                bad += 1
log(f"  q=7 (= 2^3 - 1, alternating signs, step 3): {bad} violations / {tot}")
log(f"  elapsed {time.time() - T0:.1f}s")

# =====================================================================================
# E. Section 4k: injectivity of Phi_q, its proof, and the fourth moments
# =====================================================================================
log("=" * 100)
log("E. Section 4k: Phi_q(a) = sum_i q^(n-i) 2^-D_i mod 1, D_i = a_i + ... + a_n (suffix sums)")
log("=" * 100)


def phi_scaled(path, q, S):
    """2^S * Phi_q(path) mod 2^S as an exact integer (S >= D_1)."""
    n = len(path)
    D = 0
    tot = 0
    for i in range(n - 1, -1, -1):   # i = n-1 .. 0  (1-based index i+1)
        D += path[i]
        tot += q ** (n - 1 - i) * (1 << (S - D))
    return tot % (1 << S)


# injectivity, q = 3, valuations <= 12, n <= 5
for n in range(1, 6):
    vals = set()
    cnt = 0
    for path in itertools.product(range(1, 13), repeat=n):
        vals.add(phi_scaled(path, 3, 60))
        cnt += 1
    log(f"  q=3, n={n}, valuations<=12: {cnt} paths, {len(vals)} distinct values of Phi_3 -> injective = {cnt == len(vals)}")


def v2_of_fraction(fr):
    """2-adic valuation of a nonzero rational."""
    num, den = fr.numerator, fr.denominator
    v = 0
    while num % 2 == 0:
        num //= 2
        v += 1
    while den % 2 == 0:
        den //= 2
        v -= 1
    return v


def phi_frac(path, q):
    n = len(path)
    D = 0
    tot = Fraction(0)
    for i in range(n - 1, -1, -1):
        D += path[i]
        tot += Fraction(q ** (n - 1 - i), 1 << D)
    return tot % 1


# the proof: minimal index i0 where the DEPTHS differ has the unique minimal 2-adic valuation
n = 3
paths = list(itertools.product(range(1, 6), repeat=n))
ok_depth = 0
note_first_ok = note_last_ok = 0
half = 0
tot = 0
example = None
for a in paths:
    Da = [sum(a[i:]) for i in range(n)]
    for b in paths:
        if a >= b:
            continue
        Db = [sum(b[i:]) for i in range(n)]
        tot += 1
        diff = (phi_frac(a, 3) - phi_frac(b, 3)) % 1
        v = v2_of_fraction(diff) if diff != 0 else None
        i0 = next(i for i in range(n) if Da[i] != Db[i])
        if v == -max(Da[i0], Db[i0]):
            ok_depth += 1
        if diff == Fraction(1, 2):
            half += 1
        jf = next(i for i in range(n) if a[i] != b[i])          # first index where the paths differ
        jl = max(i for i in range(n) if a[i] != b[i])           # last index where the paths differ
        if v == -max(Da[jf], Db[jf]):
            note_first_ok += 1
        if v == -max(Da[jl], Db[jl]):
            note_last_ok += 1
        elif example is None:
            example = (a, b, Da, Db, diff, v, jl)
log(f"  n=3, valuations<=5, {tot} unordered pairs: v_2(Phi(a)-Phi(b)) = -max(D_i0(a),D_i0(b)) with i0 = first index where the DEPTHS differ: {ok_depth}/{tot};"
    f" differences equal to 1/2 mod 1: {half}")
log(f"    the note's bookkeeping 'denominator 2^max(D_j(a),D_j(b)) at the first index j where the paths differ': holds for {note_first_ok}/{tot};"
    f" with j = LAST differing index (so that the later terms cancel): {note_last_ok}/{tot}")
a, b = (3, 1), (1, 2)
log(f"    explicit counterexample to the note's sketch (n=2, q=3): a={a}, b={b}: depths {[4, 1]} vs {[3, 2]}; last differing index j=2, max(D_2)=2,"
    f" but Phi(a)-Phi(b) = {phi_frac(a, 3) - phi_frac(b, 3)} (denominator 2^4 = 2^max(D_1), not 2^2); the j=2 term alone is 1/2-1/4 = 1/4")

# fourth moments: i.i.d. model (per-level multiset key) and the random-start tower (additive-energy key)
log("  fourth moments scaled by 9^n, valuations <= 8 (claims: iid 1.400,1.704,1.976,2.241; q=3 1.400,2.550,3.915,7.047):")
# sanity check of the per-level i.i.d. expectation: E e(th_p + th_q - th_r - th_s) over 2^6 bit strings = multiset indicator
Dmax = 6
allbits = np.array(list(itertools.product((0, 1), repeat=Dmax)), dtype=float)   # columns = b_1..b_6
thetas = np.zeros((len(allbits), Dmax + 1))
for d in range(1, Dmax + 1):
    thetas[:, d] = sum(allbits[:, i - 1] * 2.0 ** (i - 1 - d) for i in range(1, d + 1))
worst = 0.0
for p in range(1, Dmax + 1):
    for q_ in range(1, Dmax + 1):
        for r in range(1, Dmax + 1):
            for s in range(1, Dmax + 1):
                ex = np.mean(e(thetas[:, p] + thetas[:, q_] - thetas[:, r] - thetas[:, s]))
                ind = 1.0 if sorted((p, q_)) == sorted((r, s)) else 0.0
                worst = max(worst, abs(ex - ind))
log(f"    per-level i.i.d. check: max |E e(th_p+th_q-th_r-th_s) - 1[{{p,q}}={{r,s}}]| over depths <= 6 = {worst:.1e}")
pair2 = max(abs(np.mean(e(thetas[:, d] - thetas[:, f]))) for d in range(1, 7) for f in range(1, 7) if d != f)
log(f"    pairwise: max |E e(th_d - th_e)|, d != e, depths <= 6 = {pair2:.1e}")


def fourth_moments(n, qs, vmax=8):
    paths = np.array(list(itertools.product(range(1, vmax + 1), repeat=n)), dtype=np.int64)
    P = len(paths)
    D = np.cumsum(paths[:, ::-1], axis=1)[:, ::-1]      # D[:, i] = a_i + ... + a_n
    w = 2.0 ** (-paths.sum(axis=1))
    W = (w[:, None] * w[None, :]).ravel()
    out = {}
    key = np.zeros((P, P), dtype=np.int64)
    for i in range(n):
        mn = np.minimum(D[:, i][:, None], D[:, i][None, :])
        mx = np.maximum(D[:, i][:, None], D[:, i][None, :])
        key = key * 4096 + (mn * 64 + mx)
    _, inv = np.unique(key.ravel(), return_inverse=True)
    out["iid"] = float(np.sum(np.bincount(inv, weights=W) ** 2)) * 9.0 ** n
    del key, inv
    for q in qs:
        S = 40
        phi = np.zeros(P, dtype=np.int64)
        for i in range(n):
            phi = (phi + (q ** (n - 1 - i)) * (np.int64(1) << (S - D[:, i]))) % (1 << S)
        key = (phi[:, None] + phi[None, :]) % (1 << S)
        _, inv = np.unique(key.ravel(), return_inverse=True)
        out[q] = float(np.sum(np.bincount(inv, weights=W) ** 2)) * 9.0 ** n
        del key, inv
    return out


log("    n   iid      q=3      q=5      q=7      q=9      q=11     q=15")
for n in (1, 2, 3, 4):
    fm = fourth_moments(n, (3, 5, 7, 9, 11, 15))
    log(f"    {n}   " + "  ".join(f"{fm[k]:.3f}" for k in ("iid", 3, 5, 7, 9, 11, 15)))
log(f"  elapsed {time.time() - T0:.1f}s")

# =====================================================================================
# F. Section 4b / E4: the i.i.d.-digit model, E|f_n(k)|^2 = 3^-n; Monte Carlo at n = 10
# =====================================================================================
log("=" * 100)
log("F. i.i.d.-digit model: Monte Carlo of E[3^n |f_n(0)|^2] at n = 10 (A = 40), 2000 seeds; and window positions")
log("=" * 100)
rng = np.random.default_rng(20261004)
nseeds = 2000 if not QUICK else 300
vals0 = np.zeros(nseeds)
valsk = {1: np.zeros(nseeds), 7: np.zeros(nseeds), 40: np.zeros(nseeds)}
A = 40
for s in range(nseeds):
    rec = {10: None}
    la = window_run(1, 3, 12, A, phase_source="iid", rng=rng, record=rec)   # N = 12 so that level 10 keeps a window of 80
    g10 = rec[10]                                                           # rescaled by 3^5
    vals0[s] = abs(g10[0]) ** 2
    for k in valsk:
        valsk[k][s] = abs(g10[k]) ** 2
exact_trunc = ((1 - 4.0 ** (-A)) / 3.0) ** 10 * 3.0 ** 10
log(f"  E[3^10 |f_10(0)|^2] = {vals0.mean():.4f} +- {vals0.std(ddof=1) / math.sqrt(nseeds):.4f} (exact: {exact_trunc:.6f}); median {np.median(vals0):.3f}; max {vals0.max():.1f}")
for k in valsk:
    log(f"  window position k = -{k}: E[3^10 |f_10(k)|^2] = {valsk[k].mean():.4f} +- {valsk[k].std(ddof=1) / math.sqrt(nseeds):.4f}")
log(f"  typical value exp(E log 3^(n/2)|f_10(0)|) = {math.exp(0.5 * np.mean(np.log(vals0))):.3f}")
log(f"  elapsed {time.time() - T0:.1f}s")

# =====================================================================================
# G. Bonus: E2 (moments/tails of the dense law to n = 14) and E3 (FFT of the law at n = 11..13)
# =====================================================================================
log("=" * 100)
log("G. E2/E3 checks from the dense law (A = 40): rho_n = (2/3) 3^n mu_n on the units; moments under the uniform measure on units")
log("=" * 100)
NMAX = 14 if not QUICK else 10
law = np.array([1.0])
prevE2 = None
fft_store = {}
for m in range(1, NMAX + 1):
    qm = 3 ** m
    x = np.arange(len(law), dtype=np.int64)
    y = (3 * x + 1) % qm
    inv2 = pow(2, -1, qm)
    new = np.zeros(qm)
    for a in range(1, 41):
        y = (y * inv2) % qm
        new += 2.0 ** (-a) * np.bincount(y, weights=law, minlength=qm)
    law = new
    units = np.arange(qm) % 3 != 0
    rho = (2.0 / 3.0) * qm * law[units]
    E2 = np.mean(rho ** 2)
    inc = (E2 - prevE2) if prevE2 is not None else float("nan")
    prevE2 = E2
    line = f"  n={m:2d}: E[rho^1.5]={np.mean(rho ** 1.5):.4f} E[rho^1.9]={np.mean(rho ** 1.9):.4f} E[rho^2]={E2:.4f} (inc {inc:+.4f}) E[rho^2.5]={np.mean(rho ** 2.5):.3f} rho(-1)={(2.0 / 3.0) * qm * law[qm - 1]:.3f} rho(1)={(2.0 / 3.0) * qm * law[1]:.3f} max rho={rho.max():.3f} (=({rho.max() / 1.5 ** m:.4f})(3/2)^n)"
    if m >= 9:
        tails = [np.mean(rho > t) for t in (2, 4, 8, 16, 32, 64, 128, 256)]
        ui = [np.mean(rho * (rho > M)) for M in (4, 8, 16, 32, 64, 128)]
        line += f"\n         P(rho>t), t=2..256: {['%.3g' % p for p in tails]};  E[rho 1(rho>M)], M=4..128: {['%.4f' % v for v in ui]}"
    log(line)
    if m in (11, 12, 13):
        fft_store[m] = np.abs(np.fft.fft(law))
for m, F in fft_store.items():
    qm = 3 ** m
    t = np.arange(qm)
    prim = (t % 3 != 0)
    mass_all = np.sum(F[1:] ** 2)
    mass_prim = np.sum(F[prim] ** 2)
    log(f"  n={m}: Parseval total sum_t |mu_hat|^2 = {np.sum(F ** 2):.4f} (= 3^n sum mu^2 = (3/2)E[rho^2] check), primitive-t mass = {mass_prim:.4f}, rms over primitive t = {math.sqrt(mass_prim / prim.sum()):.2e}")
    # top coefficients among primitive t (|t| <= 3^m/2 representatives)
    idx = np.argsort(-F * prim)[:8]
    tops = [(int(i) if i <= qm // 2 else int(i) - qm, round(float(F[i]), 4)) for i in idx]
    log(f"       top primitive coefficients (t, |mu_hat|): {tops}")
    for s_ in (14, 15, 16, 17):
        if m == 13:
            log(f"       |mu_hat_13(+-2^{s_})| = {F[(2 ** s_) % qm]:.4f}, {F[(-2 ** s_) % qm]:.4f}; share of primitive mass {F[(2 ** s_) % qm] ** 2 / mass_prim * 100:.2f}%; theta = {2 ** s_ / qm:.3f}")
    if m == 13:
        log(f"       |mu_hat_13(465905)| = {F[465905]:.4f}  (465905 = -2^16 mod 3^12)")
    theta = t / qm
    bins = np.floor(theta * 64).astype(int)
    bm = np.bincount(bins[prim], weights=F[prim] ** 2, minlength=64) / mass_prim
    log(f"       64-bin mass over theta (primitive t): min {bm.min():.4f}, max {bm.max():.4f}, uniform 1/64 = {1 / 64:.4f}; max/uniform = {bm.max() * 64:.3f}")
    for mm in (3, 5, 8):
        dist = np.abs(theta * 2 ** mm - np.round(theta * 2 ** mm))   # in units of the spacing 2^-mm
        near = dist <= 0.25
        log(f"       mass within a quarter-spacing of j/2^{mm}: {np.sum(F[prim & near] ** 2) / mass_prim:.4f} (uniform 0.5)")
    for T in (100, 1000, 10000):
        sel = prim & ((t <= T) | (t >= qm - T))
        share = np.sum(F[sel] ** 2) / mass_prim
        haar = sel.sum() / prim.sum()
        sel_all = (t != 0) & ((t <= T) | (t >= qm - T))
        share_all = np.sum(F[sel_all] ** 2) / mass_all
        haar_all = sel_all.sum() / (qm - 1)
        log(f"       fixed frequencies 1<=|t|<={T}: primitive-only ratio to Haar share {share / haar:.3f}; all t!=0 (incl. multiples of 3) ratio {share_all / haar_all:.3f}")
log(f"  elapsed {time.time() - T0:.1f}s")

# =====================================================================================
# H. Longer runs: q = 3, u = 1 to N = 1500 (sections 4/4a) and three i.i.d. runs to N = 1500 (section 4/4a)
# =====================================================================================
if not QUICK:
    log("=" * 100)
    log("H. q = 3, u = 1 to N = 2500 and three i.i.d.-digit runs to N = 1500 (A = 40): typical rates, ridges")
    log("=" * 100)
    t1 = time.time()
    la = window_run(1, 3, 2500, 40)
    v = np.exp(la + 0.5 * np.arange(2501) * LN3)
    log(f"  real digits of -3^-n, u=1: rates 200..1500 {ls_rate(la, 200, 1500):.4f} (claimed 0.5691), 300..1500 {ls_rate(la, 300, 1500):.4f} (claimed 0.5689), "
        f"750..1500 {ls_rate(la, 750, 1500):.4f} (claimed 0.5695), 200..2500 {ls_rate(la, 200, 2500):.4f}, 200..800 {ls_rate(la, 200, 800):.4f}  ({time.time() - t1:.0f}s)")
    log(f"     3^(n/2)|mu_hat_n(1)| at n=200: {v[200]:.2e}, n=2500: {v[2500]:.2e} (note: 9e-2 -> 3e-15); largest past n=200: {v[201:].max():.3f} at n={201 + int(np.argmax(v[201:]))} (note: nothing above 0.08);"
        f" count of n in 201..2500 with value > 0.08: {int(np.sum(v[201:] > 0.08))}; mean per-level log increment 300..1500: {(la[1500] - la[300]) / 1200 + 0.5 * LN3:+.4f} (note -0.0133)")
    rng = np.random.default_rng(1)
    for s in range(3):
        t1 = time.time()
        la = window_run(1, 3, 1500, 40, phase_source="iid", rng=rng)
        log(f"  i.i.d. digits seed {s}: rates 300..1500 {ls_rate(la, 300, 1500):.4f}, 750..1500 {ls_rate(la, 750, 1500):.4f}, 200..1500 {ls_rate(la, 200, 1500):.4f} (note: 0.572-0.574; rms rate exactly 0.5774)  ({time.time() - t1:.0f}s)")

    # I. E5q-type ensembles: 2000 random odd 40-bit starts u (prime to q), N = 40: second and fourth moments of f_n(0)
    log("=" * 100)
    log("I. Integer-start ensembles at N = 40 (2000 random odd 40-bit units prime to q; A = 40): 3^n E|f|^2, K_n = 9^n E|f|^4, kurtosis K_n/(3^n E|f|^2)^2")
    log("=" * 100)
    # projectivity: Y_n mod q^(n-1) = Y_(n-1) in law, hence mu_hat_n(q u) = mu_hat_(n-1)(u) exactly (a multiple of q is not a cold unit)
    dev = max(abs(abs(fourier_direct(laws3_exact[n], 3 * u, 3, n)) - abs(fourier_direct(laws3_exact[n - 1], u, 3, n - 1))) for n in (4, 5, 6, 7) for u in (1, 5, 7, 11))
    dev_w = max(abs(math.exp(window_run(3 * u, 3, 7, 40)[7]) - math.exp(window_run(u, 3, 6, 40)[6])) for u in (1, 5, 7))
    dev_w2 = abs(math.exp(window_run(9, 3, 7, 40)[7]) - math.exp(window_run(1, 3, 5, 40)[5]))
    log(f"  projectivity check: max | |mu_hat_n(3u)| - |mu_hat_(n-1)(u)| | (dense law, n<=7) = {dev:.1e}; window recursion u=3u' at N=7 vs u' at N=6: {dev_w:.1e}; u=9 at N=7 vs u=1 at N=5: {dev_w2:.1e}")
    log("  so a start u with v_q(u) = j has 3^n|f_n(0)|^2 = 3^j x (a level-(n-j) unit value): sampling ALL odd 40-bit u (fraction q^-j with v_q = j) inflates the moments.")
    NS = 2000
    levels = (10, 20, 30, 40)
    for src, variant in ((3, "prime to q"), (3, "all odd u"), (5, "prime to q"), (5, "all odd u"), ("iid", "")):
        t1 = time.time()
        rng = np.random.default_rng(7)
        m2 = {n: [] for n in levels}
        m4 = {n: [] for n in levels}
        lg = {n: [] for n in levels}
        nq = 0
        for s in range(NS):
            rec = {n: None for n in levels}
            if src == "iid":
                window_run(1, 3, 40, 40, phase_source="iid", rng=rng, record=rec)
            else:
                while True:
                    uu = int(rng.integers(1 << 39, 1 << 40)) | 1
                    if variant == "all odd u" or uu % src:
                        break
                nq += (uu % src == 0)
                window_run(uu, src, 40, 40, record=rec)
            for n in levels:
                x = abs(rec[n][0]) ** 2      # = 3^n |f_n(0)|^2
                m2[n].append(x)
                m4[n].append(x * x)
                lg[n].append(0.5 * math.log(x) - 0.5 * n * LN3)
        s2 = [np.mean(m2[n]) for n in levels]
        s4 = [np.mean(m4[n]) for n in levels]
        s4e = [np.std(m4[n], ddof=1) / math.sqrt(NS) for n in levels]
        kurt = [s4[i] / s2[i] ** 2 for i in range(len(levels))]
        meanlog = np.array([np.mean(lg[n]) for n in levels])
        typ = math.exp(np.polyfit(np.array(levels), meanlog, 1)[0])
        log(f"  source {src!s:>4} {variant:>10}: 3^n E|f|^2 at n=10,20,30,40 = {['%.2f' % x for x in s2]}; K_n = {['%.0f +- %.0f' % (a, b) for a, b in zip(s4, s4e)]}; "
            f"kurtosis K_n/(3^n E|f|^2)^2 = {['%.0f' % k for k in kurt]}; typical rate 10..40 (fit of E log) {typ:.4f}; starts divisible by q: {nq}  ({time.time() - t1:.0f}s)")
    log("  (note's E5q claims: q=3: 4.19,4.04,2.94,2.86 / 1296,2880,444,920 / rate 0.5656; q=5: 1.16,1.43,1.11,1.65 / 26,244,26,331 / 0.5692; iid: 1.00,1.07,1.04,0.95 / 4.6,18.9,14.6,13.9 / 0.5684)")
log(f"  total elapsed {time.time() - T0:.1f}s")
_outf.close()
