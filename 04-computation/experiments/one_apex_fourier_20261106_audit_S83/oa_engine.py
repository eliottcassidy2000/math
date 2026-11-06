# Independent engine for the S83 audit (does not reuse the session's code).
# g(v) is computed by running the (H, A, p) automaton over the bits of v (MSB first), not by the doubling recursion.
# Phases e(theta*h) are computed in 64-bit fixed point: T = round(theta * 2^64) mod 2^64, frac(theta*h) ~ (T*h mod 2^64)/2^64,
# error <= h * 2^-65 (< 2^-36 for h <= 2^28), far below float64 theta*h (which loses ~28 bits at h ~ 2^28).
import numpy as np, numba as nb, math
from fractions import Fraction

@nb.njit(cache=True)
def g_auto(v):
    # automaton: (H, A, p) --1--> (H, 2A+1, 1); --0--> (H + A + p, 2A, 0); start (1, 0, 0); leading zeros harmless
    if v == 0:
        return 1
    nbits = 0
    t = v
    while t:
        nbits += 1
        t >>= 1
    H = 1; A = 0; p = 0
    for j in range(nbits - 1, -1, -1):
        b = (v >> j) & 1
        if b == 1:
            A = 2 * A + 1; p = 1
        else:
            H = H + A + p; A = 2 * A; p = 0
    return H

@nb.njit(parallel=True, cache=True)
def g_block(n):
    # H(x) for x in X_n, i.e. g(v) for v in [2^(n-1), 2^n), as int64
    N = 1 << (n - 1)
    out = np.empty(N, dtype=np.int64)
    for i in nb.prange(N):
        out[i] = g_auto(N + i)
    return out

@nb.njit(parallel=True, cache=True)
def expsum_fixed(vals, T, nchunks):
    # sum_i e(theta * vals[i]) with theta = T / 2^64 (uint64), chunked partial sums
    N = vals.size
    re = np.zeros(nchunks); im = np.zeros(nchunks)
    step = (N + nchunks - 1) // nchunks
    for c in nb.prange(nchunks):
        lo = c * step; hi = min(N, lo + step)
        sr = 0.0; si = 0.0
        for i in range(lo, hi):
            ph = np.uint64(T) * np.uint64(vals[i])           # wraps mod 2^64
            x = float(ph) * 5.421010862427522e-20 * 2.0 * math.pi  # 2^-64 * 2 pi
            sr += math.cos(x); si += math.sin(x)
        re[c] = sr; im[c] = si
    return re.sum(), im.sum()

def theta_fixed(theta):
    """theta as uint64 fixed point (fractional part only). Accepts float, Fraction, or mpmath mpf string via Fraction."""
    if isinstance(theta, Fraction):
        fr = theta - math.floor(theta)
        return int((fr.numerator << 64) // fr.denominator) % (1 << 64), float(fr)
    raise TypeError

def S_exact(vals, theta_frac):
    """|S_n(theta)| for theta given as a Fraction (exact rational or a 2^-80-accurate rational approximation of an irrational)."""
    T, _ = theta_fixed(theta_frac)
    re, im = expsum_fixed(vals, np.uint64(T), 256)
    return complex(re, im)

@nb.njit(cache=True)
def count_values(vals, maxval):
    c = np.zeros(maxval + 1, dtype=np.int32)
    for i in range(vals.size):
        c[vals[i]] += 1
    return c

@nb.njit(cache=True)
def _sq_stats(c):
    E = 0; D = 0; mx = 0
    for i in range(c.size):
        x = np.int64(c[i])
        E += x * x
        if x > 0:
            D += 1
        if x > mx:
            mx = x
    return E, D, mx

def energy(vals, n):
    c = count_values(vals, (1 << n) + 1)
    E, D, mx = _sq_stats(c)
    del c
    return int(E), int(D), int(mx)

# high-precision rational approximations of the five irrationals (error < 2^-100)
import mpmath
mpmath.mp.prec = 200
def frac_of(x):
    # x: mpmath mpf; return Fraction with |x - F| < 2^-120
    return Fraction(int(mpmath.floor(x * (mpmath.mpf(2) ** 130))), 1 << 130)
IRR = {
    'phi': frac_of((mpmath.sqrt(5) - 1) / 2),
    'sqrt2-1': frac_of(mpmath.sqrt(2) - 1),
    'pi-3': frac_of(mpmath.pi - 3),
    'e-2': frac_of(mpmath.e - 2),
    'sqrt3-1': frac_of(mpmath.sqrt(3) - 1),
}

def bound_A(n, q):
    """Theorem A bound / 2^(n-1), minimised over all splits n = a + L + l (a >= 2, L >= 1, l >= 1)."""
    best = float('inf'); arg = None
    for L in range(1, n - 2):
        for a in range(2, n - L):
            l = n - a - L
            if l < 1:
                continue
            R = 2 * (2.0 ** L / q + 1) * math.sqrt((3 + 2.0 ** (3 - a) * q) * (3 + 2.0 ** (2 - l) * q))
            b = 2.0 ** (-L / 2) * math.sqrt(1 + R)
            if b < best:
                best = b; arg = (a, L, l)
    return best, arg
