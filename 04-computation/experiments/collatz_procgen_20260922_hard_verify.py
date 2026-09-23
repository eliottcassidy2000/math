#!/usr/bin/env python3
"""Independent exact verifier for the hard-class lane (collatz-procgen-20260922).

Affine 2-adic shift map T(x) = (m_e x + r_e)/2 on the class x = e mod 2, given as
amap = (m0, r0, m1, r1).  For a finite word z, T^|z|(x) = (M_z x + R_z)/2^|z|.

  * phi2_mod(w, N)      Phi_T(w) mod 2^N = -R * M^{-1} mod 2^N, (M, R) of the first N letters,
                        computed by divide and conquer from the closed form
                        R(z1 z2) = M(z2) R(z1) + 2^|z1| R(z2)   (for 3x+1: R = sum_l 3^(a-l) 2^(d_l)).
  * parity_check(X, w)  independent validation: iterate T on X mod 2^(N-t) and compare parities.
  * approximant(w,a,j)  exact Phi_T(U V^inf) = P/Q for U = w[0,a), V = w[a,j):
                        Q = M_U (2^|V| - M_V),  P = 2^|U| R_V - R_U (2^|V| - M_V).
  * verify(...)         reduces P/Q, checks the isometry v_2(Phi_T(w) - P/Q) = lambda exactly,
                        and returns the exact Theorem-R gain lambda - log2(|P'| + |Q'|).
  * wang(X, N)          rational reconstruction (unique candidate of height <= 2^((N-1)/2));
    follow(c, w)        the candidate's parity vector, iterated exactly until it leaves w.
  * phi_real(w)         real value of the same series (converges for strictly supercritical words).
"""
import math
from fractions import Fraction

MAPS = {
    "3x+1": (1, 0, 3, 1),
    "3x-1": (1, 0, 3, -1),
    "5x+1": (1, 0, 5, 1),
    "7x+1": (1, 0, 7, 1),
    "19x+1": (1, 0, 19, 1),
    "mahler": (3, 0, 3, 1),
}


def mu_of(beta, amap):
    """Height exponent mu = max(1, (1-beta) log2 m0 + beta log2 m1) for frequency beta of 1s."""
    m0, _, m1, _ = amap
    return max(1.0, (1 - beta) * math.log2(abs(m0)) + beta * math.log2(abs(m1)))


def _MR(w, lo, hi, amap, cache):
    """(M, R) for the segment w[lo:hi] by divide and conquer (exact integers)."""
    n = hi - lo
    m0, r0, m1, r1 = amap
    if n <= 48:
        M, R = 1, 0
        p2 = 1
        for e in w[lo:hi]:
            if e:
                M, R = m1 * M, m1 * R + r1 * p2
            else:
                M, R = m0 * M, m0 * R + r0 * p2
            p2 <<= 1
        return M, R
    mid = (lo + hi) // 2
    M1, R1 = _MR(w, lo, mid, amap, cache)
    M2, R2 = _MR(w, mid, hi, amap, cache)
    return M1 * M2, M2 * R1 + (R2 << (mid - lo))


def MR(w, lo, hi, amap):
    return _MR(w, lo, hi, amap, None)


def phi2_mod(w, N, amap):
    M, R = MR(w, 0, N, amap)
    mod = 1 << N
    return (-R * pow(M, -1, mod)) % mod


def parity_check(X, w, N, amap):
    """True iff the parity vector of X (mod 2^N) under T equals w[0:N] (exact, T applied mod 2^(N-t))."""
    m0, r0, m1, r1 = amap
    x = X
    for t in range(N):
        b = x & 1
        if b != w[t]:
            return False
        x = (m1 * x + r1) >> 1 if b else (m0 * x + r0) >> 1
        rem = N - t - 1
        x &= (1 << rem) - 1 if rem > 0 else 0
    return True


def approximant(w, a, j, amap):
    Mu, Ru = MR(w, 0, a, amap) if a > 0 else (1, 0)
    Mv, Rv = MR(w, a, j, amap)
    p = j - a
    D = (1 << p) - Mv
    P = (Ru * -D) + (Rv << a)
    Q = Mu * D
    if Q < 0:
        P, Q = -P, -Q
    return P, Q


def v2(n):
    return (n & -n).bit_length() - 1 if n else None


def canonical_ep(w, a, j):
    """Minimal representation (a0, a0 + p0) of the eventually periodic word U V^inf, U = w[:a], V = w[a:j]:
    p0 = primitive period of V, a0 = minimal preperiod.  Same rational Phi_T, possibly smaller P, Q."""
    p = j - a
    V = w[a:j]
    p0 = p
    for d in range(1, p + 1):
        if p % d == 0 and all(V[i] == V[i % d] for i in range(p)):
            p0 = d
            break

    def z(i):
        return w[i] if i < a else w[a + (i - a) % p0]
    a0 = a
    while a0 > 0 and z(a0 - 1) == z(a0 - 1 + p0):
        a0 -= 1
    return a0, a0 + p0


def verify(w, X, Nbits, a, j, lam, amap):
    """Exact certificate for the approximant U = w[:a], V = w[a:j] with claimed common prefix lam."""
    P, Q = approximant(w, a, j, amap)
    g = math.gcd(P, Q)
    Pr, Qr = P // g, Q // g
    iso = None
    if lam < Nbits:
        mod = 1 << Nbits
        diff = (X - Pr * pow(Qr, -1, mod)) % mod
        iso = (v2(diff) == lam)
    # direct check of the common prefix length (independent of the scanner)
    p = j - a
    k = j
    while k < len(w) and w[k] == w[k - p]:
        k += 1
    lam_direct = k
    hbits = math.log2(abs(Pr) + abs(Qr))
    # canonical (minimal) representation of the same EP word: its gcd is the 'genuine' reduction
    a0, j0 = canonical_ep(w, a, j)
    if (a0, j0) != (a, j):
        P0, Q0 = approximant(w, a0, j0, amap)
        g0 = math.gcd(P0, Q0)
        same = (P0 // g0 == Pr and Q0 // g0 == Qr)
    else:
        g0, same = g, True
    # Lemma-H size of the unreduced minimal representation, for comparison
    return dict(a=a, j=j, lam=lam, lam_direct=lam_direct, gain=lam - hbits, hbits=hbits,
                gcd_bits=g.bit_length() - 1 if g > 0 else 0, canon=(a0, j0),
                canon_gcd_bits=g0.bit_length() - 1 if g0 > 0 else 0, canon_same=same, iso=iso,
                Pbits=abs(Pr).bit_length(), Qbits=Qr.bit_length())


def wang(X, N):
    """Wang rational reconstruction mod 2^N with |a|, b <= floor(sqrt(2^(N-1)))."""
    M = 1 << N
    bound = math.isqrt(M >> 1)
    r0, r1 = M, X % M
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    a, b = r1, s1
    if b == 0 or abs(b) > bound:
        return None, bound
    if b < 0:
        a, b = -a, -b
    if b % 2 == 0 or (a - b * X) % M != 0:
        return None, bound
    return Fraction(a, b), bound


def follow(cand, w, amap):
    """Length of the common prefix of the parity vector of cand (odd denominator) and w."""
    m0, r0, m1, r1 = amap
    a, b = cand.numerator, cand.denominator
    for i, e in enumerate(w):
        if (a & 1) != e:
            return i
        a = (m1 * a + r1 * b) >> 1 if e else (m0 * a + r0 * b) >> 1
    return len(w)


def recon_certificate(w, Nbits, amap):
    """('NONE', h, leave_at, cand_bits) or ('NONE', h, None, None) or ('AGREES', h, len(w), bits)."""
    X = phi2_mod(w, Nbits, amap)
    cand, bound = wang(X, Nbits)
    h = bound.bit_length() - 1
    if cand is None:
        return "NONE", h, None, None
    lam = follow(cand, w, amap)
    cb = max(abs(cand.numerator), cand.denominator).bit_length()
    if lam < len(w):
        return "NONE", h, lam, cb
    return "AGREES", h, lam, cb


def phi_real(w, amap, dps=40):
    """Real value -sum_i r_{w_i} 2^i / (m_{w_0} ... m_{w_i}) (strictly supercritical words), mpmath."""
    import mpmath
    mpmath.mp.dps = dps + 10
    m0, r0, m1, r1 = amap
    s = mpmath.mpf(0)
    M = mpmath.mpf(1)
    p2 = mpmath.mpf(1)
    tol = mpmath.mpf(10) ** (-(dps + 5))
    small = 0
    for e in w:
        M *= (m1 if e else m0)
        r = r1 if e else r0
        if r:
            term = r * p2 / M
            s += term
            small = small + 1 if abs(term) < tol else 0
            if small > 50:
                break
        p2 *= 2
    return -s


def E_tail(w, s, amap=(1, 0, 3, 1), dps=30):
    """E_s = -Phi_R(sigma^s w) for the 3x+1 map (positive for words containing a 1)."""
    return -phi_real(w[s:], amap, dps)
