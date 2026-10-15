#!/usr/bin/env python3
"""Audit C (2026-10-14) -- independent helpers for auditing THM-4653 / HYP-9276.

Written from scratch for the audit; shares no code with 04-computation/experiments/paley_k6p6_20261014/.

Fields
------
* Prime fields: chi by Euler's criterion via a discrete-log table of a primitive root (not via squares).
* Prime-power fields GF(p^r): elements are base-p digit vectors; the defining polynomial is the LAST monic
  polynomial (in reverse lexicographic order of the tail) for which the class of x generates a cyclic
  unit group of order q-1.  Finding such a generator PROVES F_p[x]/(f) is a field (every nonzero element is a
  power of x, hence a unit).  chi(x) = (-1)^log(x).  This differs from the audited script (first irreducible
  polynomial, chi from a multiplication table of squares).

M_q evaluators (M_q = - sum_{a,x} chi(a) chi(x) chi(x-1) chi(x-a) t(a) t(x), t(a) = sum_y chi(y)chi(y-1)chi(y-a))
-------------
* M_exact_direct(q): exact int64 O(q^2) via numpy.correlate / numpy.convolve (direct, no FFT), prime q.
* M_mult_fft(q): multiplicative-group FFT (discrete logs; cyclic correlations on Z/(q-1)) using
      t(x) = sum_y phi(y) phi(1-y) phi(1-xy)         (x != 0)
      h(x) = sum_y phi(y) phi(1-y) t(xy)             (x != 0)
      M    = - sum_x phi(x) phi(x-1) t(x) h(x)
  i.e. the *second* (multiplicative) forms of THM-4653 (3)-(4), not the additive FFT of the audited script.
  Rounding is checked: max distance of every rounded integer array from its float values is returned.
* M_field(F): exact O(q^2) evaluation over any GF(p^r) object F (below).
"""
import numpy as np


# ----------------------------------------------------------------------------------------------
# prime fields
# ----------------------------------------------------------------------------------------------
def factor_small(n):
    f, d = [], 2
    while d * d <= n:
        if n % d == 0:
            f.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        f.append(n)
    return f


def primitive_root(p):
    if p == 2:
        return 1
    fs = factor_small(p - 1)
    for g in range(2, p):
        if all(pow(g, (p - 1) // f, p) != 1 for f in fs):
            return g
    raise ValueError(p)


def dlog_tables(p):
    """exp[k] = g^k mod p (k = 0..p-2), log[x] for x != 0 (log[0] = -1)."""
    g = primitive_root(p)
    n = p - 1
    # vectorised powers: g^(i*B + j) = (g^B)^i * g^j
    B = int(np.ceil(np.sqrt(n)))
    small = np.empty(B, dtype=np.int64)
    small[0] = 1
    for j in range(1, B):
        small[j] = (small[j - 1] * g) % p
    gB = pow(g, B, p)
    big = np.empty(B, dtype=np.int64)
    big[0] = 1
    for i in range(1, B):
        big[i] = (big[i - 1] * gB) % p
    ex = ((big[:, None] * small[None, :]) % p).reshape(-1)[:n]
    lg = np.full(p, -1, dtype=np.int64)
    lg[ex] = np.arange(n)
    assert (lg[1:] >= 0).all() and len(set(ex.tolist())) == n
    return ex, lg


def chi_prime(p):
    """chi on Z/p as int64 array, via parity of discrete logarithm (independent of 'is a square')."""
    ex, lg = dlog_tables(p)
    c = np.zeros(p, dtype=np.int64)
    c[1:] = 1 - 2 * (lg[1:] % 2)
    return c


def M_exact_direct(p, c=None):
    """Exact M_p for a prime p by direct (non-FFT) integer correlation/convolution. Returns (M, t)."""
    if c is None:
        c = chi_prime(p)
    x = np.arange(p)
    g = c * c[(x - 1) % p]                       # g(y) = chi(y) chi(y-1)
    # t(a) = sum_y g(y) chi(y - a) = sum_z g(z + a) chi(z)
    gg = np.concatenate([g, g])
    t = np.correlate(gg, c, mode='valid')[:p]     # t[a] = sum_z gg[a + z] c[z]
    u = c * t                                     # u(a) = chi(a) t(a)
    conv = np.convolve(u, c)                      # linear convolution, length 2p - 1
    h = conv[:p].copy()
    h[:p - 1] += conv[p:]                         # cyclic fold: h(x) = sum_a u(a) chi(x - a)
    M = -int(np.sum(g * t * h, dtype=np.int64))
    return M, t, h


def M_mult_fft(p, c=None, tables=None):
    """M_p via multiplicative FFTs on Z/(p-1). Returns (M, maxdev) where maxdev is the worst rounding gap."""
    if tables is None:
        ex, lg = dlog_tables(p)
    else:
        ex, lg = tables
    if c is None:
        c = np.zeros(p, dtype=np.int64)
        c[1:] = 1 - 2 * (lg[1:] % 2)
    n = p - 1
    # k(y) = phi(y) phi(1-y) on F_p^*, indexed by log y; l(w) = phi(1-w) indexed by log w
    k = (c[ex] * c[(1 - ex) % p]).astype(np.float64)
    l = c[(1 - ex) % p].astype(np.float64)
    # t(x) = sum_y k(y) l(x y): with x = g^i, y = g^j: t_i = sum_j k_j l_{i+j}  (cross-correlation)
    K = np.fft.rfft(k)
    tf = np.fft.irfft(np.conj(K) * np.fft.rfft(l), n)
    t = np.rint(tf)
    dev1 = float(np.max(np.abs(tf - t)))
    # h(x) = sum_y k(y) t(x y)
    hf = np.fft.irfft(np.conj(K) * np.fft.rfft(t), n)
    h = np.rint(hf)
    dev2 = float(np.max(np.abs(hf - h)))
    # M = - sum_{x != 0} phi(x) phi(x-1) t(x) h(x), x = g^i
    w = c[ex] * c[(ex - 1) % p]
    ti = t.astype(np.int64)
    hi = h.astype(np.int64)
    M = -int(np.sum(w * ti * hi, dtype=np.int64))
    return M, max(dev1, dev2), ti, hi


def W_from_M(q, eps, M):
    """THM-4653 statement 2: W_q = -eps M_q + 7 + 13 eps - 2 eps q."""
    return -eps * M + 7 + 13 * eps - 2 * eps * q


# ----------------------------------------------------------------------------------------------
# prime-power fields
# ----------------------------------------------------------------------------------------------
class GF:
    """GF(p^r) with elements 0..q-1 (base-p digits, low degree first). Builds log/antilog tables."""

    def __init__(self, p, r):
        import itertools
        self.p, self.r, self.q = p, r, p ** r
        q = self.q
        pw = p ** np.arange(r)
        self.digits = np.array([[(e // p ** k) % p for k in range(r)] for e in range(q)], dtype=np.int64)
        self.pw = pw
        tails = list(itertools.product(range(p), repeat=r))[::-1]     # reverse order: differs from audited code
        found = None
        for tail in tails:
            if tail[-1] == 0:          # constant term 0 -> x is a zero divisor
                continue
            # f = x^r + tail[0] x^(r-1) + ... + tail[r-1]; reduction x^r = -(tail[0] x^(r-1) + ... + tail[r-1])
            red = [(-tail[r - 1 - j]) % p for j in range(r)]           # coefficient of x^j in x^r (low first)
            # powers of x: start at 1, multiply by x repeatedly
            cur = [1] + [0] * (r - 1)
            seen_one_at = None
            ex = [0] * (q - 1)
            ok = True
            for kk in range(q - 1):
                code = sum(cur[j] * p ** j for j in range(r))
                ex[kk] = code
                if kk > 0 and code == 1:
                    ok = False
                    break
                top = cur[r - 1]
                nxt = [0] + cur[:r - 1]
                for j in range(r):
                    nxt[j] = (nxt[j] + top * red[j]) % p
                cur = nxt
            if not ok:
                continue
            if sum(cur[j] * p ** j for j in range(r)) != 1:            # x^(q-1) must be 1
                continue
            if len(set(ex)) != q - 1 or 0 in ex:
                continue
            found = (tail, ex)
            break
        assert found is not None
        self.poly_tail = found[0]
        self.exp = np.array(found[1], dtype=np.int64)                  # exp[k] = x^k
        self.log = np.full(q, -1, dtype=np.int64)
        self.log[self.exp] = np.arange(q - 1)
        self.chi = np.zeros(q, dtype=np.int64)
        self.chi[self.exp] = 1 - 2 * (np.arange(q - 1) % 2)
        self.one = 1
        self.neg_one = self.neg(1)
        self.eps = int(self.chi[self.neg_one])

    def code(self, digs):
        return (np.asarray(digs) % self.p) @ self.pw

    def sub_vec(self, xs, a):
        """vector of (xs - a) for an index array xs and scalar a."""
        return ((self.digits[xs] - self.digits[a]) % self.p) @ self.pw

    def sub_scalar(self, a, b):
        return int(((self.digits[a] - self.digits[b]) % self.p) @ self.pw)

    def neg(self, a):
        return int(((-self.digits[a]) % self.p) @ self.pw)

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        return int(self.exp[(self.log[a] + self.log[b]) % (self.q - 1)])


def M_field(F):
    """Exact M_q over GF(p^r) object F (O(q^2)). Returns (M, t)."""
    q = F.q
    c = F.chi
    allx = np.arange(q)
    xm1 = F.sub_vec(allx, 1)
    g = c * c[xm1]                                   # chi(y) chi(y-1)
    t = np.empty(q, dtype=np.int64)
    for a in range(q):
        t[a] = int(np.dot(g, c[F.sub_vec(allx, a)]))
    gt = g * t
    M = 0
    for a in range(q):
        if c[a] == 0:
            continue
        M += int(c[a] * t[a] * np.dot(gt, c[F.sub_vec(allx, a)]))
    return -M, t
