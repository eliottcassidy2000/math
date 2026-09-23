#!/usr/bin/env python3
"""Exact word generators for the hard-class lane (collatz-procgen-20260922).

All words are exact 0/1 lists.  Irrational slopes/angles are quadratic irrationals
x = (a + b*sqrt(d))/c handled with exact integer arithmetic (math.isqrt), so every letter
is the exact letter of the real object (no floating point, no margin checks needed).

Generators
  * mechanical(alpha, rho, n, upper)    Sturmian words s_i = floor((i+1)a+r) - floor(i a + r)
  * rotation_coding(alpha, rho, arcs, n) x_i = 1 iff frac(i*alpha + rho) lies in a union of arcs
  * fixed_point(sigma, a, n)             fixed point of a substitution
  * morphic_image(word, S, prefix)       u S(w)   (quasi-Sturmian words for S(01) != S(10))
  * strip words (rational slope A/B via a finite automaton; quadratic slope via isqrt):
      random / greedy low-repetition / Sturmian-perturbed
"""
import math
import random
from fractions import Fraction


# ---------------------------------------------------------------- quadratic irrationals
class QI:
    """(a + b sqrt(d)) / c with integers, c > 0, d > 1 squarefree (shared per object)."""
    __slots__ = ("a", "b", "c", "d")

    def __init__(self, a, b, c, d):
        if c < 0:
            a, b, c = -a, -b, -c
        g = math.gcd(math.gcd(a, b), c)
        if g > 1:
            a, b, c = a // g, b // g, c // g
        self.a, self.b, self.c, self.d = a, b, c, d

    @staticmethod
    def rat(fr, d):
        fr = Fraction(fr)
        return QI(fr.numerator, 0, fr.denominator, d)

    def __add__(self, o):
        if not isinstance(o, QI):
            o = QI.rat(o, self.d)
        return QI(self.a * o.c + o.a * self.c, self.b * o.c + o.b * self.c, self.c * o.c, self.d)

    def __neg__(self):
        return QI(-self.a, -self.b, self.c, self.d)

    def __sub__(self, o):
        return self + (-o if isinstance(o, QI) else QI.rat(-Fraction(o), self.d))

    def mul_int(self, k):
        return QI(self.a * k, self.b * k, self.c, self.d)

    def sign(self):
        a, b, d = self.a, self.b, self.d
        if a >= 0 and b >= 0:
            return 0 if (a == 0 and b == 0) else 1
        if a <= 0 and b <= 0:
            return -1
        # opposite signs: compare a^2 with b^2 d  (never equal for d non-square, b != 0)
        return (1 if a > 0 else -1) if a * a > b * b * d else (1 if b > 0 else -1)

    def floor(self):
        a, b, c, d = self.a, self.b, self.c, self.d
        if b == 0:
            return a // c
        fb = math.isqrt(b * b * d)
        if b < 0:
            fb = -fb - 1           # floor(b sqrt d) for b<0 (b^2 d is not a square)
        return (a + fb) // c        # floor((n + f)/c) = floor(n/c) for integer n, f in (0,1)

    def frac(self):
        return self - self.floor()

    def __float__(self):
        return (self.a + self.b * math.sqrt(self.d)) / self.c

    def __lt__(self, o):
        return (self - o).sign() < 0

    def __ge__(self, o):
        return (self - o).sign() >= 0

    def __repr__(self):
        return f"({self.a}+{self.b}*sqrt({self.d}))/{self.c}"


def qi_golden_inv():
    """1/phi = (sqrt5 - 1)/2."""
    return QI(-1, 1, 2, 5)


def qi_from_cf(prefix, period, d_hint=None):
    """Exact value of [0; prefix, (period)] as a QI (period nonempty)."""
    # solve t = [period..., t]  ->  t = (P t + P')/(Q t + Q')  via convergent matrices
    p0, q0, p1, q1 = 1, 0, 0, 1        # matrix [[p0, p1],[q0, q1]] = identity
    for a in period:
        p0, p1 = a * p0 + p1, p0
        q0, q1 = a * q0 + q1, q0
    # t = (p0 t + p1)/(q0 t + q1)  ->  q0 t^2 + (q1 - p0) t - p1 = 0, t > 1
    A, B, C = q0, q1 - p0, -p1
    disc = B * B - 4 * A * C
    s = math.isqrt(disc)
    if s * s == disc:
        raise ValueError("rational periodic tail")
    # squarefree part
    d, f = disc, 1
    k = 2
    while k * k <= d:
        while d % (k * k) == 0:
            d //= k * k
            f *= k
        k += 1
    t = QI(-B, f, 2 * A, d)            # positive root (-B + sqrt(disc)) / (2A)
    # x = [0; prefix, t] : fold from the back
    x = t
    for a in reversed(prefix):
        x = inv_qi(x) + a
    return inv_qi(x)


def inv_qi(x):
    """1/x for a QI x = (a + b sqrt d)/c : c (a - b sqrt d)/(a^2 - b^2 d)."""
    a, b, c, d = x.a, x.b, x.c, x.d
    den = a * a - b * b * d
    return QI(c * a, -c * b, den, d)


def cf_terms_qi(x, n):
    """First n partial quotients of a QI in (0,1) (x = [0; a1, a2, ...]); returns [0, a1, ...]."""
    out = []
    for _ in range(n):
        f = x.floor()
        out.append(f)
        y = x - f
        if y.a == 0 and y.b == 0:
            break
        x = inv_qi(y)
    return out


def convergent_dens(terms):
    qs = []
    qm, q = 0, 1
    for a in terms[1:]:
        qm, q = q, a * q + qm
        qs.append(q)
    return qs


# ---------------------------------------------------------------- Sturmian / rotation words
def mechanical(alpha, rho, n, upper=False, start=0):
    """s_i = floor((i+1)alpha + rho) - floor(i alpha + rho), i = start..start+n-1 (ceil if upper)."""
    def fl(i):
        y = alpha.mul_int(i) + rho
        return -((-y).floor()) if upper else y.floor()
    out = []
    prev = fl(start)
    for i in range(start, start + n):
        cur = fl(i + 1)
        out.append(cur - prev)
        prev = cur
    return out


def rotation_coding(alpha, rho, arcs, n):
    """x_i = 1 iff frac(i alpha + rho) lies in one of the half-open arcs [lo, hi) (lo < hi in [0,1),
    or wrapping arcs given with hi < lo meaning [lo,1) u [0,hi)).  Endpoints are QI or Fractions."""
    d = alpha.d
    arcs_q = [(lo if isinstance(lo, QI) else QI.rat(lo, d), hi if isinstance(hi, QI) else QI.rat(hi, d)) for lo, hi in arcs]
    out = []
    y = rho if isinstance(rho, QI) else QI.rat(rho, d)
    for i in range(n):
        f = y.frac()
        hit = 0
        for lo, hi in arcs_q:
            if lo < hi:
                if f >= lo and f < hi:
                    hit = 1
                    break
            else:
                if f >= lo or f < hi:
                    hit = 1
                    break
        out.append(hit)
        y = y + alpha
    return out


def fixed_point(sigma, a, n):
    """Fixed point sigma^omega(a) (sigma: dict letter -> list), first n letters."""
    w = [a]
    while len(w) < n:
        nw = []
        for c in w:
            nw.extend(sigma[c])
            if len(nw) >= n:
                break
        if nw[:len(w)] != w[:min(len(w), len(nw))]:
            raise ValueError("not prolongable on this letter")
        w = nw
    return w[:n]


def pf_data(sigma):
    """Perron-Frobenius eigenvalue, left vector (lengths) and right vector (frequencies) for a
    binary substitution, from the 2x2 incidence matrix M[i][j] = |sigma(j)|_i."""
    M = [[sigma[j].count(i) for j in (0, 1)] for i in (0, 1)]
    tr = M[0][0] + M[1][1]
    det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    th = (tr + math.sqrt(tr * tr - 4 * det)) / 2
    th2 = (tr - math.sqrt(tr * tr - 4 * det)) / 2
    # right eigenvector f: M f = th f
    if M[0][1] != 0:
        f = [M[0][1], th - M[0][0]]
    else:
        f = [th - M[1][1], M[1][0]]
    s = f[0] + f[1]
    f = [f[0] / s, f[1] / s]
    # left eigenvector l: l M = th l  ->  l0 M00 + l1 M10 = th l0
    if M[1][0] != 0:
        l = [M[1][0], th - M[0][0]]
    else:
        l = [th - M[1][1], M[0][1]]
    return th, th2, l, f, M


def morphic_image(word, S, prefix=()):
    out = list(prefix)
    for c in word:
        out.extend(S[c])
    return out


# ---------------------------------------------------------------- strip words
def strip_automaton(A, B, W):
    """States k = B a_s - A s in [0, floor(B W)]; letter 1: k += B - A, letter 0: k -= A."""
    K = int(math.floor(B * W + 1e-12))
    trans = {}
    for k in range(K + 1):
        t = {}
        if 0 <= k + B - A <= K:
            t[1] = k + B - A
        if 0 <= k - A <= K:
            t[0] = k - A
        trans[k] = t
    return K, trans


def strip_matrix(A, B, W, weights=None):
    """Adjacency matrix of the strip automaton, optionally with column weights (target measures)."""
    import numpy as np
    K, trans = strip_automaton(A, B, W)
    M = np.zeros((K + 1, K + 1))
    for k, t in trans.items():
        for e, k2 in t.items():
            M[k, k2] += 1.0 if weights is None else weights[k2]
    return M


def strip_entropy(A, B, W, weights=None):
    """log2 of the spectral radius of the (weighted) strip automaton (exact eigenvalues)."""
    import numpy as np
    M = strip_matrix(A, B, W, weights)
    r = max(abs(np.linalg.eigvals(M)))
    return math.log2(r) if r > 0 else float("-inf")


def strip_word_rational(A, B, W, n, mode, seed=0, k0=None, pmax=2000):
    """Strip word of slope A/B, width W: modes 'random' (fair choice when both letters allowed),
    'greedy' (at a free choice, pick the letter minimising the best prefix-repetition ratio
    (j+LCE)/j available at the current position, periods <= pmax)."""
    K, trans = strip_automaton(A, B, W)
    rnd = random.Random(seed)
    k = K // 2 if k0 is None else k0
    w = []
    # greedy bookkeeping: for each period p, start index of the current agreement run w[i]=w[i-p]
    if mode == "greedy":
        import numpy as np
        runstart = np.zeros(pmax + 1, dtype=np.int64)   # run of w[i]==w[i-p] for i in [runstart, len)
        ps = np.arange(pmax + 1)
    for i in range(n):
        t = trans[k]
        if not t:
            raise ValueError("dead state")
        if len(t) == 1:
            e = next(iter(t))
        elif mode == "random":
            e = rnd.randrange(2)
        elif mode == "greedy":
            import numpy as np
            best_e, best_val = None, None
            L = len(w)
            for e_try in (0, 1):
                # periods p <= min(L, pmax): agreement continues iff w[L-p] == e_try
                if L == 0:
                    val = 0.0
                else:
                    pm = min(L, pmax)
                    prev = np.frombuffer(bytes(w[L - pm:L]), dtype=np.uint8)[::-1]  # prev[p-1] = w[L-p]
                    agree = (prev == e_try)
                    rs = runstart[1:pm + 1]
                    # after appending, run for period p (if agree) covers [rs, L]; the approximant
                    # U = w[0, rs-p), V = w[rs-p, rs) has |UV| = rs and lambda >= L+1
                    jj = np.where(agree, np.maximum(rs, ps[1:pm + 1]), 10 ** 12)
                    val = float(np.max((L + 1) / jj)) if agree.any() else 0.0
                if best_val is None or val < best_val or (val == best_val and rnd.random() < 0.5):
                    best_e, best_val = e_try, val
            e = best_e
        else:
            raise ValueError(mode)
        if mode == "greedy":
            import numpy as np
            L = len(w)
            pm = min(L, pmax)
            if pm > 0:
                prev = np.frombuffer(bytes(w[L - pm:L]), dtype=np.uint8)[::-1]
                agree = (prev == e)
                rs = runstart[1:pm + 1]
                rs[~agree] = L + 1
            if L + 1 <= pmax:
                runstart[L + 1] = L + 1     # period p = L+1 has no comparison yet
        w.append(e)
        k = t[e]
    return w


def strip_check_rational(w, A, B, W, k0):
    K = int(math.floor(B * W + 1e-12))
    k = k0
    for e in w:
        k += (B - A) if e else -A
        if not (0 <= k <= K):
            return False
    return True


def strip_word_quadratic(alpha, c, W, n, mode, seed=0, flip_rate=0.05):
    """Strip around a quadratic slope alpha (QI): a_s in [alpha s + c, alpha s + c + W] (c, W rational).
    modes: 'random' (fair choice when both letters allowed),
           'perturbed' (follow the lower mechanical word floor(alpha s + c + W') while allowed, and at a
           free choice flip away from it with probability flip_rate)."""
    rnd = random.Random(seed)
    d = alpha.d
    c = Fraction(c)
    W = Fraction(W)

    def lo_ok(a, s):     # a >= alpha s + c
        return (QI.rat(a - c, d) - alpha.mul_int(s)).sign() >= 0

    def hi_ok(a, s):     # a <= alpha s + c + W
        return (alpha.mul_int(s) + (c + W) - a).sign() >= 0
    # start: smallest a_0 = 0 must satisfy the strip at s = 0: need c <= 0 <= c + W
    assert c <= 0 <= c + W
    a = 0
    w = []
    mid = c + W / 2
    for s in range(n):
        opts = [e for e in (0, 1) if lo_ok(a + e, s + 1) and hi_ok(a + e, s + 1)]
        if not opts:
            raise ValueError("dead")
        if len(opts) == 1:
            e = opts[0]
        elif mode == "random":
            e = rnd.randrange(2)
        elif mode == "perturbed":
            # reference: the mechanical word through the middle of the strip
            ref = (alpha.mul_int(s + 1) + mid).floor() - (alpha.mul_int(s) + mid).floor()
            e = ref if rnd.random() >= flip_rate else 1 - ref
        else:
            raise ValueError(mode)
        w.append(e)
        a += e
    return w


def discrepancy_range(w, alpha_float):
    a, lo, hi = 0, 0.0, 0.0
    for s, e in enumerate(w, 1):
        a += e
        dlt = a - alpha_float * s
        lo, hi = min(lo, dlt), max(hi, dlt)
    return lo, hi


def complexity(w, nmax, step=1):
    out = []
    s = bytes(w)
    for n in range(1, nmax + 1, step):
        out.append((n, len({s[i:i + n] for i in range(len(s) - n + 1)})))
    return out


def write_word(path, w):
    with open(path, "w") as fh:
        fh.write("".join("1" if c else "0" for c in w))
