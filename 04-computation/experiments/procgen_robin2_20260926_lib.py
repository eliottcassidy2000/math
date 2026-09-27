#!/usr/bin/env python3
"""procgen_robin2_20260926_lib -- helpers for the robin2 lane (session collatz-procgen-20260922).

Everything is in V-coordinates: V_t = U_t - floor(c t), c = log_3 2, U_t = number of ones among the first t
letters.  delta_t = floor(c(t+1)) - floor(c t) in {0,1} (a Sturmian word of slope c).  One step from site x at
time t goes to x + eps - delta_t (eps = the letter).  S-height (the slope level of the pairpeak/robin notes):
S_t = V_t - {c t}, so for t >= 1:  S_t > 0 <=> V_t >= 1  and  S_t < M <=> V_t <= M.

Exact objects (Python integers):
  N_m(L)   Robin barrier count (sites 1..m, site m = zone, both letters go to m-1, weight 2);
  A_M(L)   hard wall count (V in [1,M] at times 1..L-1, V >= 1 at time L);
  V^{(s)}_k(x), H^{(s)}_k(x), W^{(s)}_k(x)  Dirichlet / half-line / Robin continuation counts;
  P_{s,s+n}(x,1)  confined kernel to site 1;  bridge counts with optional walls.
Analytic helpers (floats with explicit safety margins): the Hoeffding sum G(K,v), KL(p||1/2) bounds, k_b(m).
"""
import math
from fractions import Fraction

C = math.log(2) / math.log(3)
D = 1.0 - C
MU = C - 0.5


# ------------------------------------------------------------------------------------------------ floors
def floors(T):
    """F[t] = floor(c t) exactly (largest n with 3^n <= 2^t), t = 0..T"""
    out = []
    n, q, p2 = 0, 3, 1
    for _ in range(T + 1):
        while q <= p2:
            n += 1
            q *= 3
        out.append(n)
        p2 <<= 1
    return out


def deltas(T):
    F = floors(T + 1)
    return [F[t + 1] - F[t] for t in range(T + 1)]


# ------------------------------------------------------------------------------------------------ N and A
def robin_counts(Lmax, m, dl):
    """N_m(L) for L = 0..Lmax"""
    f = [0] * (m + 2)
    f[0] = 1
    out = [1]
    for t in range(Lmax):
        d = dl[t]
        g = [0] * (m + 2)
        for v in range(m + 1):
            x = f[v]
            if not x:
                continue
            if v == m:
                if d != 1:
                    raise AssertionError("zone site with delta = 0")
                g[m - 1] += 2 * x
            else:
                g[v - d] += x
                g[v + 1 - d] += x
        g[0] = 0
        f = g
        out.append(sum(f))
    return out


def dirichlet_counts(Lmax, M, dl):
    """A_M(L) for L = 0..Lmax"""
    f = [0] * (M + 3)
    f[0] = 1
    out = [1]
    for t in range(Lmax):
        d = dl[t]
        g = [0] * (M + 3)
        for v in range(M + 1):
            x = f[v]
            if x:
                g[v - d] += x
                g[v + 1 - d] += x
        g[0] = 0
        out.append(sum(g[1:]))
        g[M + 1] = 0
        g[M + 2] = 0
        f = g
    return out


# ------------------------------------------------------------------------------------------------ continuation counts
def free_counts(s, x, kmax, top, dl):
    """k -> #words of length k from site x at time s with V in [1, top] at times s+1..s+k-1 and V >= 1 at s+k.
    top = None gives the half-line count H (no top wall)."""
    cap = top if top is not None else x + kmax + 2
    f = [0] * (cap + 3)
    f[x] = 1
    out = [1]
    for i in range(kmax):
        d = dl[s + i]
        g = [0] * (cap + 3)
        for v in range(1, cap + 1):
            y = f[v]
            if y:
                g[v - d] += y
                g[v + 1 - d] += y
        g[0] = 0
        out.append(sum(g[1:]))
        g[cap + 1] = 0
        g[cap + 2] = 0
        f = g
    return out


def robin_from(s, x, kmax, m, dl):
    """k -> Robin weighted count W^{(s)}_k(x) (x in 1..m; x = m only at a zone time)"""
    f = [0] * (m + 2)
    f[x] = 1
    out = [1]
    for i in range(kmax):
        d = dl[s + i]
        g = [0] * (m + 2)
        for v in range(1, m + 1):
            y = f[v]
            if not y:
                continue
            if v == m:
                if d != 1:
                    raise AssertionError("zone site with delta = 0")
                g[m - 1] += 2 * y
            else:
                g[v - d] += y
                g[v + 1 - d] += y
        g[0] = 0
        f = g
        out.append(sum(f))
    return out


def kernel_to_site1(s, xs, n, M, dl):
    """{x: P_{s,s+n}(x,1)}: words from (s,x) to (s+n,1) with V in [1,M] at times s+1..s+n"""
    res = {}
    for x in xs:
        f = [0] * (M + 3)
        f[x] = 1
        for i in range(n):
            d = dl[s + i]
            g = [0] * (M + 3)
            for v in range(1, M + 1):
                y = f[v]
                if y:
                    g[v - d] += y
                    g[v + 1 - d] += y
            g[0] = 0
            g[M + 1] = 0
            g[M + 2] = 0
            f = g
        res[x] = f[1]
    return res


def kernel_to_site1_word(word, xs, M):
    """same as kernel_to_site1 but for an explicit delta-word (a Sturmian factor)"""
    res = {}
    for x in xs:
        f = [0] * (M + 3)
        f[x] = 1
        for d in word:
            g = [0] * (M + 3)
            for v in range(1, M + 1):
                y = f[v]
                if y:
                    g[v - d] += y
                    g[v + 1 - d] += y
            g[0] = 0
            g[M + 1] = 0
            g[M + 2] = 0
            f = g
        res[x] = f[1]
    return res


def bridge_counts(s, x0, x1, n, top, dl, bottom_on=True, top_on=True):
    """words of length n from site x0 (time s) to site x1 (time s+n); walls V>=1 at times s+1..s+n-1 (bottom) and
    V<=top at times s+1..s+n-1 (top), each optional."""
    lo = 1 if bottom_on else -(10 ** 9)
    hi = top if top_on else 10 ** 9
    f = {x0: 1}
    for i in range(n):
        d = dl[s + i]
        g = {}
        for v, y in f.items():
            for w in (v - d, v + 1 - d):
                g[w] = g.get(w, 0) + y
        if i < n - 1:
            g = {v: y for v, y in g.items() if lo <= v <= hi}
        f = g
    return f.get(x1, 0)


# ------------------------------------------------------------------------------------------------ Sturmian factors
def sturmian_factors(n, dl, scan):
    """all distinct factors of length n of delta found among start positions 0..scan-1 (dict factor -> first start)"""
    seen = {}
    for s in range(scan):
        w = tuple(dl[s:s + n])
        if w not in seen:
            seen[w] = s
    return seen


def landing_factors(n, dl, scan):
    """factors delta_s..delta_{s+n-1} with delta_{s-2} delta_{s-1} = 0 1 (s = landing time after a zone time);
    returns (dict, complete) where complete certifies that all n+3 factors of length n+2 were found"""
    allf = sturmian_factors(n + 2, dl, scan)
    complete = (len(allf) == n + 3)
    out = {}
    for w, s in allf.items():
        if w[0] == 0 and w[1] == 1:
            out[w[2:]] = s + 2
    return out, complete


# ------------------------------------------------------------------------------------------------ analytic helpers
def G_upper(K, v, T1=40000):
    """rigorous (up to float rounding, covered by the 1e-9 safety factor) upper bound for
    G(K,v) = sum_{t>=1} exp(-2 (K + v t)^2 / t)"""
    s = 0.0
    for t in range(1, T1 + 1):
        s += math.exp(-2.0 * (K + v * t) ** 2 / t)
    tail = math.exp(-4.0 * K * v - 2.0 * v * v * (T1 + 1)) / (1.0 - math.exp(-2.0 * v * v))
    return (s + tail) * (1.0 + 1e-9)


def KL_half(p):
    """KL(p || 1/2) (natural log)"""
    if p <= 0.0 or p >= 1.0:
        return math.log(2.0)
    return p * math.log(2.0 * p) + (1.0 - p) * math.log(2.0 * (1.0 - p))


def KL_half_upper(y):
    """upper bound 2y^2/(1-4y^2) for KL(1/2+y || 1/2), |y| < 1/2"""
    return 2.0 * y * y / (1.0 - 4.0 * y * y)


def k_b(m):
    """largest k >= 0 with m - 1 - mu k >= sqrt(k)  (closed-form start sqrt(k) = (-1+sqrt(1+4mu(m-1)))/(2mu),
    then exact adjustment; the condition is monotone in k)"""
    r = (-1.0 + math.sqrt(1.0 + 4.0 * MU * (m - 1))) / (2.0 * MU)
    k = max(0, int(r * r) + 2)
    while k > 0 and not (m - 1 - MU * k >= math.sqrt(k)):
        k -= 1
    while m - 1 - MU * (k + 1) >= math.sqrt(k + 1):
        k += 1
    return k


def binom_tail_exact(k, j):
    """P(Bin(k,1/2) >= j) as a Fraction"""
    if j <= 0:
        return Fraction(1)
    if j > k:
        return Fraction(0)
    s = 0
    c = math.comb(k, j)
    for e in range(j, k + 1):
        s += c
        c = c * (k - e) // (e + 1)
    return Fraction(s, 2 ** k)


def binom_tail_lower(k, j):
    """closed-form lower bound for P(Bin(k,1/2) >= j) (Robbins' Stirling bounds):
    j <= k/2: 1/2; else J consecutive terms each >= 0.6753/sqrt(k) exp(-k KL(e/k||1/2)), e <= j+J-1 <= k-1"""
    if j <= k / 2.0:
        return 0.5
    if j > k - 1:
        return 0.0
    J = min(int(math.ceil(math.sqrt(k))), k - j)
    e_last = j + J - 1
    y = e_last / k - 0.5
    return J * 0.6753 / math.sqrt(k) * math.exp(-k * KL_half_upper(y))
