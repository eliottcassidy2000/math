#!/usr/bin/env python3
"""Pascal boundary of counter-only Collatz flows, prior constraints, ascending
shadows, the leaf-section tower, finite heads, and a residue-level mass bound.

Companion note:
  05-knowledge/results/collatz_pascal_boundary_leaf_section_20261005.md

Conventions (Codex adaptive note, 2026-10-05): U(n)=oddpart(3n+1), S(n)=4n+1,
bases b have v2(3b+1) in {1,2}; a rooted n>1 with ROOT word (a_1..a_tau) has
counters L=tau-1 (base edges) and K=sum floor((a_i-1)/2) (sibling total).
  W(n)=2 K!(L+1)!/(L+K+2)!      beta(1,2) mixture (Codex P2)
  E(n)=6 (K+1)!(L+1)!/(L+K+3)!  beta(2,2) mixture (Codex B1; M=K, T=L)

No universal convergence assumption is used. Finite orbit checks certify only
the listed finite universe. All load-bearing identities use Fractions; the
residue-level bound is float64 with explicit margins (VERIFIED numerical).
Checks survive python -O (they use require, not assert).

Usage: python3 <this file> [--level-max 6] [--head-bits 20] [--json PATH]
"""
from fractions import Fraction as F
from math import comb, factorial, lgamma, exp, log
import argparse
import json
import sys
import time

import numpy as np

try:
    from numba import njit
except Exception:  # pragma: no cover
    def njit(*a, **k):
        def wrap(f):
            return f
        return wrap if not a or not callable(a[0]) else a[0]

CHECKS = 0


def require(cond, witness=None):
    global CHECKS
    CHECKS += 1
    if not cond:
        raise RuntimeError(witness)


def v2(n):
    if n == 0:
        raise ValueError(n)
    return (n & -n).bit_length() - 1


def U(n):
    t = 3 * n + 1
    return t >> v2(t)


def root_counters(n):
    """(L, K) from the actual ROOT word of a positive odd n; forward iteration."""
    if n == 1:
        return 0, 0
    L = -1
    K = 0
    while n != 1:
        t = 3 * n + 1
        a = v2(t)
        K += (a - 1) // 2
        L += 1
        n = t >> a
    return L, K


def wW(L, K):
    return F(2 * factorial(K) * factorial(L + 1), factorial(L + K + 2))


def wE(L, K):
    return F(6 * factorial(K + 1) * factorial(L + 1), factorial(L + K + 3))


def beta_frac(a, b):
    """B(a,b) for positive integers."""
    return F(factorial(a - 1) * factorial(b - 1), factorial(a + b - 1))


def polya_prob(M, T, red0=2, blue0=2):
    """Polya urn probability of a specific colour sequence with M red, T blue."""
    p = F(1)
    for i in range(M):
        p *= F(red0 + i, red0 + blue0 + i)
    for i in range(T):
        p *= F(blue0 + i, red0 + blue0 + M + i)
    return p


# ---------------------------------------------------------------------------
# S1. Decode the snippet
# ---------------------------------------------------------------------------
def section_decode():
    print("== S1. Decode: E = 6 B(M+2,T+2) is the Polya(2,2) sequence law ==")
    for M in range(0, 13):
        for T in range(0, 13):
            e = wE(T, M)  # L=T base edges, K=M sibling total
            require(e == 6 * beta_frac(M + 2, T + 2), (M, T))
            require(e == polya_prob(M, T, 2, 2), ("polya", M, T))
            require(e == wE(T + 1, M) + wE(T, M + 1), ("split", M, T))
            w = wW(T, M)
            require(w == 2 * beta_frac(M + 1, T + 2), ("W", M, T))
            require(w == polya_prob(M, T, 1, 2), ("polyaW", M, T))
    # 11 = 6*(1+1/2+1/3) = 3!*H_3 ; with P1 the beta(2,2) bound is 17/2
    H3 = F(1) + F(1, 2) + F(1, 3)
    require(6 * H3 == 11)
    # integrate 6 r(1-r) * (bound)/(1-r) exactly for polynomial bounds
    def int_poly(coeffs):  # integral_0^1 sum c_i r^i dr
        return sum(F(c, i + 1) for i, c in enumerate(coeffs))
    # Codex Green bound D/r -> 6 int r*(1+r+r^2)/r = 6 int (1+r+r^2) = 11
    require(6 * int_poly([1, 1, 1]) == 11)
    # Codex P1 bound 2+2r-r^2 -> 6 int r(2+2r-r^2) = 17/2
    require(6 * int_poly([0, 2, 2, -1]) == F(17, 2))
    # beta(1,2): 2 int (2+2r-r^2) = 16/3
    require(2 * int_poly([2, 2, -1]) == F(16, 3))
    # root rays: E(S^j 1)=6/((j+2)(j+3)) sums to 3 ; W(S^j 1)=2/((j+1)(j+2)) sums to 2
    for j in range(0, 60):
        require(wE(0, j) == F(6, (j + 2) * (j + 3)), j)
        require(wW(0, j) == F(2, (j + 1) * (j + 2)), j)
    # leaf (3|n) masses = int r/(1-r) dmu : beta(2,2) -> int 6 r^2 = 2 ; beta(1,2) -> int 2r = 1
    require(6 * int_poly([0, 0, 1]) == 2)
    require(2 * int_poly([0, 1]) == 1)
    # orders modulo 9 (the structural 6 and 3 behind 22 and 11)
    def order_mod(g, m):
        x, k = g % m, 1
        while x != 1:
            x = (x * g) % m
            k += 1
        return k
    require(order_mod(2, 9) == 6 and order_mod(4, 9) == 3 and order_mod(8, 9) == 2)
    require((1 << 6) // 3 == 21 and F(64, 3) < 22)
    print("  E(T,M)=6B(M+2,T+2)=Polya(2,2) law; split identity; 6*H_3=11;"
          " P1-integrated beta(2,2) bound 17/2; beta(1,2) bound 16/3")
    print("  root rays 3 (E) and 2 (W); leaf masses 2 (E) and 1 (W);"
          " ord_9(2)=6, ord_9(4)=3; 2^6/3 = 21.33 < 22")
    # the 27 control from the Codex note
    L27, K27 = root_counters(27)
    require(wE(L27, K27) == F(1, 21296186450), (L27, K27))
    require((L27, K27) == (40, 8), (L27, K27))
    require(root_counters(53) == (1, 3) and root_counters(113) == (1, 3))
    require(root_counters(111) == (23, 7) and root_counters(155) == (29, 8))
    print("  controls: E(27)=1/21296186450 at (L,K)=(40,8); 53,113 -> (1,3);"
          " 111 -> (23,7); 155 -> (29,8)")


# ---------------------------------------------------------------------------
# S2. Classification mechanism (Theorem A): split + positivity <-> moments
# ---------------------------------------------------------------------------
def section_classification():
    print("== S2. Theorem A mechanism: split arrays are Hausdorff moment arrays ==")
    rng = np.random.default_rng(20261005)
    for trial in range(6):
        k = int(rng.integers(1, 5))
        rs = [F(int(rng.integers(1, 97)), 97) for _ in range(k)]
        ps = [F(int(rng.integers(1, 20))) for _ in range(k)]
        tot = sum(ps)
        ps = [p / tot for p in ps]
        Lmax, Kmax = 10, 14

        def w(L, K):
            return sum(p * r ** K * (1 - r) ** L for p, r in zip(ps, rs))

        require(w(0, 0) == 1)
        for L in range(Lmax):
            for K in range(Kmax):
                require(w(L, K) > 0)
                require(w(L, K) == w(L + 1, K) + w(L, K + 1), (trial, L, K))
        # (-Delta)^L m_K with m_K=w(0,K) reproduces w(L,K)
        m = [w(0, K) for K in range(Kmax + Lmax + 2)]
        cur = list(m)
        for L in range(Lmax):
            for K in range(Kmax):
                require(cur[K] == w(L, K), ("difference", trial, L, K))
            cur = [cur[K] - cur[K + 1] for K in range(len(cur) - 1)]
        # exact fibre payment with remainder (Corollary A1)
        for L in range(0, 4):
            for K in range(0, 4):
                for J in range(1, 8):
                    s = sum(w(L + 1, K + j) for j in range(J))
                    require(s == w(L, K) - w(L, K + J), ("fibre", trial, L, K, J))
    # a split array that is NOT a moment array must have a negative entry
    m = [F(1), F(9, 10), F(1, 2)]
    require(m[0] - 2 * m[1] + m[2] < 0)  # (-Delta)^2 m_0 = w(2,0) < 0
    # an atom at r=1 leaves a root-row defect equal to the atom mass
    def w_atom(L, K):
        return F(1, 2) * F(1, 3) ** K * F(2, 3) ** L + F(1, 2) * (1 if L == 0 else 0)
    J = 40
    s = sum(w_atom(1, j) for j in range(J))
    require(w_atom(0, 0) - s - F(1, 2) == w_atom(0, J) - F(1, 2))
    require(abs(float(w_atom(0, 0) - s) - 0.5) < 1e-15)
    print("  six random discrete mixtures: positivity, split, (-Delta)^L row formula,"
          " exact fibre remainders; non-monotone row gives w(2,0)<0;"
          " atom at r=1 gives root-row defect = atom mass")


# ---------------------------------------------------------------------------
# S3. Priors: summability needs beta>1; KT prior excluded; regret table
# ---------------------------------------------------------------------------
def section_priors():
    print("== S3. Priors on the price: the root ray excludes the KT prior ==")
    # KT = beta(1/2,1/2): root-ray weight B(j+1/2,1/2)/B(1/2,1/2) = C(2j,j)/4^j
    # partial sums: sum_{j<=J} C(2j,j)/4^j = (2J+1) C(2J,J)/4^J -> infinity
    acc = F(0)
    for J in range(0, 300):
        term = F(comb(2 * J, J), 4 ** J)
        acc += term
        require(acc == F((2 * J + 1) * comb(2 * J, J), 4 ** J), J)
    require(acc > 19)  # ~ 2 sqrt(300/pi) = 19.5
    print("  KT root ray: partial sums (2J+1)C(2J,J)/4^J, at J=299: %.4f (diverges)"
          % float(acc))
    # beta(alpha,beta) mass bound 3(a+b-1)/(b-1) - b/(a+b) equals the
    # integral of (2+2r-r^2)/(1-r) = 3/(1-r) - (1-r) against beta(a,b)
    import mpmath as mp
    mp.mp.dps = 30
    for (a, b) in [(1, 2), (2, 2), (1, 3), (mp.mpf('0.5'), 2), (mp.mpf('0.5'), mp.mpf('1.5')), (3, 5)]:
        a = mp.mpf(a)
        b = mp.mpf(b)
        dens = lambda r, a=a, b=b: r ** (a - 1) * (1 - r) ** (b - 1) / mp.beta(a, b)
        val = mp.quad(lambda r: dens(r) * (3 / (1 - r) - (1 - r)), [0, 1])
        formula = 3 * (a + b - 1) / (b - 1) - b / (a + b)
        require(abs(val - formula) < mp.mpf('1e-12'), (a, b, val, formula))
    print("  beta(a,b) mass bound 3(a+b-1)/(b-1)-b/(a+b) verified against the"
          " integral for six priors; beta(1/2,2): %.4f, beta(1/2,3/2): %.4f"
          % (3 * 1.5 / 1 - 2 / 2.5, 3 * 1.0 / 0.5 - 1.5 / 2.0))
    # root-ray weight ~ j^{-beta}: summable iff beta>1 (check the exponent numerically)
    for (a, b, expo) in [(0.5, 1.5, 1.5), (1, 2, 2), (2, 2, 2), (0.5, 1.0, 1.0)]:
        j1, j2 = 2000, 4000
        w1 = mp.beta(j1 + a, b) / mp.beta(a, b)
        w2 = mp.beta(j2 + a, b) / mp.beta(a, b)
        slope = -mp.log(w2 / w1) / mp.log(2)
        require(abs(slope - expo) < 0.01, (a, b, slope))
    print("  root-ray decay exponent = beta for beta(a,b) in {(1/2,3/2),(1,2),(2,2),(1/2,1)}")
    # regret table: w_{a,b}(L,K)/M(L,K), M = best fixed price
    def M_best(L, K):
        t = L + K
        if t == 0:
            return mp.mpf(1)
        val = mp.mpf(1)
        if L:
            val *= (mp.mpf(L) / t) ** L
        if K:
            val *= (mp.mpf(K) / t) ** K
        return val
    priors = [(0.5, 0.5), (0.5, 1.5), (0.5, 2), (1, 2), (2, 2)]
    print("  regret factor w(L,K)/M(L,K) [best fixed price]; columns = priors",
          priors)
    for (L, K) in [(5, 2), (23, 7), (29, 8), (40, 15), (100, 30), (300, 100), (10, 0), (0, 10)]:
        row = []
        for (a, b) in priors:
            w = mp.beta(K + a, L + b) / mp.beta(a, b)
            row.append(float(w / M_best(L, K)))
        # the mixture never beats the best fixed price (P4 upper bound)
        require(all(x <= 1 + 1e-12 for x in row), (L, K, row))
        print("   (L,K)=(%3d,%3d): " % (L, K) + "  ".join("%.3e" % x for x in row))


# ---------------------------------------------------------------------------
# S4. Ascending shadows: Mersenne chains and the (11,7) negative cycle
# ---------------------------------------------------------------------------
def section_shadows():
    print("== S4. Ascending shadows force non-uniformity of every supersolution ==")
    for ell in range(1, 61):
        x = (1 << (ell + 1)) - 1
        for i in range(ell):
            t = 3 * x + 1
            require(v2(t) == 1, (ell, i))
            x = t >> 1
            require(x == (1 << (ell - i)) * 3 ** (i + 1) - 1, (ell, i))
        require(x == 2 * 3 ** ell - 1)
    print("  2^(l+1)-1 -> 2*3^l-1 in l U-steps of valuation 1, l<=60;"
          " ratio (2*3^l-1)/(2^(l+1)-1) -> (3/2)^l")
    # the -17 cycle: word (1,1,1,2,1,1,4), 7 odd steps, 11 halvings
    n = -17
    word = []
    for _ in range(7):
        t = 3 * n + 1
        a = v2(t)
        word.append(a)
        n = t >> a
    require(n == -17 and word == [1, 1, 1, 2, 1, 1, 4] and sum(word) == 11)
    require(3 ** 7 - 2 ** 11 == 139)
    # positive shadows n = 2^(11m) t - 17 follow the word m times and land on 3^(7m) t - 17
    for m in range(1, 4):
        for tt in range(2, 11, 2):
            n0 = (1 << (11 * m)) * tt - 17
            n = n0
            for _ in range(m):
                for a_exp in word:
                    t = 3 * n + 1
                    require(v2(t) == a_exp, (m, tt))
                    n = t >> a_exp
            require(n == 3 ** (7 * m) * tt - 17, (m, tt, n))
    print("  (-17)-cycle word (1,1,1,2,1,1,4): shadows 2^(11m)t-17 climb to 3^(7m)t-17 (t even)"
          " (m<=3, t in 2..10); multiplier 3^7/2^11 = 1.0679")


# ---------------------------------------------------------------------------
# S5. The leaf-section tower: least predecessor divisible by 3^j
# ---------------------------------------------------------------------------
def section_leaf_tower():
    print("== S5. Leaf section tower: least predecessor divisible by 3^j ==")
    table1 = {1: 6, 2: 5, 4: 4, 5: 1, 7: 2, 8: 3}
    table_unit = {1: 2, 2: 1, 4: 2, 5: 3, 7: 4, 8: 1}
    sup_ratio = {1: F(0), 2: F(0), 3: F(0)}
    arg = {}
    for v in range(1, 200001, 2):
        if v % 3 == 0:
            continue
        for j in (1, 2, 3):
            mod = 3 ** (j + 1)
            a = 1
            x = (2 * v) % mod
            while x != 1:
                x = (2 * x) % mod
                a += 1
            require(a <= 2 * 3 ** j, (v, j, a))
            z = ((1 << a) * v - 1) // 3
            require(z % 3 ** j == 0 and z % 2 == 1 and U(z) == v, (v, j))
            # least such predecessor: smaller exponents fail
            for b in range(1, a):
                if ((1 << b) * v - 1) % 3 == 0:
                    require((((1 << b) * v - 1) // 3) % 3 ** j != 0, (v, j, b))
            if j == 1:
                require(a == table1[v % 9], (v, a))
                require(z < 22 * v)
            r = F(z, v)
            if r > sup_ratio[j]:
                sup_ratio[j] = r
                arg[j] = v
        # least unit predecessor: a<=4
        a = 1
        while True:
            x = (1 << a) * v
            if x % 3 == 1 and ((x - 1) // 3) % 3 != 0:
                break
            a += 1
        require(a == table_unit[v % 9] and a <= 4, (v, a))
    for j in (1, 2, 3):
        bound = F(1 << (2 * 3 ** j), 3)
        require(sup_ratio[j] < bound, (j, sup_ratio[j], bound))
        print("  j=%d: sup_(v<2e5) rho_j(v)/v = %.6f at v=%d; bound 2^(2*3^j)/3 = %.6g"
              % (j, float(sup_ratio[j]), arg[j], float(bound)))
    # sharpness: v = 1 mod 3^(j+1) needs the full period a = 2*3^j
    for j in (1, 2, 3):
        mod = 3 ** (j + 1)
        v = 1 + mod * 2  # odd unit congruent to 1
        a = 1
        x = (2 * v) % mod
        while x != 1:
            x = (2 * x) % mod
            a += 1
        require(a == 2 * 3 ** j, (j, a))
    # the leaf predecessors of v are exactly S^(3t)(rho_1(v))
    for v in range(1, 2000, 2):
        if v % 3 == 0:
            continue
        a = table1[v % 9]
        z = ((1 << a) * v - 1) // 3
        for t in range(0, 4):
            zt = z
            for _ in range(3 * t):
                zt = 4 * zt + 1
            require(zt == ((1 << (a + 6 * t)) * v - 1) // 3 and zt % 3 == 0 and U(zt) == v)
    print("  table a_0(v mod 9) = {1:6,2:5,4:4,5:1,7:2,8:3}; unit-predecessor"
          " exponents {1:2,2:1,4:2,5:3,7:4,8:1} (bound 16/3);"
          " leaf predecessors of v are S^(3t) rho_1(v)")


# ---------------------------------------------------------------------------
# S6. Finite heads (numba) with Codex controls
# ---------------------------------------------------------------------------
@njit(cache=True)
def _counters_head(limit):
    """Return arrays L,K for odd n = 1,3,5,... < limit (index i -> n=2i+1)."""
    m = (limit + 1) // 2
    Ls = np.empty(m, dtype=np.int64)
    Ks = np.empty(m, dtype=np.int64)
    for i in range(m):
        n = 2 * i + 1
        if n == 1:
            Ls[i] = 0
            Ks[i] = 0
            continue
        L = -1
        K = 0
        while n != 1:
            t = 3 * n + 1
            a = 0
            while (t & 1) == 0:
                t >>= 1
                a += 1
            K += (a - 1) // 2
            L += 1
            n = t
        Ls[i] = L
        Ks[i] = K
    return Ls, Ks


def section_heads(head_bits):
    print("== S6. Finite heads below 2^%d (lower bounds on the masses) ==" % head_bits)
    limit = 1 << head_bits
    t0 = time.time()
    Ls, Ks = _counters_head(limit)
    # spot checks against the pure-Python counters
    for n in (1, 3, 5, 7, 27, 53, 113, 111, 155, 999, 12345, 65535, 77031):
        if n < limit:
            require((int(Ls[(n - 1) // 2]), int(Ks[(n - 1) // 2])) == root_counters(n), n)
    L = Ls.astype(np.float64)
    K = Ks.astype(np.float64)
    lg = np.vectorize(lgamma)
    logW = log(2) + lg(K + 1) + lg(L + 2) - lg(L + K + 3)
    logE = log(6) + lg(K + 2) + lg(L + 2) - lg(L + K + 4)
    Wv = np.exp(logW)
    Ev = np.exp(logE)
    n_arr = 2 * np.arange(len(Ls)) + 1
    leaf = (n_arr % 3 == 0)
    out = {}
    for bits in sorted(set([15, head_bits])):
        sel = n_arr < (1 << bits)
        Wm = float(Wv[sel].sum())
        Em = float(Ev[sel].sum())
        lamW = float(Wv[sel & leaf].sum())
        lamE = float(Ev[sel & leaf].sum())
        # fibre-complete lower bounds: root ray (2 for W, 3 for E) plus, for every
        # enumerated base b>1 (odd, not 5 mod 8), its entire sibling ray, whose
        # exact mass is w(L-1,K) by the split identity (Corollary A1)
        base = sel & (n_arr % 8 != 5) & (n_arr > 1)
        Lb = L[base]
        Kb = K[base]
        fibW = 2.0 + float(np.exp(log(2) + lg(Kb + 1) + lg(Lb + 1) - lg(Lb + Kb + 2)).sum())
        fibE = 3.0 + float(np.exp(log(6) + lg(Kb + 2) + lg(Lb + 1) - lg(Lb + Kb + 3)).sum())
        require(fibW >= Wm - 1e-9 and fibE >= Em - 1e-9, (bits, fibW, Wm, fibE, Em))
        out[bits] = dict(W=Wm, E=Em, lamW=lamW, lamE=lamE, W_fibre=fibW, E_fibre=fibE)
        print("  below 2^%d: W-head %.6f  E-head %.6f  lambda_W-head %.6f  lambda_E-head %.6f"
              % (bits, Wm, Em, lamW, lamE))
        print("             fibre-complete lower bounds: W >= %.6f  E >= %.6f" % (fibW, fibE))
    # Codex controls (exact fractions from collatz_adaptive_mixture_flow_20261005.json)
    codexW = 10824789322034673764702404645758141487821737730395371349423857 / \
        4164090077319143230070339136140045075852921948800740078822000
    codexLam = 2612987431485206015290435616122397567387944198435189 / \
        4566948047536734466280351985102680076107685597828000
    require(abs(out[15]['W'] - codexW) < 1e-9, (out[15]['W'], codexW))
    require(abs(out[15]['lamW'] - codexLam) < 1e-9, (out[15]['lamW'], codexLam))
    print("  controls: Codex finite head below 2^15: W %.9f lambda %.9f reproduced"
          % (codexW, codexLam))
    # layer masses by tau = L+1 (odd steps) for W, to see the geometric-type decay
    tau = Ls + 1
    layers = []
    for t in range(1, 12):
        sel = (tau == t)
        layers.append(float(Wv[sel].sum()))
    print("  W-head mass by odd-step layer tau=1..11:", " ".join("%.4f" % x for x in layers))
    print("  (head time %.1fs; max counters L=%d K=%d)" % (time.time() - t0, Ls.max(), Ks.max()))
    return out, layers


# ---------------------------------------------------------------------------
# S7. Residue-level upper bound on fixed-price base mass, integrated
# ---------------------------------------------------------------------------
@njit(cache=True)
def _lift_table(P, N):
    """child[c',t0] (mod N) and allowed[c',t0] for lifts c' mod P=3N, t0<P."""
    child = np.empty((P, P), dtype=np.int32)
    allowed = np.zeros((P, P), dtype=np.uint8)
    for c in range(P):
        y = c % P
        for t in range(P):
            if y % 3 != 0:
                k0 = 1 if (y % 3 == 2) else 2
                num = ((1 << k0) * y - 1) % P
                # num is divisible by 3 as a residue mod P (P divisible by 3)
                child[c, t] = (num // 3) % N
                allowed[c, t] = 1
            else:
                child[c, t] = 0
            y = (4 * y + 1) % P
    return child, allowed


@njit(cache=True)
def _accumulate(child, allowed, Fv, P, N):
    """ent[c',b] = sum_{t0 allowed, child=b} F[t0] over lifts c' mod P=3N;
    v1[b] = root's first generation (c'=1, t>=1, no self-loop)."""
    ent = np.zeros((P, N), dtype=np.float64)
    for c in range(P):
        for t in range(P):
            if allowed[c, t]:
                ent[c, child[c, t]] += Fv[t]
    v1 = np.zeros(N, dtype=np.float64)
    for t in range(1, P):
        if allowed[1, t]:
            v1[child[1, t]] += Fv[t]
    if allowed[1, 0]:
        v1[child[1, 0]] += Fv[P]  # t = P, 2P, ... (not the self-loop)
    return ent, v1


def _f_t(r, t, P):
    """f_t(r) = r^t (1-r)/(1-r^P) for 0<=r<=1, t in [0,P]."""
    if r >= 1.0:
        return 1.0 / P
    if r <= 0.0:
        return 1.0 if t == 0 else 0.0
    return exp(t * log(r)) * (1.0 - r) / (1.0 - exp(P * log(r)))


def _modes(P):
    """r*_t with m(r*)=t for t < (P-1)/2 ; m(r) = truncated geometric mean."""
    ts = np.arange(0, P + 1, dtype=np.float64)
    lo = np.zeros(P + 1)
    hi = np.ones(P + 1)
    half = (P - 1) / 2.0

    def m_of(r):
        r = np.clip(r, 1e-300, 1 - 1e-15)
        rP = np.exp(P * np.log(r))
        rP1 = np.exp((P - 1) * np.log(r))
        return r * (1 - P * rP1 + (P - 1) * rP) / ((1 - r) * (1 - rP))

    for _ in range(60):
        mid = 0.5 * (lo + hi)
        val = m_of(mid)
        up = val < ts
        lo = np.where(up, mid, lo)
        hi = np.where(up, hi, mid)
    modes = 0.5 * (lo + hi)
    modes[ts >= half] = 1.0  # increasing on [0,1]
    return modes


def _sup_vector(u, v, P, modes):
    Fv = np.empty(P + 1)
    for t in range(P + 1):
        rs = modes[t]
        cand = [_f_t(u, t, P), _f_t(v, t, P)]
        if u < rs < v:
            cand.append(_f_t(rs, t, P))
        Fv[t] = max(cand) * (1 + 1e-9)
    return Fv


def _bellman_value(ent, N, max_iter=80):
    """Adversarial-lift subtree bound: the Bellman fixed point
         V(a) = 1 + max_{lift c' of a} sum_b ent[c',b] V(b)
    found by policy iteration. Any finite fixed point V>=1 dominates every
    depth-truncated adversarial value, hence the actual subtree mass of any
    base in class a. Returns None if no transient policy is reached."""
    idx = np.arange(N)
    eye = np.eye(N)
    policy = np.zeros(N, dtype=np.int64)
    ones = np.ones(N)
    V = None
    for _ in range(max_iter):
        Mpi = ent[idx + policy * N, :]
        try:
            V = np.linalg.solve(eye - Mpi, ones)
        except np.linalg.LinAlgError:
            return None
        if not np.all(np.isfinite(V)) or np.any(V < 1 - 1e-9):
            return None
        res = np.max(np.abs((eye - Mpi) @ V - ones))
        if res > 1e-7 * (1 + np.max(V)):
            return None
        Q = (ent @ V).reshape(3, N)  # Q[i,a] = row (a + i N) of ent, dotted with V
        bell = 1 + Q.max(axis=0)
        if np.max(np.abs(bell - V)) <= 1e-9 * (1 + np.max(V)):
            return V
        new_policy = np.argmax(Q, axis=0)
        if np.all(new_policy == policy):
            return V
        policy = new_policy
    return None


def residue_level_bound(level, h=1e-3, h_tail=1e-4, tail_from=0.98,
                        point_rs=(1 / 16, 1 / 4, 1 / 2, 3 / 4, 0.9, 0.99)):
    N = 3 ** level
    P = 3 * N
    child, allowed = _lift_table(P, N)
    modes = _modes(P)
    # pointwise bounds at specific r (exact weights at that r, no interval sup)
    pointwise = {}
    for r in point_rs:
        Fv = np.array([_f_t(r, t, P) for t in range(P + 1)])
        ent, v1 = _accumulate(child, allowed, Fv, P, N)
        V = _bellman_value(ent, N)
        if V is None:
            pointwise[r] = (float('inf'), 2 + 2 * r - r * r)
        else:
            pointwise[r] = (1 + float(v1 @ V), 2 + 2 * r - r * r)
    # integrated bounds: per interval min(Bellman bound, Codex P1 bound)
    edges = [0.0]
    r = 0.0
    while r < tail_from - 1e-12:
        r = min(r + h, tail_from)
        edges.append(r)
    while r < 1.0 - 1e-12:
        r = min(r + h_tail, 1.0)
        edges.append(r)
    E_bound = W_bound = 0.0
    E_codex = W_codex = 0.0
    n_mdp = 0
    worst = (0.0, 0.0)
    for i in range(len(edges) - 1):
        u, v = edges[i], edges[i + 1]
        Bc = 2 + 2 * v - u * u  # sup of 2+2r-r^2 on [u,v]  (Codex P1, PROVED)
        B = Bc
        if u > 0:
            Fv = _sup_vector(u, v, P, modes)
            ent, v1 = _accumulate(child, allowed, Fv, P, N)
            V = _bellman_value(ent, N)
            if V is not None:
                Bm = (1 + float(v1 @ V)) * (1 + 1e-7)
                if Bm < B:
                    B = Bm
                    n_mdp += 1
        if B > worst[1]:
            worst = (u, B)
        E_bound += B * 3 * (v * v - u * u)       # 6 int_u^v r dr * B
        W_bound += B * 2 * (v - u)               # 2 int_u^v dr * B
        E_codex += Bc * 3 * (v * v - u * u)
        W_codex += Bc * 2 * (v - u)
    return dict(level=level, N=N, intervals=len(edges) - 1, intervals_bellman=n_mdp,
                E_bound=E_bound, W_bound=W_bound,
                E_codex_quadrature=E_codex, W_codex_quadrature=W_codex,
                worst_interval=worst, pointwise=pointwise)


def section_residue_bounds(level_max):
    print("== S7. Adversarial-lift (Bellman) bound on fixed-price base mass ==")
    print("  (base mass bound B(r); f_r mass = B(r)/(1-r); Codex P1: 2+2r-r^2;"
          " per interval the smaller valid bound is used)")
    results = []
    for level in range(1, level_max + 1):
        t0 = time.time()
        res = residue_level_bound(level)
        results.append(res)
        pw = res['pointwise']
        print("  level 3^%d (%d classes, %d intervals, Bellman used on %d, %.0fs):"
              " E-mass <= %.4f  W-mass <= %.4f  [P1 quadrature: E %.4f W %.4f]"
              % (level, res['N'], res['intervals'], res['intervals_bellman'],
                 time.time() - t0, res['E_bound'], res['W_bound'],
                 res['E_codex_quadrature'], res['W_codex_quadrature']))
        print("    B(r):", "  ".join("r=%.4g: %.4f (P1 %.4f)" % (r, pw[r][0], pw[r][1])
                                     for r in sorted(pw)))
        # the P1 quadrature must reproduce 17/2 and 16/3 up to interval looseness
        require(abs(res['E_codex_quadrature'] - 8.5) < 0.02, res['E_codex_quadrature'])
        require(abs(res['W_codex_quadrature'] - 16 / 3) < 0.02, res['W_codex_quadrature'])
        require(res['E_bound'] <= res['E_codex_quadrature'] + 1e-12, level)
        require(res['W_bound'] <= res['W_codex_quadrature'] + 1e-12, level)
    return results



# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--level-max', type=int, default=6)
    ap.add_argument('--head-bits', type=int, default=20)
    ap.add_argument('--json', type=str, default='')
    args = ap.parse_args()
    t0 = time.time()
    section_decode()
    section_classification()
    section_priors()
    section_shadows()
    section_leaf_tower()
    heads, layers = section_heads(args.head_bits)
    res = section_residue_bounds(args.level_max)
    best = res[-1]
    print("== Summary ==")
    print("  lower bounds (bases below 2^%d, full fibres): W >= %.6f, E >= %.6f"
          % (args.head_bits, heads[args.head_bits]['W_fibre'], heads[args.head_bits]['E_fibre']))
    print("  upper bounds (level 3^%d): W <= %.4f, E <= %.4f  (proved Codex bounds 16/3 and 11; P1 gives 17/2)"
          % (best['level'], best['W_bound'], best['E_bound']))
    print("  the leaf masses are exact: lambda_W = 1, lambda_E = 2 (Codex P5, B2)")
    print("  checks: %d, total time %.1fs" % (CHECKS, time.time() - t0))
    if args.json:
        payload = dict(
            status="PROVED scoped identities; FINITE-EXACT; VERIFIED numerical residue bounds;"
                   " universal positivity OPEN",
            checks=CHECKS, head_bits=args.head_bits, heads=heads, layers_W=layers,
            residue=[{k: (v if k != 'pointwise' else {str(kk): vv for kk, vv in v.items()})
                      for k, v in r.items()} for r in res])
        with open(args.json, 'w') as fh:
            json.dump(payload, fh, indent=1, default=str)
        print("  json written:", args.json)


if __name__ == '__main__':
    main()
