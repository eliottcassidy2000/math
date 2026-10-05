#!/usr/bin/env python3
"""Refinement floors, tower conjugation, tail cost of the atom criterion, and the
shadow theorem: backward-descent cones are 3-adic shadows of negative rational
cycle points whose 2-adic shadows (mod 2^(A+1)) are the forward-rising cylinders.

Companion note:
  05-knowledge/results/collatz_refinement_floor_shadows_20261005.md

Conventions: U(n)=oddpart(3n+1); ROOT word (a_1..a_tau) = valuations along the
odd orbit; counters L=tau-1, K=sum floor((a_i-1)/2); W(n)=2K!(L+1)!/(L+K+2)!;
fixed price f_r(n)=(1-r)^L r^K on rooted n. A word w=(a_1..a_l) is RISING iff
2^A < 3^l with A=sum a_i; its carry is c_w=sum_i 3^(l-i) 2^(a_1+..+a_(i-1)) and
its cycle point is x_w=c_w/(2^A-3^l) (negative rational, 2- and 3-adic integer).
No universal convergence assumption; finite checks certify the listed universe.

Usage: python3 <this file> [--sieve-bits 24] [--cone-depth 12] [--samples 20000]
"""
from fractions import Fraction as F
from math import factorial, lgamma, exp, log
import argparse
import json
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
    return (n & -n).bit_length() - 1


def U(n):
    t = 3 * n + 1
    return t >> v2(t)


def root_word(n, cap=10 ** 6):
    w = []
    while n != 1:
        t = 3 * n + 1
        a = v2(t)
        w.append(a)
        n = t >> a
        if len(w) > cap:
            raise RuntimeError("cap")
    return w


def counters(word):
    if not word:
        return 0, 0
    return len(word) - 1, sum((a - 1) // 2 for a in word)


def wW(L, K):
    return F(2 * factorial(K) * factorial(L + 1), factorial(L + K + 2))


def f_price(r, L, K):
    return (1 - r) ** L * r ** K


# ---------------------------------------------------------------------------
# S1. Tower conjugation along an orbit prefix (exact)
# ---------------------------------------------------------------------------
def section_tower():
    print("== S1. Tower conjugation: 2-adic cylinder of n <-> 3-adic class of U^t(n) ==")
    r = F(1, 5)
    tested = 0
    for n in (7, 27, 97, 231, 703, 871, 6171):
        w = root_word(n)
        for t in range(1, min(len(w), 9)):
            A = sum(w[:t])
            Kpre = sum((a - 1) // 2 for a in w[:t])
            img = n
            for _ in range(t):
                img = U(img)
            for s in range(0, 12):
                npr = n + (1 << (A + 1)) * s
                wp = root_word(npr)
                require(wp[:t] == w[:t], ("prefix", n, t, s))
                img2 = npr
                for _ in range(t):
                    img2 = U(img2)
                require(img2 == img + 2 * 3 ** t * s, ("image", n, t, s))
                Lp, Kp = counters(wp)
                Li, Ki = counters(root_word(img2))
                if img2 != 1:
                    require(Lp == t + Li and Kp == Kpre + Ki, ("counters", n, t, s))
                    require(f_price(r, Lp, Kp) == (1 - r) ** t * r ** Kpre * f_price(r, Li, Ki))
                tested += 1
    print("  exact on %d (n,t,s) triples: prefix fixed, U^t(n+2^(A+1)s)=U^t(n)+2*3^t s,"
          " f_r multiplicative with the prefix price" % tested)
    # floors are stopping-time bounds: W(n) <= 2/(L+2) and W(n) <= w(L,0)
    for L in range(0, 60):
        for K in range(0, 60):
            require(wW(L, K) <= F(2, L + 2))
    print("  W(n) <= 2/(L+2): an epsilon-floor forces tau(n) <= 2/epsilon - 1")


# ---------------------------------------------------------------------------
# S2. Tail cost of the atom criterion p_m >= L - R(q)
# ---------------------------------------------------------------------------
def section_tail_cost():
    print("== S2. Tail cost of the refinement criterion (Codex (10)) ==")
    # two-step family z_k=(2 S^k(1)-1)/3, k = 1 mod 9: lambda(z_k)=4/((k+1)(k+2)(k+3))
    for k in (1, 10, 19, 28):
        nk = (4 ** (k + 1) - 1) // 3  # S^k(1)
        require(nk % 9 == 5)
        zk = (2 * nk - 1) // 3
        require(zk % 3 == 0 and U(zk) == nk and U(nk) == 1)
        L, K = counters(root_word(zk))
        require((L, K) == (1, k), (k, L, K))
        require(wW(L, K) == F(4, (k + 1) * (k + 2) * (k + 3)))
    W27 = wW(*counters(root_word(27)))
    print("  lambda(27) = W(27) = %s = %.3e (27 is a leaf; 41 odd steps)" % (W27, float(W27)))
    # R(q) >= sum over the certified two-step family of z_k >= q, k=1 mod 9
    def tail_lower(q):
        k = 1
        tot = F(0)
        while True:
            nk = (4 ** (k + 1) - 1) // 3
            zk = (2 * nk - 1) // 3
            if zk >= q:
                tot += F(4, (k + 1) * (k + 2) * (k + 3))
                if k > 4000:
                    break
            k += 9
        return tot
    for e in (10, 20, 50, 100, 200):
        q = 10 ** e
        print("  q=10^%-3d: R(q) >= %.3e from the two-step family alone" % (e, float(tail_lower(q))))
    # the modulus needed for lambda(27): sum_{k>=k0} 4/k^3 ~ 2/k0^2 < W27 -> k0 > sqrt(2/W27)
    k0 = (2 / float(W27)) ** 0.5
    print("  criterion (10) at the atom 27 needs k0 ~ %.2e sibling levels, i.e. q >~ 4^(k0) = 10^%.0f"
          % (k0, k0 * log(4) / log(10)))
    # fixed price r: root ray tail beyond size q is sum_{j: S^j(1)>=q} r^j ~ q^{-log_4(1/r)}
    for r in (0.5, 0.25, 0.0625):
        L27, K27 = counters(root_word(27))
        atom = f_price(r, L27, K27)
        expo = log(1 / r) / log(4)
        q_needed = atom ** (-1 / expo)
        print("  fixed price r=%.4f: f_r(27)=%.3e, root-ray tail ~ q^(-%.3f), needs q >~ %.2e"
              % (r, atom, expo, q_needed))


# ---------------------------------------------------------------------------
# S3. Shadow theorem: rising words, cycle points, 2-adic and 3-adic shadows
# ---------------------------------------------------------------------------
def rising_words(max_len):
    """All valuation words of length l<=max_len with 2^A < 3^l."""
    out = []

    def rec(prefix, A, l):
        if l > 0 and (1 << A) < 3 ** l:
            out.append(tuple(prefix))
        if l == max_len:
            return
        # remaining letters are at least 1 each; prune when impossible to rise
        for a in range(1, 64):
            if (1 << (A + a + (max_len - l - 1))) >= 3 ** max_len and (1 << (A + a)) >= 3 ** (l + 1):
                break
            rec(prefix + [a], A + a, l + 1)

    rec([], 0, 0)
    return out


def carry(word):
    l = len(word)
    c = 0
    Apre = 0
    for i, a in enumerate(word):
        c += 3 ** (l - 1 - i) * (1 << Apre)
        Apre += a
    return c, Apre


def section_shadows(max_len):
    print("== S3. Shadow theorem (rising words of length <= %d) ==" % max_len)
    words = rising_words(max_len)
    print("  rising words: %d" % len(words))
    checked = 0
    classes = {}
    for w in words:
        l = len(w)
        c, A = carry(w)
        d = (1 << A) - 3 ** l  # negative
        require(d < 0)
        x = F(c, d)
        require(x < 0)
        r2 = (c * pow(d % (1 << (A + 1)), -1, 1 << (A + 1))) % (1 << (A + 1))  # word exact needs one more bit
        r3 = (c * pow(d % 3 ** l, -1, 3 ** l)) % 3 ** l
        require(r2 % 2 == 1)
        classes.setdefault(l, set()).add(r3)
        for s in range(0, 4):
            m = r2 + (1 << (A + 1)) * s
            v = m
            for a in w:
                t = 3 * v + 1
                require(v2(t) == a, ("forward", w, s))
                v = t >> a
            n = v
            require(n == (3 ** l * m + c) >> A, ("formula", w, s))
            require(n % 3 ** l == r3 and n % 2 == 1, ("3-adic", w, s))
            require(n > m, ("rising", w, s))
            # backward chain along w from n returns to m
            u = n
            for a in reversed(w):
                num = (1 << a) * u - 1
                require(num % 3 == 0, ("integral", w, s))
                u = num // 3
                require(u > 0)
            require(u == m, ("backward", w, s))
            checked += 1
    print("  %d exact forward/backward shadow checks; m=x_w mod 2^(A+1) -> n=x_w mod 3^l bijectively" % checked)
    # the three integer negative cycles
    for w, pt in (((1,), -1), ((1, 2), -5), ((2, 1), -7), ((1, 1, 1, 2, 1, 1, 4), -17)):
        c, A = carry(w)
        require(F(c, (1 << A) - 3 ** len(w)) == pt, (w, c))
    cw112, _ = carry((1, 1, 2))
    require(F(cw112, 16 - 27) == F(-19, 11))
    print("  cycle points: (1)->-1, (1,2)->-5, (2,1)->-7, (1,1,1,2,1,1,4)->-17, (1,1,2)->-19/11")
    # cumulative density of the descent set from cones of depth <= l
    dens = {}
    for l in range(1, max_len + 1):
        marked = set()
        for j in range(1, l + 1):
            for r in classes.get(j, ()):
                # lift class r mod 3^j to classes mod 3^l
                step = 3 ** j
                for t in range(3 ** (l - j)):
                    marked.add(r + step * t)
        dens[l] = F(len(marked), 3 ** l)
    print("  density of the backward-descent set from cones of depth <= l:")
    print("   ", "  ".join("l=%d: %.5f" % (l, float(dens[l])) for l in sorted(dens)))
    # new (primitive) cone counts per depth
    prim = {}
    for l in range(1, max_len + 1):
        base = set()
        for j in range(1, l):
            for r in classes.get(j, ()):
                step = 3 ** j
                for t in range(3 ** (l - j)):
                    base.add(r + step * t)
        prim[l] = len(classes.get(l, set()) - base)
    print("  new cone classes per depth:", prim)
    return dens, classes


# ---------------------------------------------------------------------------
# S4. Census of the descent set D (visited by a smaller odd start)
# ---------------------------------------------------------------------------
@njit(cache=True)
def _sieve(limit):
    m_odd = (limit + 1) // 2
    mark = np.zeros(m_odd, dtype=np.uint8)
    maxv = 0
    for i in range(m_odd):
        m = 2 * i + 1
        v = m
        while True:
            t = 3 * v + 1
            while (t & 1) == 0:
                t >>= 1
            v = t
            if v < m or v == 1:
                break
            if v > maxv:
                maxv = v
            if v < limit:
                mark[(v - 1) // 2] = 1
    return mark, maxv


def section_census(sieve_bits, dens_cone):
    print("== S4. Census of D = {odd n : some odd m<n has n on its orbit} ==")
    limit = 1 << sieve_bits
    t0 = time.time()
    mark, maxv = _sieve(limit)
    n = 2 * np.arange(len(mark)) + 1
    unit = (n % 3 != 0)
    leaf = ~unit
    require(not mark[leaf].any())  # leaves have no predecessors at all
    dD = mark.mean()
    dU = mark[unit].mean()
    d1 = mark[n % 3 == 1].mean()
    d2 = mark[n % 3 == 2].mean()
    print("  below 2^%d (%.1fs, max excursion %.3e): dens(D) = %.5f; among units %.5f;"
          " n=1 mod 3: %.5f; n=2 mod 3: %.5f" % (sieve_bits, time.time() - t0, maxv, dD, dU, d1, d2))
    require(d2 > 0.999999, d2)  # every n = 2 mod 3 has the smaller predecessor (2n-1)/3
    minima = [int(x) for x in n[(~mark.astype(bool)) & unit][:25]]
    print("  first unit basin minima (no smaller ancestor):", minima)
    require(minima[:4] == [1, 7, 19, 31] or minima[:3] == [1, 7, 19], minima[:5])
    # compare with cone densities: the sieve density must dominate every finite-depth cone density
    for l, dl in dens_cone.items():
        require(dD >= float(dl) - 1e-6, (l, dl, dD))
    # density of D by dyadic blocks (stability)
    blocks = []
    for b in range(10, sieve_bits):
        sel = (n >= (1 << b)) & (n < (1 << (b + 1)))
        blocks.append((b, float(mark[sel].mean())))
    print("  dens(D) by dyadic block:", " ".join("[2^%d,2^%d): %.4f" % (b, b + 1, d) for b, d in blocks[-6:]))
    return dict(density=float(dD), units=float(dU), minima=minima, blocks=blocks)


# ---------------------------------------------------------------------------
# S5. Forward versus backward descent depth on random sources
# ---------------------------------------------------------------------------
def section_forward_backward(samples, classes, max_len):
    print("== S5. Forward stopping depth versus backward descent depth (random odd n < 2^40) ==")
    rng = np.random.default_rng(2026_10_05)
    fwd = []
    bwd = []
    joint_unpaid = {}
    for _ in range(samples):
        n = int(rng.integers(1, 1 << 39)) * 2 + 1
        # forward: first odd step with U^k(n) < n
        v = n
        k = 0
        while True:
            v = U(v)
            k += 1
            if v < n:
                break
            if k > 10000:
                k = 10 ** 9
                break
        fwd.append(k)
        # backward: least l <= max_len with n mod 3^l in a rising-word class
        b = 10 ** 9
        for l in range(1, max_len + 1):
            if (n % 3 ** l) in classes.get(l, ()):
                b = l
                break
        bwd.append(b)
    fwd = np.array(fwd)
    bwd = np.array(bwd)
    print("  forward odd-step stopping depth: mean %.3f, P(>5)=%.4f, P(>10)=%.4f, P(>20)=%.4f"
          % (fwd.mean(), (fwd > 5).mean(), (fwd > 10).mean(), (fwd > 20).mean()))
    print("  backward descent depth (<= %d): P(none)=%.4f, P(=1)=%.4f, P(=2)=%.4f, P(>=3)=%.4f"
          % (max_len, (bwd > max_len).mean(), (bwd == 1).mean(), (bwd == 2).mean(),
             ((bwd >= 3) & (bwd <= max_len)).mean()))
    for kf in (3, 5, 8, 12):
        unp_f = (fwd > kf).mean()
        unp_both = ((fwd > kf) & (bwd > max_len)).mean()
        joint_unpaid[kf] = (float(unp_f), float(unp_both))
        print("  depth %2d: forward-unpaid %.4f, forward-and-backward-unpaid %.4f" % (kf, unp_f, unp_both))
    return joint_unpaid


# ---------------------------------------------------------------------------
# S6. Exact defect ledger of segment closures (general form of Codex C4/C5)
# ---------------------------------------------------------------------------
def nu(m):
    """Codex source prior nu(m)=8/(3*4^bitlength(m)) on positive odd m."""
    return F(8, 3 * 4 ** m.bit_length())


def section_segment_ledger(limit_starts, limit_targets):
    print("== S6. Segment closures: defect(y) = nu(y) - nu{starts exiting at y} ==")
    rules = {
        'stopping (exit at first value below the start)': lambda m, v: v < m,
        'single rise (exit at the first valuation >= 2)': None,  # handled inline
    }
    from collections import defaultdict
    for name in rules:
        incoming = defaultdict(F)   # sum over starts m of nu(m) * #{p in seg(m): U(p)=y}
        through = defaultdict(F)    # sum over starts m != y with y in seg(m)
        exits = defaultdict(F)      # sum over starts m with exit(m) = y
        for m in range(3, limit_starts, 2):
            seg = [m]
            v = m
            while True:
                t = 3 * v + 1
                a = v2(t)
                u = t >> a
                if name.startswith('stopping'):
                    stop = (u < m) or u == 1
                else:
                    stop = (a >= 2) or u == 1
                if stop:
                    exit_point = u
                    break
                seg.append(u)
                v = u
            w = nu(m)
            for p in seg:
                y = U(p)
                if y < limit_targets:
                    incoming[y] += w
            for y in seg[1:]:
                if y < limit_targets:
                    through[y] += w
            if exit_point < limit_targets:
                exits[exit_point] += w
        viol = []
        for y in range(3, limit_targets, 2):
            require(incoming[y] == through[y] + exits[y], (name, y))
            if exits[y] > nu(y):
                viol.append((y, exits[y] / nu(y)))
        print("  rule: %s" % name)
        print("    identity K W = (through) + (exit mass) checked at %d targets; violations (exit mass > nu(y)): %d"
              % ((limit_targets - 3) // 2, len(viol)))
        print("    first violations (y, exit/nu):", [(y, "%.3f" % float(r)) for y, r in viol[:8]])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sieve-bits', type=int, default=24)
    ap.add_argument('--cone-depth', type=int, default=12)
    ap.add_argument('--samples', type=int, default=20000)
    ap.add_argument('--json', type=str, default='')
    args = ap.parse_args()
    t0 = time.time()
    section_tower()
    section_tail_cost()
    dens, classes = section_shadows(args.cone_depth)
    census = section_census(args.sieve_bits, dens)
    joint = section_forward_backward(args.samples, classes, args.cone_depth)
    section_segment_ledger(1 << 16, 1 << 12)
    print("== Summary ==")
    print("  checks: %d, total time %.1fs" % (CHECKS, time.time() - t0))
    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(dict(checks=CHECKS, cone_density={l: str(d) for l, d in dens.items()},
                           census=census, joint_unpaid=joint,
                           status="PROVED scoped identities; FINITE-EXACT; VERIFIED census"),
                      fh, indent=1, default=str)
        print("  json written:", args.json)


if __name__ == '__main__':
    main()
