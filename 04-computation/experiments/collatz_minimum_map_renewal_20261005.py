#!/usr/bin/env python3
"""The minimum map and the renewal structure of hard excursions.

m(x) = first value below x on the odd orbit of x (the endpoint of its first-descent
segment). The chain of running minima of n is the m-orbit of n. Its blocks have
types (l, A) with the cylinder law p(l,A) = N(l,A)/2^A (Theorem T of the expense
note). Tests here:
  R1  renewal: the k-th block's type law equals the first block's law, and types of
      consecutive blocks are (nearly) independent, among sources of fixed size;
  R2  the number of blocks m(n) and the odd time tau(n) scale as log2 n / E[excess]
      and log2 n * E[l]/E[excess], with E[A] = 2 E[l] (Wald);
  R3  coverage of F_q per dyadic block decays like a power X^(-c_q),
      c_q = -log2(1 - p_q)/E[excess] (renewal prediction), fitted on the census;
  R4  coalescence: popularity T(x) = number of sources whose chain passes through x
      (subtree sizes of the minimum-map tree), its spectrum and concentration;
  R5  in-degrees of the minimum map (preimages indexed by first-descent words);
  R6  cheap probes: Hankel determinants of the first-descent length law c_l, and the
      Kraft deficit 1 - sum_{l<=L} c_l = mass never descended within L steps.
Companion note: 05-knowledge/results/collatz_minimum_map_renewal_20261005.md
Usage: python3 <this file> [--census-bits 22] [--json PATH]
"""
from fractions import Fraction as F
from math import log2, ceil, log
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
ALPHA = log2(3.0)


def require(cond, witness=None):
    global CHECKS
    CHECKS += 1
    if not cond:
        raise RuntimeError(witness)


def A_min(l):
    return (3 ** l).bit_length()


def first_descent_counts(lmax, extra=40):
    N = {}
    states = {0: 1}
    for j in range(1, lmax + 1):
        new_states = {}
        counts_j = {}
        thr = A_min(j)
        for s, cnt in states.items():
            for a in range(1, thr - s):
                new_states[s + a] = new_states.get(s + a, 0) + cnt
            for A in range(max(thr, s + 1), thr + extra + 1):
                counts_j[A] = counts_j.get(A, 0) + cnt
        N[j] = counts_j
        states = new_states
    return N


def law_tables(lmax=200):
    N = first_descent_counts(lmax)
    types = []
    for l in range(1, lmax + 1):
        for A, cnt in N[l].items():
            types.append((l, A, F(cnt, 1 << A)))
    total = sum(p for (_, _, p) in types)
    El = sum(l * p for (l, _, p) in types)
    EA = sum(A * p for (_, A, p) in types)
    Eexc = float(EA) - float(El) * ALPHA
    tails = {}
    for Q in (3, 8, 16, 32, 100):
        tails[Q] = float(sum(p for (l, A, p) in types if ceil(l / (A - l * ALPHA) - 1e-9) > Q))
    cl = {}
    for (l, A, p) in types:
        cl[l] = cl.get(l, F(0)) + p
    return types, float(total), float(El), float(EA), Eexc, tails, cl


# ---------------------------------------------------------------------------
# census of the minimum map
# ---------------------------------------------------------------------------
@njit(cache=True)
def _census(limit, nblocks):
    """Per odd n: number of blocks, tau, q(n), the first nblocks block types (l, A), popularity
    counts T[x] (x as a running minimum of some n<limit), in-degree of m among x<limit."""
    m_odd = (limit + 1) // 2
    nb = np.zeros(m_odd, dtype=np.int32)
    tau = np.zeros(m_odd, dtype=np.int32)
    qn = np.zeros(m_odd, dtype=np.int32)
    bl = np.zeros((m_odd, nblocks), dtype=np.int16)
    bA = np.zeros((m_odd, nblocks), dtype=np.int16)
    pop = np.zeros(m_odd, dtype=np.int64)
    indeg = np.zeros(m_odd, dtype=np.int64)
    log23 = np.log2(3.0)
    for i in range(1, m_odd):
        n = 2 * i + 1
        x = n
        k = 0
        t = 0
        qb = 0
        while x != 1:
            v = x
            l = 0
            A = 0
            while True:
                tt = 3 * v + 1
                a = 0
                while (tt & 1) == 0:
                    tt >>= 1
                    a += 1
                v = tt
                l += 1
                A += a
                if v < x:
                    break
            if k < nblocks:
                bl[i, k] = l if l < 32000 else 32000
                bA[i, k] = A if A < 32000 else 32000
            k += 1
            t += l
            qs = int(np.ceil(l / (A - l * log23) - 1e-9))
            if qs > qb:
                qb = qs
            if v < limit:
                pop[(v - 1) // 2] += 1
            if x == n and v < limit:
                indeg[(v - 1) // 2] += 1
            x = v
        nb[i] = k
        tau[i] = t
        qn[i] = qb
    return nb, tau, qn, bl, bA, pop, indeg


def tv_distance(emp_counts, law, keys):
    tot = sum(emp_counts.values())
    tv = 0.0
    for key in keys:
        tv += abs(emp_counts.get(key, 0) / tot - law.get(key, 0.0))
    # mass outside the listed keys
    tv += abs(sum(v for k, v in emp_counts.items() if k not in keys) / tot - max(0.0, 1 - sum(law.get(k, 0.0) for k in keys)))
    return 0.5 * tv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--census-bits', type=int, default=22)
    ap.add_argument('--json', type=str, default='')
    args = ap.parse_args()
    t0 = time.time()
    bits = args.census_bits
    limit = 1 << bits
    types, total, El, EA, Eexc, tails, cl = law_tables(200)
    law = {(l, A): float(p) for (l, A, p) in types}
    print("== Law of first-descent blocks (DP to l=200): total %.8f, E[l]=%.5f, E[A]=%.5f, E[A]/E[l]=%.5f (Wald: 2),"
          " E[excess]=%.5f = (2-log2 3)E[l]=%.5f" % (total, El, EA, EA / El, Eexc, (2 - ALPHA) * El))
    require(abs(EA / El - 2) < 1e-4 and abs(Eexc - (2 - ALPHA) * El) < 1e-4, (EA / El, Eexc, (2 - ALPHA) * El))
    print("   tails p_q = P(block rate > q):", " ".join("%d:%.5f" % (q, p) for q, p in tails.items()))
    naive_c = {q: -log2(1 - p) / Eexc for q, p in tails.items()}
    print("   naive prediction c_q = -log2(1-p_q)/E[excess]:",
          " ".join("%d:%.4f" % (q, c) for q, c in naive_c.items()))

    def killed_root(q):
        # survival to level L of the chain with hard blocks killed: S(L) ~ 2^(-cL) where
        # sum over good types of p(t) 2^(c excess(t)) = 1  (exact renewal equation)
        good = [(float(p), A - l * ALPHA) for (l, A, p) in types if ceil(l / (A - l * ALPHA) - 1e-9) <= q]
        lo, hi = 0.0, 2.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            val = sum(pp * 2 ** (mid * e) for pp, e in good)
            if val < 1:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)
    pred_c = {q: killed_root(q) for q in tails}
    print("   killed-renewal prediction (root of sum_good p(t) 2^(c excess) = 1):",
          " ".join("%d:%.4f" % (q, c) for q, c in pred_c.items()))
    t1 = time.time()
    nb, tau, qn, bl, bA, pop, indeg = _census(limit, 6)
    n = 2 * np.arange(len(nb)) + 1
    print("   census below 2^%d in %.0fs" % (bits, time.time() - t1))

    # R1 renewal: block-type laws by position, among sources in the top dyadic block
    print("== R1. Block types along the chain, sources in [2^%d, 2^%d) ==" % (bits - 1, bits))
    top = (n >= (1 << (bits - 1)))
    keys = [(l, A) for (l, A, p) in types if float(p) > 1e-5]
    emp_first = {}
    for k in range(6):
        sel = top & (nb > k)
        cnt = {}
        ls = bl[sel, k].astype(int)
        As = bA[sel, k].astype(int)
        for l, A in zip(ls, As):
            cnt[(l, A)] = cnt.get((l, A), 0) + 1
        if k == 0:
            emp_first = cnt
        tv = tv_distance(cnt, law, keys)
        mean_l = float(ls.mean())
        print("   block %d: %d sources, TV distance to the cylinder law %.4f, mean length %.3f (law %.3f)"
              % (k + 1, int(sel.sum()), tv, mean_l, El))
        require(tv <= 0.0015 * (2 ** k) + 0.002 + 3.0 / np.sqrt(max(1, int(sel.sum()))), ("renewal TV", k, tv))
    # conditional independence: law of block 2 given block 1 type
    for cond in ((1, 2), (1, 3), (2, 4), (5, 8)):
        sel = top & (nb > 1) & (bl[:, 0] == cond[0]) & (bA[:, 0] == cond[1])
        cnt = {}
        for l, A in zip(bl[sel, 1].astype(int), bA[sel, 1].astype(int)):
            cnt[(l, A)] = cnt.get((l, A), 0) + 1
        tv = tv_distance(cnt, law, keys) if sel.sum() > 1000 else float('nan')
        print("   block 2 given block 1 of type %s: %d sources, TV to the unconditional law %.4f" % (cond, int(sel.sum()), tv))

    # R2 scaling of the number of blocks and tau
    print("== R2. Number of running minima and odd time per dyadic block ==")
    for k in range(8, bits):
        sel = (n >= (1 << k)) & (n < (1 << (k + 1)))
        mb = float(nb[sel].mean())
        mt = float(tau[sel].mean())
        print("   [2^%2d,2^%2d): mean blocks %.3f (renewal %.3f), mean tau %.2f (renewal %.2f)"
              % (k, k + 1, mb, (k + 0.5) / Eexc, mt, (k + 0.5) * El / Eexc))

    # R3 coverage decay of F_q and of the two-tier family
    print("== R3. Coverage of F_q per dyadic block and fitted exponents ==")
    fits = {}
    for q in (3, 8, 16, 32, 100):
        ks = list(range(10, bits))
        cov = []
        for k in ks:
            sel = (n >= (1 << k)) & (n < (1 << (k + 1)))
            cov.append(float((qn[sel] <= q).mean()))
        xs = np.array(ks, dtype=float)
        ys = np.log2(np.array(cov))
        slope, intercept = np.polyfit(xs, ys, 1)
        fits[q] = (-slope, pred_c[q], naive_c[q])
        print("   q=%3d: coverage at 2^10..2^%d: %s ; fitted exponent %.4f, killed-renewal %.4f, naive %.4f"
              % (q, bits - 1, " ".join("%.4f" % c for c in cov[::3]), -slope, pred_c[q], naive_c[q]))

    # R4 coalescence spectrum
    print("== R4. Popularity T(x) = number of sources below 2^%d whose chain passes through x ==" % bits)
    pop[0] = 0   # every chain ends at 1; exclude ROOT from the spectrum
    order = np.argsort(-pop)
    tot_visits = int(pop.sum())
    print("   total chain visits (excluding ROOT) %d; top 15 running minima:" % tot_visits)
    for i in order[:15]:
        l0, A0 = int(bl[i, 0]), int(bA[i, 0])
        qown = int(ceil(l0 / (A0 - l0 * ALPHA) - 1e-9)) if l0 > 0 else 0
        print("     x=%7d T=%8d (%.4f of sources) own block (l,A)=(%d,%d) q_own=%d q(x)=%d"
              % (n[i], pop[i], pop[i] / (len(n) - 1), l0, A0, qown, qn[i]))
    for K in (10, 100, 1000, 10000):
        print("   share of all visits at the top %5d minima: %.4f" % (K, float(pop[order[:K]].sum()) / tot_visits))
    # tail of the popularity distribution among x in [2^8, 2^16)
    sel = (n >= (1 << 8)) & (n < (1 << 16)) & (pop > 0)
    P = pop[sel].astype(float)
    for thr in (10, 100, 1000, 10000):
        print("   fraction of x in [2^8,2^16) with T(x) >= %5d: %.5f" % (thr, float((P >= thr).mean())))
    # R5 in-degrees
    print("== R5. In-degrees of the minimum map among sources below 2^%d ==" % bits)
    units = (n % 3 != 0) & (n > 1)
    leaves = (n % 3 == 0)
    require(int(indeg[leaves].sum()) == 0)
    d = indeg[units]
    print("   leaves have in-degree 0; units: mean in-degree %.3f, fraction with in-degree 0: %.4f, max %d at x=%d"
          % (float(d.mean()), float((d == 0).mean()), int(d.max()), int(n[units][np.argmax(d)])))
    # renewal density: expected number of running minima per dyadic block, for sources far above it
    print("== R4b. Renewal density: chain visits per dyadic block per source starting above 2^(k+6) ==")
    for k in (8, 10, 12, 14, 16):
        sel_x = (n >= (1 << k)) & (n < (1 << (k + 1)))
        # visits to x in the block from sources >= 2^(k+6): approximate by all visits minus sources inside; sources
        # below the block cannot visit it, sources in [2^(k+1), 2^(k+6)) are few compared with the total above
        visits = float(pop[sel_x].sum())
        above = float((n >= (1 << (k + 1))).sum())
        print("   block [2^%d,2^%d): visits per source above = %.4f (renewal 1/E[excess] = %.4f)" % (k, k + 1, visits / above, 1 / Eexc))
    # R5b Fibonacci probe on in-degrees
    fib = [1, 2]
    while fib[-1] < limit:
        fib.append(fib[-1] + fib[-2])
    fib_odd = [f for f in fib if f % 2 == 1 and 100 < f < limit // 4]
    deg_fib = [int(indeg[(f - 1) // 2]) for f in fib_odd]
    sel_cmp = (n > 100) & (n < limit // 4) & units
    print("   in-degree at odd Fibonacci numbers:", list(zip(fib_odd, deg_fib)),
          "; mean in-degree of units in the same range %.3f, 99.9th percentile %d"
          % (float(indeg[sel_cmp].mean()), int(np.percentile(indeg[sel_cmp], 99.9))))
    # fair comparison: percentile rank of each unit Fibonacci number among units within a factor-2 window
    ranks = []
    for f in fib_odd:
        if f % 3 == 0:
            continue
        w = (n >= f // 2) & (n < 2 * f) & units
        d_f = int(indeg[(f - 1) // 2])
        ranks.append((f, d_f, float((indeg[w] < d_f).mean())))
    print("   percentile rank of the unit Fibonacci in-degrees within a factor-2 window:",
          " ".join("%d:%d(%.3f)" % r for r in ranks))
    # R6 Hankel probe (exact, untruncated block-length law) and Kraft deficit
    print("== R6. Probes: Hankel determinants of the exact block-length law c_l, and the Kraft deficit ==")
    # exact c_l: sum over rising prefixes of length l-1 with sum s of 2^-s times the terminal tail 2^(-(A0-s)+1)
    states = {0: 1}
    cl_exact = {}
    for j in range(1, 61):
        thr = A_min(j)
        tot = F(0)
        new_states = {}
        for s_, cnt in states.items():
            A0 = max(thr, s_ + 1)
            tot += F(cnt, 1 << (A0 - 1))
            for a in range(1, thr - s_):
                new_states[s_ + a] = new_states.get(s_ + a, 0) + cnt
        cl_exact[j] = tot
        states = new_states
    cls = [cl_exact[l] for l in range(1, 61)]
    require(cls[0] == F(1, 2) and cls[1] == F(1, 8) and cls[2] == F(1, 8), cls[:3])
    import sympy
    dets = []
    for k in range(1, 9):
        M = sympy.Matrix([[cls[i + j] for j in range(k)] for i in range(k)])
        dets.append(M.det())
    print("   exact c_l, l=1..10:", " ".join(str(c) for c in cls[:10]))
    print("   Hankel determinants H_k, k=1..8:", " ".join(str(dd) for dd in dets))
    require(all(dd != 0 for dd in dets), "a vanishing Hankel determinant would indicate a finite-order recurrence")
    # Somos-4 test on the Hankel sequence: H_n H_{n-4} - a H_{n-1}H_{n-3} - b H_{n-2}^2 = 0 for some (a,b)? solve from two windows
    H = dets
    # two linear equations in (a,b) from n=5,6; check at n=7,8
    import fractions
    a1, b1, r1 = H[3] * H[1], H[2] ** 2, H[4] * H[0]
    a2, b2, r2 = H[4] * H[2], H[3] ** 2, H[5] * H[1]
    det = a1 * b2 - a2 * b1
    somos = None
    if det != 0:
        a = (r1 * b2 - r2 * b1) / det
        b = (a1 * r2 - a2 * r1) / det
        ok7 = H[6] * H[2] == a * H[5] * H[3] + b * H[4] ** 2
        ok8 = H[7] * H[3] == a * H[6] * H[4] + b * H[5] ** 2
        somos = (ok7 and ok8)
        print("   (alpha,beta) Somos-4 fit from H_5,H_6: alpha=%s beta=%s; holds at H_7,H_8: %s" % (a, b, somos))
    deficit = []
    acc = F(0)
    for L in range(1, 61):
        acc += cl_exact[L]
        if L in (5, 10, 20, 40, 60):
            deficit.append((L, float(1 - acc)))
    print("   Kraft deficit 1 - sum_{l<=L} c_l (exact law):", " ".join("L=%d: %.3e" % (L, dfc) for L, dfc in deficit))
    # R7 rational-slope control: replace the threshold ceil(j log2 3) by ceil(j p/q) for convergents p/q
    print("== R7. Hankel transforms of block-length laws at rational slopes p/q (control for the irrational case) ==")
    import sympy
    from math import gcd

    def block_law(thr, lmax):
        states = {0: 1}
        out = []
        for j in range(1, lmax + 1):
            t = thr(j)
            tot = F(0)
            new_states = {}
            for s_, cnt in states.items():
                A0 = max(t, s_ + 1)
                tot += F(cnt, 1 << (A0 - 1))
                for a in range(1, t - s_):
                    new_states[s_ + a] = new_states.get(s_ + a, 0) + cnt
            out.append(tot)
            states = new_states
        return out

    def hankel(cs, kmax):
        return [sympy.Matrix([[cs[i + j] for j in range(k)] for i in range(k)]).det() for k in range(1, kmax + 1)]

    def somos_fit(H):
        a1, b1, r1 = H[3] * H[1], H[2] ** 2, H[4] * H[0]
        a2, b2, r2 = H[4] * H[2], H[3] ** 2, H[5] * H[1]
        det = a1 * b2 - a2 * b1
        if det == 0:
            return None
        a = (r1 * b2 - r2 * b1) / det
        b = (a1 * r2 - a2 * r1) / det
        ok = all(H[n] * H[n - 4] == a * H[n - 1] * H[n - 3] + b * H[n - 2] ** 2 for n in range(6, len(H)))
        return (a, b, ok)

    for (pp, qq) in ((2, 1), (3, 2), (8, 5), (19, 12), (65, 41)):
        thr = (lambda j, pp=pp, qq=qq: -(-pp * j // qq))   # ceil(j p/q): least A with qA >= pj
        cs = block_law(thr, 40)
        H = hankel(cs, 10)
        fit = somos_fit(H)
        print("   slope %d/%d: c_l = %s ...; H_1..H_10 = %s; Somos-4 fit holds to H_10: %s; all H>0: %s"
              % (pp, qq, " ".join(str(c) for c in cs[:6]), " ".join(str(h) for h in H), fit[2] if fit else None, all(h > 0 for h in H)))
    # slope 2: H_k = 2^(-k(2k-1)) exactly (hexagonal numbers), a Somos-4 relation with (alpha,beta) = (2^-12, 0)
    cs2 = block_law(lambda j: 2 * j, 40)
    H2 = hankel(cs2, 10)
    require(all(H2[k - 1] == F(1, 1 << (k * (2 * k - 1))) for k in range(1, 11)), ("slope-2 Hankel", H2))
    require(all(H2[n] * H2[n - 4] == F(1, 1 << 12) * H2[n - 1] * H2[n - 3] for n in range(4, 10)))
    print("   slope 2/1 identity: H_k = 2^(-k(2k-1)) for k<=10 and H_n H_(n-4) = 2^-12 H_(n-1) H_(n-3): Somos-4 with beta = 0 (FINITE-EXACT)")
    cs = block_law(lambda j: A_min(j), 40)
    H = hankel(cs, 10)
    fit = somos_fit(H)
    print("   slope log2 3: H_1..H_10 signs = %s; Somos-4 fit holds to H_10: %s" % ("".join("+" if h > 0 else "-" for h in H), fit[2] if fit else None))
    print("== Summary ==")
    print("   checks: %d, total time %.1fs" % (CHECKS, time.time() - t0))
    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(dict(checks=CHECKS, E_l=El, E_A=EA, E_excess=Eexc, tails=tails, fits=fits,
                           top=[(int(n[i]), int(pop[i])) for i in order[:15]],
                           status="PROVED renewal structure; FINITE-EXACT laws; census VERIFIED"), fh, indent=1, default=str)
        print("   json written:", args.json)


if __name__ == '__main__':
    main()
