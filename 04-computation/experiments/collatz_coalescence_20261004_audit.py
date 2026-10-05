#!/usr/bin/env python
# -*- coding: ascii -*-
"""
collatz_coalescence_20261004_audit.py -- independent auditor's blind re-derivation.

Written from scratch WITHOUT reading collatz_coalescence_20261004_probes.py, to check
the numerical claims of

  05-knowledge/results/collatz_coalescence_20261004_proof_strategies_and_bold_predictions.md
  05-knowledge/hypotheses/HYP-9174-x2x3-lonely-spectrum-is-the-artin-coordinate-of-2-3.md
  05-knowledge/hypotheses/HYP-9175-powers-of-three-shadow-only-cycle-points-1-or-3-mod-8-...md

PART 1  I(t) = inf_{j,k>=0} || 2^j 3^k t ||  on reduced fractions m/q.
        I(m/q) = g_q(m)/q where g_q(m) = min of D(x) = min(x, q-x) over the forward orbit
        of m under x -> 2x, x -> 3x (mod q).  Algorithm: residues are processed in
        increasing D and the value is propagated backwards along reverse edges (exact, O(q)).
        A second, per-m forward-orbit BFS and the coset formula are used as cross-checks.
PART 2  Collatz T(n) = n/2 (even), (3n+1)/2 (odd).  Bad_k = odd residues mod 2^k whose
        parity word has 3^{o_j} > 2^j for all 1 <= j <= k.  Exact DP on words + direct
        numpy enumeration of residues; 2-adic closure of <3>; Lucas/Fibonacci identity.

Usage:  python collatz_coalescence_20261004_audit.py > collatz_coalescence_20261004_audit.out
Everything is exact (Python ints / Fractions) except the floating-point DP of 2(d),
which is cross-checked against exact big-integer counts.
"""
import sys
import time
from fractions import Fraction
from math import gcd
from collections import Counter

import numpy as np

T0 = time.time()
FAILS = []


def check(name, mine, claimed):
    ok = (mine == claimed)
    print("  CHECK %-70s %s" % (name, "PASS" if ok else "FAIL"))
    if not ok:
        print("      mine    = %r" % (mine,))
        print("      claimed = %r" % (claimed,))
        FAILS.append(name)
    return ok


def elapsed():
    return "[t=%.1fs]" % (time.time() - T0)


# ----------------------------------------------------------------------------------
# PART 1 -- the x2,x3 lonely function
# ----------------------------------------------------------------------------------

def lonely_table(q):
    """g[r] = min over the forward orbit of r (under x->2x, x->3x mod q) of min(x, q-x).
    Exact.  Residue 0 has g = 0, so denominators of the form 2^a 3^b give I = 0."""
    pred = [[] for _ in range(q)]
    for r in range(q):
        pred[(2 * r) % q].append(r)
        pred[(3 * r) % q].append(r)
    g = [-1] * q
    order = [0]
    for d in range(1, q // 2 + 1):
        order.append(d)
        if q - d != d:
            order.append(q - d)
    for v in order:                      # increasing D(v)
        if g[v] >= 0:
            continue
        d = min(v, q - v)
        g[v] = d
        stack = [v]
        while stack:                     # every node that can reach v and is still
            x = stack.pop()              # unassigned gets the value d
            for u in pred[x]:
                if g[u] < 0:
                    g[u] = d
                    stack.append(u)
    return g


def lonely_direct(m, q):
    """Independent check: explicit forward orbit of m under x2, x3 mod q."""
    m %= q
    seen = {m}
    stack = [m]
    best = min(m, q - m)
    while stack:
        x = stack.pop()
        for y in ((2 * x) % q, (3 * x) % q):
            if y not in seen:
                seen.add(y)
                stack.append(y)
                best = min(best, y, q - y)
    return best


def subgroup_23(q):
    """The subgroup <2,3> of (Z/q)^x, q prime to 6 (as a set)."""
    H = {1}
    stack = [1]
    while stack:
        x = stack.pop()
        for y in ((2 * x) % q, (3 * x) % q):
            if y not in H:
                H.add(y)
                stack.append(y)
    return H


def phi(q):
    return sum(1 for m in range(1, q + 1) if gcd(m, q) == 1)


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def prime_factors(n):
    f = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            f.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        f.append(n)
    return f


def mult_order(a, q, factors_of_order):
    """Order of a in the cyclic group (Z/q)^x, q prime; factors_of_order = primes of q-1."""
    o = q - 1
    for p in factors_of_order:
        while o % p == 0 and pow(a, o // p, q) == 1:
            o //= p
    return o


def part1():
    print("=" * 100)
    print("PART 1: I(t) = inf_{j,k} ||2^j 3^k t||")
    print("=" * 100)
    ONE14 = Fraction(1, 14)

    # (a) -----------------------------------------------------------------------
    print("\n(a) reduced m/q in (0,1), q <= 600, with I(m/q) >= 1/14")
    pts = []
    for q in range(2, 601):
        g = lonely_table(q)
        for m in range(1, q):
            if gcd(m, q) == 1:
                v = Fraction(g[m], q)
                if v >= ONE14:
                    pts.append((q, m, v))
    vals = Counter(v for _, _, v in pts)
    dens = sorted(set(q for q, _, _ in pts))
    byden = Counter(q for q, _, _ in pts)
    print("  number of points           :", len(pts))
    print("  denominators               :", dens)
    print("  points per denominator     :", sorted(byden.items()))
    print("  values (value: multiplicity):", sorted(((str(v), c) for v, c in vals.items()),
                                                    key=lambda t: -Fraction(t[0])))
    print("  (denominator, value) table :")
    dv = Counter((q, str(v)) for q, _, v in pts)
    for (q, v), c in sorted(dv.items()):
        print("      q=%-4d I=%-6s  x%d" % (q, v, c))
    check("1a count = 88", len(pts), 88)
    check("1a denominators = {5,7,10,11,13,14,26,28,33,52,56}", dens,
          [5, 7, 10, 11, 13, 14, 26, 28, 33, 52, 56])
    check("1a values 1/5x4 1/7x6 1/10x4 1/11x20 1/13x30 1/14x24",
          {str(v): c for v, c in vals.items()},
          {"1/5": 4, "1/7": 6, "1/10": 4, "1/11": 20, "1/13": 30, "1/14": 24})
    # extra: point list for the record
    print("  explicit list:", ", ".join("%d/%d" % (m, q) for q, m, _ in sorted(pts)))
    print(" ", elapsed())

    # (b) -----------------------------------------------------------------------
    print("\n(b) q <= 300 prime to 6: coset formula vs direct orbit (per m) vs table algorithm")
    nq = nm = 0
    bad = []
    for q in range(5, 301):
        if q % 2 == 0 or q % 3 == 0:
            continue
        H = subgroup_23(q)
        g = lonely_table(q)
        nq += 1
        for m in range(1, q):
            if gcd(m, q) != 1:
                continue
            coset_min = min(min(r, q - r) for r in ((m * h) % q for h in H))
            direct = lonely_direct(m, q)
            nm += 1
            if not (coset_min == direct == g[m]):
                bad.append((q, m, coset_min, direct, g[m]))
    print("  denominators tested:", nq, " fractions tested:", nm, " mismatches:", len(bad))
    if bad:
        print("  first mismatches:", bad[:10])
    check("1b coset formula == direct orbit == table algorithm (q<=300 prime to 6)", bad, [])
    print(" ", elapsed())

    # (c) -----------------------------------------------------------------------
    print("\n(c) largest 12 values of I strictly below 1/14, q <= 4000")
    best_p6 = {}
    best_all = {}
    for q in range(2, 4001):
        g = lonely_table(q)
        p6 = (q % 2 != 0) and (q % 3 != 0)
        ds = set(g[m] for m in range(1, q) if gcd(m, q) == 1)
        for d in ds:
            v = Fraction(d, q)
            if v < ONE14:
                if v not in best_all:
                    best_all[v] = q
                if p6 and v not in best_p6:
                    best_p6[v] = q
    top_p6_20 = sorted(best_p6.items(), key=lambda t: -t[0])[:20]
    top_p6 = top_p6_20[:12]
    top_all = sorted(best_all.items(), key=lambda t: -t[0])[:14]
    print("  q prime to 6 (the claim's range), top 20:")
    for v, q in top_p6_20:
        idx = phi(q) // len(subgroup_23(q))
        print("      I = %-8s = %.6f   smallest q = %-5d  index [(Z/q)^x:<2,3>] = %d"
              % (v, float(v), q, idx))
    note20 = [Fraction(5, 73), Fraction(1, 17), Fraction(1, 19), Fraction(5, 97),
              Fraction(13, 259), Fraction(7, 145), Fraction(29, 601), Fraction(11, 247),
              Fraction(31, 697), Fraction(19, 431), Fraction(1, 23), Fraction(1, 25),
              Fraction(1, 29), Fraction(13, 385), Fraction(1, 31), Fraction(7, 241),
              Fraction(1, 35), Fraction(1, 37), Fraction(13, 485), Fraction(23, 865)]
    check("1c note section 6 top-20 list (q<=4000 prime to 6)", [v for v, _ in top_p6_20], note20)
    note_idx = {73: 2, 97: 2, 259: 3, 145: 2, 601: 8, 247: 3, 697: 8, 431: 10, 23: 2,
                385: 2, 241: 2, 485: 8, 865: 4}
    mine_idx = {q: phi(q) // len(subgroup_23(q)) for _, q in top_p6_20 if q in note_idx}
    check("1c indices quoted in the note for those moduli", mine_idx, note_idx)
    claimed = [Fraction(5, 73), Fraction(1, 17), Fraction(1, 19), Fraction(5, 97),
               Fraction(13, 259), Fraction(7, 145), Fraction(29, 601), Fraction(11, 247),
               Fraction(31, 697), Fraction(19, 431), Fraction(1, 23), Fraction(1, 25)]
    check("1c top-12 values below 1/14 over q<=4000 prime to 6",
          [v for v, _ in top_p6], claimed)
    check("1c smallest q attaining each value is the printed denominator",
          [q for _, q in top_p6], [v.denominator for v in claimed])
    print("  ALL q <= 4000 (including q divisible by 2 or 3) -- robustness check:")
    for v, q in top_all:
        print("      I = %-8s = %.6f   smallest q = %d" % (v, float(v), q))
    extra = [(v, q) for v, q in top_all if v not in best_p6]
    print("  values in the all-q top list that do NOT occur for q prime to 6:", extra)
    print(" ", elapsed())

    # (d) -----------------------------------------------------------------------
    print("\n(d) primes q <= 20000, q = 1 (mod 24), index exactly 2: q*max_m I(m/q) vs least QNR")
    P = [p for p in primes_upto(20000) if p >= 5]
    full = 0
    idx2_24 = []
    index_hist = Counter()
    for q in P:
        fo = prime_factors(q - 1)
        o2 = mult_order(2, q, fo)
        o3 = mult_order(3, q, fo)
        hsize = o2 * o3 // gcd(o2, o3)          # cyclic group: <2,3> has order lcm
        idx = (q - 1) // hsize
        index_hist[idx] += 1
        if idx == 1:
            full += 1
        if q % 24 == 1 and idx == 2:
            idx2_24.append(q)
    # cross-check the lcm formula for |<2,3>| against the explicit subgroup for q <= 3000
    mism = 0
    for q in P:
        if q > 3000:
            break
        fo = prime_factors(q - 1)
        o2 = mult_order(2, q, fo)
        o3 = mult_order(3, q, fo)
        if len(subgroup_23(q)) != o2 * o3 // gcd(o2, o3):
            mism += 1
    check("1d |<2,3>| = lcm(ord 2, ord 3) agrees with explicit subgroup (primes <= 3000)", mism, 0)
    print("  primes 5 <= q <= 20000:", len(P))
    print("  index histogram:", sorted(index_hist.items()))
    print("  primes with <2,3> = (Z/q)^x:", full, " fraction = %.4f" % (full / len(P)))
    check("1d number of primes 5<=q<=20000", len(P), 2260)
    check("1d primes with <2,3> = (Z/q)^x = 1598", full, 1598)
    artin2 = 1.0
    for l in primes_upto(10 ** 6):
        artin2 *= 1 - 1 / (l * l * (l - 1))
    se = (artin2 * (1 - artin2) / len(P)) ** 0.5
    print("  naive two-generator Artin product prod_l (1 - 1/(l^2 (l-1))) = %.5f ; measured %.4f ; "
          "difference %.4f = %.2f standard errors (SE %.4f)"
          % (artin2, full / len(P), full / len(P) - artin2, (full / len(P) - artin2) / se, se))
    check("1d two-generator Artin product rounds to 0.6975", round(artin2, 4), 0.6975)
    print("  primes q = 1 mod 24 with index exactly 2:", len(idx2_24))
    check("1d count of q = 1 mod 24 with index 2 = 186", len(idx2_24), 186)
    exceptions = []
    lnr_hist = Counter()
    for q in idx2_24:
        g = lonely_table(q)
        qmax = max(g[m] for m in range(1, q))     # all m in 1..q-1 are units (q prime)
        lnr = next(r for r in range(2, q) if pow(r, (q - 1) // 2, q) == q - 1)
        lnr_hist[lnr] += 1
        if qmax != lnr:
            exceptions.append((q, qmax, lnr))
    print("  least-QNR histogram over those primes:", sorted(lnr_hist.items()))
    print("  exceptions (q, q*max I, least QNR):", exceptions)
    check("1d q*max_m I(m/q) == least quadratic non-residue, no exception", exceptions, [])
    print("  first few such primes:", idx2_24[:12])
    # extra: the largest I_max over primes <= 20000 with a proper subgroup <2,3>
    print("  largest I_max(q) = max_m I(m/q) over primes q <= 20000 with index > 1:")
    imax = []
    for q in P:
        fo = prime_factors(q - 1)
        o2 = mult_order(2, q, fo)
        o3 = mult_order(3, q, fo)
        idx = (q - 1) // (o2 * o3 // gcd(o2, o3))
        if idx == 1:
            continue
        g = lonely_table(q)
        imax.append((Fraction(max(g[1:]), q), q, idx))
    imax.sort(key=lambda t: -t[0])
    for v, q, idx in imax[:12]:
        print("      I_max = %-10s = %.6f  q = %-6d index %d" % (v, float(v), q, idx))
    note_imax = [(Fraction(5, 73), 2), (Fraction(5, 97), 2), (Fraction(29, 601), 8),
                 (Fraction(19, 431), 10), (Fraction(1, 23), 2), (Fraction(197, 6563), 17),
                 (Fraction(7, 241), 2), (Fraction(5, 193), 2), (Fraction(11, 439), 2),
                 (Fraction(109, 4513), 12), (Fraction(29, 1201), 2), (Fraction(157, 6553), 56)]
    check("1d note's top-12 I_max over primes with proper <2,3> (values)",
          [v for v, _, _ in imax[:12]], [v for v, _ in note_imax])
    print("      indices quoted in the note (where given): 6563->17, 4513->12, 6553->56; mine:",
          [(q, idx) for _, q, idx in imax[:12]])
    print(" ", elapsed())


# ----------------------------------------------------------------------------------
# PART 2 -- Collatz no-descent residues, 2-adic closure of <3>, Lucas identity
# ----------------------------------------------------------------------------------

def thresholds(K):
    """thr[j] = least o with 3^o > 2^j (exact integers)."""
    thr = [0] * (K + 1)
    o, p3 = 0, 1
    for j in range(K + 1):
        p2 = 1 << j
        while p3 <= p2:
            o += 1
            p3 *= 3
        thr[j] = o
    return thr


def bad_counts_exact(K):
    """Exact word counts.  Words of length j surviving 3^{o_i} > 2^i for all i <= j.
    Only the prefixes 110 (residues 3 mod 8) and 111 (residues 7 mod 8) survive j <= 3."""
    thr = thresholds(K)
    a = [0] * (K + 2)
    b = [0] * (K + 2)
    a[2] = 1          # prefix 110: three steps, two odd
    b[3] = 1          # prefix 111: three steps, three odd
    N110 = {3: 1}
    N111 = {3: 1}
    for j in range(4, K + 1):
        t = thr[j]
        a = [0] * t + [a[o] + a[o - 1] for o in range(t, K + 2)]
        b = [0] * t + [b[o] + b[o - 1] for o in range(t, K + 2)]
        N110[j] = sum(a)
        N111[j] = sum(b)
    return N110, N111


def bad_residues_direct(k, chunk=1 << 20):
    """Direct enumeration of odd residues r mod 2^k: iterate T on the representative,
    record the parity word, test 3^{o_j} > 2^j for all j <= k.  Returns the sorted residues."""
    thr = thresholds(k)
    K = 1 << k
    out = []
    for start in range(1, K, 2 * chunk):
        r = np.arange(start, min(K, start + 2 * chunk), 2, dtype=np.int64)
        x = r.copy()
        o = np.zeros(len(r), dtype=np.int32)
        alive = np.ones(len(r), dtype=bool)
        for j in range(1, k + 1):
            odd = (x & 1) == 1
            o += odd
            alive &= (o >= thr[j])
            x = np.where(odd, (3 * x + 1) // 2, x // 2)
        out.append(r[alive])
    return np.concatenate(out)


def share_dp_float(K, checkpoints):
    """Floating-point DP for the share s_k = N110(k) / (N110(k)+N111(k)) with joint
    renormalisation (each letter has probability 1/2; only the ratio matters)."""
    thr = thresholds(K)
    a = np.zeros(K + 2)
    b = np.zeros(K + 2)
    a[2] = 1.0
    b[3] = 1.0
    res = {}
    for j in range(4, K + 1):
        a = 0.5 * (a + np.concatenate(([0.0], a[:-1])))
        b = 0.5 * (b + np.concatenate(([0.0], b[:-1])))
        t = thr[j]
        a[:t] = 0.0
        b[:t] = 0.0
        s = a.sum() + b.sum()
        a /= s
        b /= s
        if j in checkpoints:
            res[j] = a.sum() / (a.sum() + b.sum())
    return res


def part2():
    print("\n" + "=" * 100)
    print("PART 2: Collatz no-coefficient-descent residues")
    print("=" * 100)
    KMAX = 24
    N110, N111 = bad_counts_exact(KMAX)
    tot = [N110[k] + N111[k] for k in range(4, KMAX + 1)]
    n3 = [N110[k] for k in range(4, KMAX + 1)]
    print("\n(a) |Bad_k| (exact word DP), k = 4..24:")
    print("   ", tot)
    claimed_tot = [3, 4, 8, 13, 19, 38, 64, 128, 226, 367, 734, 1295, 2114, 4228, 7495,
                   14990, 27328, 46611, 93222, 168807, 286581]
    check("2a |Bad_k| k=4..24", tot, claimed_tot)
    print("\n(b) |Bad_k cap {1,3 mod 8}| = |Bad_k cap {3 mod 8}| (prefix 110), k = 4..24:")
    print("   ", n3)
    claimed_n3 = [1, 1, 2, 3, 4, 8, 13, 26, 45, 71, 142, 247, 395, 790, 1386, 2772, 5019,
                  8468, 16936, 30485, 51280]
    check("2b |Bad_k cap {1,3 mod 8}| k=4..24", n3, claimed_n3)

    # direct residue enumeration, k = 4..24 (independent of the word DP)
    print("\n  direct residue enumeration (numpy), k = 4..24:")
    direct_tot, direct_13, direct_mod8 = [], [], set()
    for k in range(4, KMAX + 1):
        R = bad_residues_direct(k)
        direct_tot.append(len(R))
        m8 = R % 8
        direct_13.append(int(np.sum((m8 == 1) | (m8 == 3))))
        direct_mod8 |= set(np.unique(m8).tolist())
        print("    k=%2d |Bad_k|=%7d  #(1 or 3 mod 8)=%6d  residues mod 8 present: %s"
              % (k, len(R), direct_13[-1], sorted(np.unique(m8).tolist())))
    check("2a direct residue enumeration == word DP, k=4..24", direct_tot, tot)
    check("2b direct #(1,3 mod 8) == word DP prefix-110 count", direct_13, n3)
    check("2b every element of Bad_k is 3 or 7 mod 8 (k=4..24)", sorted(direct_mod8), [3, 7])
    print(" ", elapsed())

    # (c) subgroup <3> of (Z/2^k)^x
    print("\n(c) <3> in (Z/2^k)^x equals {u odd : u = 1 or 3 mod 8}, and -1 not in <3>, k=3..20")
    okc = True
    for k in range(3, 21):
        M = 1 << k
        S = set()
        x = 1
        for _ in range(M // 4):
            S.add(x)
            x = (3 * x) % M
        target = set(u for u in range(1, M, 2) if u % 8 in (1, 3))
        ok = (S == target) and ((M - 1) not in S) and (len(S) == M // 4)
        okc &= ok
        if not ok:
            print("    k=%d FAIL" % k)
    check("2c <3> = {1,3 mod 8} and -1 not in <3>, k=3..20", okc, True)
    # Probe F's other half: 2^(L/2) = -1 mod 3^n, L = 2*3^(n-1), n = 2..20
    okF = all(pow(2, 3 ** (n - 1), 3 ** n) == 3 ** n - 1 for n in range(2, 21))
    check("2c (Probe F) 2^(L/2) = -1 mod 3^n for n=2..20", okF, True)

    # (d) share s_k for large k: EXACT big-integer DP to k = 6000 (sparse dict), float DP beyond
    print("\n(d) share s_k = |Bad_k cap {3 mod 8}| / |Bad_k| for large k")
    KX = 6000
    thr = thresholds(KX)
    a = {2: 1}
    b = {3: 1}
    cps_exact = [24, 50, 100, 200, 400, 800, 1600, 3200, 6000]
    ex = {}
    Wk = {}
    for j in range(4, KX + 1):
        t = thr[j]
        na, nb = {}, {}
        for o, c in a.items():
            if o >= t:
                na[o] = na.get(o, 0) + c
            if o + 1 >= t:
                na[o + 1] = na.get(o + 1, 0) + c
        for o, c in b.items():
            if o >= t:
                nb[o] = nb.get(o, 0) + c
            if o + 1 >= t:
                nb[o + 1] = nb.get(o + 1, 0) + c
        a, b = na, nb
        if j in cps_exact or j in (799, 1599, 5999):
            n1, n2 = sum(a.values()), sum(b.values())
            Wk[j] = n1 + n2
            if j in cps_exact:
                ex[j] = (n1, n1 + n2)
    cps_float = [6000, 12000, 24000]
    fl = share_dp_float(24000, set(cps_float))
    import math
    r = math.log(2) / math.log(3)
    hH = -r * math.log2(r) - (1 - r) * math.log2(1 - r)
    print("    h = H_2(log_3 2) = %.5f  (THM-4479/4495's exponent; 2^h = 2 min_theta M(theta) = %.5f)"
          % (hH, 2 ** hH))
    note_shares = {50: 0.16983, 100: 0.16504, 200: 0.16254, 400: 0.16119, 800: 0.16047,
                   1600: 0.16010, 3200: 0.15991, 6000: 0.15982}
    sk = {}
    for k in cps_exact:
        n1, tot_k = ex[k]
        sk[k] = n1 / tot_k
        Ck = math.exp(math.log(tot_k) + 1.5 * math.log(k) - hH * k * math.log(2))
        line = "    k=%5d  exact s_k = %.6f" % (k, sk[k])
        if k in note_shares:
            line += "  (note: %.5f)" % note_shares[k]
        line += "   C_k = W_k k^(3/2) 2^(-hk) = %.4f" % Ck
        print(line)
    check("2d exact shares round to the note's five-digit values at k=50..6000",
          {k: round(sk[k], 5) for k in note_shares}, note_shares)
    check("2d float DP s_6000 == exact s_6000 to 1e-12", abs(fl[6000] - sk[6000]) < 1e-12, True)
    for k in (12000, 24000):
        print("    k=%5d  float DP s_k = %.6f" % (k, fl[k]))
    extrap_1 = (6000 * sk[6000] - 1600 * sk[1600]) / (6000 - 1600)
    extrap_2 = (24000 * fl[24000] - 6000 * fl[6000]) / (24000 - 6000)
    print("    1/k extrapolation from (1600, 6000):  %.6f" % extrap_1)
    print("    1/k extrapolation from (6000, 24000): %.6f" % extrap_2)
    print("    s_6000 = %.5f (claimed 0.15982); limit claimed 0.1597 to four digits" % sk[6000])
    check("2d s_6000 rounds to 0.15982", round(sk[6000], 5), 0.15982)
    check("2d limit (1/k extrapolation) rounds to 0.1597", round(extrap_2, 4), 0.1597)
    # growth of W_k = |Bad_k| for the record (Cramer rate 2^h, polynomial factor k^(-3/2))
    print("    W_6000/W_5999 = %.6f ; (W_6000/W_3200)^(1/2800) = %.6f ; 2^h = %.6f"
          % (Wk[6000] / Wk[5999], math.exp((math.log(Wk[6000]) - math.log(Wk[3200])) / 2800), 2 ** hH))
    # HYP-9175 (iii): "resisting exponent classes ... a fraction ~ 0.16 x 2 x 2^(-(1-h)k) k^(-3/2)"
    print("    HYP-9175(iii): fraction of exponent classes a mod 2^(k-2) that resist = |Bad_k cap <3>| / 2^(k-2):")
    for k in (24, 100, 400, 1600, 6000):
        n1, tot_k = ex[k]
        frac_log2 = math.log2(n1) - (k - 2)
        claimed_log2 = math.log2(0.16 * 2) - (1 - hH) * k - 1.5 * math.log2(k)
        print("      k=%5d  log2(fraction) = %.3f ; log2 of the HYP's 0.16*2*2^(-(1-h)k) k^(-3/2) = %.3f ; ratio = %.2f"
              % (k, frac_log2, claimed_log2, 2 ** (frac_log2 - claimed_log2)))
    print(" ", elapsed())

    # (e) 2-adic digits of log_3(-5), closure of <3>, cycle points
    print("\n(e) a mod 2^(k-2) with 3^a = -5 (mod 2^k), k = 3..14")
    alist = []
    for k in range(3, 15):
        M = 1 << k
        sols = [a for a in range(M // 4) if pow(3, a, M) == (-5) % M]
        alist.append(sols[0] if len(sols) == 1 else sols)
    print("   ", alist)
    check("2e a_k = 1,3,3,11,11,11,11,11,267,267,1291,3339 (unique each k)",
          alist, [1, 3, 3, 11, 11, 11, 11, 11, 267, 267, 1291, 3339])
    print("    consistency: a_k = a_{k+1} mod 2^(k-2):",
          all(alist[i + 1] % (1 << (i + 1)) == alist[i] for i in range(len(alist) - 1)))
    print("    2-adic expansion of a (low bits first):", bin(alist[-1])[2:][::-1])
    print("    (k, modulus 2^(k-2), a):", [(k, 1 << (k - 2), alist[k - 3]) for k in range(3, 15)])
    # lift to k = 39: a mod 2^37 (the note quotes 68011179275 mod 2^37)
    ak = 1
    for k in range(4, 40):
        M = 1 << k
        cand = [ak, ak + (1 << (k - 3))]
        sols = [c for c in cand if pow(3, c, M) == (-5) % M]
        assert len(sols) == 1
        ak = sols[0]
    print("    a mod 2^37 (k = 39):", ak)
    check("2e log_3(-5) = 68011179275 mod 2^37", ak, 68011179275)
    # closure test == (c) restated: {3^a mod 2^k} = residues 1,3 mod 8, k <= 20 (done in (c))
    ints = [-1, -5, -7, -17, -25, -37, -55, -41, -61, -91, 1]
    inset = [u for u in ints if u % 8 in (1, 3)]
    outset = [u for u in ints if u % 8 not in (1, 3)]
    print("    integers = 1 or 3 mod 8 :", inset)
    print("    integers not            :", outset, " (residues mod 8:", [u % 8 for u in outset], ")")
    check("2e {-5,-7,-37,-55,-61,1} are 1 or 3 mod 8; the others are not",
          (inset, outset), ([-5, -7, -37, -55, -61, 1], [-1, -17, -25, -41, -91]))
    # also check these integers really are Collatz cycle points of T on negatives
    def T(n):
        return n // 2 if n % 2 == 0 else (3 * n + 1) // 2
    cyc = {}
    for u in ints:
        x, seen = u, []
        while x not in seen:
            seen.append(x)
            x = T(x)
        cyc[u] = (x == u)
    print("    each listed integer lies on a T-cycle:", all(cyc.values()), cyc)

    # (f) Lucas / Fibonacci identity
    print("\n(f) A = [[1,1],[1,0]]: |det(A^j - I)| = |L_j - 1 - (-1)^j|, tr(A^j) = L_j, j <= 60")
    L = [2, 1]
    for _ in range(61):
        L.append(L[-1] + L[-2])
    A = [[1, 1], [1, 0]]
    P = [[1, 0], [0, 1]]
    okf = True
    for j in range(1, 61):
        P = [[P[0][0] * A[0][0] + P[0][1] * A[1][0], P[0][0] * A[0][1] + P[0][1] * A[1][1]],
             [P[1][0] * A[0][0] + P[1][1] * A[1][0], P[1][0] * A[0][1] + P[1][1] * A[1][1]]]
        det = (P[0][0] - 1) * (P[1][1] - 1) - P[0][1] * P[1][0]
        tr = P[0][0] + P[1][1]
        okf &= (abs(det) == abs(L[j] - 1 - (-1) ** j)) and (tr == L[j])
    check("2f determinant and trace identities, j=1..60", okf, True)
    print("    |det(A^j - I)| for j=1..12:", [abs(L[j] - 1 - (-1) ** j) for j in range(1, 13)])
    print("    tr(A^j) - |det(A^j - I)| = 1 + (-1)^j: 2 for even j (the two period-2 codings), 0 for odd j")
    print(" ", elapsed())

    # (g) carry census (Probe B): words w of length A <= 18, x_0 = S_w / (2^A - 3^p) integral?
    print("\n(g) carry census: words of length A <= 18 whose carry S_w is 0 in Z/(2^A - 3^p)")
    hits = Counter()
    hits_by_shape = Counter()
    for A in range(1, 19):
        for w in range(1 << A):            # bit i of w = letter w_i (1 = odd step)
            p = bin(w).count("1")
            # S_w = sum_i w_i 2^i 3^(ones after i)
            S = 0
            ones_after = 0
            for i in range(A - 1, -1, -1):
                if (w >> i) & 1:
                    S += (1 << i) * 3 ** ones_after
                    ones_after += 1
            clock = (1 << A) - 3 ** p
            if S % clock == 0:
                x0 = S // clock
                hits[x0] += 1
                hits_by_shape[(A, p, x0)] += 1
    print("    integers x_0 hit (x_0: number of (A, word) pairs):", sorted(hits.items()))
    known = {0, 1, 2, -1, -5, -7, -10, -17, -25, -37, -55, -82, -41, -61, -91, -136, -68, -34}
    check("2g every integral carry solution with A <= 18 is a point of a known T-cycle (incl. 0 and even points)",
          set(hits) <= known, True)
    check("2g all 18 known cycle points (0, 1, 2, -1, the -5 and -17 cycles) occur", set(hits), known)
    sh18 = sorted((x0, c) for (A, p, x0), c in hits_by_shape.items() if (A, p) == (18, 9))
    print("    shape (18, 9), clock 2^18 - 3^9 = %d: hits %s" % ((1 << 18) - 3 ** 9, sh18))
    check("2g clock 242461 is hit only by the repeat of the cycle {1,2}", [x for x, _ in sh18], [1, 2])
    sh117 = sorted((x0, c) for (A, p, x0), c in hits_by_shape.items() if (A, p) == (11, 7))
    print("    shape (11, 7), clock 139: hits (x_0, count) =", sh117, " -> %d words = the 11 rotations (7 odd, 4 even points)"
          % sum(c for _, c in sh117))
    print(" ", elapsed())


if __name__ == "__main__":
    part1()
    part2()
    print("\n" + "=" * 100)
    print("SUMMARY: %d failed checks" % len(FAILS))
    for f in FAILS:
        print("  FAIL:", f)
    print(elapsed())
