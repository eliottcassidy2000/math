#!/usr/bin/env python3
"""
lrc14_threadC_support_parity_ocf_opus_s2.py    (opus-2026-06-20-S2, THREAD C)

QUESTION (THREAD C).  The tournament OCF H = signed sum over ODD-cycle
collections; EVEN cycles drop.  The LRC Weyl error is

    corr(E) = measS7(E) - iid_k = sum_{0 != n in Lambda(E)} K(n),
    Lambda(E) = { n in Z^{k-1} : sum_j n_j e_j = 0 }   (the OFFSET RELATION LATTICE),
    K(n)      = sum_{T subseteq {1..6}} (-1)^{|T|} prod_{j=1}^{k-1} chat_T(n_j),

where the e_j are the nonzero offsets and chat_T is the Fourier coefficient of the
"missed-sector" indicator (chat_T(0)=1-|T|/7, chat_T(7m)=0 = THM-503 7-vanishing,
|chat_T(n)| <= |T|/(pi|n|)).  See angleF_fourier_lattice (mac-mini-0618s7) for the
identity; it is VERIFIED there to reproduce measS7 exactly under truncation.

The SUPPORT of a relation is supp(n) = { j : n_j != 0 }.  This is the LRC analog of
a cycle's vertex/edge set: a relation sum_{j in supp} n_j e_j = 0 is an additive
"loop" among the speeds e_j (the cycle space of the speed multiset, the LRC twin of
the GF(2) cycle space of K_n that carries the OCF).

CLAIM TO TEST.  Does corr(E) stratify by support parity the way H drops even cycles?
Is the ODD-support part of corr dominant, and/or does the EVEN-support part cancel?

ENUMERATION (the key to making this exact/fast).  We do NOT brute-force a box in
Z^{d-1}.  We stratify BY SUPPORT directly: for each subset S of the nonzero speeds
with |S|=s, we enumerate the integer relations supported EXACTLY on S with bounded
coefficients (|n_j| <= N0) and sum K(n) over them.  Support 1 is impossible for a
primitive set (n_j e_j = 0 => n_j = 0).  Support 2 forces n_i e_i = -n_j e_j, i.e.
the unique primitive relation +/-(e_j/g, -e_i/g)*t -- a "2-body" cut-type term.
Support >=3 are the genuine many-body cycle terms (Schur triples a+b=c live here).

PARTS.
  1. Exact corr(E) = measS7 - iid (ground truth).
  2. Support-size + parity stratification of corr (truncated lattice sum).
  3. ODD-support sum vs EVEN-support sum (the OCF test).
  4. Per-support detail: Re K, count, 1/height weight.
  5. 7-vanishing (apex prime) filter per support layer.
Each finding is marked ANALOGY / PARTIAL-TOOL / TOOL in the write-up.
"""
from __future__ import annotations
import sys, itertools, math, cmath
from collections import defaultdict
from fractions import Fraction as F

sys.stdout.reconfigure(line_buffering=True)
TWO_PI_I = 2j * math.pi


# ----------------------------------------------------------------------------
# Ground truth.
# ----------------------------------------------------------------------------
def measS7(E):
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        if len({int(((e * xm) % 1) * 7) for e in E}) == 7:
            total += x1 - x0
    return total


def stirling2(n, k):
    return sum((-1) ** (k - j) * math.comb(k, j) * j ** n for j in range(k + 1)) // math.factorial(k)


def iid_k(k):
    return F(math.factorial(7) * stirling2(k, 7), 7 ** k)


# ----------------------------------------------------------------------------
# Fourier kernel.
# ----------------------------------------------------------------------------
def shat(n, j):
    if n == 0:
        return 1.0 / 7.0
    a = j / 7.0
    return (cmath.exp(-TWO_PI_I * n * a) - cmath.exp(-TWO_PI_I * n * (a + 1 / 7.0))) / (TWO_PI_I * n)


SUBSETS = [tuple(T) for r in range(0, 7) for T in itertools.combinations(range(1, 7), r)]
SIGN = {T: (-1) ** len(T) for T in SUBSETS}

_chat_cache: dict = {}


def chat_T(n, T):
    key = (n, T)
    v = _chat_cache.get(key)
    if v is not None:
        return v
    if n == 0:
        v = complex(1 - len(T) / 7.0, 0.0)
    elif n % 7 == 0:
        v = 0j  # THM-503 7-vanishing (apex prime)
    else:
        v = -sum(shat(n, j) for j in T)
    _chat_cache[key] = v
    return v


def K_on_support(coeffs, speeds_idx, nz):
    """K(n) where n is supported on positions speeds_idx with values coeffs.
    All other coordinates are 0.  Uses chat_T(0) on the zero coordinates."""
    d = len(nz)
    nvec = [0] * d
    for p, c in zip(speeds_idx, coeffs):
        nvec[p] = c
    K = 0j
    for T in SUBSETS:
        p = 1.0 + 0j
        ok = True
        for ni in nvec:
            c = chat_T(ni, T)
            if c == 0:
                ok = False
                break
            p *= c
        if ok:
            K += SIGN[T] * p
    return K


# ----------------------------------------------------------------------------
# Enumerate, by support, the integer relations supported EXACTLY on a chosen
# subset of speeds, with |coeff| <= N0, summing K.
# ----------------------------------------------------------------------------
def stratified_corr(E, N0, max_supp=None):
    nz = [e for e in E if e != 0]
    d = len(nz)
    if max_supp is None:
        max_supp = d
    by_K = defaultdict(complex)
    by_cnt = defaultdict(int)
    by_W = defaultdict(float)
    rng = [c for c in range(-N0, N0 + 1) if c != 0]
    for s in range(2, max_supp + 1):
        for idxs in itertools.combinations(range(d), s):
            sp = [nz[i] for i in idxs]
            last = sp[-1]
            # solve last coord: c_last = -(sum_{i<s-1} c_i sp_i)/sp[-1], must be
            # a nonzero integer with |c_last|<=N0.  Enumerate first s-1 freely.
            for coeffs0 in itertools.product(rng, repeat=s - 1):
                acc = sum(c * v for c, v in zip(coeffs0, sp[:-1]))
                if acc % last != 0:
                    continue
                cl = -acc // last
                if cl == 0 or abs(cl) > N0:
                    continue
                coeffs = coeffs0 + (cl,)
                K = K_on_support(coeffs, idxs, nz)
                by_K[s] += K
                by_cnt[s] += 1
                h = 1
                for c in coeffs:
                    h *= abs(c)
                by_W[s] += 1.0 / h
    total = sum(v for v in by_K.values())
    return by_K, by_cnt, by_W, total.real


# ----------------------------------------------------------------------------
SHAPES_K7 = [
    ("consec      {0..6}", [0, 1, 2, 3, 4, 5, 6]),
    ("near-consec {0..5,7}", [0, 1, 2, 3, 4, 5, 7]),
    ("AP step2    {0,2..12}", [0, 2, 4, 6, 8, 10, 12]),
    ("two-block   {0,1,2,3,40,41,42}", [0, 1, 2, 3, 40, 41, 42]),
    ("Schur-rich  {0,1,2,3,4,64,65}", [0, 1, 2, 3, 4, 64, 65]),
    ("dissoc      {0,1,3,7,15,31,63}", [0, 1, 3, 7, 15, 31, 63]),
]
SHAPES_K8 = [
    ("consec      {0..7}", list(range(8))),
    ("near-consec {0..6,8}", [0, 1, 2, 3, 4, 5, 6, 8]),
    ("two-block   {0,1,2,3,40,41,42,43}", [0, 1, 2, 3, 40, 41, 42, 43]),
    ("Schur-rich  {0,1,2,3,4,5,64,65}", [0, 1, 2, 3, 4, 5, 64, 65]),
    ("dissoc      {0,1,3,7,15,31,63,127}", [0, 1, 3, 7, 15, 31, 63, 127]),
]


def part1(shapes, k):
    print(f"PART 1 -- exact corr(E) = measS7 - iid_{k}   (iid_{k} = {float(iid_k(k)):.6f})")
    print(f"  {'shape':<34}{'measS7':>12}{'corr (exact)':>16}{'sign':>6}")
    for name, E in shapes:
        m = float(measS7(E))
        c = m - float(iid_k(k))
        print(f"  {name:<34}{m:>12.6f}{c:>16.6f}{('+' if c >= 0 else '-'):>6}")
    print()


def part23(shapes, k, N0, max_supp):
    print(f"PART 2/3 -- support + parity stratification of corr  (|coeff|<={N0}, supp<={max_supp})")
    print("  support 1 impossible (primitive); lowest is 2 = '2-body / cut' term.")
    print(f"  {'shape':<24}{'corr_trunc':>11}{'odd-supp':>11}{'even-supp':>11}"
          f"{'|odd|/|tot|':>11}{'lead supp':>10}")
    rows = []
    for name, E in shapes:
        bsK, bsc, bsW, tot = stratified_corr(E, N0, max_supp)
        odd = sum(v.real for s, v in bsK.items() if s % 2 == 1)
        even = sum(v.real for s, v in bsK.items() if s % 2 == 0)
        lead = max(bsK.items(), key=lambda kv: abs(kv[1])) if bsK else (None, 0j)
        oddfrac = abs(odd) / abs(tot) if abs(tot) > 1e-9 else float("nan")
        short = name.split()[0]
        print(f"  {short:<24}{tot:>11.6f}{odd:>11.6f}{even:>11.6f}"
              f"{oddfrac:>11.3f}{str(lead[0]):>10}")
        rows.append((name, E, bsK, bsc, bsW, tot, odd, even))
    print()
    return rows


def part_detail(rows):
    print("PART 4 -- per-support breakdown (Re K, count, 1/height weight)")
    for name, E, bsK, bsc, bsW, tot, odd, even in rows:
        print(f"  {name}   corr_trunc={tot:+.6f}   odd={odd:+.6f} even={even:+.6f}")
        for s in sorted(bsK):
            par = "ODD " if s % 2 == 1 else "EVEN"
            print(f"      supp={s} [{par}]  ReK={bsK[s].real:+.6f}  "
                  f"|K|={abs(bsK[s]):.6f}  count={bsc[s]:>7}  Wsum={bsW[s]:.4f}")
        print()


def part5_vanishing(shapes, N0, max_supp):
    print(f"PART 5 -- 7-vanishing (apex prime) filter per support  (|coeff|<={N0}, supp<={max_supp})")
    print("  a relation is 7-KILLED iff some coord is a nonzero multiple of 7")
    print("  (chat_T(7m)=0 for all T zeroes the whole product).")
    print(f"  {'shape':<22}{'supp':>5}{'total':>9}{'7-killed':>9}{'alive':>8}{'alive ReK':>12}")
    for name, E in shapes:
        nz = [e for e in E if e != 0]
        d = len(nz)
        per = defaultdict(lambda: [0, 0, 0.0])
        rng = [c for c in range(-N0, N0 + 1) if c != 0]
        for s in range(2, max_supp + 1):
            for idxs in itertools.combinations(range(d), s):
                sp = [nz[i] for i in idxs]
                last = sp[-1]
                for coeffs0 in itertools.product(rng, repeat=s - 1):
                    acc = sum(c * v for c, v in zip(coeffs0, sp[:-1]))
                    if acc % last != 0:
                        continue
                    cl = -acc // last
                    if cl == 0 or abs(cl) > N0:
                        continue
                    coeffs = coeffs0 + (cl,)
                    rec = per[s]
                    rec[0] += 1
                    if any(c % 7 == 0 for c in coeffs):
                        rec[1] += 1
                    else:
                        rec[2] += K_on_support(coeffs, idxs, nz).real
        short = name.split()[0]
        first = True
        for s in sorted(per):
            t, kld, alive = per[s]
            lbl = short if first else ""
            first = False
            print(f"  {lbl:<22}{s:>5}{t:>9}{kld:>9}{t - kld:>8}{alive:>12.6f}")
        print()


def main():
    print("=" * 92)
    print("THREAD C -- does the LRC Weyl error have OCF (odd-support) structure?")
    print("opus-2026-06-20-S2; exact corr; complex truncated lattice sums.")
    print("=" * 92)
    print()
    print("### k = 7  (6 nonzero speeds) ###")
    part1(SHAPES_K7, 7)
    rows7 = part23(SHAPES_K7, 7, N0=10, max_supp=6)
    part_detail(rows7)
    part5_vanishing(SHAPES_K7, N0=8, max_supp=5)
    print("### k = 8  (7 nonzero speeds) ###")
    part1(SHAPES_K8, 8)
    rows8 = part23(SHAPES_K8, 8, N0=8, max_supp=5)
    part_detail(rows8)


if __name__ == "__main__":
    main()
