#!/usr/bin/env python3
"""Independent audit script for collatz_artin_corrections_20260927 (Propositions 1
and 2, the effective-sample-size reading of the STICKY audit's 0.497, the 3-adic
marginal above the start range, the seeds' numbers, Artin's 20/19).

Everything is exact integer arithmetic except the permutation tests (numpy, own
seed) and the floating summaries.  Own code throughout; no import of the note's
script.

Run from the repo root:
  python 04-computation/experiments/collatz_artin_corrections_20260927_audit.py \
      > 05-knowledge/results/collatz_artin_corrections_20260927_audit.out
"""
from __future__ import annotations

import math
import random
import sys
import time
from collections import Counter, defaultdict

import numpy as np

T0 = time.time()


def stamp() -> str:
    return f"[{time.time() - T0:7.1f}s]"


def v2(x: int) -> int:
    return (x & -x).bit_length() - 1


def U(m: int) -> tuple[int, int]:
    x = 3 * m + 1
    v = v2(x)
    return x >> v, v


# ----------------------------------------------------------------------------
# A. forward pass: successor map and visit weights by reverse-topological sums
# ----------------------------------------------------------------------------

def forward(N: int):
    """nxt[m] = U(m) for every odd value visited by some odd start n <= N (m != 1);
    w[m] = number of odd starts n <= N whose Syracuse orbit passes through m
    (the start itself and the terminal 1 included, as in the note's script)."""
    nxt = {}
    for n in range(1, N + 1, 2):
        m = n
        while m != 1 and m not in nxt:
            x = 3 * m + 1
            m2 = x >> v2(x)
            nxt[m] = m2
            m = m2
    visited = set(nxt)
    visited.add(1)
    indeg = Counter()
    for b in nxt.values():
        indeg[b] += 1
    w = {m: (1 if m <= N else 0) for m in visited}
    stack = [m for m in visited if indeg[m] == 0]
    while stack:
        a = stack.pop()
        if a == 1:
            continue
        b = nxt[a]
        w[b] += w[a]
        indeg[b] -= 1
        if indeg[b] == 0:
            stack.append(b)
    assert all(indeg[m] == 0 for m in visited), "cycle among visited values?"
    return nxt, w


def band_table(N: int, nxt, w, min_distinct: int = 200, verbose: bool = True):
    """Per dyadic band, over visited m = 3 (mod 4) (i.e. v(m) = 1):
    R_w = sum w [U(m) = 3 mod 4] / sum w, R_1 = distinct-value version,
    n_eff = (sum w)^2 / sum w^2, plus the finite-population s.e."""
    acc = defaultdict(lambda: [0, 0, 0, 0, 0])  # sum w, sum w*both, count, count both, sum w^2
    for m, wt in w.items():
        if m == 1 or m % 4 != 3:
            continue
        i = m.bit_length() - 1
        both = 1 if nxt[m] % 4 == 3 else 0
        a = acc[i]
        a[0] += wt
        a[1] += wt * both
        a[2] += 1
        a[3] += both
        a[4] += wt * wt
    rows = {}
    for i in sorted(acc):
        sw, swb, c, cb, sw2 = acc[i]
        if c < min_distinct:
            continue
        neff = sw * sw / sw2
        Rw = swb / sw
        R1 = cb / c
        se_iid = 0.5 / math.sqrt(neff)
        # finite population (exactly half of the labels are 1, labels permuted): Var = (1/4) (n/(n-1)) (1/neff - 1/n)
        se_fp = 0.5 * math.sqrt((c / (c - 1)) * (1 / neff - 1 / c)) if c > 1 else float("nan")
        tag = "below N" if 2 ** (i + 1) <= N else ("straddles N" if 2 ** i <= N else "above N")
        z = (Rw - 0.5) / se_iid
        rows[i] = (tag, Rw, R1, neff, c, se_iid, se_fp, z, sw, swb, cb)
        if verbose:
            print(f"   i={i:>2} ({tag:>11}): R_w = {Rw:.4f}  R_1 = {R1:.6f} ({cb}/{c})  n_eff = {neff:8.0f}  sum w = {sw:9d}"
                  f"  2 s.e.(iid) = {2 * se_iid:.4f}  2 s.e.(finite pop.) = {2 * se_fp:.4f}  z = {z:+.2f}")
    return rows


def part_A(N: int):
    print(f"== A. forward pass, all odd starts n <= {N} (own code: early-stopping walk + reverse-topological visit weights) ==")
    nxt, w = forward(N)
    H = max(w)
    tot = sum(w.values())
    print(f"   distinct odd values visited: {len(w)}; total visits (sum of weights): {tot}; max weight {w[1]} at m = 1; "
          f"largest visited odd value H = {H} ({H.bit_length()} bits)")
    # Syracuse depth (number of odd steps to 1) of every visited value; the depth J needed in Proposition 1 for a value m is the
    # longest path from a start <= N to m, at most the maximal depth
    depth = {1: 0}
    for m in list(w):
        path = []
        x = m
        while x not in depth:
            path.append(x)
            x = nxt[x]
        d = depth[x]
        for y in reversed(path):
            d += 1
            depth[y] = d
    dmax = max(depth[n] for n in range(1, N + 1, 2))
    argm = max((n for n in range(1, N + 1, 2)), key=lambda n: depth[n])
    print(f"   maximal Syracuse depth (odd steps to 1) over starts n <= N: {dmax} at n = {argm}; log_3 N = {math.log(N, 3):.1f}; "
          f"3^(J+1) with J = {dmax} has {int((dmax + 1) * math.log10(3)) + 1} decimal digits")
    print(stamp())
    return nxt, w, H


# ----------------------------------------------------------------------------
# B. Proposition 1(a): the odd preimages of m are exactly {(2^k m - 1)/3 : k >= 1, 2^k m = 1 mod 3}
# ----------------------------------------------------------------------------

def part_B(N: int, nxt, w, H: int):
    print("\n== B. Proposition 1(a): preimage set and the tree recursion, checked on EVERY visited value ==")
    children = defaultdict(list)
    for a, b in nxt.items():
        children[b].append(a)
    visited = set(w)
    bad_set = 0
    bad_rec = 0
    mult3_with_children = 0
    k0_issue = 0
    checked = 0
    for m in visited:
        if m == 1:
            continue
        r = m % 3
        actual = set(children.get(m, ()))
        if r == 0:
            if actual:
                mult3_with_children += 1
            continue
        k = 2 if r == 1 else 1  # 2^k m = 1 (mod 3): k even for m = 1 (3), k odd for m = 2 (3)
        formula = set()
        while True:
            a = (2 ** k * m - 1)
            if a > 3 * H:
                break
            assert a % 3 == 0
            a //= 3
            assert a % 2 == 1
            if a in visited:
                formula.add(a)
            k += 2
        checked += 1
        if formula != actual:
            bad_set += 1
        if w[m] != (1 if m <= N else 0) + sum(w[a] for a in actual):
            bad_rec += 1
        if r == 1 and ((m - 1) // 3) % 2 == 0:
            k0_issue += 1
    print(f"   values checked (3 does not divide m): {checked}; children(m) from the forward map == formula preimages "
          f"(k >= 1, one parity class, restricted to visited values): mismatches {bad_set}")
    print(f"   recursion w(m) = [m <= N] + sum over children w(a): failures {bad_rec}")
    print(f"   multiples of 3 with a child (should be 0: multiples of 3 are leaves): {mult3_with_children}")
    print(f"   nit: for m = 1 (mod 3) the note's set {{(2^k m - 1)/3 : 2^k m = 1 mod 3}} without 'k >= 1' would include k = 0, "
          f"i.e. (m-1)/3, which is EVEN for every such m ({k0_issue} of them here), not an odd preimage; the set needs k >= 1")
    print(stamp())


# ----------------------------------------------------------------------------
# C. Proposition 1(b): the affine form of the ancestor and the bound on the carry c
# ----------------------------------------------------------------------------

def part_C(N: int, nxt, w, sample: int = 3000, seed: int = 8675309):
    print("\n== C. Proposition 1(b): ancestor = (2^K m - c)/3^j with c = sum_{l=1}^{j} 3^(l-1) 2^(K_j - K_l); the note's bound 0 < c < 2^K ==")
    rng = random.Random(seed)
    starts = [rng.randrange(1, N + 1) | 1 for _ in range(sample)]
    pairs = 0
    viol = 0
    viol_words = Counter()
    max_ratio_extremal = 0.0
    max_c_over_2Km = 0.0
    rel_shift = []
    first_viol = None
    for n in starts:
        vals = []
        m = n
        orbit = [n]
        while m != 1:
            x = 3 * m + 1
            v = v2(x)
            vals.append(v)
            m = x >> v
            orbit.append(m)
        # for each prefix n -> m_j: inverse word k_1 = v_{j-1}, ..., k_j = v_0
        for j in range(1, len(vals) + 1):
            ks = vals[:j][::-1]
            K = 0
            Kl = []
            for k in ks:
                K += k
                Kl.append(K)
            c = sum(3 ** (l - 1) * 2 ** (K - Kl[l - 1]) for l in range(1, j + 1))
            mj = orbit[j]
            num = 2 ** K * mj - c
            assert num % 3 ** j == 0 and num // 3 ** j == n, (n, j)
            pairs += 1
            if c >= 2 ** K:
                viol += 1
                if first_viol is None:
                    first_viol = (n, j, ks[:6], c, K)
                viol_words[tuple(ks[:2])] += 1
            max_ratio_extremal = max(max_ratio_extremal, c / (2 ** K * ((1.5) ** j - 1)))
            max_c_over_2Km = max(max_c_over_2Km, c / (2 ** K * mj))
            rel_shift.append(c / (2 ** K * mj))
    rel_shift.sort()
    print(f"   {sample} random odd starts, {pairs} (start, orbit value) pairs: n = (2^K m_j - c)/3^j exact in every case")
    print(f"   the note's bound c < 2^K FAILS in {viol} of {pairs} pairs ({viol / pairs:.1%}); first failure: start {first_viol[0]}, j = {first_viol[1]}, "
          f"inverse word begins {first_viol[2]}, c = {first_viol[3]} vs 2^K = {2 ** first_viol[4]}")
    print(f"   e.g. the inverse word (1, 1): ancestor (4m - 5)/9, c = 5 > 4 (m = 17: 17 <- 11 <- 7). Failures by the first two inverse letters: "
          f"{dict(sorted(viol_words.items(), key=lambda t: -t[1])[:6])}")
    print(f"   correct bounds: c/2^K = sum_l 3^(l-1) 2^(-K_l) <= (3/2)^j - 1 (equality iff all k_l = 1): max observed ratio to that bound {max_ratio_extremal:.4f};"
          f" and c < 2^K m (positivity of the ancestor): max c/(2^K m) = {max_c_over_2Km:.4f}")
    q = [rel_shift[int(f * (len(rel_shift) - 1))] for f in (0.5, 0.9, 0.99, 1.0)]
    print(f"   relative threshold shift c/(2^K m) = 1 - 3^j n/(2^K m) (the 'rounding'): median {q[0]:.2e}, 90% {q[1]:.2e}, 99% {q[2]:.2e}, max {q[3]:.4f}")
    print(f"   (the shift is a function of the word only, not of m; it is ~ sum over the path of 1/(3 a_l), large only through small intermediate values)")
    print(stamp())


# ----------------------------------------------------------------------------
# D. the table, both conventions, n_eff, sign test, permutation tests
# ----------------------------------------------------------------------------

def part_D(N: int, nxt, w, do_perm: bool = True):
    print(f"\n== D. visit-weighted vs distinct-value persistence per dyadic band, N = {N} (band by the value m with v(m) = 1; target v(U m) = 1) ==")
    rows = band_table(N, nxt, w)
    # aggregates with their own n_eff
    for tag, cond in (("m <= N", lambda m: m <= N), ("m > N", lambda m: m > N), ("all m", lambda m: True)):
        sw = swb = c = cb = sw2 = 0
        for m, wt in w.items():
            if m == 1 or m % 4 != 3 or not cond(m):
                continue
            both = 1 if nxt[m] % 4 == 3 else 0
            sw += wt
            swb += wt * both
            c += 1
            cb += both
            sw2 += wt * wt
        neff = sw * sw / sw2
        print(f"   aggregate {tag:>6}: R_w = {swb / sw:.5f}, R_1 = {cb}/{c} = {cb / c:.5f}, n_eff = {neff:.0f}, sum w = {sw}, 2 s.e.(iid) = {1 / math.sqrt(neff):.4f}, "
              f"z = {(swb / sw - 0.5) * 2 * math.sqrt(neff):+.2f}")
    # sign pattern and combined deviation over the bands with n_eff >= 100
    sel = [(i, r) for i, r in rows.items() if r[3] >= 100]
    below = sum(1 for i, r in sel if r[1] < 0.5)
    above = sum(1 for i, r in sel if r[1] > 0.5)
    num = sum(r[3] * (r[1] - 0.5) for i, r in sel)
    den = sum(r[3] for i, r in sel)
    comb = num / den
    se = 0.5 / math.sqrt(den)
    print(f"   bands with n_eff >= 100: {len(sel)}; R_w below 1/2 in {below}, above in {above}; n_eff-weighted mean deviation {comb:+.5f} "
          f"(s.e. {se:.5f}, z = {comb / se:+.2f})")
    mid = [(i, r) for i, r in sel if 14 <= i <= 23]
    below = sum(1 for i, r in mid if r[1] < 0.5)
    num = sum(r[3] * (r[1] - 0.5) for i, r in mid)
    den = sum(r[3] for i, r in mid)
    print(f"   bands 14..23 (the STICKY audit's range 10^4..10^7): {len(mid)} bands, {below} below 1/2; n_eff-weighted mean deviation {num / den:+.5f} "
          f"(s.e. {0.5 / math.sqrt(den):.5f}, z = {(num / den) / (0.5 / math.sqrt(den)):+.2f}); two-sided sign-test P(>= {max(below, len(mid) - below)} of {len(mid)} on one side) = "
          f"{2 * sum(math.comb(len(mid), t) for t in range(max(below, len(mid) - below), len(mid) + 1)) / 2 ** len(mid):.3f}")
    # z-scores summary
    zs = [r[7] for i, r in sel]
    print(f"   z-scores (iid model) of the {len(sel)} bands: max |z| = {max(abs(z) for z in zs):.2f}; mean z = {sum(zs) / len(zs):+.2f}; "
          f"sum z^2 / k = {sum(z * z for z in zs) / len(zs):.2f} (1 expected under the model)")
    if do_perm:
        rng = np.random.default_rng(8675309)
        print("   permutation test of the note's model (labels m mod 8 permuted among the m = 3 mod 4 of the band, weights fixed; 4000 permutations):")
        for i in (12, 16, 18, 22, 24, 26, 27):
            if i not in rows:
                continue
            ms = [m for m in w if m != 1 and m % 4 == 3 and m.bit_length() - 1 == i]
            wt = np.array([w[m] for m in ms], dtype=np.float64)
            lab = np.array([1.0 if m % 8 == 7 else 0.0 for m in ms])
            obs = float((wt * lab).sum() / wt.sum())
            dev = abs(obs - 0.5)
            cnt = 0
            perms = 4000
            sims = np.empty(perms)
            for t in range(perms):
                rng.shuffle(lab)
                sims[t] = (wt * lab).sum() / wt.sum()
            p = float((np.abs(sims - 0.5) >= dev - 1e-15).mean())
            print(f"      band {i:>2}: R_w = {obs:.4f}, permutation s.d. {sims.std():.4f} (model s.e. iid {rows[i][5]:.4f}, finite pop. {rows[i][6]:.4f}), two-sided p = {p:.3f}")
    print(stamp())
    return rows


def part_D2(N: int, nxt, w):
    """The STICKY audit's statistic in both banding conventions, recovered from the weights."""
    print(f"\n== D2. the STICKY audit's orbit statistic recovered from the weights (N = {N}); bands [1,10^4), [10^4,10^5), [10^5,10^6), [10^6,inf) ==")
    cuts = [10 ** 4, 10 ** 5, 10 ** 6]

    def band(x):
        b = 0
        for c in cuts:
            if x >= c:
                b += 1
        return b

    names = ["[1,10^4)", "[10^4,10^5)", "[10^5,10^6)", "[10^6,inf)"]
    # PREDECESSOR banding (= the note's convention): weight w(m) on m = 3 mod 4, band by m, target U(m) = 3 mod 4
    P = [[0, 0, 0] for _ in range(4)]
    # CURRENT banding (= the audit's headline numbers): band by x = U(m); each x = 2 mod 3 has the unique v=1 preimage m = (2x-1)/3
    C = [[0, 0, 0] for _ in range(4)]
    for m, wt in w.items():
        if m == 1 or m % 4 != 3:
            continue
        x = nxt[m]
        both = 1 if x % 4 == 3 else 0
        b = band(m)
        P[b][0] += wt
        P[b][1] += wt * both
        P[b][2] += wt * wt
        b = band(x)
        C[b][0] += wt
        C[b][1] += wt * both
        C[b][2] += wt * wt
    for label, A in (("banded by the PREDECESSOR m (note's convention)", P), ("banded by the CURRENT value U(m) (audit's headline)", C)):
        print(f"   {label}:")
        for b in range(4):
            sw, swb, sw2 = A[b]
            neff = sw * sw / sw2
            print(f"      {names[b]:>12}: P_w = {swb / sw:.4f} on {sw} steps; n_eff = {neff:.0f}; 2 s.e.(iid) = {1 / math.sqrt(neff):.4f}; z = {(swb / sw - 0.5) * 2 * math.sqrt(neff):+.2f}")
        hi = [sum(A[b][k] for b in (1, 2, 3)) for k in range(3)]
        print(f"      at or above 10^4: {hi[1] / hi[0]:.4f} ({hi[0]} steps), n_eff = {hi[0] ** 2 / hi[2]:.0f}, z = {(hi[1] / hi[0] - 0.5) * 2 * math.sqrt(hi[0] ** 2 / hi[2]):+.2f}")
    print(stamp())


# ----------------------------------------------------------------------------
# E. Proposition 2 and the 3-adic marginal above N; in-valuation laws
# ----------------------------------------------------------------------------

def part_E(N: int, nxt, w):
    print(f"\n== E. Proposition 2 (m = U(m') is 2^(-v) mod 3) and the 3-adic marginal / arrival-valuation law above N = {N} ==")
    bad = 0
    for m, x in nxt.items():
        v = v2(3 * m + 1)
        if x % 3 != (2 if v % 2 else 1):
            bad += 1
    print(f"   U(m) = 2 (mod 3) iff v odd, = 1 (mod 3) iff v even: violations on all {len(nxt)} visited edges: {bad}")
    mod3 = Counter()
    for m in w:
        if m > N:
            mod3[m % 3] += 1
    tot = sum(mod3.values())
    print(f"   distinct visited values above N: {tot}; residues mod 3 (0, 1, 2): {mod3[0] / tot:.4f}, {mod3[1] / tot:.4f}, {mod3[2] / tot:.4f} (Haar-valuation prediction 0, 1/3, 2/3)")
    # arrival edges into x > N
    preds = defaultdict(list)
    for m, x in nxt.items():
        if x > N:
            preds[x].append(m)
    parity_bad = 0
    n_first = 0
    n_noext = 0
    noext_mod2 = 0
    lawE = Counter()   # per distinct edge
    lawV = Counter()   # per visit (weight w[m])
    lawE_above = Counter()
    lawV_above = Counter()
    odd_first_frac = 0
    for x, ms in preds.items():
        ks = [v2(3 * m + 1) for m in ms]
        if len({k % 2 for k in ks}) > 1:
            parity_bad += 1
        first = any(m <= N for m in ms)
        if first:
            n_first += 1
        else:
            n_noext += 1
            if x % 3 == 2:
                noext_mod2 += 1
        for m, k in zip(ms, ks):
            kk = min(k, 9)
            lawE[kk] += 1
            lawV[kk] += w[m]
            if m > N:
                lawE_above[kk] += 1
                lawV_above[kk] += w[m]
    print(f"   distinct values above N with a visited predecessor: {len(preds)} of {tot} (every value above N is entered, none is a start); "
          f"predecessors of one value all have the same valuation parity: violations {parity_bad}")
    lo = (2 * N - 1) // 3 + 1
    lo += (3 - lo) % 4
    print(f"   values above N entered from below N (only by v = 1, m' in ((2N-1)/3, N]): {n_first} ({n_first / len(preds):.4f}); exact count of m' = 3 (4) in that range: "
          f"{len(range(lo, N + 1, 4))} (U restricted to v = 1 is injective)")
    print(f"   values above N entered ONLY from above N: {n_noext}; of these = 2 (mod 3) (odd arrival valuation): {noext_mod2 / n_noext:.4f} (Haar 2/3)")
    for label, L in (("per distinct arrival edge", lawE), ("per visit (edge weighted by w(m'))", lawV),
                     ("per distinct edge, both ends above N", lawE_above), ("per visit, both ends above N", lawV_above)):
        s = sum(L.values())
        podd = sum(c for k, c in L.items() if k % 2 == 1) / s
        print(f"   arrival valuation law {label:>40}: P(k) for k = 1..6 = {[round(L[k] / s, 4) for k in range(1, 7)]}, P(k odd) = {podd:.4f} (Haar 1/2, 1/4, ...; 2/3)")
    # density profile of distinct visited values above N by band, and the out-valuation law (2-adic, should be Haar)
    cnt = Counter()
    for m in w:
        if m > N:
            cnt[m.bit_length() - 1] += 1
    bands = sorted(cnt)
    print(f"   distinct visited values above N per dyadic band {bands[0]}..{bands[0] + 8}: {[cnt[i] for i in bands[:9]]}; ratios {[round(cnt[i + 1] / cnt[i], 3) for i in bands[:8]]}")
    # density-profile reading: if the density of visited values at size x above N is rho(x) ~ x^(-alpha), then x = 2 (3) is visited
    # iff one of its odd-k preimages (2^k x - 1)/3 (k odd) is visited, x = 1 (3) iff one of its even-k preimages is; with independent
    # visits of probability rho((2^k x)/3) the ratio of the two chances is sum_(k odd) 2^(-k alpha) / sum_(k even >= 2) 2^(-k alpha) = 2^alpha,
    # so P(x = 2 mod 3 | x visited above N) = 2^alpha/(1 + 2^alpha) (Haar alpha = 1 gives 2/3; alpha = 2 gives 4/5),
    # and the arrival law is P(k) proportional to 2^(-k alpha) within each parity class.
    lo_b, hi_b = bands[1], bands[1] + 5
    alpha_count = -math.log2(cnt[hi_b] / cnt[lo_b]) / (hi_b - lo_b)  # counts per band ~ 2^(-alpha_count i); density exponent = alpha_count + 1
    alpha = alpha_count + 1
    pred = 2 ** alpha / (1 + 2 ** alpha)
    pk = {k: (2 ** (-k * alpha)) for k in range(1, 9)}
    sodd = sum(v for k, v in pk.items() if k % 2)
    sev = sum(v for k, v in pk.items() if k % 2 == 0)
    predlaw = [round((pred * pk[k] / sodd) if k % 2 else ((1 - pred) * pk[k] / sev), 4) for k in range(1, 7)]
    print(f"   density-profile reading: counts per band fall like 2^(-{alpha_count:.3f} i) over bands {lo_b}..{hi_b}, i.e. density ~ x^(-{alpha:.3f}); "
          f"predicted P(x = 2 mod 3 | visited above N) = 2^alpha/(1 + 2^alpha) = {pred:.4f} (measured {mod3[2] / tot:.4f}; Haar alpha = 1: 0.6667; alpha = 2: 0.8000)")
    print(f"   predicted arrival law per distinct edge P(k), k = 1..6: {predlaw} (measured above); under the Haar model the Cramer exponent of the drift, "
          f"E[(3/2^v)^theta] = 1, is theta* = 1 exactly (3^theta = 2^(theta+1) - 1), so excursion heights have tail 1/h and the visited density above N falls like x^(-2)")
    lawO = Counter()
    lawOw = Counter()
    for m, wt in w.items():
        if m > N:
            k = min(v2(3 * m + 1), 8)
            lawO[k] += 1
            lawOw[k] += wt
    s1 = sum(lawO.values())
    s2 = sum(lawOw.values())
    print(f"   out-valuation law of visited m > N, distinct: {[round(lawO[k] / s1, 4) for k in range(1, 8)]}; weighted: {[round(lawOw[k] / s2, 4) for k in range(1, 8)]}; 2^-k: {[round(2.0 ** -k, 4) for k in range(1, 8)]}")
    drift1 = sum(math.log2(3) - v2(3 * m + 1) for m in w if m > N) / s1
    driftW = sum(wt * (math.log2(3) - v2(3 * m + 1)) for m, wt in w.items() if m > N) / s2
    print(f"   mean drift log2(3) - v above N: distinct {drift1:+.4f}, weighted {driftW:+.4f}, Haar {math.log2(3) - 2:+.4f}")
    joint = Counter()
    for m in w:
        if m > N:
            joint[(m % 8, m % 9)] += 1
    p8 = Counter()
    p9 = Counter()
    for (a, b), c in joint.items():
        p8[a] += c
        p9[b] += c
    mi = sum(c / tot * math.log2((c / tot) / ((p8[a] / tot) * (p9[b] / tot))) for (a, b), c in joint.items())
    df = (len(p8) - 1) * (len(p9) - 1)
    print(f"   mod 8 marginal (1,3,5,7): {[round(p8[a] / tot, 4) for a in (1, 3, 5, 7)]}; I(m mod 8; m mod 9) = {mi:.6f} bits; "
          f"chance level under independence df/(2 n ln 2) = {df / (2 * tot * math.log(2)):.6f} bits (df = {df})")
    print(stamp())


# ----------------------------------------------------------------------------
# F. (d): R_1 = 1/2 exactly? residue counts in dyadic bands and in [1, N]
# ----------------------------------------------------------------------------

def part_F(N: int):
    print("\n== F. (d) exactness of R_1 = 1/2 below N: counting m = 7 (8) among m = 3 (4) ==")
    for i in (10, 12, 14, 16, 18):
        lo, hi = 2 ** i, 2 ** (i + 1)
        n3 = len(range(lo + 3, hi, 4))
        n7 = len(range(lo + 7, hi, 8))
        print(f"   band [2^{i}, 2^{i + 1}): #m = 3 (4) = {n3}, #m = 7 (8) = {n7}, ratio {n7 / n3}")
    n3 = len(range(3, N + 1, 4))
    n7 = len(range(7, N + 1, 8))
    print(f"   [1, {N}]: {n7}/{n3} = {n7 / n3} (exact because 8 | N); general interval [a, b): ratio = 1/2 + O(1/#): e.g. [1, 10^6 + 4]: {len(range(7, N + 5, 8))}/{len(range(3, N + 5, 4))}")
    print(stamp())


# ----------------------------------------------------------------------------
# G. the seeds
# ----------------------------------------------------------------------------

def T(n: int) -> int:
    return n * (n + 1) // 2


def part_G():
    print("\n== G. section 4's numbers ==")
    # symbolic pivot identities
    try:
        import sympy as sp
        n, k = sp.symbols("n k", integer=True, positive=True)
        lin = sp.summation(n ** 2 + k, (k, 0, n)) - sp.summation(n ** 2 + k, (k, n + 1, 2 * n))
        c2 = 2 * n ** 2 + n
        sq = sp.summation((c2 + k) ** 2, (k, 0, n)) - sp.summation((c2 + k) ** 2, (k, n + 1, 2 * n))
        c = sp.symbols("c", integer=True)
        gen = sp.expand(sp.summation((c + k) ** 2, (k, 0, n)) - sp.summation((c + k) ** 2, (k, n + 1, 2 * n)))
        roots = sp.solve(sp.Eq(gen, 0), c)
        cub = sp.expand(sp.summation((c + k) ** 3, (k, 0, n)) - sp.summation((c + k) ** 3, (k, n + 1, 2 * n)))
        print(f"   symbolic: linear identity difference = {sp.simplify(lin)}; square identity (start 2n^2+n) difference = {sp.simplify(sq)}")
        print(f"   general square identity: the starts c with sum_(k<=n) (c+k)^2 = sum_(n<k<=2n) (c+k)^2 are c = {roots} (pivot c + n = 4T_n for c = 2n^2 + n)")
        print(f"   cubic analogue, difference as a polynomial in c: {sp.collect(cub, c)}")
    except Exception as e:  # pragma: no cover
        print(f"   sympy unavailable or failed ({e}); numeric checks below")
    for nn in range(1, 41):
        c = nn * nn
        assert sum(range(c, c + nn + 1)) == sum(range(c + nn + 1, c + 2 * nn + 1)) and c + nn == 2 * T(nn)
        c2 = nn * (2 * nn + 1)
        assert sum(x * x for x in range(c2, c2 + nn + 1)) == sum(x * x for x in range(c2 + nn + 1, c2 + 2 * nn + 1))
        assert c2 + nn == 4 * T(nn) and c2 == T(2 * nn)
    print("   numeric: both families and pivots 2T_n, 4T_n and start T_(2n) hold for n <= 40")
    # cubes: exact search c < 20000 for n <= 6, then ALL integers c for n <= 300 via the integer polynomial's real roots
    def cube_poly(nn: int, c: int) -> int:
        return sum((c + k) ** 3 for k in range(nn + 1)) - sum((c + k) ** 3 for k in range(nn + 1, 2 * nn + 1))
    found = {nn: [c for c in range(1, 20000) if cube_poly(nn, c) == 0] for nn in range(1, 7)}
    print(f"   cube family, c < 20000, n <= 6: {found}")
    none = True
    for nn in range(1, 301):
        # coefficients of c^3 - 3 n^2 c^2 - 3 n^2 (2n+1) c - n^2 (7n^2 + 6n + 1)/2  (constant term = 2 T_n^2 - T_(2n)^2)
        coef = [1, -3 * nn ** 2, -3 * nn ** 2 * (2 * nn + 1), -(nn ** 2 * (7 * nn ** 2 + 6 * nn + 1)) // 2]
        assert cube_poly(nn, 5) == sum(cf * 5 ** (3 - i) for i, cf in enumerate(coef))
        for r in np.roots(coef):
            if abs(r.imag) < 1e-6:
                c0 = int(round(r.real))
                for c in (c0 - 1, c0, c0 + 1):
                    if cube_poly(nn, c) == 0:
                        none = False
                        print(f"      integer solution n = {nn}, c = {c}")
    print(f"   cube family for ANY integer c (positive or negative), n <= 300: {'none' if none else 'FOUND'} (largest real root is about 3n^2 + 2n, so c < 20000 already covered n <= 6)")
    print(f"   sporadic 3^3 + 4^3 + 5^3 = 6^3: {27 + 64 + 125 == 216}")
    # square triangular numbers: brute force n <= 10^4 and the Pell orbit
    brute = [(nn, math.isqrt(T(nn))) for nn in range(1, 10 ** 4 + 1) if math.isqrt(T(nn)) ** 2 == T(nn)]
    pell = []
    x, y = 3, 2
    while len(pell) < 6:
        assert x * x - 2 * y * y == 1 and x % 2 == 1 and y % 2 == 0
        pell.append(((x - 1) // 2, y // 2))
        x, y = 3 * x + 4 * y, 2 * x + 3 * y
    print(f"   square triangular: brute force n <= 10^4: {brute}; Pell orbit of 3 + 2 sqrt 2: {pell}; T_n = {[T(a) for a, b in pell]}; agree: {brute == pell}")
    # Schur
    from itertools import product

    def sumfree_colourings(n: int, allow_equal: bool = True):
        good = []
        for col in product((0, 1), repeat=n):
            ok = True
            for a in range(1, n + 1):
                for b in range(a if allow_equal else a + 1, n + 1):
                    s = a + b
                    if s <= n and col[a - 1] == col[b - 1] == col[s - 1]:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                good.append(col)
        return good
    g4 = sumfree_colourings(4)
    g5 = sumfree_colourings(5)
    classes = [tuple(i + 1 for i in range(4) if c[i] == c[0]) for c in g4]
    print(f"   Schur (x + y = z, x = y allowed): sum-free 2-colourings of [1,4]: {len(g4)}, class of 1: {classes}; of [1,5]: {len(g5)}; "
          f"(with x != y required the counts are {len(sumfree_colourings(4, False))} and {len(sumfree_colourings(5, False))}: the '2' needs Schur's x = y convention)")
    # HKM 7825: Euclid parametrisation (all triples with c <= M) and a numpy brute force
    M = 7825
    in_triple = np.zeros(M + 1, dtype=bool)
    with7825 = set()
    ntrip = 0
    for a in range(2, math.isqrt(M) + 1):
        for b in range(1, a):
            if (a - b) % 2 == 1 and math.gcd(a, b) == 1:
                x, y, z = a * a - b * b, 2 * a * b, a * a + b * b
                if z > M:
                    continue
                kk = 1
                while kk * z <= M:
                    ntrip += 1
                    t = (kk * x, kk * y, kk * z)
                    in_triple[list(t)] = True
                    if M in t:
                        with7825.add(tuple(sorted(t)))
                    kk += 1
    print(f"   Pythagorean triples with hypotenuse <= {M}: {ntrip}; numbers in [1, {M}] in no triple: {M - int(in_triple[1:].sum())}; triples containing {M}: {len(with7825)} {sorted(with7825)}")
    # brute force cross-check
    in2 = np.zeros(M + 1, dtype=bool)
    cnt7825 = 0
    for a in range(1, M + 1):
        bmax = math.isqrt(M * M - a * a)
        if bmax < a:
            break
        b = np.arange(a, bmax + 1, dtype=np.int64)
        c2 = a * a + b * b
        c = np.sqrt(c2).astype(np.int64)
        c = np.where(c * c < c2, c + 1, c)
        ok = c * c == c2
        if ok.any():
            bb = b[ok]
            cc = c[ok]
            in2[a] = True
            in2[bb] = True
            in2[cc] = True
            cnt7825 += int((cc == M).sum()) + int((bb == M).sum()) + (int(ok.sum()) if a == M else 0)
    print(f"   brute-force cross-check: numbers in no triple: {M - int(in2[1:].sum())}; triples containing {M}: {cnt7825}; 7825 = 5^2 * 313: {7825 == 25 * 313}; "
          f"as a hypotenuse, 7825^2 = a^2 + b^2 has ((2*2+1)(2*1+1) - 1)/2 = {((5 * 3) - 1) // 2} representations (313 = 1 mod 4: {313 % 4 == 1})")
    print(stamp())


# ----------------------------------------------------------------------------
# H. Artin
# ----------------------------------------------------------------------------

def primes_upto(n: int) -> np.ndarray:
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, math.isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def part_H():
    print("\n== H. Artin's constant, the primitive-root densities for 2 and 5, and the 20/19 factor ==")
    P = primes_upto(10 ** 7)
    Pf = P.astype(np.float64)
    A = float(np.prod(1 - 1 / (Pf * (Pf - 1))))
    print(f"   A = prod_(p < 10^7) (1 - 1/(p(p-1))) = {A:.7f} (tail < 10^-7); A * 20/19 = {A * 20 / 19:.5f}; A * 19/20 = {A * 19 / 20:.5f} (the owner's '19/20' is the reciprocal factor)")
    # 20/19 from the general formula: sum over squarefree n of mu(n)/[Q(zeta_n, 5^(1/n)):Q], degree halved when 10 | n
    # extra = (1/40) * prod_(q != 2,5) (1 - 1/(q(q-1))) = A/19
    print(f"   check of the correction: extra mass (1/40) prod_(q not in {{2,5}}) (1 - 1/(q(q-1))) = A/((1/2)(19/20)(40)) = A/19: {A / ((0.5) * (19 / 20) * 40):.7f} vs A/19 = {A / 19:.7f}")
    # densities
    spf = np.zeros(10 ** 6 + 1, dtype=np.int64)
    for p in primes_upto(10 ** 6):
        sl = spf[p::p]
        sl[sl == 0] = p

    def is_prim(a: int, p: int) -> bool:
        q = p - 1
        fs = set()
        while q > 1:
            f = int(spf[q])
            fs.add(f)
            while q % f == 0:
                q //= f
        return all(pow(a, (p - 1) // f, p) != 1 for f in fs)
    for bound in (3 * 10 ** 5, 10 ** 6):
        c2 = c5 = tot = 0
        c2n = totn = 0  # the note's convention: primes 3 <= p < bound, p != 5, for both bases
        for p in primes_upto(bound - 1):
            p = int(p)
            if p == 2:
                continue
            if p != 5:
                totn += 1
                c2n += is_prim(2, p)
                c5 += is_prim(5, p)
            tot += 1
            c2 += is_prim(2, p)
        print(f"   primes < {bound}: base 2 (all odd p): {c2 / tot:.4f}; base 2 (note's convention, p != 5): {c2n / totn:.4f}; base 5 (p != 2, 5): {c5 / totn:.4f}; "
              f"A = {A:.4f}, A*20/19 = {A * 20 / 19:.4f}")
    print(stamp())


# ----------------------------------------------------------------------------

def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 6
    print(f"audit script start; N = {N}")
    nxt, w, H = part_A(N)
    part_B(N, nxt, w, H)
    part_C(N, nxt, w)
    rows = part_D(N, nxt, w)
    part_D2(N, nxt, w)
    part_E(N, nxt, w)
    part_F(N)
    # N-dependence of the band deviations (smaller N, no permutation tests)
    print("\n== D3. the same band table for smaller N (does |R_w - 1/2| track 1/sqrt(n_eff)?) ==")
    for N2 in (10 ** 5, 3 * 10 ** 5):
        nxt2, w2 = forward(N2)
        print(f"   N = {N2}:")
        rows2 = band_table(N2, nxt2, w2)
        zs = [r[7] for i, r in rows2.items() if r[3] >= 100]
        print(f"   bands with n_eff >= 100: {len(zs)}; max |z| = {max(abs(z) for z in zs):.2f}; mean z = {sum(zs) / len(zs):+.2f}; sum z^2 / k = {sum(z * z for z in zs) / len(zs):.2f}")
        print(stamp())
    part_G()
    part_H()
    print(f"\n{stamp()} DONE")


if __name__ == "__main__":
    main()
