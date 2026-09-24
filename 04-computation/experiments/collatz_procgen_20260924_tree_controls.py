#!/usr/bin/env python3
"""collatz_procgen_20260924_tree_controls.py

Lane: inverse tree mod 192 (session collatz-procgen-20260922, 2026-09-24).

Part 2: the controls behind the no-go.

  C1  SHEET: the minus sheet T_-(n) = n/2, (3n-1)/2.  Positive cycles (min <= N), basin census,
      basin densities and their residue profiles mod 3, 9, 8, 192.  Type-level statistics are
      identical to the plus sheet's after r -> -r (checked on the exceptional-class sets).
  C2  The integer points of the multiplicative exceptional set:
      Bad_+ cap [-N, N]  and  Bad_- cap [-N, N]  (Bad = parity sequences with 3^a > 2^j for every prefix).
  C3  DEFECT: two planted maps that agree with T on every class up to finitely many integers per level,
      differ only on a set of counting function O(log X), and have a divergent orbit.
  C4  Z-Collatz on [-N, N]: every integer's T_+-orbit ends in one of the five integer cycles
      {0}, {1,2}, {-1}, {-5,-7,-10}, {-17,...,-34}  (equivalently T_- on [-N, N], by negation).

Every check raises on failure.  Peak memory about 310 MB (N = 10^7 minus-sheet census); runtime about 5 s.
"""
import sys
import time
from math import log2

import numpy as np

T0 = time.time()
LOG2_3 = log2(3.0)


def check(cond, msg):
    if not cond:
        raise AssertionError("CHECK FAILED: " + msg)


def T(n, b=1):
    return n // 2 if n % 2 == 0 else (3 * n + b) // 2


def mem(tag):
    import resource
    print(f"[mem] {tag}: max RSS so far {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.0f} MB", file=sys.stderr)


def hdr(s):
    mem("before " + s[:3])
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


def cycles_up_to(N, b):
    """All T_b-cycles with minimum <= N among positive integers, by the descent census."""
    cyc = []
    seen_min = set()
    for n0 in range(1, N + 1):
        x = T(n0, b)
        steps = 1
        # follow until below n0 or back at n0
        while x > n0 and steps < 100000:
            x = T(x, b)
            steps += 1
        if x == n0:
            c = [n0]
            y = T(n0, b)
            while y != n0:
                c.append(y)
                y = T(y, b)
            if min(c) not in seen_min:
                seen_min.add(min(c))
                cyc.append(c)
        check(steps < 100000, "no runaway below the cap")
    return cyc


def basins(N, b, cycles, chunk=500_000):
    """Vectorized basin labels for 1..N (chunked in increasing n, in-place arithmetic; low memory).
    cycles: list of cycles (lists of positive ints).  A value is resolved when its orbit first
    falls below it onto an already-resolved value."""
    lab = np.full(N + 1, -1, dtype=np.int8)
    for ci, c in enumerate(cycles):
        for v in c:
            if v <= N:
                lab[v] = ci
    rounds_max = 0
    for lo in range(1, N + 1, chunk):
        hi = min(N, lo + chunk - 1)
        act = np.arange(lo, hi + 1, dtype=np.int64)
        act = act[lab[lo:hi + 1] < 0]
        x = act.copy()
        rounds = 0
        while act.size:
            rounds += 1
            odd = (x & 1).astype(bool)
            xo = x[odd]
            xo *= 3
            xo += b
            xo >>= 1
            x >>= 1
            x[odd] = xo
            del xo, odd
            below = x < act
            if below.any():
                idx = np.nonzero(below)[0]
                lv = lab[x[idx]]
                ok = lv >= 0
                if ok.any():
                    good = idx[ok]
                    lab[act[good]] = lv[ok]
                    keep = np.ones(act.size, dtype=bool)
                    keep[good] = False
                    act = act[keep]
                    x = x[keep]
            check(rounds < 5000, "basin loop terminates")
        rounds_max = max(rounds_max, rounds)
    return lab, rounds_max


# ---------------------------------------------------------------------------
hdr("C1  SHEET control: the minus sheet has three positive cycles; its type data equal the plus sheet's")
# ---------------------------------------------------------------------------
cyc_m = cycles_up_to(10 ** 5, -1)
cyc_p = cycles_up_to(10 ** 5, 1)
cyc_m.sort(key=min)
print(f"  plus sheet  T_+ cycles with minimum <= 10^5: {[sorted(c) for c in cyc_p]}")
print(f"  minus sheet T_- cycles with minimum <= 10^5: {[sorted(c) for c in cyc_m]}")
check([sorted(c) for c in cyc_p] == [[1, 2]], "plus: only {1,2}")
check([min(c) for c in cyc_m] == [1, 5, 17] and [len(c) for c in cyc_m] == [1, 3, 11], "minus: 1, 5, 17")
for c in cyc_m:
    L = len(c)
    a = sum(1 for v in c if v % 2 == 1)
    print(f"    minus cycle min {min(c):3d}: {L:2d} T-steps, {a} odd steps, multiplier 3^{a}/2^{L} = {3**a/2**L:.6f} > 1")
    check(3 ** a > 2 ** L, "minus positive cycles expand")
print("    plus cycle {1,2}: 2 T-steps, 1 odd step, multiplier 3/4 < 1 (contracting; sign law: positive plus")
print("    cycles have 2^K > 3^L, positive minus cycles have 3^L > 2^K).")

for N in (10 ** 6, 10 ** 7):
    t1 = time.time()
    lab, rounds = basins(N, -1, cyc_m)
    check((lab[1:] >= 0).all(), f"every n <= {N} lies in one of the three minus basins")
    cnt = np.bincount(lab[1:], minlength=3)
    print(f"  N = {N:>8d}: minus basins of 1, 5, 17 contain {cnt.tolist()} integers;"
          f" densities {[round(float(c) / N, 5) for c in cnt]}   ({rounds} rounds)")
    if N == 10 ** 7:
        # residue profiles of the three basins (chunked, low memory)
        for Mod in (3, 9, 8, 192):
            tab = np.zeros((3, Mod), dtype=np.int64)
            for lo in range(1, N + 1, 10 ** 6):
                hi = min(N, lo + 10 ** 6 - 1)
                nn = np.arange(lo, hi + 1, dtype=np.int64)
                tab += np.bincount(lab[lo:hi + 1].astype(np.int64) * Mod + nn % Mod,
                                   minlength=3 * Mod).reshape(3, Mod)
            tot = tab.sum(axis=0)
            prof = tab / tot
            if Mod in (3, 9, 8):
                print(f"    basin share by class mod {Mod}:")
                for ci in range(3):
                    print(f"      basin({[1,5,17][ci]:2d}): " + " ".join(f"{v:.4f}" for v in prof[ci]))
            # chi-square against the basin's global share (uniform profile), Mod-1 degrees of freedom
            chis = []
            for ci in range(3):
                pglob = tab[ci].sum() / tot.sum()
                expc = pglob * tot
                chis.append(float(((tab[ci] - expc) ** 2 / (expc * (1 - pglob))).sum()))
            print(f"    mod {Mod:3d}: chi-square of each basin's residue profile against uniform "
                  f"({Mod-1} d.o.f.): " + ", ".join(f"basin({[1,5,17][ci]}) {chis[ci]:.1f}" for ci in range(3)))
            check(all(c < 3 * (Mod - 1) + 30 for c in chis), "basins residue-equidistributed (no gross bias)")
    del lab

print()
print("  Type data are sheet-blind (PROVED in the note, Theorem 6 (transport)).  Finite check: the forward")
print("  multiplicative exceptional classes mod 2^k (no prefix with 3^a < 2^j within k steps) satisfy")
print("  Bad_k(minus) = -Bad_k(plus) as sets of residues, k = 1..18:")


def bad_classes(k, b):
    out = []
    for r in range(2 ** k):
        x = r if r > 0 else 2 ** k
        a = 0
        ok = True
        for j in range(1, k + 1):
            if x % 2:
                a += 1
            x = T(x, b)
            if 3 ** a < 2 ** j:
                ok = False
                break
        if ok:
            out.append(r)
    return out


BAD16 = {}
for k in range(1, 19):
    bp = bad_classes(k, 1)
    bm = bad_classes(k, -1)
    if k == 16:
        BAD16[1], BAD16[-1] = bp, bm
    check(sorted((-r) % 2 ** k for r in bp) == bm, f"Bad_k(minus) = -Bad_k(plus) at k={k}")
    if k in (6, 12, 18):
        print(f"    k = {k:2d}: |Bad_k| = {len(bp)} on both sheets, and the sets correspond under r -> -r")

# ---------------------------------------------------------------------------
hdr("C2  Integer points of the multiplicative exceptional set on both sheets")
# ---------------------------------------------------------------------------
print("x in Bad_b iff every prefix of the T_b-parity sequence of x has 3^a > 2^j.  For an integer whose")
print("orbit enters a cycle, membership is decided by following the orbit one full period past the entry")
print("(the cycle multiplier then repeats), so the check below is exact.")


def bad_integer_points(N, b, cycle_pts, period_max, prefilter, expanding):
    """Positive integers n <= N (on sheet b) whose T_b parity sequence stays above the line forever.
    prefilter: the residues mod 2^16 of Bad_16 on this sheet (a necessary condition).
    expanding: True if every cycle in cycle_pts has multiplier > 1.  Then an orbit that stays above
    the line for one full period after entering its cycle stays above forever (the multiplier only
    grows by the cycle factor per period).  For contracting cycles no such shortcut exists and the
    loop runs until every candidate has dropped below the line (which then must happen)."""
    tbl = np.zeros(65536, dtype=bool)
    tbl[np.array(prefilter, dtype=np.int64)] = True
    parts = []
    for lo in range(1, N + 1, 10 ** 6):
        nn = np.arange(lo, min(N, lo + 10 ** 6 - 1) + 1, dtype=np.int64)
        parts.append(nn[tbl[nn & 65535]])
    n = np.concatenate(parts)
    del parts
    x = n.copy()
    a = np.zeros(n.size, dtype=np.int64)
    since = np.full(n.size, -1, dtype=np.int64)
    cp = np.array(sorted(cycle_pts), dtype=np.int64)
    idx = np.arange(n.size)
    j = 0
    done_alive = []
    while idx.size:
        j += 1
        xa = x[idx]
        odd = (xa & 1).astype(bool)
        a[idx] += odd
        x[idx] = np.where(odd, (3 * xa + b) >> 1, xa >> 1)
        # above the line?  3^a > 2^j  <=>  a*log2(3) > j   (never equal; margin check below)
        margin = a[idx] * LOG2_3 - j
        check(np.all(np.abs(margin) > 1e-9), "no float ambiguity")
        died = margin < 0
        idx = idx[~died]
        # cycle bookkeeping: once inside a cycle, one full period decides membership
        inc = np.isin(x[idx], cp)
        newly = inc & (since[idx] < 0)
        since[idx[newly]] = j
        if expanding:
            finished = (since[idx] >= 0) & (j - since[idx] >= period_max + 1)
        else:
            finished = np.zeros(idx.size, dtype=bool)
        done_alive.extend(n[idx[finished]].tolist())
        idx = idx[~finished]
        check(j < 20000, "loop bound")
    return sorted(done_alive), int(n.size)


cp_m = sorted({v for c in cyc_m for v in c})
cp_p = [1, 2]
Nb = 10 ** 7
bad_minus_pos, pre_m = bad_integer_points(Nb, -1, cp_m, 11, BAD16[-1], expanding=True)
bad_plus_pos, pre_p = bad_integer_points(Nb, 1, cp_p, 2, BAD16[1], expanding=False)
print(f"  (prefilter: n mod 2^16 in Bad_16 leaves {pre_m} resp. {pre_p} candidates below 10^7)")
print(f"  Bad_- cap [1, 10^7] = {bad_minus_pos}      (= -(Bad_+ cap [-10^7, -1]) by negation)")
print(f"  Bad_+ cap [1, 10^7] = {bad_plus_pos}")
check(bad_minus_pos == [1, 5, 17], "Bad_- positive integers = {1, 5, 17}")
check(bad_plus_pos == [], "no positive integer in Bad_+ below 10^7")
print("  Hence (FINITE-EXACT): Bad_+ cap [-10^7, 10^7] = {-17, -5, -1}  and  Bad_- cap [-10^7, 10^7] = {1, 5, 17}.")
print("  The statement 'Bad cap Z = the points of least |x| on the three expanding integer cycles' is sheet-symmetric;")
print("  only the SIGN of those three integers tells the sheets apart.")

# ---------------------------------------------------------------------------
hdr("C3  DEFECT controls: planted density-zero modifications invisible to the type data")
# ---------------------------------------------------------------------------
print("Control P1 (one trit).  T1(n) = 2n if n = 3*2^m (m >= 1), T1(n) = T(n) otherwise.")
print("  * S1 = {3*2^m : m >= 1} is T1-invariant and T1^-1(S1) is inside S1 (T^-1 of a multiple of 3 is its")
print("    double only), so every orbit starting outside S1 is its T-orbit; S1 itself is a divergent orbit.")
print("  * T1(3*2^m) - T(3*2^m) = 9*2^(m-1): at level k (types mod 3*2^k) the type map of T1 equals that")
print("    of T except at the k-1 integers 3*2^m, m < k.  #S1 cap [1,X] = floor(log2(X/3)).")
S1 = [3 * 2 ** m for m in range(1, 40)]
for m in range(1, 39):
    s = 3 * 2 ** m
    check(2 * s - T(s) == 9 * 2 ** (m - 1), "discrepancy")
    for k in range(1, m + 1):
        check((2 * s - T(s)) % (3 * 2 ** (k - 1)) == 0, "agreement mod 3*2^(k-1) for k <= m")
# T^-1(S1) inside S1
for s in S1[:30]:
    check(s % 3 == 0, "multiple of 3: only the D-preimage")
print("  verified: discrepancy 9*2^(m-1), agreement mod 3*2^(k-1) for all k <= m (m < 39).")
print()
print("Control P2 (all 2- and 3-adic levels).  s_1 = 12 and s_(i+1) = s_i/2 + 6^i k_i, with the least k_i")
print("  making 2^(i+2) | s_(i+1) and s_(i+1) > s_i^2.  T2(s_i) = s_(i+1); T2 = T elsewhere.")
s = [12]
ks = []
for i in range(1, 7):
    si = s[-1]
    check(si % (2 ** (i + 1)) == 0 and si % 3 == 0, "invariant v2(s_i) >= i+1, 3 | s_i")
    k = 1
    # least k with 2^(i+2) | si/2 + 6^i k and si/2 + 6^i k > si^2
    base = si // 2
    kmin = max(1, (si * si - base) // (6 ** i) + 1)
    k = kmin
    while (base + 6 ** i * k) % (2 ** (i + 2)) != 0:
        k += 1
    s.append(base + 6 ** i * k)
    ks.append(k)
for i in range(1, len(s)):
    si, sn = s[i - 1], s[i]
    check((sn - T(si)) % (6 ** i) == 0, "T2(s_i) = T(s_i) mod 6^i")
    check(sn > si * si and sn % 3 == 0, "growth and 3 | s_(i+1)")
    # s_(i+1) is never a power-of-two multiple of an earlier s_j (so the D-chains stay disjoint)
    for sj in s[:i]:
        q, r = divmod(sn, sj)
        check(not (r == 0 and q & (q - 1) == 0), "S2 elements not on each other's D-chains")
print(f"  s_1..s_4 = {s[:4]}; s_5 has {len(str(s[4]))} digits, s_7 has {len(str(s[6]))} digits.")
print("  verified: T2(s_i) = T(s_i) mod 6^i (so the discrepancy tends to 0 in Z_2 AND in Z_3), s_(i+1) > s_i^2,")
print("  3 | s_i, and the D-chains {s_i 2^j} are pairwise disjoint.  The affected set is the union of the")
print("  D-chains of the s_i: O(log X * log log X) integers up to X.  T2 has the divergent orbit s_1, s_2, ...")
# statistics comparison for n <= 10^6: stopping times identical off the affected set
N = 10 ** 6
aff = set()
for si in s:
    x = si
    while x <= N:
        aff.add(x)
        x *= 2
aff |= {v for v in S1 if v <= N}
print(f"  Affected integers <= 10^6 (both controls together): {len(aff)}.  Every statement about T that")
print("  tolerates exceptions on a density-zero set of starting points (densities, log-densities, 'almost all',")
print("  tree-count lower bounds, residue-class statements with density-zero exceptions, and the one-step type")
print("  maps F_k up to finitely many exceptions per level) holds verbatim for T1 and T2; Collatz fails for both.")

# ---------------------------------------------------------------------------
hdr("C4  Z-Collatz on [-N, N]: the sheet-symmetric statement that contains Collatz")
# ---------------------------------------------------------------------------
N = 10 ** 6
labm, _ = basins(N, -1, cyc_m)
check((labm[1:] >= 0).all(), "every minus-sheet n <= 10^6 ends in 1, 5 or 17")
labp, _ = basins(N, 1, [[1, 2]])
check((labp[1:] >= 0).all(), "every plus-sheet n <= 10^6 ends in {1,2}")
print("  T_+ on [-10^6, 10^6]: 0 is fixed; every positive integer reaches {1,2}; every negative integer x")
print("  reaches -1, {-5,-7,-10} or {-17,...,-34} (checked as T_-(-x) = -T_+(x) on [1, 10^6]).")
cyc_all = {"0": [0], "{1,2}": [1, 2], "{-1}": [-1], "{-5,-7,-10}": [-5, -7, -10],
           "{-17,...}": [-17, -25, -37, -55, -82, -41, -61, -91, -136, -68, -34]}
print("  The five integer cycles of T_+ on [-10^6, 10^6], with parity words read from the listed point:")
for name, c in cyc_all.items():
    x = c[0]
    w = []
    for _ in range(len(c)):
        w.append(x % 2)
        x = T(x, 1)
    check(x == c[0], "cycle closes")
    p_, a_ = len(w), sum(w)
    kind = "contracting" if 2 ** p_ > 3 ** a_ else "expanding"
    sgn = "positive" if c[0] > 0 else ("zero" if c[0] == 0 else "negative")
    check((kind == "contracting") == (c[0] >= 0), "sign law: contracting <=> nonnegative (plus sheet)")
    print(f"    {name:14s} word {''.join(map(str, w)):12s} p={p_:2d} a={a_}  {kind:11s} {sgn}")
print("  (i) every integer orbit in [-10^6, 10^6] is bounded; (ii) the contracting integer cycles are exactly")
print("  0 and {1,2}.  Both are sheet-symmetric statements; the sign law (contracting <=> positive on the plus")
print("  sheet) turns them into Collatz (Theorem 15 of the note).")
print("  Z-Collatz (all integers, both signs, five cycles named by their parity words) is invariant under")
print("  negation: it is the same statement for T_+ and for T_-.  Z-Collatz plus the sign check '{1,2} is the")
print("  only one of the five cycles inside Z_{>0}' implies Collatz (the check is true for T_+ and false for T_-,")
print("  whose positive cycles are 1, 5, 17).  Z-Collatz is strictly stronger than Collatz: it contains the")
print("  3n-1 conjecture on the positive integers.")

print()
print(f"[tree_controls] ALL CHECKS PASSED   ({time.time() - T0:.1f} s)", file=sys.stderr)
print("[tree_controls] ALL CHECKS PASSED")
