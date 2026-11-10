#!/usr/bin/env python3
"""AUDIT check (d): analysis of the auditor's own deletion-pair tables (d_window.c).
usage: d_window_analyze.py DIR MAXORDER      (uses DIR/win_K.A.bin, .Z.bin for K = 2..9 and DIR/win_10_r.* for K = 10)
Prints, for bases of order <= 8, <= 9, ..., <= MAXORDER: realised h, exceptions to the windows [3h,5h) and [3h,9h),
the statistics quoted in section 6 of the README, and a comparison with the packet's out/window_le10.out."""
import sys, os, glob, re, statistics
D = sys.argv[1]; MAXO = int(sys.argv[2])
HM, NM = 1024, 16384; ROW = NM // 16
def load(fn):
    b = open(fn, "rb").read(); assert len(b) == (HM // 2) * ROW
    return [int.from_bytes(b[r * ROW:(r + 1) * ROW], "little") for r in range(HM // 2)]
def files(K, kind):
    f = glob.glob(os.path.join(D, f"win_{K}.{kind}.bin")) + glob.glob(os.path.join(D, f"win_{K}_*.{kind}.bin"))
    return sorted(f)
def has(T, h, n): return T[h >> 1] >> (n >> 1) & 1
def odd_range(lo, hi):          # odd integers n with lo <= n < hi (lo, hi real)
    import math
    n = math.ceil(lo); n += (n % 2 == 0)
    while n < hi:
        yield n; n += 2
A = [0] * (HM // 2); Z = [0] * (HM // 2)
here = os.path.dirname(os.path.abspath(__file__))
for K in range(2, MAXO + 1):
    fa, fz = files(K, "A"), files(K, "Z")
    if not fa: print(f"(no tables for order {K})"); continue
    for fn in fa:
        for r, x in enumerate(load(fn)): A[r] |= x
    for fn in fz:
        for r, x in enumerate(load(fn)): Z[r] |= x
    if K < 8: continue
    realised = [h for h in range(1, HM, 2) if has(A, h, h)]
    print(f"\n=== bases of order <= {K} ({len(fa)} table file(s) at order {K}) ===")
    print(f"realised odd h < 1024: {len(realised)}; unrealised below the largest realised: {[h for h in range(1, max(realised), 2) if h not in realised][:12]}{' ...' if len([h for h in range(1, max(realised), 2) if h not in realised]) > 12 else ''}; largest realised {max(realised)}")
    exc35 = {h: [n for n in odd_range(3 * h, 5 * h) if not has(A, h, n)] for h in realised}
    exc35 = {h: v for h, v in exc35.items() if v}
    exc39 = {h: [n for n in odd_range(3 * h, 9 * h) if not has(A, h, n)] for h in realised}
    exc39 = {h: v for h, v in exc39.items() if v}
    def complete_from(exc):       # least h0 such that every realised AND unrealised odd h in [h0, 1023] ... (needs all realised)
        bad = set(exc) | (set(range(1, HM, 2)) - set(realised))
        h0 = 1023
        while h0 >= 1 and h0 not in bad: h0 -= 2
        return h0 + 2
    print(f"(W1) [3h,5h): {len(exc35)} realised h with an incomplete window; largest such h = {max(exc35) if exc35 else None}; "
          f"complete (and h realised) for every odd h in [{complete_from(exc35)}, 1023]")
    print(f"(W2) [3h,9h): {len(exc39)} realised h with an incomplete window; largest such h = {max(exc39) if exc39 else None}; "
          f"complete (and h realised) for every odd h in [{complete_from(exc39)}, 1023]")
    if K == 8:
        # the auditor's own sub-statement for bases of order <= 8
        runs = []; cur = None
        for h in range(1, HM, 2):
            good = (h in realised) and (h not in exc35)
            if good and cur is None: cur = h
            if not good and cur is not None: runs.append((cur, h - 2)); cur = None
        if cur is not None: runs.append((cur, 1023))
        best = max(runs, key=lambda r: r[1] - r[0])
        print(f"   own sub-statement: with bases of order <= 8 the window [3h,5h) is complete for every odd h in [{best[0]}, {best[1]}] "
              f"({(best[1] - best[0]) // 2 + 1} consecutive odd h); all maximal runs of length >= 5: {[r for r in runs if r[1] - r[0] >= 8]}")
        small = {h: exc35[h] for h in sorted(exc35) if h <= 45}
        print(f"   incomplete windows for h <= 45 (order <= 8): {small}")
    # statistics of section 6 for 237 <= h <= 1023
    hs = [h for h in realised if h >= 237]
    if hs:
        ratios = []
        for h in hs:
            n = 5 * h - 2
            while n > h and has(A, h, n): n -= 2
            ratios.append(((n + 2) / h, h))
        print(f"full-from ratio n0(h)/h, realised 237 <= h <= 1023 ({len(hs)} values): min {min(ratios)[0]:.3f} (h={min(ratios)[1]}), median {statistics.median(r for r, _ in ratios):.3f}, max {max(ratios)[0]:.3f} (h={max(ratios)[1]})")
        least = []
        for h in hs:
            n = h + 2
            while not has(A, h, n): n += 2
            least.append((n / h, h, n))
        r, h, n = min(least)
        print(f"least extension ratio over realised h >= 237: {r:.4f} at h = {h}, n = {n}; forced bound (3h+1)/(2h) = {(3 * h + 1) / (2 * h):.4f}; "
              f"violations of n >= (3h+1)/2: {sum(1 for r_, h_, n_ in least if 2 * n_ < 3 * h_ + 1)}; pairs on the bound: {[(h_, n_) for r_, h_, n_ in least if 2 * n_ in (3 * h_ + 1, 3 * h_ + 3)][:8]}")
        def dens(T, lo, hi, hmin, hmax=1023):
            tot = pres = 0
            for h in realised:
                if hmin <= h <= hmax:
                    for n in odd_range(lo * h, hi * h): tot += 1; pres += has(T, h, n)
            return pres / tot if tot else float("nan")
        print(f"all out-sets, h >= 237: density [3h,5h) {dens(A, 3, 5, 237):.5f}; [2h,3h) {dens(A, 2, 3, 237):.5f}; [3h/2,2h) {dens(A, 1.5, 2, 237):.4f}")
        print(f"all out-sets, h >= 100: density [3h/2,2h) {dens(A, 1.5, 2, 100):.4f}")
        print(f"out-degree-one extensions, h >= 237: density [3h,5h) {dens(Z, 3, 5, 237):.4f}; [2h,3h) {dens(Z, 2, 3, 237):.4f}")
    # comparison with the packet's listing
    pk = os.path.join(here, "..", "out", "window_le10.out")
    if os.path.exists(pk):
        claimed = {}
        for l in open(pk):
            mm = re.match(r"\s*h=\s*(\d+) missing n: \[(.*)\]", l)
            if mm: claimed[int(mm.group(1))] = [int(x) for x in mm.group(2).split(",")]
        npairs = sum(len(v) for v in claimed.values())
        still = sum(1 for h, v in claimed.items() for n in v if not has(A, h, n))
        print(f"packet list of missing pairs in [3h,5h) (order <= 10): {len(claimed)} values of h, {npairs} pairs; of these, absent from the auditor's table at order <= {K}: {still} "
              f"(must be all {npairs}); auditor's exception list identical to the packet's: {exc35 == claimed}")
        if K == MAXO and exc35 != claimed:
            only_a = {h: v for h, v in exc35.items() if claimed.get(h) != v}; only_p = {h: v for h, v in claimed.items() if exc35.get(h) != v}
            print("   differences (auditor):", dict(list(only_a.items())[:10])); print("   differences (packet):", dict(list(only_p.items())[:10]))
