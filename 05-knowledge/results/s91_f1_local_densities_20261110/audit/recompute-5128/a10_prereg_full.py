# Auditor recompute-5128, task 6 (complete): independent re-evaluation of the eight pre-registered tests (279 items)
# on [2^24, 2^27), from my own census and my own predictions (rule of THM-5128 (G) / PREREGISTRATION.md).
# The tolerance rules are those written in PREREGISTRATION.md (and in the docstring of the session's window_split.py);
# the measured and predicted numbers are mine.  Also evaluated: the same tests with the standard error computed from
# the measured relative variance of r on each sheet (instead of the registered constant 0.75), and with 3 se.
import sys
import math
import time
import numpy as np
from numba import njit
from concurrent.futures import ThreadPoolExecutor
from a0_defs import psi_table32
from a6_window import predict

PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
          109, 113, 127, 131, 241, 251, 257, 263, 677, 683, 691, 2729, 2731, 2741, 8191, 43691, 65521, 65537, 65539]
SLOW = {3: 0.08, 7: 0.15, 31: 0.05, 73: 0.05, 127: 0.05}
SPECIAL = [241, 257, 683, 2731, 43691, 65537]
CONTROLS = [251, 263, 677, 691, 2729, 2741, 8191, 65521, 65539]
SIGN_HIGH = [3, 5, 7, 11, 17, 257, 65537]
SIGN_LOW = [29, 31, 37, 61, 67, 71, 73, 79, 83, 97, 101, 103, 109, 113, 127, 131]
WEIGHT_H = [0, 1, 2, 3, 8, 16, 64, 128, 1024]


@njit(cache=True)
def div_mask64(psi, all_p):
    n = psi.shape[0]
    mask = np.zeros(n, dtype=np.uint64)
    first = np.full(all_p.shape[0], -1, dtype=np.int64)
    for j in range(all_p.shape[0]):
        p = all_p[j]
        h = (p - 1) // 2
        while h < n:
            x = 2 * np.int64(psi[h]) + 1
            if x % p == 0:
                mask[h] |= np.uint64(1) << np.uint64(j)
                if first[j] < 0:
                    first[j] = h
            h += p
    return mask, first


@njit(cache=True, nogil=True)
def window_part(psi, mask, a, B, all_p, trkidx, ntrk, res, nres):
    na = all_p.shape[0]
    T = np.zeros(2, dtype=np.int64)
    dist = np.zeros((2, na), dtype=np.int64)
    anyd = np.zeros((2, na), dtype=np.int64)
    trk = np.zeros((ntrk, 2), dtype=np.int64)
    viol = 0
    lo = np.int64(1) << a
    for s in range(0, B):
        sheet = 0 if s == 0 else 1
        corr = 1 if s == 0 else 0
        Xs = (np.int64(1) << B) - 1 + corr
        h1 = np.int64(0)
        while True:
            y1 = 2 * h1 + 1
            y2crude = (4 * Xs) // ((1 + (np.int64(1) << (s + 2))) * y1)
            if (y2crude - 1) // 2 < h1:
                break
            if h1 % nres != res:
                h1 += 1
                continue
            x1 = 2 * np.int64(psi[h1]) + 1
            b = y1 << s
            y2max = (2 * Xs) // (x1 + 2 * b)
            h2max = (y2max - 1) // 2
            # lower bound for h2: M <= y2 (x1 + b) must reach lo  =>  y2 >= lo / (x1 + b)
            y2min = lo // (x1 + b)
            h2 = (y2min - 1) // 2 - 1
            if h2 < h1:
                h2 = h1
            m1 = mask[h1]
            t1 = -1
            if h1 <= 1024:
                t1 = trkidx[h1]
            while h2 <= h2max:
                M = x1 * (2 * np.int64(psi[h2]) + 1) + b * (2 * h2 + 1)
                if M <= Xs:
                    N = M - corr
                    if N >= lo:
                        w = 2 if h2 > h1 else 1
                        T[sheet] += w
                        m2 = mask[h2]
                        mm = m1 | m2
                        for j in range(na):
                            isd = (M % all_p[j] == 0)
                            if isd:
                                dist[sheet, j] += w
                            if (mm >> np.uint64(j)) & np.uint64(1):
                                anyd[sheet, j] += w
                                if not isd:
                                    viol += 1
                        if t1 >= 0:
                            trk[t1, sheet] += 1
                        if h2 != h1 and h2 <= 1024:
                            t2 = trkidx[h2]
                            if t2 >= 0:
                                trk[t2, sheet] += 1
                h2 += 1
            h1 += 1
    return T, dist, anyd, trk, viol


def main():
    a = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    B = int(sys.argv[2]) if len(sys.argv) > 2 else 27
    nthreads = 3
    t0 = time.time()
    all_p = np.array(PRIMES, dtype=np.int64)
    psi = psi_table32((1 << B) // 2 + 4)
    mask, first = div_mask64(psi, all_p)
    trkidx = -np.ones(1025, dtype=np.int64)
    for i, h in enumerate(WEIGHT_H):
        trkidx[h] = i
    with ThreadPoolExecutor(max_workers=nthreads) as ex:
        futs = [ex.submit(window_part, psi, mask, a, B, all_p, trkidx, len(WEIGHT_H), r, nthreads) for r in range(nthreads)]
        parts = [f.result() for f in futs]
    T = sum(p[0] for p in parts)
    dist = sum(p[1] for p in parts)
    anyd = sum(p[2] for p in parts)
    trk = sum(p[3] for p in parts)
    viol = sum(p[4] for p in parts)
    print(f"# census of [2^{a}, 2^{B}) done in {time.time() - t0:.1f} s; triples sheet 0 {int(T[0])}, sheet 1 {int(T[1])}; violations of confinement {viol}")
    A_, B_ = float(2 ** a), float(2 ** B)
    Wtot, Wdiv, Wtrk, nterms = predict(psi, all_p, A_, B_, B + 1, np.array(WEIGHT_H, dtype=np.int64))
    nodd = (2 ** B - 2 ** a) // 2
    # measured relative variance of r on each sheet (a7_misses.out): 0.4459, 0.3900 on [2^24, 2^27);
    # on another window: index of dispersion 2.78 / 3.08 (stable across windows) divided by the mean of r on the sheet
    if (a, B) == (24, 27):
        VREL = (0.4459, 0.3900)
    else:
        VREL = (2.78 / (int(T[0]) / nodd), 3.08 / (int(T[1]) / nodd))
    print(f"# relative variance used for the measured-variance se: {VREL[0]:.4f}, {VREL[1]:.4f}; mass per odd N: measured {int(T.sum()) / nodd:.5f}, predicted {float(Wtot.sum()) / nodd:.5f}; sheets measured {int(T[0]) / nodd:.4f}, {int(T[1]) / nodd:.4f}, predicted {Wtot[0] / nodd:.4f}, {Wtot[1] / nodd:.4f}")

    def se_reg(p, sheet):
        return math.sqrt(0.75 * p / nodd)

    def se_meas(p, sheet):
        return math.sqrt(VREL[sheet] * p / nodd)

    meas = {}
    pred = {}
    for j, p in enumerate(PRIMES):
        for sheet in (0, 1):
            d = int(dist[sheet, j])
            t = int(T[sheet])
            E = d / ((t - d) / (p - 1)) - 1
            ad = int(anyd[sheet, j])
            End = (d - ad) / ((t - d) / (p - 1)) - 1
            meas[(p, sheet)] = (E, ad / t, End, ad)
            f = Wdiv[sheet, j] / Wtot[sheet]
            D = 2 * f - 1 / p ** 2
            pred[(p, sheet)] = (p * D / (1 - D), 2 * f - f * f, p * (f * f - 1 / p ** 2) / (1 - D))

    def run_tests(se, nse, label):
        results = {}
        worst = {}

        def rec(test, name, dev, tol):
            ok = dev <= tol
            results.setdefault(test, []).append((name, dev, tol, ok))

        # A
        tot_m, tot_p = int(T.sum()), float(Wtot.sum())
        rec("A", "total", abs(tot_m / tot_p - 1), 0.002)
        for sheet in (0, 1):
            rec("A", f"sheet {sheet}", abs(int(T[sheet]) / Wtot[sheet] - 1), 0.02)
        # B
        for t, h in enumerate(WEIGHT_H):
            for sheet in (0, 1):
                rec("B", f"h1={h} sheet {sheet}", abs(int(trk[t, sheet]) / Wtrk[t, sheet] - 1), 0.005)
        # C: 31 primes up to 131 and 241, 257, 683
        for j, p in enumerate(PRIMES):
            if p <= 131 or p in (241, 257, 683):
                for sheet in (0, 1):
                    fm, fp = meas[(p, sheet)][1], pred[(p, sheet)][1]
                    cnt = meas[(p, sheet)][3]
                    tol = (0.03 if first[j] < 2 ** 12 else 0.15) + 4 / math.sqrt(max(cnt, 1))
                    rec("C", f"p={p} sheet {sheet}", abs(fm / fp - 1), tol)
        # D
        for p in PRIMES:
            if (p <= 131 and p not in SLOW) or p in SPECIAL:
                for sheet in (0, 1):
                    E, Ep = meas[(p, sheet)][0], pred[(p, sheet)][0]
                    rec("D", f"p={p} sheet {sheet}", abs(E - Ep), 0.008 + 0.08 * abs(Ep) + nse * se(p, sheet))
        # E
        for p, tol in SLOW.items():
            for sheet in (0, 1):
                rec("E", f"p={p} sheet {sheet}", abs(meas[(p, sheet)][0] - pred[(p, sheet)][0]), tol)
        # F
        for p in CONTROLS:
            for sheet in (0, 1):
                rec("F", f"p={p} sheet {sheet}", abs(meas[(p, sheet)][0] - pred[(p, sheet)][0]), 0.008 + nse * se(p, sheet))
        # G
        for p in PRIMES:
            if p <= 131 and p not in SLOW:
                for sheet in (0, 1):
                    rec("G", f"p={p} sheet {sheet}", abs(meas[(p, sheet)][2] - pred[(p, sheet)][2]), 0.012 + nse * se(p, sheet))
        # H
        for p in SIGN_HIGH:
            for sheet in (0, 1):
                E = meas[(p, sheet)][0]
                rec("H", f"p={p} sheet {sheet} (E >= 0.10)", max(0.10 - E, 0.0), 0.0)
        for p in SIGN_LOW:
            for sheet in (0, 1):
                E = meas[(p, sheet)][0]
                rec("H", f"p={p} sheet {sheet} (-0.05 <= E <= 0.02)", max(E - 0.02, -0.05 - E, 0.0), 0.0)
        print(f"== tests with {label} ==")
        total = 0
        for test in "ABCDEFGH":
            items = results[test]
            total += len(items)
            nfail = sum(1 for it in items if not it[3])
            # tightest item: largest deviation / tolerance (for H: smallest margin)
            if test != "H":
                tight = max(items, key=lambda it: it[1] / it[2])
                extra = f"tightest: {tight[0]}: deviation {tight[1]:.5f} of tolerance {tight[2]:.5f} ({100 * tight[1] / tight[2]:.0f} %)"
            else:
                margins = []
                for p in SIGN_HIGH:
                    for sheet in (0, 1):
                        margins.append((meas[(p, sheet)][0] - 0.10, f"p={p} sheet {sheet} E={meas[(p, sheet)][0]:+.4f}"))
                for p in SIGN_LOW:
                    for sheet in (0, 1):
                        E = meas[(p, sheet)][0]
                        margins.append((min(0.02 - E, E + 0.05), f"p={p} sheet {sheet} E={E:+.4f}"))
                mm = min(margins)
                extra = f"smallest margin: {mm[1]} (margin {mm[0]:.4f})"
            print(f"  {test}: {len(items)} items, {nfail} outside -> {'CONFIRMED' if nfail == 0 else 'REFUTED'};  {extra}")
            for it in items:
                if not it[3]:
                    print(f"       OUTSIDE: {it[0]}: deviation {it[1]:.5f} > tolerance {it[2]:.5f}")
        print(f"  total items: {total}")

    run_tests(se_reg, 4, "the registered rule (se = sqrt(0.75 p / #odd N), 4 se)")
    run_tests(se_meas, 4, "se from the measured relative variance (0.446 / 0.390), 4 se")
    run_tests(se_meas, 3, "se from the measured relative variance, 3 se")
    run_tests(se_meas, 2, "se from the measured relative variance, 2 se")
    print()
    print("== the Fermat / Wagstaff primes and what the D interval allows (registered rule) ==")
    for p in SPECIAL:
        for sheet in (0, 1):
            E, Ep = meas[(p, sheet)][0], pred[(p, sheet)][0]
            tol = 0.008 + 0.08 * abs(Ep) + 4 * se_reg(p, sheet)
            print(f"   p = {p:6d} sheet {sheet}: E measured {E:+.4f}, predicted {Ep:+.4f}; registered interval [{Ep - tol:+.4f}, {Ep + tol:+.4f}] ({'contains 0' if Ep - tol <= 0 else 'excludes 0'}); deviation in measured-variance se: {(E - Ep) / se_meas(p, sheet):+.2f}")
    print()
    print("== measured / predicted for the primes not in the theorem's list ==")
    for p in [47, 53, 59, 67, 71, 79, 83, 89, 97, 101, 103, 107, 109, 113, 131, 241, 683, 43691]:
        print(f"   p = {p:6d}: E measured {meas[(p, 0)][0]:+.5f} / predicted {pred[(p, 0)][0]:+.5f} (s = 0);  {meas[(p, 1)][0]:+.5f} / {pred[(p, 1)][0]:+.5f} (s >= 1);  E_nd measured {meas[(p, 0)][2]:+.5f}, {meas[(p, 1)][2]:+.5f}")


if __name__ == "__main__":
    main()
