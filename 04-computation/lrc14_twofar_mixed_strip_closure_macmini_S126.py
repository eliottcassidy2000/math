#!/usr/bin/env python3
"""Two-far mixed-strip sweep — the LAST finite piece of the all-heights
two-far gap closure (mac-mini-2026-07-19-S126, HYP-8005 consequence).

CLOSURE LOGIC (HYP-8005 + HYP-7985):
  - both fars ≥ 265                      → M ≥ 3/41 (k=2 plateau lemma);
  - y ≥ 23x (any x ≥ 14)                 → M ≥ 3/41 (k=1 lemma on S∪{x},
                                            12 speeds, settled M ≥ 1/13);
  - x < y ≤ 300                          → HYP-7985 sweep (complete);
  - REMAINING: 15 ≤ x ≤ 264, 301 ≤ y < 23x   ← THIS SCRIPT (finite).
If this strip is empty of gap families, the ENTIRE two-far stratum
(S ⊂ {1..13}, |S| = 11, any two fars, ANY height) has no M ∈ (1/14, 3/41).

Filters: proved-necessary duty/no-14-mult/rung-pair conditions (HYP-7985).
Scanner: early-exit grid q ≤ QSCAN = 200 (exit ⟹ M ≥ 3/41, sound at any
height); non-exiting families go to the exact breakpoint-moduli verifier.
"""

from fractions import Fraction as F
from itertools import combinations
from multiprocessing import Pool
import time

ONE14 = F(1, 14)
GAP_HI = F(3, 41)
QSCAN = 200
XLO, XHI = 15, 264
YLO = 301
YFACT = 23


def rung_sums(smax):
    out = set()
    for s in range(15, smax + 1):
        lo = s // 14 + 1
        if lo * 41 < 3 * s:
            out.add(s)
    return out


def endangered(removed):
    d = {q for q in removed if q >= 7}
    r = set(removed)
    if {5, 10} <= r:
        d.add(5)
    if {6, 12} <= r:
        d.add(6)
    return d


def core_tables(S, qmax):
    core = [None, None]
    for q in range(2, qmax + 1):
        row = [q] * q
        for u in S:
            um = u % q
            for k in range(q):
                r = (um * k) % q
                d = q - r if q - r < r else r
                if d < row[k]:
                    row[k] = d
        core.append(row)
    return core


def exact_M_bp(W):
    mods = set()
    for u, v in combinations(W, 2):
        mods.add(u + v)
        mods.add(v - u)
    for w in W:
        mods.add(2 * w)
    mods.discard(0)
    best = F(0)
    best_t = None
    for q in mods:
        for k in range(1, q):
            m = min(min((w * k) % q, q - (w * k) % q) for w in W)
            if F(m, q) > best:
                best, best_t = F(m, q), F(k, q)
    return best, best_t


def work(removed_tup):
    removed = set(removed_tup)
    S = [u for u in range(1, 14) if u not in removed]
    duties = sorted(endangered(removed))
    ymax = YFACT * XHI
    RS = rung_sums(13 + XHI + ymax)
    core = core_tables(S, QSCAN)
    stats = [0, 0, 0, 0]
    survivors = []
    resistant_q = 0
    for x in range(XLO, XHI + 1):
        if x % 14 == 0:
            continue
        for y in range(max(YLO, x + 1), YFACT * x):
            if y % 14 == 0:
                continue
            stats[0] += 1
            ok = all((x % q == 0) or (y % q == 0) for q in duties)
            if not ok:
                continue
            stats[1] += 1
            if not (any(u + x in RS or u + y in RS for u in S)
                    or (x + y) in RS):
                continue
            stats[2] += 1
            exited = False
            for q in range(2, QSCAN + 1):
                row = core[q]
                xm, ym = x % q, y % q
                for k in range(1, q):
                    d = row[k]
                    r = (xm * k) % q
                    dx = q - r if q - r < r else r
                    if dx < d:
                        d = dx
                    r = (ym * k) % q
                    dy = q - r if q - r < r else r
                    if dy < d:
                        d = dy
                    if 41 * d >= 3 * q:
                        exited = True
                        if q > resistant_q:
                            resistant_q = q
                        break
                if exited:
                    break
            stats[3] += 1
            if not exited:
                survivors.append(sorted(set(S) | {x, y}))
    return removed_tup, stats, resistant_q, survivors


def main():
    comps = list(combinations(range(1, 14), 2))
    t0 = time.time()
    agg = [0, 0, 0, 0]
    maxq = 0
    surv = []
    with Pool(8) as pool:
        done = 0
        for removed_tup, stats, rq, survivors in \
                pool.imap_unordered(work, comps, chunksize=1):
            for i in range(4):
                agg[i] += stats[i]
            maxq = max(maxq, rq)
            surv += [(removed_tup, W) for W in survivors]
            done += 1
            with open("05-knowledge/results/lrc14_strip_journal_S126.txt",
                      "a") as jf:
                jf.write(f"{removed_tup} stats={stats} rq={rq} "
                         f"surv={survivors}\n")
            print(f"  [{done}/78] removed {removed_tup}: "
                  f"raw {stats[0]}, scanned {stats[3]}, "
                  f"survivors {len(survivors)}, {time.time()-t0:5.1f}s",
                  flush=True)
    print(f"\nTOTALS: raw {agg[0]}, duty {agg[1]}, rung {agg[2]}, "
          f"scanned {agg[3]}, max exit-q {maxq}, "
          f"survivors {len(surv)}, {time.time()-t0:.1f}s")
    ingap = []
    for removed_tup, W in surv:
        M, t = exact_M_bp(W)
        tag = ("BELOW?!" if M < ONE14 else "TIGHT?!" if M == ONE14
               else "IN-GAP !!!" if M < GAP_HI else "ok-above")
        print(f"  survivor {W}: M = {M} at {t}  => {tag}")
        if ONE14 < M < GAP_HI:
            ingap.append((W, M))
    print("\n==== STRIP VERDICT ====")
    if not ingap:
        print("The mixed strip is EMPTY of gap families.  Together with")
        print("HYP-7985 (x<y ≤ 300) and HYP-8005 (both ≥ 265; y ≥ 23x):")
        print("THE TWO-FAR STRATUM HAS NO FAMILY WITH M ∈ (1/14, 3/41)")
        print("AT ANY HEIGHT (modulo the plateau lemma's independent read).")
    else:
        for W, M in ingap:
            print(f"  IN-GAP: {W}: M = {M}")


if __name__ == "__main__":
    main()
