#!/usr/bin/env python3
"""HYP-7990 referee — the three-far stratum sweep (mac-mini-2026-07-19-S126).

STRATUM: W = S ∪ {x,y,z}, S ⊂ {1..13}, |S| = 10, 15 ≤ x < y < z ≤ VMAX
(fars exclude multiples of 14).  Extends HYP-7985 (two-far, gap-empty to
fars ≤ 300, small-witness ceiling q ≤ 41) to |removed| = 3 — the stratum
containing GW-Hamming-2 with 24 retained, and the duty-TRADING corner.

NECESSARY FILTERS for M ∈ (1/14, 3/41) — same proofs as HYP-7985:
 (N1) duties: endangered(removed) = {q ∈ removed : q ≥ 7}
      ∪ {5 if {5,10} ⊆ removed} ∪ {6 if {6,12} ⊆ removed}
      ∪ {4 if {4,8,12} ⊆ removed}   [new triple-removal case]
      — every endangered duty needs a multiple among the fars.
 (N2) no 14-multiples (covering ⟹ M ≥ 14/183 > 3/41).
 (N3) some pair sum ∈ the in-gap rung set (integer D with 41D/3 < s < 14D).
LRC(14) itself and the tight census on this stratum are already canon:
THM-738 (every 13-speed family with ≥ 10 speeds in {1..14} satisfies
LRC(14); tights = AP, GW only) — the GAP is the only open question, so the
rung filter is sound for exactly the question asked.

SCANNER: identical soundness to HYP-7985 — candidate grid t = k/q over all
integers q ≤ 13 + 2·VMAX (superset of all breakpoint moduli), ascending q,
early exit at clearance ≥ 3/41; families finishing the grid get exact M and
an independent recheck in the parent process.

Parallelism: one task per complement (286), multiprocessing.Pool.
"""

from fractions import Fraction as F
from itertools import combinations
from multiprocessing import Pool
import time

ONE14 = F(1, 14)
GAP_HI = F(3, 41)

VMAX_MAIN = 150
VMAX_DEEP = 400          # duty-trading complements only


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
    if {4, 8, 12} <= r:
        d.add(4)
    return d


def allowed_fars(vmax):
    return [v for v in range(15, vmax + 1) if v % 14 != 0]


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


def work(args):
    removed_tup, vmax = args
    removed = set(removed_tup)
    S = [u for u in range(1, 14) if u not in removed]
    duties = sorted(endangered(removed))
    qmax = 13 + 2 * vmax
    RS = rung_sums(13 + 2 * vmax)
    fars = allowed_fars(vmax)
    nf = len(fars)
    core = core_tables(S, qmax)

    # duty masks per far + small-rung-partner flag per far
    dmask = []
    full = (1 << len(duties)) - 1
    for f in fars:
        m = 0
        for i, q in enumerate(duties):
            if f % q == 0:
                m |= 1 << i
        dmask.append(m)
    sp = [any((u + f) in RS for u in S) for f in fars]

    stats = [0, 0, 0, 0]   # raw, duty, rung, scanned
    resistant = []          # (exit_q, W)
    survivors = []          # W lists needing parent verification
    Sset = set(S)

    for i in range(nf):
        mi = dmask[i]
        for j in range(i + 1, nf):
            mij = mi | dmask[j]
            for k in range(j + 1, nf):
                stats[0] += 1
                if (mij | dmask[k]) != full:
                    continue
                stats[1] += 1
                x, y, z = fars[i], fars[j], fars[k]
                if not (sp[i] or sp[j] or sp[k] or (x + y) in RS
                        or (x + z) in RS or (y + z) in RS):
                    continue
                stats[2] += 1
                stats[3] += 1
                # inline early-exit scan
                exited = False
                best_n, best_d = 0, 1
                for q in range(2, qmax + 1):
                    row = core[q]
                    xm, ym, zm = x % q, y % q, z % q
                    for kk in range(1, q):
                        d = row[kk]
                        r = (xm * kk) % q
                        dx = q - r if q - r < r else r
                        if dx < d:
                            d = dx
                        r = (ym * kk) % q
                        dy = q - r if q - r < r else r
                        if dy < d:
                            d = dy
                        r = (zm * kk) % q
                        dz = q - r if q - r < r else r
                        if dz < d:
                            d = dz
                        if 41 * d >= 3 * q:
                            exited = True
                            if len(resistant) < 6 or q > resistant[-1][0]:
                                resistant.append((q, sorted(Sset | {x, y, z})))
                                resistant.sort(key=lambda t: -t[0])
                                del resistant[6:]
                            break
                        if d * best_d > best_n * q:
                            best_n, best_d = d, q
                    if exited:
                        break
                if not exited:
                    survivors.append(sorted(Sset | {x, y, z}))
    return removed_tup, duties, stats, resistant, survivors


def exact_M_independent(W):
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


def run_pass(tasks, label):
    t0 = time.time()
    agg = [0, 0, 0, 0]
    all_res = []
    all_surv = []
    with Pool(8) as pool:
        done = 0
        for removed_tup, duties, stats, resistant, survivors in \
                pool.imap_unordered(work, tasks, chunksize=1):
            for a in range(4):
                agg[a] += stats[a]
            all_res += resistant
            all_surv += [(removed_tup, W) for W in survivors]
            done += 1
            with open("05-knowledge/results/lrc14_threefar_journal_S126.txt",
                      "a") as jf:
                jf.write(f"{label} {removed_tup} stats={stats} "
                         f"nres={len(resistant)} nsurv={len(survivors)} "
                         f"surv={survivors}\n")
            if done % 40 == 0:
                print(f"    [{label}] {done}/{len(tasks)} complements, "
                      f"{time.time()-t0:6.1f}s, scanned {agg[3]}", flush=True)
    all_res.sort(key=lambda t: -t[0])
    print(f"  PASS {label}: raw {agg[0]}, duty-pass {agg[1]}, ", flush=True) if False else print(f"  PASS {label}: raw {agg[0]}, duty-pass {agg[1]}, "
          f"rung-pass {agg[2]}, scanned {agg[3]}, "
          f"survivors {len(all_surv)}, {time.time()-t0:.1f}s")
    print("  most exit-resistant:")
    for q, W in all_res[:6]:
        print(f"    exit q={q:>3}: {W}")
    return all_surv


def main():
    # controls (parent process, independent verifier)
    for name, W, expect in (("AP", list(range(1, 14)), F(1, 14)),
                            ("GW", list(range(1, 12)) + [13, 24], F(1, 14)),
                            ("F_3(13)", list(range(1, 12)) + [13, 36], F(3, 41))):
        Mv, tv = exact_M_independent(W)
        assert Mv == expect, (name, Mv)
        print(f"  control {name}: exact M = {Mv}  OK")

    comps = list(combinations(range(1, 14), 3))
    print(f"\n==== MAIN PASS: all {len(comps)} complements, fars ≤ {VMAX_MAIN} ====")
    surv1 = run_pass([(c, VMAX_MAIN) for c in comps], f"main-{VMAX_MAIN}")

    deep = [c for c in comps if len(endangered(set(c))) >= 3
            or {10, 12} <= set(c) or {5, 10} <= set(c) or {6, 12} <= set(c)]
    print(f"\n==== DEEP PASS: {len(deep)} duty-trading complements "
          f"(≥3 duties, or ⊇{{10,12}}/{{5,10}}/{{6,12}}), fars ≤ {VMAX_DEEP} ====")
    print("  (general 2-duty complements are complete at the main-pass window;")
    print("   scope note: their 200<fars≤600 range is not swept — recorded.)")
    surv2 = run_pass([(c, VMAX_DEEP) for c in deep], f"deep-{VMAX_DEEP}")

    print("\n==== SURVIVOR VERIFICATION ====")
    ingap = []
    for removed_tup, W in surv1 + surv2:
        Mv, tv = exact_M_independent(W)
        tag = ("BELOW FLOOR?!" if Mv < ONE14 else
               "TIGHT?!" if Mv == ONE14 else
               "IN-GAP !!!" if Mv < GAP_HI else "ok-above (slow exit)")
        print(f"  {W} (removed {removed_tup}): M = {Mv} at t={tv}  => {tag}")
        if ONE14 <= Mv < GAP_HI and Mv != ONE14:
            ingap.append((W, Mv))
    print("\n==== VERDICT ====")
    if not ingap and not (surv1 or surv2):
        print(f"NO three-far family (|S|=10, fars ≤ {VMAX_MAIN} complete; "
              f"duty-rich complements to {VMAX_DEEP}) has M in (1/14, 3/41).")
        print("No survivor even reached full scan. In particular no three-far")
        print(f"4/55 realizer with fars ≤ {VMAX_MAIN} (≤ {VMAX_DEEP} duty-rich).")
    elif not ingap:
        print("Full-scan survivors existed but all verified ≥ 3/41 or floor.")
    else:
        print("IN-GAP FAMILIES FOUND — verify in a second session before ANY")
        print("canon claim:")
        for W, Mv in ingap:
            print(f"  {W}: M = {Mv}")


if __name__ == "__main__":
    main()
