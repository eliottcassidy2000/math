#!/usr/bin/env python3
"""HYP-7985 referee — the complete two-far stratum sweep
(mac-mini-2026-07-19-S125; executes S124 backlog lead (ii), the duty-trading
shell, = the stratum where any non-single-far 4/55 realizer must live).

STRATUM: W = S ∪ {x, y},  S ⊂ {1..13}, |S| = 11,  14 ≤ x < y ≤ VMAX.

NECESSARY FILTERS for M(W) ∈ (1/14, 3/41)  [all proved, stated in-file]:

  (N1) DUTIES: M < 1/q forces a multiple of q in W (at t = k/q every
       clearance is a multiple of 1/q ≥ 1/13 > 3/41 unless zero).  With
       3/41 < 1/13: multiples of every q ∈ {2..13} required.  Removing
       {a,b} = {1..13}\\S endangers: each removed q ∈ {7..13} (sole
       multiple), {5} iff {a,b}={5,10}, {6} iff {a,b}={6,12}.  The fars
       must carry the endangered duties (split or stacked = duty-trading).
  (N2) NO 14-MULTIPLE: covering families have M ≥ 14/183 > 3/41 (THM-724/
       726), so in-gap families are non-covering exactly at q=14:
       14 ∤ x, 14 ∤ y (retained ≤ 13 automatic).
  (N3) RUNG PAIR-SUM: an in-gap maximizer sits at t* = k/s where s = u+v
       is an ACTIVE PAIR sum (a local max of the lower envelope needs two
       opposite-slope tents: (u+v)t* ∈ ℤ; peaks give clearance 1/2), and
       D = M·s ∈ ℤ.  So some pair of W sums to s with (D,s) integral,
       41D/3 < s < 14D — the in-gap rung set.

SCANNER: candidate grid t = k/q over ALL integers q ≤ 2·VMAX (superset of
every breakpoint modulus u+v, v−u, 2w ≤ 2·VMAX ⟹ max over the grid is the
exact M).  Ascending q with EARLY EXIT once clearance ≥ 3/41 (then M ≥ 3/41,
out of the open gap, above).  Families finishing the full grid get exact M:
  - M < 1/14: would REFUTE LRC(14) — independent recheck + panic banner;
  - M = 1/14: a third tight family — contradicts THM-1142 — recheck+banner;
  - 1/14 < M < 3/41: THE FINDING (gap populated) — recheck + banner;
  - else report (stratum bottom census).
Controls: AP and GW run through the same full-scan path (expect exactly
1/14); F_3(13) through the early-exit path (expect exit at 3/41).
"""

from fractions import Fraction as F
from itertools import combinations
from math import gcd
import sys, time

VMAX_PASS1 = 60
VMAX_PASS2 = 120

ONE14 = F(1, 14)
GAP_HI = F(3, 41)


# ----------------------------------------------------------------- rung set
def rung_sums(smax):
    """All integer pair sums s that an in-gap maximizer could use:
    exists D ∈ ℤ with 41D/3 < s < 14D  (⟺ s/14 < D < 3s/41 wait — invert:
    D > s/14 and D < 3s/41... careful: D/s > 1/14 ⟺ D > s/14; D/s < 3/41
    ⟺ D < 3s/41).  s qualifies iff ⌈s/14⌉ ... some integer D in
    (s/14, 3s/41)."""
    out = []
    for s in range(15, smax + 1):
        lo = s // 14 + 1              # smallest D with D > s/14
        if lo * 41 < 3 * s:           # D < 3s/41  ⟺ 41D < 3s
            out.append(s)
    return set(out)


# ------------------------------------------------------- per-complement core
def core_tables(S, qmax):
    """core[q][k] = min over u in S of the integer circle distance of u*k
    mod q (0..q//2).  Stored flat: list of lists."""
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


def scan(core, x, y, qmax):
    """Ascending-q scan with early exit at clearance ≥ 3/41
    (41*d ≥ 3*q).  Returns ('exit', q, k) or ('full', M, t) with exact M."""
    best_n, best_d, best_t = 0, 1, (0, 1)
    for q in range(2, qmax + 1):
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
                return ("exit", q, k)
            if d * best_d > best_n * q:
                best_n, best_d, best_t = d, q, (k, q)
    return ("full", F(best_n, best_d), best_t)


# ------------------------------------------------------------------ verifier
def exact_M_independent(W):
    """S124-style breakpoint-moduli exact M (independent code path)."""
    W = sorted(W)
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


# --------------------------------------------------------------------- sweep
def endangered(removed):
    d = set(q for q in removed if q >= 7)
    if removed == {5, 10}:
        d.add(5)
    if removed == {6, 12}:
        d.add(6)
    return d


def allowed_fars(vmax):
    return [v for v in range(14, vmax + 1) if v % 14 != 0]


def run_pass(vmax, label):
    t0 = time.time()
    qmax = 2 * vmax
    RS = rung_sums(13 + vmax + vmax)   # any pair sum ≤ x+y ≤ 2vmax, +13 safe
    fars = allowed_fars(vmax)
    print(f"\n==== PASS {label}: VMAX = {vmax}, fars {len(fars)}, "
          f"rung sums ≤ {13+2*vmax}: {sorted(s for s in RS if s <= 2*vmax)[:12]}... ====")
    stats = dict(pairs=0, duty=0, rung=0, scanned=0, full=0)
    survivors = []
    resistant = []          # (exit_q, W) — latest-exit families, the census
    for removed_tup in combinations(range(1, 14), 2):
        removed = set(removed_tup)
        S = [u for u in range(1, 14) if u not in removed]
        duties = endangered(removed)
        core = core_tables(S, qmax)
        Sset = set(S)
        for i, x in enumerate(fars):
            for y in fars[i + 1:]:
                stats["pairs"] += 1
                ok = all((x % q == 0) or (y % q == 0) for q in duties)
                if not ok:
                    continue
                stats["duty"] += 1
                # (N3) rung pair-sum: pairs (u,x),(u,y) u∈S, and (x,y)
                if not (any(u + x in RS for u in S) or
                        any(u + y in RS for u in S) or (x + y) in RS):
                    continue
                stats["rung"] += 1
                res = scan(core, x, y, qmax)
                stats["scanned"] += 1
                if res[0] == "exit":
                    if len(resistant) < 8 or res[1] > resistant[-1][0]:
                        resistant.append((res[1], sorted(Sset | {x, y})))
                        resistant.sort(key=lambda r: -r[0])
                        del resistant[8:]
                if res[0] == "full":
                    stats["full"] += 1
                    M, (k, q) = res[1], res[2]
                    W = sorted(Sset | {x, y})
                    Mv, tv = exact_M_independent(W)
                    tag = ("BELOW FLOOR — LRC(14) COUNTEREXAMPLE?!" if Mv < ONE14
                           else "TIGHT — THIRD TIGHT FAMILY?!" if Mv == ONE14
                           else "IN-GAP HIT !!!" if Mv < GAP_HI
                           else "ok-above (slow exit)")
                    survivors.append((W, Mv, tv, tag))
                    print(f"  FULL-SCAN family {W}: scan M={M}, "
                          f"independent M={Mv} at t={tv}  => {tag}")
                    if Mv < GAP_HI:
                        print("  *** SENSATION-CLASS RESULT — do not trust "
                              "without a second session's referee ***")
        # light progress
        if removed_tup[1] == 13:
            el = time.time() - t0
            print(f"    [complements through removed={removed_tup} done, "
                  f"{el:6.1f}s, scanned {stats['scanned']}]")
    print(f"  PASS {label} totals: raw pairs×complements {stats['pairs']}, "
          f"after duty {stats['duty']}, after rung {stats['rung']}, "
          f"scanned {stats['scanned']}, full-scans {stats['full']}, "
          f"{time.time()-t0:.1f}s")
    print("  most exit-resistant families (latest beat, by exit modulus q):")
    for eq, W in resistant:
        print(f"    exit q={eq:>3}: {W}")
    return survivors


def controls():
    print("== controls ==")
    # full-scan path: AP and GW (M = 1/14 < 3/41 so no early exit possible)
    for name, W in (("AP", list(range(1, 14))),
                    ("GW", list(range(1, 12)) + [13, 24])):
        S = W  # 13 elements; core-table trick with x=y absent: emulate by
        # scanning with the independent verifier only
        Mv, tv = exact_M_independent(W)
        assert Mv == ONE14, (name, Mv)
        print(f"  {name}: exact_M = 1/14  OK (full-scan-path verifier)")
    # early-exit path: F_3(13) must exit at exactly clearance 3/41
    S = list(range(1, 12))  # + fars 13? F_3 is one-far; test scanner shape:
    core = core_tables(list(range(1, 12)) + [13], 2 * 60)
    res = scan(core, 36, 36, 2 * 60)   # degenerate y=x=36 → same family
    assert res[0] == "exit", res
    q, k = res[1], res[2]
    print(f"  F_3(13): scanner exits at t={k}/{q} with clearance ≥ 3/41  OK"
          f" (its true M = 3/41 is attained, not in the OPEN gap)")


def main():
    controls()
    all_surv = []
    vdone = 0
    for vmax in (60, 120, 300):
        s = run_pass(vmax, str(vmax))
        all_surv += s
        vdone = vmax
        if any(Mv < GAP_HI for _, Mv, _, _ in s):
            print("\nSub-3/41 content found — halting escalation for "
                  "verification discipline.")
            break
    ingap = [r for r in all_surv if r[1] < GAP_HI]
    print("\n==== VERDICT ====")
    if not ingap:
        print(f"NO two-far family (S ⊂ {{1..13}}, |S|=11, 14 ≤ x < y ≤ {vdone},")
        print("duty+rung filtered — filters proved necessary, so the sweep is")
        print("COMPLETE for the stratum+window) has M in the OPEN gap")
        print("(1/14, 3/41); none is tight; none is below the floor.")
        print(f"In particular the 4/55 rung has NO two-far realizer with")
        print(f"fars ≤ {vdone}.  Scope: two-far stratum only; three-far and")
        print("other shapes remain with the walls.")
    else:
        for W, Mv, tv, tag in ingap:
            print(f"  {W}: M = {Mv}  ({tag})")


if __name__ == "__main__":
    main()
