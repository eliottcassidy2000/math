#!/usr/bin/env python3
"""h_unlock_table_monad_s7.py — H-UNLOCK TABLE n=3..9 (monad-compute-2026-06-04-S7).

Answers the explicit OPEN-Q-055 sub-question:
    "Is the forbidden set finite? At what n does each forbidden value 'unlock'?"

For every odd H value we compute  unlock(H) = the SMALLEST n in {3,...,9} at which
some n-vertex tournament has H(T)=H.  Built entirely from EXHAUSTIVE per-level
H-spectra (no sampling), so unlock(H) is exact for every H that unlocks by n=9:

  * n=3..7 : generated here via the validated isomorph-free orderly engine
             (h21_finite_check_v2_monad_s4, color-refinement canonical form;
             iso-class counts checked against A000568 = 2,4,12,56,456).
  * n=8    : parsed from h_spectrum_n8_exhaustive_monad.out  (full 2^28 census, S1).
  * n=9    : parsed from h_spectrum_n9_histogram_monad_s6.tsv (191,536 iso classes, S6).

Outputs the unlock CASCADE (how many NEW odd H values first appear at each n),
the explicit unlock-n for every transient gap (a value that is a gap at some
level but is achieved later), and confirms that {7,21} NEVER unlock through n=9
(the permanent-gap set).  Cross-references the SAMPLED n=10/n=11 spectra to mark
which n=9 high-end gaps are additionally seen to unlock by n<=11 (achievability
lower bound only — sampling cannot certify a gap).

Session: monad-compute-2026-06-04-S7.
"""
import sys, os, re, time
sys.stdout.reconfigure(line_buffering=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from h21_finite_check_v2_monad_s4 import (
    extend, beats_from_canon, H_count,
)

A000568 = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880, 9: 191536}
NOCAP = 10 ** 9  # no pruning -> all iso classes
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "05-knowledge", "results")


def exhaustive_Hset(maxn):
    """Return {n: set(achieved H)} for n=3..maxn via isomorph-free enumeration."""
    out = {}
    R = {1: {0}}
    for k in range(1, maxn):
        nxt = set()
        for code in R[k]:
            extend(beats_from_canon(code, k), k, NOCAP, nxt)
        R[k + 1] = nxt
        del R[k]
        n = k + 1
        exp = A000568.get(n)
        if exp is not None and exp != len(nxt):
            print(f"  ABORT: iso-class mismatch at n={n}: {len(nxt)} != {exp}")
            sys.exit(1)
        if n >= 3:
            hs = set()
            for code in nxt:
                hs.add(H_count(beats_from_canon(code, n), n))
            out[n] = hs
            print(f"  n={n}: {len(nxt):,} iso classes (A000568 OK), "
                  f"{len(hs)} distinct H, range [{min(hs)},{max(hs)}]")
    return out


def parse_n8():
    """Achieved-H set at n=8 from the exhaustive monad .out histogram."""
    path = os.path.join(RESULTS, "h_spectrum_n8_exhaustive_monad.out")
    hs = set()
    pat = re.compile(r"^\s*H=\s*(\d+):\s*(\d+)\s*$")
    with open(path) as f:
        for line in f:
            mo = pat.match(line)
            if mo:
                h, cnt = int(mo.group(1)), int(mo.group(2))
                if cnt > 0:
                    hs.add(h)
    print(f"  n=8: parsed {len(hs)} distinct H from exhaustive census, "
          f"range [{min(hs)},{max(hs)}]")
    return hs


def parse_n9():
    """Achieved-H set at n=9 from the exhaustive iso-class histogram tsv."""
    path = os.path.join(RESULTS, "h_spectrum_n9_histogram_monad_s6.tsv")
    hs = set()
    with open(path) as f:
        next(f)  # header
        for line in f:
            parts = line.split("\t")
            if len(parts) >= 2:
                h, c = int(parts[0]), int(parts[1])
                if c > 0:
                    hs.add(h)
    print(f"  n=9: parsed {len(hs)} distinct H from exhaustive iso-classes, "
          f"range [{min(hs)},{max(hs)}]")
    return hs


def parse_sampled(fname):
    """Achieved-H set from a SAMPLED spectrum .out (lower bound on achievability)."""
    path = os.path.join(RESULTS, fname)
    hs = set()
    pat = re.compile(r"^\s*(\d+)\s*\|\s*(\d+)")
    with open(path) as f:
        for line in f:
            mo = pat.match(line)
            if mo:
                h, c = int(mo.group(1)), int(mo.group(2))
                if c > 0:
                    hs.add(h)
    return hs


def main():
    print("=" * 78)
    print("  H-UNLOCK TABLE  n=3..9  (exhaustive) + n=10,11 sampled cross-check")
    print("=" * 78)
    print()
    t0 = time.time()

    print("[1] exhaustive H-sets n=3..7 (orderly engine):")
    sets = exhaustive_Hset(7)
    print()
    print("[2] exhaustive H-sets n=8,9 (parsed from saved census):")
    sets[8] = parse_n8()
    sets[9] = parse_n9()
    print(f"\n  (generation+parse done in {time.time()-t0:.1f}s)")

    # ---- unlock(H) = smallest n in 3..9 achieving H ----
    unlock = {}
    for n in range(3, 10):
        for h in sets[n]:
            if h not in unlock:
                unlock[h] = n
    maxH9 = max(sets[9])
    achieved_by9 = set(unlock)

    # ---- cascade: new values introduced at each level ----
    print("\n" + "=" * 78)
    print("  UNLOCK CASCADE — NEW odd H values first achieved at each n")
    print("=" * 78)
    print(f"  {'n':>3} {'distinct H':>11} {'maxH':>7} {'NEW @ n':>9} {'cum distinct':>13}")
    cum = set()
    for n in range(3, 10):
        new = sets[n] - cum
        cum |= sets[n]
        print(f"  {n:>3} {len(sets[n]):>11} {max(sets[n]):>7} {len(new):>9} {len(cum):>13}")

    # ---- gap analysis per level + which gaps later unlock ----
    print("\n" + "=" * 78)
    print("  GAP -> UNLOCK ANALYSIS")
    print("=" * 78)
    # A value H is a 'gap at n' if H is odd, H <= maxH(n), and H not in sets[n].
    permanent = []          # gap at every level 3..9
    unlock_n = {}           # H -> first n achieving it, for H that is a gap at some earlier level
    for h in range(1, maxH9 + 1, 2):
        # levels where h is within range but not achieved
        gap_levels = [n for n in range(3, 10)
                      if h <= max(sets[n]) and h not in sets[n]]
        if not gap_levels:
            continue
        if h in achieved_by9:
            unlock_n[h] = unlock[h]
        else:
            # never achieved through n=9, and within n=9 range -> candidate permanent
            permanent.append(h)

    print(f"\n  PERMANENT-through-n=9 gaps (odd, <= n=9 maxH={maxH9}, never achieved):")
    print(f"    {permanent[:40]}{' ...' if len(permanent) > 40 else ''}")
    print(f"    count = {len(permanent)}")
    low_perm = [g for g in permanent if g <= 609]
    print(f"    of these, LOW (<=609): {low_perm}")
    print(f"    => low permanent-through-9 gap set = {sorted(low_perm)}")

    print(f"\n  TRANSIENT gaps that UNLOCK at a later n (gap earlier, achieved by n=9):")
    print(f"    {'H':>6} {'unlock_n':>9}   (was a gap at the levels below unlock_n)")
    trans = sorted(unlock_n)
    # show all transient unlocks with unlock_n >= 4
    shown = 0
    for h in trans:
        un = unlock_n[h]
        # was it a gap before un? yes by construction of gap_levels
        print(f"    {h:>6} {un:>9}")
        shown += 1
        if shown >= 60:
            print(f"    ... ({len(trans)-shown} more)")
            break
    print(f"    total transient-unlock values: {len(trans)}")

    # ---- the famous values ----
    print("\n" + "=" * 78)
    print("  KEY VALUES")
    print("=" * 78)
    for h in (7, 21, 35, 39, 49, 63):
        u = unlock.get(h)
        status = f"unlocks at n={u}" if u else "NEVER unlocks through n=9 (PERMANENT gap)"
        print(f"    H={h:>4}: {status}")

    # ---- sampled n=10/11 cross-check on n=9 high gaps ----
    print("\n" + "=" * 78)
    print("  SAMPLED n=10 / n=11 CROSS-CHECK (achievability lower bound only)")
    print("=" * 78)
    s10 = parse_sampled("h_spectrum_n10.out")
    s11 = parse_sampled("h_spectrum_n11.out")
    print(f"  sampled n=10: {len(s10)} distinct H; sampled n=11: {len(s11)} distinct H")
    # the low permanent gaps: are they ever seen in n=10/11 samples? (should NOT be)
    for h in (7, 21):
        seen = []
        if h in s10: seen.append(10)
        if h in s11: seen.append(11)
        print(f"    H={h}: seen in sampled levels {seen if seen else '[] (absent — consistent with permanent gap)'}")
    # of the n=9 permanent-through-9 LOW gaps beyond {7,21} (should be none): report
    extra_low = [g for g in low_perm if g not in (7, 21)]
    print(f"    low permanent-through-9 gaps beyond {{7,21}}: {extra_low}")
    # how many n=9 high gaps are seen to unlock by n<=11 sampling
    n9_high_gaps = [g for g in permanent if g > 609]
    unlocked_by_sample = sorted(g for g in n9_high_gaps if g in s10 or g in s11)
    print(f"    n=9 (>609) gaps seen achieved in n=10/11 samples: "
          f"{len(unlocked_by_sample)}/{len(n9_high_gaps)}")
    if unlocked_by_sample:
        print(f"      e.g. {unlocked_by_sample[:20]}")

    # ---- save the unlock table ----
    outpath = os.path.join(RESULTS, "h_unlock_table_monad_s7.tsv")
    with open(outpath, "w") as f:
        f.write("H\tunlock_n\n")
        for h in range(1, maxH9 + 1, 2):
            u = unlock.get(h, "")  # blank = not achieved through n=9
            f.write(f"{h}\t{u}\n")
    print(f"\n  unlock table saved: {outpath}")
    print(f"    (column unlock_n blank = odd H <= {maxH9} not achieved through n=9)")
    print("\n  DONE.  total elapsed {:.1f}s".format(time.time() - t0))


if __name__ == "__main__":
    main()
