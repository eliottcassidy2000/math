#!/usr/bin/env python3
"""h_high_gap_unlock_sampling_monad_s9.py — heavy targeted n=10/11/12 sampling
to unlock the 157 "permanent-through-n=9" HIGH gaps (monad-compute-2026-06-04-S9).

Context (OPEN-Q-055 / S7 handoff):
    The exhaustive H-spectra n=3..9 leave 159 odd values <= maxH9=3357 that are
    NEVER achieved through n=9.  Of these, the LOW set is exactly {7,21} (the
    permanent gaps); the other 157 are HIGH-end sparseness in the window
    [2883, 3355] — odd values just below the n=9 maximum H=3357 that simply did
    not occur among the 191,536 iso classes at n=9.  S7's light n=10/11 sampling
    confirmed only 9/157 as achieved (transient).  A heavier, TARGETED sweep can
    convert many more from "unknown" to "achieved (transient, NOT permanent)".

Method:
    These 157 targets sit near the n=9 MAX, i.e. on the LOW side of the n>=10
    distribution.  Uniform-random tournaments at n=10 concentrate at much larger
    H, so we BIAS the sampler toward near-transitive tournaments: build T as a
    transitive order with each forward arc independently REVERSED (an "upset")
    with probability p, then sweep p over a fine grid so the achieved-H cloud
    sweeps across the [2883, 3355] window.  H is computed with the validated
    Held-Karp engine H_count() from h21_finite_check_v2_monad_s4.

    For every sampled tournament we record H in a global achieved-set and note
    the SMALLEST n at which each target value is hit.  Sampling certifies
    achievability (a concrete witness) — it is a rigorous lower bound that a
    value is TRANSIENT, never a proof of permanence.

Outputs: per-n target-hit counts, the cumulative unlocked target set, and a tsv
    (target H, first n that achieves it) for the newly-unlocked values.

Session: monad-compute-2026-06-04-S9.
"""
import sys, os, time, random

sys.stdout.reconfigure(line_buffering=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h21_finite_check_v2_monad_s4 import H_count

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "05-knowledge", "results")

SEED = 20260604
WALL_BUDGET = 1500.0          # seconds of sampling work (leave margin under timeout)
LEVELS = [10, 11, 12]         # n values to sample
WINDOW_LO, WINDOW_HI = 2883, 3355


def load_targets():
    """The 157 HIGH permanent-through-n=9 gaps, from the S7 unlock tsv."""
    path = os.path.join(RESULTS, "h_unlock_table_monad_s7.tsv")
    perm = []
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            h = int(parts[0])
            u = parts[1] if len(parts) > 1 else ""
            if u == "":
                perm.append(h)
    high = sorted(g for g in perm if g > 609)   # exclude {7,21}
    return set(high)


def random_biased_tournament(m, p, rng):
    """Transitive base (i beats j for i<j) with each forward arc reversed w.p. p.
    Returns beats[] bitmask list (beats[v] = vertices v beats)."""
    beats = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            # default i->j; with prob p flip to j->i
            if rng.random() < p:
                beats[j] |= 1 << i
            else:
                beats[i] |= 1 << j
    return beats


def main():
    t0 = time.time()
    rng = random.Random(SEED)
    targets = load_targets()
    print("=" * 78)
    print("  HIGH-GAP UNLOCK SAMPLING  (monad-compute-2026-06-04-S9)")
    print("=" * 78)
    print(f"  targets: {len(targets)} HIGH permanent-through-n=9 gaps in "
          f"[{min(targets)}, {max(targets)}]")
    print(f"  levels: n={LEVELS}   seed={SEED}   wall budget={WALL_BUDGET:.0f}s")
    print(f"  bias window of interest: H in [{WINDOW_LO}, {WINDOW_HI}]")
    print()

    # already-known unlocked (S7 sampling found 9): we just re-derive fresh.
    first_n = {}                 # target H -> smallest n that achieves it
    achieved_all = set()         # every H value ever produced (any n)
    per_n_window_hits = {n: 0 for n in LEVELS}     # samples landing in window
    per_n_samples = {n: 0 for n in LEVELS}

    # Per-n upset-probability grids, CALIBRATED (smoke tests) so the achieved-H
    # cloud sweeps the [2883,3355] window: higher n needs fewer upsets (smaller p)
    # to land at the same H.  Window-hit rate peaks near these centers.
    P_GRID = {
        10: [round(0.16 + 0.004 * k, 4) for k in range(0, 30)],   # 0.160..0.276
        11: [round(0.10 + 0.003 * k, 4) for k in range(0, 30)],   # 0.100..0.187
        12: [round(0.07 + 0.0025 * k, 4) for k in range(0, 30)],  # 0.070..0.1425
    }
    # weight time toward faster, better-covered levels (n=10 fastest).
    WEIGHT = {10: 0.42, 11: 0.34, 12: 0.24}

    for n in LEVELS:
        n_start = time.time()
        deadline = n_start + WALL_BUDGET * WEIGHT[n]
        p_grid = P_GRID[n]
        gi = 0
        batch = 400
        last_report = n_start
        while time.time() < deadline:
            p = p_grid[gi % len(p_grid)]
            gi += 1
            for _ in range(batch):
                beats = random_biased_tournament(n, p, rng)
                h = H_count(beats, n)
                achieved_all.add(h)
                per_n_samples[n] += 1
                if WINDOW_LO <= h <= WINDOW_HI:
                    per_n_window_hits[n] += 1
                if h in targets and h not in first_n:
                    first_n[h] = n
            now = time.time()
            if now - last_report >= 20.0:
                hits = len(first_n)
                print(f"  [n={n}] t={now-t0:6.0f}s  samples={per_n_samples[n]:>9,}  "
                      f"window_hits={per_n_window_hits[n]:>7,}  "
                      f"targets_unlocked={hits}/{len(targets)}")
                last_report = now
        print(f"  [n={n}] DONE  samples={per_n_samples[n]:,}  "
              f"window_hits={per_n_window_hits[n]:,}  "
              f"cum targets unlocked={len(first_n)}/{len(targets)}  "
              f"({time.time()-n_start:.0f}s)")
        print()

    # ---- report --------------------------------------------------------------
    print("=" * 78)
    print("  RESULTS")
    print("=" * 78)
    unlocked = sorted(first_n)
    still_unknown = sorted(targets - set(first_n))
    print(f"  total HIGH targets:                {len(targets)}")
    print(f"  NEWLY CONFIRMED ACHIEVED (transient): {len(unlocked)}")
    print(f"  still unknown (no witness found):     {len(still_unknown)}")
    print()
    by_level = {n: sorted(h for h, fn in first_n.items() if fn == n) for n in LEVELS}
    for n in LEVELS:
        print(f"    first achieved at n={n}: {len(by_level[n])} values")
    print()
    print(f"  unlocked target values ({len(unlocked)}):")
    for i in range(0, len(unlocked), 12):
        print("    " + " ".join(f"{h}" for h in unlocked[i:i + 12]))
    print()
    if still_unknown:
        print(f"  STILL-UNKNOWN target values ({len(still_unknown)}):")
        for i in range(0, len(still_unknown), 12):
            print("    " + " ".join(f"{h}" for h in still_unknown[i:i + 12]))
    print()

    # sanity: confirm 7 and 21 NEVER appear (consistent with permanent gap)
    print(f"  sanity — H=7 ever sampled?  {7 in achieved_all}")
    print(f"  sanity — H=21 ever sampled? {21 in achieved_all}")
    print(f"  distinct H values produced overall: {len(achieved_all)} "
          f"(range [{min(achieved_all)}, {max(achieved_all)}])")

    # ---- save ----------------------------------------------------------------
    outpath = os.path.join(RESULTS, "h_high_gap_unlock_sampling_monad_s9.tsv")
    with open(outpath, "w") as f:
        f.write("target_H\tfirst_n_achieved\n")
        for h in sorted(targets):
            f.write(f"{h}\t{first_n.get(h, '')}\n")
    print(f"\n  saved: {outpath}")
    print(f"    (blank first_n = no witness found in this sampling run)")
    print(f"\n  DONE.  total elapsed {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
