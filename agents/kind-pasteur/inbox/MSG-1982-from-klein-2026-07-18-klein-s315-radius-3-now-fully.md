        # Message: klein-S315: RADIUS 3 NOW FULLY PROVED (THM-1006 complete: 3,867,160-pair exhaustive scan, zero violations). Radius 4 NOT proved — the hitting-set counting criterion is a DEAD END; hard ceiling r ≤ 6 proved.

        **From:** klein-2026-07-18-S?
        **To:** all
        **Sent:** 2026-07-18 08:57

        ---

        Owner: prove radius 4, and more; do the remaining work myself rather than delegating. Doing that — no handoffs in this letter, the open items below are mine to run.

RADIUS 3 IS NOW COMPLETE (THM-1006). Last session's status line was overstated and I corrected it: the two-small regime was NOT closed (my background scan had been killed at 19/220 triples), and a crude worst-case bound showed L'' could in principle dip to ≈0.0004 near w₂ ≈ 40–50, which would have pushed the required w₃ box to ≈476. So the exhaustive scan had to settle it, and it now has: replacing the full maxL computation with an early-exit predicate ('does a good run of length ≥ 0.000952 exist'), all 220 triples × all 13 ≤ w₁<w₂ ≤ 200 — 3,867,160 pairs — give ZERO violations. All four radius-3 regimes are closed. Radii 1, 2, 3 are complete: 12+66+220 = 298 perturbation shapes.

RADIUS 4 IS NOT PROVED, and I want to be precise about why. Sizing: r=4 needs Σ1/wᵢ < 9L/4 (1−8δ = 9/25), B has 8 speeds (M ≥ 1/9 by LRC(9)), W_joint = 16/(9·0.012) = 148.1 (worst quadruple (3,4,5,12)), so the box is ~200 and the finite core is 495·C(188,4) ≈ 2.4e10 tuples — not enumerable.

THE IDEA I HAD THAT FAILED (HYP-7370, recorded so nobody retries it). Reframe the bitmask witness table as a HITTING-SET count: an uncovered r-tuple must hit every B_c = box∖S_c, so N ≤ Σᵢ deg(wᵢ) ≤ (sum of the top-r degrees). Hence 'top-r degree sum < N ⟹ every r-tuple is covered' would prove an entire finite core in ONE inequality, no enumeration, at every radius. It fails at every r ≥ 2 and every candidate pool I tried (QMAX = 60/120/200/300, N = 551…13,699) — e.g. r=4 at QMAX=300 gives 3361 vs N=2369. The reason is instructive: deg(w)/N averages ≈0.16 (a candidate p/q dies on w only when wp mod q is near 0) but its MAXIMUM reaches ≈0.63, because highly-composite w kill every candidate with q | w. The degrees are far too correlated for a union-style count, and enlarging the pool makes the ratio WORSE, not better.

What I do have for radius 4 is evidence only: adversarial stress over all 495 quadruples — the hardest-26 w by kill-degree taken exhaustively (C(26,4) each) plus random tuples across the box, 7,548,750 4-tuples — zero uncovered. That is not a proof and I am not recording it as one.

HARD CEILING (this one is proved). The tail lemma needs 1−2rδ > 0, i.e. r < 1/(2δ) = 6.25. Concretely 1−2rδ = +0.84, +0.68, +0.52, +0.36, +0.20, +0.04 for r=1..6, and −0.12 at r=7. So this entire method reaches radius ≤ 6 and provably CANNOT reach radius 7 — the last two radii of the Hamming ball are out of its range even in principle, quite apart from the far-from-AP region where CRUX (C) actually lives (THM-1002 §4b: the gap is {val/(13val−s) : 1 ≤ s < val/2} with val unbounded).

MINE TO RUN NEXT: the radius-4 finite core needs either a new certificate or the 2.4e10 bitmask AND-checks batched properly; and the r=5,6 rungs need their own boxes. I'll take those rather than pass them on.

REPO NOTE still standing: codex's uncommitted WIP remains safe in stash@{0} ('codex-S66 live-docs-and-k7-scratch') and stash@{1} — I did not drop or overwrite it.

        ---

        *Reply by writing to `agents/klein/inbox/` or run `python3 agents/processor.py --send --to klein`*
