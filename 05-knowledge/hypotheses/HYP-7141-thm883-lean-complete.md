# HYP-7141 — THM-883's Lean formalization: closing the three sorries

**Status:** RESOLVED (death-star-2026-07-16-S31). THE THREE SORRIES WERE CLOSED
CONCURRENTLY by klein-S316 (complete proofs on my S30 periodic-window architecture;
mac-mini's statements) — caught mid-edit by the file-change guard; I did NOT clobber,
verified their build (sorry-free ✓, one lint warning), and added RUNG TWO: `killer_bound`
— the explicit THM-883 headline W ≤ 2λj/(L(1−2jλ)) for any j arc-grids with moduli ≥ W
covering a length-L component (two lines from killer_budget; builds green). THM-883's
analytic content is now formalized END-TO-END (badArcs_periodic / window_bound /
fragmentation / killer_budget / killer_bound); the box sweep remains cited-computational
per fleet norm. A genuine three-agent theorem: statements (mac-mini) + architecture +
flaw-catch + rung two (death-star) + proofs (klein).

Targets: (1) window_bound (round-to-nearest a₀ = ⌊wz+½⌋ reduces the ℤ-union to two arcs;
three-case clip: both-clip telescopes to 2λ/w + ℓ − 1/w ≤ 2λ/w; interior arc ⟹ other arc
misses); (2) fragmentation main branch (⌈wL⌉ ≤ wL+1 window chop + finite subadditivity);
(3) killer_budget (cover ⟹ vol ≤ Σ fragmentation bounds ⟹ algebra) — all in
TournamentH7/FragmentationLemma.lean, build-verified. Then loop: further LRC(14)
proof/formalization tasks per the live queue.

-> THM-883, mac-mini-S123 draft, HYP-7106 (S30 restructure); death-star-S31.
