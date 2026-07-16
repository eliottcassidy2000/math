# HYP-7141 — THM-883's Lean formalization: closing the three sorries

**Status:** CLAIMED / IN PROGRESS (death-star-2026-07-16-S31; owner directive: finish the
formalization + loop LRC tasks). Verify-first.

Targets: (1) window_bound (round-to-nearest a₀ = ⌊wz+½⌋ reduces the ℤ-union to two arcs;
three-case clip: both-clip telescopes to 2λ/w + ℓ − 1/w ≤ 2λ/w; interior arc ⟹ other arc
misses); (2) fragmentation main branch (⌈wL⌉ ≤ wL+1 window chop + finite subadditivity);
(3) killer_budget (cover ⟹ vol ≤ Σ fragmentation bounds ⟹ algebra) — all in
TournamentH7/FragmentationLemma.lean, build-verified. Then loop: further LRC(14)
proof/formalization tasks per the live queue.

-> THM-883, mac-mini-S123 draft, HYP-7106 (S30 restructure); death-star-S31.
