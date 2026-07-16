---
id: THM-899
title: THE LATTICE LAW OF FOUR-RUNNER BOX-HIT CONSTANTS (the c_B(k) evaluation, completed and corrected) — (i) SINGLE-RELATION CLOSED FORM: for a primitive full-support relation k ⊥ v and singleton boxes, c_B(k) = −(1/(24·k₁k₂k₃k₄))·Σ_{ε∈{±1}⁴}(∏ε)·B₄({(Σεᵢkᵢ−φ)/14}), φ = Σkᵢ(2cᵢ+1) — a RATIONAL number (Bernoulli B₄ at 14ths; e.g. c_B((1,1,−1,−1),0⁴) = 11/7203; verified against direct summation to 1e-8); (ii) THE CORRECTION: the constant attaches to the RELATION LATTICE L(v) = {k : k·v = 0}, not to one vector — near-AP quadruples (n, n+4, n+1, n+3) carry a RANK-2 lattice (basis (3,1,−4,0), (2,0,−3,1); note the TRIPLE relation) and the plateau is the lattice sum R∞ = Σ_{k∈L∖0, supp≥3} ∏ĝ_Bᵢ(kᵢ) = 4.685e-3, matching the measured 4.712e-3 (dominant strata: full-support 2.95e-3 + four support-3 at 3.1–5.6e-4); (iii) THE FINAL FORM OF THM-898's law: R(v,B) = Σ_{k∈L(v), supp≥3} ∏ĝ(kᵢ) + O(1/v), every term Bernoulli-closed (B₄ full support, B₃-weighted with ĝ(0) = |B|/7 factors at support 3) — the kernel's four-runner content is EXACTLY the Bernoulli mass of the speed set's relation lattice; THM-730 counts the rank-2 (E₃-degenerate) stratum
status: (i) PROVED (∏sin expansion + Σcos(2πtx)/t⁴ = −π⁴B₄({x})/3; exact-ℚ referee); (ii)/(iii) established computationally to matching precision (lattice sum vs measured plateau; strata resolved); the B₃ closed forms at support 3 are the one-line analogues (named); K6 normalization reconciliation with codex-THM-891 remains theirs/next
source: boxeph-2026-07-16-S34 (owner: take the c_B(k) evaluation; work the two named finishers)
depends_on: [THM-898 (the stratified law), codex-S17 THM-891 (the kernel), THM-730 (counts the degenerate stratum)]
script: inline in session (referee blocks); constants: c_B((1,1,-1,-1),0000) = 11/7203
---

# THM-899 — the lattice law

The week's B₂ machinery ({x}(1−{x}) = the Q_s kernel) has its fourth-moment sibling:
box-hit relation constants are B₄-rational. The single-vector form (i) is exact and
proved; the honest completion (ii): a quadruple's scale-free remainder is the SUM over
its whole relation lattice — near-AP families are rank-2 (the E₃-richness that makes
them extremal everywhere else makes their box-hit lattice fat too), and the measured
plateau decomposes exactly into the lattice's support-strata. Final law (iii): the
residue-6 kernel's four-runner content = pair D-data + the Bernoulli mass of L(v).
What remains for the full K6 = −12 closure: codex's normalization bridge (their
462-state frame ↔ these lattice constants) — flagged to them with all constants ready.
