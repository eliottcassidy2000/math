---
id: HYP-9428
title: "After THM-4824: (1) RESOLVED, a root lattice pays a catalytic merge debt iff it is A1 or A2 (THM-4824 (G)); (2) the debt ratio kappa(c) of the (1/2,1/2)-glued family is positive and increasing for 1 < c <= 3 (window evidence; the drafted angle law (c-1)/(2c) is refuted on window evidence); (3) the aspect-2 exceptional set (REFUTED as drafted: excess 4 also at (62,0) with p = 31 = 7 mod 8, and off the axis); (4) the 13/37 coincidences are numerology, each side being an elementary identity"
status: >
  Partly RESOLVED, partly REFUTED, the rest OPEN (mac-mini-2026-10-27-slab, 2026-10-27; independently audited, corrections
  applied: 04-computation/experiments/eisenstein_two_grids_fourth_powers_20261027/audit/REPORT_THM4824.md).
  (1) PROVED as THM-4824 (G), the audit's upgrade. (2) OPEN, with window evidence only: the windowed values are upper
  bounds g <= g_W. kappa(2) >= 1/4 and kappa(3) >= 1/3 are PROVED by potentials. (3) as drafted REFUTED twice: excess 4
  at (62, 0) (p = 31 = 7 mod 8), and off the axis at (+-10, +-24). The corrected description is OPEN. (4) a typing
  verdict, recorded in the results note.
session: mac-mini-2026-10-27-slab
source: 05-knowledge/results/eisenstein_two_grids_13_37_fourth_powers_20261027.md
related:
  - 01-canon/theorems/THM-4824-catalytic-merge-debt-by-lattice-shape-square-checkerboard-and-e8-free-first-grid-collinear-eisenstein-the-critical-half-shift-glue-parity-coset-shells-zero-sum-quartics.md
  - 01-canon/theorems/THM-4820-eisenstein-merge-debt-floor-q-over-three-free-neutral-catalysts-and-the-eleven-forty-one-threshold.md
  - 05-knowledge/hypotheses/HYP-9426-slab-exactness-black-endpoint-counts-eisenstein-exceptional-set-and-catalyst-thresholds.md
---

# HYP-9428 — which lattices pay a catalytic merge debt, and the 13/37 coincidences

1. **A₂ is the only root lattice that pays — RESOLVED (THM-4824 (G)).** With the roots as unit charges, a root lattice
   (reducible or not) is free iff it is neither A₁ nor A₂. The proof uses three devices:
   - perpendicular helpers make every root multiple free;
   - support reduction handles A_n with n ≥ 3;
   - index-2 glue handles E₆ ⊃ A₅ ⊕ A₁, E₇ ⊃ A₇ and E₈ ⊃ D₈.

   The mechanism is the one anticipated here: among the ADE lattices, A₂ is the only irreducible one of rank ≥ 2 without
   two perpendicular roots.
2. **The glued family's debt ratio.** For L_c = {u ≡ v mod 2}, Q = u² + c v², put κ(c) = lim sup g/Q.
   - **Known.** κ(1) = 0, κ(2) ≥ 1/4, κ(3) ≥ 1/3, and the lattice is collinear for c > 3.
   - **Window evidence only.** Windowed values are upper bounds g ≤ g_W. At c = 3/2 the outer-half maximum is 0.1352,
     and at c = 5/2 it is 0.2794, at 4× windows. A scan from c = 11/10 to 11/4 stays below the angle law except at c = 2
     and 3, rises with c, and tends to 0 as c → 1⁺ (audit).
   - **Conjecture.** κ(c) > 0 for every c ∈ (1, 3], κ is increasing, and κ(c) → 0 as c → 1⁺.
   - **Refuted on window evidence.** The first-drafted law κ(c) = (c − 1)/(2c), equivalently −cos θ/(1 − cos θ) for the
     angle θ between the glue vectors.
   - **Open question.** Whether κ is a continuous function of the shape, or depends on the arithmetic of the form for
     rational c.
3. **The aspect-2 exceptional set (REFUTED as drafted; open in corrected form).**
   - **The draft.** It predicted excess 2 off the axis, and excess 4 on the axis at (2p, 0) exactly for primes
     p ≡ 5 (mod 8).
   - **Refuted on the axis.** (62, 0) has g = 964 = 2⌊3844/8⌋ + 4, so p = 31 ≡ 7 (mod 8) has excess 4. The full row
     k ≤ 53 has:
     - excess 4 at k = 5, 13, 16, 20, 29, 31, 35, 37, 53;
     - excess 2 at k = 23, 46, 47, 49;
     - excess 0 elsewhere.

     So the "if" direction (p ≡ 5 mod 8 gives excess 4) holds for p = 5, 13, 29, 37, 53.
   - **Refuted off the axis.** (±10, ±24) at norm 1252 has excess 4, and there are 104 off-axis excess-4 points up to norm
     11262 (audit, all certified).
   - **The wrong reason.** Primes ≡ 7 (mod 8) are also inert in Z[√−2]. And L₂ is not Z[√−2]: it is the index-2
     sublattice with the form 3x² + 2xy + 3y², the non-principal class of discriminant −32.
   - **Open.** Find the arithmetic description of the aspect-2 exceptional set. The analogue of THM-4820 (C), a seed
     criterion via the shell lemma for L₂, is the natural first step.
4. **The 13/37 coincidences are numerology.** Each statement below is an elementary IDENTITY (results note §3). The
   coincidences between them are NUMEROLOGY:
   - **Four fourth powers.** 13 of 1..50 are sums of four fourth powers, and 37 are not.
   - **Square-later.** Session69 finds 37 square-later evens and 13 failures up to 100.
   - **Lattice balls.** 13 and 37 are lattice-ball counts in both the square and the Eisenstein grid.
   - **Primes.** 13 and 37 are the two smallest primes ≡ 1 (mod 12).
   - **Triangles.** They are the primitive hypotenuses with leg 12, and the longest sides of the primitive 120° triangles in which 7 is a side adjacent to the 120° angle. (3-5-7 has 7 as its longest side.)

   **The structural part is not numerology.** opus-2026-10-27-S66's THM-4825 (concurrent, pre-audit) characterizes
   {5, 13, 37} as exactly the squarefree d ≡ 1 (mod 4) for which the owner's two grids of Z[√−d] form the whole class
   group, that is h(−4d) = 2. The NUMEROLOGY verdict covers only the coincidences between the elementary facts listed
   above, not the recurrence of 13 and 37 as such.

   **The prediction, already confirmed.** No two of these sequences keep matching. At the doubled cutoffs the square-
   later split of the evens up to 200 is 81/19, and the four-fourth-powers split of 1..100 is 21/79.
