---
id: HYP-9428
title: "After THM-4824: (1) every irreducible root lattice other than A1 and A2 has zero catalytic merge debt (with the roots as unit charges), so A2, the Eisenstein grid, is the only root lattice with a positive floor; (2) the debt ratio kappa(c) of the (1/2,1/2)-glued family is positive for 1 < c <= 3 (the drafted angle law (c-1)/(2c) is refuted); (3) the aspect-2 exceptional set: excess 2 off the axis, excess 4 on the axis at (2p, 0) for primes p = 5 (mod 8); (4) the 13/37 coincidences across the four-fourth-powers count, session69's square-later count and the two lattices are numerology, each side being an elementary identity"
status: >
  OPEN (mac-mini-2026-10-27-slab, 2026-10-27; not audited). (1) PROVED for D_n (n >= 3, including A3 = D3) and E8
  (THM-4824 (A)); FINITE-EXACT for A4 (Q <= 6) and A5 (Q <= 4); OPEN for A_n (n >= 6), E6, E7. (2) kappa(2) >= 1/4 and
  kappa(3) >= 1/3 PROVED by potentials; measured values in THM-4824 (C); the angle law REFUTED at c = 3/2 and 5/2.
  (3) FINITE-EXACT below norm 891 (THM-4824 (D)); the on-axis prime prediction holds for every prime p <= 29
  (aspect2_axis.out, exact by the window lemma), and the composite on-axis exceptions need their own law.
  (4) a typing verdict, recorded in the results note.
session: mac-mini-2026-10-27-slab
source: 05-knowledge/results/eisenstein_two_grids_13_37_fourth_powers_20261027.md
related:
  - 01-canon/theorems/THM-4824-catalytic-merge-debt-by-lattice-shape-square-checkerboard-and-e8-free-first-grid-collinear-eisenstein-the-critical-half-shift-glue-parity-coset-shells-zero-sum-quartics.md
  - 01-canon/theorems/THM-4820-eisenstein-merge-debt-floor-q-over-three-free-neutral-catalysts-and-the-eleven-forty-one-threshold.md
  - 05-knowledge/hypotheses/HYP-9426-slab-exactness-black-endpoint-counts-eisenstein-exceptional-set-and-catalyst-thresholds.md
---

# HYP-9428 — which lattices pay a catalytic merge debt, and the 13/37 coincidences

1. **A₂ is the only root lattice that pays.** With the roots as unit charges, every irreducible root lattice except A₁
   (collinear, THM-4824 (B)) and A₂ (the Eisenstein floor ⌊Q/3⌋, THM-4820) has g ≡ 0.
   - **Mechanism.** A₂ is the only irreducible root system of rank ≥ 2 without two perpendicular roots, and
     perpendicular roots are what the free constructions of THM-4824 (A) use.
   - **What is proved.** D_n for n ≥ 3, including A₃ = D₃, and E₈.
   - **What is computed.** A₄ (Q ≤ 6) and A₅ (Q ≤ 4) are free on their windows.
   - **Test.** A₆ and E₆ on windows Q ≤ 3, or a uniform construction for A_n, n ≥ 3. The natural candidate is the D_n
     induction with e_a − e_b in place of e_a ± e_b, which needs four distinct indices.
2. **The glued family's debt ratio.** For L_c = {u ≡ v mod 2}, Q = u² + c v², put κ(c) = lim sup g/Q.
   - **Known values.** κ(1) = 0, κ(2) ≥ 1/4, κ(3) ≥ 1/3; for c > 3 the lattice is collinear.
   - **Measured.** κ(3/2) ∈ [0.129, 0.141] and κ(5/2) ∈ [0.271, 0.290] (window estimates, THM-4824 (C)).
   - **Conjecture.** κ(c) > 0 for every c ∈ (1, 3], κ is increasing, and κ(c) → 0 as c → 1⁺.
   - **Refuted.** The first-drafted law κ(c) = (c − 1)/(2c), equivalently −cos θ/(1 − cos θ) for the angle θ between the
     glue vectors, fits c = 1, 2, 3 and fails at c = 3/2 and c = 5/2.
   - **Open question.** Whether κ is a continuous function of the shape, or depends on the arithmetic of the form for
     rational c.
3. **The aspect-2 exceptional set.** Exceptions below norm 891 have excess 2, except the on-axis points (±10, 0) and
   (±26, 0), which have excess 4.
   - **Prediction.** For a prime p, the point (2p, 0) has excess 4 iff p ≡ 5 (mod 8).
     - **Why.** These primes are inert in Z[√−2] (the norm form u² + 2v²), as 11 and 41 are inert in Z[ω] for the
       Eisenstein case.
     - **Check below 891.** Excess 4 occurs at p = 5 and 13, while 7 ≡ 7 (mod 8) gives no excess.
   - **Test, passed for primes ≤ 29** (`aspect2_axis.out`, exact by the window lemma).
     - **Excess 4** exactly at p = 5, 13 and 29, which is the prediction confirmed at the new point (58, 0).
     - **Excess 2** at p = 23 (≡ 7 mod 8).
     - **Excess 0** at every other prime.
   - **Composite k with excess 4.** k = 16 and k = 20, while 10, 15, 25 and 26 have none. The composite law is open.
4. **The 13/37 coincidences are numerology.** Each statement below is an elementary IDENTITY (results note §3). The
   coincidences between them are NUMEROLOGY:
   - **Four fourth powers.** 13 of 1..50 are sums of four fourth powers, and 37 are not.
   - **Square-later.** Session69 finds 37 square-later evens and 13 failures up to 100.
   - **Lattice balls.** 13 and 37 are lattice-ball counts in both the square and the Eisenstein grid.
   - **Primes.** 13 and 37 are the two smallest primes ≡ 1 (mod 12).
   - **Triangles.** They are the primitive hypotenuses with leg 12, and the primitive 120°-triangle sides with side 7.

   **The prediction, already confirmed.** No two of these sequences keep matching. At the doubled cutoffs the square-
   later split of the evens up to 200 is 81/19, and the four-fourth-powers split of 1..100 is 21/79.
