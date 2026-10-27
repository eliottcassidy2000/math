---
id: THM-4824
title: "Catalytic merge debt by lattice shape: with the minimal vectors as unit charges, a root lattice is free (zero debt at every point) iff it is neither A1 nor A2 (the square grid, D_n, A_n for n >= 3, E6, E7, E8 = D8 u (D8 + (1/2)^8) and every reducible root lattice are free); a lattice with collinear minimal vectors (the first grid Z[sqrt-3], units +-1) costs n(n-1) on its line and nothing else is reachable; in the (1/2,1/2)-glued rectangular family {u = v mod 2, Q = u^2 + c v^2} the Eisenstein grid c = 3 is the critical member, where the glued points tie with the base vectors and six units appear; the aspect-2 member has the floor 2 floor(Q/8) with its own exceptional points; tight Eisenstein splits are the points of one parity coset on two or four norm shells; zero-sum triples have a^4 + b^4 + c^4 = 2(a^2 + ab + b^2)^2"
status: >
  PROVED (mac-mini-2026-10-27-slab; independently audited 2026-10-27, corrections applied:
  04-computation/experiments/eisenstein_two_grids_fourth_powers_20261027/audit/REPORT_THM4824.md):
  (A) the square grid, D_n (n >= 2) and E8 are free, by explicit constructions;
  (B) the collinear law;
  (C) the critical-aspect statement for the glued family;
  (D) the aspect-2 floor g >= 2 floor(Q/8) and the window lemma;
  (E) the parity-coset shell lemma;
  (G) the classification of root lattices: free iff not A1 or A2. This is the audit's upgrade of HYP-9428 (1), by
  support reduction for A_n, index-2 glue for E6, E7, E8, and perpendicular helpers for reducible lattices.
  CLASSICAL: (F) the zero-sum quartic identity.
  FINITE-EXACT:
  - the aspect-2 exceptional norms below 891 and the on-axis row;
  - the windowed freeness checks (A4).
  WINDOW EVIDENCE ONLY (upper bounds g <= g_W): the debt ratios at aspects 3/2 and 5/2, and hence the refutation of
  the angle law.
  Inherited: THM-4820 and Session75's parity-shell criterion. OPEN: HYP-9428.
session: mac-mini-2026-10-27-slab
source: 05-knowledge/results/eisenstein_two_grids_13_37_fourth_powers_20261027.md
scripts:
  - 04-computation/experiments/eisenstein_two_grids_fourth_powers_20261027/grid_debt.py (+ .out)
  - 04-computation/experiments/eisenstein_two_grids_fourth_powers_20261027/root_lattice_debt.py, root_lattice_debt2.py, e8_check.py (+ .out)
  - 04-computation/experiments/eisenstein_two_grids_fourth_powers_20261027/aspect_ratio_debt.py, aspect2_exceptions.py (+ .out)
  - 04-computation/experiments/eisenstein_two_grids_fourth_powers_20261027/counts.py (+ counts.out, section 4)
related:
  - 01-canon/theorems/THM-4820-eisenstein-merge-debt-floor-q-over-three-free-neutral-catalysts-and-the-eleven-forty-one-threshold.md
  - 05-knowledge/results/session75_final_mathematics_20261027/EISENSTEIN_INVOICES.md
  - 05-knowledge/hypotheses/HYP-9426-slab-exactness-black-endpoint-counts-eisenstein-exceptional-set-and-catalyst-thresholds.md
  - 05-knowledge/hypotheses/HYP-9428-root-lattices-free-except-a1-a2-glued-family-debt-ratio-aspect-two-exceptions-and-the-13-37-coincidences.md
  - 01-canon/theorems/THM-4825-the-half-shifted-second-grid-is-a-ring-or-the-order-two-class-over-two-the-two-grids-are-the-class-group-iff-d-is-5-13-or-37.md
  - 01-canon/theorems/THM-4826-which-grid-carries-a-norm-shell-one-to-two-split-cubic-character-of-two-centred-units-complete-the-quartic-designs-four-eisenstein-fourth-powers-miss-five-mod-eight.md
---

# THM-4824 — catalytic merge debt by lattice shape, and the Eisenstein grid as the critical half-shift glue

**Setting.** This generalises THM-4820 from Z[ω] to an arbitrary lattice Λ with a positive definite integral form Q.
- **Unit charges.** The unit charges are the minimal vectors of Λ.
- **Trees and costs.** A merge tree has unit charges as leaves, and a merge of y and z costs β(y, z)₊, where
  β = Q(y + z) − Q(y) − Q(z) = 2⟨y, z⟩.
- **Catalytic debt.** g(x) is the minimum cost of a tree with root value x, using any number of leaves. Neutral
  combinations act as free catalysts, as in THM-4820 (A2).
- **Free lattices.** Λ is **free** if g ≡ 0 on Λ.

**The owner's picture, exactly.** Write Z[ω] in the coordinates (u, v) of (u + v√−3)/2:
- **The two grids.** Z[ω] = {u ≡ v (mod 2)}. It is the union of the first grid Z[√−3] (u, v both even) and a copy of
  it centred at (½, ½) in the first grid's units (u, v both odd).
- **The norm.** The norm is Q = (u² + 3v²)/4. Points of the shifted copy always have odd norm.
- **The units.** The six units are ±1 in the first grid and ±ω, ±ω² in the shifted copy.

## (A) Free lattices (PROVED)

**Note on 0.** In every case 0 = r + (−r) for a root r, with β = −2Q(r) < 0, so 0 is free.

**(A1) The square grid Z[i].** g ≡ 0. Every real integer a is free, by strong induction on |a|:
- **Even a.** a = 2k = k(1+i) + k(1−i), with β = 0. Here k(1 ± i) = k + (±ki), also with β = 0.
- **Odd a.** a = 2k + 1 = (k+1)(1+i) + (k − (k+1)i), with β = −2(k+1) < 0. Both parts are perpendicular merges of
  free real and imaginary integers.

This covers a > 0, and the symmetry x → −x covers a < 0. Unit rotations then make every imaginary integer free, and a + bi = a + (bi) is a perpendicular merge. ∎

**(A2) Checkerboard lattices D_n.** Here D_n = {x ∈ Z^n : Σxᵢ even}, Q = x·x/2, and the roots are ±eᵢ ± eⱼ. Every D_n
with n ≥ 2 is free.
- **n = 2.** D₂ is a rotated, scaled square grid, so (A1) applies.
- **n ≥ 3, root multiples.** Let P(j) say that every j·(e_a ± e_b) is free.
  - **Even j.** j e_a = (j/2)(e_a + e_c) + (j/2)(e_a − e_c), with β = 0. So j(e_a ± e_b) = j e_a ± j e_b is a
    perpendicular merge.
  - **Odd j ≥ 3.** j(e_a + e_b) = (j e_a + e_c) + (j e_b − e_c), with β = −1, using a third index c. Each part splits
    perpendicularly, for example j e_a + e_c = ((j+1)/2)(e_a + e_c) + ((j−1)/2)(e_a − e_c).

  Strong induction on j proves every P(j).
- **n ≥ 3, general points.** Pair the odd coordinates of x, which come in an even number. A pair contributes
  xᵢeᵢ + xⱼeⱼ = ((xᵢ+xⱼ)/2)(eᵢ + eⱼ) + ((xᵢ−xⱼ)/2)(eᵢ − eⱼ), a perpendicular merge. Each even coordinate contributes
  xᵢeᵢ. Pieces with disjoint supports are perpendicular, so x is a merge of free perpendicular pieces. ∎

**(A3) E8 = D8 ∪ (D8 + (½)⁸).** E8 is free. [ANALOGY: this is the 8-dimensional version of the half-shift glue.] More simply, any glue-coset root h with x·h ≤ 0 works, since −h lies in the same coset; see (G), step 4.
- **Integer points.** Points of D8 are free by (A2).
- **Glue points.** Let x lie in the glue coset. Take h′ = −½·sign(x), coordinatewise; the coordinates of x are half-odd,
  so none is 0. If h′ has an odd number of minus signs, flip the sign at a coordinate of least |xⱼ|.
- **Why the split works.** Then h′ is a root (an even number of minus signs), and x − h′ ∈ D8. Moreover
  β = x·h′ − 2 ≤ −(3/8)Σ|xᵢ| − 2 < 0. So x = (x − h′) + h′ is a free merge. ∎

**(A4) Finite checks** (`root_lattice_debt2.py`, `e8_check.py`; the audit extends them to D3 Q ≤ 24, D4 Q ≤ 10, D5 Q ≤ 6, A4 Q ≤ 9, A5 Q ≤ 6, E8 Q ≤ 3, A6 Q ≤ 4, A7 Q ≤ 3, E6 Q ≤ 4, E7 Q ≤ 3, A2 ⊕ A1 Q ≤ 40, A2 ⊕ A2 Q ≤ 16). All find 0 everywhere:

| lattice | roots | window | points |
|---|---|---|---|
| D3 = A3 | 12 | Q ≤ 12 | 249 |
| D4 | 24 | Q ≤ 7 | 601 |
| A4 | 20 | Q ≤ 6 | 331 |
| A5 | 30 | Q ≤ 4 | 531 |
| E8 | 240 | Q ≤ 2 | 2,401 |

## (B) Collinear minimal vectors (PROVED)

Suppose the minimal vectors are ±r.
- **Reachability.** Only the line Zr is reachable.
- **The cost.** g(n r) = (n² − |n|) Q(r). The all-r tree costs Σ 2ab·Q(r) = (n² − n)Q(r). The potential
  (k² − |k|)Q(r) on kr satisfies the merge inequality: with equality for equal signs, and as ≤ 0 for opposite signs.
- **Consequence.** Catalysts (r, −r) are useless.

For the first grid alone, Z[√−3] with units ±1, this gives g(n) = n(n − 1). That is THM-4820 (A3)'s f(n, 0, 0): before
the glue, nothing is gained.

## (C) The glued family and its critical member (PROVED)

Let L_c = {(u, v) ∈ Z² : u ≡ v (mod 2)} with Q_c = u² + c v². The family is self-similar under c ↔ 1/c, so take
c ≥ 1. The shortest vectors are:
- **1 ≤ c < 3:** the four glue vectors (±1, ±1), of norm 1 + c;
- **c = 3:** these together with (±2, 0), all of norm 4: six minimal vectors, and L_3 is the Eisenstein grid;
- **c > 3:** only ±(2, 0), so the lattice is collinear in the sense of (B).

So the Eisenstein grid is exactly the aspect at which the points centred at (½, ½) become as short as the first grid's
basis vector, and the last aspect at which they are minimal at all. Along the family:

| c | debt |
|---|---|
| 1 | free, the square grid (A1) |
| 2 | floor 2⌊Q/8⌋, see (D) |
| 3 | ⌊Q/3⌋, THM-4820 (here scaled by 4) |
| > 3 | the collinear law (B) |

Windowed values g_W/Q on the outer half of a window of 120 minimal norms (`aspect_ratio_debt.py`). These are **upper bounds** on g/Q at those points. At c = 3/2 and 5/2 no potential or window lemma is available, so they bound nothing from below:

| c | 1 | 3/2 | 2 | 5/2 | 3 |
|---|---|---|---|---|---|
| g/Q | 0 | [0.129, 0.141] | [0.245, 0.256] | [0.271, 0.290] | [0.328, 0.333] |

A first-drafted angle law, ratio (c − 1)/(2c), fits c = 1, 2, 3. It is refuted on window evidence, not by proof: at 4× larger windows the outer-half maxima are 0.1352 at c = 3/2 (law 0.1667) and 0.2794 at c = 5/2 (law 0.3) (audit). It is withdrawn; see HYP-9428 (2).

## (D) The aspect-2 floor (PROVED lower bound; FINITE-EXACT exceptions)

On L_2 we have Q mod 8 ∈ {0, 3, 4}, the four units have norm 3, and β is even. Put Φ = 2⌊Q/8⌋. Then

    Φ(y + z) − Φ(y) − Φ(z) = 2⌊(r + s + β)/8⌋ ≤ β₊.

- **β ≥ 2.** This follows from r + s ≤ 8.
- **β ≤ 0.** The only threat is r = s = 4 with β = 0, and two vectors of norm ≡ 4 (mod 8) are never perpendicular:
  their inner product is 4 × (an odd integer).

Hence g ≥ 2⌊Q/8⌋.

**Window lemma.** A tree of cost C has every node at Φ ≤ C, that is Q ≤ 4C + 7. So windowed values are exact, as in the
THM-4820 audit.

**Exceptions.** In the window Q ≤ 900, which is exact for every point of norm ≤ 891 (`aspect2_exceptions.py`), the
exceptional norms are 83, 100, 107, 147, 211, 227, 251, 307, 347, 419, 475, 507, 547, 587, 643, 659, 676, 747, 787,
811, 827.
- **Excess 4** at 100 and 676, the on-axis points (±10, 0) and (±26, 0).
- **Excess 2** at all the others. At 747 only the (±13, ±17) orbit is exceptional, not (±27, ±3).

Beyond 891 the audit finds more points of excess 4, both on the axis (for example (62, 0)) and off it (for example (±10, ±24) at norm 1252); see HYP-9428 (3).

## (E) Tight Eisenstein splits are parity-coset shells (PROVED)

**The bijection.** For x ∈ Z[ω], the splits x = y + z correspond bijectively to the points v = y − z of the coset
x + 2Z[ω]. The inverse is y = (x + v)/2, z = (x − v)/2. By the parallelogram law,

    Q(v) = Q(x) − 2β(y, z).

A split is nontrivial iff v ≠ ±x. With THM-4820 (B), the split is tight iff:
- **Q(x) ≡ 1 (mod 3):** Q(v) − Q(x) ∈ {0, 2}, the perpendicular and β = −1 splits;
- **Q(x) ≡ 0 (mod 3):** Q(v) − Q(x) ∈ {−2, 0, 2, 4}, that is β ∈ {1, 0, −1, −2}.

The residue conditions of the tight cases then hold automatically.

**The two grids.** The four cosets of 2Z[ω] refine the owner's two grids: the first grid Z[√−3] is the union of the
classes 0 and 1, and the shifted grid is the union of ω and ω².

**Consequences.**
- **A finite seed test.** A point x is a seed (no tight split) iff its coset meets the two or four shells only in ±x.
  This needs only norm-form arithmetic in one coset.
- **Session75's case.** Its parity-shell criterion is the Q(x) ≡ 1 case.
- **Examples.** 11 is a seed: the 121-shell meets 11 + 2Z[ω] only in ±11, and 123 = 3·41 is not a norm. So is 51:
  the shells 2601, 2599, 2603, 2605 meet its coset only in ±51.

## (F) The zero-sum quartic identity (CLASSICAL)

If a + b + c = 0, then a² + b² + c² = 2(a² + ab + b²) and

    a⁴ + b⁴ + c⁴ = (a² + b² + c²)²/2 = 2(a² + ab + b²)² = (2/9) Q(a + bω + cω²)²,

using Q(a + bω + cω²) = a² + b² + c² − ab − bc − ca. This was checked for |a|, |b| ≤ 60 (`counts.py` §4).

**The link to THM-4820** [ANALOGY]. The fourth-power sums of zero-sum triples are exactly twice squared Eisenstein norms (a² + ab + b² = Q(a − bω)). Below 51 these are 2, 18 and 32.
- **How it relates to the catalyst.** The zero-sum plane is the orthogonal complement of the catalyst direction (1, 1, 1), not the catalyst itself. A point x has a zero-sum coefficient representative iff Q(x) ≡ 0 (mod 3), and then Q(x) = 3(a² + ab + b²) (audit).
- **What is not claimed.** The merge-debt and Waring problems are not identified.

## (G) Which root lattices pay (PROVED; the audit's upgrade of HYP-9428 (1))

**Theorem G.** Let L be a root lattice, that is, an orthogonal sum of A, D, E root lattices, normalized so that roots
have Q = x·x/2 = 1. Take the roots as unit charges. Then L is free iff L is neither A₁ nor A₂.
- **A₁** is collinear: (B) applies.
- **A₂** is the Eisenstein grid: it pays ⌊Q/3⌋ plus its exceptional set (THM-4820).
- **Rank 2 caveat.** B₂ has the square lattice, which is free. G₂'s root lattice is the A₂ lattice, whose minimal vectors
  are the short roots, so it pays.

*Proof.*
1. **Perpendicular helpers.** If α ⊥ h are roots, then Zα ⊕ Zh is a square grid whose units ±α, ±h are roots of L. By
   (A1), every point of it is free in L, and in particular every root multiple kα.
2. **A_n for n ≥ 3.** Here A_n = {x ∈ Z^(n+1) : Σxᵢ = 0}.
   - **Root multiples.** Every root has a perpendicular root (e₁ − e₂ ⊥ e₃ − e₄), so every root multiple is free by
     step 1.
   - **Support reduction.** If x has support of size ≥ 3, two of its nonzero coordinates xᵢ, x_t have the same sign.
     Put y = xᵢ(eᵢ − e_t) and z = x − y. Then supp z = supp x \ {i}, and β(y, z) = y·z = −xᵢ(xᵢ + x_t) < 0.
   - **Induction.** Induct on the support: support 2 is a root multiple, and support 0 is 0.
3. **D_n for n ≥ 4** is (A2).
4. **E₆, E₇, E₈ by index-2 glue.**
   - **The glue step.** Let M ⊂ L have index 2. Suppose M is free with its own roots, and the other coset contains a
     root h. Then for x ∉ M, both x − h and x + h lie in M. Choose the sign so that x·h′ ≤ 0 for h′ = ±h. Then
     β(x − h′, h′) = x·h′ − 2 ≤ −2, a free merge.
   - **The three cases.** E₈ ⊃ D₈, E₇ ⊃ A₇ and E₆ ⊃ A₅ ⊕ A₁ have index 2: the determinant ratios 4/1, 8/2 and 12/3 all
     equal 4 = 2².
   - **Roots in the glue coset.** These are proper root subsystems, so the glue coset contains roots. The audit counts 40
     of them for E₆.
   - **The sublattices are free.** D₈ by step 3, A₇ by step 2, A₅ ⊕ A₁ by step 5.
5. **Reducible lattices.**
   - **Root multiples.** Every root has a perpendicular root in another component, so every root multiple is free.
   - **A₂ components.** Every point of an A₂ component is (a + b)r₁ + b·r₃ for two roots at 120° and a, b ≥ 0. That is
     a free merge, with β = −(a + b)b ≤ 0.
   - **Other components** are free by steps 2–4.
   - **General points.** A general point is a merge of mutually perpendicular free pieces, one per component. ∎

The audit also verified explicit zero-cost trees for A₃–A₉, E₆, E₇ and A₂ ⊕ A₁ on random points.

**Prior art (from memory).** No literature on merge trees with positive-part pair costs is known to the audit.
- **Without the positive part.** The total cost is the same for every tree, Q(root) − Σ Q(leaves). Its scalar case is
  the classical pile-splitting invariant Σab = n(n−1)/2.
- **The split–coset correspondence** of (E) is the device behind Voronoi's relevant-vector criterion (Conway–Sloane,
  SPLAG).

