# Independent audit of THM-4824, HYP-9428 and the two-grids results note

*2026-10-27. Independent adversarial audit by a subagent of mac-mini-2026-10-27-slab, of the pre-audit checkpoint
3b156ecd6d. The producer saved this report verbatim from the auditor's final message. Every number below was recomputed
with the auditor's own code (the scripts are in this directory).*

**Method.**
- Debt values come from a different algorithm (Knuth's generalised Dijkstra, `knuth_debt.py`) than the producer's
  all-pairs rounds. It was first checked against THM-4820's g(1..16).
- Each (A) construction was rebuilt as written, and the resulting trees were checked merge by merge.
- No existing file was edited and nothing was committed.

**Summary.** Most of THM-4824 holds up, but HYP-9428 (3) is refuted twice, and HYP-9428 (1) turns out to be provable.

## Per-claim verdicts

**(A) Free lattices: PASS.**
- **(A1) square grid:** correct. For a < 0 you need the symmetry x → −x, because the formula's β = −2(k+1) turns
  positive at a ≤ −3. This is presentational only.
- **(A2) D_n:** P(j) and the pairing of odd coordinates are correct. The construction as written gave zero-cost trees
  for D3–D8.
- **(A3) E8:** h′ is a root, x − h′ is in D8, and β ≤ −(3/8)Σ|xᵢ| − 2, which is attained. Simplification: any
  glue-coset root h with x·h ≤ 0 works, since −h is in the same coset.
- **(A4) table:** reproduced exactly (249, 601, 331, 531 and 2401 points, all g = 0). Larger windows also give 0: D3 to
  Q ≤ 24, D4 to 10, D5 to 6, A4 to 9, A5 to 6, E8 to Q ≤ 3 (9121 points).

**(B) Collinear law: PASS.** The potential cases are correct (opposite signs use min(a,|b|) ≤ a|b|). On Z[√−3] only
the real line is reachable and g(n) = n² − |n| for |n| ≤ 40. On L_4, g(2n, 0) = 4(n² − |n|).

**(C) Glued family: minimal vectors PASS; ratios reproduce; GAP on labels.**
- Checked exactly for c = 1/3 to 4, including 29/10 and 31/10. The c ↔ 1/c similarity and L_3 = 4 × the Eisenstein grid
  are correct.
- **The GAP:** "FINITE-EXACT: debt ratios at 3/2 and 5/2" is not justified. There is no potential or window lemma at
  those aspects, so the window values are only upper bounds (g ≤ g_W). The low ends of the stated ranges bound nothing.
- Evidence in favour: none of the 120m outer-half values changes at 240m or 480m.
- "κ(3/2) ∈ [0.129, 0.141]" presents a range of g_W/Q values as an interval for a lim sup.
- **Angle law:** "REFUTED" rests on window data, not a proof. The data are strong, though: at 4× windows the outer-half
  maximum falls to 0.1352 (c = 3/2, law 0.1667) and 0.2794 (c = 5/2, law 0.3).
- An extra scan of c = 11/10 to 11/4 stays below the law except at c = 2 and 3, rises with c, and tends to 0 as c → 1⁺.
  That supports HYP-9428 (2)'s conjectures, as window evidence only.

**(D) Aspect-2 floor: PASS.**
- Exhaustive box check (|u|, |v| ≤ 26, about 2·10⁶ pairs): Q mod 8 ∈ {0, 3, 4}, β is even, no two norm-4-mod-8 vectors
  are perpendicular, and the merge inequality holds.
- The window lemma is correct, and all 1001 points of norm ≤ 891 are certified.
- The 21 exceptional norms match exactly, with excess 4 only at 100 and 676, which are the only points of those norms.
- The on-axis row k ≤ 29 matches value for value.
- Unstated detail: at norm 747 only 4 of 12 points are exceptional, the (±13, ±17) orbit and not (±27, ±3).

**HYP-9428 (3): FAIL, both parts.** Both counterexamples are exact: the lower bound by the window lemma, the upper bound
by an explicit tree checked independently.
- **Prime law:** (62, 0) has g = 964 = 2⌊3844/8⌋ + 4, so p = 31 ≡ 7 (mod 8) has excess 4. The "only if" direction fails.
  - The "if" direction holds for every p ≡ 5 (mod 8) up to 53: 5, 13, 29, 37, 53.
  - Full row k ≤ 53: excess 4 at k = 5, 13, 16, 20, 29, 31, 35, 37, 53; excess 2 at 23, 46, 47, 49; 0 elsewhere.
- **"Excess 2 off the axis":** (±10, ±24), at norm 1252, has g = 316 = 312 + 4. There are 104 off-axis excess-4 points up
  to norm 11262, all certified.
- **The "why" does not separate the cases.** Primes p ≡ 7 (mod 8) are also inert in Z[√−2], since (−2/p) = 1 iff
  p ≡ 1, 3 (mod 8). Also, L_2 is not Z[√−2]. It is the index-2 sublattice u ≡ v (mod 2), whose form is
  3x² + 2xy + 3y², the non-principal class of discriminant −32.

**(E) Parity-coset shell lemma: PASS.**
- Brute force over every split of every x with 2 ≤ Q(x) ≤ 700 (2532 points) gives 0 mismatches with the shell
  criterion.
- Seeds up to 700 are 11, 23 and 25, matching THM-4820.
- The 11 and 51 seed tests reproduce.
- Session75's criterion is exactly the Q ≡ 1 (mod 3) case.

**(F) Zero-sum quartic identity: PASS, CLASSICAL** (a consequence of Newton's identities, checked symbolically with
sympy). The values below 51 are 2, 18, 32.
- **Typing GAP:** "the zero-sum condition is the catalyst condition read on the coefficients" is not an identity. The
  catalyst is the kernel direction (1,1,1), and the zero-sum plane is its orthogonal complement.
- x has a zero-sum representative iff 3 | a+b+c, that is iff Q(x) ≡ 0 (mod 3), and then Q = 3(a² + ab + b²). This should
  be labelled ANALOGY.

**(G) Results-note numerics: PASS for every listed count, list and table entry,** including the square-later counts and
the aliquot fibres (complete).
- **Gaps below 81:** there are five, not four. The fifth is 65–80, of length 16, which breaks the "16 − 5 + a" pattern.
  The producer's loop drops the trailing gap.
- **Side 7:** 3-5-7 is a primitive 120° triangle with 7 as its longest side. Both the note and HYP-9428 (4) need "7
  adjacent to the 120° angle".
- **50:** first only if (5,5) counts; with distinct squares the first is 65.
- **79:** "this is why 79 needs 19" is half the story. Mod 16 forces at least 15 summands; 19 comes from 79 < 3⁴, so only
  1 and 16 are available.
- **Typing:** "same simplex as THM-4820's charge vectors" and "two-grid obstruction" are untyped; they should be ANALOGY.

**(H) Typing and prior art.**
- **Root-system sentence: TRUE (classical).** A2 has no perpendicular roots; B2 has e₁ ⊥ e₂; G2 has short ⊥ long;
  rank ≥ 3 always has two non-adjacent Dynkin nodes.
  - Caveat: the unit charges are minimal vectors. G2's root lattice is the A2 lattice, whose minimal vectors are the
    short roots, and it pays. G2's perpendicular pairs are all short against long.
- **E8 as "the 8-dimensional form of the half-shift glue":** an analogy, and should be typed as one.
- **Prior art (from memory, not checked):**
  - I know of no literature on merge trees whose cost is the positive part of the pairwise interaction.
  - Without the positive part the total is the same for every tree (Σβ = Q(root) − ΣQ(leaves)). Its scalar case is the
    classical pile-splitting invariant: the sum of ab is n(n−1)/2.
  - The correspondence between splits and the coset x + 2L via v = y − z is the device behind Voronoi's
    relevant-vector criterion (Conway–Sloane, SPLAG).
  - g(4) = 19 and the 79 example are classical.
- **Minor:** `e8_check.py`'s docstring says "(E)" where it means (A3)/(A4).

## Upgrade: HYP-9428 (1) can be proved

With the roots as unit charges, a root lattice, reducible or not, is free iff it is neither A1 nor A2.
1. **Perpendicular roots:** if α ⊥ h are roots, then Zα + Zh is a square grid, so every kα is free by (A1).
2. **Support reduction in A_n:** take two coordinates xᵢ, x_t of the same sign and split off y = xᵢ(eᵢ − e_t). Then
   β = −xᵢ(xᵢ + x_t) < 0 and the support of the remainder shrinks by one. For n ≥ 3, the leftover root multiples are free
   by step 1, using e₁ − e₂ ⊥ e₃ − e₄. So every A_n with n ≥ 3 is free.
3. **Index-2 glue:** for x in the glue coset, choose ±h with x·h ≤ 0. Then x − h lies in the free sublattice and
   β ≤ −2. This handles:
   - E8 over D8;
   - E7 over A7, with E7 = {x ∈ Z⁸ ∪ (Z+½)⁸ : Σx = 0};
   - E6 over A5 ⊕ A1 (index 2, 40 roots in the glue coset, checked).
4. **Reducible lattices:** every root has a perpendicular partner in another component, so steps 1 and 2 apply.

**Evidence.**
- Explicit zero-cost trees were verified for A3–A9, E6, E7 and A2 ⊕ A1 on random points.
- Windows give g = 0 for A6 (Q ≤ 4), A7 (Q ≤ 3), E6 (Q ≤ 4), E7 (Q ≤ 3), A2 ⊕ A1 (Q ≤ 40) and A2 ⊕ A2 (Q ≤ 16).
- The FINITE-EXACT status for A4 and A5 is correct but superseded.

**Scripts.** The two refutations are `audit_p31.py` and `audit_offaxis_10_24.py`. The construction checks are
`audit_constructions.py`.
