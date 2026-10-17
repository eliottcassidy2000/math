# Independent audit of THM-4665 (EI17 over the primes of Q(zeta55))

2026-10-17. Independent referee subagent of opus-2026-10-17-S39/S40; the subagent could not write this file, so its
report is transcribed here unchanged in substance. Scripts and outputs: this folder.

**Overall verdict: CONFIRMED WITH FIXES.** Every count in statements 1, 2 and 4 was reproduced by a search that uses
none of the theorem's lemmas. The direction structure and the reduction lemma are correct. One sentence in the Reading
was wrong (MISTAKE-625). The other fixes concern citation scope, wording and provenance. The audit also identified the
"third involution" that the theorem left open.

| Item | Verdict |
|---|---|
| Direction structure (13 classes, gauge, U7..U12, 17 vertex formulas) | CONFIRMED |
| Exhaustiveness (parallelograms, elimination, degenerate branches) | CONFIRMED WITH FIXES (citation scope, wording) |
| Statement 1 | CONFIRMED |
| Statement 2 | CONFIRMED |
| Statement 3 (reduction lemma, gauge) | CONFIRMED (sharpening offered) |
| Statement 4 | CONFIRMED WITH FIXES (provenance; "mostly 8" made exact) |
| Reading: "obstruction is specific to the classes ±1 (mod 55)" | WRONG |
| Open: "third involution, not identified" | Resolved by this audit |
| No claimed exclusion of EI17 from L | CONFIRMED (theorem, note and HYP-9247 update) |

## 1. Direction structure (`a1`)

- **The graph.** PE has 17 vertices and 31 edges and is triangle-free. It is not 3-colourable but is 4-colourable
  (SAT), and it passes the Laman count. It matches the House of Graphs adjacency list in exactly 2 ways, and
  |Aut| = 2 (the automorphism fixes 17 and swaps 11 with 14).
- **The classes.** The nine 4-cycles, closed under the parallelogram identifications with signs, give 13 classes and
  no sign conflict. Independently, the 24-digit realization (refined to residual 3e-101) has 13 edge directions up to
  sign. The least gap between classes is 0.0057, and the partition and signs are the same.
- **The relations.** The cycle relations from a spanning tree have rank 6. An exact sympy solve expresses U7..U12 as
  integer combinations of the free directions. These are identical to the six formulas in the theorem, and the derived
  vertex formulas equal the scripts' VERT table exactly.
- **Degenerate branches.** The four collapses hold as identities: v5 = v11, v7 = v10, v2 = v11 and v1 = v6. Two
  branches also force v2 = v7 and v10 = v17.

## 2. Exhaustiveness

- **Parallelograms.** The forcing holds over any field of characteristic not 2, split or not.
  - Two unit circles whose centres differ by w ≠ 0 meet in at most two points: subtracting the norm equations leaves a
    line, and a line meets the irreducible conic N(x) = 1 at most twice.
  - The set {b−a, d−a} is closed under x ↦ w−x, and a fixed point would force b = d. Hence c−b = d−a.
  - **Gap:** THM-4663 (2) proves the codegree bound only for G_q, but THM-4665 cited it for "any field" and for the
    split model.
  - The tables also check codegree ≤ 2 for every c ≠ 0 at every prime run, in both models.
  - "Opposite steps equal" was loose wording: the steps satisfy u3 = −u1 and u4 = −u2.
- **Elimination.** Complete.
  - U1 and U2 are looped over all units, and U3, U6, U4 and U0 each come from a fully enumerated codegree equation.
  - U10 and U9 are checked; U7, U12, U11 and U8 are units by construction.
  - Both root formulas were re-derived.
- **Degenerate branches.** Handled correctly by both scripts: the split script enumerates them and the twisted script
  skips them, which is valid in both cases. Method C enumerates them in full and gets the same counts.

## 3. Independent recounts (`a2`, `a6`, `a7`)

- **Method B (lemma-free).** A direct search for injective embeddings in the gauge v17 = 0, v11 = −1. It uses no
  parallelograms, no direction classes, no relations and no quadratic formula.
  - It places vertices in a different order from S24 (free vertices 9 and 10 instead of 6 and 2), checking every edge
    and injectivity as it goes.
  - The twisted model uses a different representation of F_{p²} from the S39 script.
  - Apex vertices are solved from a full pair table, or by brute force over all units for p < 400.
- **Method C.** The U-system in a different elimination order.
- **Agreement.** 151 split primes and 123 twisted primes overlap with the S39 outputs, with **zero disagreements**, and
  all methods agree with each other.
- **Statement 1.** Zero at all 26 listed primes. Positive counts: 2861: 8, 3191: 16, 4621: 8, 4951: 8, 7481: 8,
  10781: 8.
- **Statement 2.** Zero at all 20 primes ≡ 54 (mod 55) below 6000.
- **Statement 4, split controls.** The audit also ran p = 3, 5 and 11, which S39 skipped. Points occur exactly at the
  11 listed primes below 700.
- **Statement 4, twisted controls.** Points occur exactly at the 14 listed primes below 600; p = 383 gives 4.
- **2441 and 2551.** The counts 2441 → 16 and 2551 → 8 are confirmed, but neither was in the committed twisted output.
  They are now cited from `a2_twisted_10mod11_table.out`.
- **Count values.** Positive counts take only the values 4, 8 and 16. All are multiples of 4, and 32 of the 33 are
  multiples of 8; the exception is p = 383 (twisted).
- **Direct checks (a6).**
  - Every twisted embedding found at 127, 149, 383, 2441 and 2551 was re-verified in pure Python as a subgraph of
    Cay(F_{p²}, μ_{p+1}), with μ built from a primitive element, and re-checked in the S39 script's own representation.
  - In the split model, the S39 script's U-tuple sets equal the audit's exactly at 103, 139, 677, 2861 and 3191.

## 4. Reduction lemma (`a4`)

- **Conjugation.** On Q(ζ55), complex conjugation is σ₋₁.
- **Split primes.** For the 32 primes ≡ 1 (mod 55) below 12100, Φ55 splits into 40 linear factors, so the residue
  field is F_p, and conjugation moves P. At all 40 primes P above p, a sample unit direction u reduces to
  (u mod P, ū mod P) = (x, 1/x).
- **Primes ≡ 54 (mod 55).** For all 20 primes, Φ55 has 20 quadratic factors. Conjugation lies in the decomposition
  group and equals Frobenius (y^p = y⁻¹, checked at 109 and 5939). Reductions of unit directions have norm 1.
- **Good reduction.** The definition is correct: for norm-one u, being a P-unit is equivalent to being a σ(P)-unit.
- **Sharpening.** At a conjugation-fixed P, v_P(u) = −v_P(ū) = −v_P(u), so every unit direction is automatically a
  P-unit. Bad reduction above the 20 twisted primes therefore always means a vertex collision.
- **Gauge.** Legitimate in both settings.
  - In L: translate, then rotate by ū5 = 1/u5, which lies in L.
  - In G_p: the affine maps z ↦ az + b with a ∈ μ_{p+1} are arc-transitive, so "not a subgraph of G_p" follows.

## 5. Typing and overclaims

The Reading's sentence is contradicted by the counts in the controls, which exclude the classes ±1 (mod 55).

- **Zero is the norm.** 108 of 119 split controls and 91 of 105 twisted controls have no point.
- **Split class 1 (mod 55).** Enriched, not depleted: 6 of 32 primes have points, against a control rate of 0.092
  (binomial tail 0.070).
- **Twisted class 54 (mod 55).** Only mildly low: 0 of 20, which has probability 0.057 at the control rate.

The note's side observation is confirmed: there are no split points at the 38 primes ≡ 10 (mod 11) below 2600
(probability 0.025). Among those primes the twisted model has points only at 2441 and 2551.

## 6. The third involution (`a8`, `a9`)

- **The map.** v1 + v14 = v2 + v11 holds identically; call the common midpoint m. The block re-placement
  ι: (U0, U4) ↦ (U9, U11) moves v4 to 2m − v4 and v6 to 2m − v6, and v13, v16 follow.
  - ι fixes U7, U10 and U12, swaps U0 with U9 and U4 with U11, and sends U8 to −U8, so it preserves the solution set.
  - ι has no fixed point, since a fixed point would force v4 = v6.
- **Even counts.** The mirror has no fixed point for p > 17: a fixed point would make every direction ±1 around an odd
  cycle. Hence faithful counts are even for p > 17 (PROVED).
- **The group.** The mirror, the graph automorphism and ι are commuting involutions generating (Z/2)³.
  - It acts freely on faithful solutions at all 29 tested primes with points, except p = 383 (twisted), where
    mirror∘automorphism fixes all four.
  - ι mapped faithful solutions to faithful ones at every tested prime; this is not proved in general.
  - The involution explains S24's two real realizations that differ only in the vertices {4, 6, 13, 16}.

## Fixes applied (2026-10-17, S40)

1. Status: parallelogram forcing by the two-circle bound in characteristic not 2 (THM-4663 (2) is its G_q case); every
   degenerate branch forces a vertex collision in every model.
2. Setting: the correct 4-cycle argument (opposite steps are negatives).
3. Statement 3: the conjugation-fixed sharpening (every unit direction is a P-unit there).
4. Statement 4: 2441 (16) and 2551 (8) as the only twisted points among the 38 primes ≡ 10 (mod 11) below 2600, cited
   to this audit's outputs; "multiples of 4; all but p = 383 multiples of 8".
5. Reading (theorem and note): the corrected base-rate statement (MISTAKE-625).
6. Open (theorem and note): the third involution ι, even counts for p > 17, and the (Z/2)³ action.
7. Credit: THM-4558 Lemma R (KNOWN: Madore arXiv:1509.07023, Prop. 3.2) cited beside THM-4628.

## MISTAKE candidate

- **MC1 (minor, interpretive) → MISTAKE-625.** The claim that "the obstruction is specific to the classes ±1
  (mod 55)", made in the theorem's Reading and the results note, is contradicted by the counts in section 5.
- **Fixes, not mistakes.** Citing THM-4663 (2) beyond G_q, the "opposite steps equal" wording, and the missing outputs
  for 2441/2551 and the ≡ 10 (mod 11) runs; all the underlying claims are true.

## Prior art (repository grep)

- **EI17.** There is no earlier finite-field treatment. Every earlier mention (S24's LLL work, the EI19/fish and EI21
  notes, THM-4661/4662, the six-colour notes, HYP-9247) is numerical or concerns field placement. S24 §5 proposed
  reduction modulo good primes, and THM-4665 credits that.
- **Embeddings in G_p.** No earlier computation embeds a given unit-distance graph in G_p. THM-4552, 4558, 4624, 4628
  and 4663 and the logic4 folder study automorphisms, colourings, spectra and walk counts of G_q.
- **Reduction mechanism.** THM-4628 (a), the same as THM-4558's Lemma R (KNOWN).

## Process notes

- The a6 run wrote numba cache files into the shared, gitignored `__pycache__`; the referee deleted only its own 12
  files and redirected the cache to `%TEMP%`.
- One quoted Bash heredoc containing apostrophes was used, without harm.
- Nothing was committed by the referee, and no tracked file was modified.
