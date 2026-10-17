# Audit report: S38 six-colour / Hadwiger–Nelson material (THM-4661, THM-4662)

Independent adversarial referee, 2026-10-16 (subagent of opus-2026-10-16-S38). Worktree `math-wt-explicit-modal-38`.
Audited at `6b7eecb78a`, then extended to `4971cf2b05` (THM-4661 statements 4(c) and 5, THM-4662) at the coordinator's
request. No existing file was edited. All scripts and outputs are in this folder, one `.out` per `.py`.

## Summary

I found no mathematical errors in THM-4661 statements 1–5 or THM-4662 statements 1–5. Everything computer-assisted
reproduces in independent implementations, including the decisive two-descent. The errors are in the remarks, the
note, the results INDEX and the probe script's LLL section. The two serious ones: the "LLL coupling lengths" are wrong
by a factor of about 1,000, and the claim that 18- and 19-vertex five-cores "exist combinatorially" is false.

**THM-4662 holds (the most important item).** My own two-descent over Q(√5) reproduces it:
- Sel² of E over Q(√5) has 4 elements and Sel² of E^δ over Q(√5) has 16, giving the same Selmer subspaces as the S38
  script.
- Every local image is maximal isotropic for the local Tate pairing, a structural check that does not use the expected
  sizes.
- A residue-symbol test at 59 primes of Q(√5) that split in F_4 shows the four image classes are independent and
  (33,33) is not in their span.
- So rank E(F_4) = 2, C3(F_4) is empty, and the fish and EI21 have no placement in Q(ζ55).

**C3 is locally solvable everywhere over F_4; justification corrected 2026-10-17.** The original sentence that −11 is a square at *every* finite place was false. At 19 the quartic defining F_4 has four simple roots 1,4,15,18, so its completions are Q_19, where −11 is nonsquare. C3 instead has finite local points there: D=10 and (z0,z1,z2)=(4,6,6) modulo 19 lift by Hensel. The complete proof is B6 of [cubic adelic coverage](../../../05-knowledge/results/cubic_adelic_20261016e.md): infinity points at the bad places 2,3,7, and the even-valuation Kummer image at every other finite place. The unchanged conclusion is a nonzero class in Sha(E/F_4)[2]. If Sha is finite, that part has dimension at least 2. See the [exact hostile control and correction lineage](../../../05-knowledge/results/anchor_deficits_cubic_hub_20261017a.md).

**THM-4661 statement 5 reproduces.** With my own pattern enumeration, encoding and solver:
- I get the same pattern counts (1, 7, 38 and 135 at n = 19, 20, 21, 22).
- Every pattern is UNSAT at n = 19, 20, 21.
- All 36 models at n = 22 satisfy the necessary conditions N1–N5, and their Alon–Tarsi coefficients reproduce exactly.

**THM-4661 4(c) and THM-4662 statement 5 reproduce** from the House of Graphs adjacency lists: Γ_9 occurs exactly at
EI19 vertex 12 and EI21 vertex 15, and in neither EI17 nor the fish.

## ERRORS

### E1 (moderate–high): the "LLL coupling relations" are artifacts
- **Where:** note §5 ("LLL-reduced integer forms of length 19,416 to 4.9 million … Powers of one conjugate do much
  better", "Coupling needs more than first powers of the conjugates"); results INDEX ("LLL length at least about 19,000"
  and "the short coupling is clock_duality_oct15's 64-term circuit"); the probe script's docstring and JSON field
  `coupling_relations_lll`; commit 00116ae780 ("≥ 19,416").
- **Cause.** `fmpz_mat.nullspace()` returns a basis of the rational null space. Its integer span is a sublattice of the
  integer relation lattice of index about 10^115, 10^109, 10^107 and 10^192 for the four triples (Gram-determinant
  ratios). LLL on that sublattice says nothing about the shortest relation.
- **Saturated lattice.** Kernel-LLL on [10^12·B | I] gives relations of L1 length 20, 70, 80 and 88. The S38 recipe,
  rerun, gives 19,416, 34,652, 37,330 and 4,922,578.
- **The length-20 relation for the triple (1,3,4)**, with ζ = ζ55 and v_r = σ_r(v):

  −ζ⁴ + v₁(ζ⁴ − ζ¹⁹ + ζ²⁹ + ζ³⁹) − v₃(1 + ζ¹¹ + 2ζ¹⁹ + ζ²² + ζ³³ + ζ³⁹) + v₄(1 + ζ⁴ + ζ¹¹ + ζ¹⁴ + ζ²² − ζ²⁴ + ζ³³ + ζ³⁴) = 0.

  It is exactly zero both in my arithmetic (basis ζ^i w^k with w = 2θ) and in the S38 script's own θ-basis arithmetic
  with its own v. It amounts to two 10-step unit walks with a common endpoint, built only from rotated first-power
  conjugates. That is shorter than C2's 64-term circuit.
- **Corrected statement.** The probes are Cartesian because their direction sets are linearly independent (Q-ranks 50,
  and 44 for probe 4), not because coupling relations are long. Rotated conjugates ζ^k v_r couple in 20 steps.
- **Effect on conclusions:** none. Walk banks built from the 20-step circuit (up to 1,279 vertices, 9 non-generator
  edges) are three-colourable, so the board's conclusion stands.
- **Ledger lesson:** a rational null-space basis is not an integer relation lattice; use kernel-LLL or HNF.

### E2 (moderate): five-cores of 18 and 19 vertices do not exist
- **Where:** THM-4661 Remarks ¶1 ("not excluded (only as failed exteriors)"); note §2 ("Five-cores of 18 or 19 vertices
  are not excluded"); note §6 ("canonical hostile … exists combinatorially").
- **Fact.** No graph on 18 or 19 vertices is triangle-free, has every pair sharing at most two neighbours, and has
  minimum degree at least 5.
  - **19 vertices.** The structure around a degree-6 vertex is forced (pair graph = octahedron K_{2,2,2}) and has no
    completion: zero by CaDiCaL enumeration and zero by plain backtracking. Dropping the root-neighbour /
    second-neighbour constraint leaves 48 labelled completions; dropping the second-neighbour constraint leaves 100.
    So both constraint families are used and the zero is not an encoding accident.
  - **18 vertices.** The graph must be 5-regular. All four rooted cases are UNSAT (SAT and backtracking). Positive
    control: the same code recovers the Clebsch graph, uniquely, at 16.
- **Consequence.** Every set of at most 19 exterior points spans a 4-degenerate unit graph (15 or fewer by counting; 16
  by the bridge note's planar Gram argument; 17 by parity; 18 and 19 by enumeration). This gives statement 1 without
  Borodin–Erdős–Rubin–Taylor. Statement 5 makes the remark moot, but the three sentences should be replaced.

### E3 (minor): "conjugation-stable" dropped
THM-4662 proves the fish/EI21 exclusion only for conjugation-stable carriers M ⊇ K (the proof builds z_i from c̄ and
d̄; for non-conjugation-stable M the compositum M·M̄ can have even degree over K). The qualifier is missing in the INDEX
(twice), the note's status line and its Γ_9 sentence, FRONTIER-CONNECTION-SIDECARS, the THM-4661 title, THM-4662
Remark 1, and commit messages 00116ae780, 4cb4cab5fc. EI19 is unaffected.

### E4 (minor): THM-4661 3(d) needs the ports in K'
"every equality tester whose old vertices lie in K' …" is vacuous; (d) needs x, y ∈ F and V(G) ⊂ M. Same gap in the
Remark "Its exterior must already pass statements 2 and 5". Fix: "every equality tester G ⊂ L' with ports x, y ∈ K' has
a failed exterior inside V(G) ∖ K', so |V(G) ∖ K'| ≥ 22".

### E5 (minor): note §4, "general remark"
"Lies … in the smallest field containing that edge" is true only for placements inside an odd-degree extension of a
conjugation-stable F containing the edge. "Rigid graphs … are the only ones this does not touch" overclaims (e.g. the
fish's genus-nine C4 is untouched).

### E6 (minor): note §5, "five complete assemblies built around that circuit"
Only rows 4 and 5 of the C2 note's table are built from the circuit; rows 1–3 are power, mixed and conjugate banks. All
five being three-chromatic is correct.

### E7 (minor): probe-script docstring
"LLL relations … explain why the probes are Cartesian": the Q-ranks explain it; the LLL lattice is not saturated.

## IMPRECISIONS

- **I1** 2(b) is the list version of Gallai's theorem (Kostochka–Stiebitz–Wirth 1996; Borodin); cite KSW.
- **I2** "the same point of C3" (note, HYP-9247, INDEX): the theorem gives existence of some F_4-point only.
- **I3** (I1), (I2) for non-conjugation-stable M are correct; worth one sentence.
- **I4** "Inherited from the cyclotomic part": probes 1–2 map to C5 and C11 via residue maps plus the sum of
  v-coordinates; probe 4 via its eleven pentagonal factors; probe 3's cyclotomic projection is three-coloured only by SAT.
- **I5** Checks script: the cubic's root mod 661 is simple (661 ∤ 592); further certificates −35 at 881, −2 at 991, 2311.
- **I6** Statement 3(b): the chain needs m ≥ 2 (automatic for d < 1/2).
- **I7** First-moment paragraph: correctly labelled and computed (17.655; average degrees 3.647, 3.684, 3.652; probes
  5.72, 5.58, 4.02, 4.04).
- **I8** Step (v) correctly counts low degree-6 vertices in W_6 (tight too); say so.
- **I9** The even-q table checks only steps (iv)–(v); steps (i)–(iii) generalize verbatim (checked).
- **I10** Board row "equality tester ⟺ χ = 7" holds only as ⟹ for a tester inside L'.
- **I11** Statement 5 typed as solver verdicts (honest); Glucose corroborates; DRAT would make it checkable.
- **I12** Note §2: "choice number at most five" refers to the 36 models found; "every other required property" should
  read "the necessary conditions N1–N5".
- **I13** Status lines say "independent audit pending"; update.
- **I14** INDEX credits clock_duality_oct15 only for the circuit; add the seventeen-vertex theorem, G7, C1.

## Requested checks on the three new statements

**(1) THM-4661 statement 5.** Root choice and degree bound correct; pattern enumeration complete (1, 7, 38, 135
reproduced by brute-force canonical forms; the two root-degree-7 patterns at n = 22 not enumerated); conditional
Gallai count correct; N5 cuts valid (none used at n ≤ 21, so the UNSAT verdicts certify N1–N4 alone); the Alon–Tarsi
sign convention is irrelevant to non-vanishing.

**(2) THM-4661 4(c).** The four 4-cycles at hub 15 and the edges 10–18, 16–17 are present in EI21 (HoG 51378);
symbolic algebra gives |P|² + |Q|² = 1, Re(P Q̄) = 0 and the three square relations of C3 with X = |v13 − v20|².

**(3) THM-4662.** Local class maps (2-adic: unit square classes mod 8, Hensel verified mod 64; split places over 11:
exact valuations by division in Z[φ]; ramified √5: exact); all local images maximal isotropic (dimensions 1, 4, 2, 2, 2
of local H¹ dimensions 2, 8, 4, 4, 4); rank E(F_4) = rank E(K0) + rank E^δ(K0) standard, dimension count 2 + 2 = 4;
F_4 square test agrees with residue symbols and the theorem's hand argument; P1, P2 exact with the claimed classes;
the covering map has zero residuals; the consequence for Q(ζ55) follows through 4(b), 4(c) and [L⁺ : F_4] = 5.

## VERDICTS

| # | Claim | Verdict |
|---|---|---|
| 1 | 2(b): Borodin/ERT, Gallai trees, (5/4)(ℓ−1), high vertex | CONFIRMED |
| 2 | Counting lemma (2); bound 1 + q + ⌈q(q−2)/2⌉ | CONFIRMED |
| 3 | 2(d): 19 excluded | CONFIRMED; stronger: no 19-vertex five-core exists |
| 4 | 2(c); (I1)–(I3) for non-conjugation-stable M | CONFIRMED |
| 5 | Statement 3 | CONFIRMED WITH FIXES (E4) |
| 6 | 4(a), certificate at 661; 4(b) | CONFIRMED (qualifier issues: E3) |
| 7 | Probe script | v, conjugates, edges, zero non-generator edges, ranks 40/50/120/44, Cartesian reading: CONFIRMED; LLL lengths and conclusions from them: REFUTED (E1) |
| 8 | Note and INDEX | CONFIRMED WITH FIXES (E1, E2, E3, E5, E6) |
| N1 | Statement 5 | CONFIRMED (independent reproduction) |
| N2 | 4(c): EI21 carries C3 | CONFIRMED |
| N3 | THM-4662, statements 1–5 | CONFIRMED (independent descent); the remark's local solvability is TRUE |

## Scripts (this folder)

audit_nineteen_vertex_exhaustive.py, audit_nineteen_vertex_backtrack.py, audit_eighteen_vertex_five_cores.py,
audit_eighteen_vertex_backtrack.py, audit_cubic_probes_independent.py, audit_short_coupling_relations.py,
audit_short_circuit_walk_banks.py, audit_field_and_counting_checks.py, audit_ei19_fish_structure.py,
audit_ei21_gamma9.py, audit_exterior_sat_independent.py, audit_alon_tarsi_independent.py,
audit_two_descent_independent.py; exploratory (inconclusive except the Gallai filter): audit_twenty_vertex_gallai_filter.py
(none of 3,000 sampled 20-vertex five-cores passes the Gallai filter), audit_twenty_vertex_degree6_count.py (SAT
timeouts), audit_twenty_vertex_cegar_fast.py (undecided within its cap). CaDiCaL 1.9.5 in this pysat build ignores
`interrupt()`, so time-limited queries used Glucose.

## Not checked

C1's valuation argument beyond its arithmetic; G7's width claims; C2's 64-term relation (cited only); the fish note's
2-descent over Q (superseded by the descent over K0); geometric realizability of any 20–22-vertex model.
