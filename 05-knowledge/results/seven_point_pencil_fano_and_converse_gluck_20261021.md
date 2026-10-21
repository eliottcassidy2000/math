# The seven-point sphere is two spheres, Fano's seventh line, and Gluck's theorem with the converse added

Session opus-2026-10-21-S49.
* **Owner prompt (verbatim):** "consider how tournament structure relates to the open mathematical questions we investigate and see how the fano plane is analogous to the below, and keep merging in a diverse set of ideas for inspiration on continued hypotheses". It came with a pasted "seven-point sphere" statement for a tetrahedron and D. Gluck, *Trivial set-stabilizers in finite permutation groups*, Can. J. Math. 35 (1983) 59–67.
* **Scripts.** All are in `04-computation/experiments/seven_point_spheres_20261021/`, each with a `.out`, each printing ALL CHECKS PASSED:

| script | contents | runtime |
|---|---|---|
| `seven_point_pencil.py` | (P1) Gram-basis proof, (P2) rational tetrahedra, (P3) `F_p`, (P4) `Q(i)`, (P5) other dimensions, (P6) Fano anchor | 7 s |
| `converse_gluck_tournaments.py` | (G1) Paley `P_q`, (G2) `P_27`, (G3) all tournaments `n ≤ 9`, (G4) `n = 15` substitution checks, (G5) all circulants `n ≤ 15` | 40 s |

* **Canon:** THM-4742 and THM-4743 (PROVED + FINITE-EXACT), HYP-9353 (OPEN). HYP-9354 was reserved and released unused.
* **Concurrent work.** The same prompt reached the session23b lane, which committed four packets on 2026-10-21 (§1). This note builds on them and states what is new relative to each.

**Status.**
* PROVED:
  * the six-label sphere `S2`;
  * the pencil;
  * the defect `Σ` and its tangency reading;
  * the field dichotomy;
  * the orbit-pairing lemma;
  * the primitive converse-Gluck classification;
  * Paley tournaments;
  * the substitution rule;
  * the infinite exceptional family.
* FINITE-EXACT: the checks above.
* ANALOGY: §4.
* OPEN: HYP-9353.

## 0. Results

1. **The pasted claim is false (session23b), and its true form is a pencil of two spheres (THM-4742).**
   * The six labels `A, B, D, E, F, H'` always lie on one sphere `S2`. Session23b's sphere `S1` carries `B, C, D, E, F`.
   * Both spheres are invariant under the inversion of radius `2R` about `O`. They meet in a circle through `B, D, E, F` in the plane `O G O'`.
   * The vertex `C = V0` misses `S2` by exactly

         Σ = |V0V1|² + |V0V2|² + |V0V3|² = −8 (V0 − O)·(G − V0).

   * So the seventh incidence holds iff the median `V0G` is tangent to the circumsphere, iff `A = C`.
2. **This is the Fano phenomenon in metric form.**
   * Six incidences are identities. The seventh is controlled by one quantity that is positive over `R`.
   * For Fano's seventh line that quantity is the determinant `2`, which vanishes only in characteristic 2. Here it is the sum of squares `Σ`.
   * `Σ` vanishes for genuine tetrahedra exactly over non-real fields: `Q(i)`, and `F_p` for `5 ≤ p ≤ 97`, in every sampled case. Even then it vanishes only with the coincidence `A = C`.
3. **Gluck's theorem with the converse added (THM-4743).** Tournaments always have Gluck subsets, since their automorphism groups have odd order (KNOWN). Admitting anti-automorphisms:
   * the exceptional tournaments (those needing a third colour) are exactly `C_3`, `Z_5{1,2}`, `P_7` when the group is primitive;
   * the Paley tournament `P_q` is exceptional iff `q = 3, 7`;
   * every `X[C_3^[k]]` is exceptional, an infinite family.
4. **Censuses (FINITE-EXACT).**
   * All tournaments with `n ≤ 9` (191,536 classes at `n = 9`): the exceptions are `C_3, Z_5{1,2}, P_7, C_3[C_3]`.
   * All circulants on `Z_n`, `n` odd `≤ 15`: exactly the family members.
   * HYP-9353 conjectures that the family is complete.
5. **Merging (§4).** Three typed bridges and one numerology guard:
   * the **sign involution** (tournament converse, `3x+1 ↔ 3x−1`, sphere inversion), which no odd-symmetry observer can break in exactly the exceptional cases;
   * the **seventh incidence** (Fano, the sphere, Collatz cycles: algebra permits and positivity forbids);
   * **`C_3` substitution** as the common grammar of rigid and of maximally symmetric self-converse tournaments.

## 1. What the session23b lane had already done today

* **`tetrahedral_inversion_spheres_20261021`:**
  * the universal claim is false, with a rational counterexample;
  * the four-point obstruction for `A, C, D, E`;
  * the real exceptional locus (`D = E`, `B = F`; five distinct points);
  * the generic sphere `S1 ⊃ {B, C, D, E, F, 4C}` with its inversion invariance;
  * `Q1(H') = 6 − 2 v·s` and `Q1(A) = μ(16 − T)`, both positive;
  * a dimension-free obstruction; `A, C, O, E` concyclic.
* **`fano_tournament_move_quotients_20261021`:**
  * a combinatorial tetrahedron→Fano bridge (six edges plus a formal point; `F7\ω = M(K4)`);
  * the real non-representability of Fano;
  * the line-flip cube `Q7 → K8` syndrome cover;
  * `P_7`'s 84 regular subsets; the order-42 and order-168 groups have none;
  * no tournament anti-automorphism preserves a single Fano plane.
* **`signed_tournament_palette_rigidity_20261021`:**
  * three colours always suffice up to converse, and are needed for `P_7`;
  * a Gluck subset leaves at most one anti-involution;
  * distinguishing numbers 2, 3, 4 for the order-21, -42, -168 actions;
  * a fixed-slack transfer formula.
* **`paley_multicoloring_and_carry_20261021`:** palette laws (not used here).

**What is new here:**
* the second sphere `S2`, the pencil, the defect `Σ` with its tangency/coincidence reading, the field dichotomy, and the three-dimensionality check;
* which tournaments need the third colour: the orbit-pairing lemma, the primitive classification, Paley, the substitution rule, the tower family, and the two censuses.

## 2. The seven-point sphere (THM-4742)

* **Labels.**
  * `A`: the second point of line `V0G` on the circumsphere.
  * `B`: the pole of the face plane.
  * `C = V0`.
  * `D = O + 4(G − O)`.
  * `E = inv(G)`.
  * `F = O + 4(O' − O)`.
  * `H'`: the orthocentre of the face `V1V2V3`.
* **Incidences** (generic tetrahedra; `ι` = inversion of radius `2R` about `O`):

| sphere | labels on it | also on it | defect of the missing labels |
|---|---|---|---|
| `S1` (session23b) | B, C, D, E, F | `4C = ι(C)` | `Pow(H') = Σ`; `Pow(A) = μ(16 − T)` |
| `S2` (new) | A, B, D, E, F, H' | `ι(A) = 4A`, `ι(H')` | `Pow(C) = Σ` |
| `S1 ∩ S2` | the circle Γ through B, D, E, F | | radical plane = plane `O G O'` |

* **Proof shape.**
  * Both spheres have constant term `4R²`; with `O = 0` that is power `4R²` at `O`. This is forced for any sphere through the inverse pair `D, E` lying on a ray from `O`.
  * Each sphere is then fixed by three dot products with the basis `v = V0`, `s = V1 + V2 + V3`, `w = O'`.
  * Every incidence becomes an identity in `v·s`, `|s|²`, `v·w` and `|w|²` ((P1)).
  * The key identities are `Σ = 6R² − 2 v·s` and `μ = Σ/|s − 3v|²`.
* **Fields.**
  * Over an ordered field `Σ > 0`, so the seventh label never joins.
  * Over `F_p` (`5 ≤ p ≤ 97`) in every one of 50,025 valid random configurations (1,388 with `Σ = 0`), and over `Q(i)` in the example `V3 = (0, i, i)`: all seven lie on one quadric iff `Σ = 0` iff `A = C`.
  * Seven distinct cospherical labels therefore never occur in the generic regime over any field.
  * Session23b's real locus has `D = E` and `B = F`, so it is not generic.
* **Dimensions** ((P5), random exact trials).
  * For 4- and 5-simplices, with factors `n + 1` and `H' = O' + n(G' − O')`, no more than `n + 1` labels are ever cospherical.
  * In the plane, `A, B, C, E` and `B, D, E, F` are concyclic.
  * The pencil is a three-dimensional phenomenon: it needs the orthocentre of a triangular face.
* **The Fano comparison.**
  * The textbook real picture of `PG(2,2)` (vertices, midpoints and centroid of a triangle) draws six lines straight and the seventh as a circle. The seventh line has determinant `−2` ((P6)).
  * The seven-label sphere likewise draws two spheres where one was claimed. The obstruction is `Σ`, a sum of squares.
  * Neither seven-label system is identified with the Fano points: session23b's §8 warns against that, and I agree. The common shape is "one positive quantity separates the real geometry from an algebraically consistent incidence".

## 3. Gluck with the converse (THM-4743, HYP-9353)

* **KNOWN.**
  * Tournament automorphism groups have odd order, so Gluck's Corollary 1 gives every tournament a vertex subset with trivial stabilizer, i.e. distinguishing number at most 2 (Albertson–Collins; Meslem–Sopena; Lozano).
  * Session23b adds that three colours always suffice when anti-automorphisms are admitted.
* **Which tournaments need the third colour?** Call them *exceptional*.

| result | statement | type |
|---|---|---|
| orbit pairing | exceptional iff an anti-automorphism fixes every regular `Aut`-orbit on subsets | PROVED |
| primitive `Aut±` | exceptional iff `C_3`, `Z_5{1,2}` or `P_7`: Gluck's Theorem 1 leaves only `S_3, D_10, F_42` among orders `≡ 2 (mod 4)` | PROVED |
| Paley | `P_q` exceptional iff `q ∈ {3, 7}`: `Aut±(P_q) ⊇ AGL(1,q)` is 2-transitive | PROVED |
| substitution | `A` prime, `B` strong and exceptional: `A[B]` exceptional iff `D(Aut±(A)) > r(B)` | PROVED |
| towers | `r(C_3^[k]) = 2`, so every `X[C_3^[k]]` (`X = C_3, Z_5{1,2}, P_7`) is exceptional | PROVED |
| census | exceptions with `n ≤ 9`: `C_3, Z_5{1,2}, P_7, C_3[C_3]`; circulants with odd `n ≤ 15`: the family members only | FINITE-EXACT |
| completeness | the family is all of them | HYP-9353 |

* **Regular-orbit counts and distinguishing numbers.**
  * `r(B)`, the number of regular `Aut(B)`-orbits on subsets, is 2, 6, 4, 2, 6 for `C_3`, `Z_5{1,2}`, `P_7`, `C_3[C_3]`, `Z_5{1,2}[C_3]`.
  * The three prime exceptions need exactly three colours.
  * So only a triangle tower, with `r = 2`, can be substituted into them without escaping. `C_3[Z_5{1,2}]` (`n = 15`) has 15,000 regular subsets, while `Z_5{1,2}[C_3]` has none.
* **Paley regular-orbit counts.** `Aut(P_q)` has 2, 4, 36, 3052, 33156, 4618100 regular orbits for `q = 3, 7, 11, 19, 23, 31`; `Aut±(P_q)` has 0, 0, 12, 1470, 16392, 2306868.
  * The closed form (G1) is `2^q − 2 − q Σ_(d | m, d > 1) (−μ(d))(2^(1+(q−1)/d) − 2)`, where `m` is the order of the centre stabilizer.
  * Brute force confirms it for `q ≤ 19`.

## 4. Tournament structure and the open questions: what merges

Each bridge names what transfers and what does not.

* **The sign involution (ANALOGY with PROVED instances).** Three objects carry an involution that preserves all the "odd" or "algebraic" structure but reverses a sign:
  * the converse of a self-converse tournament: for `P_7`, `x ↦ −x`, which swaps its two Fano planes (session23b);
  * the conjugation `x ↦ −x` taking `3x + 1` to `3x − 1` (THM-4610's remark: every Haar statement survives it, which is the KNOWN barrier);
  * the inversion `ι` of the sphere pencil, which pairs `D ↔ E` and `B ↔ F`.

  For tournaments THM-4743 decides exactly when observers (vertex subsets) fail to break the sign: only in the Gluck circulants and their triangle towers. In those cases odd symmetry is broken, but the anti-involution survives every marking.

  The Collatz reading is that Haar data are odd-symmetry observers. They cannot break `3x ± 1`, and an archimedean input is required. That part is the KNOWN barrier; no new Collatz statement is claimed.
* **The seventh incidence (ANALOGY).** In three places an algebraically consistent incidence is killed by positivity:
  * Fano's seventh line, killed by `2 ≠ 0`;
  * the sphere's seventh label, killed by `Σ > 0`;
  * nontrivial positive cycles of `3x + 1`, killed by integrality plus size.
    * A rational cycle point `c_w/(2^A − 3^l)` exists for every parity word.
    * THM-4610's (S5) shows that a 1-cycle shape has a single candidate numerator for `p = 2`, the equation Steiner's theorem excludes, but `(p−1)^k2` candidates for base `p ≥ 3`, where extra cycles do appear.

  What transfers is the shape of the proof obligation: find the positive quantity. What does not transfer is any mechanism.
* **`C_3` substitution (PROVED instances).**
  * THM-2450: rigid self-converse classes are generated by `stack` and `C_3[T1, T2, T3]`.
  * THM-4743: the exceptional classes are `X[C_3[C_3[…]]]`.
  * The same grammar has rigidity at one end and maximal sign-rigidity at the other.
  * A candidate quantitative bridge for HYP-9352's presentation counts: compare `tc(T) = H/|Aut|` with `r(T)` along the grammar. For `P_7`, `tc = 9` and `r = 4`. No inequality is claimed.
* **Numerology guard.** The exceptional orders include 7 and 21 (`P_7`, `P_7[C_3]`), and THM-1370's H-spectrum gaps are 7 and 21. No mechanism links these; recorded as NUMEROLOGY only.
* **No transfer.** As session23b notes, none of this touches LRC or the chromatic number of the plane: there is no target map.

## 5. Next steps

1. **HYP-9353:**
   * self-converse tournaments with `n = 10, 11`;
   * the 21-vertex vertex-transitive files (`04-computation/circtourn21.d6`, `othervttourn21.d6`), where `P_7[C_3]` should be the only exceptional class;
   * Cayley tournaments on `Z_3^3`;
   * unequal substitutions.
2. **Geometry.**
   * The four vertex choices give eight spheres through `D` and `E`.
   * Do the four circles `Γ_i` (planes `O G O'_i`) have a common sphere or a common radical point? A first scratch test found that `D` and `E` are not on the circumsphere of the tangential tetrahedron.
3. **Is the defect always a sum of squares?** For which label systems built from a simplex by centroids, circumcentres, poles, inversions and homotheties is every failing incidence's defect a sum of squares? This is a testable form of "real geometry forbids exactly the non-real incidences".

## 6. Status

| item | status |
|---|---|
| `S2 ⊃ {A,B,D,E,F,H'}`, pencil, `Σ`, tangency, field dichotomy (THM-4742 1–4) | PROVED |
| `S1`, the counterexample, the real exceptional locus | session23b (credited) |
| three-dimensionality; `F_p`, `Q(i)` instances | FINITE-EXACT |
| converse-Gluck: lemma, primitive classification, Paley, substitution, towers (THM-4743 1–5) | PROVED |
| censuses `n ≤ 9` (all), circulants `n ≤ 15` | FINITE-EXACT |
| completeness of the exceptional family | OPEN (HYP-9353) |
| sign involution, seventh incidence, `C_3` grammar bridges | ANALOGY |
| 7/21 coincidence with THM-1370 | NUMEROLOGY |
