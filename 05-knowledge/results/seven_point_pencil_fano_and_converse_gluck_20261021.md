# The seven-point sphere is two spheres, Fano's seventh line, and Gluck's theorem with the converse added

Session opus-2026-10-21-S49.

* **Owner prompt (verbatim):** "consider how tournament structure relates to the open mathematical questions we investigate and see how the fano plane is analogous to the below, and keep merging in a diverse set of ideas for inspiration on continued hypotheses".
* **Attachments:** a pasted "seven-point sphere" statement for a tetrahedron, and D. Gluck, *Trivial set-stabilizers in finite permutation groups*, Can. J. Math. 35 (1983) 59–67.

**Scripts.** All are in `04-computation/experiments/seven_point_spheres_20261021/`, each with a `.out`, and each prints ALL CHECKS PASSED.

| script | contents | runtime |
|---|---|---|
| `seven_point_pencil.py` | (P1) Gram-basis proof for the tetrahedron; (P1n) the same identities for every `n`, homogeneous in `R²`; (P2) rational tetrahedra; (P3) `F_p`; (P4) `Q(i)`; (P5) dimensions `n = 2…6` in exact Gram coordinates; (P6) the Fano anchor | 7 s |
| `converse_gluck_tournaments.py` | (G1) Paley `P_q`; (G2) `P_27`; (G3) all tournaments `n ≤ 9`; (G4) `n = 15` substitution checks; (G5) all circulants `n ≤ 15` | 40 s |
| `audit_S49/` | the independent audit's own code and outputs | minutes; the `n = 11` census is 15 min |

The audit code covers:
* the geometry in every dimension and over fields;
* the self-converse census for `n ≤ 11`;
* the vertex-transitive 21-tournaments and the Cayley tournaments of order 27;
* the Paley groups computed from the tournaments.

**Canon.**
* THM-4742 and THM-4743: PROVED + FINITE-EXACT, independently audited (MISTAKE-650).
* HYP-9353: PROVED by the concurrent Session24 lane (cba114930c, independently audited; [structural proof](converse_exception_classification_20261021/session24-structure/STRUCTURE.md)). It was OPEN at this note's checkpoint.
* HYP-9354 was reserved and released.

**Concurrent work.** The same prompt reached the session23b lane. It committed four packets on 2026-10-21 (§1). This note builds on them and states what is new relative to each.

**Status.**
* PROVED:
  * for every simplex dimension `n ≥ 3`: the six-label sphere `S2`, the pencil, the defect `Σ` and its tangency reading, and the field dichotomy;
  * the orbit-pairing lemma, and the fact that exceptional tournaments are strong;
  * the primitive converse-Gluck classification;
  * Paley tournaments;
  * the substitution rule and the infinite exceptional family.
* FINITE-EXACT: the checks above, plus the audit's censuses.
* ANALOGY: §4.
* PROVED (Session24): HYP-9353, completeness of the exceptional family. This note's evidence and analogies do not supply it; the Session24 proof does.

## 0. Results

1. **The pasted claim is false (session23b). Its true form is a pencil of two spheres (THM-4742).**
   * Assume `O, V0, G, O'` are not coplanar.
   * Then the six labels `A, B, D, E, F, H'` lie on one sphere `S2`. Session23b's sphere `S1` carries `B, C, D, E, F`.
   * Both spheres are invariant under the inversion of radius `√(n+1) R` about `O`. They meet in a circle through `B, D, E, F` in the plane `O G O'`.
   * The vertex `C = V0` misses `S2` by exactly

         Σ = Σ_i |V0Vi|² = −2(n+1)(V0 − O)·(G − V0).

   * So the seventh incidence holds iff the median `V0G` is tangent to the circumsphere, iff `A = C`.
   * This holds for every simplex dimension `n ≥ 3`. There `H' = O' + n(G' − O')` is the facet's own analogue of `D`: the face orthocentre for the tetrahedron. The plane is degenerate.
2. **This is the Fano phenomenon in metric form.**
   * In both, six incidences are identities, and the seventh is decided by one quantity that is nonzero over `R`.
   * For Fano's seventh line that quantity is the determinant `2`, which vanishes only in characteristic 2.
   * Here it is the sum of squares `Σ`. It vanishes nontrivially only where `−1` is a sum of at most 8 squares. That happens over `Q(i)`, and over `F_p` for `5 ≤ p ≤ 97` in every sampled case, and even then only with the coincidence `A = C`.
3. **Gluck's theorem with the converse added (THM-4743).**
   * **KNOWN.** Tournaments always have Gluck subsets, since their automorphism groups have odd order.
   * Admitting anti-automorphisms, call a tournament *exceptional* if it then needs a third colour.
   * Every exceptional tournament is strong (1b, from the audit).
   * In the primitive case the exceptions are exactly `C_3`, `Z_5{1,2}` and `P_7`.
   * The Paley tournament `P_q` is exceptional iff `q = 3, 7`.
   * Every `X[C_3^[k]]` is exceptional, an infinite family. Session24 proves that there are no others.
4. **Censuses (FINITE-EXACT).**
   * All tournaments with `n ≤ 9`, and (independent audit) all self-converse tournaments with `n ≤ 11`: the exceptions are exactly `C_3, Z_5{1,2}, P_7, C_3[C_3]`.
   * All circulants on `Z_n`, `n` odd `≤ 15`: exactly the family members.
   * The 110 vertex-transitive 21-tournaments: only `P_7[C_3]`.
   * All Cayley tournaments on `Z_3^3`, `Z_9 × Z_3` and `Z_27`: only `C_3^[3]`.
   * HYP-9353 conjectured that the family is complete, and its four pre-registered predictions passed. Session24's structural proof now establishes completeness.
5. **Merging (§4).**
   * The sign involution: tournament converse, `3x + 1 ↔ 3x − 1`, and the sphere inversion.
   * The seventh incidence: Fano (characteristic), the sphere (positivity), and Collatz cycles (open beyond `m ≤ 91`).
   * Stacking (path-rigid tournaments, THM-2477) is kept apart from triangle substitution (the exceptional family).
   * A numerology guard for 7 and 21.

## 1. What the session23b lane had already done today

* **`tetrahedral_inversion_spheres_20261021`:**
  * the universal claim is false, with a rational counterexample;
  * the four-point obstruction for `A, C, D, E`;
  * the real exceptional locus (`D = E`, `B = F`; five distinct points);
  * the generic sphere `S1 ⊃ {B, C, D, E, F, 4C}` with its inversion invariance;
  * `Q1(H') = 6 − 2 v·s > 0`, i.e. `Σ > 0` over `R`;
  * `A ≠ C` for a genuine tetrahedron;
  * `Q1(A) = μ(16 − T)`;
  * a dimension-free obstruction; `A, C, O, E` concyclic.
* **`fano_tournament_move_quotients_20261021`:**
  * the tetrahedron→Fano bridge (six edges plus a formal point; `F7\ω = M(K4)`);
  * the real non-representability of Fano;
  * the line-flip cube `Q7 → K8`;
  * `P_7`'s 84 regular subsets; the order-42 and order-168 groups have none.
* **`signed_tournament_palette_rigidity_20261021`:**
  * three colours always suffice up to converse;
  * a Gluck subset leaves at most one anti-involution;
  * distinguishing numbers 2, 3, 4 for the order-21, -42, -168 actions.
* **`paley_multicoloring_and_carry_20261021`:** palette laws (not used).

**What is new here:**
* the second sphere `S2` and the pencil with its radical plane;
* `Pow_S2(C) = Σ`, and `Σ = 0 ⟺ A = C` over general fields;
* the field condition;
* the all-dimension form;
* the answer to which tournaments need the third colour: the orbit-pairing lemma, the primitive classification, Paley, the substitution rule, the tower family, and the censuses.

## 2. The seven-point sphere (THM-4742)

* **Labels.**
  * `A`: the second point of line `V0G` on the circumsphere.
  * `B`: the pole of the facet hyperplane.
  * `C = V0`.
  * `D = O + (n+1)(G − O)`.
  * `E = inv(G)`.
  * `F = O + (n+1)(O' − O)`.
  * `H' = O' + n(G' − O')`, the orthocentre for the tetrahedron.
* **Incidences** (generic simplex: `O, V0, G, O'` not coplanar; `ι` = inversion of radius `√(n+1) R` about `O`):

| sphere | labels on it | also on it | defect of the missing labels |
|---|---|---|---|
| `S1` (session23b) | B, C, D, E, F | `ι(C)` | `Pow(H') = Σ`; `Pow(A) = μ(16 − T)` for `n = 3` |
| `S2` (new) | A, B, D, E, F, H' | `ι(A)`, `ι(H')` | `Pow(C) = Σ` |
| `S1 ∩ S2` | the circle Γ through B, D, E, F | | radical plane = plane `O G O'` |

* **Proof shape.**
  * Both spheres have power `(n+1)R²` at `O`. That is forced for any sphere through the inverse pair `D, E` on a ray from `O`.
  * Each sphere is then fixed by three dot products with the basis `v = V0 − O`, `s = Σ_i (Vi − O)`, `w = O' − O`.
  * Every incidence becomes an identity in `v·s`, `|s|²`, `v·w`, `|w|²`, `R²` and `n` ((P1), (P1n)).
* **Fields.**
  * Over an ordered field `Σ > 0`.
  * Over `F_p` (`5 ≤ p ≤ 97`; 50,025 valid configurations, 1,388 with `Σ = 0`), and over `Q(i)` with `V3 = (0, i, i)`: all seven lie on one quadric iff `Σ = 0` iff `A = C`.
  * The audit's exhaustive `F_5` run (368,640 in-regime configurations) agrees.
  * A nondegenerate zero of `Σ` needs `−1` to be a sum of at most 8 squares.
* **Dimensions.**
  * All labels lie in the 3-flat `O + span(v, s, w)`. Inside it the pencil statement holds for every `n ≥ 3` ((P1n) symbolic; (P5) exact for `n = 3…6`).
  * In the plane, `v, s, w` are always dependent, so the statement is void.
  * The first version of this note said the opposite. It tested spheres in `R^n` through `n + 1` labels, and every such system was singular. The audit caught it (MISTAKE-650).
* **The Fano comparison.**
  * The textbook real picture of `PG(2,2)` draws six lines straight and the seventh as a circle. Its determinant is `−2`, so the obstruction is the characteristic.
  * The seven-label sphere draws two spheres where one was claimed. Its obstruction is the positivity of `Σ`.
  * Neither seven-label system is identified with the Fano points (session23b §8 warns against that).
  * The common shape is "one quantity, nonzero over `R`, separates the real geometry from an algebraically consistent incidence".

## 3. Gluck with the converse (THM-4743, HYP-9353)

* **KNOWN.**
  * Tournament automorphism groups have odd order. So Gluck's Corollary 1 gives every tournament a vertex subset with trivial stabilizer, i.e. distinguishing number at most 2 (Godsil 2002, as cited by Lozano; Meslem–Sopena for cyclic tournaments).
  * Session23b adds that three colours always suffice when anti-automorphisms are admitted.
* **Which tournaments need the third colour?** Call them *exceptional*.

| result | statement | type |
|---|---|---|
| orbit pairing | exceptional iff an anti-automorphism fixes every regular `Aut`-orbit on subsets | PROVED |
| strongness | `r(T) ≥ 2` always, so non-strong tournaments are not exceptional | PROVED (audit) |
| primitive `Aut±` | exceptional iff `C_3`, `Z_5{1,2}` or `P_7`: Gluck's Theorem 1 leaves only `S_3, D_10, F_42` among orders `≡ 2 (mod 4)` | PROVED |
| Paley | `P_q` exceptional iff `q ∈ {3, 7}`: `Aut±(P_q) ⊇ AGL(1,q)` is 2-transitive | PROVED |
| substitution | `A` prime, `B` exceptional: `A[B]` exceptional iff `D(Aut±(A)) > r(B)` | PROVED |
| towers | `r(C_3^[k]) = 2`, so every `X[C_3^[k]]` (`X = C_3, Z_5{1,2}, P_7`) is exceptional | PROVED |
| census | exceptions with `n ≤ 9` (all), `n ≤ 11` (self-converse; audit): `C_3, Z_5{1,2}, P_7, C_3[C_3]`; circulants `n ≤ 15`, VT21 and Cayley 27: the family members only | FINITE-EXACT |
| completeness | the family is all of them | PROVED by Session24 (HYP-9353 closed) |

* **Substitution numbers.**
  * `r(B)` is 2, 6, 4, 2, 6 for `C_3`, `Z_5{1,2}`, `P_7`, `C_3[C_3]`, `Z_5{1,2}[C_3]`.
  * `C_3[Z_5{1,2}]` (`n = 15`) has 15,000 regular subsets, while `Z_5{1,2}[C_3]` has none.
* **Paley regular-orbit counts.**
  * The affine group has 2, 4, 36, 3052, 33156, 4618100 regular orbits for `q = 3, 7, 11, 19, 23, 31` with squares, and 0, 0, 12, 1470, 16392, 2306868 with all units.
  * The closed form for prime `q` is `2^q − 2 − q Σ_(d | m, d > 1) (−μ(d))(2^(1+(q−1)/d) − 2)`, where `m` is the order of the centre stabilizer.

## 4. Tournament structure and the open questions: what merges

Each bridge names what transfers and what does not.

* **The sign involution (ANALOGY, with PROVED instances).** Three objects carry an involution that preserves all the "odd" or "algebraic" structure but reverses a sign:
  * the converse of a self-converse tournament: for `P_7`, `x ↦ −x`, which swaps its two Fano planes (session23b);
  * the conjugation `x ↦ −x` taking `3x + 1` to `3x − 1` (THM-4610's remark: every Haar statement survives it, the KNOWN barrier);
  * the inversion `ι` of the sphere pencil, which pairs `D ↔ E` and `B ↔ F`.

  For tournaments:
  * THM-4743 proves that observers (vertex subsets) cannot break the sign in the Gluck circulants and their triangle towers, and settles the primitive case.
  * Session24's proof of HYP-9353 shows that these are the only cases.
  * The pre-audit version of this sentence asserted "only" before any proof existed (MISTAKE-650).

  The Collatz reading is that Haar data are odd-symmetry observers, which cannot break `3x ± 1`. That part is the KNOWN barrier; no new Collatz statement is claimed.
* **The seventh incidence (ANALOGY).** In three places an algebraically consistent incidence fails over the integers or the reals:
  * Fano's seventh line fails because the characteristic is not 2.
  * The sphere's seventh label fails because `Σ > 0`, a positivity obstruction.
  * Nontrivial positive cycles of `3x + 1`. A rational cycle point `c_w/(2^A − 3^l)` exists for every parity word. Cycles with `m ≤ 91` local minima are excluded by linear forms in logarithms plus computation (Steiner; Simons–de Weger; Hercher). The general case is OPEN.

  The base-`p` comparison is in `zp_rank_one_coalescence_20261008.md` (S5): a 1-cycle shape has a single candidate numerator for `p = 2` but `(p−1)^k2` for `p ≥ 3`. What transfers is the shape of the obligation: find the quantity that is nonzero over `R` or `Z`. No mechanism transfers.
* **Stacking versus triangle substitution (PROVED instances; no identification).**
  * THM-2450's grammar is historical, FINITE-EXACT through `n = 10`.
  * Current [THM-2477](../../01-canon/theorems/THM-2477-rigid-strong-classification-through-sixty-two.md) and the [all-order bypass proof](merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md#3-rigid-strong-tournaments-are-completely-classified) prove that path-rigid tournaments (`H = |Aut|`) are exactly the ordered stacks of `{1, C_3}`.
  * THM-4743 with Session24 shows that the exceptional tournaments are exactly `X[C_3[C_3[…]]]`.
  * Both lists start from `C_3`, but one is built by stacking and the other by substitution. The families must not be identified: `C_3[C_3]` is exceptional and has `tc = H/|Aut| = 3159/81 = 39`. The pre-audit text called the exceptional family "the maximally symmetric corner of the same grammar" (MISTAKE-650; caught by Session24's review).
  * A candidate bridge to HYP-9352's presentation counts: compare `tc(T)` with the palette `(q, f)` along substitution. For `P_7`, `tc = 9` and `q = 4`. Session24 counts the Hamiltonian-path orbits fixed by the converse: 3 of 9 for `P_7` and 11 of 39 for `C_3[C_3]`, although every regular subset orbit is fixed. No inequality is claimed.
* **Numerology guard.** The exceptional orders include 7 and 21 (`P_7`, `P_7[C_3]`), and THM-1370's H-spectrum gaps are 7 and 21. No mechanism links them (NUMEROLOGY). S50 studies the spectrum side.
* **No transfer.** As session23b notes, none of this touches LRC or the chromatic number of the plane.

## 5. Next steps

1. **After HYP-9353.** The exceptional family is classified (Session24). Two quantitative questions remain:
   * the palette pair `(q, f)` beyond the `C_3[·]` recursion, e.g. for `P_7[B]` and `Z_5{1,2}[B]`;
   * the relation between converse-fixed subset orbits and converse-fixed Hamiltonian-path orbits.
2. **Geometry:** the `n + 1` vertex choices give `2(n+1)` spheres through `D` and `E`; their configuration is unexplored.
3. **Is the defect always a sum of squares?** For which label systems built from a simplex by centroids, circumcentres, poles, inversions and homotheties is every failing incidence's defect a sum of squares?

## 6. Status

| item | status |
|---|---|
| `S2 ⊃ {A,B,D,E,F,H'}`, pencil, `Σ`, tangency, field dichotomy, all `n ≥ 3` (THM-4742) | PROVED (audited) |
| `S1`, the counterexample, the real exceptional locus, `Σ > 0` over `R` | session23b (credited) |
| `F_p`, `Q(i)` instances; dimensions 3–6 | FINITE-EXACT |
| converse-Gluck: lemma, strongness, primitive classification, Paley, substitution, towers (THM-4743) | PROVED (audited) |
| censuses `n ≤ 9` (all), `n ≤ 11` (self-converse), circulants `n ≤ 15`, VT21, Cayley 27 | FINITE-EXACT |
| completeness of the exceptional family | PROVED by Session24 (HYP-9353 closed) |
| sign involution and seventh-incidence bridges | ANALOGY |
| 7/21 coincidence with THM-1370 | NUMEROLOGY |

## Audit

The independent audit of 2026-10-21 (`audit_S49/`) found that THM-4743 and HYP-9353 had no mathematical error, and that THM-4742 statements 1–4 were correct.

**Corrected:**
* THM-4742 statement 5, the dimension claim, which was false: the theorem holds for every `n ≥ 3`;
* the overclaims in §4 (the "only", Collatz cycles, and the typing of Fano's obstruction);
* the regime (non-coplanarity);
* the field condition (level);
* the credits (session23b's `Σ > 0` and `A ≠ C`; Godsil; Meslem–Sopena's scope);
* the scope of THM-2450;
* an unscripted claim, now removed;
* the INDEX wording.

**Added:** the strongness lemma and the extended censuses. See MISTAKE-650.

**The concurrent Session24 review (cba114930c) independently flagged:**
* the pre-audit "only" sentence;
* the THM-2450 citation, which THM-2477 supersedes;
* the conflation of stacking with triangle substitution.

All are corrected above.
