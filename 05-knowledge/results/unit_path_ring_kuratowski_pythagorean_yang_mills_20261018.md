# Unit-path rings of the six-colour field, Kuratowski and Petersen in the finite planes, the Pythagorean tree, and a Yang–Mills reading of the energy ladder

2026-10-18, mac-mini-2026-10-18-unitpath. This note covers two owner prompts:
- "continue with another similar session, building on the ideas of other agents creatively";
- "keep going as you were and also merge in related to the Yang--Mills-LLM-experiment repository and the ternary tree
  of primitive pythagorean triples and the petersen graph, K_5, K_(3,3)".

Claims are typed as PROVED, KNOWN, CITED, FINITE-EXACT, STRUCTURAL, ANALOGY, COINCIDENCE, NUMEROLOGY, CONDITIONAL or
HEURISTIC. Canon:
- [THM-4700](../../01-canon/theorems/THM-4700-the-unit-path-ring-of-the-six-colour-field-short-reciprocal-walks-and-the-heron-census-of-coefficient-triangles.md):
  the unit-path ring, reciprocal walks, Heron census, mixed quadrilaterals.
- [THM-4701](../../01-canon/theorems/THM-4701-kuratowski-graphs-the-petersen-graph-and-the-finite-unit-circle-graphs.md):
  Kuratowski graphs, the Petersen graph and `G_q`.

Correction: MISTAKE-635. In the brief-2 note, "8 fixed primes per p" is really 8 or 24; the first 24 is at p = 2549.

## 1. What was inherited

- **clock_duality_oct15.**
  - B2 (`cubic_adelic_20261016e.md` section 3): the unit-path ring `Gamma = ⋂ O_P` over conjugation-fixed primes, once
    `O ⊆ Gamma`. Also the overring lemma, and the `t/conj(t)` pole argument.
  - The cubic recipe (`golden_contact_circuit_20261017i.md` section 4), whose recorded next target was "a coefficient
    context that connects these copies".
- **opus S41/S42.**
  - HYP-9289 (additive freeness) and its design target, a "difference-forcing layer".
  - S41 audit a4's exact ideal ladder `E_d` for `d <= 12`.
- **Mine.**
  - THM-4690(5): the 19-adic separation of the copies, valid only for step alphabets integral at the `h = 7` primes.
  - The brief-2 fixed-place census and κ table; THM-4558 (Lemma R); THM-4624 (`G_4 = Clebsch`).
  - The clock-graded ladder of the brief-2 note.

## 2. The unit-path ring of the six-colour field (THM-4700, PROVED unless typed)

**Component.** `Gamma_F = ⋂ (O_F)_P` over the conjugation-fixed primes, which lie exactly over `p = 14 (mod 15)` with
`(p/11) = -1`: 8 per `p`, or 24 when the cubic splits completely mod `p`.
- The proof needs only `zeta3` on top of clock_duality's B2 for `F' = Q(zeta5, sqrt-11, w)`. Every generator of `O_F`
  was nevertheless re-exhibited as an explicit unit sum.
- Consequences:
  - `1/n` lies in `Gamma_F` iff `n` avoids those primes. All `1/p` with `p < 29` are reachable; `1/29` is not.
  - Every coordinate factor of the cubic recipe (`A + aM`, `C + bM`, `C + cM`) has power-basis denominators only at
    2, 3 and 19, so the fourteen coefficient copies lie in the zero component of `G(F)`.
  - THM-4690(5) is therefore alphabet-relative. `sqrt-3` splits the fixed 19-primes of `F'`, and
    `delta_19 = (1 - 5 sqrt-3)/(1 + 5 sqrt-3)` has a pole at one prime of the pair and a zero at the other.

**Short walks** (exact; `delta = conj(g)/g`, `delta + delta^-1 = t^2/n - 2`):

| Steps | Points |
|---:|---|
| 2 | `1/2`, `1/5`, `1/13` |
| 4 | `1/3`, `1/4`, `1/7`, `1/19` (the 19-denominators of b, c) |
| 7 | `1/23` |
| 64 | `1/17` |
| 82 | `1/11` (clock_duality's connector) |

- **Two-step reciprocals.** These are the leg-one solutions `(2n)^2 - D y^2 = 1` for the four imaginary quadratic
  subfields `Q(sqrt-D)` of `F`, `D ∈ {3, 11, 15, 55}`. That gives `n = 2, 5, 13, 122, 181, 1985, 2521, 7562, ...`.
  `D = 55` never contributes, because its `x` is always odd.

**Heron criterion.**
- A triangle with totally positive squared sides `R_i` has vertices in the CM field `F = F+(sqrt-3)` iff
  `H = 2 Σ R_iR_j - Σ R_i^2` (= `16·area^2`) lies in `3(F+)^{*2}`.
- **Census of coefficient distances** `{1, d_a, d_b, d_c}`, where `d_k = 1/|k|`:
  - only the four equilateral triangles and the two golden ones exist (PROVED by construction);
  - the other fourteen are FINITE-EXACT exclusions, each with an explicit non-residue certificate at a degree-one place.
- **Consequences for the copies.**
  - No translation matching joins two unit-coefficient copies on the `b`- or `c`-axis.
  - Single reflection links always exist: `q' = (conj(k)/k) conj(q) (kq - s)/(conj(kq - s))`.
  - Mixed quadrilaterals (two unit sides, two `d_k` sides) always exist. The conic `sigma^2 - N t^2 = (N - 1)/3` has
    quaternion algebra `(N, -3)`, which `k = A + B sqrt-3` splits visibly; the point is `t = (1 + A)/(3B)`,
    `sigma = (A + N)/(3B)`.
  - A quadrilateral with three `d_k` sides is excluded for `b` and `c` by the house bound, `max|sigma(b)| = 3.70` and
    `max|sigma(c)| = 4.12`.
- **Reading for the six-colour design** (STRUCTURAL). Connecting copies is cheap; rigidity is not. Any coefficient
  context that couples fibre shifts on the `b`- or `c`-axis needs mixed cycles of length at least 4. This sharpens
  HYP-9289's difference-forcing target and refutes nothing.

**Walk-length geometry.**
- `L <= 1` is a norm condition, and `L <= 2` is one square class (Heron).
- `L <= 3` is the existence of a point on an explicit genus-one quartic over the degree-24 field `F+`
  (THM-4700 remarks). So three-step copy connectors are a rational-point problem. This is why finite direction screens
  (THM-4690's 450 directions) are inconclusive.
- Known lower bound: connector length `>= 3` (certified). The conjugate-modulus bound reaches 6 to 7 for the farthest
  copy pairs.

## 3. Kuratowski, Petersen and Clebsch in the finite planes (THM-4701, PROVED + FINITE-EXACT)

In `G_q = Cay(F_(q^2), mu_(q+1))`, the residue graph of Lemma R:
- **Two circles meet in at most two points** (the norm form is anisotropic). So no `G_q` contains `K_(2,3)`, and hence
  none contains either Kuratowski graph.
- **`K_4` occurs iff `q` is an odd power of 2.** The two apexes over a unit edge are at squared distance
  `N(sqrt-3) = 3`, and `3 = 1` only in characteristic 2. In odd characteristic the planar obstruction to `K_4`
  survives exactly.
- **Petersen in the golden residue plane.** `G_4` is the Clebsch graph (THM-4624). Its distance-two layer is the ten
  two-sums of the residue pentagon `mu_5`, and it induces the Petersen graph as the Kneser graph `K(5,2)`.
- **The golden Petersen realization.**
  - Outer pentagon of circumradius `R`, inner pentagram of circumradius `r` turned by 90 degrees, with
    `R^2 + r^2 = 1`, `R/r = phi`, `R r = 1/sqrt5`: a golden Pythagorean triple.
  - It lives in `Q(zeta20)`. At the prime above 2 it reduces exactly onto that Clebsch layer: `R -> w`, `r -> w^2`, so
    the image is `w mu_5 ∪ w^2 mu_5 = F_16^× ∖ mu_5`.
  - At `q = 19 (mod 20)` it reduces into `G_q`.
  - In characteristic 0 the Petersen graph is the unit-sum graph of the two-sums, since a disjoint pair sums to minus
    the fifth root.
- **Census.** For primes `q <= 113`, the Petersen graph embeds in `G_q` exactly for
  `{17, 19} ∪ {41..113} ∖ {73}`.
  - The counts through a fixed directed edge are irregular: 480, 768, 240, 3840, 480 at `q = 17, 19, 41, 43, 47`.
  - The real realization space is a smooth 2-parameter family (rigidity rank 15). So the gaps must come from Frobenius
    permuting its geometric components. Which components there are, and their fields, is OPEN.

**Connection to HYP-9289** (STRUCTURAL).
- The ports' golden pentagon has two-distance graph `K_5`. That is the 465-pattern naive block.
- Modulo 2 the same pentagon becomes the independent set `mu_5` of the Clebsch graph, whose two-sums carry the
  Petersen graph.
- So the `K_5` and the Petersen graph are two shadows of one pentagon. No colouring consequence is claimed.

## 4. The ternary tree of primitive Pythagorean triples (PROVED / KNOWN / FINITE-EXACT)

- **B2 for the rational plane** (KNOWN in substance). In `Q(i) = Q^2` the conjugation-fixed primes are `(1 + i)` and
  the inert `p = 3 (mod 4)`. So `Gamma_(Q(i)) = Z[i][1/p : p = 1 (mod 4)]`.
  - These denominators are exactly the hypotenuses of Berggren's tree (Fermat).
  - The residue colouring at the ramified prime `(1 + i)` is Woodall's `chi(Q^2) = 2`.
- **Walk lengths in the rational plane.**
  - `L(1/n) = 2` never happens, because `4n^2 - 1` is never a square.
  - The symmetric three-step walk `1/c = delta + conj(delta) ∓ 1` exists iff a leg is `(c ± 1)/2`. These are the
    near-60-degree triples, governed by `X^2 - 3Y^2 = 1` with `X = (3c ∓ 1)/2`.
  - They are the alternating branch `C, A, C, A, ...` of the tree: `(3,4,5) -> (15,8,17) -> (33,56,65) ->
    (209,120,241) -> (451,780,901) -> ...`.
  - Up to depth 9 this branch contains every near-60-degree triple (FINITE-EXACT). The Pell bijection is PROVED.
- **Branches as Pell sequences.** The classical `B, B, ...` branch (legs differing by 1) is the `sqrt2` Pell branch; the
  `C, A, ...` branch is the `sqrt3` one.
  - The six-colour field's two-step reciprocals are the `sqrt(-D)` Pell branches for `D = 3, 11, 15`.
  - Its four-step reciprocals are Eisenstein (60/120-degree) triangles: `1/7` from `(3, 5, 7)`, `1/19` from
    `(5, 16, 19)`, and the two-step `1/13` from `(7, 13, 15)` with its 60-degree angle.
  - STRUCTURAL: the same norm-form arithmetic organizes all three.
- **Collatz contrast** (ANALOGY). Berggren's tree is the inverse of a descent with a strict Lyapunov function, the
  hypotenuse. That is why every primitive triple reaches `(3,4,5)`. Collatz's inverse tree has no such function, and
  the ROOT energy programme is an attempt to replace it with an `l^2` (energy) bound.

## 5. A Yang–Mills reading (sources: the repository README only)

The repository `ethanwsussman/Yang--Mills-LLM-experiment` (owner link; treated as trusted) collects Ethan Sussman and
Levent Alpoge's attempt to have Claude complete Balaban's renormalization-group construction of Yang–Mills:
- QED3 and YM3 first, aiming at YM4;
- a self-contained "Theorem 0" on the continuum limit along Balaban's RG;
- an "RP trick" (reflection positivity) to win O(4) invariance.

The substantive files are 10 MB `.tex` sources that GitHub will not render, plus PDFs. Fetching either is a download,
which needs the owner's explicit permission, so only the README's ideas are merged here.

**(a) RG steps and the energy ladder** (ANALOGY, with exact data).
- In THM-4667's ideal ladder, depth `d` plays the role of the RG step, and the ideal `N -> infinity` model plays the
  continuum limit.
- S41's exact cross terms `E_d - E_(d-1)` for `d <= 12` decay geometrically at first, ratio about 0.54 for `d <= 4`.
  That is like a super-renormalizable flow ("UV friendlier", as the README says of QED3/YM3).
- They then cross over to a power law (`ladder_decay_exponents.out`). On `d = 6..12` a power fit, exponent 1.91 with
  rms 0.051, beats a geometric fit, ratio 0.805 with rms 0.094. The local exponents are 2.4, 2.4, 1.9, 1.8, 1.1, 1.7.
- So there is no "gap": the depth correlations decay polynomially, as in a marginal flow.
- **Positivity is the brake.** Every cross term is negative and `E_d >= 0`, so the tail after `d = 12` totals at most
  `E_12 = 1.479`. The local exponent cannot settle at or below 1 (PROVED, elementary).
- The clock-graded split (section 6) shows what the decay is made of. A non-decaying same-clock (twin) dissipation,
  about -0.10, is cancelled by growing positive odd-lag terms. The small cross term is a difference of two
  non-decaying parts, which is the familiar shape of a delicate multi-scale cancellation (HEURISTIC).

**(b) Reflections give rotations** (ANALOGY and STRUCTURAL).
- The RP trick wins continuum O(4) invariance from discrete reflection symmetry. In the unit-distance setting:
  - THM-4700(5c)'s reflection links generate new unit coefficients by reflecting `kq` across lines;
  - Berggren's tree is the orbit of `(3,4,5)` under a free monoid inside the integral Lorentz group `O(2,1; Z)`, which
    is a reflection group (Barning, Hall: classical, CITED from memory);
  - the field planes have dense rotation groups `U_F`, so "rotation invariance" is algebraic there.
- Nothing quantitative is transferred.

**(c) The additive six-colour problem as a hard-constraint `Z/6` gauge problem** (STRUCTURAL reformulation).
- HYP-9289's additive colourings `c = Σ c_i(x_i)` are 0-cochains on the coefficient factors. Their coboundaries are the
  link variables `c_i(x + kq) - c_i(x)`.
- Axis edges demand nonzero links. A contact family demands that the sum of its links over the family's axes, a
  plaquette-like flux, is nonzero.
- "Difference forcing" asks that every nonzero flux value be forced in each direction. The language suggests treating
  port extension as a flux problem; nothing is claimed beyond the reformulation.

**What would sharpen this.** The statement of Theorem 0 and the RP rotation argument (`YM_Rotation_Invariance_simplified.pdf`).
Read with the owner's download permission, they could be compared with THM-4700's reflection mechanism and with the
ladder's crossover.

## 6. Collatz: the clock-graded ladder to depth 7 or 8 (FINITE-EXACT, ideal model)

`clock_graded_ladder.py 8 64` reproduces `E_1`-`E_7` to 12 digits (S41 a4). New rows, in THM-4667 units:

| `d` | Cross term | Lag 0 (twins) | Lag 1 | Lag 2 | Lag 3 |
|---:|---:|---:|---:|---:|---:|
| 6 | -0.085089 | -0.102282 | +0.019058 | -0.004400 | +0.004011 |
| 7 | -0.059204 | -0.101833 | +0.040117 | +0.003866 | -0.000596 |
| 8 | see below | | | | |

The twin part is nearly constant from `d = 6`, and the brief-2 parity alternation fades: -0.1023, -0.1018. Lag 1 keeps
growing, and lag 2 turns positive at `d = 7`. The prediction for `d = 8` was registered before the run printed it
(`ladder_d8_prereg_20261018.json`):
- P1: the twin part lies in `[-0.106, -0.098]`;
- P2: lag 1 is positive and above +0.040;
- P3: the other lags sum to between +0.055 and +0.063.

## 7. Ledger

| Claim | Type |
|---|---|
| `Gamma_F` = intersection over fixed primes; census 8/24; reciprocal criterion; walks; Pell branches | PROVED |
| Heron criterion; six realizable coefficient triangles; reflection links; mixed quadrilaterals; house bound | PROVED |
| Fourteen excluded coefficient triangles; no 1- or 2-step copy connector | FINITE-EXACT (certified non-residues) |
| `G_q`: at most 2 common neighbours, `K_4` criterion, Clebsch/Petersen/Kneser, golden realization and its reductions | PROVED |
| Petersen census `q <= 113`; Pythagorean branch to depth 9; ladder rows `d = 7`, 8 | FINITE-EXACT |
| B2 for `Q(i)`; Woodall; Clebsch second subconstituent; Berggren as a reflection-group orbit | KNOWN / CITED |
| `K_5` vs Petersen as two shadows of the pentagon; Pell branches across the three fields; gauge reformulation | STRUCTURAL |
| RG-step reading of the ladder; RP / reflection reading; Collatz vs Berggren descent | ANALOGY |
| The cancellation shape of the ladder; persistence of the twin term beyond `d = 8` | HEURISTIC |

## 8. Next obligations

- **Three-step copy connectors.** Points on the genus-one quartics of THM-4700 for `x = b(m' - m) + s`. Local
  solubility at the places above 2, 3 and 19 would either kill them (proving connector length at least 4) or license a
  search.
- **Petersen census.** The geometric components of the Petersen realization variety, their fields of definition, and
  the resulting Chebotarev law. This includes why `q = 73` fails and why 43 has 3840 copies.
- **Ladder.** Prove that the twin part has a nonzero limit, or compute `d = 9` with a faster (C) clock-graded
  recursion.
- **Yang–Mills.** With permission, read Theorem 0 and the rotation-invariance PDF, and test whether reflection
  positivity has a usable shadow in the residue-graph Hoffman/Delsarte bounds.

## 9. Files

- `04-computation/experiments/unit_path_ring_20261018/`:
  - `explicit_inverses.py`, `unit_path_ring_six_colour.py`, `matching_criterion.py`, `heron_census.py`,
    `mixed_quadrilaterals.py`, `house_bound.py` (each with its `.out`);
  - `two_step_reciprocals.out`, `ladder_decay_exponents.py` (+ `.out`), `clock_graded_ladder_d8.out`,
    `ladder_d8_prereg_20261018.json`.
- `04-computation/experiments/kuratowski_petersen_20261018/`:
  - `residue_graph_kuratowski.py`, `petersen_census.py`, `petersen_zeta20.py`, `petersen_mod2_and_unit_sum.py`,
    `pythagorean_tree_walks.py` (each with its `.out`);
  - `petersen_census_count47.out`, `petersen_census_exist113.out`.
- `04-computation/experiments/heegner_kappa_20261018/`: the session-start κ(43) tabu probe (negative).
