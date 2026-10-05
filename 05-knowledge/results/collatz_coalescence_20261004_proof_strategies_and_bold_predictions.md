# Collatz through the whole repository: the strategy atlas (what each proof style proves, the exact wall where it stops, the barrier it cannot see), five axes along which the thread's structures coalesce (three places and the `<2,3>`-solenoid; criticality; clocks as `H^1`; the Artin coordinate; the sheet as a phase), twelve bold predictions with their cheapest tests, and four probes run today (the `x2x3` lonely spectrum below `1/14` is a quadratic-residue phenomenon, first value `5/73`; the powers of three can shadow `-5` but never `-1`, and `log_3(-5) in Z_2` is a universal resisting exponent with share `0.1597`)

**Session:** opus, `collatz-synthesis-20261004` (opus-2026-10-04-S1), 2026-10-04.
**Owner's directive:** "digest possible proof strategies for Collatz in this
repo, along with many other related and seemingly unrelated ideas, in search
of holistic synthesis and bold predictions about how structures coalesce
across domains; capture the spirit of our many novel ideas and conjecture in
that spirit, at its heart; maximum freedom in creativity."
**Status: SYNTHESIS (every inherited claim carries its original label and
ID + path) + PROVED (Proposition E, elementary; the one-line identity behind
HYP-9174; the three-place form of the clock-as-fixed-point-group reading) +
FINITE-EXACT (Probes A, B, E, F: THM-4522's 88-point census reproduced; the
value spectrum of `I` to `q <= 4000` and the prime index census to
`20000`; the carry census to `A <= 18`; `|Bad_k cap <3>|` exactly to
`k = 24` and by a dynamic programme to `k = 6000`) + CONJECTURE / SPECULATIVE
/ ANALOGY where so marked. Collatz OPEN; LRC(14) OPEN in the repo (external
claim under audit). Nothing here is a proof step for either.** New
hypothesis files: HYP-9174 (the Artin coordinate of the `x2x3` spectrum),
HYP-9175 (the powers of three against bounded-depth certificates).
Script `04-computation/experiments/collatz_coalescence_20261004_probes.py`,
output beside it (`ALL PROBES RAN in 13.1s` before the E2 extension; `~60s`
with it). Independent audit: subagent, blind re-derivation of Probes A and
E and a read of this note's typed claims (section 9).

---

## 0. The answer in one screen

1. **Every proof style in the thread is typed by one of two walls, and the
   two walls are the two halves of the conjecture.** Density, mixing,
   thin-divergence, certificate, Fourier and rank methods see the DRIFT and
   stop at the pointwise statement (the `0.95`-dimensional no-descent set
   `E_inf`, the `0.050` bits per step deficit). Diophantine, necklace,
   clock, LRC-type and rigidity methods see the SHEET and INTEGRAL walls and
   stop at bounded complexity (runs `<= 91`, shapes `A <= 22`, `n < 2^71`).
   No method in the repository or the literature sees both walls at
   unbounded complexity, and the barrier atlas's reading that this is the
   whole difficulty survives everything added since 2026-09-22 (section 2).
2. **The structures coalesce along five axes** (section 4). (I) Three
   places: the multiplier `3/2^v` is a unit of `Z[1/6]` with
   `|.|_oo |.|_2 |.|_3 = 1`, real growth is the excess of 3-adic contraction
   over 2-adic expansion, and the clock groups `Z/(2^A - 3^p)` are the
   fixed-point groups of the dual automorphisms of the three-place solenoid
   `(R x Q_2 x Q_3)/Z[1/6]`, exactly as the Fibonacci groups
   `F(2,j)^ab = Z[phi]/(phi^j - 1)` are the fixed-point groups of the cat
   map (Probe B): the golden thread is the rank-one model of a rank-two
   problem. (II) Criticality: Collatz sits exactly on the boundary of every
   classical criterion (Apery `1`, Ridout `2`, Borel-Dwork `1`, Terras
   `h* < 1` by `0.050`, Parseval `3^(-1/2)` against `log_2 3 - 1` by
   `1.3%`, the strategy-cube window `[0,1]`), and the slack is in every case
   a quantity the conjecture itself controls. (III) Clocks are `H^1`: a cycle
   is a vanishing carry class in `H^1(Z; Z[1/6]_u) = Z/(2^A - 3^p)`, a fair
   necklace split is the isotypic decomposition of that group under `Z/j`
   (THM-4515), and the quantifier dictionary of the LRC/JC thread extends:
   JC wants universal vanishing, LRC and the Collatz cycle half want
   universal nonvanishing outside an explicit finite list. (IV) The Artin
   coordinate: the discrete logarithm of `2` (and of `3`) in finite cyclic
   groups is the one coordinate behind the circulant level operator
   (THM-4520), backward-tree rigidity (THM-4523), the two-sheet clocks and
   Paley (THM-4532), the `x2x3` lonely spectrum (THM-4522; today: HYP-9174)
   and the powers of three (today: HYP-9175). (V) The sheet is a phase: the
   sheet involution is the half-turn of `<2> = (Z/3^n)^x` (THM-4521) and a
   coset swap of `<3>` in `(Z/2^n)^x` (Probe F); on the multiplicative
   spectrum it is the parity of the character, and the seed-`1` density is
   the cancellation to `O(1)` of two sums of size `0.16 (3/2)^n` (section
   4.5). Everything conjugation-invariant is sheet-blind.
3. **Twelve predictions** (section 5), typed and priced. Two were tested
   today and survived with their numbers (P3 = HYP-9174, P4 = HYP-9175);
   one is a one-line theorem (Proposition E: a cycle point is 2-adically
   approachable by powers of three iff it is `1` or `3 mod 8`); the rest
   name their cheapest decisive test.
4. **The heart, in one sentence.** Every solved relative of Collatz is a
   rank-one unit dynamics (the doubling map, the golden beta-map, `F_2[x]`,
   Mersenne, Gilbreath's additive automaton, the one-generator circulants);
   Collatz is the first rank-two one, `<2,3>` with `log_2 3` irrational, read
   at three places at once, and each coalescence found in the thread is a
   rank-one shadow of it: `<2>` cyclic modulo `3^n` gives the circulant and
   the Fourier profile, `<3>` of index two modulo `2^n` gives the sheet
   coset, `phi` gives the golden words, `2` gives Mersenne. The missing
   mechanism, in this language, is a statement about the joint
   `(log_2, log_3)` coordinate of one integer along its own orbit, with the
   archimedean sign attached; the cycle half already has one (Baker is
   rank two at the real place), the divergence half has none.

---

## 1. Inheritance pass and the concept board

**Closest proved mechanism.** THM-4476/4487/4495/4499 (thin divergence, the
dip spectrum, the Spitzer identity, the moving-barrier ballot bound): the
one all-orbits mechanism the thread owns, and its floor `o(X^(h*))` is the
method's floor (THM-4487, THM-4506). On the cycle side: THM-4484 (free and
sporadic cycles), THM-4515/4516 (fair splits, perfect-power clocks),
the eleventh note's Proposition 1 (the clock as a torsion group).
**Canonical hostile.** The `3x-1` sheet (`{-1}, {-5,-7,-10}, {-17,...}`) and
`5x+1` (THM-4474's drift control); for today's probes, `a = 7` in THM-4523's
Artin condition and `q = 73` in THM-4522's spectrum.
**Corrected near miss.** MISTAKE-553/554/555/556/559/560 (my own S15 notes):
coincidences promoted to structure, bounds read in the wrong direction,
classical results under-attributed. Today every "coalescence" carries a
typed map and a lost coordinate, and every classical name was checked
against the source list of the barrier atlas and CORE-PAPERS.
**Least-used sidecar.** The multiplicative lonely runner (THM-4522): its
`x2x3` spectrum had been computed and left as a discrete top; its mechanism
(Probe A) was not named.

**Concept board (seven live objects, compared after every probe).**

| object / representation | predicate | invariant or extremal | operation applied today | lost coordinate | cheapest next test |
|---|---|---|---|---|---|
| the `<2,3>`-solenoid `(R x Q_2 x Q_3)/Z[1/6]` and its dual lattice | cycles = carry classes vanishing in `Fix(u-hat)` | `|2^A - 3^p|`; Gersonides `1` | three-place product formula; the cat-map model (Probe B) | which fixed points carries can hit | 3x+k census vs the torsion-torsor model (P2) |
| the no-descent set `E_inf` and `Bad_k` | `Z^+ cap E_inf = empty` | dim `h* = 0.94996`; `W_k ~ 2^(hk) k^(-3/2)` | intersect with `<3>` (Probe E) | the sign | share to `k = 30` exactly (HYP-9175) |
| the 3-adic Syracuse law `mu_n` and its two spectra | mixing (2.3); H1 | Parseval `3^(-1/2)`; resonant rate `3^(h*-1)` | even/odd character split (4.5) | the sheet (phase) | fixed-frequency rate at `u = 5` (P6) |
| the Artin coordinate: `<2>` in `(Z/3^n)^x`, `<3>` in `(Z/2^n)^x`, `<2,3>` in `(Z/q)^x` | which statements are residue-expressible | cyclic / index 2 / index distribution | Probes A, E, F | sign and drift | HYP-9174 tests (i)-(iii) |
| the clock group as `H^1` and its isotypic split | nonvanishing of the carry class | `Phi_e(2^m, q^x)` factors (THM-4515) | the LRC/JC quantifier dictionary (4.3) | the integer `x` (size) | Schmidt-subspace reading of the cycle vector (D56) |
| the golden rank-one model: cat map, beta-map, `F(2,j)` | `2Theta(n) in Z[phi]` (THM-4528) | `L_j`, `|L_j - 1 - (-1)^j|` | Probe B | the second prime | negative prediction P9 |
| the coverage grammars (THM-4474/4479; codex compilers of 2026-10-03/04) | residual at reading depth `k` | sharp exponent `1 - h` | the `3^a` family through `log_3(-5)` | first-hit vs legal prefix (MISTAKES 2026-10-03) | compiler on `a = 11 mod 128` (P10) |

---

## 2. The strategy atlas

Each row: the style, what it has PROVED in the repository (ID + path or note),
the barrier it overcomes / is blind to (typing rule of the
[barrier atlas](collatz_procgen_20260922_barrier_atlas.md) section 0), and
the exact wall where it stops. "Wall" means the last true statement, not a
vague "it gets hard".

| # | style | proved in the repo (selection) | sees / blind | the exact wall |
|---|---|---|---|---|
| S1 | **Density and mixing** (Terras, Korec, Tao, GGM) | the entropy curve `E(gamma) = h(max(1/2, gamma/log_2 3))` with Terras at `gamma = 1` and Korec at `log_4 3` ([THM-4487](../../01-canon/theorems/THM-4487-dip-spectrum-entropy-curve-and-sharpness-of-thin-divergence.md)); the Mazur digest's Theorems A-C ([note](mazur_positive_density_20260928.md)): harmonic mass `3^n mu_n`, the negative cycles as resonances, the seed-1 test | DRIFT / SHEET, DEFECT, INTEGRAL, DIM | the measure-zero, dimension-`0.95` set `E_inf`; log-density only (Tao); `lim H_n(1) > 0` OPEN |
| S2 | **Thin divergence and no-descent counting** | every non-periodic orbit has `o(X^(h*))` points below `X`, reciprocal sums converge ([THM-4476](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md), [THM-4499](../../01-canon/theorems/THM-4499-thin-divergence-is-little-o-of-x-to-the-h-star.md)); `W_k = Theta(2^(hk) k^(-3/2))` by a Spitzer identity ([THM-4495](../../01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md)); landing multiplicity worst case `0.631(k-D)` ([THM-4506](../../01-canon/theorems/THM-4506-landing-multiplicity-exact-worst-case-and-recursion-saturation.md)) | DRIFT / SHEET (the bound holds on both sheets, which have cycles) | `N(X) = o(X^(h*))` is the method's floor (THM-4487); only an orbit-coupled local-time bound can lower it (THM-4506); `5x+1` out (`1/log_2 5 < 1/2`) |
| S3 | **Certificates and coverage** (residue descent, strategy cube, reroute compilers) | provability iff every parity-graph cycle has odd density `< log_3 2`, Collatz undecided at every level with fraction `0.435` ([THM-4474](../../01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md)); the price of provable descent `-> 0` at the sharp rate `2^(-(1-h)k)` ([THM-4475](../../01-canon/theorems/THM-4475-price-of-provable-descent-tends-to-zero.md), [THM-4478](../../01-canon/theorems/THM-4478-critical-tube-affine-capacity-sharp-provability-price.md), [THM-4479](../../01-canon/theorems/THM-4479-strategy-cube-distance-to-provability-sharp.md)); coefficient-descent cylinders ([THM-4512](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md)); codex's early-reroute, checkpoint and join compilers with the `3^a` hostile family ([coverage note](collatz_coverage_memory_and_quadratics_20261004.md)) | DIM only in relaxations (Applegate-Lagarias, Caraiani, `F_2[x]`) / SHEET and DRIFT | the base cycles cohere to every level ([THM-4521](../../01-canon/theorems/THM-4521-collatz-base-cycles-cohere-to-every-level-half-turn-intertwining-and-exact-energy-transfer.md)(1)): no residue certificate of either kind at any level; the residual at depth `k` always contains the shadow classes of the negative cycles (today: `log_3(-5)` for the powers of three) |
| S4 | **Diophantine cycle exclusion** (Steiner, Simons-de Weger, Hercher; Baker) | free vs sporadic cycles and Gersonides ([THM-4484](../../01-canon/theorems/THM-4484-free-and-sporadic-cycles.md)); fair splits are circulant, `2^K - q^X = prod Phi_e(2^m, q^x)` ([THM-4515](../../01-canon/theorems/THM-4515-fair-consecutive-splits-of-cycle-necklaces-are-circulant.md)); perfect-power clocks are Fermat-Catalan ([THM-4516](../../01-canon/theorems/THM-4516-perfect-power-clocks-are-fermat-catalan-identities.md)); the per-period finite check ([LRC-relations note](collatz_lrc_relations_20260930.md)); the carry census to `A <= 22` (eleventh note) | SHEET, INTEGRAL, DEFECT / DRIFT | bounded complexity: runs `<= 91` (Hercher), shapes `A <= 22`; the cycle half is of lonely-runner type (a finite check per parameter, no uniform argument) |
| S5 | **Fourier and spectral** (the 3-adic Syracuse law) | the primitive profile `M(h)` is one sequence, (2.3) forces super-polynomial decay, the maxima sit on `+-2^s` and decay at the no-descent rate ([five mirrors](collatz_five_mirrors_20260929.md)); the frequency-one reduction `c_J = 2^-s mu_hat_J(1) W` ([THM-4519](../../01-canon/theorems/THM-4519-frequency-one-coefficient-carries-the-renewal-series.md)); the level operator is a Gauss-twisted circulant with spectrum on the half circle ([THM-4520](../../01-canon/theorems/THM-4520-collatz-level-operator-is-a-gauss-twisted-circulant-with-spectrum-on-the-half-circle.md)); the character spectrum identities ([three mirrors](collatz_three_mirrors_20260929.md) Theorem 1) | DRIFT (via Tao) / SHEET (conjugation-invariant), INTEGRAL, DIM | H1 ([HYP-9166](../hypotheses/HYP-9166-h1-frequency-one-coefficient-decays-below-critical.md)): a `x2x3` square-root cancellation for one character sum; the Parseval rate is a non-normal transient (THM-4520), so no spectral-radius argument applies |
| S6 | **Rigidity and graph structure** | `Aut(N,T) = 1`, backward trees identify vertices prime to `3`, Eckmann-Hilton point `-1/5` ([THM-4523](../../01-canon/theorems/THM-4523-collatz-functional-graph-is-rigid-backward-trees-identify-vertices.md)); Collatz iff the trivial component is everything (Theorem C of [the uniqueness note](collatz_functional_uniqueness_20261001.md)); `Gamma_1` dense in `Z_2 x Z_3` ([connectivity note](collatz_connectivity_from_rigidity_20261001.md)) | neither (rigidity holds for `3n-1`, 3 components, and `5n+1`) | component membership is archimedean; "every acyclic finite view of every integer occurs in `Gamma_1`" (so no finite view decides) |
| S7 | **Conjugacies and carriers** | Bernstein-Lagarias `Phi`, Monks-Yazinski `Omega` (CITED, atlas); golden `Theta`: Collatz iff `2Theta(n) in Z[phi]` ([THM-4528](../../01-canon/theorems/THM-4528-collatz-parity-in-base-phi-golden-beta-map-and-holonomy.md)); the Hadamard identity `(1-3z)(F ⊙ G) = n + zF` and the weighted-mediant generation of every rational cycle ([pentagon note](collatz_pentagon_fixed_points_argument_styles_20260930.md)); the base-`3/2` tree ([posets note](collatz_posets_dags_zeta5_20260927.md)); the Mahler bridge ([THM-4469](../../01-canon/theorems/THM-4469-mahler-bridge-adjacent-block-pairs.md)) | none (reformulations move every barrier into the target) | Bell-Lagarias: Collatz iff the basin series is rational iff D-finite iff continues across an arc; the one Mahler equation Collatz supplies mixes roots of unity (non-rigid class); `F_n` D-finite iff rational iff periodic (Bezivin) -- no holonomic carrier |
| S8 | **Diophantine criticality** (Apery, Ridout, Borel-Dwork; the zeta(5) shape) | the orbit's own approximants have 2-adic exponent `nu_L = d_L/(L log_2 3)` with `nu_L - 1 = log_2(n C_L/m_L)/(L log_2 3)`; three-place product `eta_L/3^L`; `R_oo R_2 = 1` ([posets note](collatz_posets_dags_zeta5_20260927.md) Props 4-10) | none (diagnostic) | every slack is the descent depth, the shadow error or the valuation oscillation: "the word does not know the integer" made sharp |
| S9 | **Potentials, ranks, Lyapunov functions** | ranks are tensions, defect `log(3/2)` at `-1` ([THM-4482](../../01-canon/theorems/THM-4482-ranks-are-tensions-strategy-cube.md)); forced charges on backward trees ([THM-4483](../../01-canon/theorems/THM-4483-forced-charges-two-place-ranks.md)); no finite polynomial-valuation correction ([THM-4507](../../01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md)); no prefix rank (S13 Prop 4); the price sheet of a cycle shadow, `c_w` bits per bit, size price `K <= log_2(m + |x|)` ([pricing note](collatz_pricing_approaches_20260927.md) Theorem 1) | DIM (they say what a rank must carry) / SHEET, DRIFT | a rank must use the whole integer; the only such rank known is the stopping time; a potential pricing future regenerations is the height walk again |
| S10 | **Relaxations with choice** | the choice ladder, E-SCC `Q2` to `10^18`, `Q1` mirror, sibling dimension ladder `dim = h(log_3 2)` for every `3n+b`, positive measure for `q >= 5` ([procgen synthesis](collatz_procgen_20260922_synthesis.md)) | DIM / SHEET (the hostile points of every relaxation accumulate at the minus sheet's cycles) | "removing DIMENSION by choice exposes SHEET" (atlas section 2.4): the `-1 mod 2^j` class and the negative rationals at `-1` |
| S11 | **Undecidability and logical shape** | Kurtz-Simon `Pi^0_2`-completeness of the generalised family (CITED); completeness of tournaments is `Pi_1` while Collatz is `Pi_2` ([Camion/Busch note](camion_busch_gaps_polyhedra_collatz_20261002.md)) | n/a (binds only uniform complete proof systems) | a proof must use special properties of `3x+1` (GGM's reading); per-parameter statements are decidable |
| S12 | **Tournaments, clocks, polyhedra, motifs** | the four-vertex reading is time reversal, not the sheets ([THM-4472](../../01-canon/theorems/THM-4472-four-vertex-reading-of-3n-plus-minus-1.md)); 3-/5-node motif censuses sheet-blind (THM-4521(4)); the owner's label map and the drop multiplicities `2 + N(6d+1)` ([THM-4527](../../01-canon/theorems/THM-4527-syracuse-in-owner-labels-microcosm-and-difference-spectrum.md), [THM-4530](../../01-canon/theorems/THM-4530-syracuse-drop-multiplicities-two-copies-over-z-densities-and-pair-injectivity.md)); the Syracuse clock mod `9` and two-sheet clocks ([THM-4532](../../01-canon/theorems/THM-4532-tournament-clocks-two-sheet-clocks-all-odd-tournaments-and-the-syracuse-clock.md)); the `{1,2,4} = QR_7` coincidence explained (MISTAKE-554) | none / SHEET (every residue-defined invariant is identical for `3x+-1` after `x -> -x`, THM-4474(E)) | "drop statistics cannot constrain cycles" (THM-4530); the kinship is with the cycle half's `Pi_1` finite checks |
| S13 | **Refuted external proofs** | Kawasaki's fixed-point proof ([THM-4471](../../01-canon/theorems/THM-4471-kawasaki-fixed-point-collatz-proof-refuted.md)); the `zeta(5,7,9)` gain-sign error (S15 batch); Opfer's withdrawn Mahler route (atlas) | -- | each fails a control: `x -> x + 1`, `3n-1`, or the sign of one integral |

**What the atlas says that is new since 2026-09-22.** The procgen barrier
atlas typed the literature; the sessions since then added S2 (the method's
exact floor), S5 (the profile and the circulant), S6 (rigidity), S8 (the
criticalities) and S12 (the owner's combinatorial readings), and every one
of them lands on one of the two walls. The one statement in the whole
thread that is pointwise and not obviously of the forbidden prefix class is
the shell-revisit / orbit-coupled local-time statement (S9/S10 notes,
THM-4506), because it constrains the orbit of one integer through its own
address. That is where the honest frontier of the divergence half is.

---

## 3. The anatomy: two halves, three places, one deficit

| the half | its named target | the places it needs | the finite check | the wall |
|---|---|---|---|---|
| cycles (`forall w != (2)^k: x_w notin Z_>0`) | an `S`-unit relation with one unknown integer (THM-4490's specialisation) | `2` (the word), `3` (the carry is a unit), `oo` (size, sign) | per period / per run count (Steiner, SdW, Hercher; the repo's `A <= 22`) | uniformity in the parameter; carries come within `1/D` of multiples of the clock (pentagon note 3.2), so no metric bound |
| divergence (Lagarias's Periodicity Conjecture restricted to `Z_>0`) | `Z^+ cap E_inf = empty` (S13) | all three, with no parameter | none (Conway) | the pointwise statement; sheet-blind by necessity (PC is about all `3x+k`) |

**The exact slacks (section 3 of the posets note, S2 of this note).**

| criterion | critical value | the Collatz object | the slack = what the conjecture controls |
|---|---|---|---|
| Terras count vs integers | `2^(h*k)` vs `2^k`, `h* = 0.94996` | `Bad_k`, `W_k` | `0.050` bits per step (`1 - h*`): THM-4478/4479/4495, the localglobal note's information budget |
| 2-adic Apery exponent | `1` | `nu_L = d_L/(L log_2 3)` | `log_2(n C_L/m_L)`, the descent depth |
| Roth-Ridout three places | `2` | `r_L = 2^(d_L) m_L/3^L` against `(xi, 0, oo)` | `eta_L (n C_L)^2`, the shadow error (HYP-9164) |
| Borel-Dwork radii | `R_oo R_2 = 1` | `F_n(z) = sum 2^(d_k) z^k` | `2^(liminf - limsup) d_k/k`, the valuation oscillation; on a divergent orbit both limits are `log_2 3` and the slack is exactly `1` |
| Fourier (frequency one) | Parseval `3^(-1/2) = 0.5774` vs critical `log_2 3 - 1 = 0.5850` | `mu_hat_n(1)` | the margin `1.3%` (H1); measured rate `0.569` |
| Fourier (resonant maximum) | the no-descent rate `3^(h*-1) = 0.9465` | `M(h)` on `+-2^s` | `M(h) ~ c P_h` CONJECTURAL (five mirrors 2c) |
| strategy cube | window `[0,1]` | the base cycles `0`, `10`, `1` | maximal indecision at every level (THM-4474, THM-4521(1)) |
| drift (Tao/GGM) | `q < p^(p/(p-1))`, i.e. `3 < 4` | `log_2 3 = 1.585 < 2` | `5x+1` fails; nothing finer |
| Gersonides | `|2^a - 3^b| = 1` | the four free cycles `1, -1, {-5,-7}` and `0` | the sporadic `-17` at the clock `139 = 3^7 - 2^11` (THM-4484) |
| the lonely-runner side | `kappa(2,3) = 1/5` | the `x2x3` runners | a mod-5 obstruction; the spectrum is Artin (HYP-9174) |

Hostile control against over-reading criticality: THM-4494 (`C* <= 197/125 <
log_2 3` for AMM 12592) shows a "critical coincidence" that was a proof
artifact; the rows above are not coincidences, each is an identity or a
theorem with the slack named.

---

## 4. Five axes along which the structures coalesce

Each axis: source, target, map, preserved predicate, lost coordinate,
sidecar, cheapest decisive test (the typing discipline of
RESEARCH-PROTOCOL section 3).

### 4.1 Axis I: three places and the `<2,3>`-solenoid; the golden rank-one model

**The objects.** The Syracuse step multiplies by `3/2^v`, a unit of
`Z[1/6]`, with `|3/2^v|_oo = 3/2^v`, `|3/2^v|_2 = 2^v`, `|3/2^v|_3 = 1/3`,
product one (pricing note, Proposition 2). Over `L` odd steps real growth
is exactly `3^L/2^(d_L) = (3^(-L))^(-1) (2^(d_L))^(-1)`: the excess of 3-adic
contraction over 2-adic expansion. The 2-adic coordinate carries the word
(Terras, Bernstein-Lagarias: the dynamics is a shift), the 3-adic
coordinate carries the address (`m_l = 2^(-d_l) S_(l-1) mod 3^l`, the
crossings note), and the archimedean coordinate carries the size, the sign
(the sheet) and the drift.

**The fixed-point-group reading, at three places (PROVED, one line).** Let
`Sigma = (R x Q_2 x Q_3)/Z[1/6]`, the compact dual of the discrete group
`Z[1/6]`, and for a shape `(A,p)` let `u = 3^p/2^A in Z[1/6]^x` act on
`Sigma` through the dual of multiplication by `u`. Then
`Fix(u-hat) = dual of Z[1/6]/(u - 1)Z[1/6] = dual of Z[1/6]/(3^p - 2^A) ≅ Z/|2^A - 3^p|`,
because the clock is prime to `6`. The eleventh note's Proposition 1
([monotile/discrepancy note](collatz_lucas_monotile_discrepancy_20260930.md))
proved this on the 2-adic solenoid `Sigma_2 = dual(Z[1/2])`; the same proof
holds for the three-place solenoid, and it is the right space because the
sheet and the drift live on the factor `R`. The carry class `S_w mod
(2^A - 3^p)` is a point of `Fix(u-hat)`, and an integer cycle of the shape
exists iff some word's carry is the zero point (the Bohm-Sontacchi condition,
THM-4484). Probe B re-checks the census to `A <= 18`: the zero class is hit
exactly at `1` (shape `(2,1)` and its repeats), `-1` (`(1,1)`), `-7, -5`
(`(3,2)`), the seven points of the `-17` cycle (`(11,7)`), and nothing else.

**The golden model.** For the Fibonacci cat map `A = [[1,1],[1,0]]` on
`T^2 = R^2/Z^2`, `Fix(A^j) ≅ Z^2/(A^j - I)Z^2 ≅ Z[phi]/(phi^j - 1)`, of order
`|det(A^j - I)| = |L_j - 1 - (-1)^j| = |F(2,j)^ab|`, the Spectre's monodromy
tower (the owner's `phi^8 + 1 = 7 phi^4`, `phi^10 = 1 + 11 phi^5` are its
levels `4, 5`, eleventh note); while the golden shift (the parity language of
the standard Collatz map, THM-4528) has `tr(A^j) = L_j` points of period
dividing `j`. Probe B checks both to `j = 60`: the cat map loses exactly the
shift's two points `0^inf` and `(01)^inf` when `j` is even. So: **the
Fibonacci groups are to the golden beta-map what the clock groups are to
Collatz** -- the fixed-point groups of the natural extension, and the
"lattice cycles" (`(1/2)Z[phi]` for the beta-map, exactly two, THM-4528; `Z`
for Collatz, conjecturally four) are the periodic points that land on a
dense line. The golden thread is the rank-one model (one unit `phi`, two
real places); Collatz is rank two (two units `2, 3`, three places), which is
why every "golden" coincidence in the thread (the Pillai clocks `5, 11, 29`;
the Sturmian thin word; the AMM golden thresholds
[THM-3017](../../01-canon/theorems/THM-3017-golden-threshold-for-the-checkpoint-closure-capacity-criterion.md)/[THM-3027](../../01-canon/theorems/THM-3027-capacity-threshold-is-log-sqrt5-phi.md))
stops at a finite census (THM-4520's note: sixteen Lucas/Fibonacci gaps,
all with `K <= 8`).

Typing: source = Collatz cycle equation; target = periodic points of a
compact-group automorphism; map = Pontryagin duality on `Z[1/6]`; preserved
= the group `Z/|clock|` and the rotation symmetry; lost = which classes the
carries hit (the content) and the integer `x` itself; sidecar = the carry
cocycle (Axis III); decisive test = none beyond the census (this axis is a
dictionary, EXACT, not a theorem about Collatz).

### 4.2 Axis II: criticality everywhere, and the shape shared with transcendence

The table of section 3 is the evidence. The reading: Collatz is not a hard
instance of some criterion, it is the *boundary case* of every criterion the
thread has tried, and the amount by which each criterion misses is a
quantity equivalent to the conjecture (descent depth, shadow error,
valuation oscillation, no-descent probability). This is why "the word does
not know the integer" (posets note 3.6) is the sharp statement and why
pure counting cannot finish: counting is exactly critical.

The same shape appears in my exp-integral thread: `Phi_Q(z) = int_0^1
e^(zQ(t)) dt` is an `E`-function in the parameter
([S4 reformulation](../../07-reflections/exp-integral-is-an-E-function-in-the-parameter-general-reformulation-and-rigidity-opus-S4.md)),
its transcendence at `z = 1` is proved in degree `2`
([S4 theorem](../../07-reflections/degree-2-exp-integral-is-transcendental-rigorous-via-linear-beukers-shidlovskii-opus-S4.md))
and for pure powers in every degree
([S4 theorem](../../07-reflections/pure-power-exp-integral-transcendental-all-degrees-via-1F1-one-over-d-opus-S4.md)),
and the general case reduces to a non-resonance of critical values
([S4 reduction](../../07-reflections/degree-3-exp-integral-rigidity-is-critical-value-nonresonance-computational-opus-S4.md)):
Siegel-Shidlovskii delivers everything except a lower bound at one critical
point, which the problem's own structure makes exactly critical. H1 (the
frequency-one coefficient) is the Collatz instance of the same shape: the
Gauss-sum spectrum is flat (THM-4519/4520), Parseval gives the typical
coefficient, and the needed statement is a lower-order gain at one frequency
that the recursion's own transients put on the edge. The zeta(5) record is a
third instance: the public proofs are 2-adic Apery-type approximants with
exponent above `1` and arithmetic holonomy bounds, and on a Collatz orbit
the Apery exponent is exactly `1` (section 3). Typing: ANALOGY (the shape of
the hard direction), with the lost coordinate named: the exp-integral has a
linear differential equation and Collatz's `F_n` has none (Proposition 11
of the posets note, Bezivin), so nothing transfers as a method. What
transfers is the diagnosis: **the missing input is a non-vanishing at a
critical point, never a better average.**

### 4.3 Axis III: clocks are `H^1`; the quantifier dictionary

With the twelfth note's cocycle (`S_(uv) = 3^(p_v) S_u + 2^(A_u) S_v`,
`D_(uv) = 2^(A_u) D_v + 3^(p_v) D_u`; the weighted-mediant law generates
every rational cycle, pentagon note 3.2), the clock group is
`H^1(Z; Z[1/6]_u) = Z[1/6]/(u-1) = Z/(2^A - 3^p)` (the group cohomology of
the cyclic group generated by one period acting by the multiplier), the
carry `S_w` is a 1-cocycle, and a cycle is a coboundary (twelfth note:
"cycles = coboundaries"). Fair necklace splits (THM-4515) are the isotypic
decomposition of this group under the rotation `Z/j`: the DFT diagonalises
the cycle equation and the clock factors as `prod_(e|j) Phi_e(2^m, q^x)` --
the same cyclotomic factorisation as `Res(t^j - 1, t^2 - t - 1) = prod_(e|j)
Res(Phi_e, t^2 - t - 1)` for the Fibonacci groups (Probe B's identity read
factor by factor). Zsigmondy (eleventh note, Proposition 3) is the
statement that repeated words multiply the carry by the cyclotomic cofactor.

**The dictionary.** The D5 thread (HYP-9031, THM-2473, THM-2542) found that
JC(2) wants universal `H^1`-vanishing (the abelianised class is a Kummer
class) while LRC(14) wants universal nonvanishing (the LRC(14) frontier is
written in Cech `H^1` language, THM-2542), with anti-parallel quantifiers.
Collatz splits across the dictionary: the cycle half is a universal
nonvanishing of the carry class, with an explicit finite list of exceptions
(the Gersonides free words and `-17`), exactly LRC's shape (the sixteenth
note's table: the exceptional objects satisfy short relations, Beck-Everett's
harmful relations have a parity layer, the Collatz clock has `(-1)^A mod 6`
and an odd carry), while the divergence half has no clock, no group, no
parameter: `H^1` of the free monoid is not finite, and that is the exact
sense in which it is "global". Typing: EXACT for the identifications of the
groups; ANALOGY for the transfer of methods (Beck-Everett's window-product
Fourier argument needs a product over independent factors, while the Collatz
word sum is a transfer-matrix product over levels, D54).

### 4.4 Axis IV: the Artin coordinate

The discrete logarithm of `2` (or `3`, or the pair) in a finite cyclic group
is one coordinate that the thread keeps rediscovering:

| where | the group | the fact used | what it gives |
|---|---|---|---|
| THM-4520 | `(Z/3^n)^x = <2>`, cyclic of order `L_n = 2 3^(n-1)` | `2` generates | the level operator is a circulant; the Fourier profile is one sequence `m_n(k) = mu_hat_n(2^k)`; the half-turn `2^(L/2) = -1` is the sheet (THM-4521) |
| THM-4523 | `(Z/a)^x` | `2` must generate for backward separation (`a = 7` fails: half the units bare) | backward separation; A-injectivity for every non-Wieferich `a` (the audit's sketch, MISTAKE-553 entry) |
| THM-4532 | `(Z/p)^x`, `(Z/9)^x` | discrete-log coordinates; the sheet is the parity of the exponent | Paley minus a vertex is a two-sheet clock; `log_2 S(A) = 2(A mod 3) - v mod 6` |
| THM-4522 + HYP-9174 (today) | `(Z/q)^x` and its subgroup `<2,3>` | `I(m/q) = (least absolute residue of the coset)/q` | the `x2x3` lonely spectrum: top `{1/5, ..., 1/14}` where `2` is a primitive root, then `5/73` where `<2,3>` is the squares and `5` the least non-residue; index distribution and the two-generator Artin fraction `0.707` (Probe A) |
| HYP-9175 (today) | `(Z/2^k)^x` and its subgroup `<3> = {1, 3 mod 8}` | `-1 notin <3>`; closure of `<3>` in `Z_2^x` | which cycle points the powers of three can shadow; the resisting exponent `log_3(-5) in Z_2` (Probe E) |

Typing of the axis: source = each of the five problems; target = the
cyclic group with its generator; map = discrete logarithm; preserved = every
residue-expressible predicate; lost = the archimedean sign and the drift
(the two things residue methods never see: THM-4474(E), the barrier atlas);
sidecar = positivity; decisive test = whether a proposed argument changes
when `3x+1` is replaced by `3x-1` (if it does not, it is residue-only and
cannot finish). **The owner's "local versus global" thesis, in this
coordinate:** the local problem is the Artin coordinate of `<2,3>` at one
modulus (and is trivial or Artin-type, Probe A), the global problem is the
joint coordinate of one integer at all moduli with the sign attached.

### 4.5 Axis V: the sheet is a phase; the observer

**Three exact facts.** (a) On the 3-adic side the sheet involution is inside
the cyclic group: `u -> -u` is `k -> k + L/2` on `<2>` and conjugates the
Fourier coefficients, `mu_hat_n(-u) = conj mu_hat_n(u)` (THM-4521(2); Probe
F). (b) On the 2-adic side it is outside: `-1 notin <3>` in `(Z/2^n)^x`,
`n >= 3`, so the sheet swaps the cosets `3^Z` and `-3^Z` (Probe F). (c) On
the multiplicative spectrum it is the parity of the character:
`psi_j(-1) = (-1)^j`, and Theorem 1(ii) of the three-mirrors note gives
`E_n + O_n = rho_n(1) - rho_(n-1)(1)` and `E_n - O_n = rho_n(-1) - rho_(n-1)(-1)`
for `E_n = sum_(psi even, prim) E psi(Y_n)` and `O_n` the odd sum.

**The exhibit.** From the S19 data (`rho_n(-1) = 0.9748 (3/2)^n` at `n = 18`,
increasing; `rho_n(1) = 0.637, 0.425, 0.3155, 0.2736` at `n = 6, 10, 14, 18`):
`E_n ≈ -O_n ≈ (0.975/6)(3/2)^n`, i.e. `E_18 ≈ +240`, `O_18 ≈ -240`, while
`E_18 + O_18 = rho_18(1) - rho_17(1)` is of size `10^(-2)`. **The seed-`1`
density is the cancellation, to four digits at level `18`, of a
sheet-symmetric and a sheet-antisymmetric sum each of size `0.16 (3/2)^n`.**
Everything conjugation-invariant (H1, the mixing estimate (2.3), the
Parseval mass, the motif censuses, the profile `M(h)`) sees `E_n` only, i.e.
the `-1` resonance; all the information about the positive integers is in
`O_n`. This is the Fourier-side form of the barrier atlas's "every
mechanism that overcomes DRIFT is blind to SHEET". (DERIVED from Theorem
1(ii), PROVED, and the OBSERVED spike law; the numbers are S19's.)

**The observer.** THM-4523 is an observer statement: every vertex prime to
`3` is identified by what it sees looking backward, and the only point where
the two backward branches coincide (where an observer could not tell the two
copies apart) is `-1/5`, the fixed point of `6x + 1`, a 2-adically odd
3-adic number that no rational graph contains. The owner's "two copies" is
exact in THM-4530: every drop `d` is realised `2 + N(6d+1)` times, the two
unit preimages `8d + 1` and `-4d - 1` being one on each sheet. The three
appearances of `5` -- `kappa(2,3) = 1/5` (THM-4522), the Eckmann-Hilton point
`-1/5`, the Monks-Yazinski conjugate `Omega{-5,-7,-10} = {4/5, 2/5, 1/5}` --
have two mechanisms, not one: `5 = 2 3 - 1` for the second and third (the
clock of the one-step mixed map, and `3 = 2^(-1) mod 5` which collapses the
`<2,3>`-orbit of `1/5` to one `<2>`-orbit), and `5 = ` the least prime
outside `{2,3}` for the first; typed EXPLAINED COINCIDENCE (the lesson of
MISTAKE-554 applied before writing).

---

## 5. Bold predictions, typed and priced

Each: the prediction; type; the structure it asserts; the cheapest decisive
test; status today.

**P1 (the three-place divergence proof). SPECULATIVE.** The divergence half
will be proved, if at all, by an argument that reads the archimedean
coordinate and both adic coordinates of one orbit together, and it will
prove Lagarias's Periodicity Conjecture for every `3x+k` at once, because
every correct divergence argument is sheet-blind (the atlas; PC is a
statement about all sheets). Its margin will be the `0.050` bits per step
deficit, and its form will be an orbit-coupled local-time bound (THM-4506's
conclusion) rather than any count. Diagnostic (cheap): any proposed
divergence argument that changes when `3x+1` is replaced by `3x-1` is wrong;
any that uses only residues is incomplete (Axis IV).

**P2 (the torsion-torsor model of rational cycles). MODEL, testable.** The
number of integer cycles of `3x+k` (`k = +-1 mod 6`) of shape `(A,p)` is
distributed as the number of zero hits of `C(A-1,p-1)/p`-many near-uniform
carry classes in `Z/|2^A - 3^p|` (a Poisson law with mean the eleventh
note's sum, restricted to the convergent shapes), and the sporadic `-17` is
the one hit the model owes `3x+1`. Test: THM-4484's free/sporadic census for
`k <= 1000` against the Poisson means per shape; a systematic excess at the
convergent shapes would refute the model and locate a structure (Lagarias
1990's `k^(1-eps)` counts are the thing to reproduce first).

**P3 = HYP-9174 (the Artin spectrum). CONJECTURE; tested today.** Below
`1/14` the `x2x3` lonely spectrum is the Artin coordinate of `<2,3>`: values
`n_q/q` with `n_q` the least absolute residue of a proper coset, first `5/73`,
then `1/17, 1/19, 5/97, 13/259, ...`; the index-2 primes `q = 1 mod 24` give
exactly the least quadratic non-residue over `q` (186 of 186 to `20000`); the
two-generator Artin fraction is `0.707 +- 0.010` against the heuristic
`0.6975`; only accumulation point `0`. Survived (Probe A).

**P4 = HYP-9175 (the powers of three). PROVED + FINITE-EXACT + CONJECTURE;
tested today.** A cycle point is a 2-adic limit of powers of three iff it is
`1` or `3 mod 8` (Proposition E, section 6), so `3^a` shadows `-5, -7, -37,
-55, -61` and never `-1, -17, -25, -41, -91`; the class `a = log_3(-5) mod
2^(k-2)` (`1, 3, 3, 11, 11, 11, 11, 11, 267, ...`) resists every depth-`k`
certificate, and codex's residual `a = 3 mod 8` of 2026-10-04 is its first
three digits; the share of `Bad_k` inside `<3>` converges to `0.1597`,
strictly below the plain tilt `1 - log_3 2 = 0.369`, by ballot conditioning.
Survived (Probe E).

**P5 (the seed-`1` density is a sheet cancellation). DERIVED + OBSERVED.**
`E_n = -O_n + O(1)` with both of size `(0.975/6)(3/2)^n`; the limit
`lim_n sum_(m <= n) (E_m + O_m) = lim rho_n(1)` exists and is positive iff
Mazur's seed-1 test passes (S19 Theorem C: `limsup n^(1/6) H_n(1) > 0` is
forced; the limit is OPEN). Prediction: the first sheet-sensitive quantity
ever bounded from below in this thread will be an odd-character sum `O_n`,
and no bound on `|mu_hat|`, on the mixing distances or on the profile can
replace it. Test: the character-spectrum output of S23 already contains
`E_n` and `O_n` to level 16; print them and the cancellation ratio.

**P6 (fixed frequency against moving frequency). CONJECTURE.** The decay
rate of `|mu_hat_n(u)|` for a *fixed* unit `u` (`u = 1, 5, 7, 11, ...`) is
one universal constant `rho_1 = 0.569 +- 0.003 < 3^(-1/2)`, while the
resonant maximum `M(h)` is a moving-frequency phenomenon (`+-2^s`,
`s ≈ h log_2 3 - 6`) at the no-descent rate `3^(h*-1)`; the gap between
`rho_1` and the Parseval rate is a genuine exponential deficit, not a
transient (the rescaled data `3^(n/2)|mu_hat_n(1)|` fall from `9e-2` to
`3e-15` over `n = 200..2500`, i.e. `0.9866` per level). Test: mac-mini's
rescaled recursion (`collatz_h1_20260929_*`) run at the frequency
`u = 5 = 2^(log_2 5)` to `n = 2500`; equality of block rates with those of
`u = 1` confirms, a different rate refutes.

**P7 (the `L^p` criticality of the Syracuse law). CONJECTURE.** The 3-adic
Syracuse law is absolutely continuous (conditional on (2.3), Mazur digest)
with density in `L^p` exactly for `p < 2`: `E[rho_n^2]` grows linearly
(`0.31` per level, S19), which is `P(rho > t) ~ c t^(-2)` with `c ≈ 0.38`,
fed by the forward closures of the negative-cycle spikes. Test: `E[rho_n^p]`
for `p = 1.5, 1.9, 2.1` at levels `<= 18` from the existing FFT law;
boundedness for `1.9`, linear growth at `2`, exponential at `2.1`.

**P8 (the two uniformities). SPECULATIVE.** If arXiv:2609.02604 (LRC for 14
and 15 runners; audit in progress, OPEN-QUESTIONS) stands, every LRC(`k`)
will fall to the same two-layer architecture (a sieve local at every prime
plus one archimedean bound), and the Collatz cycle half with `<= m` runs is
the same architecture with the places swapped (a 2-adic sieve on words plus
Baker at the real place; Simons-de Weger, Hercher). Neither architecture
extends to the uniform statement, and the uniform statements need the same
kind of new input: a relation-lattice covering-radius bound for LRC and a
Schmidt-subspace-type statement on the cycle vector
`(2^(d_0), ..., 2^(d_(p-1)), x 2^A, x 3^p)` for Collatz (D56). No cheap
test; the prediction is about the shape of the next theorem on either side.

**P9 (a negative golden prediction). FALSIFIABLE.** No golden constant will
enter the conclusion of a Collatz theorem except through the parity language
(the golden-mean shift, THM-4528) and its periodic-point counts `L_j`; every
other golden appearance (the Pillai clocks `5, 11, 29`, the AMM thresholds,
the Spectre tower, Zeckendorf/Wythoff seeds) is the rank-one model of Axis I
showing through, and stops at a finite census. A single Collatz theorem with
`phi` in a rate or a constant that is not `L_j`-combinatorial refutes this.

**P10 (what the coverage residual is). CONJECTURE, test for the codex
lanes.** The residual of every reroute grammar whose rules read a bounded
number `k` of steps of the source is, up to the common-future joins, the
union of the shadow classes of the rational cycle points admissible to the
family: for the powers of three, exactly the `s_k W_k` classes of HYP-9175
(`51280` of `2^22` exponent classes at `k = 24`). Test: run the checkpoint
and early-reroute compilers of 2026-10-04 on `a = 11 mod 128` and `a = 267
mod 512` and count the certified members; a count above `0` with reading
depth `<= 7` resp. `<= 9` would show that the joins read deeper than their
declared depth (and would be worth a theorem).

**P11 (Collatz is Rudolph's world). SPECULATIVE-STRUCTURAL.** The
absolute continuity of the Syracuse law, the Poisson law of the deep 3-adic
zero lines `u 2^Q = -+1 mod 3^d` (three mirrors, section 3), Furstenberg's
theorem behind THM-4522(S) and H1's "`x2x3` square-root cancellation" are
four faces of the measure rigidity of the joint `x2, x3` action (Rudolph,
Johnson; effective versions by Bourgain-Lindenstrauss-Michel-Venkatesh, via
Z. Wang, CITED in THM-4522): the Syracuse law is the stationary measure of
the `x3` contraction with `x2^(-1)`-geometric weights, and the only
mechanism that can produce the exponential gain of H1 beyond Parseval is an
entropy input of that kind. Test: none cheap; the statement that would make
it precise is an identification of `mu` as a `<2,3>`-quasi-invariant measure
on the 3-adic solenoid with positive entropy, and that is the first thing to
try to write down.

**P12 (the rank-two principle, the heart). SPECULATIVE, falsifiable in the
negative.** Every relative of Collatz that is solved is rank one (one unit:
the doubling map, the beta-map, `F_2[x]`'s `x+1` and `x`, Mersenne,
Gilbreath's additive automaton, the one-generator circulants), or is rank
two with the ratio of logarithms rational (then the clock degenerates and the
dynamics is a one-generator system in disguise); Collatz is the first rank-two
system with `log_2 3` irrational that is asked a question at three places. The
prediction: **no argument whose only arithmetic input is rank one at the adic
places (Terras's bijection, 3-adic mixing of `<2>`, residue certificates,
the circulant, the profile) will prove divergence, and the first genuinely new
divergence result will be a rank-two adic statement about one orbit: a bound
on the joint `(2-adic word, 3-adic address)` coincidences
`2^(d_l) m_l = 3^l n + S_(l-1)` along a single orbit (the orbit-level Poisson
law for the ridge seeds of S22/S23), which is pointwise and therefore not in
the class excluded by Proposition 4 of the directions note.** Negative test:
exhibit a published or in-house divergence statement for `3x+1` that is not
sheet-blind and uses only one adic place -- none is known (the atlas).

---

## 6. The probes (FINITE-EXACT) and Proposition E (PROVED)

All in `collatz_coalescence_20261004_probes.py`; output beside it.

**Probe A (the Artin coordinate of the `x2x3` spectrum).** Positive control:
over all reduced `m/q` with `q <= 600`, exactly `88` points have `I >= 1/14`,
with denominators `5, 7, 10, 11, 13, 14, 26, 28, 33, 52, 56` and values
`1/5 (4), 1/7 (6), 1/10 (4), 1/11 (20: q = 11, 33), 1/13 (30: q = 13, 26, 52),
1/14 (24: q = 14, 28, 56)` -- THM-4522 Theorem S reproduced. The identity
`I(m/q) = n_q(m)/q` for `q` prime to `6` is one line (the orbit is a coset);
for `q = 2^e 3^f q'` the script takes the minimum over the intermediate
denominators. Values below `1/14`, moduli prime to `6` up to `4000`:
`5/73 (index 2), 1/17, 1/19, 5/97 (2), 13/259 (3), 7/145 (2), 29/601 (8),
11/247 (3), 31/697 (8), 19/431 (10), 1/23 (2; -1 is a non-residue), 1/25,
1/29, 13/385 (2), 1/31, 7/241 (2), 1/35, 1/37, 13/485 (8), 23/865 (4)`.
Primes to `20000`: index distribution `{1: 1598, 2: 465, 3: 86, 4: 42, 5: 16,
6: 19, 7: 5, 8: 13, 9: 3, 10: 4, 12: 3, 14: 1, 15: 1, 16: 1, 17: 1, 24: 1,
56: 1}`, two-generator Artin fraction `1598/2260 = 0.7071`; for the 186
primes `q = 1 mod 24` of index exactly `2`, `q I_max(q)` is the least
quadratic non-residue in every case; the largest `I_max` over primes with a
proper subgroup: `5/73, 5/97, 29/601, 19/431, 1/23, 197/6563 (index 17),
7/241, 5/193, 11/439, 109/4513 (12), 29/1201, 157/6553 (56)`. Why `73`: it is
the least prime `1 mod 24` (both `2` and `3` squares, `-1` a square), and
`5` is the least non-residue; no smaller modulus prime to `6` has a proper
`<2,3>` whose extremal coset misses `+-1` (`23` and `47` have index `2` but
`-1` in the other coset).

**Probe B (clocks as fixed-point counts).** `|det(A^j - I)| = |L_j - 1 -
(-1)^j|` and `tr(A^j) = L_j` for `j <= 60`; the carry census of every shape
`(A,p)`, `A <= 18`, hits the zero class of the clock only at the known
cycles (`1`, `-1`, `-7/-5`, the `-17` cycle at `(11,7)` with its seven
rotations) and their repeats; `242461 = 2^18 - 3^9` is hit only by the
repeat of `1`.

**Probe E (the powers of three) and Proposition E.**

*Proposition E (PROVED).* `Z_2^x = {+-1} x (1 + 4Z_2)`, and `-3`, being
`5 mod 8`, topologically generates `1 + 4Z_2`. So the closure of
`{3^a : a in Z}` in `Z_2^x` is `{(-1)^a (-3)^a : a in Z_2} = (1 + 8Z_2) cup
(3 + 8Z_2)`: a 2-adic unit is a limit of powers of three iff it is `1` or
`3 mod 8`. Among the integer cycle points of `3x+1` this admits `-5, -7`
(`3, 1`), `-37, -55, -61` (`3, 1, 3`) and `1`, and excludes `-1, -17, -25,
-41` (`7`) and `-91` (`5`). For `x` admissible, `log_3(x) in Z_2` is the
unique 2-adic integer with `3^(log_3 x) = x`, and `3^a = x mod 2^k` for every
`a = log_3(x) mod 2^(k-2)`; such `3^a` follow the word of `x` for `k` steps,
so for `x = -5` (word `(110)^inf`, in `E_inf`) they have no coefficient
descent within `k` steps. ∎ Digits: `log_3(-5) = 1 mod 2, 3 mod 4, 3 mod 8,
11 mod 16, ..., 11 mod 128, 267 mod 512, 1291 mod 1024, 3339 mod 2048, ...,
68011179275 mod 2^37`.

*Counts.* `Bad_k subset {3, 7 mod 8}` (the second letter of a no-descent word
is `1`), and `Bad_k cap <3> = Bad_k cap {3 mod 8}`. Exact residue scans:
`|Bad_k| = 3, 4, 8, 13, 19, 38, 64, 128, 226, 367, 734, 1295, 2114, 4228, 7495,
14990, 27328, 46611, 93222, 168807, 286581` for `k = 4..24` (THM-4479's
table extended), and `|Bad_k cap <3>| = 1, 1, 2, 3, 4, 8, 13, 26, 45, 71, 142,
247, 395, 790, 1386, 2772, 5019, 8468, 16936, 30485, 51280`; shares `0.333 ->
0.1789`. The `(j, o)` dynamic programme (exact condition, float probabilities;
it reproduces `51280/2^24` at `k = 24`) gives shares `0.16983, 0.16504,
0.16254, 0.16119, 0.16047, 0.16010, 0.15991, 0.15982` at `k = 50, 100, 200,
400, 800, 1600, 3200, 6000` and `1/k`-extrapolations `0.15975, 0.15973,
0.15972, 0.15972`. The limit is a ballot-conditioned constant: under the
critical tilt a third letter `0` has probability `1 - log_3 2 = 0.369`, but the
survivors are conditioned to stay positive, which favours the higher start
`111` over `110`; with `h` the harmonic function of the killed tilted walk
the limit is `(1-r) h(0.170) / ((1-r) h(0.170) + r h(1.755))`, `r = log_3 2`,
and the data give `h(0.170)/h(1.755) = 0.3249`.

**Probe F (half-turn and coset swap).** `2^(L/2) = -1 mod 3^n` for
`2 <= n <= 20`; `-1 notin <3> mod 2^n` for `3 <= n <= 20`.

---

## 7. Verdicts

| claim | status |
|---|---|
| strategy atlas rows S1-S13 | SYNTHESIS of PROVED/CITED/FINITE-EXACT results, each with ID + path; no status changed |
| the two-wall reading survives everything since 2026-09-22 | SYNTHESIS (the honest frontier is the orbit-coupled local-time statement) |
| three-place fixed-point-group reading of the clock | PROVED (one line; the eleventh note's Proposition 1 at three places) |
| cat-map / golden-shift / Fibonacci-group identities | PROVED (classical; Probe B to `j = 60`) |
| `I(m/q) = n_q(m)/q` | PROVED (one line) |
| 88-point census; the value spectrum to `4000`; prime index census to `20000`; least-non-residue law on 186 primes | FINITE-EXACT |
| two-generator Artin fraction `0.707` vs `0.6975` | FINITE-EXACT vs HEURISTIC (within one standard error) |
| Proposition E (closure of `<3>`; shadows; the resisting exponent) | PROVED |
| the share limit `0.1597` | FINITE-EXACT to `k = 24`, dynamic programme to `6000`, extrapolated; the harmonic-function formula CONJECTURED |
| the seed-`1` cancellation `E_n ≈ -O_n` | DERIVED from S23 Theorem 1(ii) (PROVED) and S19's spike law (OBSERVED) |
| P1, P8, P11, P12 | SPECULATIVE (typed, with the negative tests named) |
| P2, P6, P7, P10 | CONJECTURE with a cheap test named |
| P9 | negative prediction, FALSIFIABLE |
| the exp-integral / zeta(5) / Collatz shape | ANALOGY (lost coordinate: no differential equation for `F_n`) |
| Collatz, LRC(14) | OPEN |

**What changes for the repo.** Two hypothesis files (HYP-9174, HYP-9175);
the clock-as-`H^1`-as-fixed-points dictionary is now stated at three places;
the `x2x3` spectrum has a mechanism; the `3^a` hostile family of the
coverage lane has its residual explained (the `-5` shadow through
`log_3(-5)`) and a concrete test (P10); the Fourier side has a stated
reason why conjugation-invariant statements cannot see seed `1`. A
META-PATTERNS candidate card is proposed at the end of this note (not
promoted here: it has four threads of evidence but one session's authorship).

---

## 8. Directions and obligations

* **D62.** HYP-9174 tests (i)-(iii); in particular the GRH-conditional
  two-generator Artin density (Pappalardi-type) against `0.707`.
* **D63.** HYP-9175: `h` numerically (Lemma M of THM-4499), the formula
  against `0.15972`; the residue scan to `k = 30`; the compiler test P10.
* **D64.** P5/P6/P7 from the existing S19/S23 outputs and mac-mini's
  rescaled recursion (three cheap runs, one afternoon).
* **D65.** P2: the torsion-torsor model against THM-4484's census for
  `k <= 1000`.
* **D66.** P11: write the Syracuse law as a `<2,3>`-quasi-invariant measure
  on the 3-adic solenoid and compute its entropy; this is the one direction
  here that could turn a SPECULATIVE prediction into a theorem statement.
* **Obligations.** Independent audit of Probes A and E and of Proposition E
  (subagent, section 9); the owner's three-place reading of Axis I should be
  cross-checked against the eleventh note's resultant column; nothing here
  changes any canon status.

**META-PATTERNS candidate (not promoted): "Read the Artin coordinate before
the residue."** Trigger: a finite cyclic unit group in which one of the
problem's multipliers is (or fails to be) a generator. Action: express the
statement in the discrete logarithm of that multiplier; what is expressible
there is residue-only and sheet-blind; what is not (sign, size, drift) is
the content. Counterindication: when the multiplier generates a proper
subgroup (`a = 7` in THM-4523, `q = 1 mod 24` in HYP-9174) the coordinate
splits into cosets and the coset is the first new datum. Evidence from
distinct threads: THM-4520 (Collatz Fourier), THM-4522/HYP-9174 (lonely
runner), THM-4523 (rigidity), THM-4532 (tournaments), HYP-9175 (coverage).

---

## 9. Reproduction and audit record

```bash
python 04-computation/experiments/collatz_coalescence_20261004_probes.py
# writes collatz_coalescence_20261004_probes.out; ~60 s (Probe A 6 s; the E2 dynamic programme to k = 6000 ~40 s)
```

Inputs cited from the repository: THM-4522's census (Theorem S), S19's
tables (`mazur_positive_density_20260928.md`, section 6 and "Global
statistics"), S23's Theorem 1, the eleventh note's Proposition 1, the
pricing note's Theorem 1 and Proposition 2, and the identification
`|Bad_k| = W_k` of THM-4479/THM-4495 (our exact scan gives `3, 4, 8, 13, 19,
38, 64, 128, 226, 367, 734, ..., 27328` for `k = 4..20`, equal to `W_4..W_20`
of the [Spitzer note](collatz_nodescent_order_20260926_spitzer_ballot.md)'s
table; the blind audit re-derives the counts with its own code).

Audit record: see the end of this file (appended by the session after the
subagent's verdict).
