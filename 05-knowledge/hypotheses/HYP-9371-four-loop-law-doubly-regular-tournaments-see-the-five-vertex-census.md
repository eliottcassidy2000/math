---
id: HYP-9371
title: "Four-loop law for doubly regular tournaments. (1) The four-loop (excess-3) term X_3 of the linear-forest loop expansion (111 cores) is computed by a validated engine. (2) For a doubly regular T, X_3 is an affine function of the 5-vertex census (PROVED), so across the doubly regular 23-tournaments it is exactly affine in the number TT5 of transitive 5-sets, Paley included. (3) Exact forest numbers give N_8(D) - N_8(D') = 7168 (TT5(D) - TT5(D')) at n = 19 and 23, which is the four-loop z^8 coefficient. (4) The strict four-loop prediction of H(D)/H(P_23) - 1 has the right sign and order but 38% of the exact TT5 slope (1.006e-7 against 2.665e-7). It has no Paley offset: the offset and the scatter first appear at forest length 10, where only five-loop clusters see the six-vertex census. Conjecture: the 5-census of a doubly regular tournament is affine in TT5 with an n-independent direction, hence N_8 = U_8(n) + 7168 TT5 for every order. (5) Five loops (section 8): with the 1076 excess-4 cores and the excess-5 cores with 5 vertices, Delta c_10 = delta_10 + Delta X_4^(10) + Delta X_5^(10) holds exactly for every class, so the six-vertex census enters at forest length 10 through five loops (the L = 10 part of Paley's offset); loops <= 5 predict the L <= 10 forest-length slices, 75-83% of every H difference, with no fitted parameters."
status: >
  OPEN (mac-mini-2026-10-22-loop4); independently audited 2026-10-22
  (04-computation/experiments/sierpinski_ocf_20261022/audit_loop4/REPORT.md). The audit reproduced the forest numbers of
  all 13 classes with its own DP, checked the engine on two further tournaments to z^16 and confirmed the Proposition.
  Its three errors and six minor items are corrected here.
  Section 8 (five loops, same session, later): independently audited 2026-10-22
  (04-computation/experiments/sierpinski_ocf_20261022/audit_loop5/REPORT.md; corrections applied, MISTAKE-667). Engine
  VERIFIED (equal to decomp.c on all excess <= 3 cores and, in the audit, on 499 excess-4 cores; X_4 equal to exact
  blow-up scaling on six tournaments). The z^10 identity on 13 classes is FINITE-EXACT: the audit recomputed it in
  exact rational arithmetic through an independent kernel and integer homomorphism numbers.
  Engine: VERIFIED. Exact blow-up scaling and exact set partitions agree to 2.3e-12 at excess 3, and the z^8
  coefficients of the 13 doubly regular 23-tournament differences agree exactly with exact forest numbers.
  Census reduction (2): PROVED below from THM-4757 (2).
  Forest numbers, the DRT(23)/DRT(19) tables, ΔN_8 = 7168 ΔTT5, the anatomy by forest length and the out-of-sample
  doubled-Paley DRT(23) class: FINITE-EXACT.
  Conjecture (C1) and its consequence (C2) for all n, and the asymptotic reading (C3): OPEN.
session: mac-mini-2026-10-22-loop4
source: 05-knowledge/hypotheses/HYP-9367-three-loop-law-the-loop-expansion-is-asymptotic.md
scripts:
  - 04-computation/experiments/sierpinski_ocf_20261022/cores_enum.py (+ cores_e1/2/3.json): cores of any excess, checked against the configuration-model count
  - 04-computation/experiments/sierpinski_ocf_20261022/decomp.c, gen_decomp_e3.py (+ .log): run decompositions (C port of loop_engine.decompositions_fast), 5-matrix weight tensors
  - 04-computation/experiments/sierpinski_ocf_20261022/loop_engine4.py: DP-optimal contraction, Cauchy z-series, two-matrix reduction for doubly regular T
  - 04-computation/experiments/sierpinski_ocf_20261022/excess3_validate.py (+ .out), drt_reduction_check.py (+ .out), paley_x3_scaling.py (+ .out)
  - 04-computation/experiments/sierpinski_ocf_20261022/drt23_classes.py (+ .json, .out), drt23_four_loop.py (+ .json, .out)
  - 04-computation/experiments/sierpinski_ocf_20261022/forest_numbers.c, drt23_forest_exact.py (+ .json, .out), drt23_anatomy.py (+ .json, .out), drt19_check.py (+ .out), drt_census5_direction.py (+ .out)
  - 04-computation/experiments/sierpinski_ocf_20261022/drt_c1_test.py (+ .out), drt23_doubled_check.py (+ .out, drt23_doubled.json): (C1), (C2) and the TT5 law on DRTs built by skew-Hadamard doubling
  - 04-computation/experiments/sierpinski_ocf_20261022/loop_engine5.py, cores_e4.json, cores_e5_v5.json, sign_network_basis5_check.py (+ .out), excess4_validate.py (+ .out), excess4_validate_R7_contour.py (+ .out): five-loop sign-network engine and its validation
  - 04-computation/experiments/sierpinski_ocf_20261022/drt23_five_loop.py (+ .out, .json), drt23_five_loop_analysis.py (+ .out, .json), paley_x4_scaling.py (+ .out): five loops on DRT(23)
related:
  - 01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md
  - 05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md
  - 05-knowledge/hypotheses/HYP-9367-three-loop-law-the-loop-expansion-is-asymptotic.md
  - 05-knowledge/hypotheses/HYP-9366-cosine-twin-and-the-alternating-end-of-the-forest-series.md
---

# HYP-9371 — four loops see the five-vertex census

Notation as in THM-4757:
* `Ŝ = S/n` and `κ = H 2^(n−1)/n! = Σ_L φ_L g(L)`, with `φ_L = N_L/n^L` and `g(L) = n^L/(n)_L`;
* `log Φ_T(z) = Σ_e X_e(z)`.

A tournament is *doubly regular* (DRT) if `S Sᵀ = nI − J`. `TT5` is the number of transitive 5-subtournaments.

## 1. The four-loop engine (VERIFIED)

**Cores.** The excess-3 cores are the connected multigraphs with loops, minimum degree ≥ 3 and `E = V + 3`. There are 111 of them:

| `V` | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| cores | 1 | 7 | 20 | 36 | 30 | 17 |

* `Σ 1/|G| = 1661/240`, where `|G|` is the automorphism count of `loop_engine.automorphisms`.
* Completeness check: `exp(Σ_cores x^V y^E/|G|)`, summed over excess ≤ 3, equals the configuration-model count `(1/V!) Σ_d (2E−1)!!/Π d_v!` exactly, at all 12 `(V, E)` points (`cores_enum.py`, `cores_enum.out`). The connected cubic part, `1105/1152`, is the Janson–Knuth–Łuczak–Pittel value `85085/82944` with the disconnected kernels removed.

**Decompositions.**
* `decomp.c` is a line-by-line C port of `loop_engine.decompositions_fast`. Its descriptor maps equal the Python ones key for key on all 3 + 15 cores of excess 1 and 2, and on the 8 excess-3 cores with `V ≤ 2`.
* It writes the weight tensor in the 5-matrix basis of `loop_engine_fast`, checked equal to `CoreB.Wb` at excess 1 and 2.
* A 6-vertex core has `4⁶ · 2⁹ = 2.1·10⁶` configurations. All 111 tensors take 209 s. They are sparse: 545 nonzeros for the 4-loop bouquet, and 12405–39285 of `5⁹` (0.6–2%) for the 6-vertex cores.

**Contraction.**
* An exact dynamic program over subsets of the ≤ 10 operands (the weight tensor and the edge propagators) minimises the multiply-add count under a cap on intermediate size.
* Each pairwise step is a batched matmul.
* Worst core at `n = 7`: `≈ 4·10⁸` multiply-adds (core 105).

**Validation** (`excess3_validate.out`). The z-series coefficients of `X_3`, from Cauchy sums of the engine, were compared with two exact rational computations that use no core theory: the blow-up scaling `log Ψ(Ŝ ⊗ J_m/m) = Σ_e m^(−e) X_e` (audit 2's `b_mscaling.py`) and set-partition enumeration (`b_partitions.py`).

| tournament | orders checked | max error |
|---|---|---|
| `R_7`, `P_7`, the non-circulant regular 7-tournament | `z⁶ … z¹²` | 2.3e-12 |
| `R_9` and two non-vertex-transitive regular 9-tournaments | `z⁶ … z¹⁰` (m-scaling), `z⁶, z⁸` (set partitions, equal to m-scaling as rationals) | 1.4e-12 |
| `P_11`, a non-vertex-transitive regular 11-tournament | `z⁶, z⁸` | 3.6e-16 |

* The overall maximum is `2.3·10^(−12)`, attained at the vanishing `z², z⁴` coefficients (Cauchy-sum roundoff). At `z⁶ … z¹²` the errors are ≤ `6·10^(−13)`.
* The audit added a non-vertex-transitive regular 9-tournament (gentourng class 8, exact `X_3 = −4304/59049, −1413992/4782969, −113969696/387420489` at `L = 6, 8, 10`) and a further random regular 11-tournament, agreeing to `10^(−15)`. It also pushed the non-circulant 7-tournament to `z¹⁶`, agreeing to `4·10^(−13)`; there the 6-vertex cores contribute about `10^(−4)`.
* `X_1` and `X_2` from the same engine equal the THM-4757 closed form and the HYP-9367 engine to 13 digits, at real and complex `z`.

**Doubly regular reduction.**
* For a DRT, `Ŝ² = −P/n` with `P = I − J/n`, and every propagator vanishes on the all-ones vector. So each propagator is `αP + γŜ`, where `α` and `γ` are the even and odd parts of the scalar propagator at `λ = ±i n^(−1/2)`.
* A core contributes `⟨W ×_e C2(z), Γ2⟩/|G|`, where `Γ2` is the `2^E` tensor of `{P, Ŝ}` homomorphism numbers of the core. It does not depend on `z`.
* On the doubly regular `P_7` and `P_11`, this two-matrix evaluation equals the general engine to 13–14 digits for `X_1, X_2, X_3` (`drt_reduction_check.out`).

**Power counting** (`paley_x3_scaling.out`, `q = 19 … 83`).
* `q X_1 → −7/2` and `q² X_2 → 29` (both reproduced).
* `q³ X_3(P_q)` = −370.3, −363.5, …, −345.5, extrapolating to `≈ −341 ± 1`. The fit residual is 0.03, not 10^(−7), because `TT5(P_q)` carries arithmetic fluctuations. Adding `TT5(P_q)` as a regressor drops the residual to `4·10^(−5)`, with limit −340.7 (audit).
* The loop coefficients grow like 3.5, 29, 341 (ratios 8.3 and 11.8). This is the asymptotic-series signature of HYP-9367.

## 2. Four loops see exactly the five-vertex census (PROVED)

**Proposition.** Let `T` be doubly regular of order `n`. Then

    X_3(T; z) = A_n(z) + Σ_τ B_(n,τ)(z) c_τ(T),

where `c_τ` counts the 5-subtournaments of type `τ`, and the coefficients `A_n, B_(n,τ)` depend only on `n`.

*Proof.* By THM-4757 (2), `X_3 = Σ μ(π) Z_π / p!` over connected clusters of excess 3. Here `Z_π = Σ_f Π_arcs (zŜ)_(f(·) f(·))` is a homomorphism sum, over all maps `f`, of the quotient multigraph. That multigraph is a subdivided core with ≤ 6 branch vertices.

1. *Collapse onto the core.* Summing over the images of the degree-2 blocks of a core path gives a product of factors `±zŜ`, i.e. `±(zŜ)^λ`. For a DRT, `Ŝ² = −P/n`, so this lies in `span{P, Ŝ}`. Odd `λ` gives a multiple of `Ŝ`, even `λ ≥ 2` a multiple of `P`. Hence `Z_π = Σ_(f: core → [n]) Π_e (a_e P + c_e Ŝ)` with scalar coefficients.
2. *Parity.* In a given cluster each core path has a definite length `λ_e`. So the cluster's term carries one labelling `t ∈ {P, Ŝ}^E`, with `Ŝ` on the odd paths, and its total length is `L = Σ λ_e ≡ #Ŝ(t) (mod 2)`. Reversing an odd run is a sign-reversing involution on clusters: it maps `π` to its mirror image and preserves the excess (THM-4757 §1). So `X_3` has only even powers of `z`: the clusters with odd `#Ŝ` cancel, and `X_3` is the sum over clusters with even `#Ŝ`.
3. *Reduction to five vertices.* Split each homomorphism sum by the kernel of `f` into injective sums. Inside an injective sum, `P` between distinct points is `−1/n`, `P` on a point is `1 − 1/n`, and an `Ŝ`-loop is 0.
   * A 6-vertex core has 9 edges, so an even `#Ŝ` is at most 8. Then some vertex `w` has `Ŝ`-degree ≤ 2.
   * Summing `f(w)` over the points not used by the other five vertices gives, by `Ŝ1 = 0`, `Ŝ_xx = 0` and `ŜŜᵀ = (nI − J)/n²`, an `n`-dependent combination of injective sums on ≤ 5 vertices.
   * *Step 2 is not needed for the conclusion* (audit). An all-`Ŝ` 9-edge labelling of a cubic 6-vertex core never sees the 6-census: loops give 0, a double edge reduces to fewer vertices, and `K_(3,3)` and the prism each have an automorphism reversing an odd number of edges, so the injective sum vanishes.
4. *Census.* An injective `Ŝ`-digraph sum on `k ≤ 5` vertices is a fixed linear combination of the `k`-vertex census counts, times `n^(−#arcs)`. For DRTs the census on ≤ 4 vertices is determined by `n` (HYP-9365 (d)). Every out- and in-neighbourhood is a regular tournament of order `(n−1)/2`, which fixes the numbers of `TT4`, of 3-cycles with a source and of 3-cycles with a sink, hence all four 4-vertex counts. ∎

**Hierarchy.** The same count says what the next order sees. At excess 4:
* a 6-vertex core has 10 edges, an even number, so all can carry `Ŝ` with every degree ≥ 3, giving genuinely 6-vertex counts;
* a 7-vertex core has 11 edges, so at most 10 carry `Ŝ`, and it reduces to 6-vertex counts by the argument of step 3;
* an 8-vertex cubic core has 12 edges, all of which can carry `Ŝ`, giving 8-vertex structure.

So the first dependence of `log Φ` on more than the 5-census is at five loops.

**Numerics** (`drt23_four_loop.out`).
* Over the 13 classes, `X_1`, `X_2` and the `V ≤ 4` part of `X_3` coincide exactly; the series differ by at most `8.8·10^(−22)`.
* `ΔX_3` is affine in `TT5` with residual `3·10^(−18)`, and `P_23` lies on the line.
* The genuinely 6-vertex part of `ΔX_3` (the 6-vertex cores minus their 15 single vertex merges) is nonzero, `≈ 0.36%` of `ΔX_3`. It is also exactly affine in `TT5`, as step 3 says it must be.
* Hostile control (audit): a 10-edge statistic (prism plus a chord) separates `D23_1` and `D23_2`, which share their 5-census, as 71920 against 67824. The 6-vertex part of `ΔX_3` is nevertheless identical for the two. So it is not the case that the 6-census happens to be affine in `TT5`.

## 3. Prediction against the exact doubly regular 23-tournaments (FINITE-EXACT)

**Data.** The 15 codes of `drt23_codes_all.txt` fall into 12 nauty classes; with `P_23` that makes 13. Class order as in `drt_census_hierarchy.out`.
* `H` was recomputed for all 13 (`hk64.c`). It reproduces the four values of `drt_results.out`.
* The classes 6/10 and 7/11 are converse pairs (`drt23_classes.json`).

**Universal parts at `n = 23`.** `X_1 = −0.1428639`, `X_2 = +0.0537243`, `X_3(P_23) = −0.0298764`, so `n³X_3 = −363.5`.

**Columns** (`exact` is `log H(D) − log H(P_23)`).
* `strict` is `ΔX_3(1)`, the complete non-universal `n^(−3)` term. Through three loops all doubly regular tournaments of one order agree (HYP-9365 (d)), and the measure enters only at the next order.
* `dressed` adds to `ΔX_3` its own `n^(−4)` finite-size response, `Σ_j [1 + (j D_1 + C(j,2))/n] [z^j] ΔX_3` with `D_1 = 2n/(n+1)`. This is NOT the complete `n^(−4)` term: `ΔX_4(1)` is missing.

| class | TT5 | exact | strict `ΔX_3(1)` | dressed | exact/strict |
|---|---|---|---|---|---|
| D23_0 | 1478 | −1.6106e-5 | −4.024e-6 | −1.199e-5 | 4.00 |
| D23_1 | 1466 | −1.9373e-5 | −5.232e-6 | −1.559e-5 | 3.70 |
| D23_2 | 1466 | −1.9345e-5 | −5.232e-6 | −1.559e-5 | 3.70 |
| D23_3 | 1470 | −1.8104e-5 | −4.829e-6 | −1.439e-5 | 3.75 |
| D23_4 | 1492 | −1.2372e-5 | −2.616e-6 | −7.793e-6 | 4.73 |
| D23_5 | 1494 | −1.2077e-5 | −2.415e-6 | −7.193e-6 | 5.00 |
| D23_6 = D23_10ᵀ | 1484 | −1.4512e-5 | −3.421e-6 | −1.019e-5 | 4.24 |
| D23_7 = D23_11ᵀ | 1482 | −1.5276e-5 | −3.622e-6 | −1.079e-5 | 4.22 |
| D23_8 | 1490 | −1.2560e-5 | −2.817e-6 | −8.392e-6 | 4.46 |
| D23_9 | 1486 | −1.3949e-5 | −3.220e-6 | −9.591e-6 | 4.33 |

**TT5 slopes** over the 12 non-Paley classes, per transitive 5-set, in units of `H/H(P_23)`:

| | slope | `× H(P_23)` | line at Paley (`TT5 = 1518`) | max residual |
|---|---|---|---|---|
| exact | `2.6648·10^(−7)` | `4.1999·10⁹` (the HYP-9367 law) | `−5.47·10^(−6)` | `3.7·10^(−7)` |
| strict | `1.0061·10^(−7)` (38%) | `1.586·10⁹` | 0 | `3·10^(−18)` |
| dressed | `2.9972·10^(−7)` (112%) | `4.724·10⁹` | 0 | `2·10^(−17)` |

**Order 19** (`drt19_check.out`). For `D_19` against `P_19`, with `ΔTT5 = 18`:
* exact `+3.173·10^(−5)`;
* strict `+8.54·10^(−6)` (exact/strict = 3.72);
* dressed `+2.92·10^(−5)` (92%).

**Out-of-sample class** (`drt23_doubled_check.out`). Normalising the doubled Paley skew-Hadamard matrix of order 24, `[[K, K], [−Kᵀ, Kᵀ]]` with `K` from `P_11` (`drt_c1_test.py`), gives a new DRT(23) class with `TT5 = 1430`. That is `ΔTT5 = −88`, far outside the annealed range −52 … −24.
* Exact `H = 15759752857737405`, so `H/H(P_23) − 1 = −2.8814·10^(−5)`.
* The empirical `TT5` law of the 12 annealed classes predicts its `H` to `−1.1·10^(−7)`, inside the in-sample residual `3.7·10^(−7)`. The line, with its offset from Paley, extrapolates; Paley is the outlier, not the annealed family.
* Four loops: strict `−8.85·10^(−6)` (exact/strict = 3.25), dressed `−2.64·10^(−5)` (92% of exact).
* Exact `ΔN_8 = 7168 ΔTT5` again. By forest length: `L = 8` 43%, 10 39%, 12 14%, 14 3%.

**What the four-loop term gets right.**
* The sign and the ordering: exact `H` is monotone in `TT5`, and `D23_1`, `D23_2` share their 5-census and differ by `2.8·10^(−8)`.
* The exact affine dependence on the 5-census.

**What it does not give.**
* The magnitude: the strict term is 20–27% of each annealed-class difference and 38% of the slope; the dressed value is 60–81%. For the doubled-Paley class and for `n = 19` the strict term gives 31% and 27%, the dressed 92%.
* Paley's offset `+5.47·10^(−6)` from the non-Paley line, and the `3.7·10^(−7)` scatter. By §2, four loops cannot depend on anything beyond `TT5` here.

## 4. Anatomy from the exact forest numbers (FINITE-EXACT)

**Forest numbers.** `F(T, x) = Σ_σ x^(#forward steps)` was computed modulo three primes at `x = 0 … 22` by a Held–Karp DP (`forest_numbers.c`, 0.8 GB, 1.5 s per point), then interpolated and lifted by CRT. THM-4757 §5 then gives `N_L` exactly.

All 13 classes pass these checks: `f_0 = f_22 = H`, `Σ f_k = 23!`, `f` palindromic, `N_L = 0` for odd `L`, and `Σ_L N_L/(n)_L = H 2^22/23!` exactly.

**`N_8` is the four-loop term.**
* `N_2, N_4, N_6` are universal: 506, 91080, 7832880.
* `N_8(D) − N_8(P_23) = 7168 ΔTT5` for all 12 annealed classes and the doubled-Paley class.
* `ΔN_8 = 7168 ΔTT5` also at `n = 19`: `N_8 = 46577664` and `46706688`, `ΔTT5 = 18`.
* `7168 = 7·2^10`. The coefficient of `z⁸` in `Δ log Φ` equals the engine's `δ_8 = [z⁸]ΔX_3` exactly. No other excess is non-universal at `z⁸`: excess ≥ 4 has ≤ 4 blocks there, and excess ≤ 2 reduces to ≤ 4-vertex cores.

**Higher coefficients are not four-loop.** `D23_0`, exact `Δc_L` against four-loop `δ_L`:

| `L` | 8 | 10 | 12 | 14 | 16 |
|---|---|---|---|---|---|
| exact `Δc_L` | −3.661e-6 | +1.711e-6 | −6.61e-7 | +2.61e-7 | −1.08e-7 |
| four-loop `δ_L` | −3.661e-6 | −3.18e-7 | −4.2e-8 | −3.0e-9 | −2.6e-10 |

* The exact coefficients alternate. The successive ratios −0.47, −0.39, −0.39, −0.41 approach `1/u_0 = −0.448` per `z²`, where `u_0` is the zero of the forest polynomial `Φ(u)` (`u = z²`) nearest 0: −2.23181 for `P_23`, −2.23167 for `D23_0`. So `log Φ_T` has radius `|z| = 1.49`.
* Every individual `X_e` is singular only at `z² = ±n`, so no finite loop order sees this zero. The zero is on the negative `z²` axis (`z = ±1.494i`), like the nearest forest-polynomial zeros of almost-regular tournaments in the S54 lane: `arg z² = π`, `|z|` 1.26–1.37 for `n = 8–14`, growing with `n` (`04-computation/experiments/loop_nonregular_20261022/phi_zeros_model.out`).
* At `z¹⁰` the higher-loop terms `Δc_10 − δ_10` are −0.38 to −0.61 times the four-loop `δ_8`, against `δ_10/δ_8 = +0.087`.

**Where the `H` difference lives.** Exactly, `ΔH/H = Σ_L g(L) Δφ_L / κ`:

| | `L = 8` | 10 | 12 | 14 | ≥ 16 |
|---|---|---|---|---|---|
| D23_0 | 35% | 43% | 18% | 3.4% | 0.2% |
| D23_1 | 38% | 42% | 16% | 3.3% | 0.2% |
| D23_4 | 30% | 46% | 20% | 3.7% | 0.06% |
| `D_19` vs `P_19` | 53% | 43% | 0.4% | 3.5% | 0.2% |

* The four-loop term gives the `L = 8` slice exactly, and nothing else exactly.
* Propagated through the exact `Φ_(P_23)`, the four-loop `ΔX_3` overshoots every `L ≥ 10` slice. Its sum, `−8.4·10^(−4)` for D23_0, is 50 times too large, because the exact response weights grow to `w_22 ≈ 1.4·10⁷` at the top of the forest series. The same happens with the loop measure `exp(X_0 + X_1 + X_2)`, whose `Σ_L g φ_L` is 214 instead of `κ = 2.557` (negative control in `drt23_four_loop.out`).
* Paley's offset, by forest length: `L = 10` 62%, `L = 12` 34%, `L = 14` 5%, −0.7% at `L ≥ 16`, none at `L = 8`. At `L = 10` only five-loop clusters (6 blocks) can carry the 6-census. At `L ≥ 12` excess 4–6 clusters with 6–8 vertices can contribute.
* `D23_1` and `D23_2` (equal 5-census) first differ at `L = 10`, consistent with the six-vertex census entering at five loops.

## 5. Why four loops do not account for the differences, and what they do give

* **Structurally,** the four-loop term is a pure five-vertex effect (§2). The doubly regular 23-tournaments span a one-parameter 5-census, so four loops predict an exact line in `TT5` through Paley. The observed offset of Paley and the `3.7·10^(−7)` scatter first appear at `L = 10`, where only five-loop clusters bring in the 6-vertex census. Their `L ≥ 12` parts may also involve higher loops.
* **Quantitatively,** at `n = 19` and `23` the `n^(−4)` correction is larger than the `n^(−3)` term, so the strict truncation is not quantitatively useful at these orders.
  * The strict four-loop term is 20–27% of each `n = 23` difference and 38% of the exact `TT5` slope.
  * Its own `n^(−4)` finite-size response is twice itself: `dressed/strict = 2.98`, with factor `1 + (8D_1 + 28)/n = 2.88` on the dominant `z⁸` coefficient. It brings the slope to 112% of exact at `n = 23` and the `n = 19` difference to 92%.
  * The remaining `n^(−4)` pieces (`ΔX_4`) are not computed.
* **Exactly,** four loops fix the leading forest-number difference, `ΔN_8 = 7168 ΔTT5`. Its slice of the `H` differences is 28–38% at `n = 23` and 53% at `n = 19`.

## 6. Conjectures

**(C1)** The 5-vertex census of a doubly regular tournament of order `n` is affine in `TT5`, with an `n`-independent direction. Per transitive 5-set, the change in each of the 12 isomorphism types of 5-tournaments is as follows. Types are labelled by scores and number of Hamiltonian paths; the two types with key `(1,1,2,3,3)|9` are split by whom the score-2 vertex beats.

| type | change |
|---|---|
| `(0,1,2,3,4)\|1` (transitive) | +1 |
| `(0,1,3,3,3)\|3` | −1 |
| `(0,2,2,2,4)\|3` | −1 |
| `(0,2,2,3,3)\|5` | +1 |
| `(1,1,1,3,4)\|3` | −1 |
| `(1,1,2,2,4)\|5` | +1 |
| `(1,1,2,3,3)\|9`, score-2 vertex beats both score-1 vertices | −3 |
| `(1,1,2,3,3)\|9`, score-2 vertex beats a score-1 and a score-3 vertex | +5 |
| `(1,2,2,2,3)\|11` | +1 |
| `(1,2,2,2,3)\|13` | −7 |
| `(1,2,2,2,3)\|15` | +3 |
| `(2,2,2,2,2)\|15` | +1 |

The census is by complete isomorphism type (`drt_census5_direction.py`, canonical forms). The earlier key `(scores | Hamiltonian paths)`, used in `drt_census_hierarchy.py` and in HYP-9367's "rank 1 across the 11 census types", merges the two `(1,1,2,3,3)|9` types into a net +2. With all 12 types the differences still have rank 1 (audit).

FINITE-EXACT, on 16 non-Paley classes (14 up to converse) at three orders, each against the Paley tournament of its order:
* `n = 19`: `D_19` (`drt_census5_direction.out`);
* `n = 23`: the 12 annealed classes, plus the doubled-Paley class with `ΔTT5 = −88`, which was not used to find the vector;
* `n = 31`: two classes from the twice-doubled Paley skew-Hadamard matrix of order 32, not converse to each other, both with `ΔTT5 = −300` and hence one census data point (`drt_c1_test.out`).

The Paley switching families of `P_27` and `P_31` are single classes.

**(C2)** For every doubly regular `T`, `N_8(T) = U_8(n) + 7168 TT5(T)`. This follows from (C1). The non-universal part of `[z⁸] X_3` involves only single-arc paths on 5-vertex cores, so the non-universal part of `n⁸ [z⁸] X_3` is an `n`-independent linear combination of 5-census counts.
* FINITE-EXACT from exact forest numbers at `n = 19` (1 pair) and `n = 23` (13 pairs).
* VERIFIED by the engine at `n = 31` (2 pairs): `n⁸ Δδ_8/ΔTT5 = 7168.000000`.

**(C3)** The asymptotic reading of HYP-9367 for differences: `log H(D) − log H(D′) = ΔX_3(1)(1 + O(1/n)) + O(n^(−4))`. The `O(n^(−4))` term carries the 6-vertex census and the 8-vertex cubic structure of excess 4; the 6-census comes first, at `L = 10`. At `n = 19–23` the `O(1/n)` factor is ≈ 2.6–3.7: 2.65 for the slope, 3.25 for the doubled-Paley class, 3.72 at `n = 19`. That most of it is the first-order finite-size dressing is a HEURISTIC reading: the first-order weights are far from the exact response weights (`w_8 = 2.88` against 19.3).

## 7. Refutation targets and next steps

**Refutation targets.**
* A doubly regular pair at any order whose 5-census difference is not proportional to the vector of (C1).
* A pair with `ΔN_8 ≠ 7168 ΔTT5`.

**Next steps.**
* **Five loops.** `ΔX_4` needs the excess-4 cores (up to 8 vertices, 12 edges) and their decompositions (`4⁸ · 2^12` configurations for a cubic 8-vertex core). Testable target: the exact `Δc_10 − δ_10`, e.g. `+2.029·10^(−6)` for D23_0, which is `ΔX_4^(10) + ΔX_5^(10)`. `K_5` enters at excess 5, `z¹⁰`.
* **Third-order doubly regular constant.** The engine makes the universal `n^(−3)` coefficient of `κ_DRT = e(1 − 3/(2n) + 25/(8n²) + c_3/n³ + …)` computable from `q³X_3(P_q) → ≈ −341` plus cross and finite-size terms.

## 8. Five loops (same session; independently audited, `audit_loop5/REPORT.md`)

**Cores.**
* Excess 4 has 1076 cores. By vertex count `V = 1 … 8`: 1, 10, 48, 153, 277, 323, 193, 71. The 71 with 8 vertices are the connected cubic multigraphs with loops on 8 vertices.
* For `z¹⁰` the excess-5 cores with `V ≤ 5` are also needed: 2046 of them, 1451 with `V = 5`.
* Canonical forms are by individualisation–refinement (`cores_enum.py`; it reproduces excess 1–3). The configuration-model identity holds exactly at all 20 `(V, E)` points up to excess 4, and at all 25 up to excess 5 (top level `V ≤ 5`).

**Engine.** Explicit enumeration does not scale: a cubic 8-vertex core has `4⁸ · 2¹² = 2.7·10⁸` configurations, and a dense weight tensor `5¹²` entries. `loop_engine5.py` instead writes the decomposition sum of a core as a tensor network over half-edge orientation signs `o(x) = ±1`.
* A vertex carries `Σ_π mob(d − k − |π|) Π_(pairs) [o_x = −o_y]` over partial matchings `π` of its free half-edges, where `k` counts pairs already fixed.
* An edge carries an m1 term `R_(o_t) [o_t = −o_h]` and an m2 term `Q_(o_t, −o_h)`.
* The contraction also counts closed runs, cycles of m1 edges joined by pass-throughs. Inclusion–exclusion over sets `A` of half-edge-disjoint closed trails removes them: `Valid = Σ_A (−1)^|A| N(A)`. In `N(A)`, the runs of `A` are forced, with one orientation sign per run, and everything else is free.
* For a doubly regular `T` the edge weights are Ising-like:
  * m1: `P` coefficient `a`, `Ŝ` coefficient `g·o_t`;
  * m2: `P` coefficient `((b − a) + (b + a) o_t o_h)/2`, `Ŝ` coefficient `(h/2)(o_t − o_h)`;
  * with `u² = −z²/n`: `a = u²/(1 − u²)`, `b = u²/(1 + u²)`, `g = z/(1 − u²)`, `h = −2u²z/((1 + u²)(1 − u²))`.

**Validation.**
* The 5-basis weight tensors equal `decomp.c`'s integer tensors on all 129 cores of excess ≤ 3, in 4.6 s against 209 s (`sign_network_basis5_check.out`).
* In matrix mode, `X_1 … X_3` equal `loop_engine4` to `10^(−16)` on a non-vertex-transitive regular 11-tournament (reproduced by the audit: `audit_loop5/a5_task6_modes.out`).
* The `X_4` z-series equals the exact blow-up scaling (`excess4_validate.out`):

| tournament | orders | max relative error |
|---|---|---|
| `P_7` | `z⁶ … z¹⁴` | `4.5·10^(−14)` |
| non-circulant 7-tournament | `z⁶ … z¹⁴` | `4.8·10^(−11)` |
| non-vertex-transitive 9-tournament | `z⁶ … z¹⁰` | `3.9·10^(−11)` |
| `R_7`, `R_9` | `z⁶ … z¹⁴`, `z⁶ … z¹⁰` | `2.1·10^(−9)`, `2.6·10^(−9)` |

  The `R_7`, `R_9` errors are contour aliasing: their largest eigenvalues are the biggest. With 64 points on `|z| = 0.35`, `R_7` agrees to ≤ `3·10^(−14)` at every order (`excess4_validate_R7_contour.out`), and `R_9` to `1.8·10^(−15)` (audit).
* The audit compared the engine's 5-basis tensors with `decomp.c` exactly on all 489 excess-4 cores with `V ≤ 5` and ten with `V = 6`. It checked the 7- and 8-vertex cores against a first-principles kernel (20514 values, 0 mismatches). It also matched `X_4` against exact blow-up scaling on a further regular 9-tournament up to `z¹²` (`8.4·10^(−15)`).
* For the DRT(23) classes, the cores with ≤ 4 vertices are universal at excess 4 (class differences ≤ `2·10^(−18)`).

**Power counting** (`paley_x4_scaling.out`). `q⁴ X_4(P_q)` = 6341, 6055, 5739, 5513 at `q = 19, 23, 31, 43`, extrapolating to about `5.0·10³` (4926 from the last two points, 5038 from a quadratic fit in `1/q`). The universal loop terms of `P_q` alternate in sign and grow: about `−3.5/q`, `+29/q²`, `−341/q³`, `+5·10³/q⁴`. At `q = 23` they are −0.143, +0.054, −0.030, +0.022.

**The `z¹⁰` identity (FINITE-EXACT).** At `L = 10` only excess 3, 4 and 5 can differ between doubly regular tournaments of one order: excess ≥ 6 has ≤ 4 blocks there, and excess-5 cores with ≥ 6 vertices have ≥ 11 edges. So, class by class,

    Δc_10 (exact, forest numbers) = δ_10 + ΔX_4^(10) + ΔX_5^(10).

It holds to ≤ `3·10^(−18)` for all 13 non-Paley classes (11 up to converse; `drt23_five_loop.out`). The audit recomputed it exactly, as rationals, through an independent kernel and integer homomorphism numbers (`audit_loop5/a5_task5_z10_*.out`). For D23_0: `+1.7109136619·10^(−6) = −3.183739·10^(−7) + 2.74117952·10^(−6) − 7.11891954·10^(−7)`.
* The non-TT5 part of `Δc_10` is reproduced exactly. That is Paley's offset `+8.694·10^(−7)`, the scatter, and the `6.33·10^(−9)` split of D23_1 / D23_2, whose 5-censuses are equal. The six-vertex census enters through excess-4 clusters with 6 blocks.
* `ΔX_4` and `ΔX_5` have no `z⁸` term (≤ `2·10^(−18)`), so `Δc_8 = δ_8` is again exact.
* At `z¹⁰` the loop orders do not decrease at `n = 23`. For D23_0, `ΔX_4^(10)` is 8.6 times `δ_10` and of opposite sign, and `ΔX_5^(10)` is a quarter of `ΔX_4^(10)`, again of opposite sign. Across the classes these ratios range over 6.6–10.1 and 0.22–0.34.

**Predictions** (`drt23_five_loop_analysis.out`; 13 non-Paley classes including the doubled-Paley one; slopes and offsets in this section are 13-class fits, while sections 3–4 quote the 12-class values `2.6648·10^(−7)` and `+5.47·10^(−6)`).

| class | exact | strict `n^(−3)` | strict `n^(−4)` (complete) | loops ≤ 5, slices `L ≤ 10` |
|---|---|---|---|---|
| D23_0 | −1.6106e-5 | −4.02e-6 (25%) | −8.44e-6 (52%) | −1.266e-5 (79%) |
| D23_1 | −1.9373e-5 | −5.23e-6 (27%) | −1.071e-5 (55%) | −1.552e-5 (80%) |
| D23_4 | −1.2372e-5 | −2.62e-6 (21%) | −5.83e-6 (47%) | −9.42e-6 (76%) |
| D23_5 | −1.2077e-5 | −2.41e-6 (20%) | −5.49e-6 (45%) | −9.09e-6 (75%) |
| D23_dbl | −2.8815e-5 | −8.85e-6 (31%) | −1.744e-5 (61%) | −2.384e-5 (83%) |

* **Strict `n^(−4)`.** It is `ΔX_3(1) + ΔX_4(1)` plus the first-order finite-size response of `ΔX_3`, the complete `n^(−4)` difference.
  * It gives 45–61% of each difference and a `TT5` slope of `1.87·10^(−7)`, 71% of the exact `2.650·10^(−7)`.
  * Paley's offset comes out `+9.6·10^(−7)`, 17% of the exact `+5.53·10^(−6)`. The non-`TT5` residuals correlate with the exact ones at `+0.994`, at about 18% of their size.
  * `ΔX_4(1)` itself (`+3.5·10^(−6)` for D23_0) has the sign opposite to the differences. The `n^(−4)` gain comes from the finite-size response of `ΔX_3`.
* **By forest length.** Loops ≤ 5 fix `Δφ_8 = δ_8` and `Δφ_10 = Δc_10 + φ_2 Δc_8` exactly. With `g(L) = n^L/(n)_L` and the universal `κ = e(1 − 3/(2n) + 25/(8n²))`, with no exact input, they give:
  * 75–83% of every difference;
  * a `TT5` slope of `2.32·10^(−7)` (88%);
  * a Paley offset of `+3.39·10^(−6)` (61%).

  These slices equal the exact ones: to `1.6·10^(−12)` with the exact `κ`, and to `3.1·10^(−5)` relative with the universal one. The remainder sits at `L ≥ 12`. At `L = 12` every excess from 3 to 7 contributes. Excess 3 and 7 enter only through the 5-vertex census; excess 7 gives `ΔX_7^(12) = +3.853·10^(−9)` per transitive 5-set (audit). Structure on 6–8 vertices comes from excess 4–6.

**What this settles.**
* Paley's offset, the scatter and the equal-5-census split go beyond the 5-vertex census. Their `L = 10` parts are exactly six-vertex census effects: 61% of the offset and 89% of the D23_1/D23_2 split. They first appear there, through five loops, and the engine reproduces them exactly. The `L ≥ 12` parts may also carry 7- and 8-vertex structure.
* The strict `1/n` expansion improves slowly at `n = 23`: 25% at `n^(−3)`, 52% at `n^(−4)` for D23_0.
* Organised by forest length, the expansion is exact slice by slice. Through five loops it fixes about four fifths of every difference.
* Next step: `z¹²` needs excess 5 with ≤ 7 vertices, excess 6 with ≤ 6 vertices and excess 7 with 5 vertices. The exact target is in `drt23_anatomy.out`.

