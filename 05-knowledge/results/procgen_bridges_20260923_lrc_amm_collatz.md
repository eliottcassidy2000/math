# LRC(14), AMM 12592 and the Collatz work: typed bridges, where the golden ratio comes from, and the owner's inequality principle made concrete

**Status.** SYNTHESIS with FINITE-EXACT computations. No HYP or THM file was created, and nothing here changes a canonical status. **LRC(14) is OPEN. Collatz, the Periodicity Conjecture (PC) and both halves of E-SCC are OPEN. AMM 12592 has `C*` in `[11/8, C_*]`, with HYP-9128 and HYP-9129 OPEN.**
* **AUDIT (orchestrator, 2026-09-23).** The LRC row `{1..12, 5460}` was re-computed independently by an exact union of rational bad arcs: the lonely measure is `P(M=0) = 301/10296 = 0.0292347`. This equals the lane's excised-Bonferroni value. The AM–GM reading of the Collatz drift is elementary and was re-checked: for `x/2, (qx+1)/2` the arithmetic-mean step factor is `(q+1)/4`, which equals 1 iff `q = 3`, and the geometric mean `sqrt(3)/2` gives the drift `-log(2/sqrt3)` per step. The 'Wang–Massey' citation for the `+-1` Hankel determinants of `sum w^(2^j)` is from memory and UNVERIFIED; the determinants were checked exactly only for `n <= 100` (and mod 2 for `n <= 300`).
* **PROVED here:** only elementary identities. These are the eigen-data identifications in section 2, the fold identity `p = sin^2(pi x)` in B2, and the pointwise Bonferroni/owner-count defect identities in section 3.
* **Exact rational arithmetic:** every LRC measure below, including the two "excised Bonferroni" loneliness certificates in section 4.1.
* **CITED, not re-read:** Wang–Massey (1986) and Bugeaud–Kim (via Corollary D1 of the hard-class note).
* **Bridge labels:** each bridge is labelled REAL, ANALOGY or NUMEROLOGY, and each label comes with the test that decided it.

Session `collatz-procgen-20260922`, bridges lane (mac-mini), 2026-09-23. It answers the owner's question: "think how k=13 lonely runner and amm12592 relate to our collatz work", together with the owner's inequality principle.

* **Scripts** (`04-computation/experiments/`):
  * `procgen_bridges_20260923_phi.py`
  * `..._adelic.py`
  * `..._defects.py`
  * `..._twoblock.py`
  * `..._run.py` (runner)
* **Output:** [`procgen_bridges_20260923.out`](procgen_bridges_20260923.out).

---

## 0. Answers in brief

1. **One genuine shared engine.** AMM and Collatz close their arguments with the product formula: a nonzero integer `N` has `prod_v |N|_v = 1`. They differ in which place carries the analysis.
   * **AMM:** THM-4467 is Pólya's Hankel-determinant theorem. The analytic place is `inf`, and the gain comes from continuation.
   * **Collatz:** Theorems R, D and Y, and the theta lane's H1 (Bézivin's Hankel method), use the same engine with the places exchanged. The analytic place is 2 and `inf` is the arithmetic place. The Hankel valuations grow like `n^2` for Pólya and like `n^3` for q-series tails.
   * **LRC:** LRC uses the other half of the pair, Siegel/Dirichlet pigeonhole (THM-2052) and transference (THM-4009), to produce short relations. Nothing in the LRC program plays Liouville's role.
   * **Label:** REAL for AMM–Collatz, ANALOGY for LRC.
2. **The adelic-capacity (Pólya–Bertrandias) form cannot carry Collatz.** The natural function forms of the Bernstein series (both parametrizations) and of the theta series have `sum_v log R_v = liminf - limsup <= 0`. Its terms are `{2,3}`-units, and the value is exactly 0 when a density exists. The AMM lacunary class is also at 0. Theorem A's gain is continuation beyond the circle of convergence, which lacunary Collatz series do not have.
   * At the prime 2, AMM is inert as well. Every shifted Hankel determinant `H_n^(1)`, `H_n^(2)` of the parity skeleton `sum w^(2^j)` equals `+-1`, checked exactly for `n <= 100`. It is odd for `n <= 300` (checked) and for all `n` by Wang–Massey.
   * By AMM Lemma P, the same Hankel determinants of `phi` are odd for **every complement-symmetric fair extractor**. So no 2-adic Bertrandias gain exists.
3. **The golden ratio.** The four occurrences are:
   * Collatz Theorem S: `mu(mu-1) < phi`;
   * Collatz Lemma P (Tschakaloff): `mu_bar < phi`;
   * AMM THM-3027/THM-3009: `tau* = phi^-2` and `C_* = 1 + log_5 phi^2`;
   * AMM Theorem B / Long's point: `w = -1`, `p = -1/phi`.

   All four are eigen-data of the smallest hyperbolic unit `[[1,1],[1,0]]` of `GL_2(Z)` or of its square. The co-variation tests (PHI2) show that four unrelated knobs move them.
   * Across programs the shared `phi` is **NUMEROLOGY**.
   * Inside AMM it is **one mechanism with three faces (REAL)**.
   * The genuine common pattern is a meta-analogy: **`phi` bounds a one-scale or one-parameter method, and coupling scales beats it.**
     * Theorem S's 1.8668 becomes Bugeaud–Kim's 2.50994 (PROVED).
     * Collatz Lemma P's `phi` becomes 7/4, 28/11 and 2.878 (theta lane H1–H3, PROVED today).
     * AMM's golden `C_*` falls to about 1.570 with super-blocks (NUMERICAL, HYP-9128).
   * **LRC(14) has no golden threshold.** Its frontier constants are rationals built from sevenths, the clocks and pigeonhole.
4. **The dyadic structure is shared in name only.**
   * **AMM:** the parity skeleton is automatic, of Artin–Schreier type, and has a perfect GF(2) linear-complexity profile.
   * **Collatz:** the theta and cube skeletons have GF(2)-random profiles.
   * **LRC:** the "clock two" is the same Z/2 quotient as AMM's `p <-> 1-p`, under `p = sin^2(pi x)`. That identity is exact, but it transports AMM's complex continuation to imaginary times.
   * **Label:** ANALOGY.
5. **S596 "same two-block question".** The THM-4447 capacity rule, re-run for other runner counts `n` and tail counts `r`, shows how the LRC(14) small clocks `{2,3,4}` arise.
   * With `r = 3` tails they are `{2,3,4}` for every `n` in `{12,14,16,20,28}`, whether or not `2n-1` is a power of 3. So they are pigeonhole clocks.
   * No logarithmic linear form occurs in the LRC residual.
   * HYP-9122's finite form lives on the record denominators of `log_2 3`.
   * **Verdict:** the shared shape (easy by density, hard on a thin arithmetic residual) is ANALOGY. The "prime-2 versus prime-3 seam" is NUMEROLOGY for LRC(14). The Baker import to LRC has no target.
6. **The owner's principle.** Every inequality has an exact defect identity. The computed ledger (section 3) finds three different situations:
   * **LRC: the defect is localized.** The Bonferroni defect `E[C(M-1,5) 1{M>=1}]` sits almost entirely in the neighbourhoods of the row's pack clocks:
     * more than 99% for AP13, for `{1..12, 5460}` and for random rows;
     * 87% at clocks 13 and 26 for `26*{1..12} u {339}`.
   * **AMM: tight at one point.** On the lacunary class the Bernstein-box triangle inequality is tight only at the golden point `w = -1`. Per level it is exponentially sharp everywhere, since the ratio lies in `[2/pi, 1]`.
   * **Collatz: spread.** The drift is exactly an AM–GM defect, `log(2/sqrt3)` per step (`q = 3` is the unique odd `q` with arithmetic mean 1). Its equality set has dimension 0.95 (PROVED), so it cannot be excised; choice (E-SCC) is what localizes it.
7. **Levers (section 4).** The LRC lever was tested here and works. Removing the pack-clock neighbourhoods and applying degree-5 Bonferroni on the complement gives these values, each a rigorous lower bound for `P(M=0)` and an exact loneliness certificate:

   | row | raw BONF5 | excised BONF5 | exact `P(M=0)` |
   |---|---|---|---|
   | `{1..12, 5460}` | −6.73 | +0.029235 | 0.029235 |
   | `26*{1..12} u {339}` | −6.73 | +0.0078 | 0.0292 |
   | random 13-sets | about −1.6 to −2.0 | within 0–19% of `P(M=0)` | — |

   For `{1..12, 5460}` the excised bound equals `P(M=0)` exactly. For the AMM and Collatz levers see sections 4.2 and 4.3.

---

## 1. Bridge table

The notation is `source -> target`, then the map, the preserved predicate (PRES), the destroyed information (DESTR), the sidecar, the cheapest decisive test, its result, and the label.

| # | bridge | map / PRES / DESTR / sidecar | test and result (`.out` section) | label |
|---|---|---|---|---|
| B1 | **Product-formula engine.** AMM THM-4467 (Pólya + Fatou–Gauss) -> Collatz Theorem R, D, Y and theta H1 | **Map:** "integer `N != 0`, small at one place => `N = 0`". **PRES:** the vanishing mechanism. Both closing integers are Hankel determinants: of `phi` (AMM), and of the tails of `Theta` (H1). **DESTR:** which place is analytic (`inf` against 2); function against value. **Sidecar:** the continuation domain `W_gamma` (AMM); the valuation law `v_2(T_n) = L n(n+1)(2n-1)/2` (H1). | AD1: the function forms of the Bernstein series have `sum_v log R_v <= 0`, with `= 0` exactly when a density exists (six words, table). AMM's gain `Lambda(gamma) > 0` is continuation. The Hankel valuation laws grow like `n^2` (capacity) against `n^3` (q-series Vandermonde). | **REAL** (one engine, exchanged places). Function-level capacity is **refuted** as a Collatz route. |
| B1' | LRC's integrality engine | THM-2052 (pigeonhole: `(Q+1)^s > QB+1` forces short relations), THM-4009 (Banaszczyk transference), and the Sungkawichai–Trakulthongchai finite-field shortcut, which works for prime `k+1` and fails at the composite `14 = 2*7`. This is the Siegel half; LRC has no Liouville half. | Literature and canon read. | **ANALOGY** |
| B2 | **The Z/2 fold.** AMM `p <-> 1-p` and quotient `w = p(1-p)` -> LRC half-translate `x -> x+1/2` and clock-two quotient `y = 2x` -> Collatz parity split | **Map:** `p = sin^2(pi x)`: the lifts `y/2`, `(y+1)/2` go to `p`, `1-p`; `4w = sin^2(pi y)`; Catalan branch `z(w) = sin^2(pi y/2)`. 2-adically, `L(p) = p(1-p)/2` is isometrically conjugate to `T`. **PRES:** the Z/2 fibre structure. **DESTR:** AMM's complex continuation, which becomes imaginary time (the golden point is `x = i asinh(phi^-1/2)/pi`); and the affine arithmetic (`L`-orbits of rationals have doubling heights, `T`-orbits cycle). | AD3: identities to `4.7e-14`; isometry on 20,000 pairs; heights `1.6 -> 101` bits under `L`, while `1/3, 5/7, -1/5` cycle under `T`. | **REAL** as an identity; **ANALOGY** as a transfer |
| B3 | **The golden ratio.** Theorem S / Collatz Lemma P <-> AMM THM-3027 / Theorem B | All four are eigen-data of `[[1,1],[1,0]]` (PHI1). **PRES:** only the algebra. | PHI2: four knobs, four responses. Section 2. | **NUMEROLOGY** across programs; **REAL** within AMM |
| B4 | **Dyadic structure.** AMM (Lucas/Catalan parity), Collatz (2-adic parity vectors), LRC (clock two; `14 = 2*7` lift wall) | Reduce mod 2 and read the GF(2) Hankel profile. | DF5, AD2: AMM has a perfect profile (all Hankel determinants `+-1` in range). Squares and cubes: jump-size spectra match i.i.d. bits (`n = 6000`). | **ANALOGY** (opposite GF(2) behaviour) |
| B5 | **S596 two-block, in its HYP-9122 finite form.** LRC small-clock residual <-> Collatz loops at the record clocks of `log_2 3` | **PRES:** "a finite covering statement at each clock". **DESTR:** the LRC clocks are pigeonhole-determined and finite; the Collatz clocks are Diophantine and infinite. | BR1: THM-4447's rule reproduces `{2,3,4}` and the exact signatures of THM-4447 (13); it gives `{2,3,4}` for `r = 3` at every `n` tested and moves with `r` (`{2}` for `r = 2`). BR2: the records `3, 5, 17, 29, 41, 94, ...` match HYP-9122. | **ANALOGY** (shape); "2 vs 3" is **NUMEROLOGY**; the Baker import has **no target** |
| B6 | **Mechanical words.** AMM ARCH profile `floor(gamma*(m+k))`, Collatz Sturmian parity words and Ostrowski splicing, LRC THM-536/778 (Sturmian seven-sector walk, centred Beatty ranks) | Shared toolbox (three-distance, Ostrowski, Beatty). The slopes are log-ratios in both AMM and Collatz. | BR3: `gamma* = [0;1,1,2,19,2,9,...]`; the AMM block optima `14/9, 25/16, 53/34, 83/53, 157/100` are not convergents of `C_*`. No Diophantine obstruction is attached to `gamma*`. | **ANALOGY** |
| B7 | **Resonance localization of the discarded structure** (section 3) | LRC: pack clocks, which are roots of unity `e^(2 pi i a/q)`. AMM: `w = -1` and the cyclotomic handoff states `Phi_2^2 Phi_4 Phi_8`. Collatz theta: the KRVZ cyclotomic divisors. | DF1, DF3 and the theta lane. | **ANALOGY**, with a REAL instance in each program |
| B8 | **Two-point integrality.** AMM Theorem A's fold (Green interaction `g_U(0,1)`, which takes the reach from 0.0955 to 0.3775) <-> Collatz Proposition T (two places) <-> LRC's two-lift failure | "Two" enters as two points, two places and two lifts respectively. | No shared predicate was found. | **ANALOGY** |

**What the bridges say about the owner's question.**
* AMM and Collatz share a real engine. It is literally the same Hankel-determinant argument with the analytic place exchanged. The theta lane's H1 was, in effect, THM-4467's engine transplanted to `Q_2`.
* LRC and Collatz share only a shape: easy by density, hard on a thin residual. In the current LRC frontier that residual is finite pigeonhole data together with a linear relation lattice. It is not an exponential Diophantine problem.
* LRC and AMM share an exact Z/2 fold and nothing downstream of it.

---

## 2. The golden ratio: four occurrences, four mechanisms

PHI1 (exact, sympy). Let `M = [[1,1],[1,0]]`.

| occurrence | equation | eigen-datum | origin (mechanism) | sharp for | beaten by |
|---|---|---|---|---|---|
| Collatz Theorem S: `mu(mu-1) < phi` | `limsup q_(n+1)/q_n >= phi` | spectral radius of the CF matrix `[[a,1],[1,0]]` at `a = 1` | Hurwitz/Fibonacci slowest denominator growth, inside a one-scale minimax over the defect position `j0` | nothing: it is a worst case of the split | Theorem D + Bugeaud–Kim, `eta < 2.50994` (a `sqrt10` constant) |
| Collatz Lemma P: `mu_bar < phi` | `max_y (y^2+2y)/(y^2+1) = phi` at `y = (n+m)/n = phi` | top eigenvalue of `[[0,1],[1,1]]` (a Rayleigh quotient; value = argmax since `phi = 1 + 1/phi`) | the valuation-versus-height quadratic-form pencil of ONE Padé family | the one-parameter Zudilin family | H1 `7/4`, H2 `28/11`, H3 `2.878` (theta lane) |
| AMM THM-3027 / THM-3009 | `(1-tau)^2 = tau` | `tau* = phi^-2` = smaller eigenvalue of `M^2` (Zudilin's `(3-sqrt5)/2` is the same number) | double root (tangency) of the entropy rate function | separately balanced blocks | super-blocks, about 1.570 (NUMERICAL) |
| AMM Theorem B / Long's point | `p(1-p) = -1` | `p, 1-p` = the eigenvalues `-1/phi, phi` of `M` (trace 1 = the fairness involution, det -1 = the natural boundary `|w| = 1`); `S(p) = sqrt5 = lambda_max - lambda_min` | the continuation domain first touches the lacunary natural boundary | same class | same |

**Co-variation (PHI2).** Each occurrence has its own knob, and they respond differently.
* **AMM capacity boost `b`.** Our independent re-solution of the THM-3027 rate problem gives `gamma*(b) = log phi/log((b+phi)/phi)` to 7 digits for `b = 2..5`, with `tau* = phi^-2` for every `b`. So the golden value is **universal** under this knob.
* **AMM alphabet `q` (THM-3009).** Here `delta^q = (1-delta)^(q-1)`, so the golden value is **special** to `q = 2`.
* **Theorem S.** The threshold is set by the **word class** through `L`: `mu_S = 1.8668, 2, 2.132, 2.303`. The map does not move it.
* **Collatz Lemma P.** The threshold is set by the **construction**: 3/2 at `x = 0` and `phi` at `x = 1/phi`, for every `rho`.

No knob moves an AMM `phi` and a Collatz `phi` together. The appearance of `phi` in both programs is the algebra of the smallest hyperbolic element of `GL_2(Z)`, which is reached whenever a one-parameter two-scale extremum is taken.

Inside AMM, THM-3027 (tangency), THM-3009 (ARCH), Theorem B (`w = -1`) and Long's evaluation point `x = -phi^-2` are one saddle. PHI1 records the supporting checks:
* `1+x = 1/phi` and `1+2x = phi^-3` at Long's point;
* the Catalan S-fraction at `w = -1` has convergents `-F_n/F_(n+1)`.

**The Theorem S golden threshold is a single-scale artifact (PHI3).** On Sturmian words of six slopes and four intercepts (`N = 30000`), we compared the prefix repetition exponent delivered at each scale with the single-scale minimax `mu_S(q_(n+1)/q_n)`.
* The delivered exponent exceeds the minimax at every scale, by at least `+0.118` overall and by at least `+0.557` on golden-slope words.
* The characteristic golden word delivers `2.618 = 1+phi` at every scale (Berthé–Holton–Zamboni, as cited in the hard-class note).
* Single windows can dip below Bugeaud–Kim's 2.50994 (for example 2.442), because a window maximum is not a limsup.

So the golden number in Theorem S is the price of treating `j0` and `q_(n+1)/q_n` as independent adversaries. The Ostrowski coupling in Corollary D1 removes it.

**LRC(14) (PHI4).**
* The LRC(14) section of CURRENT-FRONTIER contains no golden or Fibonacci token.
* Among the 279 theorem files named `lrc`/`lonely`, none mentions golden, Fibonacci or `sqrt5`. Every occurrence of "phi" there is a variable name.
* Of the 7 Fibonacci-titled carrier theorems, 3 state in their headers that no LRC transfer follows.

---

## 3. The owner's principle: defect ledger

The principle: every Cauchy–Schwarz, AM–GM, Jensen, triangle or union bound has an exact defect that can be computed from the values being averaged. A defect that is large *and* structured is a lever.

### 3.1 The ledger

| program | inequality (where) | exact defect identity | value on the extremal object | verdict |
|---|---|---|---|---|
| AMM | box majorant `|W_m(p)| <= S(p)^(d_m)` (THM-4467 step 1), per level | the zonotope support function: `max |E_m| / S^d = max_theta E|cos(theta - arg z_k)|` | in `[2/pi, 1]` for `d` up to 2000 at every tested point; `= 1` at `p = -1/phi` (DF1a) | **exponentially sharp**: no per-level lever |
| AMM | the same, on the lacunary (golden) class, per block | `delta(p) = gamma log S - log|1-p| = -log|p(1-p)|` on `|p| S^gamma = 1` | `min delta = 0` exactly at `arg p = pi` when `gamma = gamma*`; negative below `gamma*` (`-0.086` at 1/2, `-0.204` at `gamma_1`); `< 0.05` on 27.5% of the boundary angles at `gamma*` (DF1b) | **tight at one point (`w = -1`)**; sub-golden schemes must create a zero there (they do: for example the `N = 8` optimal ratio-4 handoff state `-Phi_2^2 Phi_4 Phi_8`, and the golden-zero families) |
| AMM | Pólya's Hankel bound (Hadamard/Fekete) | Gram–Schmidt: `det^2 = prod ||a_i^perp||^2` | `|H_n| = 1` against a Hadamard bound of `2^40.5` at `n = 100` (AD2) | huge but **useless**: Pólya needs only `|H| >= 1`, and capacity is sharp for general series |
| AMM | Liouville at `p = 2` (adelic) | product formula; `v_2(H_n(phi))` | `= 0` for every `n` (perfect profile), for every complement-symmetric fair extractor | **inert**: a decisive negative for a Bertrandias lever |
| AMM | one-point against two-point Pólya (Theorem A0 -> A) | Green interaction `g_U(0,1)` | reach `0.0955 -> 0.3775` | **harvested** (THM-4467) |
| AMM | THM-3027 Stirling (max term against sum) | `<= log R` | subexponential | none |
| Collatz | drift (Haar AM against GM) | `log(AM/GM)`; AM `(q+1)/4`, GM `sqrt q/2` | `q = 3`: AM = 1, defect `log(2/sqrt3) = 0.14384`/step, i.e. the whole descent; equality set = words with empirical GM `>= 1`, of dimension `h(log_3 2) = 0.94996` (PROVED, ladder); identical on `3x-1` (DF2) | **large but spread and SHEET-blind**: cannot be excised |
| Collatz | Theorem S's split (`max(g_A, g_B)` minimax in `j0`) | realized `e_n - mu_S(q_(n+1)/q_n)` | `>= +0.557` on golden words (PHI3) | **harvested** (Theorem D + Bugeaud–Kim) |
| Collatz | Lemma G Liouville, `|N|_2 >= 1/|N|` | the forced odd divisor (product formula) | squares: KRVZ cyclotomic divisors, exponent ladder `1.4265 -> 0.951 -> 0.882 -> 0.815 -> 0.560 -> 0.496`; cubes: **zero** (theta lane D1, D2, K7) | **harvested** (squares); **none** (cubes) |
| Collatz | Lemma H' heights (triangle on `R_z`) | canonical gcd | at most 28.6 bits residual after canonicalization (theta lane D2) | **sharp** |
| Collatz | Q2 budget `exp(0.1144 D)` | worst chain | PROVED sharp | **sharp** |
| Collatz | Pólya–Bertrandias, function level | `sum_v log R_v` | `= 0` exactly (AD1) | **critical** |
| LRC | union bound / Bonferroni (THM-2051 BONF5) | `P(M=0) - BONF_j = E[C(M-1,j) 1{M>=1}]`, `j` odd (checked exactly, `j = 3, 5`) | raw BONF5 = `-9.72` (AP13), `-6.73` (both frontier controls), `-1.56` to `-2.00` (random); the defect sits more than 99% within `|t - a/q| < 1/(28q)` of the pack clocks, or 87% at `q = 13, 26` for `26*{1..12} u {339}` (DF3) | **large, structured, localized: LEVER** (section 4.1) |
| LRC | pair energy against union mass (THM-4449 (13)) | `E(T) - mu(F_T) = E[(|K_1||K_2| - 1) 1_F]` | `16/847` at `(1,11,121)` and `8/1449` at `(1,9,23)`, both recomputed; `0` at `(1,7,11)`; 21% of `E` at `(1,9,23,37)` (DF4) | **harvested** at three tails (sharp caps); grows with the tail count |
| LRC | fibre union bound (THM-4447 (6)) | overlaps of the killed label sets | zero exactly on the hostile partitions (THM-4447 (9)) | the residual clocks `{2,3,4}` are **exactly where this defect can vanish** |
| LRC | second moment / Paley–Zygmund (THM-660/661, on `W`) | `P - PZ = P CV^2/(1+CV^2)` | degree 2 fails at `k = 8,9,10` and degree 4 clears them (THM-661) | **harvested** |
| LRC | integrality floor of `Var(M)` | `Var(M) >= {m}(1-{m}) = 6/49` | actual 2.81 (AP13), 1.63–1.90 (random) (DF3) | not a lever alone (the covering slack, T1535) |

### 3.2 Where the discarded structure lives

In each program the defect is concentrated on a resonance set, and the program's architecture treats that set separately.

| program | resonance set | how the architecture treats it |
|---|---|---|
| LRC | pack clocks (rational times `a/q`, where `q` divides the gcd of a large sub-pack) | THM-4447 closes clocks `>= 5`; clocks `{2,3,4}` remain |
| AMM | the golden point `w = -1`, where the box inequality is tight for the lacunary class | the Corollary in the uniform-frontier note: sub-golden schemes must be analytic at `w = -1` |
| Collatz, relaxed | `-1` and `1/2` | the E-SCC escape lemmas |
| Collatz, proper | the drift's defect set has dimension 0.95 | none: no finite family of resonance neighbourhoods covers it |

This is the owner's principle in structural form. **A defect is a lever when it localizes.** It localizes in LRC and in AMM. In Collatz it does not. Choice is exactly the operation that localizes it, which is why the relaxation "almost closes".

---

## 4. Levers, one per program

### 4.1 LRC(14): clock-excised Bonferroni (TESTED: works on the test rows)

The pointwise inequality `1{M=0} >= sum_(k<=5) (-1)^k C(M,k)` holds at every `t`. So its integral over any set `E` bounds `P(M=0 and t in E)` from below. Take `E` to be the complement of the neighbourhoods `|t - a/q| < c/(14q)` of the row's pack clocks. The DF3 results (exact rational sweep):

| row | raw BONF5 | clock set, `c` | excised BONF5 | `P(M=0)` |
|---|---|---|---|---|
| AP13 (tight) | `-9.72` | `{1,2}`, `c = 1` | `0` | `0` (a tight row cannot be certified) |
| `{1..12, 5460}` | `-6.73` | `{1,2}`, `c = 1` | **`+0.029235`** | `0.029235` (sharp) |
| `26*{1..12} u {339}` | `-6.73` | `{1,2}` or `{1,2,3,4}` | still `-5.4` or `-4.9` | — |
| `26*{1..12} u {339}` | `-6.73` | `{1,2,13,26}`, `c = 1/2` | **`+0.0078`** | `0.0292` |
| random 13-sets from `[1,80]` | `-1.56` to `-2.00` | `{1,2}` | `+0.0867`, `+0.1062`, `+0.1280` | `0.107`, `0.106`, `0.140` |

Every positive entry is an exact loneliness certificate. The lever is the owner's principle applied literally: the union-bound defect is `E[C(M-1,5) 1{M>=1}]`, it is carried by the clocks, so excise them and handle them by the clock-capacity theorems (THM-4447). This matches the corpus architecture and quantifies why that architecture is right.

**Next decisive test (not run).** For a proof it is not enough to do this row by row.
* **What is needed:** a uniform lower bound for the excised BONF5 in terms of the relation code. Fourier-expand the degree-`<= 5` terms restricted to a smooth version of `E`, and show that clock-induced relations are exactly the ones the excision kills (a THM-2051/2052 upgrade).
* **The finite test:** run the excised BONF5 on the THM-2923 terminal rows and on the THM-4447 (9) hostile phases lifted to 13-speed rows.
  * If it stays positive off the tight AP family, the lever is live.
  * If it fails on a row that is not clock-tight, record that row as the obstruction.

### 4.2 AMM 12592: the golden-point zero (tests run: per-level sharpness, 2-adic inertness)

**Findings, all computed.** Two of the tempting levers are closed:
* **The per-level box bound.** It is exponentially sharp (DF1a), so there is no per-level lever.
* **The prime 2.** It is inert: all Hankel determinants of the parity skeleton are units (AD2, DF5), so there is no Bertrandias lever at 2.

So any uniform improvement past `gamma_1 = 0.3775` must come from **cross-level phase structure**. DF1b shows where that structure is binding for the lacunary class: exactly at `w = -1`, the point where every sub-golden construction must, and does, place cyclotomic zeros.

**Proposed test (not run).** For ratio-4 super-blocks, run CP-SAT for the minimal block deadline ratio with the handoff state constrained to vanish at `w = -1` to order at least `nu`. Do this for `nu = 0..4` at `N = 8, 16`, reusing the uniform-frontier CP-SAT model.
* **A positive deadline cost per unit order of vanishing** would make "sub-golden implies zeros at `-1`" quantitative. That would give a lower bound for the super-block class (HYP-9128/9129).
* **Zero cost** would mean the golden-point lever is empty, and the gap `[1.3775, 1.570]` is a sup-norm phenomenon away from `w = -1`.

### 4.3 Collatz: no unharvested inequality lever; the lever is where the defect localizes

The theta-beyond-phi lane already executed the owner's principle on the value side. It proved HYP-9131 and moved the height exponent from 1.4265 to 0.496. It also proved that for cubes the Liouville/Hankel defect is zero (K7).

This lane adds three things:
* the drift is the AM–GM defect `log(2/sqrt3)` exactly, and it is sheet-blind;
* its equality set has dimension 0.95 and cannot be excised;
* Theorem S's split was a single-scale artifact.

So the Collatz analogue of the LRC lever is not another inequality refinement. It is passage to a relaxation in which the defect localizes. In E-SCC Q2 it localizes on one link type, the `1/2`-transfer. The lever is then the exact escape price there, `2c(k-1)/3` in `(1,2)`: the calculable defect of the see-saw.

**Proposed test (not run).** Tabulate the per-digit cost of the canonical `1/2`-transfer chains by precision `k` against the Ostrowski digits of `log_2 3`. Use the endgame lane's chain data (no new search).
* **If the sharp rate 0.1144 is attained only at a sparse set of precisions**, a resonance-excised budget would bound the chains that avoid them. HYP-9124 would then reduce to those resonant links.
* **If the worst links are dense**, this route is closed. The missing mechanism is then only the HARD-class transversality of foundry v4 (for cubes, HC-TH5: a determinant with a quartic forced divisor).

---

## 5. Controls, scope and what is not claimed

* **Not claimed:**
  * no progress on LRC(14), Collatz, PC or `C*`;
  * no new theorem beyond elementary identities.
* **The LRC certificates in 4.1 are row-level.** The rows are already known to be lonely; the point is the defect anatomy.
* **Finite-window Diophantine estimates (PHI3)** are not limsups, and are used only as per-scale lower evidence.
* **The Hankel `+-1` statement is FINITE-EXACT** (`n <= 100` exactly, `n <= 300` mod 2).
  * The all-`n` parity follows from Wang–Massey (CITED from memory, not re-read).
  * We did not search the literature for the integer `+-1` pattern; no priority is claimed.
* **The logistic fold (B2)** is an exact identity, but it is not a transfer of results.
* **The THM-4447 generalization (BR1)** is a parametric re-run of its capacity rule. It assumes the relevant LRC base case for other `n`, and it is used only to show that the open clocks track the tail count.
* **Negative results worth keeping:**
  * the adelic-capacity route for PC is exactly critical;
  * the prime 2 is inert for AMM;
  * the per-level Bernstein box is exponentially sharp;
  * the S596 Baker import has no target in the current LRC residual.

---

## 6. Reproduction

```bash
python3 04-computation/experiments/procgen_bridges_20260923_run.py      # all four parts, sequentially
# or: ..._phi.py, ..._adelic.py, ..._defects.py, ..._twoblock.py  (each accepts --quick)
```

* **Resources:** about 45 s total; one process at a time; each part peaks below 100 MB RSS.
* **Output:** deterministic (fixed seeds), with script SHA-256 hashes in the header of [`procgen_bridges_20260923.out`](procgen_bridges_20260923.out).
* **Exactness:** LRC measures are exact `Fraction` sweeps. The golden identities use sympy. The Hankel determinants are exact (Bareiss) and GF(2).

## Sources

* **Repository notes:**
  * [amm12592_procgen_20260923_uniform_frontier](amm12592_procgen_20260923_uniform_frontier.md) (Theorem A/B, Lemma P (parity), Lemma S, CP-SAT table)
  * [collatz_procgen_20260922_transversality_foundry](collatz_procgen_20260922_transversality_foundry.md) (Theorem R, S, Lemma J/H)
  * [collatz_procgen_20260922_hard_class](collatz_procgen_20260922_hard_class.md) (Theorem D, Corollary D1, Lemma P (Tschakaloff), Theorem Y)
  * [collatz_procgen_20260923_theta_beyond_phi](collatz_procgen_20260923_theta_beyond_phi.md) (H1–H3, K7, defect ladder)
  * [collatz_procgen_20260922_synthesis](collatz_procgen_20260922_synthesis.md) (§2c, §2d, §4)
  * [S596 reflection](../../07-reflections/lrc-collatz-the-same-two-block-question-s596.md)
  * HYP-9122, HYP-9131
* **Canon:** THM-4467, THM-3027, THM-3009, THM-3342, THM-4447, THM-4449, THM-2051, THM-2052, THM-4009, THM-661, THM-663, THM-4107 (a prior "2:3" LRC-blind carrier).
* **External (CITED via the corpus; not re-read this session):**
  * Sungkawichai–Trakulthongchai (LRC through 13 runners);
  * Bugeaud–Kim (`Dio >= 2.50994`);
  * Berthé–Holton–Zamboni;
  * Wang–Massey (perfect linear-complexity profiles);
  * Pólya (1928), Fatou (1906), Bézivin (1998), KRVZ.
