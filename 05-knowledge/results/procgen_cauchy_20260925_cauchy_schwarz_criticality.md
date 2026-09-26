# Cauchy–Schwarz criticality: the square loss in the provability price is the AM-fairness `g(2) = 1`; a proved exponent 0.1445 (was 0.7737), the conditional 0.1001, and a harmonic bound that sees past moments

**Status.**
- **PROVED** (hand proofs below; every identity is also checked in exact arithmetic by the scripts):
  1. **Theorem CS.** Every member of `P_L` has upper flip density `δ̄ >= ρ_L^2 / M2(L)`, where `M2(L) = E_Haar[W_L^2]` and `W_L(v) = sum_(k<L) sum_(T^k n = v) 3^(a(n,k))/2^k` (THM-4475 §3 bounds this weight by its maximum). The depth-`k` backward tree of `v`, with its weights, is a function of `v mod 3^k` exactly, so `W_L` is a function of `v mod 3^(L-1)`: for `v >= 1` the integer backward tree is the formal 3-adic tree.
  2. **Theorem T (transfer operator; criticality).**
     - `Z_k := sum_(T^k n = v) 3^a/2^k` equals `L^k 1`, the density of the forward 3-adic walk `X -> X/2` or `(3X+1)/2` started from Haar measure.
     - For all `f, g` in `L^2(Z_3)`: `<Lf, Lg> = <f, g> + (1/4)(<f∘τ, g> + <f, g∘τ>)`, with `τ(u) = 3u + 1`.
     - The diagonal coefficient is `g(2) = 1`; for `qx+1` it is `g_q(2) = (1+q)/4`.
     - Hence `||Z_(k+1)||^2 = ||Z_k||^2 + γ_k/2` exactly. All of the sibling dependence sits in the cross term `γ_k = <Z_k∘τ, Z_k>`.
     - Ladder expansion: `γ_k = 3 sum_j 4^-j <Z_(k-2j)∘S^j, Z_k> + 4^(-⌊k/2⌋) c_k`, with `S(x) = 4x + 1`.
  3. **Theorem M (computer-assisted).**
     - `||Z_k||^2 <= 253.82 · 1.0312629^(k-14)` for all `k >= 14`. This comes from a certified super-eigenvector of a class-mass majorant mod `3^15`.
     - Hence `δ_L >= 2^(-(0.1445+o(1))L)`, improving THM-4475's `0.7737`.
     - Explicit: `δ_40 >= 2.65·10^-9` (THM-4475: `6.1·10^-12`) and `δ_100 >= 2.1·10^-13` (THM-4475: `2.1·10^-26`).
  4. **Theorem H (the moment-method limit).**
     - (a) For every Hölder exponent `p != 2` the Hölder bound has exponent `> 2(1-h) = 0.1001`.
     - (b) More strongly, the best bound obtainable from the law of `W_L` alone is `δ*_L(ρ_L) = min{μ(B) : ∫_B W_L >= ρ_L}`, and `δ*_L(ρ_L) = O(ρ_L^2 / L)`. So no distribution-only argument beats the exponent `2(1-h)`.
  5. **AM-fairness is the quadratic loss.**
     - `g_q(s) = E_fwd[w^(s-1)]`, where the forward step factor `w` is `1/2` or `q/2`.
     - The identity `(3/4) log2 3 - 1 = 1 - h(3/4)`, which makes 4(b) sharp, is `g(2) = 1`, i.e. THM-4470's `(q+1)/4 = 1`.
     - The loss factor of any moment method is `κ/(κ-1)`, where `κ` is the non-trivial root of `g = 1`. `κ = 2` iff `E_fwd[w] = 1`.
  6. **Theorem P (the harmonic bound).**
     - Every member of `P_L` has `δ̄ >= H_L := E[1_Bad_L(ω) / max_(j<L) W_L(X_j)]`. Here `ω` is a uniform parity word and `X_j = T^j(n) mod 3^(L-1)` is the walk driven by `ω` from a Haar start.
     - `H_L >= max(ρ_L^2/M2(L), ρ_L/max W_L)`.
  7. **The fixed point.**
     - The forward walk contracts the Wasserstein-1 distance on `Z_3` by `2/3`. Its unique fixed point is the Syracuse law `π`, the law of `sum_t 3^(t-1) 2^-(G_1+...+G_t)` with `G_i` i.i.d. and `P(G = g) = 2^-g`.
     - Exact: `||π||_n^2 = sum_(lev ξ <= n-1) P_(1/4)(ξ) |ν̂(ξ)|^2`, with the Poisson kernel `P_(1/4)(ξ) = 15/(17 - 8 cos 2π{ξ})` in the coordinates `Λ = log_4(1+3Y)`.
     - Equivalently `||π||_(n+1)^2 = sum_(j in Z) 4^-|j| C_j^(n)`, with the ladder autocorrelations `C_j^(n)` of §7.
     - The level energies `E_n = ||π||_n^2 - ||π||_(n-1)^2` satisfy `E_1 = 2/3` and `E_2 = 10/21`.
     - The `L^2` mass on the classes 1 and 2 mod 3 is exactly `1 : 4` at every level.
  8. **Controls.**
     - **SHEET.** `Z^-_k(v) = Z_k(-v)`, so `M2`, `ρ_L`, `H_L` and every bound here coincide on the `3n−1` sheet.
     - **DRIFT (`5x+1`).**
       - `||Z^(5)_k||^2 >= (3/2)^k`.
       - The Cauchy–Schwarz bound is `<= β_L^2 (2/3)^(L-1)`, and even the harmonic bound is `<= 0.9905^(L-1)`.
       - Both tend to 0 although `β_L -> 0.176`. The observed `5x+1` price floor is invisible to single-hit covering.
  9. **The mean-field model.**
     - The i.i.d. smoothing transform `W = W'/2 + (3/2) B W''` (with `B ~ Bernoulli(1/3)`) contracts the Zolotarev metric `ζ_s` with constant `g(s)`. It is a strict contraction exactly for `1 < s < 2`; the constant is `1` at `s = 2`.
     - `E[Z_k^2] = 1 + k/2`.
     - In this model the harmonic bound is `>= ρ_L 2^(-O(√L))`, so HYP-9137 holds there. The asymptotic step uses the classical meander scaling (UNVERIFIED citation); the strip data are FINITE-EXACT to `L = 800`.
- **FINITE-EXACT.**
  - `M2(L)` as exact rationals for `L <= 16`; for example `M2(8) = 501975/4096`, `M2(16) = 762.553`.
  - `||Z_k||^2` and `γ_k` for `k <= 15`. The increment is `0.3565` at `k = 15` and still slowly decreasing; independent siblings would give `0.5`.
  - `H_L` exactly for `L <= 13`, and by Monte Carlo to `L = 16`:
    - `H_13 = 3.176·10^-4`, which is 8.8 times THM-4475's bound and 69 times the Cauchy–Schwarz bound;
    - `H_L ≈ 3.2 ρ_L / M2(L)` for `11 <= L <= 16`.
  - Bad orbits carry `2.0–2.6` times the average hub weight, where Cauchy–Schwarz must allow `1/ρ_L = 31`.
  - `δ*_L ≈ 30 ρ_L^2/M2(L)` for `L <= 16`.
  - Syracuse level energies to `n = 15` (`E_15 = 0.4708`, slowly rising towards `≈ 0.476`); ladder correlations; tails `x^2 P(dπ/dμ > x) ≈ 0.9` on `[6, 32]`.
  - The certified `θ_r` for `r <= 15`.
- **CITED.**
  - In-repo: [THM-4475](../../01-canon/theorems/THM-4475-price-of-provable-descent-tends-to-zero.md) (§3, and (B) `ρ_L >= 2^(-(1-h)L)/poly(L)`), [THM-4470](../../01-canon/theorems/THM-4470-collatz-pairing-ladder-am-fair-and-defect-blind.md) §1, and the [inverse-tree note](collatz_procgen_20260924_inverse_tree_mod192.md), Props. 9–10.
  - arXiv abstracts read this session: Tao 1909.03562; Jelenković–Olvera-Cravioto 1012.2165; Alsmeyer–Biggins–Meiners 0906.3133; Buraczewski–Damek–Zienkiewicz 1504.03144.
  - Everything else in §9 is marked UNVERIFIED.
- **OPEN.**
  - **(O1)** `M2(L) = poly(L)`. Equivalent in practice: the cross terms `γ_k` stay bounded, or the level energies `E_n` of `π` stay bounded. With (O1), Theorem CS gives the coordinator's exponent `0.1001`.
  - **(O2)** HYP-9137 through Theorem P. It suffices that typical bad orbits are hub-biased by at most `2^(o(L))` (condition (SD), §5).
  - **(O3)** Absolute continuity of `π`, and the tail `P(dπ/dμ > x) ~ c x^-2`.
- **CORRECTED.**
  - The heuristic "second moments grow only polynomially because `g(2) = 1`" treats the two children of a node as independent. In the tree they are `3x+1` and `x` (with `x = E(u)`), deterministically linked. `g(2) = 1` is exactly the neutral diagonal; the growth rate is decided by the cross term, which is not controlled by `g`.
  - So `0.1001` is conditional on (O1), and `0.1445` is what is proved.
  - `W = W'/2 + 1_legal (3/2) W''` holds for the 3-adic tree only in the mean. In law, the exact object is the functional equation `p = Lp` on `Z_3` (§7).

Session `collatz-procgen-20260922`, Cauchy–Schwarz lane, 2026-09-25.
- Scripts: `04-computation/experiments/procgen_cauchy_20260925_{moments,syracuse,majorant,harmonic,meanfield,run}.py`.
- Output: [procgen_cauchy_20260925.out](procgen_cauchy_20260925.out). It takes about 60 s; every process stays below 500 MB footprint and 700 MB maximum RSS.
- No HYP or THM file was created.

## 0. Setting and the numbers at a glance

**Notation.** `T(n) = n/2` for even `n` and `(3n+1)/2` for odd `n`. `η := 1 - h = 0.050044`, with `h = h(log_3 2)`.
- `ρ_L = |Bad_L|/2^L` is the density of the classes mod `2^L` whose parity word keeps `3^(a_k) > 2^k` for all `k <= L`. `ρ_L = 2^(-ηL + O(log L))` (THM-4475 (B)).
- `P_L` and `δ_L` are as in THM-4475: the members of the pairing family in which every `n >= 3` descends within `L` steps, and the infimum of their flip densities.

**Backward tree.**
- `D(v) = 2v`, with weight `1/2`.
- `E(v) = (2v-1)/3`, with weight `3/2`, legal iff `v ≡ 2 (mod 3)`.
- `Z_k(v)` is the weighted count of the depth-`k` preimages. `W_L = sum_(k<L) Z_k`, and `M2(L) = E_Haar[W_L^2]`.

**Transfer operator.** `L f(v) = f(2v)/2 + (3/2)[v ≡ 2] f(E v)`, so that `Z_k = L^k 1`.

| lower bound on `δ_L` | exponent `-(1/L) log2` | status |
|---|---|---|
| THM-4475 sup bound `ρ_L / sum g_k` | 0.7737 | PROVED (THM-4475) |
| Cauchy–Schwarz with the certified `θ_15` (Theorems CS + M) | **0.1445** | **PROVED** |
| Cauchy–Schwarz if `M2` is polynomial | 0.1001 `= 2η` | conditional on (O1) |
| any distribution-only argument | `>= 0.1001` | **PROVED** (Theorem H) |
| harmonic bound `H_L` | `H_L ≈ 3.2 ρ_L/M2(L)` on `L <= 16`, i.e. only polynomial loss against `ρ_L` | FINITE-EXACT; asymptotics OPEN (O2) |
| sharp (HYP-9137) | `η = 0.0500` | OPEN; upper half is THM-4475 |

## 1. Theorem CS (PROVED)

**Lemma 1 (what the tree depends on).**
- A backward word with `e` `E`-steps is legal at `v` iff `v` lies in one residue class mod `3^e`. By induction on the `E`-steps: after a prefix with `e'` of them the node is `(2^(k') v - c)/3^(e')`, and the next `E`-step needs this `≡ 2 (mod 3)`, which is one class mod `3^(e'+1)` inside the previous class (cf. inverse-tree note, Prop. 9).
- Hence the depth-`k` tree of `v` (its legal words and their weights `3^e/2^k`) is a function of `v mod 3^k`. This is sharp: `E^k` is legal exactly on `-1 + 3^k Z_3`.
- For an integer `v >= 1`, every legal `E`-preimage `(2v-1)/3` is a positive odd integer, and `D(v) = 2v`. So the integer backward tree is the formal one.
- Consequently `W_L` is `3^(L-1)`-periodic on the integers, and its logarithmic averages are Haar averages.

**Theorem CS.** For every `M` in `P_L` (any threshold), `δ̄(M) >= ρ_L^2 / M2(L)`.

*Proof.* This is THM-4475 §3 with one extra Cauchy–Schwarz step.
- Let `A` be the set of members of flipped pairs, so `d̄(A) = δ̄`.
- For `n in Bad_L` with `n >= n_0`, some `T^k n` with `k < L` lies in `A`; otherwise `M` follows `T` and `n` does not descend. Let `v(n)` be the first such point.
- Given `ε > 0`, there is `n_1` with `v/n <= (1+ε) 3^a/2^k` for `n >= n_1`. With `Y = 2^L X`:

  `sum_(n in Bad_L ∩ [n_1, X]) 1/n <= (1+ε) sum_(v in A, v <= Y) W_L(v)/v <= (1+ε) (sum_(v in A, v <= Y) 1/v)^(1/2) (sum_(v <= Y) W_L(v)^2/v)^(1/2)`.

- The left side is `ρ_L ln X + O_L(1)`, since `Bad_L` is periodic mod `2^L`.
- The last factor is `(M2(L) + o(1)) ln Y`, by Lemma 1.
- `sum_(v in A, v <= Y) 1/v <= (δ̄ + o(1)) ln Y`, since upper logarithmic density is at most upper density.
- Let `X -> ∞` and then `ε -> 0` to get `ρ_L <= (δ̄ M2(L))^(1/2)`. ∎

**Data (FINITE-EXACT, exact rationals).**

| `L` | 8 | 10 | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|---|
| `M2(L)` | 122.55 | 218.95 | 353.92 | 437.66 | 533.21 | 641.27 | 762.55 |
| `M2/L^3` | 0.239 | 0.219 | 0.205 | 0.199 | 0.194 | 0.190 | 0.186 |
| i.i.d.-sibling value | 134.0 | 242.5 | 397.0 | 494.0 | 605.5 | 732.5 | 876.0 |
| `ρ_L^2/M2` | `4.49e-5` | `1.78e-5` | `8.60e-6` | `4.59e-6` | `3.76e-6` | `2.44e-6` | `1.36e-6` |
| THM-4475 `ρ_L/sum g` | `7.47e-4` | `2.28e-4` | `7.34e-5` | `3.61e-5` | `2.18e-5` | `1.17e-5` | `5.76e-6` |

- At these `L` the sup bound is still ahead: the ratio CS/sup rises from 0.06 (`L = 8`) to 0.24 (`L = 16`).
- With the rigorous majorant of §3, the proved CS bound overtakes it at `L ≈ 22`.
- `M2/L^3` still decreases; the data fit `M2 ≈ 0.06 L^3` plus lower-order terms.

## 2. The exact second moment: where `g(2) = 1` enters, and where it does not (PROVED + FINITE-EXACT)

**Theorem T.**
- (i) `L` is the Perron–Frobenius operator of the walk `X -> X/2` or `(3X+1)/2` (probability 1/2 each) on `Z_3`. The map `X -> (3X+1)/2` pushes density `f` to `3 f(E v) [v ≡ 2]`. So `Z_k = L^k 1` is the density of `X_k` from a Haar start, and `∫ Z_k = 1` (checked exactly for `k <= 15`).
- (ii) For real `f, g`:

  `<Lf, Lg> = (1/4 + 3/4) <f, g> + (1/4)(<f∘τ, g> + <f, g∘τ>)`,   `τ(u) = 3u+1`.

  *Proof.* Expand `(Lf)(Lg)`. On `{v ≡ 2}` the substitution `w = E v` has `dv = dw/3` and `2v = 3w + 1`. So `∫_(v≡2) f(2v) g(Ev) dv = <f∘τ, g>/3`, and `∫_(v≡2) f(Ev) g(Ev) dv = <f, g>/3`. ∎
- For `qx+1` (`E_q v = (2v-1)/q`, weight `q/2`) the same computation gives `||L_q f||^2 = ((1+q)/4) ||f||^2 + (1/2) <f∘τ_q, f>`, with `τ_q(u) = qu + 1`. The diagonal coefficient is exactly `g_q(2)`.
- (iii) Hence `||Z_(k+1)||^2 = ||Z_k||^2 + γ_k/2` exactly. `g(2) = 1` makes the diagonal **neutral**; it says nothing about `γ_k`.
- (iv) **Ladder expansion.**
  - `3x+1 ≡ 1 (mod 3)` has only the `D`-child `6x+2`, whose children are `3S(x)+1` and `S(x)`, with `S(x) = 4x+1`. This is the sibling ladder `E D^2 = S E` of the inverse-tree note, Prop. 11.
  - Hence `Z_a(3x+1) = sum_(j=1)^(J) 3·4^-j Z_(a-2j)(S^j x) + 4^-J Z_(a-2J)(3 S^J x + 1)`, with `J = ⌊a/2⌋`.
  - So `γ_k` is a `4^-j`-weighted sum of **ladder correlations** `<Z_(k-2j)∘S^j, Z_k>`. The script checks this exactly for `k <= 12`.

**Why independence fails.**
- The two children of a legal node `u` are `2u = 3x + 1` and `x = E(u)`: the `D`-subtree is a deterministic function of the `E`-subtree's root.
- An i.i.d. cascade would have `γ_k = 1`, i.e. `E Z_k^2 = 1 + k/2`. The exact values are smaller.

| `k` | 3 | 6 | 9 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|---|
| `||Z_k||^2` | 9/4 | 3.5356 | 4.7393 | 5.8401 | 6.1988 | 6.5567 | 6.9132 |
| increment `γ_(k-1)/2` | 0.5 | 0.4067 | 0.3817 | 0.3615 | 0.3587 | 0.3579 | 0.3565 |
| i.i.d. `1 + k/2` | 2.5 | 4 | 5.5 | 7 | 7.5 | 8 | 8.5 |

- The increments decrease slowly (`γ_15 = 0.7127`).
- **Near-martingale property (FINITE-EXACT).** `<Z_(k-1), Z_k>/||Z_(k-1)||^2 = 1 ± 3·10^-4` for `k >= 8`, although no exact martingale structure is available in the digit filtration.

## 3. Theorem M: a certified growth rate (PROVED, computer-assisted)

**The majorant.** For a class `c mod 3^r`, let `u_k(c) = ∫_(c+3^r Z_3) Z_k^2`, and let `U_k` be the sums over the three lifts, at level `r-1`. The squared recursion splits into three terms:
- `∫_c Z_k(2v)^2 = u_k(2c)`;
- `∫_c [v≡2] Z_k(Ev)^2 = U_k(Ec)/3`;
- the cross term, which equals `(1/3) ∫_(Ec + 3^(r-1)Z_3) Z_k(3w+1) Z_k(w) dw`. Cauchy–Schwarz bounds it by `(1/3) (3 u_k(2c) U_k(Ec))^(1/2)`, because `w -> 3w+1` maps `Ec + 3^(r-1) Z_3` onto `2c + 3^r Z_3` and divides measure by 3.

Hence:

`u_(k+1)(c) <= Φ(u_k)(c) := u_k(2c)/4 + [c ≡ 2 mod 3] ( (3/4) U_k(Ec) + (√3/2) (u_k(2c) U_k(Ec))^(1/2) )`.

**Certification.**
- `Φ` is monotone and positively 1-homogeneous. So `e > 0` with `Φ(e) <= θ e` and `u_14 <= C e` give `u_k <= C θ^(k-14) e` for all `k >= 14`.
- Classes `c ≡ 0 (mod 3)` carry exactly `u_k = 4^-k 3^-r`, because `Z_k = 2^-k` there. They receive a tiny `e(c) = ε`, chosen so that `u_14 <= C e` holds with the same `C`.
- `e` comes from power iteration. The ratio `max_c Φ(e)(c)/e(c)` is evaluated in IEEE double precision: a few flops and one correctly rounded square root per entry, relative error `< 10^-14`. It is then inflated by `1 + 10^-9`.

| `r` | 1 | 5 | 10 | 12 | 13 | 14 | **15** |
|---|---|---|---|---|---|---|---|
| certified `θ_r` | 1.34307 | 1.09396 | 1.04981 | 1.04071 | 1.03713 | 1.03402 | **1.0312629** |
| `log2 θ_r` | 0.4255 | 0.1296 | 0.0701 | 0.0576 | 0.0526 | 0.0483 | **0.0444** |

**Result.**
- With `C Σe = 253.82`: `||Z_k||^2 <= 253.82 · 1.0312629^(k-14)`.
- So `M2(L) <= (sum_(k<L) ||Z_k||)^2 = O(1.0312629^L)`.
- With `ρ_L >= 2^(-ηL)/poly(L)`: **`δ_L >= 2^(-(2η + 0.04441)L - O(log L)) = 2^(-(0.1445+o(1))L)`.**
- A second run of the majorant, iterated directly from the exact `u_14`, gives sharper constants for moderate `k`:
  - rigorous `||Z_20||^2 <= 8.700` (the exact trend extrapolates to about `8.69`);
  - explicit bounds `δ_22 >= 2.1·10^-7`, `δ_40 >= 2.65·10^-9`, `δ_100 >= 2.09·10^-13`.

**Limits.**
- `θ_r - 1` falls only by about 8% per level. A class-level Cauchy–Schwarz cannot see that the true cross term stays bounded, `γ_k ≈ 0.71`, while `||Z_k||^2 -> ∞`.
- Reaching `θ = 1` needs (O1).

## 4. The moment-method limit: the square is exactly the Cauchy–Schwarz loss, forced by `g(2) = 1` (PROVED)

**Theorem H.**

(a) **Hölder.**
- For `1 < p < ∞`, Hölder in place of Cauchy–Schwarz gives `δ̄ >= ρ_L^(p/(p-1)) / (E W_L^p)^(1/(p-1))`.
- `(sum w_π)^p >= sum w_π^p` and `μ(C_π) = 3^-e` give, exactly, `E[Z_k^p] >= g(p)^k`, so `E W_L^p >= g(p)^(L-1)`. Also `E W^p >= 1`.
- Hence the exponent is at least `e(p) = [pη + max(0, log2 g(p))]/(p-1)`.
  - For `p < 2`: `e(p) >= pη/(p-1) > 2η`.
  - For `p > 2`: `log2 g` is convex with `(log2 g)'(2) = ln(27/16)/(4 ln 2) = 0.18872 > η`, so `log2 g(p) > η(p-2)` and `e(p) > 2η`.
- Only `p = 2` reaches `2η`, and only if `M2` is subexponential. Values: `e(1.5) >= 0.150`, `e(1.9) >= 0.106`, `e(2) >= 0.1001`, `e(2.1) >= 0.114`, `e(3) >= 0.236`.
- The exact `E[Z_15^p]/g(p)^15` is 2–17 over `p = 1.25..4`: only polynomial factors above the lower bound.

(b) **Any distribution-only argument.**
- The strongest consequence of "each bad orbit meets a flip `v`, which serves log-weight `<= W_L(v)/v`" together with the law of `W_L` is `δ̄ >= δ*_L(ρ_L)`, where `δ*_L(ρ) = min{μ(B) : ∫_B W_L >= ρ}`: the adversary puts `A` on the heaviest classes.
- **Claim:** `δ*_L(ρ_L) = O(ρ_L^2/L)`.
- *Proof.* Let `B` be the union of the legality classes `C_π` of backward words of length `j` with at least `e_0 = ⌈3j/4⌉` `E`-steps.
  - The walk's own reversed last `j` steps form such a word with probability `P_j = P(Bin(j, 1/2) >= e_0)`. So `∫_B Z_k >= P_j` for every `k >= j`, and `∫_B W_L >= (L-j) P_j`.
  - `μ(B) <= sum_(e >= e_0) C(j,e) 3^-e <= P_j / x_j`, with `x_j = 3^(e_0)/2^j`.
  - `P_j x_j >= C(j, e_0) 3^(e_0)/4^j >= c/(j+1)`, because `sum_e C(j,e) 3^e = 4^j` peaks at `e = 3j/4`. This is `g(2)^j = 1`, or the identity `(3/4) log2 3 - 1 = 1 - h(3/4)`.
  - Choose the largest `j` with `(L-j) P_j >= ρ_L`. Then `j ≈ 0.265 L`, `P_(j+1) >= P_j/2`, and `μ(B) <= (j+1) P_j^2/c = O(ρ_L^2/L)`. ∎
- **Data.** `δ*_L ≈ 30 ρ_L^2/M2(L)` for `8 <= L <= 16`: Cauchy–Schwarz is tight up to a constant factor among distribution-only bounds.

**Why the square, and why exactly 3n+1 (the owner's thesis made exact).**
- `g_q(s) = 2^-s (1 + q^(s-1)) = E_fwd[w^(s-1)]`: backward legality has probability `1/(2w)` for a forward factor `w`.
- So:
  - `g(1) = 1` is conservation of mass;
  - `-g'(1) = -E log w = log(2/√3)` is the drift. The geometric mean `√3/2 < 1` governs typical orbits and the tree density `R = 1/log(2/√3) = 6.95212` (Prop. 10);
  - `g(2) = E w = (1+q)/4` is the **arithmetic** mean.
- The non-trivial root `κ` of `g = 1` is the tail exponent of the cascade, and the best moment-type loss is `ρ -> ρ^(κ/(κ-1))`. This is a heuristic for general `q`; for `q = 3` it is Theorem H(b).
- `κ = 2` iff `E w = 1` iff `q = 3`. That is THM-4470's AM-fairness `(q+1)/4 = 1`, the same equation as the pair-sum identity `f(2i-1) + f(2i) = (q+1)i + ...`.
- In Theorem CS the square is lost in the step `E[1_Bad Ŵ] <= E[Ŵ]` (§5). The inequality `ρ^2 <= E[1_Bad Ŵ] E[1_Bad/Ŵ]` is AM–HM for the hub weight `Ŵ` along bad orbits. What Cauchy–Schwarz abstracts away is **which orbit meets which hub**.
- The data show that the abstraction is costly. Bad orbits carry `E[Ŵ | Bad]/E[Ŵ] = 2.0, 2.2, 2.2, 2.4, 2.6` times the average hub weight (`L = 12..16`). Cauchy–Schwarz must allow `1/ρ_16 = 31`.

## 5. Beyond moments: the harmonic bound (Theorem P PROVED; its exponent OPEN)

**Theorem P.** For every `M` in `P_L`:

`δ̄ >= H_L = E[1_Bad(ω) / Ŵ],   Ŵ := max_(j<L) W_L(X_j)`.

*Proof.*
- With `v(n)` as in §1, `Ŵ(n) >= W_L(v(n))`. So

  `sum_(n in Bad_L ∩ [n_1, X]) 1/(n Ŵ(n)) <= sum_(v in A, v <= Y) (1/W_L(v)) sum_(T^k n = v, k<L) 1/n <= (1+ε) sum_(v in A, v <= Y) 1/v`.

- `1_Bad(n)/Ŵ(n)` is periodic mod `2^L 3^(L-1)`: `T^j n mod 3^(L-1) = (3^(a_j) n + c_j) 2^-j` depends only on `n mod 3^(L-1)` and the word.
- By the Chinese remainder theorem its mean over a period is `H_L`, with `X_0 = n mod 3^(L-1)` independent of `ω`. ∎

**Comparison with the other bounds.**
- `H_L >= ρ_L/max W_L` is trivial.
- `H_L >= ρ_L^2/M2(L)`: by Cauchy–Schwarz, `ρ_L^2 <= E[1_Bad Ŵ] H_L`, and `E[1_Bad Ŵ] <= sum_(j<L) E W_L(X_j) = sum_j ∫ Z_j W_L = M2(L)`.

**Data.** Exact for `L <= 13` (DFS over bad prefixes, all `X_0 mod 3^(L-1)`); Monte Carlo with `4·10^5` samples for `L = 14..16`, which reproduces the exact values at `L = 12, 13`.

| `L` | 8 | 10 | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|---|
| `H_L` | `1.658e-3` | `8.34e-4` | `4.872e-4` | `3.176e-4` | `2.728e-4` | `2.016e-4` | `1.356e-4` |
| `H_L/ρ_L` | 0.0223 | 0.0133 | 0.00883 | 0.00709 | 0.00609 | 0.00510 | 0.00420 |
| `H_L M2(L)/ρ_L` | 2.74 | 2.92 | 3.13 | 3.10 | 3.25 | 3.27 | 3.21 |
| `H_L/δ*_L` | 1.21 | 1.64 | 2.05 | 2.25 | 2.51 | 2.77 | 3.02 |

- `H_L` exceeds every distribution-only bound, by a factor that grows with `L`.
- Over the whole range `H_L ≈ 3.2 ρ_L/M2(L)`. That is only a polynomial loss against `ρ_L`, whereas Cauchy–Schwarz loses the factor `ρ_L`.

**What a proof of HYP-9137 through Theorem P needs.** Two Jensen steps give `H_L >= ρ_L 2^(-E_(ω~Bad) log2 E_(X_0)[Ŵ | ω])`. So HYP-9137 follows from:

**(SD) quenched sibling decorrelation:** `E_(ω~Bad_L) [log2 E_(X_0)[max_(j<L) W_L(X_j) | ω]] = o(L)`.

- **Data.** This quantity is `5.59, 6.36, 7.01, 7.36` at `L = 8, 10, 12, 13`, about `log2 M2(L) - 1.4`.
- **The annealed version fails.** Under the conditioning on `Bad_L` the own-path weight `2^(S_j)` (with `S_j = a_j log2 3 - j`) has expectation `>= 2^(c j)`, because rare bad words with heavy prefixes dominate.
  - Exact DP: `E[2^(S_100) | Bad_200] = 5.8·10^4 = 2^(0.158·100)`.
  - The large-deviation heuristic (tilt the prefix to odd frequency `3/4`, then survive from the raised height) gives `≈ 2^(0.14 j)` for `j <= L/2`.
  - So `E[Ŵ | Bad]` grows exponentially, while typical bad words keep `max_j S_j = O(√L)`.
  - The strip counts confirm this: `P(max S <= 2√L | Bad_L) >= 0.998` for `50 <= L <= 800`.
  - The quenched log in (SD) is essential.
- **(SD) holds in the mean-field model.** Replace the sibling subtrees along the orbit by independent cascades. Each sibling carries expected weight `<= L/2`, and the own path from `X_j` back to `X_i` weighs `2^(S_j - S_i)`.
  - So `E[Ŵ | ω] <= L^3 2^(max_j S_j)` on bad words, and `H^mf_L >= ρ_L 2^(-2√L)/L^3 · 0.998` for `50 <= L <= 800` (FINITE-EXACT strips).
  - Asymptotically this follows from meander scaling (UNVERIFIED citation).
- **Conclusion.** In the mean-field model the sharp exponent `η` of HYP-9137 holds. The only obstruction in the true tree is the 3-adic correlation between an orbit's sibling subtrees and its own (2-adic) badness.
- **The structure beyond moments is therefore this:** the pointwise covering relation (Theorem P), plus quenched independence of 2-adic badness and 3-adic hub weight along orbits. "Hubs are rare" alone is exactly what the moments see, and it gives `2η` (Theorem H).

## 6. Controls

**SHEET (PROVED).**
- `U(n) = (3n-1)/2 = -T(-n)`: its backward `E`-child is `(2v+1)/3` (legal iff `v ≡ 1 mod 3`), and `Z^-_k(v) = Z_k(-v)`. The script checks this for all `v mod 3^k`, `k <= 12`.
- Negation preserves Haar measure, `ρ_L` (Terras bijection, THM-4475 (D)) and the harmonic functional.
- So every bound here holds verbatim for the `3n−1` pairing family: the Cauchy–Schwarz and harmonic bounds are sheet-blind, as THM-4475's construction is.

**DRIFT, `5x+1` (PROVED + FINITE-EXACT).**
- `g_5(2) = 3/2` is the exact diagonal coefficient, so `||Z^(5)_k||^2 >= (3/2)^k`. Exactly, the ratio is `1.535 · 1.5^k` at `k = 9`.
- Hence the Cauchy–Schwarz bound is `<= β_L^2 (2/3)^(L-1) -> 0`, while `β_L -> 0.1761`.
- Moreover `g_5'(1) = (1/2) log(5/4) > 0`: typical orbits expand. The own path gives `Ŵ >= 5^(a_(L-1))/2^(L-1)`, so `H^(5)_L <= min_t ((2^t + (2/5)^t)/2)^(L-1) = 0.9905^(L-1)`.
- **Reading.**
  - For `3x+1`, whose drift is contracting and which is AM-fair, one flip per bad orbit suffices; THM-4475's construction uses exactly one. The single-hit relaxation is then informative, and all its losses are polynomial or `2^(o(L))`, conjecturally.
  - For `5x+1`, one flip cannot create descent. A positive-drift orbit needs its odd-step frequency pushed from `1/2` below `log_5 2 = 0.43`, i.e. `Θ(L)` flips (heuristic).
  - So the observed floor (price about `0.64 β_L`, THM-4475 §6) is a multi-flip phenomenon, invisible to Theorem CS and Theorem P alike.
  - The degradation of Cauchy–Schwarz (`g_5(2) > 1`) is a symptom, not the mechanism.

## 7. The fixed point: smoothing transform versus the exact 3-adic equation

**The exact 3-adic fixed point (PROVED).**
- The adjoint `L*` is the Markov operator of the forward walk. `x -> x/2` is a 3-adic isometry and `x -> (3x+1)/2` contracts by `1/3`, so the synchronous coupling gives `W_1(μL*, νL*) <= (2/3) W_1(μ, ν)` in the 3-adic metric.
- By Banach there is a unique fixed point `π`. Reversing time gives `π = law of Y = sum_t 3^(t-1) 2^-(G_1+...+G_t)` (up to conventions, the 3-adic Syracuse variable of Tao, arXiv:1909.03562; his definition UNVERIFIED).
- **The owner's "fixed point chain growth" is this equation:** `p = Lp`, i.e. `p(v) = p(2v)/2 + (3/2)[v≡2] p(Ev)`, for `p = dπ/dμ`, and `Z_k = L^k 1 -> π` weakly.
- Its class averages are `(0, 1, 2)` on the classes `0, 1, 2 mod 3`, the stationary law of the class chain. These are the coefficients `0, R, 2R` of Prop. 10 divided by `R`.
- In the carry-free averaged model the tree density of a root `a` is `R p(a)/a` (heuristic; §3.2 of the inverse-tree note is its class-level form). So `p` is the full 3-adic profile of the tree's large-scale shape.

**`L^2` structure (PROVED identities, FINITE-EXACT values).**
- Write `Y = 2^-G (1 + 3Ȳ)`, `Λ = ℓ(Y) := log_4(1 + 3Y)` (an isometry of `Z_3`, with inverse `t -> (4^t - 1)/3 = S^t(0)`), and `κ_± = sum_(g even/odd) 2^-g δ_(-g/2)`.
- Then `ν = law(Λ)` satisfies `ν|_(ε+3Z_3) = (Ψ_ε)_*(κ_ε * ν)`, with `Ψ_ε` a 3-adic similarity of ratio `1/3`. Hence:
  - `||π||_n^2 = 3 (||κ_+ * ν||_(n-1)^2 + ||κ_- * ν||_(n-1)^2) = sum_(lev ξ <= n-1) P_(1/4)(ξ) |ν̂(ξ)|^2`;
  - `P_(1/4)(ξ) = 3(|κ̂_+|^2 + |κ̂_-|^2) = 15/(17 - 8 cos 2π{ξ})`. It averages to `(1 + 4^-N)/(1 - 4^-N)` over any full level: criticality in Fourier form. Checked to `10^-15` for `n <= 11`;
  - equivalently `||π||_(n+1)^2 = sum_j 4^-|j| C_j^(n)` with `C_j^(n) = 3^n sum_y π_n(y) π_n(S^j y)`, the density of `Λ - Λ'` at `j`. Checked to `10^-14`;
  - `|κ̂_+|^2 : |κ̂_-|^2 = 1 : 4` pointwise, so the class `L^2` masses are exactly `1/5 : 4/5` at every level.
- Level energies: `E_1 = 2/3`, `E_2 = 10/21`, then `0.4616, 0.4642, ..., E_12 = 0.4689, E_15 = 0.4708`. The increments shrink (`6.4·10^-4`, then `5.7·10^-4`), so `||π||_n^2 ≈ 0.47 n`.
- So (O1) at the fixed point reads: the off-diagonal ladder correlations `C_j` stay bounded. At `n = 14`:
  - `C_1..C_4 = 0.728, 0.473, 1.237, 0.860`;
  - `C_9 = 2.11`, `C_27 = 2.36` (growth only with `v_3(j)`: a logarithmic singularity of the law of `Λ - Λ'` at `0`);
  - `C_1` still creeps up by about `0.0012` per level.
- In Fourier terms (O1) is an `ℓ^2` equidistribution statement: the level-`m` Fourier mass of `ν` does not pile up at archimedean-small frequencies `{ξ}`, where `P > 1`.
- Tao's method (per his abstract) estimates the characteristic function of a skew random walk on a 3-adic cyclic group at high frequencies. (O1) needs square-root cancellation on average over each level. We did not check how the two statements relate.
- `max_y π(y + 3^n Z_3) = π(-1 + 3^n Z_3) ≈ 1.462 · 2^-n`. The `L^∞` dimension is `log_3 2 = 0.631 > 1/2`, so no single point carries a divergent local `L^2` mass. The divergence of `||π||_n^2` is a bulk (multifractal) effect of the `s = 2` tilt, with odd frequency `3/4`.

**Tails.**
- 3-adic: `x^2 P(p_15 > x) = 0.85, 0.92, 0.94, 0.99, 0.89, 0.96` at `x = 6, 8, 12, 16, 24, 32`, a plateau that the finite resolution cuts off beyond. The finite layers show the same: `x^2 P(Z_14 > x) ≈ 0.6–0.9` on `[4, 64]`.
- Mean field: `x^2 P(W > x) ≈ 1.6–1.7` on `[32, 128]`.
- Both are `~ c x^-2` with different constants, since the 3-adic second moment grows with slope 0.357 against 0.5.
- **The mean-field theorem.** For `W = W'/2 + (3/2) B W''` (i.i.d., `E W = 1`), `P(W > x) ~ c x^-κ` with `κ = 2`. This is the Kesten–Goldie / Guivarc'h / Liu theory (UNVERIFIED as to exact hypotheses). The abstract of Jelenković–Olvera-Cravioto, arXiv:1012.2165, describes an implicit renewal theorem for power tails of `R = sum C_i R_i + Q` with general weights. Non-lattice holds because `log 3/log 2` is irrational.
- For the 3-adic `p` the `x^-2` tail is (O3).

**Contraction (mean field; PROVED from the two ideal properties of `ζ_s`, each a one-line consequence of its definition as a supremum over functions with `(s-1)`-Hölder derivative).**
- For laws with equal means and `1 < s <= 2`: `ζ_s(Sμ, Sν) <= E[sum_i A_i^s] ζ_s(μ, ν) = g(s) ζ_s(μ, ν)`. This uses `ζ_s(cX, cY) = c^s ζ_s(X, Y)`, `ζ_s(X+Z, Y+Z) <= ζ_s(X, Y)` for independent `Z`, and one swap per summand.
- `g(s) = 0.987, 0.974, 0.966, 0.975, 0.988, 1.000` at `s = 1.1, 1.25, 1.5, 1.75, 1.9, 2`.
- So **the Banach argument works exactly on `1 < s < 2` and fails at `s = 2`, the same place where Cauchy–Schwarz loses the square**.
- The mean-field fixed point has moments of every order `< 2` (Biggins/Liu theory, UNVERIFIED as cited). It provably has no second moment: if `E W^2 < ∞`, the fixed-point equation would give `E W^2 = g(2) E W^2 + E[sum_(i≠j) A_i A_j] = E W^2 + 1/2`.
- For the 3-adic density the same identity (Theorem T) says that `p ∈ L^2` would force `<p∘τ, p> = 0`. Numerically `||π||_n^2` diverges linearly (§7).
- In the 3-adic model the analogue of the martingale-limit question is `L^1` convergence of `Z_k`. The increments `||Z_(k+1) - Z_k||_1` fall slowly (ratios `0.905 -> 0.949`), as in the mean field (`E|Z_(k+1) - Z_k| = 0.225, 0.124, 0.076, 0.053` at `k = 10, 20, 30, 39`). This is consistent with (O3), not a proof.

**Mean-1 martingale and drift (PROVED, elementary).**
- `E Z_k = 1` exactly for all `k`: this is `s = 1`.
- `-g'(1) = log(2/√3) = 0.143841` is the drift. `R = 1/(-g'(1))` is the tree-density constant of Prop. 10.
- `g'(2) = ln(27/16)/4 = 0.1308 > 0` is the Kesten–Goldie slope at `κ = 2`.

## 8. What would close the gaps (for the coordinator to number; no files created)

- **(Q1) Polynomial second moment.** `sup_k γ_k < ∞`; equivalently in practice, the ladder correlations `<Z_a∘S^j, Z_b>` (`j >= 1`) are bounded uniformly, or `sup_n E_n < ∞` for the Syracuse law.
  - Consequence: `M2(L) = O(L^3)` and `δ_L >= 2^(-2ηL - O(log L))`, the coordinator's `0.1001`.
  - Evidence: `k <= 15`, `n <= 15`.
  - The exact reduction is §2 (iv) and §7.
  - Known weaker results: the certified `θ_15` (Theorem M).
- **(Q2) (SD) ⟹ HYP-9137.** Quenched hub weight along bad orbits `2^(o(L))`.
  - Evidence: `H_L M2/ρ_L ≈ 3.2` stable for `L <= 16`; true in the mean-field model.
  - This isolates the one missing ingredient: 3-adic sibling decorrelation conditional on 2-adic badness.
- **(Q3)** `π` is absolutely continuous and `P(dπ/dμ > x) ~ c x^-2`.

## 9. Citations

**Read this session.** arXiv abstracts, fetched through the arXiv API with a generic user agent on 2026-09-25:
- Tao, arXiv:1909.03562 (almost all Collatz orbits attain almost bounded values; characteristic function of a skew random walk on a 3-adic cyclic group at high frequencies).
- Jelenković and Olvera-Cravioto, arXiv:1012.2165 (implicit renewal theorem for trees with general weights).
- Alsmeyer, Biggins and Meiners, arXiv:0906.3133 (fixed points of the smoothing transform).
- Buraczewski, Damek and Zienkiewicz, arXiv:1504.03144 (positivity of the tail constant for `R = sum_(i<=N) A_i R_i + B` with i.i.d. `A_i` and `E|A|^α = 1/N`). Not directly applicable: our two weights are not identically distributed.

**UNVERIFIED** (cited from memory for orientation only; not read):
- Kahane–Peyrière (1976) and Biggins (1977): non-degeneracy of the mean-1 martingale when `g'(1) < 0`.
- Durrett–Liggett (1983): fixed points of the smoothing transformation.
- Guivarc'h (1990) and Liu (2000): `x^-κ` tails.
- Rösler (1992), Zolotarev (1976), Neininger–Rüschendorf (2004): contraction method and ideal metrics.
- von Bahr–Esseen (1965).
- Iglehart (1974) and Bolthausen (1976): random-walk meanders.
- The precise form of Tao's Fourier decay proposition.

**In-repo.** THM-4470, THM-4475, HYP-9137, and the inverse-tree note (Props. 9–11).

## 10. Reproduction

`python3 04-computation/experiments/procgen_cauchy_20260925_run.py > 05-knowledge/results/procgen_cauchy_20260925.out`

Five parts run one after another as separate processes, in about 60 s. Each stays below 500 MB footprint (the Syracuse part peaks near 490 MB) and 700 MB maximum RSS. Seeds are fixed; only timing and memory lines vary between runs.

| part | script | content |
|---|---|---|
| 1 | `procgen_cauchy_20260925_moments.py` | exact tables and Gram matrix (`k <= 15`), identity and ladder checks, `M2(L)` for `L <= 16`, Hölder, `δ*_L`, tails, SHEET, DRIFT |
| 2 | `procgen_cauchy_20260925_syracuse.py` | `π mod 3^n` (`n <= 15`), level energies, `1 : 4`, ladder `C_j`, Poisson identity by FFT, density tail |
| 3 | `procgen_cauchy_20260925_majorant.py` | certified `θ_r` (`r <= 15`), explicit rigorous bounds |
| 4 | `procgen_cauchy_20260925_harmonic.py` | `H_L` exactly (`L <= 13`) and by Monte Carlo (`L <= 16`), quenched vs annealed, hub bias |
| 5 | `procgen_cauchy_20260925_meanfield.py` | `g_q`, Zolotarev constants, mean-field population dynamics, `L^1` increments, strip confinement |
