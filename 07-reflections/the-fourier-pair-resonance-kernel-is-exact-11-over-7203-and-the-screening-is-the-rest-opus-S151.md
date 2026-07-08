---
source: opus-2026-07-08-S151
status: the Fourier pair-resonance kernel derived EXACTLY (c_pair = 11/7203 = theta^3(2/3-theta),
  via Bernoulli B4), giving THM-641's pair overlap mass a closed value and the R2-linear leading
  term Var_resonance = (R2/2)c_pair; the screening (kps-S81's ~96% cancellation) quantified
  exactly and shown to be the k-dependent multiple-overlap content. Does NOT close Var<=c*R2 with
  the tight constant -- the triple/quad masses (the screening) remain, as kps-S81 flagged.
tags:
  - lrc14
  - covering-floor
  - variance
  - fourier
  - additive-energy
  - resonance-lemma
  - bernoulli
---

# The Fourier pair-resonance kernel is exactly 11/7203, and the screening is the rest

**opus-2026-07-08-S151 (HYP-5407).** Owner: prove the resonance lemma `Var(W) <= c*R2` via the
Fourier pair-resonance kernel. The pair kernel comes out in exact closed form; the tight
constant does not — it is screened by the multiple-overlap (triple/quad) terms, exactly the
target kps-S81 and mac-mini isolated concurrently. This note pins both.

## The pair-resonance kernel, exactly

Arc indicator `a_i(x,y) = chi(e_i x - y)`, `chi = 1_{(-theta,0]}`, `theta = 1/7`; Fourier
`hhat(n) = (1 - e(n theta))/(2 pi i n)`, `|hhat(n)|^2 = sin^2(pi n theta)/(pi^2 n^2)`. The
overlap of two arcs at phase-offset `delta` is the tent `t(delta) = (theta - ||delta||)_+`, the
autocorrelation of the arc, so its Fourier coefficient is `|hhat(n)|^2` and its variance is
`sum_{n!=0}|hhat(n)|^4`. Matched-difference pairs (`e_i - e_j = e_k - e_l`, the additive energy)
carry the *same* frequency `e(n(e_i-e_j)x)`, so they add coherently and the `R2`-linear leading
term of `Var(W)` is

> **`Var_resonance(W) = (R2/2) * c_pair`,  `c_pair = sum_{n!=0}|hhat(n)|^4`.**

The kernel has an **exact closed form** (the point of this session):

> **`c_pair = sum_{n!=0}|hhat(n)|^4 = theta^3 (2/3 - theta) = 11/7203`**  (`theta = 1/7`).

Three independent confirmations: (i) the direct `sum sin^4(pi n/7)/(pi^4 n^4)` numeric; (ii) the
`k=2` tent variance `2 theta^3/3 - theta^4` (a two-arc family has exactly one difference, so
`Var(W) = c_pair` at `k=2` — verified equal); (iii) the Bernoulli evaluation
`sum cos(2 pi n a)/n^4 = -(pi^4/3) B_4(a)` with `sin^4 = (3 - 4 cos 2u + cos 4u)/8`:
`c_pair = (1/4)[1/30 + (4/3)B_4(1/7) - (1/3)B_4(2/7)]`, `B_4(x)=x^4-2x^3+x^2-1/30`. This gives
klein's THM-641 pair overlap mass its exact value, in Fourier form.

## The screening: why the tight constant is not the pair kernel

`Var_resonance = (R2/2)c_pair` is the leading term, exact at `k=2` (`s=1`), but at `k=11` the
true `Var(W)` is `~8%` of it — the multiple-overlap **screening** (kps-S81's `~96% cancellation`;
mac-mini's `pair-part 27x`). Exact screening factor `s(k) = Var(block_k)/Var_resonance(block_k)`:

| k | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|----|----|----|----|
| `s(k)` | 1.00 | 0.80 | 0.64 | 0.52 | 0.39 | 0.31 | 0.22 | 0.16 | 0.11 | **0.080** | 0.058 | 0.043 |
| `Var(W)` | .0015 | .0061 | .0136 | .0236 | .0328 | .0424 | .0464 | **.0483** | .0478 | .0470 | .0448 | .0429 |

Two facts kill any hope of a universal `Var(W) = c*R2`:
- **`s(k)` falls with coverage `k*theta`** — each added arc screens the pair variance by roughly
  `5/7` (`= 1 - 2 theta`, the disjoint fraction) at large `k`. The tight `k=11` constant
  `c = s(11) c_pair/2 = 6.1e-5` is the *screened* value, not the pair kernel `c_pair/2 = 7.6e-4`.
- **`Var(W)` is NON-monotone in `k`** — it peaks at `k~9` (`k*theta ~ 1.3`, arcs cover once over)
  and falls thereafter (as `k -> inf` the circle is fully covered, `W -> 0`, `Var -> 0`), while
  `R2` grows monotonically as `(k-1)k(2k-1)/3`. So `Var ~ c*R2` is at best a **local** relation in
  the `k=11` window, and the constant is genuinely `k`-dependent.

## What this means for the resonance lemma (honest)

- **Delivered:** the pair-resonance kernel `c_pair = 11/7203 = theta^3(2/3-theta)`, exact
  (Fourier/Bernoulli); the `R2`-linear leading term `Var_resonance = (R2/2)c_pair`; the exact
  screening `s(k)` (the `~96%` cancellation at `k=11`) and the non-monotonicity of `Var(W)`.
- **Not delivered:** the tight `Var(W) <= c*R2` with `c <= 7e-5`. The pair kernel gives only the
  leading term / a `~12x`-loose bound; the tight constant is the *screened* value, and the
  screening = the triple- and quad-arc overlap masses — precisely the concrete target kps-S81
  named (the `THM-641`-analog for `|S|=3,4`, derivable by the same Bernoulli-periodization the
  pair kernel used here). A rigorous `Var <= (R2/2)c_pair` (the loose pair bound) would need the
  net higher-order terms to be `<= 0`; the sign is mixed (`+` triple-triple, `-` pair-triple), so
  even the loose bound is not immediate.
- **The clean route forward** (Fourier): the triple mass is the Fourier coefficient of the
  three-arc overlap, a product of three `hhat`'s constrained to frequency-sum zero; its
  `R2`-and-higher-energy-linear leading term has the same closed-form structure as `c_pair`
  (a `B_6`/`B_4` Bernoulli combination). Summing the pair kernel `+` the triple/quad kernels is
  the Fourier form of the screening — the same object as kps's real-space overlap-mass laws,
  and the place the tight `c` lives.

## Ledger

- EXACT/NEW: `c_pair = 11/7203 = theta^3(2/3-theta) = sum_{n!=0}|hhat(n)|^4` (Fourier/Bernoulli,
  3-way verified); `Var_resonance = (R2/2)c_pair`; the exact screening `s(k)` (1.00 -> 0.043,
  k=2..13); `Var(W)` non-monotone (peaks k~9).
- REDUCTION: the tight `c` = the screened value = pair kernel times `s(k)`; `s(k)` = the
  triple/quad overlap masses (kps-S81's target), derivable by the same Bernoulli method.
- HONEST: does NOT close `Var <= c*R2` with the tight constant; pins the pair kernel exactly and
  quantifies the screening precisely.
- Files: `lrc14_fourier_pair_kernel_opus_S151.py` (+out).
- Builds on / cites: kps-S81 (the exact `W = sum_S (-1)^|S| L_S` decomposition, the `96%`
  cancellation, the triple/quad target), mac-mini LEM-006 (`far<=E[W]^2 <=> Var<=near`,
  pair-part `27x`), klein THM-641 (pair overlap mass) + THM-656 (`Var(F)=R2*V1`), opus-S150
  (the unification). External: Bernoulli polynomials / Clausen `sum cos(2 pi n a)/n^4`.
