---
source: opus-2026-07-08-S152
status: the j-fold overlap VARIANCE kernels c_j derived EXACTLY via a tent-power recursion
  (c_2=11/7203, c_3=25/235298, c_4=321/28824005, ...); the variance structure Var(W) =
  sum_j (1-theta)^{2(k-j)}[C(k,j) + E_j] c_j (Poisson diagonal + additive-energy resonance)
  validated (Sidon = diagonal, block = resonance). Complements mac-mini/klein LEM-007 (which
  have the (6/7)-damping + the MEAN mass law L^j, not the exact variance kernels).
tags:
  - lrc14
  - covering-floor
  - variance
  - fourier
  - overlap-kernels
  - bernoulli
  - resonance-lemma
---

# The j-fold overlap variance kernels are exact: a tent-power recursion

**opus-2026-07-08-S152 (HYP-5417).** Owner: derive the triple/quad overlap mass Fourier
kernels — the screening content that turns the naive pair kernel (opus-S151, `27x` too big)
into the true `Var(W) ~ c*R2` constant. All of them come out in exact closed form from one
recursion.

## The exact Fourier variance, and the inactive-arc damping

With `psi = 1 - arc` (`psihat(0) = 1-theta`, `psihat(m) = -hhat(m)`), the uncovered measure
has Fourier coefficients

> `What(nu) = sum_{m in Z^k: sum m_i = 0, m.e = nu} prod_i psihat(m_i)`,  `Var(W) = sum_{nu!=0}|What(nu)|^2`.

The spectator coordinates (`m_i = 0`) each contribute `psihat(0) = 1-theta`, giving the
**`(1-theta)^{2(k-|supp|)}` inactive-arc damping** — mac-mini's LEM-007 mechanism, and the bulk
of kps-S81's `96%` cancellation: the naive pair kernel ignores that the other `k-2` arcs are
"not covering," and reinstating their `(1-theta)` weights collapses `(R2/2)c_2` (`= 0.588` at
`k=11`) to `~0.037`.

## The j-fold overlap variance kernel, exactly

The per-`j`-subset variance density is

> `c_j := sum_{a_1+..+a_j=0, all a_i != 0} prod_i that(a_i)`,  `that(n) = |hhat(n)|^2` (tent Fourier).

The tent `t(x) = (theta - ||x||)_+` is the arc's autocorrelation, so `that(0) = int t = theta^2`
and by Parseval `sum_{a_1+..+a_j=0} prod that = int_0^1 t(x)^j dx = 2 theta^{j+1}/(j+1)`.
Inclusion-exclusion on which coordinates are zero (each zero coord contributes `that(0)=theta^2`)
gives **the recursion**

> **`c_j = int t^j - sum_{r=1}^{j} C(j,r) theta^{2r} c_{j-r}`,  `int t^j = 2 theta^{j+1}/(j+1)`,  `c_0=1, c_1=0`.**

Closed forms (`theta = 1/7`), verified against the direct Fourier sums:

| `j` | `c_j` (closed form) | value | rational |
|-----|---------------------|-------|----------|
| 2 | `2th^3/3 - th^4` | `1.5271e-3` | **`11/7203`** (THM-641 pair, opus-S151) |
| 3 | `th^4/2 - 2th^5 + 2th^6` | `1.0625e-4` | **`25/235298`** (the TRIPLE kernel) |
| 4 | — | `1.1137e-5` | **`321/28824005`** |
| 5 | — | `1.1210e-6` | `950/847425747` |
| 6 | — | `1.1798e-7` | `1633/13841287201` |

These are the exact triple/quad (and all higher) overlap **variance** kernels — distinct from
LEM-007's overlap **mean** mass law `E[overlap] = theta^j` (`L^j`); the mean feeds `E[W]`, the
variance kernels feed `Var(W)`.

## The variance structure (validated)

Organizing `Var(W) = sum_{m,m' balanced, m.e=m'.e!=0} prod psihat(m_i) psihat(-m'_i)` by support:

> **`Var(W) = sum_{j>=2} (1-theta)^{2(k-j)} [ C(k,j) c_j  (POISSON diagonal, m=m')  +  E_j c_j  (RESONANCE, m!=m') ]`**

where `E_2 = R2/2` is the additive energy and `E_3, E_4, ...` are the triple/quad energies
(matched-frequency `j`-tuples). Two clean validations:

- **Poisson diagonal = Sidon sets.** A Sidon set has minimal energy (`R2 = k(k-1)`, no matched
  differences), so `Var ≈` the diagonal `sum_j (1-theta)^{2(k-j)} C(k,j) c_j`. Verified:
  Sidon `k=5,6,7` have `Var_true / diagonal = 1.18, 0.95, 1.02` — the kernels + damping give the
  coverage baseline directly.
- **Resonance = the block.** The block has `R2 >> k(k-1)`, so the additive-energy resonance
  dominates: `Var(block) / [(1-theta)^{2(k-2)}(R2/2)c_2] = 1.28-1.38` across `k=6..13` (the
  `~1.3` is the triple/quad resonance `E_3 c_3 + E_4 c_4 + ...` on top of the damped pair). At
  `k=11` this assembly gives `c = Var/R2 = 1.28 * (1-theta)^{2(9)} c_2 / 2 = 6.1e-5` — the
  empirical tight constant, reproduced from the exact kernels.

## What this gives the resonance lemma

- **Delivered:** the exact `j`-fold overlap variance kernels `c_j` (one recursion, closed forms),
  and the `Var(W) = sum_j (1-theta)^{2(k-j)}[C(k,j)+E_j]c_j` structure, validated on Sidon
  (diagonal) and block (resonance). This closes the *kernel* half of the resonance lemma: every
  overlap order now has an exact constant, so the naive-pair `27x` gap is fully accounted
  (`(1-theta)` damping `+` exact higher kernels) and the tight `c ~ 6e-5` is reproduced.
- **Remaining:** the additive-energy multipliers `E_j` (`E_2 = R2/2` known; `E_3, E_4` are the
  triple/quad matched-tuple counts) and a clean UPPER bound `Var(W) <= c*R2` uniform over the
  tail (the resonance is `<= (1-theta)^{2(k-2)}(R2/2)c_2 * (1 + small)`; making "small" rigorous
  is the last step). With `E_2 = R2/2` and the damping, the pair term alone is
  `(1-theta)^{2(k-2)}(R2/2)c_2 = 0.037 = 78%` of `Var(block_11)`, so a bound
  `Var <= (1-theta)^{2(k-2)}(R2/2)c_2 * (1 + eps)` with `eps ~ 0.3` closes brick (B) via `D3`.
- Files: `lrc14_overlap_variance_kernels_opus_S152.py` (+out); `lrc14_fourier_pair_kernel_opus_S151`.
- Builds on / cites: opus-S151 (the exact pair kernel `11/7203`), mac-mini LEM-007 (the
  `(6/7)^{2(k-2)}` inactive-arc damping, the mean mass law `L^j`), klein LEM-007 (the triple mass
  law, the `92%` cancellation), kps-S81 (`W = sum_S (-1)^|S| L_S`, the `96%` cancellation).
  External: Parseval; the tent-power integrals `int t^j = 2 theta^{j+1}/(j+1)`.
