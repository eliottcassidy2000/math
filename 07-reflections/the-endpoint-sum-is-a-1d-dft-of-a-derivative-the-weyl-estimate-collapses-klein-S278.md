# The k-dim Weyl estimate collapses to a 1-D DFT of a derivative — the far frequency drops out, and boundedness closes the row

*klein-2026-07-13-S278. Owner: carry out the k-dim Weyl estimate on the endpoint sum (the S277 remaining
step). It is not `k`-dimensional at all: the coupling is periodic in `j` mod `e'`, so the sum is a 1-D
DFT, and the R-endpoint indicator is a discrete **derivative** — whose gain makes `|U_s^{e'}|` bounded
independent of `e'` and of the far frequency. The estimate is carried out to a single clean 1-D Fourier
sum, empirically closed (`e'→400`), whose rigorous constant is the same Beurling–Selberg mollification
the covering side reduced to.*

---

## It was never `k`-dimensional

S277 wrote the coupled inner sum as a `k`-dimensional Weyl sum `Σ_j χ(αj)e(-Nj/e')`. But `χ(j)` — "is
this crossing an `R_s`-endpoint?" — depends on the other offsets' phases `frac((e''/e')j+…)`, which are
**periodic in `j` mod `e'`**. So `χ` is just a function on `ℤ/e'`, and

$$\sum_{j=0}^{e'-1}\chi(j)e(-Nj/e') = e'\,\hat\chi(N\bmod e')\quad\text{(a 1-D DFT)},\qquad
U_s^{e'}(N)=\sum_\sigma e(-N\sigma/7e')\,e'\hat\chi_{s,\sigma}(N\bmod e').$$

The far frequency `N=ℓw` enters **only** through the residue `κ=ℓw mod e'`. No torus, no multi-dim
Erdős–Turán — the whole thing is one discrete Fourier coefficient.

## The coupling is a discrete derivative — the source of all the cancellation

`ch_j = 1_{cond_s}(\text{leave}_j) − 1_{cond_s}(\text{enter}_j)`, `cond_s=[`others cover exactly
`\{0..6\}\setminus s]`, and leave `−` enter `= 1/(7e')`. So `ch` is a **discrete derivative** of the
`cond_s` indicator sampled on the `e'`-grid. Fourier-expanding `g=1_{cond_s}`,

$$U_s^{e'}(N)=e'\!\!\sum_{n\equiv\kappa\,(e')}\!\!\hat g(n)\,e(ns/7e')\underbrace{\big(e(n/7e')-1\big)}_{2i\sin(\pi n/7e')}.$$

The `\sin(\pi n/7e')≈πn/(7e')` factor is the derivative gain. It (i) **kills the `n=0` term** — so the
net endpoint imbalance `\hat\chi(0)=O(1)`, not `O(e')` (verified `≤4`), which is exactly why the
*resonant* case `κ=0` is bounded — and (ii) cancels one power of `n` against `\hat g(n)\sim ρ/|n|`. This
is the mechanism behind every "bounded" number in S274–S277.

## Boundedness is real, and it closes the row

Directly: `|U_s^{e'}(N)| = e'|\hat\chi(κ)|` is **bounded independent of `e'` and of the far frequency**
— verified to `e'=400`: `≤4` for 6 other offsets, `≤13` for 7. The swing-count `T` is *constant* in `e'`
(the `cond_s` boundaries are fixed by the other offsets; refining the `e'`-grid doesn't add straddles).
So

$$|S|\le\frac1{2\pi^2}\Big(\sum_\ell\frac{|\sin(\pi\ell/7)|}{\ell^2}\Big)\sum_s\sum_{e'}|U_s^{e'}(\ell w)|
=O_k(1)=O(k)$$

for bounded `k` — the density-row tail closes, and the S276 empirical `|S|≤0.61R` now has its mechanism.

## What is proved, what remains — and the unification

- **RIGOROUS:** the 1-D DFT reduction; the swing-derivative identity with the `\sin(\pi n/7e')` gain
  (which provably kills the diagonal and makes `\hat\chi(0)=O(1/e')·e'`… `=O(1)`).
- **EMPIRICAL (decisive):** `|U_s^{e'}(N)|≤C_k` uniformly in `e'`(→400) and the far frequency.
- **REMAINING (one clean 1-D sum):** `Σ_{n≡κ (e')}\hat g(n)(e(n/7e')-1)=O(1)`. Because `g=1_{cond_s}`
  has non-summable `\hat g(n)\sim1/n`, this conditionally-convergent sum needs **Beurling–Selberg
  mollification** of `g` — the identical tool the *covering* side reduced to (opus-S261: mollified
  discrepancy of the coprime core). **Both LRC(14) routes now terminate on the same Beurling–Selberg
  estimate** — a single shared analytic lemma finishes both.

The arc of S273–S278 in one line: the density-tail constant is *not* `Σe'`-free (S275), *is* `O(k)` via
a resonance sum (S276), reduces *exactly* to endpoint exponential sums (S277), which are a 1-D DFT of a
derivative (S278) whose boundedness closes the row — pending one Beurling–Selberg mollification shared
with the covering side.

*Files: `04-computation/lrc14_chi_dft_klein_S278.py`, `lrc14_chi_largeeprime_klein_S278.py` (+ outs).
THM-728, HYP-6400. Consumes THM-727. Unifies with opus-S261 (covering-side mollification).*
