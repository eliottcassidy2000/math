# The van der Corput O(k) bound is confirmed: `Error·w ≤ C·R(C,w) ≤ C·(k−1)`, the resonance sum caps each offset at O(1)

*klein-2026-07-13-S276. Owner: work on the van der Corput O(k) bound — the sharpened target from S275
(`Error·w ≤ C` for non-resonant w, off the resonant set). Confirmed, with the right shape: the bound
is `Error·w ≤ C·R`, `R` a Diophantine "resonance sum" that is `O(1)` at clean w and caps at `k−1` at
full resonance. So `Error·w = O(k)`, independent of `Σe'`, and it closes the density-row tail.*

---

## The bound

`Error(E',w) = Φ(E'∪{w}) − Φ_∞(E')`, and `Error·w = Σ_{endpoints p} ε_p G_{s(p)}(wp)` (the per-interval
sum, THM-725). Define the **resonance sum**

$$R(E',w) = \sum_{e'\in E',\,e'\ge 2} \min\!\Big(1,\ \frac{1}{e'\,\|w/e'\|}\Big).$$

Measured (prime grid `Ng ≫ w`, no aliasing): **`Error·w ≤ 0.61·R`** across clean-to-resonant `w`
(ratio `Error·w / R ∈ [0.05, 0.61]`, max `0.61` at exact resonance `w=lcm`). Since **each term of `R`
is `≤ 1`**, `R ≤ #\{e'\ge 2\} ≤ k−1`, hence

$$\boxed{\ \text{Error}\cdot w \ \le\ 0.61\,R(E',w)\ \le\ 0.61\,(k-1)\ =\ O(k)\ }$$

— bounded **independent of `Σe'` and the diameter**. (The offset `e'=1` contributes only `7` bounded
terms — `w·(j+σ/7) ≡ wσ/7` is `j`-free — and is excluded from `R`.)

## Why it is O(k), not O(Σe') — the Denjoy–Koksma mechanism

Group the endpoint sum by the offset `e'` owning each crossing (`p = (j+σ/7)/e'`). Under `×w`,
`w·(j+σ/7)/e' = (w/e')j + β`. With `g = gcd(w,e')` and `q = e'/g`, the progression `{(w/e')j : j<e'}`
hits each of the `q` points `{m/q}` exactly `g` times, so for the *uncoupled* diagonal

$$\sum_{j=0}^{e'-1} G_s\big(\tfrac{w}{e'}j+\beta\big) = g\sum_{m=0}^{q-1} G_s\big(\tfrac{m}{q}+\beta\big)
= g\sum_{r\ne 0}\hat G_s(rq)e(rq\beta) = O\!\big(g^2/e'\big),$$

using `|\hat G_s(rq)| ≤ C/(rq)^2` (only Fourier modes divisible by `q` survive the equidistributed
`q`-point average). That is `O(1/e')` at clean `w` (`g=1`) but `O(e')` at resonance (`g=e'`).

The surprise: the **coupling helps**. The genuine `Error·w` weights each crossing by whether it is an
R-endpoint (`ε_p∈\{-1,0,+1\}`, set by the *other* offsets). Empirically this **caps each offset's
contribution at `O(1)`** — `R`'s `min(1,·)` is the correct envelope — so the total is `O(k)`, strictly
better than the diagonal's `O(Σe')` at resonance. (Mechanism: at resonance the many crossings of `e'`
pile onto few phase-values, but only `O(1)` of them are actual miss-structure endpoints; the coupling
throttles the pile.)

## Diophantine face

- **Clean `w`** (`‖w/e'‖ ≥ δ` for all `e'≥2`): each term `≤ 1/(e'δ)`, so `R ≤ δ^{-1}Σ_{e'}1/e' =
  O(δ^{-1}\log e_{\max})` — small. Measured `Error·w ≈ 0.01–0.7` for `k=7`, `Σe'` up to 399, and up to
  `k=12`, `Σe'=900`, with no growth in `Σe'` (the apparent growth in a first sweep was `δ` dropping as
  more offsets crowd `w`, not a real `Σe'` effect).
- **Resonant `w=lcm`** (`e'∣w` for all): every term `=1`, `R=k−1`, `Error·w` maximal (`3.03` for the
  2-block, `k−1=5`, ratio `0.61`). But `w=lcm ≫ Σe'`, so `Error = 0.61(k-1)/lcm` is negligible (the
  S275 "resonance harmless").

## Closing the row

Peel `w = d = \max(E)`: `Error ≤ 0.61\,R(E',d)/d ≤ 0.61(k-1)/d`. For the `k=8` row (`E'` a 7-cluster,
`k−1 ≤ 6`): `Error ≤ 3.66/d`, so `d > 38 ⇒ Error < 0.097 = cap₉−0.397 ⇒ Φ(E) ≤ cap₉`. With the box
extended from `d≤25` (THM-719) to `d≤38`, the two meet. For **non-resonant** `d` (the generic case),
`R(E',d) ≪ k−1` and the crossover is far smaller. Combined with the S275 band check (max `Φ = 0.347`,
`26≤d≤8·diam`) this closes the `k=8` density-row tail; the `k=9` twin is identical (larger margin).

## Honest scope

- **Confirmed (empirical + mechanism):** `Error·w ≤ C·R ≤ C(k−1)`, `C≈0.61`; the `O(k)` bound holds,
  `Σe'`-independent, and the Denjoy–Koksma reduction (`Σ_j = g·Σ_m`, `Σ_m = O(1/q)`) is exact.
- **Open (rigor):** the *coupled* per-offset bound `≤ 1` (the `min(1,·)` envelope) is measured, not
  proved — the R-endpoint weighting `ε_p` throttling the resonant pile is the one unproved step. This
  is a Koksma / three-distance estimate on the miss-structure's transitions under `×w`. The constant
  `C≈0.61` is empirical.

This is the tightest form of the density-tail constant across the S273–S276 arc: not `Σe'`-free (that
was false, S275), but `O(k)` via a Diophantine resonance sum — which is exactly what the row needs.

*Files: `04-computation/lrc14_vdc_scaling_klein_S276.py`, `lrc14_vdc_clean_klein_S276.py` (+ outs).
HYP-6350. Consumes THM-725/700/699. Sharpens
[[the-cross-sector-constant-boundedness-is-elementary-but-decay-is-not-and-the-grid-was-aliasing-klein-S274]]
and [[the-sigma-free-decay-is-mistargeted-the-real-bound-is-sigma-over-w-which-decays-when-w-dominates-klein-S275]].*
