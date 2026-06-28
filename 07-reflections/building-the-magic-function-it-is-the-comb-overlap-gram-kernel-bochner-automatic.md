# Building the cyclotomic magic function: it is the comb-overlap GRAM KERNEL (Bochner-automatic), not a spatial minorant

*mac-mini-2026-06-28-S75. The owner asked to build out the cyclotomic Delsarte/Beurling-Selberg magic function
(the analytic resolution S74/kps-S31an converged on). I built it — and the key lesson is WHERE it lives: not a
spatial trig-minorant (slow, not Bochner), but the **comb-overlap Gram kernel** in the resonance domain, where
positivity is automatic and the cap is exact. Builds on THM-576 (cap = pairwise avoidance), kps HYP-3212/3213
(Delsarte/Chebyshev), [[the-one-obstruction-worst-case-algebra-vs-analytic-equidistribution-apex7-is-the-lee-yang-zero]].*

## The wrong build (honest negative): the spatial minorant
First I built the obvious thing: the optimal trig-polynomial minorant `F(x) ≤ 1_lonely(x)`, maximize `∫F`
(`lrc_cyclotomic_magic_function_macmini_S75.py`, LP). It **converges slowly** (`∫F` for `{1,13}`: gap 0.13
even at degree 56) and the optimal `F` is **NOT Bochner** (`f̂` has many negative coefficients). So the magic
function does **not** live in the spatial domain — minorizing the complicated lonely set blindly is the wrong
move. (This is the spatial shadow of S74: a config-blind certificate is weak.)

## The right build: the magic function IS the comb-overlap Gram kernel
The cap is the inclusion-exclusion (THM-576): `cap_P = Σ_{S⊆P} (−1)^|S| meas(∩_{p∈S} D_p)`,
`D_p = {x : ‖px‖ < 1/14}` (a width-`1/7` comb of `p` teeth). The **order-2 core** is the kernel
```
   K(p,q) = meas(D_p ∩ D_q) = ⟨ 1_{D_p}, 1_{D_q} ⟩   (inner product of comb indicators)
```
- **It is BOCHNER-positive AUTOMATICALLY** — `K` is a Gram matrix of the indicator vectors `1_{D_p}`, so it is
  PSD by construction. (Verified: kernel-matrix eigenvalues on `{1,5,7,8,9}` all `> 0`.) This is exactly the
  `f̂ ≥ 0` "magic" condition kps named — and here it costs nothing, because the kernel is a Gram matrix. **The
  positivity is structural, not imposed.**
- **Clean closed forms** (`lrc_magic_function_resonance_domain_macmini_S75.py`):
  - `K(1,q) = 1/(7q)` exactly — minimized at the antipode `q=13 ≡ −1 (mod 14)` giving `1/91`. The cap minimizer
    `{1}∪top` is the **least-resonant** comb pair. This is the cyclotomic core.
  - `K(7,q) = 1/49 = 1/7²` for **every** `q` — the apex prime 7 saturates (its comb is the full `1/7`-grid).
  - The general `K(p,q)` is NOT `gcd²/(7pq)` (that guess FAILED — tooth alignment matters); only the
    `p=1` and `p=7` rows are clean. Honest: the full kernel closed form is open.

## The cap reconstruction (the magic function works, exact)
With `K` the order-2 magic function and higher orders the corrections (VERIFIED):
```
  j   P            order-2          full = cap        correction beyond order-2
  2   {1,13}       66/91 EXACT      66/91             0          (pairwise magic fn is EXACT)
  3   {1,12,13}    0.6154           60440/100000      needs order-3
  4   {1,11,12,13} 0.5281           1979/4004=cap_9   −37/1092   (clean rational)
  5   {1,5,7,8,9}  0.4852           2243/5880=cap_8   −61/588    (clean rational)
```
So: **the order-2 Gram-kernel magic function is EXACT for `j=2`** and the dominant term throughout; the binding
rows `j=4,5` (= `k=9,8`) need an **order-3 correction that is a clean rational** (`−37/1092`, `−61/588`). The
"entire analytic difficulty" (THM-576's two constants) = these order-3 triple-overlaps `meas(D_p∩D_q∩D_r)`.

## The single-arc peeling LEMMA (PROVED) — the magic function closes on speed 1
Pushing the kernel to all orders gave a clean, **elementary-proved** lemma:
> **LEMMA (single-arc dominance).** For any `S` with `1 ∈ S` and all speeds `≤ 13`:
> `meas(∩_{p∈S} D_p) = 1/(7·max(S))`.
> **Proof.** `D_1 = (−1/14, 1/14)` is a single arc. For `p ≤ 13`, `1/p > 1/14`, so the only tooth of `D_p`
> meeting `D_1` is the central one `(−1/(14p), 1/(14p))`. The intersection over `S` is the smallest central
> tooth `(−1/(14 max S), 1/(14 max S))`, of measure `1/(7 max S)`. ∎
Verified: the formula holds for **all** `S ∋ 1` with speeds `≤ 13`; the **only** failures are at speeds
`≥ 14` (the apex modulus `14`, where `1/p ≤ 1/14` lets a second tooth enter `D_1`). The apex `14` is exactly
the boundary of the lemma.

Summing the inclusion-exclusion over subsets containing `1` (the `(−1)^{|S|}/(7 max S)` telescopes, leaving
only the `min`) gives a **PROVED peeling recursion**:
> **`cap(P) = cap(P∖{1}) − (1/7)·(1 − 1/min(P∖{1}))`**  for `1 ∈ P`, speeds `≤ 13`.
Checks: `cap({1,13}) = 6/7 − (1/7)(12/13) = 66/91` ✓; `cap({1,5,7,8,9}) = 2243/5880 = cap_8` ✓ (the binding
row, reconstructed exactly). So the magic function **closes on the speed-1 coordinate in closed form** — the
hardest-looking row `k=8` included. This generalizes THM-576's elementary `j=2` proof to peel speed 1 from
**any** `P`. The remaining work is the same peeling for `P∖{1}` (no single arc — `D_q` has `q` teeth), i.e. a
closed form for the top-speed sub-cap; that is the residual order structure.

## What this build delivers, and what remains
- **BUILT:** the magic function = the comb-overlap **Gram kernel** `K(p,q)=⟨1_{D_p},1_{D_q}⟩`, **Bochner-positive
  by construction**, reproducing the cap via inclusion-exclusion; `K(1,q)=1/(7q)`, `K(7,q)=1/49`. The Delsarte
  positivity (`f̂≥0`), the Toeplitz PSD margin, and the Perron mode are all just the **Gram/PSD structure of
  this kernel** — kps's "4 faces of one magic function" = the one Gram kernel.
- **REMAINING (the order-3 magic function):** close the triple-overlap term `meas(D_p∩D_q∩D_r)` in closed form
  (the binding constants `−37/1092`, `−61/588`); this is the order-3 / associator / non-SOS residue of S74. The
  even/order-2 (Gram, Bochner) half is built and exact; the odd/order-3 half is the remaining construction.
- **HONEST:** this builds the certificate and pins its positivity to the Gram structure; it is NOT a proof of
  LRC(14). The general kernel closed form (beyond `p∈{1,7}`) and the order-3 term are open. But the magic
  function now has a concrete, correct form — the comb-overlap Gram kernel — and the remaining work is the
  order-3 triple-overlap, two explicit rationals.

## Convergence with kps S31ao (HYP-3214): the Fourier side = the Fejér kernel
kps independently built the magic function on the FOURIER side and identified it exactly:
**it is the Fejér kernel `F_7`** — `(de Moivre cubic)²(2cos t) = (V_7(u)−2)/(u−2) = sin²(7t/2)/sin²(t/2) =
F_7(t)`, with `F̂_7(n) = (7−|n|)_+` (Fejér weights ≥ 0 = Delsarte/Bochner), `F_7(0) = 49 = 7²`, and **double
zeros at the de Moivre angles** (LP-sharpness, Cohn-Elkies/Viazovska). This is the **Fourier dual of my spatial
comb-overlap Gram kernel**: `K(p,q) = ⟨1_{D_p},1_{D_q}⟩` is the autocorrelation whose Fourier transform is the
Fejér weights; my `K(7,q) = 1/49` is the spatial echo of kps's `F_7(0) = 49`. The two builds are the same object
from two sides:
- **kps (Fourier):** `F_7`, `F̂≥0`, double zeros at de Moivre, the modular `Γ_0(7)` Eisenstein version, `F_k`
  = AP autocorrelation.
- **me (spatial/measure):** the comb-overlap Gram kernel, the PROVED single-arc lemma `meas(∩_S)=1/(7 max S)`
  for `S∋1` (speeds ≤13), and the **peeling recursion** `cap(P)=cap(P∖{1})−(1/7)(1−1/min(P∖{1}))`.
My Verblunsky/OPUC angle (S73d) = kps's "OPUC Christoffel-Darboux" facet. The peeling recursion + single-arc
lemma are the **proved spatial content** that complements kps's Fejér-kernel Fourier construction.

Related: HYP-3227 (this), HYP-3221 (one obstruction), THM-576 (cap=pairwise), kps HYP-3214 (Fejér kernel),
HYP-3212/3213 (Delsarte/Chebyshev), the scripts above, OPEN-Q-108.
