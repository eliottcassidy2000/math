# The cyclotomic ideal: Joukowski, de Moivre, and the dip as a rational-approximation defect

*kind-pasteur-2026-06-27-S31ak. The owner's connection: the ideal uniform PGF `1+z+···+z⁶` (perfect
7-fold symmetry) has roots at the 7th roots of unity, which Joukowski `w=z+1/z` sends to the **de Moivre
angles** `{2cos(2πj/7)} = {−1.8019, −0.4450, 1.2470}` (roots of the cyclotomic cubic `x³+x²−2x−1`). So
**cap = the 7th-cyclotomic ideal; dip = deviation from 7-fold symmetry = Im(w) = the real-rootedness
defect.** Plus: rational approximations of irrationals/transcendentals; merge the ferromagnetic picture
and creative niche sets. This linearizes my Lee-Yang circle (HYP-3099) and reframes the dip arithmetically.*

## The verified geometry (the lens is real)
- **Ideal PGF `1+z+···+z⁶` = (z⁷−1)/(z−1)`** — roots = the 6 non-trivial 7th roots of unity, on `|z|=1`,
  perfect 7-fold symmetry.
- **Joukowski `w=z+1/z`** sends `e^{2πij/7} ↦ 2cos(2πj/7)`, collapsing the 6 roots (conjugate pairs `j`,
  `7−j`) to the **3 de Moivre angles** `{−1.8019, −0.4450, 1.2470}` = the roots of `x³+x²−2x−1` (the real
  cyclotomic cubic, the maximal real subfield `ℚ(cos 2π/7)` of `ℚ(ζ₇)`). VERIFIED exactly.
- My Lee-Yang result: the **coverage PGF** `Q(z)=Σ q_t z^t` has 6 zeros in 3 conjugate pairs near a
  circle `|z|=R`. The **R-normalized Joukowski** `w_j = z_j/R + R/z_j` sends them to 3 points near the
  real axis — the **de Moivre lattice**. So the owner's map *is* the linearization of my circle.

## What the computation shows (`lrc_joukowski_cyclotomic_dip_kps.py`)
- **The coverage approaches the cyclotomic ideal as `k→∞`:** the Joukowski `Re(w)` distance to the de
  Moivre angles **decreases** `0.380 → 0.212 → 0.169` for `k=8,9,10`. The **dominant mode matches almost
  exactly** already at `k=8`: `Re(w)₁ = −1.804` vs de Moivre `−1.8019` (the most-negative angle, the
  "ground" mode). So the deepest cyclotomic mode is locked; the deviation lives in the sub-dominant modes.
- **The real-rootedness defect `Σ|Im(w)|` is genuinely present** (consec_8: `0.174`) and is the
  Joukowski avatar of my off-circle `λ`. **Honest caveat:** it is *not* numerically equal to the dip
  (`dip_8=0.0141`), and consec does *not* minimize the raw defect (rank 38/401) — a low-coverage config
  can sit exactly on a circle (defect 0). So "**dip = Im(w)**" holds as a *qualitative identification*
  (deviation-from-7-fold = off-circle = real-rootedness defect), not a clean equality. The precise
  statement is the next item.

## The dip as a RATIONAL-APPROXIMATION defect (the owner's lead, made precise)
The de Moivre angles are **degree-3 irrationals** (`ℚ(cos 2π/7)`, the heptagon is non-constructible). The
**cap `C(k+1,2)/91` is RATIONAL** (THM-577). So:
> **The cap is a RATIONAL APPROXIMATION of the cyclotomic (irrational) ideal, and the dip is the
> approximation defect.** The LRC is, at bottom, simultaneous **Diophantine approximation**
> (`M(S)=max_t min_i ‖s_i t‖` is "how badly can all `s_i t` be approximated by integers"), and the
> three-gap structure (THM-565) is the continued-fraction skeleton. The dip `= cap − (cyclotomic ideal)`
> is the **error of the best rational (binomial/Pascal) approximant to the 7-fold-symmetric algebraic
> value** — and it is rational *because* the De Moivre solvability of the `k=8` resolvent collapses to `ℚ`
> (perfect-square discriminant 9, HYP-3132). The two facts fuse: **the resolvent is rational (the dip is a
> rational number) precisely because the 7-fold cyclotomic obstruction is "tamed" at the apex** — the
> heptagon's degree-3 irrationality is folded away by the `s↦6−s` reflection into the solvable biquadratic.
> The *residual* irrationality (the genuine real-rootedness defect, `Im(w)` in the sub-dominant modes) is
> the **odd/Worpitzky/non-associative half** (S31ai) that does *not* rationalize.

## Merge with the FERROMAGNETIC picture (HYP-3161)
The de Moivre ground angle `−1.8019` is the **most-ordered cyclotomic mode**, and the dominant coverage
zero locks onto it — the **ferromagnetic ground state's dominant correlation IS the cyclotomic ground
mode.** The ferromagnetic transition (`Σκ₂` sign-flip at `k=5→6`) is the onset of 7-fold cyclotomic order:
below it the orbit is too sparse to feel the 7 sectors (antiferro/disordered), above it the cyclotomic
ideal organizes the correlations (ferro/ordered). So: **antiferro (k≤5) = pre-cyclotomic; ferro (k≥6) =
cyclotomic-ordered; the de Moivre angles = the ordered-phase normal modes; the dip = the residual
real-rootedness defect of the ordered phase (the φ⁴ vertex / the heptagon's irreducible irrationality).**

## Creative niche sets (exploration discretion)
- **Beraha number `B₇ = 2+2cos(2π/7) = 3.2470`** — the `n=7` member of the Beraha sequence (chromatic-root
  / Tutte-polynomial accumulation points, `B₅=φ+1`). The coverage zeros' radius/angle may track `B₇`; the
  Beraha numbers are exactly the `2+2cos(2π/n)` cyclotomic family. LEAD: is the cap a Tutte/chromatic
  value of a 7-structured graph (the Fano/Heawood)?
- **The heptagon & the Lagrange/Markoff spectrum.** `2cos(2π/7)`'s continued fraction is eventually
  periodic (it's a cubic irrational — **not** a quadratic, so NOT periodic by Lagrange; its CF is a genuine
  open problem, like all cubics). The LRC's "badly approximable" tight configs sit in the **Lagrange
  spectrum**; the apex-7 worst case = a cubic (heptagonal) obstruction, beyond the quadratic
  (golden-ratio/`B₅`) cases that are classically tractable. **This is why 14=2·7 is the first open case:
  the cubic (degree-3, non-periodic-CF) cyclotomic obstruction first appears at 7.**
- **`PSL(2,7)` / the Fano plane / the Klein quartic** (order 168, the 7-point geometry) — the symmetry
  group of the 7-fold ideal; the de Moivre cubic is its `ℚ`-character field.
- **Salem/Pisot connection:** the resolvent's root tower `{2,4,8,16}` (HYP-3132) is dyadic (Pisot `2`);
  the cyclotomic `2cos(2π/7)` is the *unit* side. The dip lives between the Pisot (dyadic, even,
  rational-collapsing) and the Salem/cyclotomic (heptagonal, odd, irrational) towers.

## Net
- **VERIFIED:** the coverage PGF's Joukowski image approaches the de Moivre angles (cyclotomic ideal) as
  `k→∞`, dominant mode locked at `k=8`; the lens is structurally correct.
- **REFRAMED:** dip = the rational-approximation defect of the binomial/Pascal cap against the 7-fold
  cyclotomic (cubic-irrational) ideal; rational *because* the apex resolvent De-Moivre-collapses to `ℚ`,
  with the residual `Im(w)` defect = the odd/non-associative half.
- **MERGED:** ferromagnetic order (k≥6) = cyclotomic order; de Moivre angles = ordered normal modes;
  transition = cyclotomic onset. **`14=2·7` is first-open because 7 is the first apex whose cyclotomic
  ideal is a CUBIC irrational (non-periodic CF, heptagon non-constructible)** — the Diophantine wall.
- **HONEST:** "dip = Im(w)" is qualitative not exact; consec doesn't minimize the raw defect; the niche
  sets (Beraha/Markoff/PSL(2,7)/Salem) are flagged leads, not theorems.

→ HYP-3162 (this), HYP-3099 (Lee-Yang circle), HYP-3132 (biquadratic resolvent, ℚ-collapse), HYP-3161
(ferromagnetic transition), THM-577 (rational cap), THM-565 (three-gap = CF skeleton), the heptagon/
de Moivre cubic, Beraha numbers, the Lagrange spectrum, [[lrc14-thread]].
