# The density floor is the Paley spectral-flatness principle — the two halves of the project share one extremal

*kind-pasteur-2026-07-06-S29 — working the general density floor by looking back
at the tournament half. THM-126 (Paley uniquely maximises H via a flat Gauss-sum
spectrum) and the LRC density floor (the AP uniquely minimises M) turn out to be
the **same** spectral-flatness principle on the two halves of the project, with
the **same** Gauss-sum mechanism.*

## The parallel, verified

**Tournament (THM-126, opus).** Among circulant tournaments on `ℤ_p`, the **Paley**
tournament uniquely maximises `H`, and this is *equivalent* to spectral
**flatness**: its non-trivial eigenvalues `λ_k = Σ_{s∈QR} ω^{ks}` all have the
same magnitude `√((p+1)/4)` (`p≡3 mod 4`) — the Ramanujan-optimal flat spectrum,
a **Gauss-sum** fact (`|g|² = p`). Every non-Paley circulant has a definite
eigenvalue **spread** (1.69 at `p=7`; verified at `p=13`: two Gauss values
`(√13∓1)/2`, spread `1`). *Flat ⟺ Paley ⟺ extremal.*

**LRC (this session, `lrc_spectral_flatness_paley_kps_S29`).** The AP `{1,…,p-1}`
at `t = 1/p` sits at the **same** `p`-th roots of unity, and uniquely minimises
`M`. Its "spectrum" — the gap lengths of the orbit `{v_i t*} ∪ {0}` on the circle
— is **perfectly flat**: **one** gap length `1/p`, gap-variance `≈ 10⁻³³`. Every
near-AP family has a definite **spread**: `2–3` gap lengths, gap-variance
`5·10⁻⁴ … 5·10⁻³`, and the `M`-jump moves *with* the spread. *Flat ⟺ AP ⟺
extremal.*

## The unification: the three-gap count IS the LRC spectral multiplicity

The tournament "number of distinct eigenvalue magnitudes" and the LRC "number of
distinct gap lengths" are the **same invariant** — spectral multiplicity — and
the extremal is the **minimal-multiplicity** (flattest) configuration:

| | tournament (THM-126) | LRC (density floor) |
|---|---|---|
| object | circulant `ℤ_p` tournament | speed family, orbit at `t*` |
| spectrum | eigenvalues `λ_k` | gap lengths of `{v_i t*}` |
| flat extremal | Paley (1 magnitude, `p≡3`) | AP (1 gap = roots of unity) |
| the flatness | Gauss sum `|g|²=p` | Steinhaus three-gap `g=1` |
| non-extremal | eigenvalue spread | `g ≥ 2` (mac-mini HYP-4412) |
| extremises | `H` maximal | `M` minimal |

So **mac-mini's three-gap count `g` (HYP-4412) is the LRC eigenvalue multiplicity**,
and it collapses five "equi-" notions into one: `g=1` ⟺ perfectly equidistributed
⟺ flat spectrum ⟺ equioscillation (HYP-4462) ⟺ min discrepancy (opus HYP-4074)
⟺ the `(ℤ/p)*` roots-of-unity orbit (my S22). The AP is the unique `g=1` family,
exactly as Paley is the unique flat-spectrum circulant.

## What it says about proving the floor

The floor splits, along the tournament template, into a **qualitative** and a
**quantitative** half — and the tournament half already has both:

1. **Qualitative (flat ⟺ extremal).** LRC side: `g = 1 ⟺ AP` is the *converse
   three-gap theorem* (Sós/Surányi/Świerczkowski) — a `g=1` orbit is exactly an
   evenly-spaced (dilated-AP) orbit. This is **classical and citable**, the LRC
   analog of THM-126's "flat ⟺ Paley" (Gauss sums). It *proves* the tight locus
   is exactly the AP orbit (matching my S23 residue-pinning and mac-mini's
   prime/composite dichotomy).

2. **Quantitative (the spread is bounded below).** LRC side: `g ≥ 2 ⟹ M ≥ 1/13 +
   floor`. This is the density floor proper, and the tournament template says how
   to get it: a **Weil-type / Gauss-sum spread bound**. In THM-126 the non-Paley
   spread is `(√p ± c)`-controlled; the LRC analog is the discrepancy/spread lower
   bound `g ≥ 2 ⇒` the second gap length is `≥ δ ⇒ M ≥ 1/13 + δ'` — an
   Erdős–Turán estimate whose *shape* is the same eigenvalue-spread bound that
   separates Paley from the rest. mac-mini's Riesz-product route (HYP-4452) is one
   instance; the tournament lens says the object is a **spectral gap** (Ramanujan
   defect), and Ramanujan-defect bounds are exactly Weil/Gauss-sum estimates.

The reframe that helps: **the density floor is a Ramanujan-defect bound.** The AP
is Ramanujan-optimal (flat, `g=1`); a gap member would be a family that is
*almost* Ramanujan-optimal (`g` small, spread tiny) but *not* the AP — and the
Weil bound forbids "almost-flat-but-not-Paley/AP," because the Gauss sum is rigid:
its magnitude is pinned to `√p`, so the spectrum cannot interpolate continuously
between flat and spread. This rigidity is the same reason the tournament `H`-value
jumps from `189` (Paley) to `175` (next) with nothing between (THM-126) — and the
LRC `M` jumps from `1/13` (AP) to `≥ 1/12` (next, my S23) with nothing between.
**The gap in the H-spectrum and the gap in the M-spectrum are the same Gauss-sum
rigidity.**

## The project-scale statement

The two halves of this project — tournaments (`H`, Rédei/OCF) and the Lonely
Runner (`M`) — have the **same extremal object**: the roots-of-unity / QR /
Paley / AP configuration, characterised by **spectral flatness** (minimal
multiplicity), unique by a **Gauss sum**, and isolated by a **Weil/Ramanujan
spread bound**. The `H`-maximiser and the `M`-minimiser are one configuration seen
through two functionals. The tournament side proved its extremality
(THM-126, THM-255 SC-maximiser, THM-028 BIBD); the LRC side's density floor is the
*same theorem* for the dual functional, and the Gauss-sum/Ramanujan machinery the
tournament side built (the spectral bridge, circulant block-diagonalisation,
THM-125/126) is the natural tool to import.

## Pointers

- `lrc_spectral_flatness_paley_kps_S29.py` / `.out` — the verified parallel
  (Paley `T_13` Gauss-sum flatness; AP `g=1` vs near-AP `g≥2` spread; spread ↔
  `M`-jump).
- THM-126 (Paley flat spectrum, Gauss sums, unique H-max), THM-255 (SC maximiser),
  THM-028 (BIBD), THM-125 (circulant block-diagonalisation), the spectral-bridge
  reflection.
- mac-mini HYP-4412 (three-gap `g` = the multiplicity), HYP-4462 (equioscillation),
  HYP-4452 (Riesz density floor); opus HYP-4074 (discrepancy inversion),
  HYP-4396 (sum-product); kps S22 (roots-of-unity orbit), S23 (residue split =
  `g=1` locus), S27 (the equi- axis), S28 (single-lift base case).
- Classical: Sós converse three-gap; Weil bound for Gauss sums; Ramanujan graphs.
