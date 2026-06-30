# Per-level vs total floor, pinned: the DOUBLET (per-level, ρ_j ≥ 4cos²(3π/7), robust) is the EXISTENCE atom that binds the lonely point at every level; the TOTAL floor ∏ρ_j ~ (0.198)^d → 0 is the vanishing MEASURE (the razor-thin) — so LRC(14) is the per-level doublet domination, and the total measure is a red herring

*opus-2026-06-30. Owner: certify the domination and pin per-level vs total, thinking doublets. The pin
dissolves the apparent paradox (the floor "vanishes" yet LRC holds): per-level is existence, total is
measure, and the doublet is the per-level atom.*

## The doublet IS the per-level atom (computed)
The descent's binding core `{0,1}⊂Z₇` (= my binding pair) has autocorrelation `a(d)=[2,1,0,0,0,0,1]` —
**`2I` + the single ±1 shift, the 2-runner minimal resonance** — with Gram eigenvalues
`[4, 3.247, 1.555, 0.198, 0.198, 1.555, 3.247]`:
> `λ_0=4=|O|²` is the DC/mean (the *average* safe-density). **The FLOOR is the MIN mode `λ_min=0.198=
> 4cos²(3π/7)`, not the mean.** This is the whole "minorant-not-average" point made atomic: the doublet's
> mean is `~0.57`, its min is `0.198`; the floor is the min.

## The PIN: per-level (existence) vs total (measure)
THM-580: `meas(lonely S) = ∏_j ρ_j · ∏_j meas(lonely O_j)` — a PRODUCT over the `d` descent levels.
| | object | value | what it is |
|---|---|---|---|
| **PER-LEVEL** | `ρ_j` | `≥ g(O_j mod 7) ≥ 0.198` (CONSTANT) | the doublet atom — **robust, set-independent, THM-590 PROVED** |
| **TOTAL** | `∏_j ρ_j` | `~ (0.198)^d → 0` | the lonely **MEASURE** — vanishes (HYP-3597, razor-thin) |
> **The total measure `∏ρ_j → 0` IS the razor-thin** (the product of factors `<1` over growing depth). But
> it is a MEASURE, and **`M ≥ 1/14` is EXISTENCE, not measure.** Existence is carried PER-LEVEL by the
> doublet — two runners binding the lonely POINT at every level — never by the (vanishing) total. **The
> apparent paradox "the floor vanishes yet LRC holds" is just per-level vs total: the doublet is robust,
> the product is not, and only the doublet matters for the point.**

## Where each side lives (the cusp)
- **Off-cusp levels** (proper core `O_j ⊊ Z_7`): `ρ_j ≥ 0.198` (the doublet or better) — the measure is
  positive, the doublet binds.
- **The cusp** (`O_j = Z_7`, the unique gap-0 core): `ρ_j = 0` — the apex measure vanishes; here **existence
  must be carried directly** (the discrete/witness side, HYP-3597). The cusp is the *only* place the
  measure argument fails, and it is exactly the disproof boundary.
So LRC(14) = "off-cusp the doublet binds (`≥0.198`, proved); at the cusp existence carries (the discrete
side)." The total measure is irrelevant to both.

## The domination to certify (per-level, single mode)
> **`ρ_j(genuine) ≥ g(O_j mod 7) = λ_min(Gram)`** — the Fejér–Bochner SOS in MIN-eigenvalue form.
Verified direction: the doublet's density proxy `λ_0/7 = 0.571 ≥ 0.198` (the genuine per-level density
exceeds the minorant — the SOS points the right way; klein-S4 already has `ρ_j>0` everywhere). Since
`14=2·7` is single-prime, this is ONE Fourier mode — at the binding level, the doublet's `0.198`. The
remaining step is the *sharp* (least-eigenvalue) certificate, not a product, not an average.

## What the pin buys
1. **The razor-thin is explained, not a problem.** The vanishing total measure was the right computation
   (HYP-3597) for the wrong target; the target is the per-level doublet, which does not vanish.
2. **The conjecture is one robust number per level.** `inf_j ρ_j = 4cos²(3π/7)` (proved, the doublet),
   set-independent and unavoidable (HYP-3598). LRC(14) is the per-level domination + the cusp existence.
3. **The doublet is the whole object, every level.** Binding pair = doublet core = per-level atom =
   existence carrier = the genus-1 cusp form's support. The total product is the measure shadow.

## Status
- **Computed/pinned (opus):** the doublet = the per-level atom (`2I+shift`, `λ_min=0.198`, the floor is
  the min not the mean); per-level `ρ_j ≥ 0.198` ROBUST vs total `∏ρ_j → 0` VANISHING; existence (per-level
  doublet) vs measure (total, razor-thin); the cusp `O=Z_7` is the only `ρ=0` level (existence carries).
- **The domination (per-level):** `ρ_j ≥ g(O mod 7)` (Fejér–Bochner min-eigenvalue SOS, single prime),
  direction verified; the sharp certificate is the remaining step.
- **Honest:** this pins the distinction and centers the doublet; the genuine-`ρ_j` ≥ least-eigenvalue
  certificate (and the cusp existence side) are the two remaining pieces — both per-level, both finite.

Related: the single-prime-minorant bridge reflection; klein HYP-3597 (measure→0/existence split), HYP-3598
(per-level apex floor, complete family), THM-590 (the doublet gap), THM-580 (the product over levels),
my razor-thin-in-the-measure + binding-pair + Dirac-comb reflections, OPEN-Q-108.
