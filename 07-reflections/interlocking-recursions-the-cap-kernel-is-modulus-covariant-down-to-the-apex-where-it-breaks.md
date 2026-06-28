# Interlocking recursions: the cap kernel is modulus-covariant (pure scale) down to the apex n/2=7, where it breaks — the three modes concentrate at the Eisenstein fold

*mac-mini-2026-06-28-S75c. The owner: think of interlocking recursions, the tiling/half-tiling models, and
moment-order depth, and leverage them for LRC proofs. Integrating the three-mode atlas (HYP-2900/2901), the
half-tiling DP (HYP-2690), the moment-depth dichotomy (HYP-2900 corrections), and my three-gap kernel (S75b)
into ONE interlocking picture — and a tested finding that locates the entire difficulty at the apex fold.
Builds on [[how-the-route-sharpened-and-the-three-gap-stern-brocot-recursion-in-the-cap-kernel]], HYP-2900/2901, HYP-2690.*

## Three ORTHOGONAL interlocking recursions (one per coordinate of the cap)
The cap is a triple-recursive object; the project's three "modes" (HYP-2900/2901) are recursions on three
independent coordinates, and my recent work is the Legendre leg:

| mode | recursion on | character | object | my work |
|---|---|---|---|---|
| **Möbius** | the moment ORDER `j` (inclusion-exclusion depth) | `μ` (`+++---+`) | the IE skeleton, all orders | the cap = IE (S75) |
| **Eisenstein** | the MODULUS `n` (even fold `2q→q`) | `χ₃` (`++-`) | `14→7→2` (exposes apex 7) | the apex-14 boundary |
| **Legendre** | the SPEED ratio `a/b` (three-gap / continued fraction) | `χ₇` (`++-+--+`) | the odd 3-set Venn at 7 | the kernel `K(a,b)` (S75b) |

Composition (HYP-2901): `14 = 2·7 = Eisenstein(14→7) ∘ Legendre(7)` on the Möbius floor; computed by the
**half-tiling `(β,τ)` DP** (HYP-2690), with the moment-depth **dichotomy** (HYP-2900): tight→Bonferroni-3,
slack→equidistribution `(6/7)^r`. The single-arc peeling (S75) is the depth-collapse where the Legendre
three-gap is trivial (speed 1).

## The tested interlock: the kernel is MODULUS-COVARIANT, and the apex is where it BREAKS
(`the modulus-fold scan`)
- **Legendre is modulus-covariant for `n ≥ 2·max(speed)`:** `K^{(n)}(a,b) = (2/n)·h(a,b)`, with `h(a,b)` the
  **modulus-free three-gap invariant** (gap-counts identical at `n=14,28`). So the threshold `1/n` is a pure
  SCALE; the **Eisenstein fold `n→2n` is exactly `×½`** (verified ratio `K^7/K^14 = K^14/K^28 = 2`).
- **The covariance BREAKS at the apex `n/2`:** at `n=7` the fold ratio is `8/3` (not 2) for `(3,5)`, and the
  gap-counts deviate (`b=5`: `1,2,4,6` vs the stable `1,2,3,4`). The break is at **speeds `> n/2 = 7`** — exactly
  the antipode region where my `K(a,13)=(2a−1)/(91a)` deviates and the binding constants `cap_8, cap_9` live.

> **The interlock, precisely:** for speeds `≤ 7` (`= n/2`) the kernel is the clean modulus-free three-gap (low
> moment-depth, Eisenstein-trivial); for speeds `8..13` (`> n/2`, the antipode half) the kernel deviates — and
> THAT deviation is the entire LRC(14) difficulty. **The three recursions all concentrate at the apex `7 = n/2`:**
> the Eisenstein fold lands there, the Legendre three-gap breaks there, and the Möbius moment-depth rises there
> (order-2 → order-3 corrections).

## Moment-order depth interlocks with speed-vs-apex
The moment-order depth (the IE order where the cap closes) is governed by the speed size relative to the apex:
- **speeds `≤ 7`:** clean three-gap, `K=(2/n)h`, low depth — the cap closes at low order (pairwise/Möbius).
- **speeds `8..13` (antipode):** the second tooth catches (the modulus-covariance break), generating the
  order-3+ deviations — higher moment-depth, the binding constants.
- **speed 1 (the extreme small):** the single-arc lemma collapses ALL orders to `1/(7 max)` (the peeling) — the
  deepest depth-collapse. So the **moment-depth interlocks monotonically with the speed ordering**: smallest
  speed → total collapse (peel), middle speeds (`≤7`) → low depth, large speeds (`>7`) → the binding depth.

## kps S31aq convergence: the MOMENT-ORDER DEPTH is `(p+1)/2` (the family law)
kps independently traced the same sharpening and **quantified the moment-order depth** — the owner's exact
concept (HYP-3216): for `LRC(2p)`, the apex moment-order **DEPTH = (p+1)/2**, a VERIFIED family law:
```
   n=6 (p=3): depth 2    n=10 (p=5): depth 3    n=14 (p=7): depth 4 = the cubic wall, FIRST hard
```
via TWO interlocking recursions = `14 = 2·7`: the **7-recursion** = the cyclotomic moment-order ladder
`cap_k = cap_{k-1} + k/91` (Faulhaber, VERIFIED; cyclotomic degree `(p-1)/2`, triangular widths `3,2,1`), and
the **2-recursion** = the 2-adic reflection fold (the apex biquadratic `u⁴−5u²+4 → degree 2 in v=u²`, halving via
`s↔6−s`). This is exactly the user's "moment-order depth," and it **interlocks with my modulus-covariance break:**
my break is at the apex `n/2`, and kps's depth `(p+1)/2=4` is reached precisely there (the antipode half). So:
**moment-order depth `(p+1)/2` (kps, the order axis) = the speed-vs-apex break at `n/2` (me, the speed axis) = the
Eisenstein fold image (the modulus axis).** Three axes, one concentration point, one depth `(p+1)/2`. And kps
credits my S75b three-gap kernel as the **third (Diophantine) recursion** alongside their cyclotomic + 2-adic two.

## The leverage (an interlocking-recursion proof architecture)
1. **Easy half (speeds ≤ 7, modulus-covariant):** the three-gap kernel is a pure scale `(2/n)h`; the cap
   contribution is the clean Legendre three-gap, low moment-depth — handled by the half-tiling DP.
2. **Apex half (speeds 8..13, the break):** the antipode deviations `(2a−1)/(91a)` etc. — the binding
   constants, the order-3 Möbius corrections, the genuine difficulty. This is the **Eisenstein fold 14→7
   image** where the structure breaks.
3. **Analytic dichotomy:** tight (binding far cluster) → Bonferroni-3 at the apex; slack/large →
   equidistribution `(6/7)^r`. The apex break is the binding/tight side.
4. **Base:** the Eisenstein descent `14→7→2` bottoms at LRC(≤6) (proven). The interlocking is the prime-tower
   descent carrying the three recursions together.

## Honest status
- **TESTED:** modulus-covariance `K^{(n)}=(2/n)h` for `n≥2·max(speed)`; Eisenstein fold `=×½`; the break at the
  apex `n/2=7` (the `(3,5)` `8/3` ratio, the `b=5,6` gap-count deviations); the moment-depth ↔ speed-vs-apex
  interlock; the single-arc total collapse at speed 1.
- **SYNTHESIS:** the three modes are recursions on three orthogonal coordinates (order/modulus/speed); they
  interlock and **concentrate at the apex `7 = n/2`**, which is precisely where the kernel's modulus-covariance
  breaks. The LRC(14) difficulty is the antipode half (speeds `>7`), the Eisenstein-fold image of the apex.
- **NOT a proof.** This precisely *locates* the difficulty (the apex-half deviation, bounded moment-depth) and
  organizes the three recursions; it does not close the binding constants. But the recursive structure is now
  described as a triple interlock with a single concentration point (the apex `n/2`). LRC(14) open.

Related: HYP-3232 (this), HYP-3230 (the three-gap kernel), HYP-2900/2901 (the three modes), HYP-2690 (half-tiling DP),
HYP-3227 (single-arc peeling), OPEN-Q-108.
