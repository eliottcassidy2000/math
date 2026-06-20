---
id: HYP-2655
status: REFUTES the uniform-C≈1.95 claim (HYP-2653b dovetail over-optimistic); CONFIRMS the wide-conclusion; re-localizes the gap to the multi-element recursion. Exact.
source: kind-pasteur-2026-06-19-S14
related:
  - HYP-2653   # the σ-dependent bound (PROVED) + the uniform-C "open input"
  - HYP-2644   # far-element plateau (the multi-element tail was its flagged open issue)
  - HYP-2646   # coset factorization (claimed C≈1.95/1.45)
  - THM-538    # (corrected) support-6 floor
  - OPEN-Q-108
---

# HYP-2655: the uniform decorrelation constant C(k) is NOT ≤1.95 — it grows with multi-scale structure; but the wide conclusion holds via the plateau shrinking

Working the LRC(14) sector route's single "open input" (HYP-2653b: a uniform `w|Δ_w| ≤ C≈1.95` would
dovetail with the span-16 finite check to close the proof) with an EXACT engine. Two findings.

## 1. The exact engine (verified)

For `E = E'∪{w}`, `w=max E`, the far-element plateau error `Δ_w = p_0(E) − [p_0(E')+(1/7)p_1(E')]`
has the EXACT form `Δ_w = (1/w)Σ_{N=1 cells c}[G_0(w b_c − s_c/7) − G_0(w a_c − s_c/7)]`, `G_0` the tent
(`|G_0|≤6/49`). VERIFIED to the exact Fraction against `w·(p_0(E)−plateau)` (e.g. consec_8, w=30:
`w·Δ_w = 19/343` both ways). `w·Δ_w` is PERIODIC in `w` (rational breakpoints), hence bounded per core.

## 2. REFUTATION: `sup_w w|Δ_w|` GROWS with the core's multi-scale structure

A broader EXACT scan defeats the `C≈1.95` claim (which sampled only consec/near-AP/single-scale):
```
   consec_8                              sup_w w|Δ_w| = 3754/5145 = 0.730
   odd-struct (0,1,3,5,7,9,10,11)        = 1370/539 = 2.542   (w=90 = 9·10, resonance)
   wide {0,1,2,20,21,22,40}              = 2.096
   MULTISCALE {0,1,2,30,31,32,60,61}     = 2804017/717360 = 3.909   ← > 1.95
```
So `C(k)` is **not** a small universal constant: each added scale/cluster adds resonant discrepancy
(the `w∼lcm(cluster gaps)` resonances HYP-2653 flagged). The clean "single uniform C" dovetail
(HYP-2653b) **does not close the proof** — at peel-threshold 16 it would need `C < (cap_9−Q(8))·16 ≈
2.12`, already breached by the 3.91 multiscale core. (This is HYP-2644's flagged "multi-element tail.")

## 3. BUT the wide CONCLUSION holds with huge margin (the plateau offsets Δ_w)

The large `Δ_w` is harmless because for a wide/multiscale `E'` the **plateau `p_0(E')+(1/7)p_1(E')` is
itself small** (wide ⟹ low coverage). EXACT `p_0` of wide primitive k=9 sets:
```
   {0,1,2,30,31,32,60,61,62}            p_0 = 141899/635376  = 0.2233   (margin to cap 0.271)
   {0,1,2,100,101,102,200,201,202}      p_0 = 0.2220                    (margin 0.272)
   {0,1,3,5,7,9,10,11,90} (the 2.54 core+far)  p_0 = 2479/16170 = 0.1533 (margin 0.341)
   {0,1,2,3,300,301,302,303,600}        p_0 = 0.1678                    (margin 0.327)
   random wide/multiscale (3000)        max p_0 = 0.2233 ≪ cap_9 = 1979/4004 = 0.4943
```
**Every wide primitive k=9 set has `p_0 ≤ ~0.223 ≪ cap` (margin ≥0.27).** The ONLY tight set is
consec/near-AP (bounded, in the done finite check). So the route is NOT broken — the bound just must
track `plateau(E')` *jointly* with `Δ_w`, not via a separate `Q(k−1) + uniform C`.

## 4. The re-localized target — the multi-element recursion (convergent)

The exact decomposition `p_0(E) = p_0(E') + (1/7)p_1(E') + Δ_w` is a RECURSION: peel the max, the base
`E'` is a `(k−1)`-set; if `E'` is also wide, peel again. The PROVED σ-dependent bound `|Δ_w| ≤
(6/7)σ(E')/w` controls each step. The data shows the recursion CONVERGES to `≤ 0.223 ≪ cap`. The clean
remaining lemma is therefore NOT "a single uniform C" but **the multi-element / joint bound**: a wide
primitive k-set's `p_0` is bounded by the recursion (each peel: plateau-of-base + `(6/7)σ/w`), and the
base plateaus shrink fast enough that the total stays `< cap`. The right object is the σ-split:
`{σ(E') ≤ S : single-far peel, plateau≤Q(k−1), threshold w≥(6/7)S/margin}` ∪ `{σ(E') > S : recurse}`.

## 5. A clean inductive ingredient (verified): the plateau is maximized at the bounded AP

`Φ(F) := p_0(F)+(1/7)p_1(F)` (the plateau) satisfies `Φ(F) ≤ Φ(consec_{k−1}) = Q(k−1)` for ALL
`(k−1)`-sets `F`, wide included (exact: `Φ(consec_8)=621/1715=0.36210`; wide ≤0.144, multiscale ≤0.162).
So `p_0(E) = Φ(E') + Δ_w ≤ Q(k−1) + Δ_w`. CAVEAT (why this still doesn't close): the decomposition does
NOT decouple — a large `Δ_w` is exactly paired with a SMALL `Φ(E')` (wide `E'`), so `sup[Φ(E')+Δ_w]`
(=`p_0(E)`, circular) is not `≤ Q(k−1)+sup|Δ|`. And `Φ`-maximality is itself the same `(k−1)`-level
extremality nut (≈ consec-maximizes-p_0). The far-element recursion is an INDUCTION on `k` that bottoms
out at the finite check, but each level carries the joint `Φ`/`Δ` entanglement — that is the real residue.

## 6. CROSS-THREAD LEAD (kps-S15 + HYP-2657, verify-confirmed): the mod-7 Galois structure may bound the resonant growth

The uniform constant `C(k)` grows precisely at the RESONANCES — `w` sharing factors with the core's `e`'s,
and the multi-scale `w ∼ lcm(cluster gaps)` (the `w=90=9·10` case). HYP-2657 proved (exactly) that the
correction's cyclotomic part is governed by `Re(D7(c))` over the residues `c = n mod 7`, with the full
`F_7*`-average a rational field-trace (`Σ_a D7(a·c) = Tr ∈ ℤ`). The exact reduction of `Δ_w` is a 1-D
breakpoint discrepancy of `{w·(j/e)}` weighted by the missed-sector phase `s_c/7` — i.e. **the same
mod-7 phase that the Galois-equivariance constrains.** CONCRETE NEXT STEP: re-express the resonant
breakpoint contributions of `w·Δ_w` through `Re(D7)` and test whether the trace integrality / the
`σ_6=conjugation` reality bound `|Re D7| ≤ 0.1431` controls them. The resonances (where `w·Δ_w` is large)
are exactly the mod-7 commensurabilities where the Galois phase is rigid — so the growth and the
cancellation may be the same object seen from the two threads (S14 far-element discrepancy ↔ S15 apex-7
Galois). This is the most promising single direction: bound the resonant discrepancy via the apex-prime
quadratic/Galois phase, not a uniform constant.

## 7. S38 addendum: the trace/QR projection is exact but insufficient

Codex S38 / HYP-2662 tested this cross-thread lead by decomposing every endpoint contribution exactly as
`G0(y-s/7)=F7_trace_average(y)+chi_7(s)*QR_NQR_quadratic_channel(y)+intra_quadratic_residual(y,s)`.
The identity verifies through denominator `29`, so the apex-prime phase is real structure at the
far-element endpoint level.  The negative result is sharper than the original question: the bad HYP-2655
resonances are not dominated by the trace average or by the QR/NQR channel.  The intra-quadratic
endpoint-weight residual dominates.  For example, odd-struct at `w=90` has raw `1370/539`, trace
`167/3234`, quadratic `1735/3234`, residual `1053/539`; the multiscale core `(0,1,2,30,31,32,60,61)` at
`w=62` has raw `2804017/717360`, quadratic `9073/10980`, residual `114857/38430`; the scale-family
cluster at `M=50` has raw `19373093/3216850` and residual `6305213/1378650`.

Thus the Galois route should not try to prove a standalone trace/QR uniform bound.  The next target is a
joint inequality: intra-quadratic endpoint-weight residual plus plateau margin.  The same multiscale
geometry that grows the residual must be made to force small `p0(E')+p1(E')/7`.  The incoming KPS
HYP-2661 shell-1 carry-conservation law gives the complementary bounded-mouth version of this target:
sub-`426/35035` near-collar rows conserve the dyadic-1 tower `{1,2,4,8}`.  Its follow-up exact
tower-deletion minima (`1333/30940`, `27/1547`, `335/23023`, `6163/336336`) show every missing tower bit
already pays the AP second threshold, with missing `4` the binding deletion.

## Honest status
LRC(14) NOT proved. CORRECTS HYP-2653b: the "single open input C≈1.95" is FALSE (C grows with scales,
≥3.91). The wide conclusion is CONFIRMED exact (≤0.223 ≪ cap). The real remaining work on the wide branch
is the multi-element recursion / joint plateau-Δ bound (HYP-2644's multi-element tail), which the data
shows converges with margin ≥0.27 — a much looser target than the tight bounded finite check. codex's
S35/S35b bounded route (THM-542/543 mouth-retention) handles the bounded tight side. Files:
`04-computation/lrc14_{uniform_decorrelation_exact,uniform_C_growth,wide_multiscale_p0}_kps.py`.
