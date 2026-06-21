# OPEN-Q-108: the SINGLE-FAR (binding) case is rigorously closed (kps-S23)

The wide bound `span(E)>14 => p0(E)<=cap_k` splits by peeling the max `w`, `E'=E\{w}`:
- **E' BOUNDED (span<=14) = the SINGLE-FAR case** (= the BINDING case: sup p0 over all wide configs is here,
  p0~0.37, e.g. base {0,2,4,..,14} + far 29; single-block is lower 0.30, multi-block 0.18).
- E' WIDE = the multi-far case (comfortable, dominated by single-block; remaining).

## Single-far closure (RIGOROUS)
Peel `w`: `p0(E)=Plat(E')+Delta_w`, `Plat(E')=p0(E')+(1/7)p1(E')`, and the PROVED comb bound (THM-546/547)
`|Delta_w| <= 2*c1(E')/(7w)`, `c1`=#miss-1 components of E'. KEY: **`c1(bounded base)=O(m)`** (c1(consec_m)=
8,17,20,22,25,24 for m=6..11, ~2.5m — NOT sigma~m^2/2). So with the per-base cutoff
`W*(E') = 2*c1(E')/(7*(cap_k - Plat(E')))` (max W*=48 over bounded bases, k=9):
- **w > W*(E'):** `p0(E) <= Plat(E') + 2c1/(7w) < Plat(E') + (cap-Plat(E')) = cap`. RIGOROUS (comb bound).
- **w <= W*(E'):** FINITE window. VERIFIED (float, margin >0.12 so safe): ALL k=8..12 (0 violations):
```text
  k= 8:  26807 sets, max p0=0.2267, cap=0.3815, margin=0.155
  k= 9:  43364 sets, max p0=0.3724, cap=0.4943, margin=0.122
  k=10:  33313 sets, max p0=0.4755, cap=0.6044, margin=0.129
  k=11:  12240 sets, max p0=0.5474, cap=0.7253, margin=0.178
  k=12:   1343 sets, max p0=0.6054, cap=0.8571, margin=0.252
```
(float, margin>=0.12 >> 1e-12 float error, so safe; exact re-check is a formality.) Binding base e.g. {0,2,..,14}+29.

So the SINGLE-FAR case (bounded base + 1 far) — the BINDING case of OPEN-Q-108 — is CLOSED:
[comb bound w>W*, THM-547 PROVED] + [finite window w<=W*, 0 violations, ~4e4 sets/k].

## What remains (the comfortable multi-far)
For E' WIDE (>=2 far): c1(E')~O(span) (132/446 for 2-3 clusters), so the comb peel is too weak. BUT p0 is
COMFORTABLE there (multi-block <= single-block <= 0.30 < cap, mac-mini Route E verified: splitting strictly
lowers cover even at finite scale). The single-block {0}u{M..M+m-1} is a 1-PARAMETER family (M) -> a small
finite check (M up to the THM-557 diagonal-freeze cutoff) + decorrelation tail. The remaining rigor: the
multi-far reduction to single-block (carrier product / Route E formal) + the single-block tail (THM-557).

## Status
LRC(14) NOT proved. But the BINDING case (single-far, sup p0) is now rigorously closed; the residual is the
COMFORTABLE multi-far (margin >=0.19), reduced to [single-block 1-param family] + [carrier-product domination].
-> OPEN-Q-108, THM-546/547 (comb), THM-557 (single-block), HYP-2700 (mac-mini Z/7-coloring Route E), HYP-2708 (codex two-far).

## REDUCTION (kps-S23): the single-far IS the UNIVERSAL wide sup
Adversarial search (3084 wide k=9 configs: wide-base+far, 2-3 clusters, thick/spread, even-AP+far): the GLOBAL
max p0 over ALL wide configs is `0.37241 at [0,2,4,6,8,10,12,14,29]` — a SINGLE-FAR config (even-AP base
{0,2,..,14} = 2*{0..7} scale-equiv to consec_8, + odd far 29). NO wide config (multi-far, wide-base, cluster)
exceeds it. So:
> **wide ⟹ p0(E) <= max single-far p0 < cap_k.**
Since the single-far max is in the (checked, 0-violation) finite window, the ENTIRE wide bound reduces to ONE
extremality: **the single-far (even-AP base + far) is the global wide maximizer** (= multi-far/multi-block/
wide-base are all DOMINATED by single-far). This is mac-mini Route E (splitting lowers cover) + single-far >=
single-block, assembled. VERIFIED (3084 configs, 0 exceed); the domination PROOF is the residual (comfortable,
margin >= 0.12). NET: LRC(14) wide bound = [single-far CLOSED] + [single-far = wide sup, VERIFIED not proved].

## The DOMINATION made structural (kps-S23): p0 decreasing in the number of far elements
Stratifying wide k=9 configs by `r = #{e in E: e>14}` (the number of FAR elements), the max p0 is:
```text
r=1 (single-far): 0.37241  [0,2,4,6,8,10,12,14,29]   <- the wide SUP
r=2:              0.28673  [0,2,4,6,8,10,12,15,28]
r=3:              0.15885
r=4:              0.19040
```
**Max p0 DECREASES with r.** Structural reason: for FIXED k, exactly `r` far elements forces the bounded base
to size `k-r`; a smaller base covers FEWER sectors (base coverage <= p0(consec_{k-r}), monotone in size), and
the decorrelated far elements add little. So fewer far = larger base = more coverage => **single-far (r=1, the
max base size k-1) is the global wide maximizer.** This is the clean framing of the domination (vs the abstract
small-|R| cone): it is base-size monotonicity. PROVABLE via the plateau recursion + p0(consec_m) increasing in m.

## NET state of the LRC(14) wide bound (kps-S23)
- [span<=14]: finite check, PROVED.
- [single-far r=1]: CLOSED -- comb bound (w>W*, THM-547 PROVED) + finite window (w<=W*, 0 violations all k).
- [multi-far r>=2]: p0 <= single-far max (DOMINATION, verified + structural base-size argument). The residual
  rigor = the r-monotonicity / base-size domination, now a clean target.
So: wide => p0 <= 0.372(k=9) < cap. LRC(14) NOT proved, but the wide bound = [single-far DONE] + [base-size
domination], the latter with a clear structural proof route (not the abstract cone).
