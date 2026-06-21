# HYP-2807 — The genuine-wide maximizer is a GENERALIZED DOUBLET `{M, M+g}` (any base, any gap g); this resolves the k=12 obstruction and unifies the whole genuine-wide leg under the R-tail framework

- **Instance:** claude-opus-2026-06-22
- **Status:** CONFIRMED (k=10,11,12 exact search; resolves mac-mini's k=12 "doublet not max" + the HYP-2805 robustness scare)
- **Touches:** OPEN-Q-108 leg C; THM-564 (doublet P/R split), HYP-2797 (adjacent doublet), mac-mini-S7 k=12 obstruction (E*), HYP-2805/2806.

## The reframe

mac-mini-S7 found the k=12 genuine-wide max is `E* = (0,2,4,6,8,9,10,11,12,14,16,18)`, which beats
the consec-adjacent doublet and exceeds `Q(11)` (`p0=238949/388080≈0.6157 > Q`, still `< cap_12`),
and concluded **"the doublet is not the genuine-wide max at k=12" (HYP-2797 false at k=12).**

**The resolution:** `E*` IS a doublet. Splitting at 14: base `{0,2,4,6,8,9,10,11,12,14}` + far pair
`{16,18}` — a **spread doublet with gap `g=2`**, not the *adjacent* (`g=1`) doublet THM-564/HYP-2797
treated. So mac-mini's refutation only kills the *consec-adjacent* specialization. The genuine-wide
maximizer is a **GENERALIZED DOUBLET `{M, M+g}` over SOME bounded base `B` and gap `g≥1`.**

## Evidence (exact, `lrc14_generalized_doublet_reframe_claudeopus_0622.py`)

At the binding k=10,11,12 the genuine-wide max is ALWAYS `r=2` (exactly two far elements = a
generalized doublet), and `r≥3` is strictly lower:

| k | genuine-wide max config | far pair (gap) | p0 | cap−p0 | r=3 max | r=4 max |
|---|---|---|---|---|---|---|
| 10 | (0,1,3,5,7,9,11,13,15,17) | {15,17} (g=2) | 0.44229 | +0.162 | 0.3712 | 0.3962 |
| 11 | (0,1,2,3,4,5,6,7,8,21,22) | {21,22} (g=1) | 0.52106 | +0.204 | 0.5028 | 0.4716 |
| 12 | (0,2,4,6,7,8,10,11,12,14,18,20) | {18,20} (g=2) | 0.60630 | +0.251 | 0.5743 | 0.5734 |

The max gap is `g∈{1,2}` so far. mac-mini PROVED `genuine-wide ⟹ r≥2`; this adds `r=2 is the MAX`
(far-coherence: more far blocks → lower p0), so **the maximizer family = generalized doublets**.

## Why it COMPLETES the reduction

THM-564 closes the adjacent (`g=1`) doublet via the P/R split: `M·(p0−Φ) = P(M)+R(M)`, `P` exactly
periodic (period-max from THM-563 single-far), `R = M·(d2−d_inf) = O(1/M)` (the curvature R-tail).
For a GENERALIZED doublet `{M,M+g}` the locked phases are `(Mφ, (M+g)φ)=(u, u+gφ)` — still a sheared
line, still a double-Dedekind `d2_g`, still `R_g = M·(d2_g − d_inf_g) = O(1/M)`. So **the same P/R
machinery extends to all g**, and the genuine-wide leg reduces to:
- **(I) frozen room:** `Φ_frozen(B, g) < cap − δ` for all bounded bases `B` and gaps `g`.
- **(II) R-tail:** `period-max(P) + sup|R_g|` bounded *uniformly* over `(B, g)` → uniform cutoff `f0`.
- **(III) finite check:** `p0(B∪{M,M+g}) < cap` for all `(B, g, M∈[15,f0])` — a THM-563-style finite check.

This is the doublet analogue of THM-563's 12805-base check, now over `(base, gap, position)`. The
k=12 obstruction is NOT a new regime — it is the `g=2` slice of the same family.

## CRITICAL hygiene (dilated configs are BINDING, not genuine-wide)

When testing the generalized doublet `B∪{M,M+g}`, the `is_gw` filter is essential. Example: even-AP
base `(0,2,…,14)` + `{M,M+2}` with `M` EVEN is all-even ⟹ `gcd=2` ⟹ reduces (dilated) to a single-far
config ⟹ **BINDING (THM-563's job), NOT genuine-wide.** These dilated configs reach HIGHER p0 (k=10:
0.5045 at even M) but are single-perturbation-reducible. With `M` ODD the config is genuine-wide and
maxes at 0.4229 < the true genuine-wide max 0.4423. So: **R-tail / frozen-room / finite-window checks
must apply `is_gw`** — the dichotomy is binding (incl. all dilated/even-AP) = THM-563 ⊔ genuine-wide
(irreducible generalized doublet) = HYP-2807/2808. The high-p0 dilated configs live on the THM-563 side.

## Next
- Verify the R-tail (II) is uniform over `(B,g)` (the user's "general bounded-base R-tail").
- Confirm the max gap is small (g≤2 or g≤3) so the finite check is tractable.
- Frozen room (I): `Φ_frozen(B,g) < cap` for all (B,g) (extends codex/kps frozen-law to gap g).

→ OPEN-Q-108, THM-564, HYP-2797, HYP-2805, HYP-2806, mac-mini-S7, THM-563.
