# The bottom of the n=14 spectrum is two ladders rooted at the AP — with plateaus

**Instance:** klein-2026-07-19-S321 (owner: "look for more structural surprises,
synthesize with bleeding edge from multiple pulls as it changes").
**Evidence:** exhaustive bottom-spectrum census (interval (0, 1/13), heights ≤ 36,
harness v2, 7 families total) + exact-M verification sweep, all predictions OK
(`lrc14_bottom_spectrum_ladders_klein_S321.py` + frozen out). Everything below is
exact rational arithmetic on named families; the census part is exhaustive.

## 1. The complete sub-1/13 spectrum at heights ≤ 36 is SEVEN families

| M | families (w_max) |
|---|---|
| 1/14 | {1..13} (13), {1..11,13,24} (24) |
| 2/27 | {1..9,11,12,13,20} (20), {1..9,11,13,20,24} (24), {1..12,26} (26), {1..9,11,13,20,36} (36) |
| 3/41 | {1..11,13,36} (36) |

Nothing else. No 3/40, no 4/53, no k ≥ 3 stratum value — through height 36 the
attained bottom spectrum is exactly {1/14, 2/27, 3/41}.

## 2. Every one of these families is a rung of a parametrized ladder

- **The 12m-ladder** L_B(m) = B ∪ {12m} over a 12-base B:
  - B₁ = {1..11,13}: L₁(1) = **{1..13} = the AP itself**; L₁(2) = {1..11,13,24}
    (the second floor realizer); L₁(3) = {1..11,13,36} (the 3/41 witness);
    L₁(4) = 4/53; L₁(5) = 1/13; then 6/77, 7/89, … (THM-1230's ladder, re-rooted).
  - B₂ = {1..9,11,13,20}: L₂(1..3) = the three plateau 2/27 realizers;
    L₂(4) = 4/53; L₂(5) = 1/13; L₂(6..) = 6/77, 7/89, … — **the same tail**.
- **The 13s-ladder** K(s) = {1..12} ∪ {13s}: K(1) = **{1..13} again**;
  K(2) = {1..12,26} = 2/27; K(3) = {1..12,39} = 3/40; K(4) = 4/53; K(5) = 5/66;
  K(6) = 6/79 — the Kravitz rungs s/(13s+1), verified exact through s = 6.

**Both ladders are rooted at the AP** (m = 1 and s = 1 both give {1..13}), and they
carry between them every family the exhaustive census found.

## 3. The plateau law (the mechanism, verified 18/18)

> **M(B ∪ {12m}) = max( m/(12m+5), c_B )**, where c_B is the base's best
> far-element-immune pinch (c_{B₁} = 1/14 via the floor; c_{B₂} = 2/27 via the
> base-internal pair (7,20) at q = 27).

Mechanics: the far element 12m **kills the base's q = 12 maximizer instantly**
(12m·(a/12) ≡ 0), collapsing M from M(B) = 1/12 (both bases — verified) to a
two-peak competition: the universal far pinch on the active pair **(5, 12m)** worth
m/(12m+5) (D = m, s = 12m+5 — THM-1269's D = M·s and THM-1291's CF active-leg law
in action; 5 is the shared convergent leg), versus the base's own c_B. Where the
formula would dip below c_B — including where it would violate LRC(14) itself
(m ≤ 2 gives m/(12m+5) < 1/14) — the value **plateaus at c_B**. The ladder
accumulates at 1/12 = M(B) exactly: a concrete instance of G-K Theorem 1.3's
from-below accumulation, with the limit value visibly equal to the deleted-element
spectrum value it converges to.

## 4. The surprises, stated plainly

1. **Height order inverts value order across strata.** h_min(2/27) = 20 <
   h_min(3/41) ≤ 36 < h_min(3/40) ∈ (36, 39] although 1/14 < 3/41 < 2/27 < 3/40.
   The k = 2 mediant 3/41 is realized at a LOWER height than the k = 1 Kravitz
   rung 3/40 despite being the smaller (harder) value. Realization height is
   governed by ladder position (12m vs 13s), not by the value's size.
2. **Realizer multiplicity = plateau length.** 2/27's four realizers are not
   sporadic: three are the L₂ plateau (m = 1, 2, 3) and one is K(2). Any 12-base
   whose internal pinch equals a rung value contributes a plateau of realizers.
   Multiplicity is a mechanical quantity — and if some base's plateau were
   unbounded, the value would be an infinite-multiplicity **density point** in
   G-K's den() sense; here the plateau ends at m = 4 (the formula overtakes), so
   the mechanism produces finite plateaus with computable lengths.
3. **The floor's second realizer is a plateau artifact.** {1..11,13,24} = L₁(2)
   sits at 1/14 only because 2/29 < 1/14 is LRC-impossible — the conjecture's own
   floor truncates the formula. The AP's "sporadic twin" is nothing but the m = 2
   rung clipped by the floor.
4. **The two-ladder picture explains the census wipeout shape.** Everything the
   exhaustive census found below 1/13 through height 36 lies on these two ladders;
   the first genuinely off-ladder object (if any) must appear at height ≥ 37 —
   the B = 45 completion (in flight) decides the next window.

## 5. Bleeding-edge convergence (same day, three agents)

- **mac-mini S126** (three-far stratum, zero survivors): "exit ceiling q = 41
  AGAIN" — my THM-1290 spectroscopy has the exceptions: exactly 50 families
  escape past 41 into witnesses q ∈ [43, 48] at heights ≤ 55, and NONE at
  56–64. The survivor-identity rerun (in flight) hands kps's HYP-8040 escape
  atlas its missing hard cases.
- **kps S128c95** (HYP-8040, small-witness law ¾ proved): the "all-heights law =
  named CRT intersection" is exactly what the census's filter-stack joint
  unsatisfiability (THM-1290 S320 corollary) instantiates at bounded height.
- **opus S410** (rung-theory manifesto, in flight): the two-ladder root-at-AP
  picture + the plateau law + the (5, 12m) universal active leg are direct
  manifesto inputs — the "rung" is not a value but a LADDER POSITION, and values
  are shared between ladders (4/53 = L(4) = K(4); 1/13 = L(5); the AP = both
  roots).

## 6. Cross-links

THM-1290 (+S320 extension) · THM-1230 (the L₁ ladder, now re-rooted at the AP) ·
THM-1235/1268/1269 (slack/D-coordinates; the plateau law refines "slack = D − k") ·
THM-1291 (CF active-leg; the (5,12m) universal leg) · THM-1289/HYP-7930-UPDATE
(from-below accumulation; density points) · kps HYP-8040 (escape atlas) ·
mac-mini HYP-7990 (q=41 exit ceiling) · opus S410 (rung manifesto) ·
kind-pasteur ladder-realization lrc14_ladder_realization_crossN (K(2) mechanism).
