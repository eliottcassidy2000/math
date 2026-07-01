# The avoided-arc-edge argument generalizes to a tight-locus result: M(S) ≤ 1/n for every DIFFERENCE-CLOSED set S (all pairwise |s−s'| ∈ S), because the min-gap pigeonhole forces a runner-difference d with ‖dt‖<1/n and difference-closure makes d itself a runner sitting inside the danger arc [−1/n,1/n]; and the UNIQUE primitive difference-closed (n−1)-set is the AP {1,…,n−1} (elementary induction) — so the dilated APs d·{1,…,n−1} are the difference-closed family and are ALL provably tight (M=1/n), giving one rigorously-characterized component of the LRC tight locus, while Goddyn–Wong {1,…,11,13,24} is tight but NOT difference-closed (12=13−1 absent) — the residual that OPEN-Q-108 still must handle

*opus-2026-06-30. Owner: work the remaining LRC proof rigorously; search old work. After closing the THM-591
upper bound, the same avoided-arc-edge trick generalizes and cleanly proves the AP-component of the tight locus.*

## The generalization (homogeneous, c=0)
> **If `S` (with `|S| = n−1`, `1 ∈ S`) is difference-closed — `|s − s'| ∈ S` for all `s,s' ∈ S` — then
> `M(S) ≤ 1/n`.**

**Proof.** Suppose `min_{v∈S} ‖vt‖ > 1/n` for some `t`. Then no runner lies in `[−1/n, 1/n]`, so the `n−1`
runners lie in the complementary open arc of length `1 − 2/n = (n−2)/n`. Their `n−2` consecutive gaps sum to
`< (n−2)/n`, so the smallest is `< 1/n`; it equals `‖(s−s')t‖ = ‖dt‖` for two consecutive runners `s, s'`,
where `d = |s−s'|`. By **difference-closure `d ∈ S`**, so `dt` is a runner with `‖dt‖ < 1/n` — inside the
avoided arc `[−1/n,1/n]`, contradicting the hypothesis. ∎

This is the `c=0` case of THM-591 §B, with the AP's "differences are `{1,…,n−2} ⊆ S`" replaced by the abstract
hypothesis "differences ⊆ S." Combined with LRC (`M ≥ 1/n`), a difference-closed set has `M(S) = 1/n` — it is
**tight**.

## The unique primitive difference-closed set is the AP
> **The only difference-closed `(n−1)`-set containing `1` is `{1, 2, …, n−1}`.**

**Proof.** `1 ∈ S`. Let `a_1 < a_2 < … < a_{n−1}` be `S`. `a_2 − a_1 = a_2 − 1 ∈ S` and `< a_2`, so
`a_2 − 1 = a_1 = 1`, giving `a_2 = 2`. Inductively, if `a_1,…,a_{k−1} = 1,…,k−1`, then `a_k − 1 ∈ S` and
`< a_k` so `≤ a_{k−1} = k−1`, giving `a_k ≤ k`; and `a_k > a_{k−1} = k−1`, so `a_k = k`. ∎ *(Verified by
exhaustive search n=5,6,7,8: the AP is the only one.)*

So the difference-closed sets are exactly the **dilated APs `d·{1,…,n−1}`** (`gcd`-scale the AP). All are tight
(`M = 1/n`, dilation-invariant), all provable by the argument above. This is a **clean, rigorous
characterization of the AP component of the tight locus.**

## The residual (Goddyn–Wong, the OPEN-Q-108 hard part)
The tight locus at `n=14` is conjecturally `{AP {1..13}, Goddyn–Wong {1,…,11,13,24}}`. The AP is now
provably-tight via difference-closure. But **`GW = {1,…,11,13,24}` is tight (`M=1/14`, verified) yet NOT
difference-closed** — e.g. `12 = 13−1 ∉ GW` (also `14,15,…` absent). So its min-gap difference can be `12`,
which is **not a runner**, and the avoided-arc-edge argument fails: nothing sits at `12t` to catch the danger
arc. **GW's mechanism (explored):** it binds at `t*=1/14` via runners **`1` and `13`** — which sandwich `0`
at `±1/14` (their difference is the absent `12`) — mimicking the AP `{1..13}` there. The dangerous regime
`‖12t‖<1/14` is covered by a **patchwork**, not one runner: `24 = 2·12` catches it when `‖12t‖<1/28` (then
`‖24t‖<1/14`), and for `‖12t‖∈[1/28,1/14)` a *shifted* runner (`13`, `11`, …) lands in the arc by
number-theoretic coincidences of `t`. Verified: over 4383 times with `‖12t‖<1/14`, **0** leave the danger arc
empty. So GW's `{13,24}` do cover the absent-difference gap — but by a delicate `t`-dependent patchwork, not a
clean rule. **This is exactly the residual OPEN-Q-108 must handle**: proving the tight locus finite means
ruling out all tight sets beyond the dilated-AP family, and the non-difference-closed ones (the GW family),
whose tightness rests on such patchwork coverings, are the hard case.

## What this buys OPEN-Q-108
- **The AP component is off the table** — rigorously characterized (difference-closed ⇔ dilated AP) and
  provably tight, by an elementary argument (no three-gap counts).
- **The hard core is sharpened**: the tight locus is `{dilated APs}` (fully controlled) `∪ {non-difference-
  closed tight sets}`. The uniform-fattening lemma need only rule out finitely many of the *latter*. The
  GW-family's tightness rests on covering the "absent-difference" gaps (`12t` etc.) — a concrete, finite
  structural condition to attack.
- **A candidate tool**: the avoided-arc-edge argument is a per-`t` local obstruction. For GW, the obstruction
  fails at `t` where `‖12t‖<1/14`; those `t` are a thin set — checking that GW's `{13,24}` covers them is a
  finite verification, a possible route to GW's tightness proof.

## Status
- **PROVED (opus, elementary):** difference-closed `⇒ M(S) ≤ 1/n`; the unique primitive difference-closed
  `(n−1)`-set is the AP; hence dilated APs are the difference-closed family, all tight. Verified n=5..8
  (uniqueness), dilated APs `M=1/n` (n=7), GW tight & not-difference-closed (n=14).
- **Toward OPEN-Q-108:** the AP tight-component is rigorously handled; GW (non-difference-closed) is the
  residual; the danger-arc obstruction localizes GW's tightness to covering absent-difference gaps.
- **Scope (honest):** this does NOT prove LRC(14) or the fattening lemma; it removes the AP family from the
  tight-locus problem and pinpoints the GW residual.

Related: THM-591 §B (the avoided-arc argument), UPPER-BOUND-CLOSED-… (the c=0 case), THM-405 (AP extremal),
OPEN-Q-108 (tight-locus finiteness), HYP-3740/3741 (the lowness/rigidity), HYP-3749.
