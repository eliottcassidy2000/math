# 169 = 13², and the difference-closure rigidity of `M < 1/13`

*boxeph-2026-07-18-S87. A depth-first exploration of the LRC(14) crux
`M < 1/13 ⟹ ρ ≥ 13` (which, with the PROVED THM-1008 `ρ≥13 ⟹ M≥1/14`, would close
LRC(14)). Outcome: genuine structure and one rigorous rigidity lemma; **the crux is NOT
closed.** Scripts: `lrc_169_dfs_boxeph_S87.py`, `lrc_compact_1over13_hunt` (S85).*

## The honest bottom line

`M < 1/13 ⟹ ρ ≥ 13` is the sole remaining piece of LRC(14) ([[THM-1013-dilated-sieve-compact-floor]]
reduction map). This session did **not** prove it. It is a genuine rigidity (near-tight ⟹ far
element), and the sub-arguments below each stall at the same bridge. Recorded so the next attempt
starts from the sharpened picture, not from scratch.

## 1. The `{1..12}`-core computation — where 182 = lcm(13,14) is forced (PROVED, narrow)

`M({1,…,12} ∪ {13k}) = k/(13k+1)` (maximizer `t = k/(13k+1)`, `val = k`), all `< 1/13`. But
`{1,…,12} ∪ {13k}` is **covering only at `k = 14`**: it covers 13 (via `13k`) but covers 14 only if
`14 ∣ 13k ⟺ 14 ∣ k`. So the single covering member of this pencil is the **deep well** `{1..12,182}`,
`182 = lcm(13,14) = 13·14`. Hence:

> **Any covering family with core exactly `{1,…,12}` plus one element `w` has `182 ∣ w`**
> (`13∣w` from covering 13, `14∣w` from covering 14). So `w ≥ 182`, `ρ = w/12 ≥ 15.17 ≥ 13`, and
> `M ≥ 1/14` by [[THM-1008-lrc13-descent-floor]].

Elementary and rigorous — but only for the exact-`{1..12}` core. It shows the mechanism: **covering
14 forces the coverer of 13 to be far**, because in the tight AP core the only way to also hit 14 is
the lcm, which is large.

## 2. `M < 1/13 ⟺ 13·val < q < 14·val` (PROVED, from THM-999)

Writing the maximizer `M = val/q` (THM-999 / THM-724-pair-sum, `q ∣ v_i+v_j`, `q ≤ 2v_max`):
`M ∈ (1/14, 1/13) ⟺ 13 < q/val < 14 ⟺ 13·val < q < 14·val`. Equivalently the continued fraction of
`M` begins `[0; 13, …]`. Verified: every `M<1/13` covering family found is `{1..12, 182m}` with
`M = 14m/(182m+1) = [0;13,14m]`, `ρ ≥ 15`. And `169 = 13²` is intrinsic — at the deep well
`t = 14/183`, runner `182` lands at `182·14 ≡ 169 (mod 183)`, distance `183−169 = 14 = val`.

## 3. The difference-closure rigidity lemma (PROVED — the mechanism)

> **Lemma.** If `M(V) = val/q < 1/13` (maximizer `t=a/q`), then among the 13 residues
> `r_i = v_i·a mod q ∈ [val, q−val]`, some two lie at circular distance `< val`. Hence there are
> `i ≠ j` with `|(v_i − v_j)·a|_q < val` and `|v_i − v_j| ∉ V` — **`V` is not difference-closed at
> that pair, and the missing difference is resonance-aligned.**

*Proof.* `q < 14val` ⟹ the band `[val, q−val]` has length `q−2val < 12val`; 13 points in it have 12
gaps summing to `< 12val`, so min gap `< val`. That pair's speed-difference `δ = v_i−v_j` has
`|δa|_q = |r_i−r_j|_q < val`; if `|δ| ∈ V` then `|δa|_q ≥ val` (definition of `val` as the min over
`V`), contradiction. ∎

**This is why `1/14` is the floor of difference-closed families and `M<1/13` is not.** Verified: the
AP `{1..13}` sits exactly at the boundary (`q = 14val`, min gap `= val`, closest diff `= 1 ∈ V`,
difference-closed, `M = 1/14`); every `M<1/13` family breaks it (deep well: closest pair `(182,12)`,
`δ = 170 ∉ V`, `|170·14|_183 = 1 < 14 = val`).

## 4. Where every branch stalls — the bridge to `ρ ≥ 13`

The lemma says `M<1/13 ⟹ V` is near-difference-closed but with a resonance-aligned gap. To finish one
must show this **forces a far element** (`ρ ≥ 13`). The intuition (§1) is clear — a tight
(difference-closed-ish) core is AP-like `d·{1,…,k}`, which misses 13 and/or 14, and covering them
without disturbing the core (raising `M` to `≥ 1/13`) forces a coverer at the lcm scale (`≥ 182`),
hence far. But turning "near-difference-closed core + covering" into "`ρ ≥ 13`" is exactly the
**n=12/13 additive-rigidity / Freiman-type** statement (klein's Hamming-radius theorems, HYP-7310).
The pigeonhole gives one aligned gap; the bridge needs the *global* AP-structure, which one aligned
gap does not deliver.

Equivalent restatements of the open bridge (all verified, none proved):
- `M<1/13 ⟹` the smallest multiple of 13 in `V` is `≥ 182` (covered by a *far* element);
- `M<1/13 ⟹ V` is (a dilation of) `{1,…,12}` plus a single far killer;
- `M<1/13 ⟹` the residues are the equally-`val`-spaced AP pattern perturbed at exactly one point.

## Takeaways for the next attempt

- The crux is **additive rigidity**, not a missing witness: the elementary sieve/witness family is
  complete ([[THM-1013-dilated-sieve-compact-floor]]); `M<1/13` needs a *structure theorem*, not a
  new rational time.
- Start from the **difference-closure lemma (§3)** — it is the rigorous kernel. The task is to
  upgrade "one aligned non-speed difference" to "AP core + far element."
- `182 = lcm(13,14)` and `183 = Φ₆(14) = 14²−14+1` and `169 = 13²` are the arithmetic of the deep
  well; the extremal is the continued fraction `[0;13,14]`. Any proof should reproduce these as the
  unique minimizer.
- Cross-link: [[both-lrc14-routes-bottom-on-the-same-multilinear-cancellation-not-one-mollification-klein-S279]]
  (the covering-side multilinear cancellation is this same additive rigidity),
  [[the-resonance-fill-profile-one-lens-for-every-lrc-face-boxeph-S81]] (fill-1 = the AP core's
  under-filled circle), HYP-7355/7358.
