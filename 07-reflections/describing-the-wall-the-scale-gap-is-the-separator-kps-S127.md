# Describing the wall — the scale gap is the separator, and it is the dilated-interval extremum

*kind-pasteur-2026-07-10-S127. Owner: "formalize the strict rounding identity; work on describing the
wall better." The identity is formalized (`LRCStrictRuler`, kernel-pure). This note is the description —
and it opens with a correction of the thing I was about to write.*

---

## The wall, stated

After `LRCStrictRuler` the remaining content of LRC(14) is **one integer statement**:

> **`StrictlyLiveSupply`** — every residual family `v` admits a modulus `q` and a multiplier `p`,
> `0 < p < q`, with all thirteen residues strictly inside the band:
> `q < 14·((vᵢ·p) mod q) < 13q` for every `i`.

`lrc14_of_strictlyLiveSupply : LRCUpTo13 → StrictlyLiveSupply → LRC14Statement`, axioms
`[propext, Classical.choice, Quot.sound]`. No measure theory, no continuum, no Fourier. The integrality
of the residues supplies the margin for free: `q < 14r` forces `q + 1 ≤ 14r`, so `ε = 1/(14q)` works
uniformly, with no minimum over `i`.

## The correction I nearly published

I was about to describe the wall as *"covering ⟹ μ(S) > 0"* — the covering hypothesis excluding the
tight locus. **That is false.** Exact rational computation (`lrc14_tight_locus_anatomy_kps_S127`):

```
2·{1,…,13} = {2,4,6,…,26}      covering: YES      μ: 0      primitive: NO      ratio: exactly 13
```

The dilate of the tight AP **is covering** and has **measure-zero safe set**. Covering does not exclude
the tight locus. (klein-S206's statement is safe — it says *primitive* covering — but the loose version
is wrong, and I would have written the loose version.)

## What actually separates them: the scale gap

Every tight family — every dilate `c·{1,…,13}` — has

$$\frac{\max v}{\min v} = \frac{13c}{c} = 13, \qquad \text{exactly}.$$

And the residual's `GapFamily` hypothesis is precisely `¬(∀ i j, |vᵢ| ≤ 13|vⱼ|)`, i.e. **ratio strictly
greater than 13**. So the tight locus is excluded from the residual not by covering, not by primitivity,
but by the **scale gap** — and it is excluded on a knife-edge: the tight families sit at ratio `= 13`,
and `spread13` dispatches everything at ratio `≤ 13`. The boundary of the `spread13` branch *is* the
tight locus.

That is not a coincidence of the assembly's design. It is the geometry of the problem.

## The three faces of the tight locus are one object

| face | characterization |
|---|---|
| **measure** | `μ(S) = 0` — the safe set is isolated points, not intervals |
| **scale** | `max = 13·min` — ratio exactly 13 |
| **additive** | `S` is a **dilated interval** `c·{1,…,13}` — the `E₃`-extremal |

The third is mine: `LRCSchurRigidity.schurCount_eq_choose_iff_dilated` says `E₃(S) = C(13,2)` iff `S` is a
dilated interval, and `LRCE3Budget.dilated_max_eq_card_mul_min` says a dilated interval has
`max = card · min = 13·min`. So the `E₃`-extremum, the ratio-13 locus, and the `μ = 0` locus are the same
set. The scale gap `GapFamily` is exactly "strictly off the `E₃` extremum," and
`LRCE3Budget.E3_lt_choose_of_gap` already proves the qualitative half of that in Lean:

> scale gap ⟹ `E₃(S) < C(13,2)` — strictly off the extremum.

## So the wall is a quantitative stability statement

The qualitative statement is proved. The wall is its quantitative shadow:

> **Crossing the scale gap strictly forces `μ(S) > 0`.**
> Equivalently: off the dilated-interval extremum, the safe set has interior, not just points.

This is exactly the gap I named in S126 on the Freiman ladder: *deficit `> 0` is easy; how much deficit
buys how much structure is hard.* Here the deficit is `C(13,2) − E₃` and the structure bought is `μ`.
The wall is the stability constant relating them.

## Why it resists

- **No absolute bound exists.** The off-line/character sum is irreducibly signed (klein-S222, nine
  independent confirmations; my own S124 correction). Anything that takes absolute values diverges.
- **The exact expansion does not truncate.** klein's Möbius reconstruction
  `μ = (6/7)¹³ + (6/7)¹¹Σc₂ + (6/7)¹⁰Σc₃ + R≥₄` is *alternating with order-one terms* on
  relation-stacked families — for the deep well, layer₃ is `−0.50` against a total of `+0.024`. You
  cannot bound the deficit by bounding a prefix.
- **klein's transfer is not a way around it.** `|LM(q) − q·μ(S)| ≤ K(S) ≤ Σvₗ` closes the *modulus*
  side, converting `μ > 0` into explicit live rulers at every `q > Σv/μ`. But it consumes `μ > 0`. It
  cannot produce it.

## What the wall does *not* require

Two things that look like they should be needed, and are not:

- **No uniform floor.** `μ` on covering families gets small (measured minimum `0.0367` on `[1,22]`;
  the deep well is `0.0239`), and a uniform `μ₀` over the infinite residual class may not exist. It is not
  needed: klein's Corollary 1 makes the certification *per family* a priori complete — any `μ(S) > 0`
  gives live rulers at all `q > Σv/μ(S)`, and the bank below that threshold is finite and explicit.
- **No measure theory.** `LRCStrictWitnessFloor` shows a single strict-margin point yields the floor by
  the reverse triangle inequality, and `LRCStrictRuler` shows an integer ruler yields the point. The
  whole chain from an integer certificate to LRC(14) is now kernel-pure and elementary.

## The wall, in one sentence

> **Does a speed set whose ratio exceeds 13 — one strictly off the dilated-interval extremum — necessarily
> admit a modulus at which all thirteen residues sit strictly inside `(q/14, 13q/14)`?**

Everything else about LRC(14) is proved, kernel-pure, in Lean. This is the open case of the Lonely Runner
Conjecture, and it is a stability question about the extremum I spent this session characterizing.

*Files: `LRCStrictRuler.lean` (the strict rounding identity, sorry-free, kernel-pure),
`lrc14_tight_locus_anatomy_kps_S127.py`/`.out` (the exact-`μ` sweep and the `2·{1,…,13}` correction).
Builds on klein THM-685, `LRCSchurRigidity`/`LRCE3Budget` (kps), `LRCStrictWitnessFloor` (kps). See
[[the-residual-is-one-measure-floor-kps-S127]] and
[[the-two-axes-share-a-threshold-e3-peel-ladder-kps-S126]].*
