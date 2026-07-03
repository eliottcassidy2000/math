# The c=7 Trichotomy — the wall breaks three different ways

*kind-pasteur-2026-07-02-S24. A reflection on why the seven-runner "wall" in LRC(14)
does not fall to a single mechanism, and why that is the correct picture rather than
a failure of technique.*

## The wall, restated

Seven danger arcs, each of circle-measure `1/7`, can tile the circle: `7 · 1/7 = 1`.
So no per-runner union bound crosses `c = 7` — the union of seven `1/7`-arcs can be
the whole circle, leaving no lonely point. Every agent who has attacked the compressed
`≥ 7`-far leg of `lrc14_of_farcut_split` has hit this same `7 × 1/7 ≥ 1` wall.

The path-Bonferroni ("Hunter") correction says: the union is only the whole circle if
the arcs are *disjoint*; any overlap between consecutive arcs is a credit that reopens
a lonely gap. So the wall falls **iff the seven arcs are forced to overlap**. The whole
question becomes: *are they?*

The answer depends on the **speed scale of the block**, and it splits three ways.

## Regime A — clustered-huge (`w₁ ≥ 7392`): the arcs are frozen, and can't pack

If the seven far speeds are seven **distinct** integers within the drift budget
`1232·(w₇ − w₁) ≤ w₁`, then because they are distinct integers `w₇ − w₁ ≥ 6`, so
`w₁ ≥ 1232·6 = 7392`. Over a window of two base periods (`2/w₁` wide, tiny) the seven
runners' phases are **frozen** at nearly-constant offsets. Frozen seven `1/7`-arcs
cover the circle only by *perfect packing* — offsets exactly `{0, 1/7, …, 6/7}` — and
perfect packing pins the offset-sum to `3 ≡ 0`. Citing the **sum-combo**
`C₀ = Σ(wⱼ − w₁)` off the origin (one extra citation slot) contradicts that pin
quantitatively: some spacing exceeds `1/7` by `θ/22`, and its midpoint clears every
runner by `1/14 + θ/44`. `seven_gap_deficit` → `cluster_sweep_step` →
`cite_cluster7_lonely`. **Closed, unconditional, kernel-pure.**

The margin arithmetic is exact and beautiful: the bounded runners drift only
`B·2/w₁ ≤ 1/182` over the sweep window, so the `1/13` citation margin lands on
`1/13 − 1/182 = 13/182 = 1/14` — the sharp band, with nothing to spare and nothing
lacking.

## Regime B — spread (some consecutive `D·L ≥ 2`): the arcs are forced to overlap

If the block is *not* tightly clustered, some consecutive pair `(w₁, w₂)` has spacing
`D = w₂ − w₁` large enough that, across the window, the `w₂`-teeth sweep through a full
`w₁`-tooth: the residue `r ≡ −jD (mod w₂)` wraps, and each wrap deposits trapezoid
overlap. That overlap is the Hunter credit. klein's `LRCSpreadPairFloor`
(`clipLen_tooth_tooth` → `per_tooth_ge_trap` → `walk_one_wrap`) is exactly this: the
first pristine-window pair credit, `≥ L/49 − err`.

**The one-pair simplification (this session).** For `c = 7` the singles density is
`7 · (L/7) = L`, which cancels the window length *exactly*. So the Hunter ledger
`good ≥ L − singles + credits ≥ credits − F` needs only **one** pair credit to be
positive — the first, on the pristine window. The depletion-transported later pairs
(`pairs 2..6`) are not needed; they only add nonnegative mass (`pairCredits_ge_first`).
`cite_hunter_c7_onepair` reduces the entire spread-7 case to klein's single
first-pair floor. **Reduced to one named obligation, nearly discharged.**

This is worth pausing on: the general ledger asks for `(c−1)` pair credits totalling
`(c−1)·L/49`. At the exact wall `c = 7`, six of those seven credits are surplus. The
wall is thinnest precisely where we feared it thickest.

## Regime C — near-equal-small (`w₁ < 7392`, tight, e.g. 23…29): open

Seven far integers that are tight *but small* — consecutive-ish, like `23, …, 29` —
satisfy **neither** closer. `cite_cluster7` fails: the sweep window `2/w₁ ≈ 0.087` is
so wide that the offsets drift `~(w₇−w₁)·L ≈ 0.5`, destroying the frozen-arc picture
(the `1/616` deficit surplus is swamped). `cite_hunter_c7_onepair` fails: with
`D = 1` and `L ≈ 0.01`, the residue walk `D·L ≈ 0.01 < 2` never wraps inside the
window, so the first-pair floor can be zero.

This is the genuine near-equal regime the project has circled for many sessions — the
drifting pair-floor (`min pairmass/L ∈ [0.929, 1.000]·1/49`, exactly `1/49` at `D = 1`
but requiring the *right* window), mac-mini's `JointRateCore` per-cell obligation,
klein's `seven_commensuration`. It is the one-integer-wide slice of the wall that is
still standing.

## The meta-point: scale selects mechanism

The lesson is not "we need a cleverer single argument." It is that **the correlation
structure of seven danger arcs is a function of the block's speed scale**, and three
scales demand three mechanisms:

| regime | scale | why arcs overlap (or don't) | mechanism |
|---|---|---|---|
| A clustered-huge | `w₁ ≥ 7392` | frozen; can't perfectly pack | sum-combo gap-deficit (sweep) |
| B spread | some `D·L ≥ 2` | walk wraps → trapezoid credit | first pristine pair-floor |
| C near-equal-small | `w₁ < 7392`, tight | frozen-ish but window too coarse | **open** (drifting floor) |

klein's meta-observation dovetails: union-bound failure and pair-credit surplus are
the *same* phenomenon (danger-arc correlation) with opposite signs — the wall breaks
worst exactly in the tight regime that supplies the overlap that rescues it. Regime C
is where that duality is most delicate: the arcs *are* correlated (tight block), but
the window we are allowed (set by the citation margin `δ ~ 1/B`) is too coarse to
*see* the correlation as measure. Closing C is a statement about choosing the window
commensurate with the block — which is precisely what a *drifting* (window-adaptive)
pair-floor does.

Seven is not one wall. It is three, and two of them are down.

---
*Linked: [[everything-is-the-triangle]] (the 1/7 arc as the hypotenuse scale),
MISTAKE-071 (the vacuous integer-shift form that forced this split), HYP-3980/3981,
klein HYP-4022/4023, mac-mini HYP-3876.*
