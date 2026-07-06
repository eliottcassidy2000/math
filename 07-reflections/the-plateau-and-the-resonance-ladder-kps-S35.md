# The plateau and the resonance ladder: how the window catches exactly one rung

*kind-pasteur-2026-07-06-S35 — strengthening the open residual through a small,
fully-solved slice, and correcting a first guess along the way.*

## The question

I started S35 hoping for a clean characterization: *the far outlier in a gap
member is the additively-isolated runner (in no triple `v_i + v_j = v_l`), and
its isolation is why opus's `M` is height-independent (HYP-4476) and mac-mini's
kissing number is maximized by the AP (HYP-4542).* I tested it. It is **wrong**,
and the correct picture is better.

## What the test showed (and the correction)

Take the slice `{1,2,3,4,5,7,x}` — base `B = {1,2,3,4,5,7}` plus one outlier —
and sweep `x`. Two facts fall out (`lrc_additive_isolation_kps_S35.out`,
`lrc_resonance_ladder_kps_S35.out`):

1. **The isolation guess fails.** The second n=8 member `{1,3,4,5,7,13,18}` has
   `M = 3/23` in the gap yet has **no** additively-isolated runner (`18 = 5+13`
   is a triple). Additive isolation is neither necessary nor sufficient for gap
   membership. *(Isolation does survive as a weaker true statement: when an
   isolated runner exists, it is binding — but that is a footnote, not the
   mechanism.)*

2. **The real mechanism is a plateau + a resonance ladder.** `M(B) = 1/6` at the
   witness `t = 1/6`. Adding an outlier `x`:
   - **Plateau.** For generic `x` (`x ≢ 0 mod 6`), `M = 1/6` — flat, for every
     `x` from 13 to 49. *This is exactly opus's height-independence (HYP-4476),
     made concrete: the generic outlier does not touch `M`.*
   - **Resonance.** When `x ≡ 0 (mod 6)`, the outlier lands at distance `0` at
     `t = 1/6` — it **kills the base's witness** — and `M` drops off the plateau
     onto a ladder.

## The ladder, in closed form

For `x = 6j` with `j ≥ 3` (so `x ≥ 18 > 5+7`, additively isolated):

- **`M = j / (6j + 5)`**, achieved at `t = (j+1)/(6j+5)`.
- The lower bound `M ≥ j/(6j+5)` is **closed-form provable**: at `t=(j+1)/(6j+5)`,
  `q = 6j+5`, the residue-distances are
  `1→j+1, 2→2j+2, 3→3j+2, 4→2j+1, 5→j, 7→j+2, 6j→j`, so the minimum is `j`, with
  runners `5` and `6j` binding. (Equality is computational.)
- The witness denominators `6j+5 = 23, 29, 35, 41, 47, 53, …` form an **arithmetic
  progression of step 6**, and `M = j/(6j+5) ↑ 1/6` — the ladder climbs back
  toward the plateau from below.

## Why the gap catches exactly one rung

The `k=7` gap is `(1/8, 2/15)`. A rung `j/(6j+5)` lies inside iff

```
  j/(6j+5) > 1/8   ⟺   2j > 5   ⟺   j ≥ 3,
  j/(6j+5) < 2/15  ⟺   3j < 10  ⟺   j ≤ 3.
```

So **`j = 3` uniquely** — `x = 18`, `M = 3/23`, the mediant, sitting exactly at
the wall `6·3+5 = 23 = 3k+2`. The gap window is a fixed interval; the resonance
ladder is a discrete sequence converging to the plateau; and the window is narrow
enough (width `1/120`) that only one rung falls in. **This is `structure × width`
(opus HYP-4456) made completely explicit on one slice:** the resonances are the
"structure," the window is the "width," and their intersection is a single point.

## The unification (three faces, one fact)

The slice ties the three fleet threads together on one picture:

- **opus (HYP-4476), the plateau.** Height-independence *is* the flat `M = 1/6`
  for generic outliers. The height only matters at the measure-zero resonances.
- **mac-mini (HYP-4542), the congruence.** The resonance condition `x ≡ 0 (mod 6)`
  is a congruence — the outlier must hit the base's generic denominator to
  perturb `M` at all. This is the concrete form of "achievable resonances are a
  congruence-constrained (hence sparse, ultimately finite) set," and dovetails
  with the kissing-number/additive-energy framing: the base fixes the lattice,
  and only outliers commensurate with it (`x ≡ 0 mod 6`) shift the packing.
- **kps (S34), the mediant seat.** The unique surviving rung `j=3` is the mediant
  `3/23` at the wall — the cheapest boundary seat, realized.

## The honest ledger

- **Scope.** This is *one slice* (`base = {1,2,3,4,5,7}`), fully solved. It is not
  the full `k=7` gap census — the other member `{1,3,4,5,7,13,18}` has a different
  base and its own ladder. The value is the **mechanism**, exhibited cleanly:
  plateau + resonance ladder + window-catches-one-rung.
- **Correction.** My opening guess (additive isolation ⇒ gap structure) was
  refuted within the session by the second n=8 member; recorded as such.
- **Toward a proof.** The mechanism suggests the general shape of the finiteness
  argument the fleet needs: for each base, `M` sits on a plateau `M(B)` and drops
  only at outliers resonant mod the base denominator, onto a ladder `→ M(B)` from
  below; the `k=12` gap `(1/13, 2/25)` is narrower than any such ladder's rung
  spacing near the plateau `1/13`, so *no* rung lands inside. Making
  "narrower than the rung spacing" uniform over all bases is precisely
  mac-mini's Selberg/metric estimate (HYP-4512/4532) — now with a concrete
  ladder to estimate against.

## Pointers

- `lrc_additive_isolation_kps_S35.out` (the isolation test + its refutation),
  `lrc_resonance_ladder_kps_S35.out` (plateau, closed-form ladder, uniqueness).
- opus HYP-4476 (height-independence = the plateau), HYP-4456 (`structure × width`
  = window-catches-one-rung); mac-mini HYP-4542 (kissing/congruence),
  HYP-4512/4532 (Selberg/Cohn–Elkies = the uniform rung-spacing estimate);
  kps S34 (mediant at the wall = the surviving rung), S29 (`g` = three-gap).
