# The mod-2g lift closes only its base case — the rest is a 2-adic descent

*kind-pasteur-2026-07-10-S127. Owner: "prove the (2,2) mod-2g lift, the last detuned residual." I proved the
part that closes and report honestly that the rest does not — it descends to the open core. This note is
the boundary.*

---

## The residual, and why the branch count is powerless there

The last detuned residual is the congruent half-harmonic pair: `v = g·H ∪ {δ₁, δ₂}` with `q₁ = q₂ = 2`, so
`δ₁ ≡ δ₂ ≡ g/2 (mod g)`. A `q = 2` coordinate contributes `badCount = g/2` — exactly half the branch
interval `[0,g)` — so two of them tile `[0,g)` completely and no single branch clears both. The count is not
merely inconvenient here; it is provably tight (`badCount₁ + badCount₂ = g`), which is why this case sits
outside the generic dispatch.

## The lift is a descent, not a door

monad's THM-678 names the escape: **lift to modulus `2g`.** There the two detuned become `q = 4`
(`gcd(δᵢ, 2g) = g/2`), and `¼ + ¼ < 1` — the generic count would fire. But doubling the modulus has a price,
and the price is the whole difficulty:

> every **odd** harmonic multiplier `m` becomes a *new* half-harmonic of `2g` — `g·m ≡ g (mod 2g)`, i.e.
> `q = 2` at `2g`.

So the lift trades the two original half-harmonics for however many odd-multiplier ones. I verified the
arithmetic exactly (`lrc14_two_detuned_lift`): with `k` odd multipliers, the count at `2g` is
`½ + k·½` — `< 1` only when `k = 0`. For `k ≥ 1` the lift fails at `2g` too, and one must lift again to `4g`,
where the same thing recurs. The obstruction **descends the 2-adic tower**; it does not terminate at a
citation. Its true home is the open core — mac-mini's THM-676 descent-burden (the `≥ 11` half-sum moduli,
Freiman-rigid at APs) and klein's pair-sum / dead-zone route, both of which attack it from the measure side.

**I did not fake a closure.** The general `(2,2)` residual is not independently dischargeable from
LRC(≤13), and the canon says as much.

## What the lift does close — the terminating base case

The descent halts in one step exactly when `k = 0`: **no odd harmonic multiplier**, i.e. the entire
harmonic part is divisible by `2g`, not merely `g`. Then at `2g` the only non-multiples are the two detuned
speeds, both `q = 4`, and opus's *already-proved* generic `d = 2` dispatch fires. That is
`lonely14_of_two_detuned_lift2g`, kernel-pure:

```lean
theorem lonely14_of_two_detuned_lift2g (cite) (v) (hv) (g) (hg : 2 ≤ g) (i₁ i₂) (hne)
    (hlift : ∀ j ≠ i₁, i₂, 2*g ∣ v j)  (hδ1 : ¬ g ∣ v i₁) (hδ2 : ¬ g ∣ v i₂) :
    ∃ t, Lonely 14 v t
  := lonely14_of_two_detuned' cite v hv (2*g) … i₁ i₂ …
```

The pretty part is how little it needs. Opus's dispatch at `2g` wants `¬ 2g ∣ δᵢ` and `(q₁,q₂) ≠ (2,2)` at
`2g` — and **both follow from `¬ g ∣ δᵢ` alone**, no gcd-valuation arithmetic: `2g ∣ δ ⟹ g ∣ δ`, and
`2g / gcd(δ,2g) = 2 ⟹ gcd(δ,2g) = g ⟹ g ∣ δ`. So the only real hypothesis is the `k = 0` divisibility of
the harmonics. The lift is opus's theorem viewed at the doubled scale; the content is recognizing *when* the
view is still a `d = 2` picture.

## The boundary, stated cleanly

The `(2,2)` residual splits along `k = #odd multipliers`:

- **`k = 0` (base case):** CLOSED — `lonely14_of_two_detuned_lift2g`, kernel-pure.
- **`k ≥ 1` (descent):** the open-core 2-adic tower — cited, not closed; the fleet attacks it from the
  measure side (THM-676, klein's dead-zone).

This is the honest shape of a residual that is *genuinely* residual. Two sessions of counting closed the
generic bulk of the detuned dispatch to this pair; the pair closes to its base case; and the base case's
complement is not a smaller detuned problem but the same open core the whole endgame reduces to. The
difficulty was never eliminated — it was driven, step by step, to exactly where the measure floor already
lives. That convergence is the signal: the detuned peel and the measure floor are two faces of one wall.

*Files: `LRCTwoDetunedLift.lean` (`lonely14_of_two_detuned_lift2g`), `lrc14_two_detuned_lift_kps_S127.py`/`.out`
(the `k`-sweep verifying the split). Builds on opus's `LRCTwoDetunedClearing`. Closes the base case of the
residual left open by [[the-detuned-citation-shrinks-to-its-half-harmonic-residual-kps-S127]]; the descent
is monad's THM-678, the open core.*
