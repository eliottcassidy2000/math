# The residual is one measure floor — the LRC(14) formalization's final shape

*kind-pasteur-2026-07-10-S127. Owner: "close out the residual obligation; get the LRC(14) proof
formalization into its best possible state." I could not close the residual — it is the open analytic
core. What I could do is make the formalization say, exactly and kernel-purely, what remains. This note
records the resulting shape.*

---

## The statement, in one line

```
lrc14_of_measureFloor : LRCUpTo13 → SafeMeasureFloor → LRC14.LRC14Statement
    axioms: [propext, Classical.choice, Quot.sound]
```

**LRC(14) follows, with foundational axioms only, from the LRC(≤13) citation plus a single measure
floor.** Nothing else remains: every branch of the grand assembly is discharged, and the last
`native_decide` (the `C(22,13) = 497 420` window census) was removed earlier today.

## Why the residual *is* a measure floor, and nothing more

`Lonely 14 v t` unfolds to `∀ i m, 1/14 ≤ |vᵢ·t − m|`. That is not a statement *about* the safe set —
it *is* membership in it. So the grand assembly's `ResidualObligation` says precisely: *every residual
family has a nonempty safe set.*

That framing is tautological and therefore useless — until klein's **THM-685** (today). Its transfer
theorem `|LM(q) − q·μ(S)| ≤ K(S) ≤ Σvₗ` is elementary (band-rounding identity + a crossing count +
per-interval sampling; `q` arbitrary, primality nowhere), and it closes the **modulus side** by counting
alone. What it leaves is klein's own verdict: *the remaining analytic content of the covering case is
exactly the measure floors.*

So the useful reduction is not to *nonemptiness* (a point — which is what you are trying to prove) but to
**positive measure** (a quantity — which is what the analytic machinery actually bounds). That asymmetry
is the whole content of this file:

- `lonely_of_safePeriod_measure_pos` : `0 < volume (safePeriod v) → ∃ t, Lonely 14 v t`.
  A positive-measure set is nonempty, and any of its points is a witness.
- `residualObligation_of_measureFloor` : `SafeMeasureFloor → ResidualObligation`.
- `lrc14_of_measureFloor` : the citation plus the floor gives LRC(14).

The floor is *strictly stronger* than the obligation — and that is exactly why it is the right hypothesis.
You cannot exhibit the lonely point; you can bound the measure of the set of them.

## Two forms, because the machinery produces intervals

The analytic stack does not hand you a measure; it hands you an interval. mac-mini's
`LRCIntervalBridge.Ico_subset_safeSet_of_bounds` produces `Ico a b ⊆ safeSet`; boxeph's `μ_L` evaluator
produces rational interval certificates. So the file carries both entry points:

```
SafeIntervalFloor  →  SafeMeasureFloor  →  ResidualObligation  →  LRC14Statement
   (interval)          (measure)            (a point)             (the theorem)
```

with `safePeriod_measure_pos_of_Icc_subset` the only real step (`volume (Icc a b) = ofReal (b−a) > 0`,
then monotonicity). `lrc14_of_intervalFloor` is the version a per-family interval certificate discharges
directly.

## What the day's arithmetic actually was

The endgame did not close by one grand argument. It closed by *subtraction* — each branch peeled until
one quantity was left standing:

| branch | discharged by |
|---|---|
| non-covering | the `1/q` sieve |
| ratio ≤ 13 | `spread13_lonely` |
| window ≤ 22 | LEM-024 six-witness pigeonhole (opus-S202, on the kps-S127 witnesses) — **native_decide removed today** |
| dominant / repeated / detuned | peels and dispatches |
| common-residue | THM-682(a), `M ≥ 8/17` |
| modulus side of the covering case | THM-685 transfer (klein, today) |
| **everything else** | **`SafeMeasureFloor`** |

The pigeonhole and the transfer landed within hours of each other, and together they collapsed the two
things that had looked like walls — a half-million-family census and an infinite family of rulers — into
counting arguments. What survived is not a wall. It is a number: `μ(S) > 0`.

## The honest ledger

- **Proved, kernel-pure:** the entire assembly, top theorem included. `[propext, Classical.choice,
  Quot.sound]`, no `native_decide`, no `sorry`.
- **Cited (owner directive):** LRC(≤13).
- **Open:** `SafeMeasureFloor` — for every covering, scale-gapped, compressed, distinct,
  difference-primitive, non-near-AP, no-detune family reaching past the window, the safe set has positive
  measure. This is the density/witness floor: mac-mini is one measure-theory brick from it on the
  interval route; klein's Corollary 1 turns any floor into explicit rational certificates at every
  `q > Σv/μ₀`, with the small-`q` bank a priori finite.

I did not close the residual. I made the formalization state it — once, precisely, and with nothing else
attached. That is the best available state, and the next agent who proves a measure floor for the
residual class finishes LRC(14) by supplying one hypothesis.

*Files: `LRCResidualMeasureFloor.lean` (sorry-free, kernel-pure), `LRC14GrandAssembly.lean` (the swap).
Builds on: klein THM-685 (transfer), opus LEM-024 + kps-S127 (the six witnesses and the swap),
mac-mini `LRCIntervalBridge`, boxeph `μ_L`. See
[[the-final-rung-lives-on-the-diagonal-dyadic-carrier-free-triangles-kps-S127]].*
