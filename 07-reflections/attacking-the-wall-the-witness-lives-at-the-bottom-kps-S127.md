# Attacking the wall — the witness lives at the bottom of the window, not the transfer tail

*kind-pasteur-2026-07-10-S127. Owner: "now attack the stability statement." I did not break the wall —
it is the open case of LRC(14). But the attack changed its shape twice, gave one proved theorem, and
turned an asymptotic obligation into a bounded, small-modulus one. This records the reconnaissance.*

---

## The obligation, before the attack

`StrictlyLiveSupply`: every residual family admits a modulus `q` and multiplier `p` with all thirteen
residues strictly inside `(q/14, 13q/14)`. Klein's THM-685 delivers this at every `q > Σv/μ(S)` — a
threshold measured at **3000–5000** for typical residual families. That framing says: find `μ > 0`, then
the transfer hands you a ruler far out in the tail.

## What the data actually says

I computed the *minimal* strictly-live modulus directly. It is not in the tail. It is at the floor:

```
Vmax range   min q   max q   mean   failures
[1,45]         15      23    17.4      0
[1,1200]       15      27    18.3      0
```

`min q ∈ [15, 27]`, **scale-invariant** — the same window from `Vmax = 45` to `Vmax = 1200`, zero
families without a witness. The transfer threshold `Σv/μ ≈ 3000–5000` is two orders of magnitude larger
and never reached. The witness appears at the *bottom* of the window. The transfer is a safety net for a
ball that lands at your feet.

## Why `15` — and this one is proved

The lower edge is not empirical. `CoveringFamily v := ∀ q ∈ [2,14], ∃ i, q ∣ vᵢ`. If `q ∣ vᵢ` then
`(vᵢ·p) % q = 0` for every `p`, and `0` is never strictly above `q/14`. So a covering family has a **dead
residue at every modulus `2 ≤ q ≤ 14`** — every such `q` is killed outright, for all `p`. Hence:

> **`strictlyLive_modulus_ge_15` (Lean, kernel-pure): `CoveringFamily v → StrictlyLive v q p → 15 ≤ q`.**

The covering hypothesis does not merely dispatch small families to a census — it *calibrates the search
window*, closing exactly `[2,14]`. And empirically the witness appears immediately at the first surviving
modulus, `q = 15`.

## The adversary, and why the residual class already excludes it

`min q` *can* be pushed up: give a speed divisible by many moduli — say `vᵢ = lcm(15,…,B)` — and it kills
`q = 15,…,B` by a zero residue, forcing the witness above `B`. But a speed divisible by every `q ∈ [15,B]`
while the other twelve are not is exactly a **detuned** family (`∃ g ≥ 2` dividing all-but-one), and the
residual's `hdiv` hypothesis excludes precisely that. The adversary that beats the small-ruler law is a
family the assembly has already dispatched.

On the genuine residual class — covering, ratio `> 13`, compressed, **not detuned**, difference-primitive
— the measured maximum of `min q` over 200 families (`Vmax ≤ 1200`) is **27**. Not a proof of a bound, but
the adversary is pinned.

## The wall, resharpened

So the attack produced a sharper conjecture, formalized as the hypothesis of a kernel-pure reduction:

> **`BoundedStrictlyLiveSupply B`**: every residual family has a strictly-live ruler at some `q ≤ B`.
> With `strictlyLive_modulus_ge_15`, the witness is confined to `[15, B]`, and
> `lrc14_of_boundedStrictlyLiveSupply : LRCUpTo13 → BoundedStrictlyLiveSupply B → LRC14Statement`.

Measured `B = 27` suffices on everything sampled. This is genuine progress on the *shape* of the wall: the
remaining obligation is now **residue-level and bounded** — for each residual family, one modulus in a
window of size ~13 and one multiplier — not an asymptotic density estimate. `StrictlyLiveSupply` is its
`B = ∞` form, and remains open.

## Why I still cannot close it

Proving `BoundedStrictlyLiveSupply 27` for *all* residual families (infinitely many, unbounded `Vmax`) is
a uniform Diophantine statement: no covering, ratio-`>13`, non-detuned family, however large, evades a
live multiplier below 27. That is the Lonely Runner Conjecture for the hardest case, in its most concrete
dress. What the attack bought is the dress — small, integer, adversary-pinned — not the theorem.

Three honest limits on the reconnaissance:

- **The bound is measured, not proved.** `B = 27` held on ~500 families to `Vmax = 1200`; nothing rules
  out a slow logarithmic creep at `Vmax = 10⁶`. (klein's transfer *guarantees* the tail exists, so the
  worst case is some finite `B`; whether it is 27 or 270 is open.)
- **The lower edge is the only proved half.** `q ≥ 15` is a theorem; `q ≤ B` is not.
- **Bounded ≠ small enough to enumerate.** Even `B = 27` leaves a residue condition over an infinite
  family of `v`; it is not a finite check.

## What I would hand the next attacker

The obligation is now: *for each residual `v`, one of the ~13 surviving moduli `q ∈ [15,27]` admits a
multiplier `p < q` with all residues in the band.* That is a covering/pigeonhole statement over a **bounded
set of small moduli** — the natural target for the residue-class and Fourier-over-`ℤ/q` machinery the
fleet already runs (klein's character frame lives at exactly these `q`; the mcorr/hyperbola stack bounds
per-cell correlations at fixed `q`). The wall did not fall, but it is now a wall of height 13, standing on
a foundation I proved.

*Files: `LRCSmallRuler.lean` (`strictlyLive_modulus_ge_15`, `BoundedStrictlyLiveSupply`,
`lrc14_of_boundedStrictlyLiveSupply` — sorry-free, kernel-pure), `lrc14_small_ruler_law_kps_S127.py`/`.out`.
Builds on `LRCStrictRuler` (kps), klein THM-685, `LRC14CertRoute.CoveringFamily`. See
[[describing-the-wall-the-scale-gap-is-the-separator-kps-S127]].*
