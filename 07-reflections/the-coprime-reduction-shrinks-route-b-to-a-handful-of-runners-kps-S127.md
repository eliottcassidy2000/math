# The coprime reduction: Route B's anti-concentration is on a handful of runners, not thirteen

*kind-pasteur-2026-07-11-S127 cont.45. Owner: "work the sharpest remaining math, pulling often." The sharpest
remaining math (klein-S258 finish map) is Route B — spread divisor-complete families clear at some non-14
`q ∈ [15,29]`. opus-S238 found no small-prime shortcut because at a prime the whole family is coprime. This
session shows that at a **composite** modulus the effective family collapses to a handful, and that collapse
is exactly what divisor-completeness forces.*

---

## The lemma

At modulus `q` with danger arc `{0, ±1}` (safe band `[2, q-2]`, the regime `q ∈ [15,28]`), take a **unit**
multiplier `p` (`gcd(p,q)=1`). For a runner `v` write `g = gcd(v,q)`:

- `g = 1` (coprime): `v·p` is a unit, in danger iff `p ≡ ±v⁻¹` — 2 bad multipliers;
- `1 < g < q` (proper common factor): **always safe**. `v·p ≡ 0 (mod g)` so it cannot be `±1`; and
  `v·p ≡ 0 (mod q)` would need `(q/g) ∣ p`, impossible for a unit `p`. Concretely `r = v·p mod q` is a
  *nonzero multiple of g*, so `g ≤ r ≤ q - g`, hence `2 ≤ r ≤ q - 2`;
- `g = q` (`q ∣ v`): stuck at `0` for every `p` — always in danger.

So **at a unit multiplier the only runners that can be unsafe are those coprime to `q`** (plus any runner
divisible by `q`, which must simply be absent). Clearing via a unit multiplier is therefore:

> `q ∤ vᵢ` for all `i`, **and** the coprime-to-`q` runners miss a unit fold-class.

And a clean **sufficient** condition: if `q ∤ vᵢ` for all `i` and `#{i : gcd(vᵢ,q)=1} ≤ φ(q)/2 − 1`, then the
coprime runners occupy fewer than the `φ(q)/2` fold-classes, so one is empty — pick a unit `p` there and the
family clears.

## Why this is the right handle — and why the prime analysis stalled

opus-S238 analysed the *prime* moduli `17, 19, 23` and found no small shortcut: at a prime, *every* speed is
coprime, so the anti-concentration is on all 13 runners, and an adversary can occupy every fold-class. The
composite moduli are the opposite. At `q = 18 = 2·3²`, `24 = 2³·3`, `20 = 2²·5`, `26 = 2·13`, a runner
divisible by any prime factor of `q` is automatically safe — and **divisor-completeness forces exactly those
runners to exist** (multiples of `2, 3, …, 14`). So the coprime-to-`q` sub-family is small: over 3500 random
spread divisor-complete families, the minimum over the window of `#coprime-to-q` has median **3** — the
13-runner problem shrinks to a 3-runner problem. The very structure that makes divisor-complete families the
hard core (they contain multiples of everything) is what empties the effective family at composite moduli.

**Verified:** the elementary sufficient condition (`q ∤ vᵢ` and `#coprime ≤ φ(q)/2−1` at *some* window `q`)
provably clears **91.5%** (3201/3500) of spread divisor-complete families. The remaining 8.5% all clear too —
but by the genuine fold-class miss (coprime runners still miss a class even when there are `≥ φ(q)/2` of
them), now a covering question on `≤ 8` runners instead of 13.

This is the quantitative mechanism behind opus-S239's "spread = bad coverer": a bad coverer leaves many
multipliers uncovered *because* most of its runners are structurally safe at composite moduli, so few runners
are left to do the covering.

## The error the check caught

My first statement of the lemma said "every non-coprime runner is safe at a unit `p`" — and the machine
check flagged 40154 violations. They were all the `g = q` case: a runner *divisible by* `q` sits at `0` for
every multiplier, unit or not. The corrected lemma excludes multiples of `q` (which is why the guarantee
carries the `q ∤ vᵢ` clause). A one-line computational assertion turned a false lemma into a true one before
it reached the write-up — the same discipline that has repeatedly saved this project (exhaust, don't assume).

## What is formalized

`LRCThreeGapConsecutive.lean` (kernel-pure `[propext, Classical.choice, Quot.sound]`):
`inBand_of_proper_common_factor` — a runner with a proper common factor `g` (`2 ≤ g`, `g ∣ vᵢ`, `g ∣ q`) and
`q ∤ vᵢ·p` is `inBand` (safe), via `g ≤ r ≤ q - g`. This is the machine-checked core of the coprime
reduction; the fold-class-miss clearing step (pigeonhole on the coprime sub-family) is the remaining Lean
work, now on a small family.

## Honest scope

This does not close Route B — opus-S238's "no small shortcut" stands: the 91.5% is a union over the window,
no single `q` suffices. What it adds is a **reduction of the effective problem size**: the spread-bulk
anti-concentration is really a fold-class covering statement on the *coprime-to-q* sub-family (median 3, at
most ~8), not on all 13 runners. The residual 8.5% is the genuine remaining core, and it is now small enough
that a per-`q` finite fold-class argument — or the LEM-010 diameter-bounded census — is within reach. The
right next move is to attack the residual as a bounded covering problem on the coprime sub-family.

## Complement to klein's exact prime formula (THM-718, same day)

klein-S259 derived the exact clearing-count at a **prime** `q` with `q ∤ vᵢ`:
`(q−1) − |{±j·vᵢ mod q : 1 ≤ j ≤ m}|`, `m = ⌈q/14⌉−1` — a fold-class covering number on the *full* speed set
(at a prime every speed is coprime). The coprime reduction is the **composite-modulus** companion: at composite
`q` the same covering number is on the *coprime-to-q sub-family only*, and unit multipliers replace all
multipliers. klein's prime formula quantifies why no single prime shortcut exists (full family, adversary
covers); the composite reduction shows why the *union* over composite moduli is favorable (shrunken family,
few coverers). Two moduli regimes, one covering-number picture.

*Files: lrc14_coprime_reduction_kps_S127.py/out, LRCThreeGapConsecutive.lean (`inBand_of_proper_common_factor`).
HYP-6100. Sharpens opus-S238 (primes → full family) and opus-S239 (bad-coverer mechanism); the composite
companion to klein-S259/THM-718 (exact prime clearing-count); feeds klein-S258 Route B, opus-S235 band-edge.
Extends [[the-three-gap-disjunction-collapses-to-consecutive-for-the-ap-corner-kps-S127]].*
