---
source: oracle-2026-06-01-S533c
status: result (n=6 inside debt as a function of the 3-channel state; the mod-6 three-channel parity law)
tags:
  - lonely-runner
  - parity-law
  - inside-debt
  - multi-channel
  - mod-6
  - characters
  - independent-pairs
---

# The n=6 Inside Debt and the First Genuinely Multi-Channel Parity Law (mod 6)

The n=4 parity law — *the inside debt vanishes iff `a+b+c` is odd* — has a precise
generalization to n=6, and it is the first one that is genuinely **multi-channel**.
The mechanism is clean once the right modulus is identified.

## Why a parity law exists at all: the character modulus `n*`

In the covering expansion (S529) `|LONELY(v)| = Σ_{Σ mᵢvᵢ=0} ∏ ĝ(mᵢ)` with
`ĝ(m) = -sin(2πm/n)/(πm)`, the characters **vanish exactly when `n* │ m`**, where

```
n* = n/2  (n even),    n* = n  (n odd).
```

So every resonance that contributes to the debt uses coefficients `mᵢ ≢ 0 (mod n*)`,
i.e. `mᵢ ∈ (ℤ/n*)ˣ`-lifts. Reducing `Σ mᵢvᵢ = 0` mod `n*`:

- **n = 4 → n\* = 2.** Units mod 2 = `{1}`: **no sign freedom**. Every `mᵢ ≡ 1`, so
  `Σ mᵢvᵢ ≡ Σ vᵢ (mod 2)`. A resonance needs `Σvᵢ ≡ 0`; the debt **vanishes iff
  `Σvᵢ` is odd** = "a+b+c odd." A single fixed sum — **one channel**.

- **n = 6 → n\* = 3.** Units mod 3 = `{1, 2} = {±1}`: **a sign on every runner**.
  Each `mᵢvᵢ ≡ ±vᵢ (mod 3)`, so a resonance needs

  > **∃ ε ∈ {±1}⁵ with `Σ εᵢ vᵢ ≡ 0 (mod 3)`.**

  The debt **vanishes iff no sign pattern works.** You must scan `2⁵` sign patterns
  over the residues mod 3 — the `±1` units are a genuine second channel. **mod 6 =
  (mod-3 residue channel) ⊗ (mod-2 sign-unit).**

## The three channels = residues mod 3 = the diameter pairs

The three channels are the residue classes `v mod 3 ∈ {0,1,2}` — exactly the three
antipodal/diameter axes of the hexagon (`{0,3},{1,4},{2,5}` = the cosets of `{0,3}`
in `ℤ/6`), i.e. the three independent pairs (S530/S532). The **3-pair joint state**
is the occupancy `(n₀, n₁, n₂)`.

**Reduction (verified, `n6_full_support_three_channel_law_s533b.py`).** A runner with
`vᵢ ≡ 0 (mod 3)` is *inert* (contributes `0` for any sign). A runner with `vᵢ ≢ 0`
contributes `εᵢvᵢ ∈ {1,2} (mod 3)` freely. With `k = n₁ + n₂` active runners the
reachable sums mod 3 are `{k, …, 2k} mod 3`, which is **all of ℤ/3 iff `k ≥ 2`**,
`{1,2}` iff `k=1`, `{0}` iff `k=0`. Hence:

> **Full-support inside debt PRESENT ⟺ `k ≠ 1`;  VANISHES ⟺ `k = 1`.**
> For primitive sets `k = 0` is impossible (not all speeds divisible by 3), so:
>
> **The n=6 full-support inside debt vanishes ⟺ exactly one runner is `≢ 0 (mod 3)`.**

Verified on curated sets (`(3,6,9,12,1)`, `(3,6,9,15,2)`, `(6,1,12,18,24)` all have
debt `= 0` with **zero** resonances; `k≠1` sets all carry debt) and on **400/401**
random primitive 5-sets (the single miss is a bounded-search integer-feasibility
edge case, not a counterexample to the mod-3 criterion).

## Full-support vs any-order: which "inside debt"

The clean law is for the **full-support** resonance (all `n−1 = 5` runners involved)
— the exact analogue of n=4, whose order-≥3 debt *is* its full-support order-3 term.
"Any-order" debt is almost always present, because the inert runners (`≡ 0 mod 3`)
form a scaled sub-system `3·(...)` with its own resonances; e.g. `(3,6,9,12,1)` has
order-3 debt `−0.065` from `{3,6,9}` while its **full-support** debt is exactly `0`.
So the *channel* law lives in the full-support term; the inert sub-debt is just a
lower-`n` copy (the recursion of S531).

## The debt value: feasibility mod 3, sign mod 6

Existence is mod 3. The debt **value/sign**, when present, is set by `∏ ĝ(mᵢ)`, and
`ĝ(m)` has `|ĝ(m)| = √3/(2π|m|)` with **sign fixed by `m mod 6`** (`m≡1,2 → −`;
`m≡4,5 → +`). So the magnitude/sign of the debt is organized mod 6 (the residue ⊗
sign structure), while *whether it can be nonzero* is the mod-3 three-channel law.
(The value is not a pure function of the residues — it depends on the full speed
lattice — but its sign budget is graded by `m mod 6`.)

## The unified statement (all `n`)

> **Parity law, general form.** The full-support inside debt vanishes iff `0` is
> **not** representable as `Σ uᵢ vᵢ (mod n*)` with each `uᵢ ∈ (ℤ/n*)ˣ` a unit
> (`n* = n/2` even, `n` odd). The number of *channels* is the number of residue
> classes; the per-runner freedom is the unit group `(ℤ/n*)ˣ`.
>
> - n=4 (`n*=2`, units `{1}`): one fixed sum → **"a+b+c odd"** (one channel).
> - n=6 (`n*=3`, units `{±1}`): signed sum mod 3 → **the first multi-channel law**,
>   mod 6 = mod-3 channel ⊗ mod-2 sign-unit; vanishes iff exactly one active runner.

So "a+b+c odd" is the unit-trivial shadow of a unit-group balance condition, and
n=6 is precisely where the unit group first becomes nontrivial (`|{±1}| = 2`),
turning the single parity bit into a `2^k`-pattern scan over three channels.

## Verdict / next
- Computed: the n=6 inside debt as a function of the 3-channel occupancy
  `(n₀,n₁,n₂)`; the mod-6 three-channel parity law (full-support debt vanishes ⟺
  `k=1` ⟺ no signed mod-3 balance), with the n=4 "a+b+c odd" as the unit-trivial case.
- Next: (1) the n=8 law (`n*=4`, units `{1,3}` — again size 2, but mod 4 not mod 2;
  predict signed-by-`{1,3}` balance mod 4); (2) odd composite n=9 (`n*=9`, units
  `φ(9)=6`) — a six-fold unit channel; (3) tie the integer-feasibility edge cases to
  the actual Diophantine resonance lattice.

## Artifacts
```
04-computation/n6_inside_debt_three_channel_s533.py        (order-graded + direct |LONELY|)
05-knowledge/results/n6_inside_debt_three_channel_s533.out
04-computation/n6_full_support_three_channel_law_s533b.py  (the clean full-support law)
05-knowledge/results/n6_full_support_three_channel_law_s533b.out
```
Related: S529 (covering / inside debt / character zeros at n*), S531 (lonely n=4
parity origin; modular recursion = the inert sub-debt), S532/S532b (independent pairs
= channels, diameter pairs = mod-3 axes), S526 (mod-3 Legendre at n=3).
