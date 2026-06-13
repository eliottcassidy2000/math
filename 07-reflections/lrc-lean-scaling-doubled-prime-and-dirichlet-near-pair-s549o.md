---
source: oracle-2026-06-01-S549o
status: formalization (Lean): scaling/repeated-addition, doubled-prime sieve, Dirichlet near-pair pigeonhole
tags:
  - lonely-runner
  - lean
  - formalization
  - repeated-addition
  - doubled-prime
  - pigeonhole
  - dirichlet
---

# Formalizing the Recent LRC Ideas in Lean: Scaling, the Doubled-Prime Sieve, and the Dirichlet Near-Pair

Three of the recent conceptual threads are now machine-checked, extending
`LonelyRunner.lean` (THM-369). Each `#print axioms` reports only the Mathlib
foundations `[propext, Classical.choice, Quot.sound]` — no `sorry`, no project axiom.

## `lonely_scale` — scaling invariance = the repeated-addition reduction (S548)

```lean
theorem lonely_scale (n : ℕ) (v : ι → ℤ) (t : ℝ) (c : ℤ) (hc : c ≠ 0) :
    Lonely n (fun i => c * v i) (t / (c : ℝ)) ↔ Lonely n v t
```

The machine-checked form of "multiplication = repeated addition" at the heart of LRC:
since `(c·v_i)·(t/c) = v_i·t`, loneliness is invariant under scaling the speeds by a
nonzero `c` and the time by `1/c`. So an arithmetic-progression family `c·(v_i)`
reduces to `(v_i)` at the scaled time `c·t` — the runner position `v_i·t` is `t`
*repeated-added* `v_i` times. This is the formal version of the S548 hyperoperation
reduction (the AP/repeated-addition orbit is the level-2 object) and the gcd/scaling
normalization used throughout the small-`n` proofs.

## `lonely_doubled` — the doubled-prime / `n*=n/2` sieve (S546)

```lean
theorem lonely_doubled (p : ℕ) (hp : 0 < p) (v : ι → ℤ)
    (hdiv : ∀ i, ¬ ((p : ℤ) ∣ v i)) : Lonely (2 * p) v ((1 : ℝ) / p)
```

For the *doubled* dimension `n = 2p`, if no speed is divisible by `p` then `t = 1/p`
(clearance `2/n = 1/p`) is lonely — a one-line corollary of the master sieve at
`q = p ≤ n`. This is the Lean witness of the S546 finding that doubled-prime
dimensions `n = 2p` inherit the clean `mod p` channel (`n* = p`); e.g. `n = 14 = 2·7`
is lonely at `t = 1/7` whenever no speed is a multiple of `7`.

## `near_pair` — the Dirichlet box pigeonhole ("always a near pair", S536/S539)

```lean
theorem near_pair (n : ℕ) (hn : 0 < n) (x : Fin (n + 1) → ℝ) :
    ∃ i j : Fin (n + 1), i ≠ j ∧ |Int.fract (x i) - Int.fract (x j)| < 1 / n
```

The genuinely new ingredient: among any `n + 1` reals, two have fractional parts
within `1/n`. Proof = Dirichlet's box principle — sort the points into the `n` boxes
`⌊n·fract⌋ ∈ {0,…,n-1}` (via `Int.floor_nonneg`, `Int.floor_lt`), apply the Finset
pigeonhole `Finset.exists_ne_map_eq_of_card_lt_of_maps_to` (`n+1 > n`), and convert
equal boxes to `|fract i − fract j| < 1/n` via
`Int.abs_sub_lt_one_of_floor_eq_floor`. In the runner picture (observer + `n-1`
runners = `n` points; the extra point gives the `n+1`) this is the machine-checked
form of two repeatedly-encountered facts: **the half-turn relation always carries a
tie** (the LRC trienerment is never a pure tournament, S539) and **some gap is `≤ 1/n`
(the apex pigeonhole, S530/S536)**. It is the first pigeonhole result in the LRC Lean
development.

### A formalization note (recorded so it isn't re-hit)

The first build left `near_pair` depending on `sorryAx`: an `omega` could not see
through a `set box := fun i => …` binding (the hypothesis `box i = box j` stayed
un-beta-reduced, so `omega` saw an opaque atom). The fix was to **drop `set` and pass
the explicit `(⌊n·fract (x i)⌋).toNat` to the pigeonhole**, so the returned equality
is *directly* the `toNat` equality and `omega` (with the two `0 ≤ ⌊·⌋` facts) closes
the floor-equality goal. Lemma: prefer explicit functions over `set` when a later
tactic must reason through the definition.

## Verdict / next
- Three new axiom-clean lemmas in `LonelyRunner.lean`: `lonely_scale` (repeated-addition
  reduction, S548), `lonely_doubled` (doubled-prime sieve, S546), `near_pair`
  (Dirichlet near-pair pigeonhole, S536/S539).
- Concrete next Lean targets: (1) the apex/largest-gap statement (some gap `≥ 1/n`)
  as the dual of `near_pair`; (2) the AP-orbit at `t = a/n` as the tight witness via
  `lonely_scale` + `initial_segment_unit_lonely`; (3) `Dirichlet`-style "some multiple
  of `t` is within `1/(N+1)` of `0`" (the anti-loneliness obstruction).

## Artifacts
```
04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean
  (Cases section: lonely_scale, lonely_doubled, near_pair)
```
Related: THM-369 (the Lean sieve), S548 (repeated addition / hyperoperation),
S546 (doubled primes), S536/S539 (pigeonhole / trienerment ties), S530 (apex).
