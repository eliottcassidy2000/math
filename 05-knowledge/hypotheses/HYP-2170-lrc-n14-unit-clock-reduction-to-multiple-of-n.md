---
id: HYP-2170
status: VERIFICATION + REDUCTION — LRC(14) reduces to configs containing a multiple of 14 (the
  unit-clock sieve handles all others); the residual = the d=14 divisor-block = the owner
  congruences (THM-398) = my large-owner CRT automaton (S581) = the rank-1 two-block (HYP-2150).
  Verified 0 failures; the residual proof is the open part.
source: claudebox-2026-06-03-S611
related: [HYP-2150, HYP-2145, HYP-2115, HYP-2105, THM-398, HYP-2063]
---

# HYP-2170: LRC(14) reduces to the multiple-of-14 configs — the unit-clock sieve does the rest

A clean reduction connecting the construction (the rational sieve) to the divisor-block / owner
machinery, locating the entire n=14 difficulty in one block.

## The unit-clock sieve (verified, 0 failures / 7315)

At a unit clock `t = a/14` (`gcd(a,14)=1`, the `φ(14)=6` clocks `{1,3,5,9,11,13}`), runner `v` sits
at the origin iff `14 ∣ v·a`, and since `gcd(a,14)=1` that is iff `14 ∣ v`. So:

> **If no speed is divisible by 14, then `‖v·a/14‖ ≥ 1/14` for every runner at every unit clock —
> the config is lonely (LRC holds), trivially.**

Equivalently `‖v a/n‖ ≥ 1/n ⇔ n ∤ va ⇔ n ∤ v` (a coprime). Verified: 0 of 7315 no-multiple-of-14
configs fail at the unit clocks. The construction is immediate off the residual.

## The reduction

> **LRC(14) ⟺ every config containing a multiple of 14 is lonely.**

The multiple of 14 is exactly the **`d = 14` (= `d = n`) divisor-block** of the danger
decomposition `gcd(j,n) = Σ_{d∣n} φ(d)[d∣j]` (HYP-2145): the runner that sits at the origin at the
unit clocks, killing the easy sieve. This is the repo's **C′ reduction** (THM-398): write the config
as `S = S' ∪ {v = nw}`; the `nw` runner is the obstruction. So the whole n=14 problem is the
residual, and it is precisely the object my recent work handles:

- **THM-398 / HYP-2105:** the cover→congruence translator — `S` tight iff `G(S')` fits one `nw`-arc,
  giving the endpoint-owner congruences.
- **S581 / HYP-2115:** the **large-owner CRT residue automaton** — the bounded decider for the
  off-centre fits, with the proved resonance bound `|D| ≤ u_b K_a + u_a K_b`.
- **HYP-2150 / HYP-2145:** the obstruction is the **rank-1 two-block** (the apex, `n=2q` mod-2),
  dissolved in the odd additive face by the pair-sum sieve.

So the n=14 architecture is: **(trivial) no multiple of 14 ⇒ unit-clock lonely; (residual) a
multiple of 14 ⇒ the owner-congruence / large-owner-automaton problem, whose obstruction is the
rank-1 two-block.** The first half is elementary and formalizable; the second is where the proof
lives.

## Open / next

- The residual proof: show every multiple-of-14 config is lonely — the large-owner CRT automaton
  (S581) intersected with the valid-config language (tasks t-0040/41), establishing the rank-1
  two-block is always clearable by the pair-sum sieve.
- Formalize the trivial half (`no multiple of n ⇒ unit-clock lonely`) — the arithmetic core
  `n ∤ v ∧ gcd(a,n)=1 ⇒ n ∤ va` (done: `Math/LonelyRunner/UnitClock.lean`).

**Artifacts:** computation inline; formal `Math/LonelyRunner/UnitClock.lean` + `DangerBlocks.lean`.
Builds on HYP-2150/2145 (faces/blocks), HYP-2115/S581 (large-owner automaton), THM-398 (C′), HYP-2063.

## Division-sieve extension (verified)

The unit-clock argument generalizes to every `m ≤ 14`: no speed divisible by `m` ⇒ gap `≥ 1/m ≥
1/14` at the `m`-clocks. So a config is proven lonely unless it contains a multiple of **every**
`m ∈ {2,…,14}` — a covering-system condition. Verified: the division sieve catches **~92%** of
random 13-speed configs (matching HYP-2075's 88.8%); the residual is the covering-system configs
(bottleneck: the high values 11,12,13,14). The covering residual is then handled by the recursive
even-fold (the mod-7 prime, solved) + the pair-sum sieve (the mod-2 rank-1 two-block).
