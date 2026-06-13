---
source: oracle-2026-06-01-S515
status: exploratory (wacky proof attempts at LRC n=18; honest partial results)
tags:
  - lonely-runner
  - n18
  - prime-square
  - triadic-ladder
  - bounded-ansatz
  - parity-split
  - wacky
---

# Wacky Attempts at LRC for n = 18

Same whimsy as the n=14 pass (S514), now at `n = 18 = 2·3²` — which is *richer*
than `14 = 2·7`: it has a prime **square** `9 = 3²`, and since `16 ≤ 18` its sieve
even forces the n=16 "sedenion" dyadic gate. Screened in
`lrc_n18_wacky_attempts_s515.py`.

## A — anti-concentration (same wall)

`N(t) = #{i : ‖v_i t‖ < 1/18}`, `E[N] = 17·2/18 = 17/9 ≈ 1.889`, Poisson
`P(N=0) ≈ e^{-17/9} ≈ 0.151`. Measured lonely measure:

```
initial 1..17 (tight) 0.00000   ← only exact resonance reaches 0 (boundary)
drop17 add18          0.01592
drop16 add18          0.01501
drop9  add18          0.01414
primes + 18           0.06482
```

Identical shape to n=14: generic ≈ Poisson, resonance concentrates `N` toward 0,
only the exact initial segment touches 0. The 2nd moment frames it, can't beat it.

## B — bounded ansatz (and where 18 is harder than 14)

Witness at `t = j/(18 s)`. Cofactors found:

```
initial 1..17   s=1      drop16 add18  s=3       primes+18   s=2
drop17 add18    s=2      drop9  add18  s=10  ←!   random sets s ≤ 4 (speeds<800)
```

Most sets use a *small* cofactor (`s ∈ {1,2,3}`), and over random sets the
cofactor is stable at `s ≤ 4` even up to speed 800 (no blow-up with speed size).
**But the prime-square perturbation `drop9 add18` needs `s = 10`** — markedly
larger than n=14's worst (`s = 7`). Removing the `9 = 3²` gate and replacing it
with `18` forces a finer witness (`18·10 = 180`). So:

> the prime square `9 = 3²` is exactly what makes n=18's ansatz *coarser* than
> n=14's — the richer 3-adic structure demands a larger cofactor on the
> resonant-perturbation sets.

This is the honest 18-vs-14 difference: same bounded-ansatz *shape*, but a worse
constant, and the worsening is localized at the `3²` gate.

## C — parity split: 18 = 9 doubled → a PROVED LRC@9 half

Odd `O` + even `E = 2W`; since `|O|+|E| = 17`, the smaller side has `≤ 8` speeds —
a **proved** LRC@9 instance (Trakulthongchai et al.). Both halves individually
have loneliness witnesses at threshold `1/9 ≥ 1/18`. The obstruction is the same
as at n=14: **coupling** `t` (odds) with `2t` (evens). "18 = 9 doubled," exactly
parallel to "14 = 7 doubled."

## D — triadic split: the ternary ladder (NEW vs n=14)

n=14 had two *distinct* primes; n=18's `3²` adds a genuinely new axis. Split
speeds by residue mod 3; speeds `≡ 0 mod 3` are `3w` and see `‖w(3t)‖`, folding
via `t ↔ 3t`. The chain `3 | 9` is a **ternary endpoint-debt ladder** — the exact
analog of codex's dyadic `2 | 4 | 8 | 16` ladder for n=16, one prime over. So
n=18 carries a *baby triadic version* of the n=16 dyadic debt: the `3²` gate
leaks `3`-adically just as the `2⁴` gate leaks `2`-adically.

## E — sieve gates: n=18 sits above BOTH n=16 and n=9

A counterexample must contain a speed divisible by every `m ∈ {2,…,18}`; the
maximal prime-power gates are

```
16 = 2⁴,  9 = 3²,  5, 7, 11, 13, 17.
```

So an n=18 counterexample **inherits the n=16 sedenion dyadic 16-gate** (because
`16 ≤ 18`) **and** adds a triadic `9`-gate. In the sieve tower, n=18 dominates
both n=16 (dyadic) and n=9 (triadic): its obstruction literally contains theirs.
This is the cleanest sense in which n=18 is "harder" — it is the first level that
forces *both* a 2⁴-gate and a 3²-gate at once.

## Verdict and the 18-specific program

- The three n=14 phenomena (anti-concentration wall, bounded ansatz, parity →
  proved sub-case) all recur.
- **New at 18:** the `3²` creates a *triadic ladder* (D) and *coarsens the ansatz*
  (B, `s = 10`), and the sieve forces the inherited `16`-gate (E).
- So the natural program is a **two-prime descent**: handle the `2`-part by the
  parity split (→ LRC@9, proved) and the `3`-part by the triadic `3|9` ladder,
  then CRT-couple via `t ↔ 6t` (since `6 = 2·3`). n=18 is the first frontier where
  *both* a dyadic and a triadic descent are needed simultaneously — a genuine
  two-adic-place problem, vs n=14's single odd prime.

Whimsical one-liner: **n=18 = 9 doubled, wearing n=16's dyadic gate as a hand-me-
down and a 3² triadic ladder of its own.**

## Next
1. Build the triadic `3|9` endpoint-debt ledger explicitly (mirror codex's n=16
   dyadic ladder) and compare leak depths.
2. Test the ansatz cofactor on adversarial `3²`-resonant large-speed sets — does
   `s` stay near 10 or grow?
3. Attempt the `t ↔ 6t` CRT coupling of the parity (LRC@9) and triadic halves.

## Artifacts
```
04-computation/lrc_n18_wacky_attempts_s515.py
05-knowledge/results/lrc_n18_wacky_attempts_s515.out
```
