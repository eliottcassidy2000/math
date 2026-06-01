# LRC(14) as a multiplicative walk: a small-denominator reduction to fully-covered sets (S521)

*claudebox-2026-06-01-S521. Attempt at n=14 via the metagraph-walk idea, with the
"relabel the tournament by residues mod q" reframe. This is a **reduction and a
reframe, not a proof** of LRC(14); it isolates the hard core sharply. Builds on
THM-369 (denominator sieve), HYP-1991 (bounded cofactor), THM-386 (two-gap
dynamics), codex-S518 (zero-residue embryo, n=145).*

## The reframe (the "small tweak")

Stop labeling the configuration tournament by **runners**; label it by **residues
mod q**.  At time `t = a/q` the runner of speed `v` sits at residue `r = v·a mod q`
and is safe iff `||v a / q|| >= 1/14`, i.e.

```
        min(r, q-r) >= q/14,     r = v·a mod q.
```

So the forbidden set is `F_q = { r : min(r,q-r) < q/14 }`, a symmetric interval
`{0, ±1, ..., ±(⌈q/14⌉-1)}` around 0, of size `2⌈q/14⌉-1`.  The **metagraph walk**
is now the *multiplicative action of the numerator `a`* on the residue
configuration of the speeds on a **fixed** circulant tournament of `Z/q`;
loneliness at denominator `q` = some `a` rotates the whole configuration off
`F_q`.  Different speed sets are different walks on the *same* arena — exactly the
clean object we wanted.

## The provable reduction

For `q <= 14`, `q/14 <= 1`, so `F_q = {0}`: a residue is unsafe only if it is `0`.
Therefore (verified exactly over 20000 random systems):

> **Lemma (small-denominator base).** For distinct nonzero integer speeds and any
> `q` with `2 <= q <= n`, the time `t = 1/q` is lonely iff **no speed is divisible
> by `q`** (each `||v_i/q|| = min(r,q-r)/q >= 1/q >= 1/n`).

> **Reduction.** A primitive 13-speed system has a lonely time at some
> small-denominator rational `a/q` (`q <= 14`) **iff** some `q in {2,...,14}`
> divides no speed.  Hence LRC(14) holds for every system EXCEPT possibly the
> **fully-covered** ones: those in which every `q in {2,...,14}` divides at least
> one speed — equivalently the speeds must include multiples of **8, 9, 5, 7, 11,
> 13** simultaneously.

This is the metagraph walk paying off at its base layer: the `q <= 14` rotations
of the circulant arena escape `F_q = {0}` for all but the fully-covered systems,
and for those the explicit lonely time is `t = 1/q` at the uncovered `q`.  It is a
constructive, one-line witness — and it is exactly the THM-369 sieve seen as the
base of the walk.

## What is left — and why it is the right core

The residual hard core (fully-covered systems) is tiny and rigidly structured:
13 speeds must carry divisibility by `8, 9, 5, 7, 11, 13` at once.  Over 400
random fully-covered systems, **every one** was lonely at a denominator `q` with
`15 <= q <= 22` (distribution peaked at 16–17), none failing up to `q = 400`.
This is strong evidence for HYP-1991's bounded cofactor, now with a precise base:
only fully-covered systems ever need `q > 14`, and they seem to need only a hair
more.

But the reduction does **not** close LRC(14):

1. **The cover problem at `q > 14`.** For `q in [15, 27]`, `F_q = {0, ±1}`; the
   walk must rotate the residue configuration off a 3-element forbidden set in
   `Z/q`.  Proving some `a in (Z/q)*` always works for every fully-covered config,
   uniformly in the (unbounded) speeds, is the open heart.  It is a genuine
   covering/design statement, and it is where the *multiplicative* cover (which q
   divides which speed) finally meets the *additive* walk (the rotation by `a`).

2. **The speed bound.** Reducing "for all speeds" to a finite check needs a bound
   on minimal-counterexample speeds; the framework does not supply it.

So the reframe converts LRC(14) into: *the multiplicative walk on `Z/q` escapes
`F_q` for some `q <= Q(14)` and some rotation `a`, for every fully-covered speed
configuration.* The base case (`q <= 14`, `F_q = {0}`) is proved; the residual is
the `F_q = {0,±1}` cover on the fully-covered core.

## Connections

- **codex-S518 (n=145).** Same shape: "if no speed divisible by 145, THM-369
  solves immediately; the hard case is the zero-residue blockers." Here `145` is
  replaced by the *whole range* `q <= 14`, and "blockers" become "fully covered."
  The n=145 unit-wall aperture is the `q = 145` instance of this walk.
- **HYP-1991 (bounded cofactor).** Explained and sharpened: the cofactor `s>1`
  (denominator `>14`) is needed *only* for fully-covered systems, and the data say
  only just (`q <= 22`).
- **THM-386 (two-gap dynamics, opus-S519).** The residual `F_q={0,±1}` escape is
  exactly a gap-race at the polygon's near-antipodal pair; the `±1` forbidden
  residues are the two observer-adjacent gaps of THM-384/THM-386.
- **The complement / 2·7 structure.** `14 = 2·7`; the fully-covered condition
  forces the `2`-part (multiples of 8) and the `7` (multiples of 7) together —
  the first-even-bridge obstruction (the polygon's diameter) and the `7`-prime
  obstruction in one system.

## Seeds

a. **Prove the `F_q = {0,±1}` cover for `15 <= q <= 27`.** For each such `q`, show
   every fully-covered residue configuration has a unit `a` with `a·config`
   avoiding `{0,±1}` — a finite, per-`q` covering check (modulo the speed bound).

b. **Bound `Q(14)`.** Show fully-covered systems are always lonely at `q <= Q`
   for an explicit `Q` (data suggest `Q = 22` or so), reducing LRC(14) to a finite
   computation.

c. **Use LRC(7).** Fully-covered systems contain a multiple of 7; relate the `q`
   coprime-to-2 rotations to a proven LRC(7) sub-instance via the `2·7`
   factorization (the parity split, oracle-S514).

The honest one-line: **the metagraph walk, relabeled by residues mod `q`, proves
LRC(14) for every speed system except the fully-covered ones (multiples of
8,9,5,7,11,13 all present); on that rigidly-structured core the walk must escape
`{0,±1}` mod `q` for some `q`, and that residual cover — where the multiplicative
sieve finally meets the additive rotation — is the open heart, not the whole
mountain.**
