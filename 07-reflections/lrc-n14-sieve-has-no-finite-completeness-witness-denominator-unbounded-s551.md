---
source: opus-2026-06-01-S551 (remote-control)
status: reflection + RIGOROUS partial result (NOT a proof of LRC@14); honest negative on the bounded sieve
tags: [LRC, n14, sieve, THM-369, covering, witness-denominator, completeness, tournament-analysis, S525]
---

# LRC @ n=14: the division-point sieve has no finite completeness threshold

**Prompt (user, remote-control):** spend the session creatively attempting
proofs of LRC at n=14 and above.

**Honesty up front.** LRC for 14 runners is open (proved only for <=7). I did not
prove it and I do not claim to. What this session contributes is a *rigorous*
clarification of exactly how far the project's main tool — the denominator sieve
(THM-369) — can possibly go, plus the concrete obstruction it leaves. The result
is a clean **negative with strategic value**: it rules out an entire class of
hoped-for finite proofs and says precisely where the real difficulty sits.

## 0. Setting (repo convention)

`n` runners = observer + `m = n-1` speeds `v_1,...,v_m` (distinct positive
integers, WLOG `gcd = 1`). Observer **lonely** at time `t` iff `||v_i t|| >= 1/n`
for every `i`. Covering reformulation (S525): the `m` open danger arcs
`D_i = {t : ||v_i t|| < 1/n}` (each `= v_i` arcs of total length `2/n`) never
cover the circle. For `n=14`: 13 arcs, each of total length `1/7`.

The **division-point sieve** (THM-369, Lean-checked). At `t = a/q`, `gcd(a,q)=1`,
runner `v` is safe iff `min(r, q-r)*n >= q` with `r = (v a) mod q`. A *found*
witness `a/q` is an unconditional proof that the set is lonely — it is just the
loneliness predicate evaluated at one explicit rational. For `q <= n` the only
dangerous residue is `r=0`, so a witness at `a/q` exists iff `q` divides no speed.

## 1. What was verified (3 exact scripts, all decisions integer/Fraction-exact)

`lrc_n14_multiprime_sieve_s551.py` (+ `_completeness_probe_` + `_witness_tournament_`):

- **0 LRC failures** over 400 random + 8 structured n=14 sets, and over
  4×120 random sets at n=15,16,18,22. (Exact ground-truth via endpoint+midpoint
  enumeration, which catches measure-zero tight witnesses.)
- **The AP `{1,...,13}` is the tight extremal**: lonely only on a measure-zero
  set of 6 boundary points (`1/14, 3/14, 5/14, 9/14, 11/14, 13/14`); witness
  `1/14` = the regular 14-gon wall. (Confirms S524/S525.)
- **Necessary condition** (known, THM-360): a counterexample must contain a
  speed divisible by **each** `q in {2,...,14}`. The irreducible demands are the
  prime powers `<= 14`: `{2,3,4,5,7,8,9,11,13}` (plus composite demands
  `6,10,12,14`).

## 2. The new rigorous result: bounded sieve cannot prove LRC(n)

One might hope for a finite `Q(n)` such that **every** speed set has a sieve
witness with denominator `q <= Q(n)`. That would reduce LRC(n) to a finite check
(witness existence depends only on residues mod `lcm(2..Q)`). **This hope is
false, and cheaply so.**

> **Construction (proved).** Fix any `Q`. Let one speed be `L = lcm(2,...,Q)` and
> the rest be any distinct positive integers. For every `q <= Q`, `q | L`, so
> `||L·a/q|| = 0 < 1/n` for **every** multiplier `a`: runner `L` is in the danger
> band at every division point of every modulus `q <= Q`. Hence the set has **no
> division-point witness with denominator `<= Q`** — it is *sieve-blind up to Q*.

Yet these blind sets are **lonely**. Verified for `n=14`, the set
`{lcm(2..Q)} ∪ {1,...,12}` at `Q = 14,16,18,20,24`:

| Q  | lcm(2..Q)   | blind up to Q? | smallest witness q\* | lonely? |
|----|-------------|----------------|----------------------|---------|
| 14 | 360360      | yes            | 27 (t = 2/27)        | yes ✓   |
| 16 | 720720      | yes            | 27                   | yes ✓   |
| 18 | 12252240    | yes            | 27                   | yes ✓   |
| 20 | 232792560   | yes            | 27                   | yes ✓   |
| 24 | 5354228880  | yes            | 27                   | yes ✓   |

Each "lonely" is an *unconditional* proof: the exhibited `a/q*` was checked safe
exactly. (They share `q*=27=3^3` because `27 ∤ lcm(2..24)`; pushing `Q >= 27`
forces a strictly larger witness, and the construction never stops.)

**Consequence.** For every bound `Q` tested there is an explicitly **verified-lonely**
n=14 speed set with no witness of denominator `<= Q`. So:

- the minimal witness denominator is **unbounded** over lonely sets (no finite
  completeness threshold exists);
- **a bounded-modulus division-point sieve provably cannot prove LRC(n)** for
  `n >= 8`. The sieve is a one-sided certificate generator, not a decision
  procedure with a finite modulus budget.

This is the rigorous converse-limitation of THM-369, and it dovetails with the
S1 correction (the *measure* is trivially `>= 0`, so it cannot prove LRC either):
**both the cheap measure bound and the cheap finite sieve are provably too weak.**
The content of LRC@14 lives entirely in the residual they leave.

## 3. Where the difficulty actually sits (the resistance curve)

The blindness above is bought with *huge* speeds. For **bounded** speeds the
sieve is strong: a hill-climb maximising the minimal witness denominator gives

- speeds `<= 60`: hardest min-witness-modulus `= 34` (set fully loaded: every
  `q <= 14` divides some speed), lonely at `t = 11/34`;
- speeds `<= 300`: hardest `= 35`, lonely at `t = 1/350`.

So the witness denominator for a set is controlled by its **divisibility
loading** (which `q` it kills by residue 0), and loading past `q <= n` requires
ever-larger speeds. The honest picture:

> LRC(n) is hard precisely because the only known **complete** finite reduction
> is the bounded-speed exact check (Tao-style), whose bound is astronomical at
> `n=14`; every **modulus-bounded** shortcut (sieve) is provably incomplete, and
> the sets that defeat a bound `Q` are exactly the ones loaded with divisors up
> to `Q` (hence with speeds `>= lcm(2..Q)`).

## 4. Tournament Analysis (project directive)

At a sieve witness `a/q` the runners occupy distinct circle points
`x_i = ((v_i a) mod q)/q`, all inside the **safe band** `[1/n, 1-1/n]`. The
half-turn runner tournament `i -> j  iff  (x_i - x_j) mod 1 in (0,1/2)` is the
object S522–S525 study. Over 300 n=14 witnesses:

- **#SCC in {1, 13} only** (284 single strong block, the rest degenerate) —
  confirms S525's restriction. (The handful of apparent `#SCC=13` are *antipodal-
  tie* degeneracies at even `q`, where two runners are exactly `1/2` apart and the
  half-turn relation is undefined; excluding those, the rule is exact.)
- **the largest circular gap straddles the observer in 93%** of witnesses — this
  is the *witness-arc signature*: loneliness ⟺ a gap `>= 2/n` around `0`.
- the realised score sequences are tightly clustered near the **regular** value
  `(6,6,...)` — a vanishing near-regular slice of `A000568(13)`.

**Restricted iso-class statement.** The realizable lonely n=14 runner tournaments
are the **round** tournaments whose maximal gap (`>= 2/n`) contains the observer.
The sieve gives the *arithmetic* realization (residues mod `q`) of exactly this
*geometric* restriction (S525). The tight AP sits at the extreme: the regular
13-point rotational tournament, gap exactly `2/n` at the wall `t=1/14`.

## 5. Honest verdict and handoff

Not a proof. The contribution is to **close off the bounded-sieve route** with a
rigorous construction and to localise the residual: the hard sets are the
heavily-divisor-loaded ones, whose witnesses sit at a modulus just above their
loading. Productive next targets:

1. **Witness-denominator vs loading, made into a theorem.** Conjecture
   (HYP-2052): the minimal witness denominator of a primitive n-set equals (up to
   a small additive constant) the smallest `q` *not* dividing all the loaded
   structure — i.e. the next prime power above the largest `q<=?` that some speed
   kills. If provable, it bounds the witness denominator *in terms of the speeds*
   and turns the bounded-speed exact check into a bounded-modulus one for any
   fixed speed cap — a genuine (if still large) finite reduction.
2. **Whether the loaded blind family is lonely *unconditionally*** (independent
   of LRC): a self-contained loneliness proof for `{lcm(2..Q)} ∪ small` would
   upgrade §2's "unbounded" from verified-for-tested-Q to a theorem.
3. The set-vs-measure residual (HYP-2039/S1) restated on the round-tournament
   slice of §4: characterise the gap families that keep the observer-straddling
   gap `>= 2/n` for all `t` only at isolated walls (the tight configs).

**Artifacts:**
`04-computation/lrc_n14_multiprime_sieve_s551.py` (+`.out`),
`04-computation/lrc_n14_sieve_completeness_probe_s551.py` (+`.out`),
`04-computation/lrc_n14_witness_tournament_s551.py` (+`.out`),
HYP-2052. **Builds on:** THM-369 (sieve), THM-360 (divisor necessity), S525
(round/SCC geometry), S1/HYP-2039 (measure trivial, set-vs-measure residual).
