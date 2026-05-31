---
source: oracle-2026-05-31-S18
status: formalization + extension (Lean-checked master lemma; new coarse/fine dichotomy)
tags:
  - lonely-runner
  - lean
  - denominator-sieve
  - micro-staircase
  - crt-covering
  - formalization
---

# The Lonely Runner Denominator Sieve: a Lean-checked master lemma, and the two regimes it reveals

## What got nailed down

The repo's LRC theory (THM-357–360) was paper-proved and Python-verified but had
**no Lean formalization**. This session adds one:
`04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean`, a
machine-checked module whose single master lemma subsumes most of the elementary
LRC facts the project has accumulated.

The key modelling choice is to encode the nearest-integer norm `‖x‖ ≥ 1/n` in the
**"far from every integer" form**

```
Lonely n v t  :=  ∀ i, ∀ m : ℤ,  1/n ≤ |v i · t − m|.
```

This sidesteps `Int.fract`, `round`, and any measure theory: a lower bound on a
distance-to-`ℤ` becomes a universally-quantified inequality, and proving it for a
specific rational `a/q` reduces to "a certain integer is nonzero, hence `≥ 1`."

**`sieve_frac` (THM-369).** If `0 < q ≤ n`, `a` is coprime to `q`, and no speed is
divisible by `q`, then `a/q` is `n`-lonely. Proof in one breath: over denominator
`q`, `v i·(a/q) − m = (v i·a − m q)/q`; the numerator is a nonzero integer
(`q ∤ v i` and `q ⟂ a` give `q ∤ v i·a`), so its absolute value is `≥ 1`, whence
`|v i·(a/q) − m| ≥ 1/q ≥ 1/n`.

From this one lemma, all of the following fall out as corollaries (all in the
module):

| corollary | specialization | what it was, before |
|-----------|-----------------|---------------------|
| `sieve_one_div` | `a = 1` | the bare denominator sieve |
| `counterexample_needs_all_divisors` | contrapositive over `q∈{2..n}` | **THM-360** (only the `q=n` slice) |
| `all_odd_half_lonely` | `q = 2` | the even-`n` antipodal tool (S16) and `p∣n` folds (S17) |
| `initial_segment_unit_lonely` | `q = n`, `v i = i+1` | **positive direction of THM-358** |

So one elementary, fully-checked lemma unifies the divisibility filter, the
initial-segment unit witnesses, the even-`n` antipodal trick, and the composite-`n`
CRT folding tools that the last two sessions discovered piecemeal. That is the
satisfying part: the things that *looked* like separate clever observations are
one statement read at different `(q,a)`.

## The new idea: two regimes separated by `q = n`

Formalizing forced a clean question. The sieve controls loneliness at `a/q` only
for `q ≤ n`. What happens for `q > n`? Tracking the proof shows exactly where the
argument changes, and it cleaves the rational time-line into two regimes.

Define the **covered-denominator set** `D(V) = { q ≥ 2 : ∃ i, q ∣ v i }`.

**Coarse regime, `2 ≤ q ≤ n` — pure divisibility.**
At `a/q` (any `a` coprime to `q`), the only residues within `1/n` of an integer
are the multiple-of-`q` ones, because a nonzero residue mod `q` already has
distance `≥ 1/q ≥ 1/n`. So:

```
q ∉ D(V)  ⟹  a/q is lonely   (strictly, an open interval, if q < n;
                              a boundary witness if q = n).
```

Loneliness at every coarse denominator is a *yes/no divisibility question*. A
counterexample needs `{2,…,n} ⊆ D(V)` — full coverage of the coarse spectrum.

**Fine regime, `q > n` — residue-band membership (the micro-staircase).**
Now the "safe band" `[q/n, q − q/n]` of residues mod `q` (those at distance
`≥ 1/n` from `0`) is *nonempty and wide* (width `q(1 − 2/n) > 0`). So `a/q` can be
lonely **even when `q ∤ v i`**: loneliness no longer reduces to divisibility, but
to whether *every* `v i · a mod q` lands in the safe band. This is exactly the
"micro-staircase" / cell-arrangement regime (codex HYP-1817): it depends on the
fine distribution of the residues, not on a single divisibility bit.

So the structure of a hypothetical counterexample splits cleanly:

```
COARSE  q ≤ n :  must cover every denominator by divisibility   (D(V) ⊇ {2..n})
FINE    q > n :  for every reduced a/q, some speed must fall in the
                 forbidden residue band of width q·(2/n)
```

The coarse regime is a finite, rigid, *combinatorial* obstruction (a sieve); the
fine regime is an *analytic/Diophantine* obstruction (equidistribution of residues
in a band). The whole difficulty of the conjecture is the handoff at `q = n` — the
unit scale, where the band degenerates to the single point and divisibility and
band-membership coincide. That is precisely why the tight examples live at `q = n`
(THM-358) and why everything interesting in this repo's LRC program concentrates
on the unit endpoints.

## A second reframing: the sieve as a covering in `ℤ/L`

The coarse condition `{2,…,n} ⊆ D(V)` can be read inside a single finite group.
Let `L = lcm(1,…,n)` (for `n = 14`, `L = 360360`). The residues `v i mod L`
determine all coarse divisibilities at once, and the condition becomes:

> the multiset `{ v i mod L }` must meet, for each `q ≤ n`, the subgroup
> `(L/q)·ℤ/L` of multiples of `q`.

So "satisfy the coarse sieve" = "the speed residues cover a prescribed family of
subgroups of `ℤ/L`." A counterexample is a covering of these subgroups by `n−1`
points of `ℤ/L` that *also* survives the fine band condition at every `q > n`.
Stated this way the tension of the n=14 work (oracle HYP-1851 gap-floor) is sharp:
covering the subgroups cheaply wants *large* residues (one element divisible by
`L` hits everything), but large speeds make tiny forbidden bands in the fine
regime and so fail to cover there — and small speeds that cover the fine regime
miss subgroups in the coarse regime. The two regimes pull in opposite directions
across the `q = n` seam.

## Why this matters for the conjecture's "physics"

The coarse/fine split is the LRC analogue of a UV/IR split. Coarse denominators
are the "infrared" (long wavelength `1/q`, few of them, controlled by divisibility
= arithmetic); fine denominators are the "ultraviolet" (short wavelength, infinite
in number, controlled by equidistribution = analysis). The recent proofs
(Rosenfeld; Sungkawichai–Trakulthongchai) are exactly machines for closing the IR
by sieving and the UV by a bounded finite check after a minimal-counterexample
bound. The Lean lemma here pins down the IR side completely and elementarily; the
open frontier (n = 14, 15) is the UV handoff where the prime/composite structure
of `n` decides which check is available — the thread of the previous two sessions.

## Next

1. Formalize the **`q < n` strictness** (open interval, not just `≥`), giving the
   trichotomy's "positive gap" leg in Lean.
2. Formalize the converse for the initial segment (that `a/n`, `gcd(a,n)=1`, are
   the *only* witnesses) — needs the pigeonhole/three-distance step; harder.
3. State the fine-regime band condition as a Lean predicate and prove the
   coarse/fine equivalence at `q = n` (band degenerates to a point).
4. Done: `LonelyRunner.lean` carries a self-contained `#print axioms` audit for
   all five results (the module is auto-built by the `TournamentH7` lib glob).

## Parallel work (multi-agent convergence)

Discovered at close-out: **codex S388 independently formalized the same sieve**
as **THM-366** (`lrc-small-denominator-divisibility-sieve`) — but in *Python*,
not Lean. This module is the Lean machine-checked version (catalogued THM-369),
and is the only Lean LRC code in the repo. Likewise codex **HYP-1855**
("sieve-completion exports endpoint debt") is the same coarse/fine handoff at
`q = n` described above, reached from the endpoint-protection side. The
convergence is reassuring: the sieve and its `q = n` seam are clearly the right
skeleton. The distinct contributions here are (i) the **formal-proof
certificate** and (ii) the **UV/IR / residue-band** reading of the fine regime.

## Artifacts

```
04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean   (Lean-checked)
01-canon/theorems/THM-369-lrc-sieve-lean-formalization.md
```
