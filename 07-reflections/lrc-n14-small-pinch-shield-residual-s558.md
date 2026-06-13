---
source: codex-2026-06-02-S558
status: proof-search; one lemma proved, direct N3 route refuted
tags:
  - lonely-runner
  - n14
  - hyp-2060
  - pinch
  - blocker-cover
  - tournament-analysis
---

# Small pinches, shields, and the remaining HYP-2060 proof residual

The HYP-2060 handoff asked for a contradiction among the necessary conditions.
The nearest promising interface is HYP-2059:

```text
counterexample => optimal binding pair has reduced sum s >= 15.
```

So if one could show that every 13-set has a reduced-sum-`<=14` pair whose
pinch time clears all runners, LRC@14 would follow.  That target is too strong.

## The failed target

Exact probes found sets with no reduced-sum-`<=14` clearing pinch.  The cleanest
HYP-2060-style example from this pass is

```text
V = (1, 2, 9, 26, 110, 153, 166, 170, 178, 190, 192, 196, 201).
```

It satisfies the checkable sieve/normalization/short-resonance conditions used
by `lrc_counterexample_necessary_conditions_s554.py`, but no small pair-cell
clears.  It is still very lonely: the exact optimum is `M=23/79` at a large
pair-pinch, so this is not a counterexample.  It only refutes the simple N3
route.

## The proved replacement lemma

The useful exact fact is THM-396.  Let `D=a+b` and
`s=D/gcd(a,b) <= 14`.  At pinch times `m/D`, the pair `(a,b)` is safe exactly
when `m` is not a multiple of `s`.  A single third speed `c` can be dangerous at
all such pair-safe residues only if

```text
D | c.
```

So universal blockers of small pinches are exactly sum-multiple shields.

This changes the proof search.  A small-pinch failure is no longer mysterious:
for each short pair-cell, either there is a shield runner divisible by the pair
sum, or the remaining runners form a finite cover of the pair-safe residues.

## Why the lemma is not enough

The stronger dream is false.  The pair `(3,12)` has `D=15`, `s=5`, and its
pair-safe residues can be covered by the five non-shield residues

```text
{14, 8, 10, 11, 13} mod 15.
```

No one of those residues is a universal blocker, but together they cover every
safe pinch residue.  Thus the residual proof problem is a blocker-cover problem,
not only a shield-existence problem.

## The HYP-2060 route after S558

The promising contradiction is now more concrete:

```text
A1-A4 sieve coverage
+ C3 short resonance
+ B1 no zero-near time
=> every short pair-cell is shielded or cover-blocked.
```

If the shield branch dominates, the speed set inherits many additive divisibility
relations `c |? a+b` in the stronger form `a+b | c`.  That looks incompatible
with primitive minimality only after a descent argument, not immediately.

If the non-shield cover branch dominates, the blocker residues must distribute
around small denominators like the `(3,12)` / `15` cover.  That should feed the
high-energy and return conditions C2/C4, but the exact contradiction is still
open.

## Handoff

Do not try to prove N3 as stated.  Try instead:

1. Prove a shield graph descent: repeated sum-multiple shields force a common
   divisor or an AP-scaled wall, contradicting F1/G4.
2. Classify minimal non-shield covers for `s<=14`, then show HYP-2060's
   mod-7 singleton coupling cannot realize all required covers at once.
3. Lift the pair-cell Tournament Analysis from the S558 script into a bounded
   exact search for minimal shield-cover certificates.
