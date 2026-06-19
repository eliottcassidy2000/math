# LRC14 perturbed unit walls: the first branch is exact, later branches are real, and the hard core lines up with unit-complete rows

**Source:** codex-2026-06-17-S1. Prompt: long session on LRC14 built around the
idea that a witness near `a/14` should survive whenever the non-14 part misses a
unit residue.

## The clean proved part

For a one-multiple row

```text
S = A U {14m},
```

the first branch off a unit wall is completely rigid.

Write

```text
t = a/14 + delta
```

with `a` a unit mod `14` and `0 < delta < 1/(28m)`. Then for each `v in A`,
letting `s = a v mod 14` in `{1,...,13}`,

```text
||v t|| > 1/14   exactly while   delta < (13 - s)/(14v).
```

For the 14-multiple itself,

```text
||14m t|| = 14m delta > 1/14   exactly while   delta > 1/(196m).
```

So the first branch is a literal open interval problem:

```text
1/(196m) < delta < min_v (13 - a v mod 14)/(14v).
```

The negative branch is the same with `(s-1)/(14v)`. This is THM-540.

The striking simplification is that the residue effect is not just on the
boundary classes `±1`. Positive `delta` helps every `s in {1,...,6}` and hurts
every `s in {7,...,13}`; negative `delta` reverses that. The user's boundary
picture was the right seed, but the whole first branch is one residue-slope
table.

## The assumption that broke

The tempting overstatement was:

```text
missing unit residue  =>  the a/14-perturbation works immediately.
```

That is false as a **first-branch theorem**.

Example:

```text
A = (1,2,3,4,5,6,7,8,9,10,11,26),   S = A U {14}.
```

`A` misses unit residue `13`, but THM-540 certifies no first branch at any unit
wall.

The obstruction is not the missing unit class; it is the **delta window**. The
26-runner makes the first positive window near `1/14` too short:

```text
1/364 < delta < 1/196
```

has no solution. So "missing a unit residue" must be paired with a quantitative
window statement; residue information alone is not enough.

## The surprise: later branches rescue the failed local rows

The same counterexample still has a unit-wall witness:

```text
t = 657/8008 = 1/14 + 85/8008.
```

So the branch picture is richer than the first linear window. A speed can become
unsafe on the first branch, wrap, and then return to safety on a later branch.

This is the key structural update from the session:

```text
first branch = theorem
later branches = real phenomenon
```

The exact script confirms this boundedly:

- among all `1372` rows `A subset [1,16]\\{14}`, `|A|=12`, missing a unit residue,
  `m <= 4`, every single row has **some** unit-wall branch witness;
- among all `15666` rows `A subset [1,18]\\{14}`, `|A|=12`, missing a unit residue,
  `m <= 3`, again every single row has **some** unit-wall branch witness.

That is HYP-2564: the user’s stronger criterion survives computationally, but now
as an all-branches statement, not a first-order perturbation lemma.

## Why this matters for the real frontier

The new exact low-L story (THM-522 / HYP-2561) says the current extremizer is

```text
{1,2,...,11,13,36},
```

and this row is **unit-complete**. So the two independent proof programs are now
pointing at the same residual:

```text
unit-incomplete rows -> unit-wall branch witnesses
unit-complete rows   -> the genuine hard core
```

That is the real payoff. The singular-series/measure extremizers are not hiding
among the easy missing-residue rows. They already sit in the unit-complete
residual the perturbation lemma does not touch.

## Assumption challenge

I explicitly considered three quotients before settling on the right one:

1. **Unit residues.**
   Preserves which `a/14` walls are even available, but destroys the actual
   branch widths.
2. **First-branch residue slopes.**
   Preserves the exact local inequality; this is THM-540. Destroys later-wrap
   rescues.
3. **Full unit-wall branches.**
   Preserves the witnessed phenomenon in bounded data, but no clean closed-form
   theorem fell out in this session.

The challenged assumption was:

```text
missing a unit residue should already solve the local branch.
```

The correction is:

```text
missing a unit residue + a nonempty delta window solves the first branch;
missing a unit residue alone may still solve a later branch.
```

That is a materially better statement than where the session started.
