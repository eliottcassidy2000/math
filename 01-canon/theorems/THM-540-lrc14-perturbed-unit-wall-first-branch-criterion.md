---
id: THM-540
title: LRC(14) one-14-multiple first-branch unit-wall criterion — exact delta windows at t = a/14 +/- delta, with a clean residue-miss corollary
status: PROVED (elementary first-branch calculation); VERIFIED exactly on bounded exhaustives (`lrc14_perturbed_unit_wall_codex.py`): the constructed witness has 0 mismatches on 1365 certified rows with A subset [1,16]\\{14}, m<=4 and on 15113 certified rows with A subset [1,18]\\{14}, m<=3
source: codex-2026-06-17-S1
depends_on:
  - THM-398   # C'(14) frame and the role of a 14-multiple
related:
  - HYP-2564  # the stronger all-branches unit-wall hypothesis
  - THM-522   # low-L extremizers now live in the unit-complete residual
  - reflection: lrc14-perturbed-unit-wall-first-branch-and-later-branches
external: lonely runner conjecture, n = 14
---

# THM-540 — the first branch off a unit wall is an exact one-inequality problem

Consider a one-multiple row

```text
S = A U {14m},
```

where `A` is a 12-set of positive integers, none divisible by `14`, and `m >= 1`.
Fix a unit `a in (Z/14Z)^* = {1,3,5,9,11,13}` and look near the unit wall `a/14`.

## Statement

Write

```text
t = a/14 + delta
```

with `delta > 0` and `delta < 1/(28m)`. For each `v in A`, let `s_v` be the residue
`a v mod 14`, taken in `{1,...,13}`. Then:

```text
||v t|| > 1/14   for every v in A
```

holds whenever

```text
delta < min_{v in A} (13 - s_v)/(14 v).
```

At the same time

```text
||14m t|| = ||14m delta|| = 14m delta > 1/14
```

holds exactly when

```text
delta > 1/(196m).
```

Therefore, if

```text
1/(196m) < min( 1/(28m), min_{v in A} (13 - s_v)/(14 v) ),
```

then every `delta` in that open interval gives a strict witness

```text
t = a/14 + delta
```

for `S`.

Symmetrically, for

```text
t = a/14 - delta,   delta > 0,
```

the exact criterion is

```text
1/(196m) < min( 1/(28m), min_{v in A} (s_v - 1)/(14 v) ).
```

## Proof

Fix `a` and `delta > 0`. For any `v in A`, write

```text
v a = 14 k_v + s_v,   s_v in {1,...,13}.
```

Then

```text
v t = k_v + s_v/14 + v delta.
```

Because `delta < 1/(28m) <= 1/28`, we are still on the first local branch off `a/14`
for the `14m`-runner, so

```text
||14m t|| = 14m delta.
```

Thus `||14m t|| > 1/14` is exactly `delta > 1/(196m)`.

Now fix `v in A`.

- If `1 <= s_v <= 6`, then `s_v/14 < 1/2`, so while we stay on this branch the
  distance to the nearest integer is

  ```text
  ||v t|| = s_v/14 + v delta.
  ```

  This remains `> 1/14` until the point where the branch reaches the opposite
  `1/14`-boundary, namely when

  ```text
  s_v/14 + v delta = 13/14,
  ```

  i.e. `delta = (13 - s_v)/(14v)`.

- If `s_v = 7`, then `s_v/14 = 1/2`, and

  ```text
  ||v t|| = 1/2 - v delta,
  ```

  so `||v t|| > 1/14` exactly while `delta < 6/(14v) = (13 - s_v)/(14v)`.

- If `8 <= s_v <= 13`, then we are already on the descending side, so

  ```text
  ||v t|| = (14 - s_v)/14 - v delta,
  ```

  and `||v t|| > 1/14` exactly while

  ```text
  delta < (13 - s_v)/(14v).
  ```

Intersecting these inequalities over all `v in A` gives the positive-branch criterion.
The negative-branch formula is identical after replacing `delta` by `-delta`; the
descending and ascending sides swap, and the endpoint becomes `(s_v - 1)/(14v)`. ∎

## Immediate consequences

1. **Residue miss is necessary on the chosen branch.**
   On the positive branch, if some `v in A` satisfies `a v == 13 (mod 14)`, then its
   upper bound is `0`, so the first positive branch fails immediately. On the negative
   branch, the obstructing class is `a v == 1 (mod 14)`.

2. **Clean corollary.**
   If `A` misses the positive-branch hurt class `-a^{-1} (mod 14)` and

   ```text
   14m > max_{v in A} v,
   ```

   then the positive first branch works. Indeed `13 - s_v >= 1` for every
   `v in A`, so

   ```text
   min_{v in A} (13 - s_v)/(14v) >= 1/(14 max A) > 1/(196m).
   ```

3. **AP one-drop corollary.**
   For `{1,...,13}\\{j} U {14m}`, the first branch discharges exactly the odd-unit drops
   `j in {1,3,5,9,11,13}` already at `m = 1`; the even drops and `j = 7` are exactly the
   unit-complete residual. This is verified exactly in the companion script/output.

## Scope / honesty

THM-540 is a **first-branch** theorem: it controls the initial branch off `a/14` and
gives an exact open interval for `delta`. It is not a complete classification of all
unit-wall witnesses of the form `a/14 + delta`. In fact, the companion script exhibits
an explicit row

```text
A = (1,2,3,4,5,6,7,8,9,10,11,26),   m = 1,
```

where the first branch fails but a later branch near `1/14` still gives a witness
`t = 657/8008 = 1/14 + 85/8008`. So the raw user intuition "missing a unit residue
should force a unit-wall witness" survives as a plausible **stronger hypothesis**
(HYP-2564), but only the first-branch version is proved here.

What is proved:

- the exact positive and negative first-branch inequalities;
- the clean corollary `14m > max A` on a branch that misses its hurt residue;
- exact computational verification of the constructed witness on all bounded rows where
  the theorem certifies a branch (`0` mismatches on `1365` certified rows with
  `A subset [1,16]\\{14}, m<=4`, and on `15113` certified rows with
  `A subset [1,18]\\{14}, m<=3`).

What remains open is the stronger global branch statement HYP-2564. The bounded
evidence is unusually strong: every tested row with a missing unit residue has **some**
unit-wall branch witness, even when the first branch fails.
