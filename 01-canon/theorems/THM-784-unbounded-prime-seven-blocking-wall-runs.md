---
id: THM-784
title: Fast-owner refinement makes prime-seven blocking-wall runs unbounded
status: PROVED (elementary exact chamber calculation; independently rediscovered and checked by two agents) + VERIFIED (exact-Fraction script)
source: codex-2026-07-14-S10 (correction to the raw-wall conjectures in THM-779 and THM-783)
depends_on:
  - THM-773   # prime-seven token formula
  - THM-779   # blocking-wall criterion
related: [THM-778, THM-781, THM-783, HYP-6840, HYP-6845]
verification: 04-computation/lrc14_unbounded_blocking_runs_codex_S10.py
  (+ 05-knowledge/results/lrc14_unbounded_blocking_runs_codex_S10.out)
---

# THM-784 — fast-owner refinement makes blocking-wall runs unbounded

## Theorem

For every integer `N >= 1`, set

`A = {1,2,3,4,5,8,10}`, `f_N = 560N+1`, and `W_N = A union {f_N}`.

At the prime-seven lens, the whole open interval

`J = (5/16, 7/20)`

is blocked by the seven owners in `A` alone.  Moreover, `J` contains exactly
`21N` walls of `f_N`, no walls of `A`, and all of those fast walls are covered.
They are therefore `21N` consecutive covered walls in the global wall schedule
of `W_N`.

Consequently there is **no absolute constant** bounding the number of walls in
a prime-seven `r=8` blocking run.  In particular, the empirical constants 5 and
6 in THM-779 and THM-783 cannot be promoted to a universal raw-wall bound.

## Proof

On `J`, nearest-integer rounding is constant for the seven slow owners:

| `w` | range of `wx` on `J` | `round(wx)` | `w^{-1} (mod 7)` | token `-w^{-1}round(wx)` |
|---:|---:|---:|---:|---:|
| 1 | `(5/16,7/20)` | 0 | 1 | 0 |
| 2 | `(5/8,7/10)` | 1 | 4 | 3 |
| 3 | `(15/16,21/20)` | 1 | 5 | 2 |
| 4 | `(5/4,7/5)` | 1 | 2 | 5 |
| 5 | `(25/16,7/4)` | 2 | 3 | 1 |
| 8 | `(5/2,14/5)` | 3 | 1 | 4 |
| 10 | `(25/8,7/2)` | 3 | 5 | 6 |

Thus the token tuple, in the displayed owner order, is

`(0,3,2,5,1,4,6)`,

which is exactly `F_7`.  The two endpoint half-integers belong to owners 8 and
10, but the interval is open; the constancy of all seven rounding values also
shows that no slow wall lies inside `J`.  Hence these seven owners cover all
seven sheets throughout `J`.

Now `f_N = 560N+1` is distinct from the slow owners and is `1 (mod 7)`.  Its
walls in `J` are the points `(m+1/2)/f_N` satisfying

`5/16 < (m+1/2)/(560N+1) < 7/20`.

The left inequality is equivalent to `m >= 175N`; the right is equivalent to
`m <= 196N-1`.  There are therefore exactly

`(196N-1) - 175N + 1 = 21N`

such walls.  No slow wall separates them.  At each fast wall the fast owner is
absent, while the seven slow tokens remain the perfect rainbow above, so the
wall is covered by THM-779's criterion.  Between fast walls the slow rainbow
already covers every sheet.  This proves every assertion.  ∎

## What this corrects — and what survives

This theorem refutes the **abstract raw-wall-count** conclusion in THM-779 §4
and the absolute-wall conjecture in THM-783 §7.  It also refutes the claimed
equivalence between a wall-count bound and a bound of the form `extent < C/g`:
arbitrarily many walls can be inserted into one fixed interval by increasing
only the fastest speed.

It does **not** refute THM-783's conditional extent theorem.  Here the
second-fastest speed is `g=10` and

`|J| = 3/80 < 1/10 + 2/(560N+1)`.

Nor does this construction alone exhibit a five-speed core whose closed safe
component is exactly `J`; its rigorous scope is the eight-owner deck and wall
schedule.  A valid LRC14 exit theorem must therefore control **metric extent or
core incidence**, not count arbitrarily fine wall events.

## Tournament and information-boundary reading

This family separates two vertex choices that had been conflated.

- **Runner vertices.**  Use the pairwise speed comparison, gauged by increasing
  speed.  The tournament is transitive, its score histogram is
  `{0,1,...,7}`, it has no directed cycle, eight singleton SCCs, and one
  Hamiltonian path.  Changing `N` changes only the metric on the fastest
  vertex; the tournament cannot see the `21N` refinement.
- **Wall-event vertices.**  Use chronological order in `J`, with the same time
  gauge.  This is a transitive tournament on `21N` vertices, again with one
  Hamiltonian path and no cycles, but now its size records the refinement.  It
  still forgets the decisive fact that the seven non-eventing owner tokens form
  the fixed rainbow `(0,3,2,5,1,4,6)`.

So neither bare tournament quotient preserves deck blocking.  The correct
object is a **metric, owner-labelled wall-event stalk over a token chamber**:
chronological Hamiltonian path plus owner labels, mesh/interval lengths, and
the seven-token fibre.  The challenged assumption is that wall count is a
complexity measure.  Fast refinement shows it is only a choice of temporal
resolution.

