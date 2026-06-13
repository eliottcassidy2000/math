---
source: codex-2026-06-03-S591
status: support-calculus formalization extending THM-401/HYP-2132; exact named-row audit; HYP-2135
tags: [LRC, sumset-support, THM-401, HYP-2132, pair-sum, support-sheaf, nonunit-holes, denominator-shields, tournament-analysis]
---

# Sumset support after THM-401: the missing labels are the proof object

The user asked to work on sumset support and to formalize/extend THM-401 +
HYP-2132.  The clean realization is that THM-401 gives the correct base
modulus `C=2n-1`, but "support" is not one thing.

There are at least five support layers:

```text
speed shells        which antipodal shells of Z/C are occupied by speeds;
pair shells         which antipodal shells are hit by pair sums;
denominators        which actual pinch denominators D=v_i+v_j appear;
shields             which D are killed by a speed divisible by D;
visibility/lifts    which holes are unit-invertible, nonunit, or outside C.
```

This changes the slogan from:

```text
study V+V
```

to:

```text
study the labelled projection V -> shells <- pair sums, with divisibility and unit action retained.
```

That is more rigid and much less glamorous, which is usually a good sign in
this repo.

## What S591 Shows

`04-computation/lrc_sumset_support_calculus_s591.py` is intentionally small.
It audits AP rows, the S573 open-gap rows, and n=14 boundary rows.

The AP identity is now crisp:

```text
AP_n={1,...,n-1}
speed shells: perfect transversal
pair shells:  everything except shell 1
low D<=n:     3,4,...,n
shields:      every D<n is shielded by speed D
defect:       D=n is unshielded
```

So AP is not "maximal pair support."  It is "perfect speed support plus one
precisely unshielded floor denominator."

The lifted rows explain why the earlier floor-only classification broke:

```text
open_gap_n7_S573
  speed shells perfect
  pair shells complete
  exact witness denominator 33
  M=5/33 in (1/7,2/13)

nonunit_hole_n8_A/B
  no unit holes
  missing speed shell 6 has gcd(6,15)=3
  exact witness denominator 23
  M=3/23 in (1/8,2/15)
```

The `C`-support ledger sees why unit inverse clocks do not close these rows,
but it does not yet see the final witness until we lift to the exact pair
denominator and endpoint/DUN labels.

For n=14:

```text
AP_n14:
  perfect speed support, pair shell 1 missing, D=14 unshielded, M=1/14.

Vstar_n14:
  misses nonunit shell 12, doubles nonunit shell 3, pair shells complete,
  D=14 unshielded, M=1/14.

doubled_apex_edge_n14:
  misses unit shell 13, so t=25/27 gives the exact edge M=2/27.

unit_shift_AP_n14:
  misses unit shell 1, but speed 14 shields D=14; the unit inverse gives only
  the 2/27 edge lower bound, while exact M=1/8.
```

That last pair is important.  A unit hole guarantees an edge witness, but it
does not decide exact looseness after denominator shields change the local
pinch geometry.

## The `n=11,12,13` Efficiency Angle

The shell strata around the paper frontier are:

```text
n=11 C=21: nonunit strata gcd 3 and gcd 7.
n=12 C=23: prime, every shell is unit-visible.
n=13 C=25: nonunit stratum gcd 5.
n=14 C=27: nonunit strata gcd 3 and gcd 9.
```

This is the arithmetic rotation from earlier sessions in a sharper form.
The paper's exact finite sieve likes primes in a different place; the repo's
certificate ledger likes prime `C=2n-1`.  So total `n=12` should be unusually
cheap for this support calculus: any missing shell is immediately unit-visible.

For `n=11,13,14`, the only rows not killed by unit-hole clocks must hide in
nonunit fibers or lift denominators.  That is a much smaller search target
than raw speed tuples.

## New Proof Picture

The support calculus wants a sheaf:

```text
base:
  antipodal speed shells of Z/(2n-1)

stalk over a shell:
  speeds occupying it, pair denominators touching it, shields, endpoint owners,
  D/U/N obligations, exact witness denominators

gluing:
  unit multiplication, reflection, CRT decomposition, pair-denominator lifts

defect:
  nonunit ramification or monodromy after an endpoint label is forgotten
```

Then the proof algorithm becomes:

```text
1. If a unit shell is missing, use the inverse C-clock for a 2/C witness.
2. Otherwise, if all speed shells are perfect, inspect the unshielded low D layer.
3. If D=n survives, route to AP/Vstar-like floor pincer.
4. If D=n is shielded or an exact witness denominator lifts beyond C, attach
   endpoint-owner and D/U/N labels and discharge by HYP-2088/HYP-2122.
5. If only nonunit holes remain, descend by gcd/CRT or lift the ramified fiber.
```

The wild part: raw `V+V` may be the wrong object; `V -> shell(V)` and
`pairs(V) -> Den(V) -> shell(Den(V))` form a pushout diagram.  The information
is not in either side alone, but in the failure of the two supports to glue
after shields and units act.

## Tournament Analysis Note

S591 uses proof lenses as tournament vertices.  The ranking is:

```text
low_denominator_shields
> unit_visible_holes
> witness_denominator
> speed_shell_transversal
> nonunit_holes
> gcd_strata
> pair_sum_shell_holes
> raw_sumset_size
```

The tournament is transitive.  The edge flips versus cost-only ranking say the
same thing the examples say: cheap support size is not proof power.  The
payload lives in whether the support carries the observer-coupled labels:
denominator, shield, unit inverse, endpoint owner, or lift.

## What To Try Next

1. Build a 64-row n=14 fixed-boundary support-sheaf table: speed shell,
   pair-shell, low denominator, shield, unit/nonunit, D/U/N, endpoint owner.
2. For `n=12`, use prime `C=23` to prove the support calculus kills every
   non-perfect row by unit inverse clocks before any paper-style prime sieve.
3. For `n=11,13`, isolate only nonunit-hole and lift-denominator residuals,
   then compare their counts to S580's certificate compression estimates.
4. Make the lift denominator a first-class automaton state: `C-clock`,
   `n-clock`, `D=v_i+v_j clock`, `endpoint/Phi clock`.

The moral: after THM-401, the object is not the sumset.  It is the labelled
support-residue left after projecting the sumset through addition,
multiplication by units, and divisibility.

**Artifacts:** HYP-2135;
`04-computation/lrc_sumset_support_calculus_s591.py`;
`05-knowledge/results/lrc_sumset_support_calculus_s591.out`.
