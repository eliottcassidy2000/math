# LRC n=16 Antipodal/Fan Proof Search

codex-2026-05-31 S392

This pass attacked the `n=16` Lonely Runner case from a different angle than
the previous sedenion/debt runs.  The goal was not another speed-set search,
but a new certificate that a counterexample would have to satisfy.

Two reframes survived.

## Antipodal Pair Quotient

At `n=16`, split the speeds into even and odd.

Under the half-turn `t -> t+1/2`, even speeds are unchanged.  Odd speeds move
by `1/2 mod 1`, so a single odd speed cannot protect both sides of an
antipodal pair.

Therefore any open cover must cover every `x in [0,1/2]` by

```text
even-forbidden(x)
OR
(odd-forbidden(x) AND odd-forbidden(x+1/2)).
```

This is an exact quotient certificate.  It is weaker than the original
problem, but any counterexample must pass it.

The initial segment `{1,...,15}` passes the quotient only as a boundary case:

```text
pair_measure=1/2
pair_gap=0
pair_boundary_witnesses=4
```

That matches the original LRC behavior: the initial segment is tight but not
an open cover.

The gated rows fail more strongly.  In the structured audits:

```text
drop 15 add 16:              pair_gap/th=0.038462
best 8-ladder:               pair_gap/th=0.007576
odd units plus gates:        pair_gap/th=0.071429
fan64 plus odd breakers:     pair_gap/th=0.202083
fan128 plus small breakers:  pair_gap/th=0.101042
```

The one-gate antipodal scan checked all

```text
{1,...,15} - {r} + {16q}, 1<=r<=15, 1<=q<=32.
```

All `480` rows had positive antipodal pair gaps.  The closest rows had
`pair_gap/th=0.027344`.

So the half-turn quotient is not just a slogan.  It gives an independent
obstruction before the full endpoint system is even built.

## Maximal Dyadic Fan

The previous endpoint-debt idea suggested choosing a largest speed `v` and
arguing that its endpoints cannot all be protected by lower speeds.

That naive proof is false in a very structured way.

For dyadic `v`, the endpoints of speed `v` are all covered by the nine-speed
fan

```text
v/2, (v/32)*1, (v/32)*3, ..., (v/32)*15.
```

The exact audit:

```text
v=16   covers   32/32 endpoints, gcd=1
v=32   covers   64/64 endpoints, gcd=1
v=64   covers  128/128 endpoints, gcd=2
v=128  covers  256/256 endpoints, gcd=4
v=256  covers  512/512 endpoints, gcd=8
v=512  covers 1024/1024 endpoints, gcd=16
```

For `v>=32`, the normalized fan is always

```text
(16,1,3,5,7,9,11,13,15).
```

This kills one possible proof route but opens a better one: a large fan is
imprimitive.  It can close the top endpoints only by concentrating into a
scaled odd residue fan.  The remaining six speeds must break the gcd and also
handle the fan's descendant endpoint debt.

The stress examples did not get close to a cover.  Fan rows with primitive
breakers still had ordinary positive gaps, positive pair gaps, unprotected
endpoints, and `coreE=0`.

## New Proof Target

A counterexample now has to satisfy three simultaneous constraints:

```text
1. it has a 16-gate, by the unit-skeleton lemma;
2. it has no antipodal pair gap;
3. every large maximal dyadic closure either contains the scaled fan or pays
   more endpoint debt elsewhere.
```

The next proof attempt should try to prove that constraints 2 and 3 are
incompatible.

Odd double-coverage wants two odd lanes on every pair not already covered by
evens.  The maximal fan wants a scaled odd residue fan plus a half-gate.  Once
primitivity forces gcd-breaker speeds into the set, those two demands appear
to pull in opposite quotient directions.

The compact slogan:

```text
n=16 LRC is an antipodal fan incompatibility problem.
```
