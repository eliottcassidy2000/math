# LRC n=16 Sedenion-Row Proof Search

codex-2026-05-31 S390

I tried to make `n=16` behave like the pure Cayley-Dickson/sedenion row:
not by importing algebra, but by asking what breaks when the denominator is a
pure doubling.

The clean branch is now completely visible.

For denominator `16`, if no speed is divisible by `16`, then every odd unit
point

```text
a/16, a odd
```

survives as a lonely boundary witness.  Multiplication by any nonzero residue
mod `16` sends odd residues to nonzero residues, so the distance is at least
`1/16`.  Only residue `0 mod 16` kills the unit skeleton.

So a counterexample must contain a `16`-gate.

That sounds like progress, but it also tells us where the difficulty lives:
the proof has to show that a `16`-gate cannot be inserted without creating
endpoint debt somewhere else.

## What the Exact Run Tried

The script `lrc_sixteen_sedenion_proof_search_s390.py` tried six angles:

```text
1. unit-skeleton gate lemma
2. half-turn parity
3. structured exact audits
4. one-gate replacement scan
5. dyadic endpoint-debt cascade
6. forced-gate random stress
```

The structured audits are the most revealing.

```text
initial segment:
  boundary_only, unprotected=8, layer=1

best 2-ladder:
  gap/th=0.030303, first layer=2

best 4-ladder:
  gap/th=0.015152, first layer=4

best 8-ladder:
  gap/th=0.007576, first layer=8
```

That is almost too neat.  The pure doubling row exposes debt at exactly the
half-layer you ask it to use.  The `8`-ladder is:

```text
(1, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 120)
```

It is not close to a disproof in the endpoint sense: positive gap, `140`
unprotected endpoints, terminal `coreE=0`.  But it is exactly the right lab.

## Failed Disproof Attempts

The one-gate replacement scan tested all

```text
{1,...,15} - {r} + {16q}
```

for `1<=r<=15`, `1<=q<=64`.  All `960` were positive-gap.  The closest rows
had `gap/th=0.013672`; no boundary-only or open-cover candidates appeared.

The random forced-gate stress was also dull in the best possible way:
`48/48` primitive forced-gate sets were positive-gap, and the closest endpoint
audits still peeled to `coreE=0`.

So speed-set search is not giving a counterexample.  It keeps paying debt.

## The New Proof Shape

The promising object is the endpoint of a gate.

For speed `v`, endpoints are

```text
(16m +/- 1)/(16v).
```

A speed `p` protects such an endpoint iff

```text
|p(16m +/- 1) - a*16v| < v
```

for some integer `a`.

The residue audit says:

```text
p = 0 mod 16v      protects all endpoints of v
half-gate residues protect the largest non-super chunks
lower dyadic pieces protect fragments
```

One local fact is especially sharp.  If `v=16` is treated as a maximum-speed
branch and only lower protectors `p<v` are allowed, then all `32` endpoints of
`v=16` require exactly nine lower residues:

```text
(1, 3, 5, 7, 8, 9, 11, 13, 15)
```

For `v=2,4,8`, lower residues cannot cover every endpoint at all.  So the
dyadic debt is not just poetic; it has a hard local set-cover profile.

If `v` is the largest speed in the set, there is no super-gate available.
Thus the largest gate's endpoint debt must descend to lower dyadic layers.
Those lower layers then need protection too.  Eventually the debt reaches the
unit layer, which asks for a `16`-gate again.

That is the sedenion analogy made useful: a gate closes one layer and leaks
the obstruction into the next lower layer, like zero-divisor debt.

## The Proof I Would Try Next

Assume a primitive `n=16` counterexample.

1. By the unit lemma, choose at least one `16`-gate.
2. Choose a speed that is maximal in the endpoint-debt order, not necessarily
   by integer size alone.
3. Its endpoints cannot be protected by a super-gate.
4. Charge its protected endpoints to lower dyadic layers.
5. Use half-turn parity: odd protectors cannot protect both sides of an
   antipodal pair.
6. Prove the resulting debt-flow has positive divergence unless some endpoint
   remains unprotected.

The missing lemma is a sharp inequality for this flow.  If it exists, it
should prove the pure `n=16` case before touching the messier `n=18` torsion
payload.

The slogan:

```text
n=16 LRC is a dyadic endpoint-debt theorem.
```
