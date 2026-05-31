---
source: codex-2026-05-31-S399
status: exploratory reflection
tags:
  - lonely-runner
  - n16
  - bruhat-tits
  - dyadic-tree
  - harmonic-flow
---

# LRC and the Bruhat-Tits Tree

The useful correction is small but sharp:

```text
The n=16 dyadic endpoint law is not just "2-adic."
It is a conservative kernel on the 2-adic tree.
```

For the Bruhat-Tits tree of `PGL_2(Q_2)`, boundary points are 2-adic
directions.  The pure dyadic LRC endpoints

```text
(16m +/- 1)/(16u),    u=2^k
```

are finite boundary cylinders.  Owner `u=2^k` is a horosphere.  Protector
`p=2^j q` is a move whose tree drop is `L=k-j`, and whose boundary direction is
the odd residue `q mod 16`.

THM-367 then becomes a finite Hecke-kernel statement.

## The Kernel

After normalizing by the `2u` endpoints in the owner horosphere, S399 gives:

```text
L = 1:  active q = {1,15},                 capacity per active q = 1/2
L = 2:  active q = {1,3,13,15},            capacity per active q = 1/4
L >= 3: active q = all eight odd classes,  capacity per active q = 1/8
```

In all three regimes:

```text
number of active residues * capacity per residue = 1.
```

That is the new small answer.  THM-367 is a mass-preserving transition rule.
The shallow layers are ramified: they use fewer boundary apertures with larger
capacity.  The deep layers are spherical: all odd residues are active, each
with one eighth of the endpoint horosphere.

So the local theorem is already a Markov/Hecke operator hiding inside endpoint
arithmetic.

## The Lower Truncation

A maximal speed cannot use protectors above itself.  This replaces the full
symmetric kernel by the one-sided lower-protector view:

```text
drop 1: {1}
drop 2: {1,3}
drop 3: {1,3,5,7}
drop >=4: all odd residues
```

This is why the proof feels recursive.  At shallow drops the maximal owner sees
only half the apartment.  By drop `4`, it sees the full odd boundary.  The first
place where a lower cover can be exact is `u=16`, precisely when that full odd
boundary becomes available.

This also explains the shape of the exact lower cover:

```text
(1,3,5,7,8,9,11,13,15).
```

It is all eight odd boundary directions plus the one radial half-turn `8`.

## The Fan

For `u>=32`, the constructive nine-cover is:

```text
u/2 + (u/32)*{1,3,5,7,9,11,13,15}.
```

After dividing by its gcd, it is always:

```text
(16,1,3,5,7,9,11,13,15).
```

So the fan is one radial direction plus the complete eight-direction boundary
sphere five levels below.  Its total endpoint incidence is always `3u`, hence
average indegree `3/2`.

The fan has two apex endpoint cylinders hit by all nine speeds.  Those apexes
look like finite Bruhat-Tits vertices where every branch of the star meets.  It
also has degree-2 bands and degree-1 leaves.  This is a tree current, not a
random set cover.

The catch is the gcd.  For `u>=64`, the fan is imprimitive.  That is exactly
where the proof should bite: local harmonicity is easy, primitive harmonicity
is the hard part.

## Why Private Leaves Stopped Working

At `u=16`, private endpoints relative to all lower protectors force the exact
nine residues.  S399 counts `24` such private endpoints.

At `u>=32`, there are no private endpoints relative to all lower protectors.
The all-lower incidence graph has the full nonempty core `(2u,u-1)`.

So the private-leaf proof is not false.  It is the base case.  Higher dyadic
layers replace private leaves by harmonic currents.  The missing theorem should
not say:

```text
every dyadic horosphere has private endpoints.
```

It should say:

```text
every primitive finite harmonic current has positive divergence somewhere.
```

That divergence can show up as an actual interval gap, an endpoint leaf after
global restrictions, or an Egyptian/q defect in the S397 sense.

## Link Back to Tournament Work

This mirrors several old tournament patterns in the repo.

Endpoint transfer work kept finding that local private witnesses prove the
first layer, while higher layers require rank, core, or flow invariants.  The
n=16 LRC row has the same transition:

```text
private endpoints at u=16
core/current structure at u>=32
```

The Bruhat-Tits tree is also the valuation skeleton of the multiplication
mode graph on natural numbers.  Addition `X+Y=Z` collapses to an order shadow.
Multiplication `X*Y=Z` collapses to divisibility.  Taking `v_2` of the
multiplicative shadow gives the tree depth; keeping odd residue `q mod 16`
keeps the boundary direction.  LRC endpoint protection is a metric thickening
of this natural-number mode graph.

That puts S397 and S399 together:

```text
Egyptian fractions: reciprocal mass splits by divisor residue.
Bruhat-Tits view: dyadic mass flows through boundary residue directions.
LRC: interval coverage demands both split mass and primitive tree flow.
```

The product-sum equations `X+Y=XY`, `X+Y+Z=XYZ`, and so on are the equality
version of the same balance problem.  LRC is the inequality version.

## New Search Shape

A disproof construction for `n=16` should not be imagined as a lucky dense
interval cover.  It would have to be a finite harmonic current on a truncated
2-adic tree, with:

1. all endpoint horospheres receiving enough incoming mass;
2. lower-protector truncation respected at maximal speeds;
3. primitive gcd equal to `1`;
4. no positive archimedean gap left over.

The existing fan passes the local current test and fails the primitive/global
test.  That is the right toy obstruction.

The next concrete experiment should assign every speed a Bruhat-Tits state:

```text
(v2(speed), odd_part mod 16, gcd class, fan-apex load)
```

and add a divergence score to the n=16 branch-and-bound searches.  The goal is
to prove that every 15-speed primitive state vector has positive divergence.

That would turn the lonely runner proof at `n=16` into a no-harmonic-current
theorem on a finite 2-adic tree.
