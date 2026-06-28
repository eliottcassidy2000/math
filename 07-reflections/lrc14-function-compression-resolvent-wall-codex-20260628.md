# LRC14 Function-Compression Resolvent Wall Reflection

Source: codex-2026-06-28  
Anchor: HYP-3150 / T1215 / LTI-276 / LTT-174

## One Sentence

The useful abstraction is not "small tournaments look like groups"; it is:
every proof move is a function, and every compression must prove that the next
observable factors through it or pay for the lost coordinate as a sidecar.

## Map 1: Compression Tower

```text
ordered pair (a,b)
  q = forget order
  survives: a+b, a*b
  lost: a^b, b^a unless orientation sidecar is kept

K3 edge cube
  q = score class C/T
  survives: edge-flip Markov kernel, stationary split, -1/3 mode
  lost: minority edge, source-sink order, Worpitzky descent, state PGF curve

K4 fixed-path cube
  q = x=a OR c, y=b OR c
  survives: T,+,-,S class
  lost: flip word, S fiber PGF z+3z^2+z^3, c-canary status, deletion stability

k=8 resolvent
  q = u -> v=u^2
  survives: u^4-5u^2+4, quadratic v^2-5v+4
  lost: sign/odd coordinate, reflection-resurrection data, boundary leakage
```

The tower suggests a proof style: compress only after writing the precise
observable being pushed forward.  The same quotient can be legal for class and
illegal for PGF, legal for value and illegal for derivative, legal for even
coordinate and illegal for odd sidecar.

## Map 2: Sidecar Periodic Table

```text
orientation sidecar:
  ordered endpoints, a^b versus b^a, tip/tail commutator

role sidecar:
  K3 minority edge, exit edge, edge-flip gate

curve sidecar:
  state-level PGF/root curve, Lee-Yang locus, log derivative

fiber sidecar:
  S fiber mass, weight PGF, lower-order leakage

canary sidecar:
  fixed filler coordinate, collapse slice, exact transversal

deletion sidecar:
  live core, deletable coordinate, delete-one-stable representative

odd sidecar:
  sign after u -> u^2, antisymmetric block, boundary leakage

degree sidecar:
  effective degree, algebraic wall, named reason degree <=4 is enough
```

The recurrent structure is that "single values" are almost always too cheap.
The proof wants typed pushforwards: value, curve, root, fiber, derivative,
orientation, and deletion status are different functions.

## Creative Hypotheses

1. The HYP-3140 `Rprime` fiber-PGF theorem should be recast as a
   factor-through theorem for the next observable, not as a scalar lower
   bound.  Ask which coefficient curve survives the quotient.
2. The k=8 U4 sidecar is a quartic only in raw form.  Its proof-facing object
   is a two-sheet square map over a quadratic base, plus odd-coordinate debt.
3. The K4 `S` fiber is the baby version of bounded-core slack: high
   representation multiplicity is useful only when the canary/deletion
   coordinate says why it cannot be deleted.
4. The K3 `-1/3` mode is a local signature of legal compression with one
   missing oscillatory coordinate.  Search for `-1/3`-like contractions in
   HYP-3129 signed SPEC low modes.
5. Worpitzky data is best viewed as a fiber derivative.  It measures how the
   transitive fiber splits by order after the class quotient has already
   compressed too much.
6. A future genuine quintic object should trigger a route alarm.  The current
   LRC14 evidence says the hard branch is solvable because every live
   compression either has degree <=4 or has an even/filler/fiber sidecar that
   drops the effective degree.

## Mainline Cross-Signal

Incoming S71 sharpened the boundary after this scout was reserved:

```text
score sequence determines tournament iso class for n=3,4
score sequence fails to determine iso class at n=5 (#iso=12, #score=9)
```

That is the same as saying the commutative face works through the K4 table and
breaks precisely when order-sensitive data becomes unavoidable.  S71 then
matches the break to the LRC cap dip:

```text
|P|<=3: dip=0
|P|=4: tiny dip
|P|=5: large k=8 binding dip
```

The parity split is even more suggestive:

```text
even +6S4 = symmetric / biquadratic / sum-product
odd  -9S3 = antisymmetric / Worpitzky / ordered orientation
|odd|/|even| = 3.15
```

So the proof route should not over-credit the even biquadratic side.  The odd
Worpitzky/orientation side is the larger local correction and probably the
place where the missing endpoint/order sidecar has to be paid.

Incoming kps-S31ah adds the root-curve version: the AP/consecutive rows are
most circular for k=8,9 in the tested bank, and off-circle variance lambda
correlates positively with the coverage gap.  That makes "whole PGF curve"
literal: root circularity is not a decorative signal, it is one of the
coordinates that prevents illegal scalar compression.

## Working Rule

Before a quotient enters the LRC14 proof packet, add:

```text
compression_map
observable_factors_through
fiber_collision_class
ordered_sidecar_required
fiber_pgf_curve_status
canary_deletion_status
even_resolvent_variable
effective_degree
terminal_exit_or_named_debt
```

This rule is stricter than the metaphor and more useful than the metaphor.
It makes "below Abel-Ruffini" a checkable finite-packet discipline rather than
a slogan.
