---
source: codex-2026-06-03-S590
status: synthesis + finite AP-clock orbit-functor audit
tags: [LRC, rigidity, orbits, functor, monodromy, CRT, tournaments, HYP-2134]
---

# Orbit Functors: Rigidity Is Which Action Preserves the Label

The local/global split was a good doorway, but it is too binary for what the
repo is showing.  A local fixed point is only useful when some action transports
it without losing the LRC predicate.  The right unit of thought is therefore:

```text
orbit functor + retained label + closure rule.
```

The AP witness set is the clean example.  HYP-2124 says the witnesses are
exactly the unit clock points `j/n`.  But the same set behaves differently
depending on the action:

```text
(Z/n)*       one static orbit, always
x -> 2x      internal only for odd n
j -> -j      antipodal pairs
CRT          prime-power channels
quotient     empty, rigid, leaking, or transporting fibers
```

So `n=14` is not simply rigid or nonrigid.  It is rigid under the full unit
action, rigidly paired under reflection, and broken under doubling because the
mod-2 coordinate collapses.  That is a much sharper sentence than "even cases
are hard."

## Vertex-Transitive Trienerment

The slogan

```text
vertex-transitive trienerment <=> regular polygon point-set
```

is now visibly a statement about one functor: cyclic sharp transitivity.  If
the functor is dihedral, the object can be an imprimitive bracelet.  If the
functor is nonabelian Cayley, a local root star can propagate through relators
without any order-`n` perimeter cycle.  If the functor is source deletion, the
rigid object is not a polygon at all; it is a marked source whose deletion is
collision-free.

This is the same correction as the LRC observer story.  The observer-blind
shape can look symmetric while the observer-coupled labels are doing all the
work.

## Types Of Rigidity Beyond Local/Global

Static rigidity: a finite witness orbit is pinned by a group action.  AP under
`(Z/n)*` is the model.

Dynamical rigidity: the proof move itself is an endomorphism of the witness
orbit.  Doubling has this property iff `n` is odd.

Reflective rigidity: a proof obligation is paired with its antipode.  This is
where pincer jaws naturally live.

Factor rigidity: CRT projections localize a defect to one prime-power block.
For `n=14`, the odd block is ordinary and the 2-block is the leak.

Quotient rigidity: a fiber map is useful only if fibers are empty, pure, or
have controlled transport rows.  HYP-1783 is the vocabulary for this.

Monodromy rigidity: a quotient loop lifts to labels.  Trivial monodromy means
the labels return unchanged; nontrivial monodromy is exactly how information
can be hidden while every local shadow looks legal.

Isostatic rigidity: the active obligations have exactly enough constraints.
This is a way to turn "worry-set" into a finite constraint count rather than
an aesthetic statement about symmetry.

Spectral rigidity: character blocks diagonalize orbit defects.  In the `2q`
case, the bet is that the non-principal leak is concentrated in the 2-adic
block.

Fiber stiffness: a class can be locally hard because boundary-label flips cost
more than inner flips.  The old antiferromagnet/fiber-bundle notes suggest
that ensemble isotropy can hide class-level anisotropy.

Source rigidity: adding or deleting a marked source is exact.  This is the
cleanest tournament analogue of an LRC observer-source certificate.

## Multiplication, Addition, Odd, Even

Multiplication and addition are different rigidity languages.

Multiplication by a unit moves along the AP witness clock without changing the
row.  Multiplication by `2` is a legal witness-dynamics move only when `2` is a
unit.  That is why odd/even matters: at even `n`, doubling stops being an
automorphism and becomes a collapse.

Addition creates folds and denominators.  A pair `a+b=D` defines a pinch clock
and a shield condition `D|v`.  This is observer-coupled data.  It cannot be
replaced by balanced additive energy, because balanced energy forgets the
observer and therefore misses the hard fold labels.

The synthesis is:

```text
odd/multiplicative = orbit propagation
even/doubling      = projection defect
additive/fold      = denominator pincer
CRT                = where the projection defect lives
```

## What To Do Next

The next useful checker should take a residual object and report:

```text
1. witness residue orbit under the available actions,
2. reflection pincer pairs,
3. CRT block where each action fails or stays internal,
4. quotient-fiber purity/leakage,
5. monodromy of retained labels around quotient loops,
6. boundary/inner stiffness ratio for observer-coupled labels.
```

For `n=14`, that means starting from the six unit witnesses, pairing them into
three antipodal jaws, localizing the only dynamic break to mod 2, and refusing
to trust any unmarked quotient that does not carry source, denominator,
endpoint, and pincer labels.

## Anchor

`04-computation/orbit_rigidity_functor_atlas_s590.py` (+ `.out`);
HYP-2134.
