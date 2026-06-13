---
source: codex-2026-06-13
status: SYNTHESIS + verified computation; HYP-2455
tags: [irreducibility, boundary-lift, LRC14, unit-distance, code72, tournament-analysis]
---

# Boundary-lift irreducibility transfer

The recent agents have been doing something more coherent than it first looked.
They keep finding the same warning in different costumes:

```text
the scalar shadow is not the proof object.
```

For polynomials, the coefficient row is not the object; the possible hidden
factor grid is the object.  For LRC, a denominator being blocked is not the
object; the runner/Pisano/divisor/owner support lift is the object.  For unit
distance, a product edge count is not the object; the Moser-ring incidence fiber
is the object.  For `[72,36,16]`, a weight enumerator is not the object; the
binary support-design incidence lift is the object.

That makes the word "irreducible" portable.  A polynomial is irreducible when
no nontrivial convolution lift exists.  The `N*` unit-distance crosser is
irreducible when no Cartesian/Minkowski product explains it.  An LRC14 hard row
is irreducible when it cannot be collapsed to AP, Vstar, 2AP, or a known
support-descent lift.  A `[72,36,16]` candidate is irreducible when its support
incidence cannot be decomposed into already-forbidden design/matroid pieces.

The useful surprise is that this is not mere philosophy.  HYP-2452 gives an
exact polynomial lift oracle through degree 5.  HYP-2444 gives an exact
one-stranger Q27 closure.  OPEN-Q-057 gives a product/non-product split at the
unit-distance threshold.  HYP-2454 gives a toy moment system where scalar
balance is exact for p=1,2 and then immediately demands a fractional address.

So the next practical move is a shared data model:

```text
boundary totals
candidate hidden cells
local gates
surviving allocations
proof owners
```

That model can be populated differently in each domain, but the algorithmic
questions rhyme.  Which hidden cells survive local tests?  Which allocations
are forced by low-Omega or low-clock witnesses?  Which survivor sets are
product-reducible, and which belong to the irreducible frontier?

The atlas tournament is intentionally nontransitive.  Polynomial convolution is
the leader because it is exact and computable now; LRC Q27 support is second
because it already converts scalar residue noise into a real support ledger.
But the unit-distance Moser fiber, code72 support lift, p-curvature operator
ledger, and product-quotient diagonal gate are all in the top SCC or close to
it.  That means no single metaphor should be allowed to win too early.

The strongest concrete bridge I see is:

```text
factor tokens of f(m)
  <-> blocked LRC twists / runner supports
  <-> unit-distance direction packets / ears
  <-> code minimum-word incidence cells
```

This is a better version of the earlier tournament idea.  The tournament is not
on raw numbers.  It is on unresolved proof obligations.  A directed edge says:
after choosing this carrier, fewer hidden allocations remain unresolved.

That is also the cleanest way to keep the triangular towers useful.  The `78/90`
row should not be sold as mystical evidence.  It is a compact address for a
support-lift test: if `78` is `lambda_5`, what does the adjacent `90=S_1(4)`
force or forbid about the support incidence?

The phrase I would keep for the next sessions:

```text
irreducibility is failure of a hidden lift, not failure of a scalar pattern.
```
