# Reciprocal Laurent transport as boundary flux

Date: 2026-09-03

Status: **STRUCTURAL SYNTHESIS FROM THM-4377; ONE PRECISE NEXT TEST.** This
reflection is non-canonical. It does not prove bracket factorization, packet
entry, a seam theorem, `JC(2)`, or `DC(2)`.

## Trigger

THM-4375 showed that reciprocal boundary types can have unequal numbers of
polynomial source monomials. The tempting conclusion was that reciprocity had
simply disappeared at source level. THM-4377 retains the full exponent tuple
and shows something sharper: reciprocity remains a deterministic translation
on the ambient Laurent lattice, while polynomiality cuts that translation at
three coordinate walls.

The distinction is useful because it turns an unexplained dimension mismatch
into an exact boundary term.

## The transport picture

For a bilateral type `(w+d,w)<->(w,w+d)`, the map is

```text
T_d(a,b,c,e)=(a+2d,b-2d,c+d,e).
```

It preserves the packet address and the `e` grade. On polynomial exponent
tuples it becomes a partial bijection. Its entire failure is supported on

```text
source loss:   b<2d,
target loss:   a<2d or c<d.
```

This is formally analogous to a deterministic Markov transport stopped at an
absorbing boundary. The analogy is deliberately limited: there are no
probabilities here, and arbitrary source coefficients are not transported
unless an additional coefficient sidecar is supplied. What survives exactly
is the basis-level coupling and its boundary flux.

At `w>=ceil(ell/2)` no target mass is lost. The free source module therefore
splits as

```text
larger fibre = transported reciprocal fibre + b-boundary jet.
```

At `w>=ceil(ell/2)+d-1`, that jet has the fixed coordinate inventory
`b=0,...,2d-1`. This is the stable object that raw fibre cardinalities were
hiding.

## Connection contract

```text
source:      polynomial source-monomial basis of a bilateral packet
target:      the reciprocal packet basis plus boundary defects
map:         partial Laurent translation T_d
preserved:   packet pair, e grade, integral basis elements on the matched core
destroyed:   polynomial nonnegativity at a, b, and c walls
sidecar:     wall label and, in the stable range, the finite b-jet coordinate
hostile:     ell=3,w=1,d=1 has equal fibre sizes but no matched basis element
```

## Comparison with the live concept board

1. **Pascal circuit quotient (THM-4369).** The circuit theorem quotients packet
   types after applying the trace map. Laurent transport lives one level
   earlier, on individual source monomials. The quotient forgets precisely the
   wall label needed here.

2. **Reciprocal packet orbits (THM-4375).** Orbit geometry supplies the
   `(w,d)` coordinates; the mutation supplies the missing lift and explains
   every imbalance.

3. **Reciprocal eigenlattice (THM-4378 candidate).** The signed involution is
   lawful after quotienting by Pascal circuits, whereas the source lift is
   only partial before quotienting. Their discrepancy is a concrete place to
   search for the sidecar carried through the quotient.

4. **Source-normal depth hierarchy (THM-4376 and THM-4380 candidate).** The
   depth hierarchy can be bracket-blind. Boundary jets offer a finite source
   coordinate on which a missing bracket functional might concentrate, but no
   such factorization is presently known.

5. **Exceptional conductor fibres (THM-4034 and THM-4381 candidate).** Both
   mechanisms repair a lossy quotient by retaining a finite boundary object:
   conductor fibre labels there, coordinate-wall labels here. The common move
   is structural; the algebraic objects and consequences are different.

## The next decisive test

Let `Br` denote the relevant bracket residual map on a stable bilateral packet
prefix, once its source coefficients have been typed. In the monomial splitting

```text
Z[E_+] = T_(-d) Z[E_-] direct_sum Z[B_+],
```

compute the two blocks of `Br` separately. The cheapest useful question is
not whether reciprocity is a bracket symmetry globally. It is:

```text
Does Br(T_(-d)y) reduce to a signed reciprocal transform of Br(y),
leaving a residual map supported only on Z[B_+]?
```

A positive answer would reduce the stable asymmetric packet calculation to a
`2d`-dimensional boundary operator. A negative answer should be recorded by
the first basis vector in the transported core where intertwining fails and by
the missing coefficient coordinate. Either outcome is structurally useful.

The first hostile should be the smallest stable asymmetric control
`ell=10,w=8,d=4`, followed by the transient mixed-wall control
`ell=10,w=4,d=1`. Passing only the stable case is not enough: the transient
case distinguishes genuine wall localization from accidental large-scale
cancellation.

## Stopping certificate

The basis transport, wall decomposition, and stable jet are exact in THM-4377.
No current theorem identifies the bracket residual with an operator that
intertwines `T_d`. Until that operator is written on the same monomial basis,
claims that the boundary jet carries the bracket obstruction are hypotheses,
not consequences.
