# Clebsch and Bruhat are carriers, not decorations

The three web prompts resolved into a useful split.

The Clebsch graph is the folded 5-cube.  For LRC residual masks, that is exactly
what appears after choosing one tangent sector and folding the other five
sectors by complement.  The quotient has 16 classes, but every class mixes
missed depths.  So it cannot prove a scalar `q_t` inequality.  Its value is
different: closed neighborhoods form a pair-balanced `2-(16,6,2)` design, so it
is a covariance carrier.

That is the right object after S31l and HYP-2890.  The residual is signed; the
proof needs a Jensen/Schur-convex labelled functional.  Clebsch gives a finite
way to keep pair labels while quotienting the 64 residual masks.

The truncated octahedral graph is the Bruhat graph on `S4`: adjacent swaps of a
four-slot word.  It has 24 vertices, 36 edges, 6 commutation squares, and 8
braid hexagons.  Its edge orbits split into outer swaps and middle swaps.  This
is exactly the HYP-2889 lesson in graph form: local compression is not a
single unlabeled monotone operation.  Middle swaps carry different information.

The unital prompt contributed the discipline rather than a literal object:
choose a tangent, distinguish tangent/secant incidence, and only then
scalarize.  In the LRC proof this means:

```text
tangent sector -> Clebsch covariance class -> Bruhat compression face -> scalar comparison
```

The proposed next move is concrete.  Take the HYP-2890 residual leak, decompose
it over tangent Clebsch classes, then classify compression moves by Bruhat
square versus braid-hexagon faces.  If square faces are design-balanced and
nonpositive, the near-violations should concentrate in braid hexagons, which
are finite AP/Freiman packets.

No LRC14 proof is claimed, but this gives a sharper local shape for the
remaining residual-leak lemma.
