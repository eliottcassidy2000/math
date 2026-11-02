# Twenty objects: a small exact model and its missing information

Status: **PROVED finite relation isomorphism**, with a full explicit mapping
in twenty_controls.json. This is not a new square-sum path census or an
isomorphism of the underlying labelled graphs.

Let X be the20 directed Hamilton paths of the square-sum graph on1..25.
The [inherited joint path kernel](../joint_incidence_clocked_frames_20261102/PATH_KERNEL.md)
gives its compatibility graph4K2+12K1. Let R reverse a path, and let B
exchange the two compatible paths when a partner exists. R is a free
involution; it preserves compatibility and commutes with B on its8-point
domain. No compatible pair consists of a path and its reverse. Thus those
eight paths form two four-point orbits under R and B; the other twelve
form six reversal pairs. All these assertions are directly enumerated.

## A K4 reference model that preserves more than the count

Let Y be the20 three-edge subsets of K4's six labelled edges. Their
isomorphism types have multiplicities

    four stars + four triangles-with-isolate + twelve spanning paths.

Let C complement the three edges within K4, and let pi=(01)(23) permute
the four labels. Define R_Y=C pi. Define B_Y=C only on the star/triangle
sector. Then B_Y has four pairs and twelve isolated objects, and R_Y is
free and commutes with B_Y. On the active sector there are two four-point
orbits; on the inactive sector there are six two-point orbits.

Matching those orbits gives a bijection f:X->Y satisfying

    f R=R_Y f,   domain(B_Y)=f(domain(B)),   f B=B_Y f.

The supplied script chooses lexicographic representatives and records every
path and its K4 label. This preserves reversal AND compatibility, not just
twenty objects. The choice of pi and orbit representatives is a marking;
no canonical S4 action on the integer square-sum graph is asserted.

The original signed-arc incidence remains necessary for other consumers.
The checker reproduces the [prior rigidity certificate](../lower_frontier_complement_fibonacci_20261102/COMPLEMENT_CARRIER.md): starting with
degree colors, refinement by the sorted neighboring colors yields12,23,25
classes. Every automorphism preserves each refinement, so25 singleton classes
force the identity. Its abstract compatibility symmetries are therefore not
vertex relabellings.
The reference model does not compute actual endpoints or future larger
square-sum graphs without the stored path table.

For the specified square25 bank, each path cube has2^8 guard orientations
and each compatible intersection has2^6. Carrying these two weights gives

    2^32+20*2^8(z-1)+4*2^6(z-1)^2
      =4294962432+4608z+256z^2.

The remaining268 tournament edges are free. This recovers the already proved
polynomial2^268 times the display; the new contribution is the explicit
K4 model of its two compatible involutions. The inherited minimal33-state
complement observer can use the same labels, with its returning/closed
leaf distinction still retained. The twenty-object model alone is not that
closed transition system.

## A useful failed shortcut: rotate the six edge slots

The same Y is the weight-three slice of six bits, of size binom(6,3)=20.
Rotating the six slots has four necklaces with lengths2,6,6,6. In a chosen
normal basis of F64, that rotation is Frobenius, so this is an exact
necklace/field correspondence for that marked operation.

It does NOT preserve B_Y. In lexicographic edge order
(01,02,03,12,13,23), the star{01,02,03} rotates to{02,03,12}, a spanning
path with no B_Y partner. The rotation loses an enabled operation. Hence
Frobenius compression cannot replace the K4 carrier for this language.
This hostile identifies the missing coordinate: the incidence of the six
slots as edges of four vertices, rather than six exchangeable bit positions.

## Related15,24,26 counts, with their exact operation

Under Frobenius, the nonzero elements of F16 divide into one degree1
element, two degree2 elements, and twelve degree4 elements:

    15=1+2+3*4.

In F27, the prime field contributes two nonzero elements, and the other24
elements form eight degree3 orbits:

    26=2+8*3,  24=3^3-3.

These follow by counting fixed fields, and are elementary recovered field
facts, not novel results. In F8 the nonzero decomposition is7=1+2*3, whereas
1+2+4 denotes the separate multiplicative exponent orbit modulo7.
The [marked124 packet](../aliquot_power_124_20261026/tournaments/TOURNAMENTS_124.md)
already proves that one such orbit has linear rank2 and another rank3.
An orbit length therefore does not determine linear information or a normal
basis. This is the same defect exposed by rotating K4's edge slots.

The common mechanism is to name both the carrier and its allowed operation
before quotienting. Equal counts can suggest this model, but the commuting
involutions and the rotation hostile are what make it mathematically useful.

The [native26 test](NATIVE_26.md) now supplies a stronger consumer boundary:
two actual orientation sources in the same old complement state differ at
only6--10, but the identical10->26->23 attachment gives counts1 and2.
On the complete active universe the minimally refined terminal-extension
observer has38 states. This is an extension of the orientation consumer,
not a failure of the bijection on the twenty individual paths.
