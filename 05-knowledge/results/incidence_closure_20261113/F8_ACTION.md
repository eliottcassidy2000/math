# The order-21 core is a semilinear F8 action, with a lost phase

**Status: PROVED IN NOTES after independent review; see [review](audit/F8_REVIEW.md).**
The Cayley tournament and its non-self-converse certificate come from incoming
Session319. The group identification below is elementary finite-field theory;
no external priority claim or implication for abc, F1, or the 7/21 Hamiltonian
exclusions is asserted.

Let F8=F2[theta]/(theta^3+theta+1). The cubic has no F2 root and is therefore
irreducible. Its nonzero elements form a group of prime order7, and theta is
neither zero nor one, so theta generates that group. Frobenius x->x^2 has
order3: its orbit of theta is theta,theta^2,theta^4, the literal 1-2-4 exponent
cycle.

For a modulo7 and b modulo3 define the invertible F2-linear map

    L_(a,b)(x)=theta^a x^(2^b).

Composition gives

    L_(a,b) L_(c,d) = L_(a+2^b c, b+d),

which is exactly the group law of the 21-vertex core. The maps are all distinct:
evaluation at1 determines a, and then evaluation at theta determines b. Thus
these are the 21 semilinear transformations GammaL_1(F8), realized faithfully
as F2-linear maps on the three-dimensional field. The Cayley tournament is on
these **operators**, not on the seven nonzero field points. Multiplication by
an element of its specified connection set gives the intrinsic orientation.

## The natural seven-point quotient loses an actual orientation

The map L->L(1)=theta^a has fibres of size3; it forgets b. It is the natural
transitive action on nonzero field elements, with Frobenius as the stabilizer
of1. However these fibres are not modules of the actual Cayley tournament,
so the tournament does not descend through this map.

Using the incoming connection set S, take

    u=(0,0), u'=(0,1), v=(1,1).

Both u and u' evaluate to1. Since (1,1) belongs to S, u->v. But
u'^(-1)v=(4,0) does not belong to S, so v->u'. More completely, the nine
arcs between the evaluation fibres a=0 and a=1 split eight in one direction
and one in the other. There is no uniform orientation between these two
quotient points. Retaining the ternary Frobenius phase b recovers every
operator, its group multiplication and the original tournament incidence.
No tie is artificially turned into a tournament edge.

This is an exact source/target map and an exact obstruction: a group action
quotient need not be a tournament quotient. The phase lost by evaluation is
analogous in role to the sheet signs lost by merging, but the groups, carriers
and reconstruction formulas differ. The current signed operator D and this
Frobenius phase are not identified with one another.

In the [joint-shadow family](JOINT_SHADOW_BOUNDARY.md), using U and U^op as
ordered blocks loses independent orientations even after retaining Q, weighted
arc usage and endpoints. Here F8 explains the actual symmetry mechanism that
makes those endpoints flat. The number21 is the order/vertex count of the
core, not a claim that H(U)=21; transitivity only implies 21 divides H(U).
Thus this bridge respects the forbidden-H21 boundary rather than confusing
it with a group order. It supplies no prime-infinitude or abc implication.

[The finite checker](audit/f8_action.py) verifies all21 maps on all8 field
elements, every pairwise composition, F2-linearity and invertibility, every
Cayley orientation and the explicit quotient hostile. Normal and optimized
receipts agree. The displayed algebra proves the identification; the finite
checks independently certify this concrete representation.
