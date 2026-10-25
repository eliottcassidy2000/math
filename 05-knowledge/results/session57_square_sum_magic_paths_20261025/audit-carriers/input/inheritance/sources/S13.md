# Two-sheet repair: converse is a polarity, not a point collineation

Status: PROVED transported incidence structure, with exact finite verification. This is a NEW relation on two typed copies of the boundary quotient and blocker set. It is not claimed to be native arc-intersection incidence.

## Base object and relation

Let K=F64, F=F4, and Tr(z)=z+z^4+z^16. Write boundary quotient labels e modulo21 using README.md's adjusted exponent map. Take two sheets:

- Point sheet e represents the one-dimensional F-subspace [g^e].
- Line sheet l represents L_(g^-l)={x:Tr(g^-l*x)=0}.

Declare point e incident with line l exactly when

    Tr(g^(e-l))=0.

Changing either exponent by21 multiplies the trace argument by a nonzero F4 scalar, so incidence is well-defined. The trace pairing (a,x)->Tr(ax) is nondegenerate: for a nonzero choose x=a^-1, giving Tr(1)=1. Its21 coefficient classes therefore parametrize all21 projective lines. There are42 typed objects,105 incidences, and degree5 at every object.

## Genuine collineations and converse polarity

Each G21 map e->16^k*e+9b acts on both sheets by the same label formula. It is induced on K by the F4-linear map x->g^(9b)*x^(16^k). On a line coefficient g^-l the corresponding transformed coefficient is g^(-16^k*l-9b), giving exactly the asserted line-label action. Thus these really are collineations, not just arbitrary point permutations.

Converse acts by

    (e,point) -> (-e,line),
    (l,line)  -> (-l,point).

The new incidence test is unchanged, since (-l)-(-e)=e-l. This is precisely the trace polarity [x]->{y:Tr(xy)=0}; it squares to the identity, reverses the two sheets, conjugates multiplication by g^9 to its inverse, and commutes with the16th-power Frobenius. The five absolute points are the trace-zero projective points. The exact enumeration finds labels {3,6,7,12,14}.

There is no conflict with the earlier negative result: point-only inversion e->-e sends no projective line to a projective line. The typed map switches points and lines; that additional type coordinate is the essential repair.

## Lifting through the actual native boundary

On two typed copies of the63 blockers, define

    (B,point) incident with (C,line)
    iff Tr(Psi(B)/Psi(C))=0.

Because each native-boundary fiber is an F4× coset, this relation depends only on partial(B) and partial(C). Every quotient point and line has exactly three blockers over it. Each of the105 base incidences lifts to a K3,3, producing126 vertices,945 bipartite edges, and degree15. This is a threefold vertex blow-up, not a graph covering in the locally bijective sense. The checker verifies that the equal-neighborhood classes in each sheet are exactly the original three-element boundary fibers.

The actual blocker converse B->iota(B), together with a sheet switch, realizes the polarity because Psi(iota(B))=Psi(B)^-1. Actual G21 acts within each sheet. Thus the larger typed universe preserves a genuine relation while retaining both the exact blocker objects and their actual symmetries.

## Native-relation hostile

This newly transported relation cannot be inferred merely from the number of deleted arcs two blockers share. In the checker bit convention, blocker pairs(63,28679) and(63,455) both share three arcs, but only the first pair is incident in the new point/line relation. The full intersection-count contingency table is included in two-sheet.json. No equality with native arc-intersection, HP counts, or native XOR is asserted.

## Verification

`two_sheet.py` checks all441 base pairs under all21 affine maps and polarity, the105 flags and regular degrees, all3969 lifted pairs, all945 edges, twin/fiber equivalence, absolute points, and the explicit native-intersection hostile. Normal and optimized outputs are identical and included in portable replay.
