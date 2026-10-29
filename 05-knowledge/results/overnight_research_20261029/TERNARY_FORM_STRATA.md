# All native ternary forms, and the necessary extra symplectic plane

2026-10-29. **FINITE-EXACT complete form census; PROVED necessary completion
reduction.** Same-agent verification using independent integer/local Smith
calculations. The index3^10 attainment question remains OPEN.

The source is the saved nine-coordinate Paley-seven three-run module,
V=C3^2+C9^6+C81, with every original port projection and all126 native
clocks retained. The nearest mechanism is
[TERNARY_EXTENSION_OBSTRUCTION T1--T3](TERNARY_EXTENSION_OBSTRUCTION.md).
The hostile remains81/27; the corrected near miss is assuming that the
one saved optimal form represents every possible restriction. The
least-used sidecar recovered here is the COMPLETE invariant-form module,
including its nonunit coefficients.

## S1. The complete form module is C3^4+C9

Write each alternating form with denominator81. Its upper-triangle entry
(i,j) is divisible by3^(4-min(e_i,e_j)). The36 uniform variables over
Z/81 have a representation kernel of order3^87. The full native clock
invariance equations have kernel order3^93. Their quotient, the group
of actual invariant alternating forms, therefore has order3^6=729.

The checker reruns the full local elimination, then extracts five explicit
independent form generators of orders3,3,3,3,9. All729 combinations are
different and exhaust the invariant-form group. For each form, two
different elimination paths compute the radical: integer Smith factors
of the map to the actual dual, and local unit-pivot/division layers.

| Radical order | Number of invariant forms |
|---:|---:|
|3^8|324|
|3^12|324|
|3^14|70|
|3^16|10|
|3^18|1|

In particular, **no invariant source form has radical order3^10**.
This closes the second branch listed in the previous handoff. Every
index3^10 completion must start with one of the324 forms of radical3^8.
Their radicals all equal the previously proved common radical R, since
they contain R and have the same order. The first-step obstruction is
unchanged.

## S2. Twenty sufficient source-form representatives

Solve the full heterogeneous endomorphism equations A C=C A for all
native clocks C. The image of this endomorphism ring on V/3V has
dimension6; all729 residue endomorphisms and all324 invertible residues
are enumerated and lifted. The checker also retains explicit invertible
maps1+tE from the integral endomorphism generators. An endomorphism is
invertible exactly when its action on V/3V is invertible.

Their certified subgroup acts on forms by b->A^T b A. Its optimal-form
orbits have sizes

    four orbits of size9, and sixteen orbits of size18.

Thus **20 supplied forms suffice for a complete existence search**.
We do not claim that these are20 inequivalent classes under the full
clock centralizer: an ungenerated automorphism could merge them. The
receipt supplies every chosen matrix and a spanning tree taking each
of the324 optimal forms to its representative.

The reduction retains the original labels constructively. If
b'=A^T b A and i embeds(V,b), then iA embeds(V,b'). For the original
210-port projection P, the transported port images are iAP. Products
along the saved spanning tree recover the complete map. These source
module automorphisms are not asserted to be native vertex re-markings.
This is an optimization of the abstract embedding problem, with the
port maps carried explicitly, not permission to forget them in a native
attachment or arithmetic operation language.

## S3. Every index3^10 completion has one of six defect-plane actions

Suppose a perfect paired W contains V equivariantly with index3^10.
Its restriction has radical R of order3^8 by S1. Put K=V-perp. Then

    |K|=3^10,   rad(K)=K intersect V=R,
    P=K/R is a perfect alternating group of order9.

Consequently P is an elementary symplectic plane F3^2. Furthermore

    L=V+K=R-perp,
    L/V is P,
    L/R=(V/R) orthogonal-sum P.                       (S1)

All statements are equivariant. The quotient P is the two-dimensional
paired defect that any index3^10 candidate must supply.

Let F be the native vertex order-three clock, R_run the run cycle, and
S the run swap on P. Since an involution in SL2(F3) is either I or-I,
S R_run S=R_run^-1 forces R_run=I. The remaining F can be the identity
or a nontrivial order-three unipotent. Its two nontrivial symplectic
conjugacy classes are represented by

    [[1,1],[0,1]] and [[1,2],[0,1]].

Together with S=+I or-I this gives exactly **six** plane-action types.
The checker enumerates all24 elements of SL2(F3), all18 literal triples
satisfying the group relations, and their six conjugacy classes of sizes
1,1,4,4,4,4. The shear with S=diag(1,-1) is a hostile: it satisfies the
S3 group relation but S reverses the alternating form, so it is not
an admissible symplectic plane action.

The scalar-involution argument holds for every odd characteristic.
Characteristic two has a different involution boundary; no import from
the fifteen F2 Lagrangians is being made.

## S4. The reduced next test

Choose one of20 source forms and one of six necessary plane actions:
120 sufficient starting cases, before their extension carries are chosen.
First construct the index-nine intermediate L in S1, retaining source
injection, pairing and every clock. Then seek a coisotropic perfect
completion of L of index3^8. Equivalently, a length-ten root chain from V
has one radical-increasing and nine decreasing steps. Its first two steps
inside L have the same S sign and elementary quotient C3^2; their joint
run-cycle action on that quotient is trivial. The vertex clock may have
either of the two unipotent classes or be trivial.

These necessary constraints do not construct L or W. In particular,
enumerating the six quotient actions does not solve their lifting and
commuting-carry equations. Neither group order nor a separately extendible
clock supplies those equations for free.

This refines the paired-group/clock bridge, not an arithmetic Selmer map.
The full native minimum still has exponent a in{10,12,14,16,18} at3,
and the already proved dyadic minimum remains2^40. Sha finiteness,
Collatz, and LRC(14) remain open.

Reproduce with `python -B ternary_form_strata.py` and its -O variant.
All729 radicals, clock-isometry edges, plane actions, and81/27 controls
are checked with explicit exceptions. Normal/optimized receipts agree.
