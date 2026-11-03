# A marked counting refinement of the inherited fixed-H core series

**PROVED counting identities / FINITE-EXACT atom and finite checks.**
The ordinary fixed-H rational series is inherited from
[Session24 multiplicities](../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md)
and [Session153 species](../session153_gadgets_20261103/species/SPECIES_TROPICAL.md).
This packet spells out the one-port and ordered-distinct-two-port refinement
and independently computes its small atoms. No external priority claim.

Fix an ordered word W of nontrivial strong tournament classes. Let its
length be k, total order m, and Hamiltonian count h. For its i-th component,
let r_i count vertex orbits under its automorphism group and q_i count
orbits of ordered DISTINCT vertex pairs. Set

    R=sum r_i,   Q=sum q_i + 2 sum_(i<j) r_i r_j.

All automorphisms preserve the ordered strong components. The full group
is the product of their groups. Thus R and Q are respectively the one-port
and two-port orbit counts of the core. Inserting j singleton components
gives binom(j+k,k) ordinary classes, each with

    one-port orbits = R+j,
    ordered-distinct-two-port orbits = Q+2Rj+j(j-1).

These formulas retain the same-component joint data q_i; it cannot be
reconstructed as r_i(r_i-1). For the directed triangle, r=1 but q=2.

Consequently this core contributes the following three ordinary generating
functions, with z recording total vertex order:

    F_W = z^m/(1-z)^(k+1),
    R_W = z^m [R/(1-z)^(k+1) + (k+1)z/(1-z)^(k+2)],
    Q_W = z^m [Q/(1-z)^(k+1) + 2R(k+1)z/(1-z)^(k+2)
               + (k+1)(k+2)z^2/(1-z)^(k+3)].

Sum over actual ordered core words of product h. For h<=21 this sum is
complete using the fourteen strong atoms independently verified in
`audit.json`: one each at H3,H5,H11,H13, two each at H9,H17,H19, and four
at H15 (two of order5, two of order6). The empty core gives h=1, including
the order0 unmarked tournament. Negative singleton counts contribute zero.

For example, for every n>=3 the H3 fibre consists of exactly n-2 ordinary
classes, and has

    (n-2)^2 one-port classes,
    (n-2)[(n-3)(n-2)+2] ordered-distinct-two-port classes.

The checker compares all three channels against EVERY ordinary class
through order7. Their all-order validity follows from component uniqueness,
not from extrapolating the finite comparison. They count marked graph
classes; they are not a smallest observer for arbitrary future attachments.
In particular, equal H, R and Q need not identify a literal source, clock
or port response. The prior81/27 and55->111 controls remain relevant to
that stronger claim and are not replaced by these generating functions.
