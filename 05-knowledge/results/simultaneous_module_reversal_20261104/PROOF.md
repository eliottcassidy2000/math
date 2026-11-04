# Simultaneous disjoint module reversals cannot cancel

**PROVED S1-S3; same-session logical review and independent exact algorithms.**
Named statements S1-S3 below concern pairwise
disjoint modules in one fixed quotient. Nested or overlapping reversal
families require the separate scope discussion at the end.

## Definitions and inherited inputs

Let Q be a tournament and let B_i be pairwise disjoint, nonempty named
vertex sets. Replace its vertex i by a tournament M_i on B_i, obtaining
T=Q[M_1,...,M_k]. Unselected exterior vertices can be singleton blocks.
Let I be a set of indices with |B_i|>=2. For epsilon in {0,1}^I set

    T_epsilon=Q[M_1^(epsilon_1),...,M_k^(epsilon_k)],
    M_i^0=M_i, M_i^1=M_i^op.

All names and all arcs between different blocks are fixed. A port query
adds x->y with a_u=[u->x],c_u=[u->y]. No coordinates are sorted or
identified under a vertex permutation. The full response F_T has bidegree
at most(2,2), hence total degree at most four. Its arbitrary difference
may have nonzero cubic **and quartic** terms.

The closest proved mechanism is **E3/E4**,
[all-order single-converse separation](../two_port_parity_separation_20261104/PROOF.md),
which gives an actual core query separating M_i and M_i^op by2 modulo4.
The other input is the complete ordered path-cover profile

    q_j(M)=number of ordered spanning covers by j nonempty directed paths.

Reversing each path preserves its ordered slot, so q_j(M)=q_j(M^op)
exactly. The occupied-run expansion, most recently recovered in
[S188's skeleton transform](../session188_structural_bridges_20261104/skeleton/PROOF.md),
states H(R[M])=sum_j b_j(R,z)q_j(M), with integer nonnegative coefficients.
Its arbitrary-uniform-port consequence is already **G2** in
[S186 gate](../session186_proof_frontiers_20261104/gate/PROOF.md).

## S1. Exact isolation of one module in a reversal cube

Fix i in I and any query on B_i. Lift it to the whole named tournament
by giving each port the root i's incidences to every exterior vertex:
if i->j in Q, both ports beat every vertex of B_j; if j->i, every vertex
of B_j beats both ports. These masks depend on Q and i, not epsilon.

**Theorem S1.** The actual integer value F_(T_epsilon)(lift_i(a,c))
depends on epsilon only through epsilon_i.

**Proof.** For j!=i, the entire block B_j remains a module after both
ports are adjoined: all its other exterior incidences are uniform by
substitution, and the ports see it uniformly by construction. Apply the
occupied-run expansion to that actual extended tournament at B_j.
Reversing M_j preserves every q_l, so it leaves this actual Hamilton
count unchanged over the integers. This argument is valid for every
choice of the other block orientations. Reverse them one at a time
until they are all in their original orientation. Only epsilon_i remains.
All operations concern disjoint blocks, so their incidences and module
status remain valid throughout. QED.

The equality is exact, not merely modulo four. It does not assert
additivity of response functions at arbitrary attachments, and does not
discard any interaction or quartic coefficient there.

## S2. Every reversal bit has a fixed named query coordinate

For each i choose, by E3, a core query (a_i,c_i) satisfying

    F_(M_i)(a_i,c_i)-F_(M_i^op)(a_i,c_i)=2 modulo4.

Let q_i be its lift from S1, and set

    A_i=F_(T_0)(q_i),
    d_i=F_(T_(unit_i))(q_i)-A_i.

**Theorem S2.** For every epsilon and every i,

    F_(T_epsilon)(q_i)=A_i+epsilon_i*d_i,
    d_i=2 modulo4.                                      (1)

**Proof.** The exact affine two-value formula is S1. For its nonzero
residue, contract B_i while retaining all other original blocks, giving
an actual quotient R_i. The port-extended B_i is a module, so **T1**,
[ambient residue transfer](../two_port_twelve_principles_20261103/PROOF.md), gives

    d_i=H(R_i)(F_(M_i^op)(a_i,c_i)-F_(M_i)(a_i,c_i)) mod4.

H(R_i) is odd, so d_i=2 modulo4. No denominator, response normalization,
or cancellation of R_i over the integers is involved. QED.

Thus the fixed |I| actual queries q_i recover the entire reversal vector:

    epsilon_i=(F_(T_epsilon)(q_i)-A_i)/2 modulo2.

All2^|I| named response functions are distinct. In particular **no
nonempty simultaneous set of disjoint module reversals can cancel** in
the full response, even modulo four. This holds at every module size,
number of modules and ambient order, for strong or nonstrong quotients.
It excludes collisions between any two reversal vectors, not only
between a vector and the all-zero state.

The query count is a constructive upper bound within this known family;
no optimality claim is made for general integer-valued measurements.

## S3. Profile-preserving changes have the same isolation property

S1 also holds if each B_j can be filled from any prescribed family of
same-named tournaments with the same complete q profile. For two global
states, choose a differing block i whose core responses differ modulo4.
A separating core query lifts to a global separator: all other blocks
can be replaced without changing its integer value, and T1 preserves
the residue at i. Conversely a block with identical full core response
is not declared distinguishable by this theorem.

This is a restricted product-injectivity theorem. Equality of H alone
does not meet its hypotheses; all q_j are required to remove exterior
changes from the actual extended count.

## Scope boundary: nested and overlapping operations

Pairwise disjointness is part of S1-S3. An outer module reversal together
with reversals of all its children changes the quotient orientation
while restoring the children's internal orientations. The ports cannot
in general be made nonuniform on the outer module and simultaneously
uniform on every changed child. S1 is therefore not a proof for arbitrary
nested families or for arbitrary legal sequences of module reversals.

For an arbitrary family of modules of the *original* tournament, the
phrase simultaneous also needs an edge convention. The natural XOR
convention reverses an arc when an odd number of selected modules contain
its two endpoints. For overlapping modules this need not be a sequence
of legal module reversals: in TT3, {0,1} and {1,2} are modules initially,
but flipping both yields a directed triangle and H changes from1 to3.

Identical final graphs from redundant operations are tautological equality,
not an exact collision of distinct named tournaments. The broad nested/
overlapping question is not declared closed by S1-S3. Retained complete
two-level controls and their uniform-probe hostile are described in the
certificate and review. No exact full-response collision is claimed.

## Exact controls: interactions exist, but disappear on the decoder

**FINITE-EXACT.** Take Q=C3 with rows (2,4,1), and three TT2 blocks
on names {0,1},{2,3},{4,5}, respectively. The base tournament has rows

    T=(14,12,56,48,35,3).

Write F_e for the response when the blocks indicated by the binary mask e
are reversed, with Q fixed. The complete mixed interaction is nonzero:

    [a_0 a_2 c_1 c_4](F_0-F_1-F_2+F_3)=-2.             (2)

Thus even modulo four, adding single-flip differences does not give the
full simultaneous difference. The entire 484-coordinate interaction is
retained, with 128 nonzero coordinates. A literal sixteen-assignment finite
difference of four actual Hamilton counts independently verifies (2).

Nevertheless the three decoder queries (a,c)=(48,49),(3,7),(12,28)
have value393 if their selected block is unchanged and483 if it is
reversed, independent of the other two bits. Their difference90 is2
modulo4. These are integer bit masks, with physical name u at bit2^u.

For the nested boundary, reverse Q alone, restoring all child interiors:

    U=(50,48,11,3,44,12).

This is the result of reversing the whole T and all three child pairs.
All64 attachments uniform on the three pairs have F_T=F_U modulo4.
Exactly24 of those assignments differ over the integers, so this is not
an exact uniform-response collision. A nonuniform actual query gives

    F_T(a=1,c=5)=251,     F_U(a=1,c=5)=249.

The full difference also has nonzero quartics. This is a hostile to
uniform-modulo-four quotient contraction, not to S1, and not a full-response
collision. Both quotient orientations and all child reversals yield16
distinct complete responses in this example. These finite facts do not
establish the unrestricted nested theorem.
