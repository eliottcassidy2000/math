# Independent audit: cyclic lattices and full-group kernels

Status: **PASS**. Source manifest SHA256:
731636f7675a95316861a994130ec2545147c09594f8a16b4dfb4a181da2b084.
All five source files and all three pinned inherited sources match. No source
file was edited. Source normal and optimized replays pass 160423 checks and
byte-match their frozen receipts.

## Mathematical review

1. A genuine C_m x C_n action gives a nonempty transporter t+H. The inverse
   image of H in Z^2 is full rank because it contains both period vectors.
   Taking the least horizontal period h, least positive projected vertical
   coordinate v, and reduced lift s proves the stated HNF basis. Membership
   requires v dividing delta_b and h dividing delta_a-s(delta_b/v). The
   division-free second congruence has modulus hv, exactly as stated; replacing
   it with modulus h would generally erase information. Covolume hv, relation
   size mn/(hv), and both one-sided kernel counts follow.
2. The arbitrary-enumeration C3^2 hostile is valid. The identity-block
   transporter is a graph of the digit translation, and its first three
   values 1,2,0 cannot be an affine C9 graph. This rejects the encoding, not
   the underlying subgroup-coset description in its correct coordinates.
3. The global lattice count uses M=diag(m_i)Z^V contained in Lambda. Therefore
   |Lambda/M|=[Z^V:M]/[Z^V:Lambda]. Solvability by the augmented integer matrix
   and character constancy on a lattice basis are exact. Neither assertion
   remains a single coset after arbitrary non-coset filters; the source
   explicitly retains that boundary.
4. W/image(D) correctly classifies fixed-local-type attachment torsors under
   shared cyclic charts. The cardinality formula is orbit-stabilizer for the
   homomorphism D. It is not an unrestricted tournament-orbit formula.
5. The full-group result is precisely the subgroup projection/kernel
   construction of Goursat's lemma. The quotient relation is functional, but
   its original kernel choices remain present. After fixed source-target
   references, all simultaneous actual transporters are a coset of the
   intersection of source stabilizers. Path inventories or unspecified
   choices can instead produce unions; the source does not collapse them.
6. Visible one-sided kernels act independently in their own cells, embedding
   their product in the cell-preserving source automorphism group. Failure of
   any local owner set to be a base leaves a nonidentity global symmetry.
   This proves the sum-of-base-sizes lower bound, not sufficiency.
7. For R9=C3[C3,C3,C3], two distinct vertices have three common outneighbors
   exactly when they belong to the same inner triple; otherwise they have
   one. Thus the three modules are intrinsically identifiable. Every
   automorphism permutes them by an outer C3 rotation and rotates each inner
   C3 independently, giving exactly 81 automorphisms. Equality-of-outer-index
   blocks retain a one-sided C3^3 kernel of order27; the declared C9 subgroup
   retains its order3 translation kernel. Their base sizes are respectively
   three and one. For an isolated cell the group changes to the full internal
   group, but those respective base sizes remain three and one. The caveat in
   the source is correct and necessary.
8. Lifting congruences to bounded integer equalities is legitimate only with
   shared variables, actual guards, and finite bounds on their auxiliaries.
   The inherited conformally primitive-move criterion then applies to that
   compiled model. No claim of complete native-port representation follows.

## Independent controls

The independent checker does not reuse the source's block-shift or HNF code.
It enumerates all 247 subgroups of C_m x C_n for 1<=m,n<=6, derives each
lattice basis, and checks every offset, kernel size and covolume. It also
checks all 60 subgroups of S3 x S3 for projection counts, normal kernels and
quotient functionality. R9 module detection and both isolated-cell base
sizes are checked separately. Normal and optimized runs each pass 5981
explicit checks and have identical output.

The test programs are evidence for bounded controls; the displayed arguments
establish the all-size results within the declared chart models. A compact
complete native tournament observer, unrestricted electrical-guard linearity,
and efficient full-group computation are not claimed.
