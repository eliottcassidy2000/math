# Independent audit of the integral mutation and triangle interface

**VERIFIED mathematical audit / FINITE-EXACT independent controls.**
The audited files are [matrix/PROOF.md](../matrix/PROOF.md) and
[matrix/check.py](../matrix/check.py). No mathematical correction is
required. This auditor imports neither matrix file as executable code.

## The congruence has the correct signs and domain

Writing x=b_ik, y=b_kj, the off-pivot entry of E B E^t is

    b_ij + [x]_+ y + x[-y]_+.

The extra term is xy if x,y are positive, -xy if both are negative,
and zero if their signs differ; zero entries also agree. These are exactly
the four cases of the ordinary Fomin–Zelevinsky formula. Entries incident
to k reverse sign because E_kk=-1, and b_kk=0 prevents a quadratic extra
term. This works for arbitrary integral skew matrices, not only tournaments.

Only column k differs from the identity. Its diagonal coefficient is -1,
so det E=-1; applying the same matrix twice cancels each off-diagonal
column coefficient and gives E^2=I. Hence its inverse is integral.
Left and right unimodular multiplication preserves integer Smith factors
and determinant. Over any field the rank is preserved. For every modulus
q, E^(-t) gives the actual isomorphism between the old and new kernel
modules. No assertion of a residue-only mutation function is needed.
The Pfaffian transformation for even dimension also has the asserted sign.

## The triangle basis is integrally unimodular

For the displayed cycle a->b->c->a, the columns
e_a,e_b,z=e_a+e_b+e_c,e_r+e_a-e_b form a triangular change of basis with
diagonal1. Pairing e_a and e_b gives1. Their pairings with z vanish.
For f_r=e_r+e_a-e_b the two mixed pairings are

    <e_a,f_r>=1-1=0,     <e_b,f_r>=1-1=0.

The row sum of the triangle form is zero, so <z,f_r>=3. In the residual
pairing all correction terms cancel, leaving <f_r,f_s>=(B_R)_rs.
This proves the claimed J2 direct sum D3(R) over the integers.

Modulo3 the isolated z coordinate contributes exactly one additional
kernel dimension, while J2 is invertible. Scaling the source row and
column of B_(1->R) by3 gives D3(R), hence the factor9 in the determinant.
The odd-dimension zero-determinant boundary, including empty R, is
correctly stated. The independent checker additionally verifies the basis
on132 small skew residual matrices, including zeros and weights±2.

## The six-vertex hostile excludes every ordinal-triangle endpoint

Directly from the literal masks,

    U=(52,33,34,39,46,0)

has a global sink at5 and indegree-one vertices3 and4. The first five
vertices form the supplied strong tournament. V=(58,60,57,16,32,8)
is C3->C3 and has indegree-one vertices0,1,2. Independent literal
enumeration of all720 vertex orders gives H(U)=H(V)=9.
The determinant sums over all720 permutations give25 and81.

Every six-vertex tournament with an ordinal C3 factor can be written
L->C3->R with |L|+|R|=3. If both L and R are nonempty, each has at
most two vertices, hence is transitive. Mutating its successive global
sources moves the left block behind the graph, making the triangle
extreme without changing determinant. If the factor is already extreme,
the three-vertex residual is transitive or cyclic. Formula M2 then gives
determinant9 or81. Thus excluding the determinant25 source uses no
all-length H-invariance assumption.

The independent checker covers this classification with all40 canonical
block templates: all four values of |L|, both triangle orientations, and
every tournament on L and R. Relabelling preserves determinant, so these
templates cover every labelled endpoint in question. Their determinant
set is exactly{9,81}. This is a complete small interface check, not an
enumeration of ambient mutation classes.

## Prime-count addendum: an actual17 return with a lattice obstruction

The supplied H5 graph is G=(10,12,1,4), marked at1. Its literal rooted
counts are(H,q,s,F)=(5,3,2,6), so the root has defect F-H=1.
Adding a vertex4 beaten only by1 gives(10,28,1,4,13); adding vertex5
beaten only by4 gives exactly

    A=(10,28,1,4,45,15).

Thus A is an actual3H+2 child at prime count17. Its two indegree-one
roots are4 and5. The comparison chart is

    B=(58,60,1,52,36,4),

also at count17, with indegree-one roots0 and1. Independently enumerating
all paths and ordered two-bin path covers gives:

| Source/root | H | H(delete root) | Paths starting at root | E(root,2) | Defect |
|---|---:|---:|---:|---:|---:|
| A/4 | 17 | 5 | 5 | 17 | 0 |
| A/5 | 17 | 9 | 6 | 20 | 3 |
| B/0 | 17 | 1 | 1 | 17 | 0 |
| B/1 | 17 | 9 | 8 | 18 | 1 |

The defect is checked both as F-H and as q-s, rather than assuming their
equality in the counting code. Direct permutation determinant sums give
det B_A=49 and det B_B=1. Therefore no ordinary FZ mutation word of any
length carries A to B, regardless of intermediate weights or relabelling.
This obstructs this proposed same-order prime-count rearming step.
It does not exclude every other ready source at count17, a different
vertex order, or operations including birth/deletion.

## The residue hostile is valid

The two skew matrices with (b_01,b_02,b_12)=(1,1,1) and(-2,1,1)
agree modulo3. At pivot1, the first has a positive two-step contribution1,
while the second has mixed signs and contribution0. Their new b_02
entries are therefore2 and1, distinct modulo3. This refutes a global
residue-only update rule. It does not claim that these two matrices
are in one specified tournament's mutation class.

The warning is compatible with integral invertibility: each individual
step is a congruence once its integer sign-dependent E is fixed, while
residues alone do not select that E. The Collatz comparison keeps exactly
this information boundary and does not claim a dynamical conjugacy.

## Reproduction

Run `python3 -B audit_matrix.py` and `python3 -B -O audit_matrix.py`.
The controls are375 weighted three-vertex mutation checks,132 triangle
basis checks,40 ordinal-triangle templates, two literal H9 sources, and
the H5-to-H17 birth plus two fully counted H17 rooted sources.
Normal and optimized outputs agree. Universal conclusions rest on the
identities and the structural argument above.
