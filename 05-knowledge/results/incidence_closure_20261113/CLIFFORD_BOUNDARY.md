# Converse, a Clifford generator, and the self-converse dimension defect

**Status: PROVED IN NOTES, C1–C3, after independent review; see [AUDIT.md](AUDIT.md).**
This is a new exact representation map, with an explicit obstruction to a
stronger interpretation. It does not assert Bott periodicity for tournament
counts. The signed-incidence decomposition is inherited from the independently
audited candidate `pairings/SIGNED_LIFT.md`, S4.

## C1. A genuine Clifford representation, with a grading guard

Let X be the ordinary tournament classes at one order, R their converse
permutation, s the number of self-converse classes, and t the number of
two-element converse orbits. On V=R^X, R^2=I, so R defines an **ungraded**
representation of Cl_(1,0)=R[e]/(e^2-1). It also gives the decomposition

    V=V_+ direct-sum V_-,   dim V_+=s+t,   dim V_-=t.

Each fixed class contributes one positive vector; each pair contributes its
sum and difference. Thus trace(R)=dim V_+-dim V_-=s. This is an ordinary
graded dimension difference, without a claim that it is a topological index.
R itself is even relative to this decomposition: calling it an odd generator
of a graded Clifford module would be a type error.

An invertible linear J satisfying JR=-RJ exists iff s=0. Necessity follows
because J maps V_+ isomorphically to V_-. If the dimensions agree, choose
an isomorphism B:V_+ -> V_- and set

    J = [ 0  -B^(-1) ; B  0 ].

Then J^2=-I and JR=-RJ. Together R and J give an ungraded Cl_(1,1)
representation. This proves the iff and the negative-square extension at once.

Every tournament order has a transitive self-converse class, so the original
ordinary-class representation never meets s=0. If new positive and negative
dimensions are permitted, balancing requires k_- - k_+=s; the least total
addition is s, achieved by adding s purely negative dimensions. These are
formal vector dimensions, not new tournament classes. Tensoring with a
trivially graded R^16 multiplies the defect by 16; it does not balance it.

## C2. Preserving incidence adds a spectral obstruction

The symmetric ordinary dart matrix C commutes with R and splits into C_+
and D, where D is the signed block of `pairings/SIGNED_LIFT.md`. Even on a
balanced space, requiring [J,C]=0 imposes

    B C_+ = D B.

Consequently C_+ and D must be similar; conversely any such isomorphism B
gives the displayed commuting J. A balanced two-dimensional hostile has
R=diag(1,-1) and C=diag(1,2). Equal sector dimensions do not give an
incidence-preserving extension because the scalar operators 1 and 2 differ.
If exactly s new negative dimensions carry an operator E, the precise condition
becomes similarity of C_+ and D direct-sum E. An appropriate E is not supplied
by the dimension count alone. These are linear tests, separate from the
integrality and tournament-realization guards for actual darts.

## C3. Why this is not another class pairing

In the ordinary-class basis, R is a nonnegative permutation matrix. No nonzero
nonnegative matrix J can anticommute with it: JR and RJ are nonnegative, so
JR+RJ=0 forces JR=0, whence J=0 since R is invertible. In particular the
Clifford extension cannot be realized by another permutation of classes.
The involutions enumerated in the signed-lift note instead **commute** with R
and preserve C; these are different algebraic conditions.

At orders 5,6,7 the dimensions (s+t,t) are respectively (10,2), (34,22),
(272,184). The dimension defects are 8,12,88. The same s is the unavoidable
fixed-class exponent for the compatible involution polynomials when the
merged incidence graph is rigid. This is an exact common representation
source for those numbers, not a proof of a periodic pattern in them.

The map is ordinary classes -> their real permutation representation with R.
It preserves the converse character and, when C is included, the full signed
incidence block. It forgets the preferred class basis unless that basis is
retained. Cl_(1,1) extension requires balancing and a chosen intertwiner;
Bott/graded Morita periodicity needs additional Clifford generators and the
appropriate category. No eight-order shift or scalar tournament recurrence
follows. The flat-transport versus underlying-bundle counterexample in
`bridges/PROOF.md`, B5, is a separate warning about a different forgetful map.
