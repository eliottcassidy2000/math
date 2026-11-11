# High degree is not a new topological type

**PROVED elementary normalization statements; inherited Keller monoid facts;
OPEN all-degree Jacobian conclusion.** Work over C throughout.

Write `pdeg(F)=max(deg P,deg Q)` for the total polynomial degree of
`F=(P,Q)`, and `d(F)=[C(x,y):C(P,Q)]` for its generic fibre degree. These
are different quantities. A Keller map has constant nonzero Jacobian.

## 1. An exact degree-inflation theorem

For any Keller map F, both P and Q are nonconstant. For an integer N with
`N deg Q > deg P`, put

`A_N(u,v)=(u+v^N,v)` and `F_N=A_N o F=(P+Q^N,Q)`.

Then:

1. `det JF_N=det JF` and `pdeg(F_N)=N deg Q`.
2. F_N is a polynomial automorphism if and only if F is.
3. `d(F_N)=d(F)`.
4. Their nonproperness sets satisfy `S_(F_N)=A_N(S_F)`; the covers over
   their complements are isomorphic after this target coordinate change.

Proof. The inverse of A_N is `(u-v^N,v)` and its Jacobian is one. The
strict degree inequality prevents cancellation of the leading term of
Q^N, proving (1). Composition with an automorphism proves (2), and
`C(P+Q^N,Q)=C(P,Q)` proves (3). If x_j escapes every compact set and
F(x_j) tends to y, then F_N(x_j) tends to A_N(y). Conversely, convergence
of F_N(x_j) implies convergence of F(x_j) by the continuous inverse of
A_N. This proves (4), including the converse inclusion. QED.

More generally, source and target polynomial automorphisms A,B preserve
invertibility, generic degree, and properness, with

`S_(A o F o B)=A(S_F)`.

A homeomorphism and its inverse take compact sets to compact sets, so B
preserves escaping sequences. This also gives a direct proof of the last
formula, without assuming a general nonproperness composition theorem.

**Consequence.** If any planar Keller counterexample exists, its
two-sided automorphism orbit contains counterexamples of arbitrarily
large total degree, with the same covering and boundary topology.
Therefore a high-degree representative does not require a new topological
anomaly. The real issue is whether a noninvertible orbit exists at all,
and whether its *least-complexity representative* can be controlled.

The positive control `(x+y^N,y)` is an automorphism of degree N with
inverse `(x-y^N,y)`. Neither high degree nor a long recursive expression
by itself measures a covering defect. This is analogous to adding
combinatorially harmless padding: the analogy concerns normalization,
not equality of tournament and Keller invariants.

## 2. What a degree bound would have to mean

Define

`mu(F)=min_{A,B in Aut(C²)} pdeg(A o F o B)`.

The minimum exists because it is the minimum of a nonempty set of positive
integers. This definition gives no algorithm for finding it. Invertible F
has mu(F)=1; conversely mu(F)=1 implies invertibility. Thus merely saying
that every mu is one is a reformulation of JC(2), not an intermediate
theorem.

Useful non-tautological targets must constrain a *specified computation*,
for example a degree-reducing source/target move with a verified stopping
rank, or a faithful bounded boundary grammar whose repeated states admit
an actual polynomial contraction. A uniform bound mu(F)≤D for nonunit
composition atoms would reduce the task to a bounded-degree family of
polynomial systems, but would still require proving that family has no
solutions. A bounded degree is not a finite list of complex coefficients.

## 3. Composition termination is inherited, not enough

THM-1330, `keller-monoid-exact-picture-inverse-jelonek-cusp-rule.md`, proves
the monoid facts used here: generic degree is multiplicative, and generic
degree one Keller maps are automorphisms. Its finite factorization into
nonunit composition atoms is inherited, not newly claimed here.

An additional elementary observation clarifies its domain. If a dominant
polynomial map F factors as G o H with G,H maps C²→C² and det JF a nonzero
constant, the chain rule gives `(det JG)(H) det JH=c`. The units of
C[x,y] are constants, so det JH is constant nonzero. H is therefore
dominant, and injectivity of its coordinate-ring pullback implies det JG
is also constant nonzero. Thus an actual polynomial factorization
automatically stays inside the Keller monoid.

This does not apply to an abstract block system of monodromy: its
intermediate normalization need not be C². Nor does every recursively
drawn boundary tree yield a polynomial factorization. The affine-source
recognition condition is the important additional structure in
[the boundary recursion bridge](../topology/boundary_recursion_bridge.md).

## 4. Scope of the geometric question

JC(2) here concerns polynomial maps of complex dimension two, hence real
dimension four. Their infinity divisors are complex curves and their
complements carry nontrivial monodromy. Classification of real flat
surfaces does not classify these polynomial maps or their boundaries.
The promising replacement for an unspecified exotic anomaly is a precise
failure of chart entry, realizability, or termination. The other notes in
this packet expose those failures as separate mathematical questions.

## Evidence

The all-N proof is above. `python3 degree/check.py` checks explicit
compositions and inverse formulas at finite parameters; it is a control,
not an all-degree census. THM-1330 was read at starting main
`1432c0ac0d987f1f7d8219a04c5326b2826e552d`. Its unrelated higher-dimensional
claims are not used in this note. No global classification is asserted.
