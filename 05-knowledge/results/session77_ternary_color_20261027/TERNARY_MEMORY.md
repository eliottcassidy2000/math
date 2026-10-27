# Ternary composition remembers more than fifteen Hamiltonian paths

**Status: elementary proof and exact finite checks.** The substitution theorem
below is inherited from the repository's audited THM-3121, not claimed as a
new discovery. The explicit small witness and its independent checks make the
loss of a scalar summary concrete. No assertion about physical QCD follows.

## 1. What a ternary node must retain

For a tournament T, let pc_T(r) be the number of spanning covers by r unordered,
vertex-disjoint directed paths. Singletons are allowed. In particular,
pc_T(1)=H(T), the number of Hamiltonian paths.

For a quotient tournament Q with adjacency matrix A_Q, put X=diag(x_i). Its
nonempty directed-walk content series is

\[
 W_Q=\mathbf1^T X(I-A_QX)^{-1}\mathbf1.
\]

For nonempty factors S_i, the exact substitution law is

\[
 \operatorname{pc}_{Q[S_1,\ldots,S_q]}(d)
 =\sum_{c_i\ge1}
 \frac{\prod_i c_i!}{d!}[x^c]W_Q^d
 \prod_i\operatorname{pc}_{S_i}(c_i).
\]

The sums terminate at c_i=|S_i|. This says the complete path-cover profile is
a closed counting carrier for tournament substitution; it does not say that
computing such a profile is easy or that this is a globally minimal encoding.

### Why the factorials are necessary

Split every global path at changes of factor block. Within block i, the runs
form a spanning c_i-path cover. The block names along each global path form a
quotient walk. Conversely, choose each factor cover and an ordered d-tuple of
quotient walks. Its c_i occurrences of block i can be assigned the distinct
actual factor paths in c_i! ways. Now divide by d!: the assembled global
paths use disjoint nonempty vertex sets, so permuting them acts freely.

Bare quotient words can repeat. Division by d! **before** actual-path
decoration is not justified. For example, the ordered pair of singleton
quotient words (1),(1) has coefficient one in W_Q²; dividing this bare count
by 2! gives 1/2. Assigning two distinct actual singleton components first
contributes 2!, and the correct unordered cover count is one. The same issue
occurs with every factor present: [x²yz]W_C3⁴=12, so bare division by 4!
again gives 1/2, while decorating the two occurrences of the first block
contributes 2! and restores the integer one.

This ownership distinction is the precise memory issue. It is not enough to
remember an abstract word or the scalar number of full paths.

## 2. The directed-triangle law

For C3 oriented 1→2→3→1, a walk is determined by its starting port and length.
Therefore

\[
 W_{C3}=\frac{x+y+z+xy+yz+zx+3xyz}{1-xyz}.
\]

For positive a,b,c its coefficient is 3 when a=b=c, 1 when
max(a,b,c)−min(a,b,c)=1, and zero otherwise. Define
F_i(r)=r! pc_{S_i}(r), zero outside the factor's range. Then

\[
 H(C3[S_1,S_2,S_3])
 =3\sum_r F_1(r)F_2(r)F_3(r)
 +\sum_r\sum_i F_i(r+1)\prod_{j\ne i}F_j(r)
 +\sum_r\sum_i F_i(r)\prod_{j\ne i}F_j(r+1).
\]

The three-port transfer sees the diagonal and its two nearest balanced shells.
It is not the product of the three scalar H values.

The triangle itself has pc=(3,3,1), hence F=(3,6,6). The law gives

\[
 H(C3[C3,C3,C3])=1377+810+972=3159,
\]

independently checked by enumerating all vertex permutations. The complete
nine-vertex profile is

    (3159, 21303, 42201, 36450, 16119, 3906, 522, 36, 1).

## 3. Two strong five-vertex factors with the same H=15

Let A have the following ten arcs on vertices 0,1,2,3,4:

    0→1, 0→3, 1→4, 2→0, 2→1,
    3→1, 3→2, 4→0, 4→2, 4→3.

The set {0,2,3} is a cyclic module, with

    {0,2,3} → 1 → 4 → {0,2,3}.

Thus A is C3[C3,K1,K1], up to the displayed names. Let B have arcs

    0→3, 0→4, 1→0, 1→4, 2→0,
    2→1, 3→1, 3→2, 4→2, 4→3.

B is the regular cyclic tournament R5. Both are strong, but their outdegree
lists are (2,1,2,2,3) and (2,2,2,2,2), so they are not isomorphic. Exact path
cover counts are

\[
 \operatorname{pc}_A=(15,39,33,10,1),\qquad
 \operatorname{pc}_B=(15,45,35,10,1).
\]

Both have fifteen Hamiltonian paths, but they differ on two-path covers.
Their automorphism orders are 3 and 5, respectively; thus H/|Aut| is 5 versus
3. This is a collision of H, **not** a collision of the normalized count.

The checker independently enumerates every labelled tournament through order
four (75 inputs in total), using two different counting procedures. None has
the same H with different complete profiles. Orders one and two are included.
Consequently order five is the smallest order for this particular kind of
scalar/profile collision.

## 4. One fixed ternary attachment observes the missing memory

Attach the same two named singleton factors to either T=A or T=B, in the
same cyclic quotient orientation:

    T → u → v → T.

The factors u and v each have profile (1). In the triangle coefficient law,
the only allowed contents are (1,1,1) and (2,1,1). Therefore, for every
nonempty tournament T,

\[
 \boxed{H(C3[T,K1,K1])=3H(T)+2\operatorname{pc}_T(2).}
\]

It follows that

\[
 H(C3[A,K1,K1])=123,\qquad H(C3[B,K1,K1])=135.
\]

Both seven-vertex outputs are strong. Their full profiles are

    A-context: (123, 555, 825, 530, 158, 21, 1),
    B-context: (135, 615, 885, 550, 160, 21, 1).

Thus an equal scalar H at a factor is not preserved as equality after the
same lawful ternary attachment. The missing two-run count is observable.
This is exactly a failure of compositional sufficiency, rather than a
numerical resemblance between fifteen and another occurrence of fifteen.

## 5. Relation to typed color contractions, and limits

The analogy with typed tensor composition is precise at the level of what
must be preserved: graph substitution needs factor-path runs and their
assignment to quotient-word occurrences, while color tensor contractions
need boundary representation types, port order and intertwiner channels.
The operations themselves are not identified. Path covers are nonnegative
integer counts of disjoint combinatorial objects; color tensors have signed
or complex coefficients and linear relations. Factorials in this theorem
must not be substituted for SU(3) recoupling coefficients.

Likewise a quotient C3 has three named factor ports. It is not automatically
a full rooted ternary tree, a trivalent color-web vertex or a residue-refinement
tree. Repeated composition is legitimate only with the appropriate carrier
and its actual operation retained.

## 6. Reproduction and exactly what was checked

The companion `check_ternary_memory.py` is a standalone standard-library
script. Run it with either ordinary Python or optimized Python. It uses
explicit exception checks, not assertions that disappear in optimized mode.

Two independent methods compute profiles:

1. Held–Karp subset Hamiltonian counts, then a set-partition recurrence that
   assigns the smallest remaining vertex to a unique path block.
2. Vertex permutations with all legal cut sets, followed by division by d!
   for the order of the actual, disjoint global paths.

The script verifies both complete factor and seven-vertex context profiles,
75 small tournaments, 216 positive C3 kernel cells, 64 complete C3-substitution
profiles, 64 complete general-quotient profiles, 4,350 integral kernel cells,
and H=3159 by independent permutations. It checks strongness, automorphism
orders and the concrete C3/R5 identifications. Three deliberate false claims
are rejected: scalar sufficiency, wrong diagonal coefficient and division
before decoration. Ordinary and optimized outputs agree byte-for-byte.

Only this explanatory note and the companion source are intended as the lean
public research package; audit outputs and run bookkeeping are not needed to
read the mathematics or run the checker.

## Source and credit

The exact composition theorem, rational walk kernel and ownership repair are
inherited from
[THM-3121: the path-cover walk-content substitution kernel](https://github.com/eliottcassidy2000/math/blob/02d1218a9a48986564693cae1ae65a050fb9f0e6/01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md).
That source credits THM-1975 for earlier scalar-H hostile examples. This note
makes no claim to discovery priority for the obstruction or its formula.
