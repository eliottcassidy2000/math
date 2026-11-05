# Physical cyclic counts as positive Gram kernels and rooted-forest determinants

Session197. **PROVED all-order finite-kernel theorem.** Fix any nonempty tournament V. The counts H(C3[T,U,V]) form a positive-definite bilinear kernel on finitely many path-cover coordinates of T and U. Its coefficient matrix is a totally nonnegative Jacobi matrix after positive factorial scaling, and its determinant is an explicit weighted rooted-forest sum on a path. This is not a TN assertion about every single-module substitution operator.

The literal fragment principle is inherited from the pinned S195 path-cover algebra: https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/pathcover-algebra/PROOF.md . The tridiagonal kernel and forest interpretation are derived below without a literature-priority claim.

## 1. Exact physical tridiagonal kernel

Let m=|V|>=1 and L=m+1. Write c_r(A) for the number of unordered r-path covers of a tournament A, and put

    f_r=r! c_r(V) for1<=r<=m,
    f_0=f_(m+1)=f_(m+2)=0.

All f_r for1<=r<=m are positive integers. Indeed every partition into r nonempty vertex sets has at least one Hamilton path in each block. For any nonempty T and U use their ordered-cover vectors

    x_r=r! c_r(T), y_r=r! c_r(U), 1<=r<=L.

Then the physical cyclic quotient T->U->V->T satisfies

    H(C3[T,U,V])=x^T J_V y,                              (1)

where J_V is symmetric tridiagonal with

    (J_V)_(r,r)=f_(r-1)+3f_r+f_(r+1),
    (J_V)_(r,r+1)=(J_V)_(r+1,r)=f_r+f_(r+1).            (2)

Coordinates with more than L runs cannot occur.

### Proof and normalization

Cut a Hamilton path at all interblock edges. The maximal run colors follow T,U,V,T cyclically. Given their counts i,j,k, the number N_(i,j,k) of legal color words is the number of their three possible starting colors that give those counts. Assigning the distinct local cover fragments to their color occurrences gives i!j!k! possibilities.

A legal word has run counts differing by at most1. If all three counts are equal, N=3. Otherwise, when their maximum and minimum differ by1, exactly one starting color gives the specified counts. Hence for i=j=r the possible V-run counts are r-1,r,r+1 with weights1,3,1, omitting0 because V is nonempty. For j=i+1 the possible V-run counts are i,i+1, each with weight1. All other pairs i,j have coefficient0. Summing the actual V covers gives precisely (2). Cutting and concatenating fragments are inverse, and no adjacent same-color fragments or closed fragment cycles are permitted.

The matrix on unordered cover vectors is

    K_V=D J_V D, D=diag(1!,2!,...,L!).                    (3)

This distinction prevents an ordered/unordered factorial error.

## 2. A positive weighted-square decomposition

Set w_r=f_r+f_(r+1)>0 for1<=r<=m. Directly expanding gives, for every real vector x,

    x^T J_V x = 2f_1 x_1^2
                + sum_(r=2)^m f_r x_r^2
                + sum_(r=1)^m w_r (x_r+x_(r+1))^2.       (4)

The empty middle sum when m=1 is intentional. If (4) vanishes, the first term forces x_1=0; the edge terms successively force x_2,...,x_L=0. Thus J_V is positive definite. This proof uses actual positive cover counts of V and does not assume positivity survives an arbitrary algebraic change.

Since J_V is symmetric tridiagonal with positive off-diagonal entries and positive definite, its LDL^T factorization has positive diagonal pivots and a unit lower bidiagonal factor with positive subdiagonal entries. Each such factor is TN, so finite Cauchy-Binet proves J_V is totally nonnegative. Positive diagonal scaling in (3) preserves TN, giving the same conclusion for K_V. Strict positivity of every minor is not asserted; tridiagonal sparsity forces some zeros.

## 3. Exact weighted rooted-forest determinant

Conjugate J_V by diag(1,-1,1,-1,...). This does not change its determinant and makes the off-diagonals -w_r. The resulting matrix is the weighted path Laplacian plus diag(lambda), where the path has vertices1,...,L and

    lambda_1=2f_1,
    lambda_r=f_r for2<=r<=m,
    lambda_L=0.

Consequently

    det J_V = sum_(F subset of path edges)
                product_(e in F) w_e
                product_(components C of F) sum_(v in C) lambda_v. (5)

Equivalently, choose one marked root in every forest component and multiply its vertex weight with all selected edge weights.

For completeness, write the conjugated matrix as B diag(w) B^T+diag(lambda), with B an oriented incidence matrix of the path. Apply Cauchy-Binet to the combined edge-incidence and weighted coordinate columns. An L-column selection has a nonzero squared determinant exactly when its edge set is a forest with one selected coordinate/root in each component. The unsigned determinant is1: remove a leaf of a rooted tree component successively until only its selected coordinate remains. Other selections are singular by a rootless component or the column count. Summing the root choices gives (5). Zero vertex weights simply remove terms. The full path with a root at1 gives a positive term, consistent with positive definiteness.

Finally

    det K_V = (product_(r=1)^L (r!)^2) det J_V.            (6)

This is a finite exact determinant, not a limiting continued-fraction claim.

## 4. Physical Gram inequalities and exact rank

Reversal preserves every unordered path-cover count, so H(C3[T,T,V]) and H(C3[T,T^op,V]) are equal numerically. From (1) and positive definiteness,

    H(C3[T,U,V])^2
      <= H(C3[T,T^op,V]) H(C3[U,U^op,V]).                 (7)

Equality holds exactly when the first L ordered-cover vectors of T and U are proportional. It does not require or imply graph isomorphism. For a list of inputs, its pairwise cyclic-count matrix is positive semidefinite. It need not be totally nonnegative in an arbitrary ordering of those inputs.

The Gram rank is exactly L when inputs Tr1,...,Tr_L are allowed. Their ordered-cover matrix is triangular with entries r!S(n,r) and diagonal n!, so its determinant is nonzero. The corresponding physical Gram determinant is precisely det K_V, by (6). Longer lists have rank at most L for this fixed third block. This fixed-context finite rank is compatible with S195's unbounded-rank theorem when the third blocks range over unbounded sizes.

If V is self-converse with a specified reversing map, then C3[T,T^op,V] is self-converse by swapping the first two blocks and applying that map to V. In particular V=Tr_m has this property. Self-converseness is not asserted for an arbitrary V merely from the algebraic Gram identity. Every cyclic graph here is strong because all three blocks are nonempty.

## 5. The first two physical kernels and a genuine arithmetic distinction

For V=K1, f_1=1, and the unordered kernel is

    K=[3 2
       2 4].

If h=H(T) and u=c_2(T), its paired count is

    3h^2+4hu+4u^2=2h^2+(h+2u)^2.                        (8)

This family has no output7p with p prime. If7 divides (8), then both h and h+2u are divisible by7, because -2 is a quadratic nonresidue modulo7. Hence49 divides the output. A7p output would have to be49, but h>=7 would already make2h^2>=98. This argument treats the only possible exceptional prime p=7 explicitly.

For V=Tr2, f_1=1,f_2=2. Formula(2) gives

    J=[5 3 0      K=[5  6  0
       3 7 2         6 28 24
       0 2 2],       0 24 72].

The rooted-forest sum is det J=32, and det K=4608. This is exactly the separately frozen three-cover cyclic gate. Unlike (8), it has a physical7p witness: the input with cover vector(23,273,611) gives37050461=7*5292923. That witness and its complete restricted inverse are proved in `THREE_COVER_CYCLIC_GATE.md`; they are mentioned here as the demonstrated contrast, not re-certified by the present determinant argument.

The comparison says that adding one vertex to this particular exterior transitive block removes the universal norm obstruction. It does not imply that the larger kernel covers all prime or7p targets, nor that arbitrary quadratic-form vectors are realized by tournaments.

## 6. Reproduction and named dependencies

Run `python3 verify_cyclic_gram_forests.py` and its `-O` form. The checker enumerates legal run-color words, checks the Jacobi formula on actual input triples, verifies the weighted-square and rooted-forest determinant identities for arbitrary positive coefficient controls and actual third blocks, checks finite TN minors and the physical Cauchy-Schwarz inequality, and reproduces the first two kernels. All computations are exact.

The inherited path-cover gluing theorem is the pinned S195 source linked above. The computational helper is `verify.py`, SHA256 `10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb`. The arithmetic contrast cites frozen `THREE_COVER_CYCLIC_GATE.md`, SHA256 `0630c4ff7f9a52e69cbce4f18cdac8678609d4565a3e42f0afec8d588afe37cf`; the main kernel, positivity, forest and Gram proofs do not depend on that companion theorem.
