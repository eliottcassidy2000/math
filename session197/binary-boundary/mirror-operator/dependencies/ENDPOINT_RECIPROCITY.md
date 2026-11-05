# Endpoint reciprocity from a squarefree path-kernel identity

Session 197. **PROVED exact signed-to-physical endpoint identity.** Let S be a tournament and a one of its vertices. Then

    sum_(directed paths P ending at a) (-1)^|P| H(S-P)
       = -H_start_at_a(S).                               (1)

All paths in the sum are nonempty ordered vertex words. The complement is induced, H(empty)=1, and the right side counts Hamilton paths whose first vertex is a. In particular the alternating sum is always nonpositive. This is an exact combinatorial count, not positivity inferred from arbitrary algebraic manipulation.

## 1. A finite commutative ring of labelled subsets

Work over the squarefree ring

    R=Q[x_v : v in V(S)] / (x_v^2 : v in V(S)).

For a subset U write x^U=product_(v in U) x_v. The variables commute. A product with a repeated vertex vanishes, so disjoint-subset convolution is ordinary multiplication in R. Every element with constant term 1 is invertible by a finite geometric series.

Let A be the zero-diagonal adjacency matrix of S and X=diag(x_v). Define

    K(x)=X(I-A X)^(-1),
    F(x)=1+1^T K(x) 1.

The inverse is a finite series: every sufficiently long walk repeats a labelled variable. The entry K_ij is exactly the sum of x^P over directed vertex-simple paths starting at i and ending at j. Repeated-vertex walks vanish. Consequently

    F(x)=sum_(U subset V(S)) H(S[U]) x^U.                 (2)

In particular K_ii=x_i, and the empty subset contributes 1 only through the initial constant in F.

## 2. The exact rank-one resolvent relation

The tournament relation is

    A^T+I=J-A, J=1*1^T.

Adding diagonal entries to an adjacency matrix does not change the kernel X(I-BX)^(-1) in this ring: every walk using a self-loop contains x_v twice. Therefore

    -K(-x)^T
       = X(I+A^T X)^(-1)
       = X(I+(A^T+I)X)^(-1)
       = X(I-A X+J X)^(-1).

Put M=I-A X and x^T=1^T X. The last inverse is (M+1*x^T)^(-1). The Sherman-Morrison identity holds over this commutative ring because M and the denominator 1+x^T M^(-1)1=F are invertible. It gives

    -K(-x)^T = K(x) - (K(x)1)(1^T K(x))/F(x).             (3)

One can also verify this identity by multiplying both sides by M+1*x^T, so no analytic convergence or division by a zero polynomial is involved. A denominator-free form convenient for exact checking is

    F(x) [K(x)+K(-x)^T] = (K(x)1)(1^T K(x)).             (4)

The diagonal-removal step above applies to the path kernel with its initial X. It is not a claim that diagonal entries may be removed from a determinant or arbitrary matrix inverse.

## 3. Endpoint extraction and the sign

Multiply (3) on the right by 1 and use 1^T K1=F-1. This gives

    F(x) K(-x)^T 1 = -K(x)1.                            (5)

The a-th coordinate of K(-x)^T1 enumerates paths ending at a, with sign (-1)^|P|. In the product with F, repeated vertices vanish, and extracting the coefficient of x^(V(S)) forces the other factor to be the full complement of P. The left coefficient is therefore exactly the sum in (1). The a-th coordinate of K1 enumerates paths starting at a; its full-vertex coefficient is H_start_at_a(S). This proves (1), including its minus sign and the singleton boundary.

Transposing the argument also gives the symmetric identity with paths starting at a on the left and Hamilton paths ending at a on the right. The nonpositivity conclusion applies to these specific alternating endpoint sums. It does not make all coefficients of the rational-looking expressions in (3) nonnegative.

Summing (5) over a yields F(x)F(-x)=1. This familiar unrooted tournament path reciprocity is recovered as a consequence; the proof here retains the endpoint information lost by that scalar identity.

## 4. Boundary interpretation and examples

For a one-vertex S, the left side is -1, agreeing with its unique Hamilton path starting at a. In a transitive tournament, the left side is -1 when a is the unique source and 0 at every other vertex.

More generally, if a separate construction has a correction of the form -d^2 with

    d=sum_(P ending at a) (-1)^|P| H(S-P),

then (1) identifies that correction as -H_start_at_a(S)^2. This is a conditional substitution into the stated correction formula; this packet does not independently assert or prove a particular graph construction's correction formula.

The actual values need not be only 0 or -1. The checker records a five-vertex rooted tournament for which H_start_at_a=5 and hence the alternating response is -5.

## 5. Reproduction and scope

Run `python3 verify_endpoint_reciprocity.py` and `python3 -O verify_endpoint_reciprocity.py`. The checker enumerates literal directed path words and independently computes induced-subset Hamilton-path totals. It verifies (1) for all 5,405 rooted labelled tournaments through order five. It separately checks the entire squarefree matrix identity (4), the vector identity (5), and F(x)F(-x)=1 on all labelled tournaments through order four and additional deterministic order-five and order-six controls. These exact finite controls supplement the all-order matrix proof.

The verifier imports only the basic physical graph utilities and subset Hamilton-path counter in `verify.py`, SHA256 `10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb`. It implements squarefree multiplication and the kernel checks directly. No unverified positivity, real-rootedness, spectral approximation or probabilistic identity test is used.
