# Sharp stability on the H=1 locus

Status: PROVED, independently audited. All orders are covered.

## Inheritance and scope

A is a real weighted tournament: zero diagonal, A_ij+A_ji=1 and entries in
[0,1]. H is its weighted directed Hamilton-path sum, u=det(I+A),
v=(-1)^n det(A), V=sum_{i<j} A_ij A_ji, and Delta=H-u-|v|.
The Session137 global theorem gives Delta>=0, and its polynomial stability
coefficient beta_n is sharp at n=3,4 but not claimed optimal thereafter.
The extremal hostile to dimension-free stability is an ordinal sum of fair
2-vertex blocks, with at most one singleton.

Live concept board: binary rounding / H=1 support; pair-block product /
variance; orientation orthants / exact Bernstein coefficients; signed
root recurrence / triangle margin; dimension parity / perfect matchings.
The present note closes the first two concepts, not the full cube.

## The H=1 locus

**Theorem.** H(A)=1 iff, after relabeling, A is an ordinal sum of singleton
blocks and genuinely fractional two-vertex blocks.

Independently orient each edge i,j with probability A_ij. Multiaffinity of
the Hamilton-path sum gives H(A)=E H(T). Every binary tournament has at
least one Hamilton path and has exactly one iff it is transitive. Thus
H(A)=1 iff every completion with positive probability is transitive.

Two fractional edges sharing a vertex, say ij and jk, allow either
orientation of each. Whichever orientation of ik has positive probability,
one can choose ij,jk to complete a directed triangle. Therefore fractional
edges form a matching. If an outsider sees the two endpoints of a
fractional pair differently, choosing the appropriate orientation of that
pair again creates a supported triangle. Hence every outsider faces a pair
uniformly; in particular the deterministic edges between two fractional
pairs all face the same way. Contract the pairs. The quotient is a binary
tournament with no directed triangle, and hence is transitive. This gives
the claimed ordinal decomposition. Conversely each completion of such an
ordinal sum is transitive, so H(A)=1.

The binary unique-Hamilton-path characterization follows directly from the
inherited odd-cycle formula: H=1 plus nonnegative odd-cycle-family terms,
and every directed triangle contributes 2. A nontransitive binary
tournament contains a directed triangle, whereas a transitive tournament
has only its ordered Hamilton path.

## Exact scalar minimization

Write r for the number of fractional pairs, t for the number of singletons,
q_i=p_i(1-p_i) in (0,1/4], and S=sum_i q_i. Block triangularity gives

    H=1, u=prod_i(1-q_i), V=S,
    |v|=prod_i q_i if t=0, and |v|=0 if t>0.

Put epsilon=1_{t=0}. With all q except one fixed, write
C=prod_other(1-q_i), D=prod_other q_i, S0=sum_other q_i. Then

    Delta/V = [1-C+(C-epsilon D)q]/[S0+q],

whose derivative has the sign of

    C(1+S0)-1-epsilon D S0 <= 0.

Indeed prod_i(1-q_i)<=1/(1+sum_i q_i): multiply the elementary inequalities
1/(1-q_i)>=1+q_i and then expand the product on the right. For at least
one other positive q, this inequality is strict. Thus for r>=2 the ratio
strictly decreases in every q_i, and is minimized iff all pairs are fair.
For r=1,t>0 it is identically 1; for r=1,t=0 (n=2) it is identically 0.

For t>0 the minimum at fixed r is

    a_r = (4/r)(1-(3/4)^r).

The sequence a_r decreases strictly because it is the arithmetic mean of
1,(3/4),...,(3/4)^(r-1). When t=0 the minimum is

    b_r = (4/r)(1-(3/4)^r-(1/4)^r) < a_r.

For n>=3, allowing all possible pair/singleton decompositions, the sharp
H=1 coefficient is therefore

    k_(2r)   = (4/r)(1-(3/4)^r-(1/4)^r), r>=2;
    k_(2r+1) = (4/r)(1-(3/4)^r),           r>=1.

For even n, any decomposition with singletons has fewer than r pairs, so
its ratio is at least a_(r-1)>a_r>b_r. For odd n the maximal r is optimal.
Equality with V>0 at n>=4 is exactly a maximal number floor(n/2) of fair
pairs, with no singleton when n is even and one singleton when n is odd.
At n=3 every H=1 matrix with V>0 has equality, in agreement with the global
identity Delta=V.

The constants at n=3..10 are 1,3/4,7/8,3/4,37/48,87/128,175/256,39/64.
They are asymptotic to 8/n. They are global UPPER bounds on the optimal
stability constant, not yet all-order global LOWER bounds.

## Separate-edge interpolation identity

For any edge e={i,j} of variance q_e, let A^0,A^1 be its two endpoint
orientations and R=A-{i,j}. Direct determinant expansion gives

    u(A)=(1-p)u(A^0)+p u(A^1)-q_e u(R),
    v(A)=(1-p)v(A^0)+p v(A^1)-q_e v(R),

with the endpoint convention adjusted to the definition of p. H is affine
in that edge. Hence each signed margin G_s=H-u+s v obeys

    G_s(A)=(1-p)G_s(A^0)+pG_s(A^1)+q_e(u(R)-s v(R)).

This is a usable mechanism; it does not alone prove the sharp coefficient,
because u(R)-s v(R) can lie below that coefficient.
