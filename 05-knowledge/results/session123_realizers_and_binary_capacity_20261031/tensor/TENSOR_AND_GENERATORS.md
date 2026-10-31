# Binary tensor, unimodular capacity Gram matrix, and minimal generators

**PROVED deductions from the all-budget C5 capacity theorem.** The kernel and
Boolean incidence identities below are elementary/classical; no priority
claim is made. These are statements about the specified composition class's
capacity, not unrestricted Ramsey numbers or graph adjacency spectra.

Write

    F(A,B)=sum_(0<=i<A,0<=j<B) 2^popcount(i & j),
    N_s(n)=#{0<=i<n : every one-bit of s is a one-bit of i}.

A,B are positive integers, extended to zero by F(0,B)=F(A,0)=0. Let
L=2^ceil(log2 A), with L=1 when A=1.

## 1. Exact Boolean Gram and tensor representations

Every subset of the common one-bits contributes once, so

    2^popcount(i & j)=sum_s [s subset_bits i][s subset_bits j],
    F(A,B)=sum_s N_s(A)N_s(B).                             (1)

Only s<min(A,B) contributes. On indices0,...,2^d-1 the density kernel is

    K_d = [[1,1],[1,2]] tensor_power d.                    (2)

This is a literal Kronecker tensor identity, not merely shared terminology.
The tensor separates independent binary coordinates; taking rectangle sums
then supplies actual graph capacity by the separately proved C5 theorem.
Neither a rectangle indicator nor a general composition tree is claimed to
factor as a pure tensor.

For any N>=1 let Z_(i,s)=[s subset_bits i], 0<=i,s<N, and T be the lower
triangular matrix of ones. Z is unit lower triangular and

    (F(a,b))_(1<=a,b<=N) = (TZ)(TZ)^transpose.              (3)

Thus every such capacity matrix is positive definite, has determinant1,
and has an integer inverse. Its inverse is

    (Z^(-1)T^(-1))^transpose (Z^(-1)T^(-1)),

where T^(-1) is first difference and
Z^(-1)_(i,s)=(-1)^(popcount(i)-popcount(s))[s subset_bits i].
These identities follow from Boolean inclusion-exclusion, including at
non-power-of-two truncations because every subset of i is <=i<N.
Consequently F(A,B)^2<F(A,A)F(B,B) when A!=B; equality holds when A=B.
For A<B, the feature N_(B-1) is zero at A and one at B, while N_0 is
positive at both, so the feature vectors are not proportional.

The incidence Z is also the parity of the Pascal matrix. Over F2,
(1+x)^i=product_(j:bit_j(i)=1)(1+x^(2^j)); the coefficient of x^s is one
exactly when s is a binary subset of i. This proves the required binary
binomial identity directly, without importing an unproved number-theory bridge.

## 2. Every binary-subset counting function is superadditive

For all s,x,y>=0,

    N_s(x+y)>=N_s(x)+N_s(y).                               (4)

Proof by induction on s. At s=0 this is equality. Write s=d+t, where d is
its highest power of two and0<=t<d. Put c=N_t(d). The indicator defining
N_s has period2d; one period has c hits. Therefore it suffices to prove(4)
for residues0<=x,y<2d. In this range

    N_s(x)=0 for x<=d, and N_s(d+u)=N_t(u) for0<=u<d.

If both residues are <=d, the right side is zero. If precisely one exceeds
d, write it d+u. Without wrap, monotonicity of N_t proves the assertion.
With wrap, the new residue is <=d and the whole period contributes c,
which is at least N_t(u). If both exceed d, write d+u,d+v. The left side
is c+N_s(u+v). If u+v<=d, the induction hypothesis bounds the right side
by N_t(u+v)<=c. If u+v>d, it bounds the right side by
N_t(u+v)=c+N_t(u+v-d), using the period d of N_t. This is precisely the
left side. All residues and wrap boundaries are covered.

## 3. Complete additive-generator classification

For fixed A, call a positive budget b irreducible when

    F(A,b)>F(A,x)+F(A,b-x) for every0<x<b.

Budget1 is irreducible vacuously. The irreducible budgets are exactly

    {1,2,...,L}.                                          (5)

In particular the unique necessary irreducible item set is this complete
initial interval, and the least possible prefix bound generating the whole
capacity row by disjoint-union knapsack is M_A=L.

For0<x,y and x+y<=L, assume y>=x and let d be the highest power of two
not exceeding y. If x+y<=2d, then x<=d, and

    N_d(x+y)-N_d(x)-N_d(y)=x>0.

If x+y>2d, both x and y are less than2d and

    N_(2d)(x+y)-N_(2d)(x)-N_(2d)(y)=x+y-2d>0.

The selected power of two is always at most L/2<A (the case A=1 has no
such x,y). Its coefficient N_s(A) in(1) is therefore positive. Every other
term is nonnegative by(4). This proves strict superadditivity for every
split with total at most L.

For b>L, the exact period-from-origin law gives

    F(A,b)=F(A,L)+F(A,b-L),

so b is reducible. Any exact finite additive generating set must include each
irreducible budget: no sum of smaller budgets can attain it, by repeated
strictness. Conversely the first L values generate every b by the displayed
translation. This classifies capacity items, not every extremal graph or
its isomorphism/decomposition type.

## 4. Categorical boundary

Equation(2) is a genuine symmetric tensor computation in integer matrices.
Equation(3) retains every prefix feature instead of replacing it by a slope.
The exact capacity theorem supplies the separate realization map to graphs.
These numerical identities do not by themselves provide a realizability
PCA, an effective topos, a strong monoidal equivalence of graph categories,
or an unrestricted circuit lower bound.
