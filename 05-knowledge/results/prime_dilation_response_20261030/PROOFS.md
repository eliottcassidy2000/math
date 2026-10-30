# Joint prime-dilation response algebra

2026-10-30. **PROVED scoped**, same-agent proof review; independently computed
finite controls are listed in AUDIT.md. Local labels D1--D8 belong to this file.
No claim of priority for general Hankel, Smith, or observability methods.

Throughout, s(n)=sigma(n)-n for positive integers n. Fix a nonempty ordered
list P=(p_1,...,p_r) of distinct primes. Write

    n = u product_i p_i^a_i,   gcd(u, product P)=1,
    C = sigma(u),             R_p(a)=(p^(a+1)-1)/(p-1),
    d = 2^r,                 f(t)=s(n product_i p_i^t_i).

The legal letters are multiplication by these named primes, with
t in N^r. An aliquot step is NOT one of these letters. The observer is f,
over Z or a specified quotient ring. All source states here are actual
positive integers; formal linear combinations enter only in D5.

## D1. Factorized response and an integral closed cube

For every n,P,t,

    f(t)=C product_i R_(p_i)(a_i+t_i)-n product_i p_i^t_i.       (1)

The commuting shifts T_i satisfy

    (T_i-1)(T_i-p_i) f=0,
    f(t+2e_i)=(p_i+1)f(t+e_i)-p_i f(t).                        (2)

Thus the d entries z_e=f(e), e in {0,1}^r, are an integer carrier closed
under every legal word, also after reduction modulo ANY M>=1. On each
axis the update matrix is

    M_p = [[0,1],[-p,p+1]].                                   (3)

No division by p or p-1 is used in this update. Multiplicativity of sigma
proves (1); the two local sequences have roots 1,p, proving (2). The fully
mixed carry is particularly simple:

    product_i (T_i-p_i) f(t)=C,
    C=sum_(S subset P) [product_(p outside S)(-p)] s(n product S). (4)

Indeed R_p(a+1)-p R_p(a)=1, and any factor annihilates the second product
in (1). This extracts the divisor sum of the P-free core, not the core u.

For r>=2 the response tensor on the binary cube has tensor rank exactly
two over Q. It is a difference of two pure product tensors by (1).
The two local vectors (R_p(a),R_p(a+1)) and (1,p) have determinant -1,
so every nontrivial bipartition flattening has rank two. This low tensor
rank concerns the initial response tensor; it is not the dimension of
the linear space of all its continuations.

## D2. A nonlinear evaluator from inherited exact reconstruction

The source inverse is INHERITED from Session78, ALIQUOT.md, "any known
prime dilation recovers the source", and Session81's prime-power extension.
Given promised valid exact readings x=s(n), y=s(pn) for a known prime p,

    a=v_p(y-x),  c=y-px,  n=R_p(a)c-x,  sigma(n)=n+x.           (5)

For ANY supplied factored multiplier m=product_q q^b_q, let e_q=v_q(n).
Then

    sigma(mn)=sigma(n) product_(q|m)
                       (q^(e_q+b_q+1)-1)/(q^(e_q+1)-1),
    s(mn)=sigma(mn)-mn.                                      (6)

This composes the inherited inverse with the usual prime-exponent update.
For each distinct q the division replaces its sigma factor and is exact;
the unknown cofactor of n need not be factored. Local valuations suffice.
The pair must be valid, or Session78's divisor-sum realization guard must
be checked. Factoring an unprovided multiplier, obtaining the initial
observations, or testing arbitrary candidate observations is not free.
This is an efficient sufficient evaluator, not a minimal number-of-integers
claim: with unrestricted divisor computation, the integer n alone suffices.

## D3. Full integral Hankel determinant

For e,h in {0,1}^r define the actual d by d response matrix

    H[e,h]=f(e+h).

Put X=product_i p_i^(a_i+1), A=product_i(p_i-1), Gamma=C X-n A.
Since C>=u and product_i p_i/(p_i-1)>1, Gamma>0. Then

    det H=(-1)^(r 2^(r-1)) C^(d-1) Gamma
          product_i p_i^[(a_i+1)(2^(r-1)-1)].                (7)

In particular every actual source has rational Hankel rank exactly d.

Here is an integral proof that also retains the lattice lost by a rational
eigenbasis. For one axis put R=R_p(a), b=p^(a+1), R'=pR+1, and

    K=[[R,R'],[R',pR'+1]],
    U=[[-p,1],[R',-R]],       V=[[1,-1],[0,1]],
    v=(1,p)^T.

Direct multiplication gives det U=-1, det V=1 and

    U K V=diag(1,b),   Uv=(0,1)^T,   v^T V=(1,p-1).

Consequently, in subset order with the full subset last,

    (tensor U_i) H (tensor V_i)
       =diag(C product_(i in S) p_i^(a_i+1)) - n e_all w^T,
    w_S=product_(i in S)(p_i-1).                             (8)

This has only one off-diagonal row. Its last diagonal is Gamma.
Multiplying the diagonal entries and the signs of tensor U_i proves (7).
Both changes of basis are unimodular over Z, including at ramified or
coincident-root primes. Formula (8) is a full integer normal presentation,
not just a determinant computation.

Any rational linear weighted model with fixed commuting prime-update
operators, one source vector and a linear observation factors H through
its state space. Its dimension is therefore at least d. D1 attains d.
This proves linear minimality even though D1 has a two-product formula
and D2 has a much smaller nonlinear exact description.

## D4. Complete Smith form for every finite prime set

Write the positive Smith factors of H as h_1|...|h_d. For r=1,

    (h_1,h_2)=(gcd(n,C), C Gamma/gcd(n,C)).                   (9)

For r>=2, the following gives every primary layer. At an arbitrary prime
ell write c=v_ell(C), e=v_ell(n), g=v_ell(Gamma).

* If ell is NOT in P, the valuations in increasing order are

      min(c,e), c [d-2 copies], c+g-min(c,e).                 (10)

* If ell IS in P, sort the following multiset:

      c [d/2 copies], e, (c+e+1) [d/2-2 copies], c+1+g.       (11)

Multiply the resulting prime powers positionwise. Only primes dividing
C Gamma product P need occur. These formulas include common content,
singular reductions, and every size of ramification; over Z, H is nonsingular.

Proof. Work locally at ell using (8). The first d-1 diagonal entries
have valuations alpha_1<=...<=alpha_(d-1). Every nonzero k-minor is either
a product of k such diagonals, Gamma times k-1 of them, or one last-row
off-diagonal entry times k-1 diagonals excluding that entry's column.
These descriptions exhaust row and column choices, because every other
row has only its diagonal entry. The valuation of the gcd of all minors
is therefore the minimum of these possibilities; no cancellation is hidden.

If ell is outside P, all alpha_j=c and some off-diagonal entry has
valuation e (the empty-subset column). Also g>=min(c,e). For k<d the
k-th determinantal divisor has valuation (k-1)c+min(c,e); at k=d it is
(d-1)c+g. Consecutive differences prove (10).

If ell=p_i is in P, then e=a_i. There are d/2 low diagonals of valuation
c and d/2-1 high diagonals of valuation c+e+1. All off-diagonal valuations
are at least e. The proper column S={i} is high and attains e because
p_i-1 is a unit. For k<d one may exclude this high column and still
choose the k-1 cheapest diagonals. Since g>=e, putting
S_j=alpha_1+...+alpha_j gives the exact determinantal valuation

    v_ell(Delta_k)=S_(k-1)+min(alpha_k,e)  (k<d),
    v_ell(Delta_d)=S_(d-1)+g.

Taking differences gives (11); if e<c the first difference is e,
otherwise e occurs after the d/2 low differences. For r=1 the minimum
entry valuation is min(c,e), and the determinant is C Gamma, proving (9).

Example: n=1, P=(2,3), in bit order with 2 first,

    H = [[0,1,1,6], [1,3,6,16], [1,6,4,21], [6,16,21,55]],
    Smith(H)=(1,1,1,24).

## D5. Smallest closed linear carrier at every modulus

The universal integral algebra of these shifts is

    A_P=Z[X_1,...,X_r]/((X_i-1)(X_i-p_i) : 1<=i<=r).        (12a)

It is free of rank d on the squarefree monomials. The functional
lambda_n(X^e)=f(e) has pairing lambda_n(ab) with Gram matrix H.
For any commutative ring R and finite free commutative R-algebra A with
linear functional lambda, the future-invisible set

    I_lambda={a : lambda(ab)=0 for every b in A}

is an ideal: lambda((ca)b)=lambda(a(cb))=0. Consequently A/I_lambda is
the smallest quotient algebra preserving every weighted continuation.
This elementary statement is the general mechanism; D3--D4 explicitly
solve its pairing and integral precision layers for A_P. A history product
means concatenation/convolution of histories, not multiplication of their
observed scalar values. Commutativity is a hypothesis, not an analogy.

Fix n,P and M>=1, and allow arbitrary Z/M-linear mixtures of prefix
histories with scalar future observations. Recurrence (2) reduces all
prefixes and all suffixes to the binary cube. A cube coefficient vector
is invisible under every suffix precisely when it lies in ker H.
This kernel is stable under each prefix shift: appending a letter to a
prefix has the same response as prepending it to the suffix, and the
letters commute. Thus the minimal observable linear module is

    (Z/M)^d / ker H  ~= im H,
    im H ~= direct_sum_i Z/(M/gcd(M,h_i)),
    |im H|=product_i M/gcd(M,h_i).                            (12)

The isomorphism is of additive modules, not a diagonalization of every
shift at once. Necessity follows because distinct response functions
must remain distinct in any exact model; sufficiency is the stable kernel
quotient. The integer cube and the unimodular maps supply a constructive
presentation with the legal operators retained.

This answers linear-mixture minimality. Actual individual-source reachable
states form a generally much smaller subset. They cannot be freely replaced
by every element of im H; that would drop the realization constraint.

## D6. Colliding roots, nilpotent memory, and exact actual clocks

Suppose ell is prime, every p_i=1 mod ell, and ell does not divide C.
Then Gamma=C X mod ell is a unit, so H is invertible modulo ell^k for
every k>=1. In characteristic ell the faithful operator algebra is

    F_ell[epsilon_1,...,epsilon_r]/(epsilon_i^2),
    T_i=1+epsilon_i.                                        (13)

Its dimension is 2^r: the cube shifts of the source span because H is
invertible. Thus repeated eigenvalue 1 does NOT mean one-dimensional
memory. The top nilpotent monomial remains observable. Frobenius sends
every element to its constant term, since ell>=2; it is not an automorphism
of this algebra. Each time shift 1+epsilon_i is an automorphism.

The exact order L_p(k) of M_p modulo ell^k is

    ell^k                            if ell is odd,
    2^k                              if ell=2 and p=1 mod4,
    2^max(1,k-v_2(p+1)+1)             if ell=2 and p=3 mod4.   (14)

Indeed (M_p-I)^2=(p-1)(M_p-I), so

    M_p^t=I+(1+p+...+p^(t-1))(M_p-I).

An off-diagonal entry of M_p-I is 1, hence the order is the least t>0
with the geometric sum divisible by ell^k. For odd ell, binomial expansion
of (1+ell^a v)^ell increases the valuation of the difference from 1 by
exactly one, while raising to an exponent prime to ell preserves it.
This gives v_ell(sum)=v_ell(t). For ell=2 and p=1 mod4 the same doubling
argument gives v_2(sum)=v_2(t). For p=3 mod4 an odd t gives an odd sum;
for even t, factoring p^t-1 initially as (p-1)(p+1) and then doubling
gives v_2(sum)=v_2(p+1)+v_2(t)-1. These prove (14).

The smallest deterministic carrier for the actual source's future response
has exactly

    product_i L_(p_i)(k) states.                             (15)

Proof of no hidden joint clock relation: the cube shifts of the initial
vector form the columns of H, hence a basis modulo ell^k. If a product
of shifts fixes that vector, it fixes every cube shift, and is the identity
matrix. A tensor product of invertible matrices over Z/ell^k is the identity
only when each factor is scalar. To see this for two factors A,B, write
A tensor B as blocks a_ij B. Off-diagonal blocks force a_ij=0 by multiplying
by B^-1; a diagonal block equal to I forces a_ii a unit and B scalar.
Induction handles all factors. A power of M_p is scalar only when its
geometric-sum off-diagonal entry vanishes, when it is I. Therefore every
coordinate exponent is divisible by its separate period. All clock states
are reachable, and distinct response vectors are distinguished by a cube
suffix. This proves deterministic minimality, not just matrix order.

For odd ell the actual carrier has ell^(kr) states whereas the linear
mixture module has ell^(k 2^r) elements. These are different questions.
For P=(7,11,13), ell=2, C odd, the actual state counts at k=1,...,5 are
8,16,64,512,4096; the three periods are supplied individually by (14).

## D7. All-order actual hostiles and the 1001 chain

For any finite set P of distinct odd primes and an outside odd prime q,
the actual sources n=product P and q have

    s(n product S)=s(q product S)=1 mod2   for every S proper subset P,
    s(n product P)=0 mod2,   s(q product P)=1 mod2.            (16)

An odd integer has odd sigma exactly when every prime exponent is even.
For a proper S, n product S has a remaining odd exponent; q product S
always retains the odd exponent of q. At S=P, n product P=n^2 is square.
Subtracting the odd source gives (16). This is an all-size family showing
that ALL proper-subset responses can miss the full joint response.

More generally for odd P, sigma(u) is odd iff the odd part of u is square.
If C is odd, D6 gives binary Hankel rank d and d actual states. If C is
even, (1) is the constant -n mod2 for every dilation, giving one actual
state and binary Hankel rank 1 if n is odd, 0 if n is even.

For P=(7,11,13), the pair in (16) is 1001 versus 3. Along the inherited
actual edges 1001 ->343 ->57, the P-free cores are 1,1,57 and their sigma
values 1,1,80. The binary ranks are 8,8,1 and the actual state counts
8,8,1. This is a precise loss mechanism under an actual aliquot step,
not an identification of the three sources' divisor lattices.

At any precision ell^k dividing C, (1) reduces to the single multiplicative
character -n product p_i^t_i. If all p_i are units, the reachable states
are the orbit of -n under the subgroup generated by the p_i modulo ell^k;
its size equals that subgroup's image size modulo ell^(k-min(k,v_ell(n))).
For n=57, P=(7,11,13), this applies through k=4 at ell=2. The generated
unit group is full: 13 generates 1 mod4 by the preceding valuation proof,
and 7 supplies the other mod4 coset. State counts are 1,2,4,8. At k=5,
the C term becomes visible; FINITE-EXACT enumeration gives 32 states.
The general nonunit-C actual orbit classification at all precisions remains
OPEN here, although D4--D5 already determine its linear mixture module.

## D8. Actual aliquot steps do not descend to the dilation quotient

For EVERY finite set of odd dilation primes, 25 and 49 have identical
complete binary dilation-response traces: multiplying an odd square by
m has square status exactly when m does. Nevertheless

    s(25)=6, s(49)=8, s(6)=6=0 mod2, s(8)=7=1 mod2.

Thus a single additional actual aliquot step separates two states identified
by the full dilation observer. No carrier retaining only that equivalence
class can be closed under both operations with the same observation.
The failed implication is extension of a quotient beyond its certified
operation language. This is a decisive obstruction, not a conjecture about
aliquot termination or a lower bound for every richer representation.
