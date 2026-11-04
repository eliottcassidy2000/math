# A positive Seidel energy detects every quotient at a dyadic cover fixed point

**Status: structural proof with exact rational and independent physical/run controls.** For every tournament
quotient Q of order m>=2, a finite dyadic modulus can be computed for which
all sufficiently dyadic transitive child substitutions have a child-uniform
two-port separator through every fixed ambient. There is no upper restriction on
quotient order or primality. Order one is excluded: it is self-converse and
cannot be separated from itself. Arbitrary child profiles are **not** covered.

The proof uses an exact rational sum of squares, not a finite quotient census.
The underlying matrix inverse identities are standard linear algebra, and no
literature-priority claim is made. The point is their combination with the
physical run-word and ordered-cover interfaces.

## 1. Physical observer, profiles, and the geometric functional

For a tournament A, q_k(A) counts ordered spanning covers by k nonempty
directed paths; H(A)=q_1(A) is odd. Quotient children retain disjoint physical
name sets. Ports x,y have edge x->y. Query bits a_i,c_i mean every name of
child i points to x,y, respectively. Target Q^op[B_i] uses the same child
sets and the same complete profiles as source Q[A_i].

The inherited exact cut-aware run identity says

    q_h(Q[A_i]+x+y) = sum_k K_h(Q,a,c;k) product_i q_(k_i)(A_i).

Each child is visited k_i>=1 times, each port once, and h-1 cuts split the
word into ordered paths. An uncut adjacent pair must be an arc and cannot
repeat a child. A cut allows any adjacent symbols, including equal children.
Filling the runs uses one ordered child cover, retaining all physical names.
This is the same interface explicitly proved in Session189's prime-seven
modulo-sixteen packet; it does not use any seven-name specialization.

Consider the formal rational/2-adic profile

    g_k = (-2)^(k-1),       k>=1.

It is not asserted to be a complete profile of any finite tournament.
Applying g to the run identity defines a convergent 2-adic functional:
a term with total run excess sum_i(k_i-1)=e has 2-adic valuation at least e.
For each modulus only finitely many run shapes can contribute. Keep a cut
variable t, so coordinate q_h is the coefficient of t^(h-1).

## 2. Resolvent calculation on each induced subset

Let B be Q's 0/1 adjacency matrix, let S=B-B^T be its skew Seidel matrix,
and write J=11^T. For an induced subset I, use B_I,S_I. The empty-subset two-port walk functional is W_empty(t)=1+2t;
its Hamilton coordinate at t=0 is one.

Give each old run occurrence weight -2. The walk adjacency with cuts is
B_I+tJ: its diagonal t permits a repeated child precisely across a cut.
Set v=1+2t and

    K_I(t) = (I+2B_I+2tJ)^(-1) = (S_I+vJ)^(-1),
    G_I(t) = -2 K_I(t).

These inverses exist as rational functions and as 2-adic formal series in t,
since their constant matrices are I modulo two. They count arbitrary old run
segments, including repeated child visits. The ports themselves occur once.

Write r=K_I(0)1, s=K_I(0)^T1 and h=1^T K_I(0)1. From
K+K^T=2(K1)(K1)^T, one gets h=h^2. Its parity is |I|, so h=0 for even I and
h=1 for odd I. Multiplying the same identity by 1 gives s=-r in the even
case and s=r in the odd case.

For even nonempty I, S_I is invertible: modulo two it is I+J, whose square
is I when |I| is even. Define

    r_I = S_I^(-1)1.

Then 1^T r_I=0 and

    K_I(t)=S_I^(-1)+v r_I r_I^T,
    K_I(t)1=r_I,       K_I(t)^T1=-r_I.                 (1)

For odd I, put r_I=K_I(0)1. Then 1^T r_I=1, S_I r_I=0, and the rank-one
inverse formula gives

    K_I(t)1=K_I(t)^T1=r_I/v,
    K_I(t)+K_I(t)^T=2 r_I r_I^T/v.                    (2)

All coordinates of every r_I have odd denominators. Indeed the defining
integer matrices are invertible modulo two.

Let A=r_I dot a and C=r_I dot c, restricted to I. The total walk functional
using both ports once, allowing all old walks on I, is

    W_I(t) = L_a M_xy R_c + L_c M_yx R_a,             (3)

where, with K=K_I(t),

    L_a=1-2 1^T K(a+t1),
    R_a=1-2((1+t)1-a)^T K1,
    M_xy=1+t-2((1+t)1-a)^T K(c+t1),
    M_yx=t-2((1+t)1-c)^T K(a+t1).

The terms 1+t and t are the direct x-to-y and y-to-x transitions: the first
can be an arc or a cut, while the second requires a cut. Prefix and suffix
old segments may be empty. This retains all port orders and all cuts.

Substitution of (1)-(2) into (3) gives the exact formulas

    W_I(t) = v (1+2A)^2(1+2C)^2,         |I| even,
    W_I(t) = -v^(-3)(1-2A)^2(1-2C)^2,    |I| odd.     (4)

For example, in the even case L_a=R_a=1+2A and
M_xy+M_yx=v(1+2A)(1+2C). In the odd case
L_a=-R_a=(1-2A)/v and
M_xy+M_yx=(1-2A)(1-2C)/v. These verify (4) directly.

Converse replaces K_I by its transpose. It sends r_I to -r_I on even I and
fixes r_I on odd I.

## 3. Inclusion-exclusion, a fixed point, and a Gram kernel

To require every child name at least once, use inclusion-exclusion over I.
The factor (-2)^(-m) converts each child's visit weight (-2)^k to g_k.
Consequently the geometric two-port functional is exactly

    F_Q^g(a,c;t)
       = 2^(-m) [ v sum_(I even) f_I^+(a)f_I^+(c)
                    +v^(-3) sum_(I odd) f_I^-(a)f_I^-(c) ],         (5)

where

    f_I^+(a)=(1+2 r_I dot a)^2,
    f_I^-(a)=(1-2 r_I dot a)^2,

and f_empty=1. Although (5) has a displayed power-of-two denominator,
its Boolean coefficients are 2-adic integers by the original convergent
run functional. Equation (5) is an identity of rational functions embedded
in the 2-adics, not a passage to a real limit of large Hamilton counts.

At t=0, (5) is a positive semidefinite Gram kernel on real query vectors:
it is a positive sum of products of the same squared linear functions.
This real interpretation is used only to establish nonzero rational values.
The approximation to actual child profiles below is 2-adic.

Without ports, the same calculation is particularly simple. The walk
functional is -2*1^T K_I(t)1, which is zero on even I and -2/v on odd I.
Inclusion-exclusion therefore gives

    q_h(Q[g,...,g]) = (-2)^(h-1),       every h>=1.                 (6)

Thus g is an all-order cover-profile fixed point for **every** nonempty
quotient. All scalar ordered-cover coordinates at this formal point forget
the quotient orientation. The two-port observer does not.

## 4. Exact positive energy in the shared named cubics

Let

    Delta_Q^g(a,c;t)=F_Q^g(a,c;t)-F_(Q^op)^g(a,c;t).

The odd-subset terms cancel. For even I,

    (1+2A)^2(1+2C)^2-(1-2A)^2(1-2C)^2
       =8(A+C)+32(A^2 C+A C^2).

Hence

    Delta_Q^g(a,c;t)=(1+2t) Delta_Q^g(a,c;0).          (7)

This is a complete all-cover identity, not merely Hamilton nonvanishing.
The contrast has total Boolean degree at most three. All quartic coordinates
are retained and equal zero, as is also expected from uniform whole-converse
profile duality.

For distinct quotient names i,j define the shared named cubic

    K_ij=[a_i a_j c_i] Delta_Q^g(a,c;0).

Boolean reduction keeps the shared physical child name i across the two
ports. Since a_i,a_j are distinct variables, (5) yields

    K_ij=2^(6-m) sum_(I even, i,j in I) r_(I,i)^2 r_(I,j).           (8)

**Theorem 1 (positive Seidel energy).** For every tournament Q on m>=2 names,

    E(Q):=sum_(i!=j) S_ij K_ij
        =2^(6-m) sum_(nonempty even I) sum_(i in I) r_(I,i)^2
        >0.                                                       (9)

**Proof.** Interchange the finite sums in (8). For each even I,

    sum_(i,j in I) S_ij r_i^2 r_j
        =sum_i r_i^2 (S_I r)_i
        =sum_i r_i^2.

Every term is a nonnegative rational square, and every two-name subset has
a nonzero vector. Thus the sum is strictly positive. QED.

In particular at least one complete named cubic is nonzero. No endpoint
parity, second endpoint invariant, primality, or finite representative list
is needed. This does not prove injectivity between unrelated quotients.

## 5. Actual dyadic transitive children and arbitrary ambient lift

**Lemma 2 (uniform transitive approximation).** If L>=1 and N is any positive
multiple of 2^L, then, for every k>=1,

    q_k(TT_N) = k! S(N,k) = (-2)^(k-1) (mod 2^L).                  (10)

**Proof.** Use the exact surjection formula

    k!S(N,k)=sum_(j=0)^k (-1)^(k-j) binom(k,j) j^N.

Every even j contributes zero modulo 2^L because N>=2^L>=L. Every odd j
satisfies j^N=1 modulo 2^L: for L>=2 this follows from the elementary unit
group exponent bound dividing 2^(L-1), and L=1 is immediate. The remaining
sum over odd j is (-1)^(k-1)2^(k-1). This also covers k>N, where both sides
vanish modulo the stated modulus. QED.

Because (9) is a nonzero rational 2-adic integer, its valuation v_2(E(Q))
is a finite nonnegative integer. Compute

    L(Q)=1+v_2(E(Q)).

This is a constructive finite threshold: enumerate the nonempty even subsets,
solve each integer Seidel system exactly over the rationals, sum (9), and
reduce E(Q)=p/q to lowest terms. The reduced denominator q is odd and p>0.
Then L(Q) is one plus the number of factors of two dividing p. No existence
oracle or unspecified asymptotic threshold is needed; the construction is
not claimed to run in polynomial time.

**Theorem 3 (all-order dyadic-transitive noncancellation).** Let Q have m>=2
names. For each i take TT_(N_i) with N_i a positive multiple of 2^L(Q).
More generally, allow any actual child A_i whose complete normalized profile
satisfies q_k(A_i)/H(A_i)=g_k modulo 2^L(Q) for all k. Compare Q[A_i] with
Q^op[B_i], where each corresponding complete profile is preserved.
Then a child-uniform physical two-port query separates them through every
fixed ambient root substitution.

A separating query can be selected from the union, over i!=j, of the eight
corners on variables {a_i,a_j,c_i}, all other query bits zero. This finite
battery is independent of child sizes and of the ambient.

**Proof.** Let Lambda=product_i H(A_i), an odd common factor. The run identity
and profile congruences give the actual complete profile contrast

    Delta q_1=Lambda Delta_Q^g(a,c;0),
    Delta q_2=2 Lambda Delta_Q^g(a,c;0),
    Delta q_h=0 for h>=3                         (mod 2^L(Q)).      (11)

The convergence is uniform in the needed cover coordinates: modulo 2^L only
run shapes with total excess below L survive, hence there are finitely many
words and cuts to check. Any larger shape has a profile-product factor
2^L. Thus no infinite unsupported summation or unbounded ambient tail is
being discarded.

The actual cubic energy is Lambda E(Q) modulo 2^L(Q), which is nonzero.
Therefore some K_ij is nonzero modulo this modulus. Its eight-term finite
difference cannot be nonzero if all eight actual Hamilton contrasts vanish.
Choose a corner with nonzero contrast delta.

Lift both ports with the core root's incidences to all exterior names. The
augmented core remains a common module. The ambient Hamilton difference is
sum_h c_h Delta q_h, with nonnegative integer coefficients determined by
the fixed ambient. Here c_1=H(R) and c_1+2c_2=H(R[TT2]), which is odd.
By (11), the lifted contrast is (c_1+2c_2)delta modulo 2^L(Q), and so remains
nonzero. This works for every fixed ambient and the same selected core query.
QED.

L(Q) is a sufficient modulus, not a sharp child-size threshold. Individual
nonzero cubic coefficients can give smaller sufficient moduli. The theorem
also applies to unequal child orders satisfying (10). Arbitrary small or
unmatched child profiles are not covered, and no conclusion about every
remaining prime-eight original-module case is inferred.

## 6. Two structural extensions of the actual profile class

**Substitution closure.** For a fixed L, the class of actual tournaments with
normalized complete profile g modulo 2^L is closed under arbitrary quotient
substitution. Indeed (6) and multilinearity give q_h(Q[A_i])=Lambda*g_h
modulo 2^L and H(Q[A_i])=Lambda modulo 2^L. Dividing by the latter odd unit
restores g_h. This covers recursively substituted dyadic transitive building
blocks, not just a literal transitive child. It does not assert that every
tournament belongs to this class.

**Profile-dominating targets.** Theorem 3 also gives a fixed source-relative
separator against every Q^op[B_i] with q_h(B_i)>=q_h(A_i) over the integers,
once a baseline query is paired with J(a,c)=(1-c,1-a). Complete augmented
profile duality gives opposite baseline integer contrasts at z and Jz,
including after the stated ambient lift. Choose the negative one. Uniform
run substitution has nonnegative coefficients, so raising target child
profiles preserves its strict negative sign. A battery containing both
members works without first knowing the baseline sign.

In particular, when A_i=TT_(N_i), the target interiors B_i may be arbitrary
tournaments on the same physical child sets. For every tournament B on N
names, the exact unordered-partition formula and H(B[S])>=1 give

    q_h(B) = h! sum_(unordered h-partitions P) product_(S in P) H(B[S])
           >= h! S(N,h) = q_h(TT_N).

Thus the transitive-source result does not require profile-preserving target
interiors. For a prime quotient it includes source-relative original-module
families containing that quotient root oddly: proper original modules lie
inside children. It does not cover arbitrary root-family quotient changes
when Q is composite, nor imply pairwise reconstruction of arbitrary outputs.

## 7. A fixed-modulus theorem for every odd carousel

Let n=2k+1>=3 and let Q be the carousel with i->i+1,...,i+k modulo n.
Use the original physical quotient labels 0,...,n-1. Then

    [a_0 a_k c_0] Delta_Q^g = -8(n-4).                (12)

This is a closed all-order identity, not an extrapolation from small n.
Since n-4 is odd, its residue is always eight modulo sixteen.

**Proof.** Switch the Seidel signs on the in-neighbors k+1,...,2k of vertex0.
The switched tournament is transitive in the order

    0, k+1, 1, k+2, 2, ..., 2k, k.

Indeed each of the two original blocks is transitively ordered; after the
switch, in-name k+s precedes out-name r exactly when s<=r. Vertex0 precedes
all of them. In this displayed order the original Seidel matrix is D T D,
where T is the transitive skew matrix (positive above its diagonal) and the
signs of D alternate. A global sign change in D does not change D T D, so
index positions by 1,...,n and use d_j=(-1)^j.

For an even selected subset with ordered positions j_1<...<j_(2r), the
transitive skew inverse is

    (T_(2r)^(-1))_(ab)=(-1)^(a+b) (T_(2r))_(ab),

as direct multiplication and the finite alternating sums verify. Put
w_a=(-1)^a d_(j_a). The induced inverse vector therefore has coordinate

    r_(j_a)=w_a (sum_(b>a) w_b - sum_(b<a) w_b).       (13)

Sum (13)'s cubic r_1^2 r_n over all even subsets containing positions1 and n,
which are original labels0 and k. Encode any subset by

    u_l=(-1)^l (-1)^(number of selected positions <=l),  u_0=1.

Even cardinality forces u_n=-1. Selecting positions1 and n forces u_1=1
and u_(n-1)=-1. The remaining u_2,...,u_(n-2) are arbitrary independent signs,
giving a bijection with the 2^(n-3) allowed subsets. At a selected position,
w_a=u_(j_a); more generally its contribution to the sum of w's is
(u_l+u_(l-1))/2. Thus the total sum of w's is

    Z = u_2+...+u_(n-2),

and (13) gives r_1=Z-1 and r_n=Z+1. Averaging over the independent signs,
E[Z]=E[Z^3]=0 and E[Z^2]=n-3. Consequently

    sum_(even I containing 1,n) r_1^2 r_n
       =2^(n-3) E[(Z-1)^2(Z+1)]
       =2^(n-3)(4-n).

Multiplication by the factor 2^(6-n) from (8) proves (12), including n=3
where Z is the empty sum. QED.

**Theorem 4 (uniform carousel family).** For every odd carousel order n>=3,
if all actual child normalized profiles equal g modulo sixteen, a query in
the eight-corner cube {a_0,a_k,c_0} separates the quotient converse in every
fixed ambient. In particular this holds for arbitrary unequal transitive
child orders divisible by four, and for recursively substituted children
in the same profile class. The target profiles may be preserved or dominated
as in Section6; for a transitive source, arbitrary target child interiors
are allowed after pairing the cube with its converse-complement cube.

To check the transitive-child assertion, every odd integer j has j^4=1
modulo sixteen, while every even j has j^N=0 modulo sixteen for positive
N divisible by four. The surjection proof of (10) therefore gives g at that
fixed modulus. The complete contrast relation (7) supplies the arbitrary
ambient lift; it is not inferred from the cubic alone. This theorem imposes
a profile condition and is not an all-profile carousel theorem.

## 8. Scope, retained data, and reproduction

The structural proofs are the exact resolvent identities, inclusion-exclusion,
positive energy, and transitive profile congruence. The producer verifier
checks all 75 labeled quotients through order four, 25,188 direct resolvent
identities, 885 transitive-profile congruences, and the independent-sign
carousel controls at every odd order from three through nineteen.

`rational_banks.json` retains nine complete quotient-named block-uniform
coefficient banks through carousel order eleven: 10,098 coordinates on each
side, including zero, cubic and quartic coordinates, plus the common odd-subset
bank and every rational induced-subset vector. The common odd-subset bank
alone can have powers of two in its denominators; the complete response and
cover coefficients have odd reduced denominators. All numerator/denominator
pairs are exact. Complete eight-corner source and target values are retained
for each named example. The contrast's quartic zero entries are not omitted.

These are complete *block-uniform* banks, not full physical banks for the
much larger inflated tournaments. Independently written audit packets retain
separate literal physical named banks, full query vectors, direct cut-aware
run counts and ambient-word controls. They distinguish those literal counts
from the rational geometric functional. The geometric profile is never
represented as an actual finite tournament.

Run `python3 verify.py > normal.json`, `python3 -O verify.py > optimized.json`,
and compare the two receipts. The verifier uses only the Python standard
library, keeps checks active under optimization, and deterministically
regenerates its three data files. Frozen producer files are PROOF.md,
README.md, verify.py, rational_banks.json, small_energy_bank.json,
carousel_controls.json, normal.json, and optimized.json. Exploratory code is
not a dependency.

This result supplies an unbounded noncancellation family and an exact scalar
fixed-point boundary. It is not a proof for arbitrary child profiles, an
exact response collision, unrestricted reconstruction, or an asymptotic
circuit lower bound. Order one is not a noncancellation case. The subsequent
source-thickening application, if included, has its own stated integration
inputs and review; it is not used to prove this packet.
