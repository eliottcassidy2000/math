# Finite centroid addresses, rational cycles, and a local discriminant11 return

2026-10-05. **PROVED:** finite-centroid word recovery, exact
denominator/depth, total rational membership, the displayed period3 cycles,
and the exact lattice/cubic/local-polynomial identities. **FINITE-EXACT:**
the declared finite graph censuses and independent controls.
**INHERITED PROVED:** marked fan-matrix decoding and the controller
serialization framework. No universal Collatz coverage, global modular
identification, or conductor derivation is claimed.

The denominator9 cycle is real, but it cycles a rational **point**.
The inherited finite decoder consumes a marked **matrix** with a
determinant-specified depth. A useful stronger result separates them:
the special point A_w(1,1,1)/3 itself recovers every finite fan word w.
Its reduced denominator is exactly3^(length(w)+1). General rational
points need not belong to this finite-centroid language.

There is also a concrete arithmetic bridge. The cycle's transverse
return has polynomial z²+8z+27 and discriminant−44. A signed cubic
root, on an explicitly retained index2 lattice, has polynomial z²+z+3,
the local polynomial obtained from the level11 eta product's coefficient
a3=−1. This is an exact matrix identity with stated sign and lattice,
not an identification of the full dynamics.

Artifacts: [script](../../04-computation/experiments/centroid_membership_20261005.py)
and [saved output](centroid_membership_20261005.out).

## 1. Inheritance and the missing stopping coordinate

The closest proved mechanism is
[small-port refinement](small_port_refinement_20261004.md), sections2--4.
For c in{0,1,2}, set a=c+1,b=c+2 modulo3 and

    A_c=[e_a,e_b,(1,1,1)/3],
    T_c=A_c^(-1),
    T_c(x)=(x_a-x_c,x_b-x_c,3x_c).

The chart image has x_c minimal; its interior has a unique minimum.
Words are written outside first:
A_(c0...cd)=A_c0...A_cd. The inherited matrix decoder obtains d from
|det A_w|=3^(-d), strips exactly d charts by the image centroid's
unique minimum, then requires a terminal permutation matrix.
It does not iterate an arbitrary rational point until something
unproved happens.

The canonical hostile below has denominator9 and never reaches the
center. The corrected near miss is treating rationality, or even a
power-of-three denominator, as proof of finite address membership.
The least-used sidecar is the terminal marker and its exact denominator.
Our board is **word / marker / denominator / point orbit / transverse
operator / lattice / controller guard**.

Use fan charts here. The six barycentric flag charts from the same note
have a different stripping map and different periodic points. Chart
names cannot be omitted from a period statement.

## 2. A single chosen point retains every finite fan word

Let o=(1,1,1)/3 and define E(w)=A_w o.

**Finite-centroid theorem.** For every finite word w of length d:

1. E(w) has reduced denominator exactly3^(d+1).
2. Its three primitive positive numerators are each1 modulo3.
3. Unique-minimum stripping recovers w, stopping at o after exactly d
   steps. In particular E is injective, including across different depths.

**Proof.** Start with the primitive triple(1,1,1), denominator3. Applying
an outer chart replaces a numerator triple(u,v,z), up to the chart's
coordinate placement, by

    (3u+z,3v+z,z),

and multiplies the denominator by3. If all input numerators are1 mod3,
all output numerators are1 mod3. Since the denominator is a power of3,
the point is still reduced. Induction proves statements1--2.
The chart slot is the unique minimum: the other two numerators exceed
it by positive multiples of3. Stripping recovers the inner point.
The center has no unique minimum, so it cannot masquerade as a
nonempty address. Denominator and successive slots prove statement3.

This yields a short total membership algorithm for any exact rational
point in the closed simplex. Reject a zero coordinate. Normalize to
primitive numerators and denominator D. Unless D=3^k with k>=1, reject.
Strip exactly k−1 unique-minimum charts, requiring each normalized
denominator to drop by3, and accept exactly when the remaining primitive
triple is(1,1,1). Reversing the accepted steps proves sufficiency, not only
necessity. The input center is the empty word and is accepted before
any tie rejection.

There is an exact counting consequence. At reduced denominator D=3^k,
the language has precisely3^(k−1)=D/3 points. The total number of
primitive positive triples summing to D is

    binom(D−1,2)−binom(D/3−1,2)
      =4*3^(2k−2)−3^k.

For k=1 the second binomial term is interpreted as zero. The language
therefore occupies exactly3/(4D−9) of that finite primitive grid.
This is a fixed-denominator counting statement, not a Lebesgue-density
or source-coverage statement.

The terminal port frame is still lost: E(w)=A_w P o for every permutation
matrix P. The point recovers w and hence its canonical chart matrix;
P must remain a separate label if geometric interface composition needs
it. The inherited marked matrix recovers both.

## 3. A total rational orbit classifier and its exact boundaries

There is also an independent finite-state decision procedure, useful for
understanding nonmembers. Write a normalized rational point as
(u,v,z)/D, with primitive nonnegative integer numerators summing to D.
Keep this initial D fixed during stripping.

If a coordinate is zero, reject as a boundary point. If all coordinates
are equal, accept the center. If the minimum occurs twice but the point
is not the center, stop at a noncentral tie. Otherwise strip the unique
minimum. Its integer numerator triple is still positive and still sums
to D. There are only binom(D−1,2) such positive triples, so a trajectory
must reach one of those stopping states or repeat a state. At most that
many transitions and one final inspection suffice.

A repeated state determines an exact eventually periodic orbit: retain
the initial state, preperiod, cyclic states, and their slot labels.
It is not an accepted finite word. Conversely a center hit reconstructs
a finite word by reversing the actual invertible chart steps, proving
agreement with the denominator decoder.

Unique-minimum stripping from a strictly positive point remains strictly
positive. Thus boundary zeros do not appear unnoticed along this branch;
they are rejected at entry. Noncentral ties are stopped rather than
resolved by an arbitrary orientation choice. If one wanted a boundary
address convention, that would be a different, explicitly labelled language.

There is a useful arithmetic constraint. For a primitive input triple,
the gcd of the stripped numerators is either1 or3. Indeed a common
divisor g divides3u,3v,3z, since it divides both differences and3 times
the minimum. Primitivity gives g|3. Consequently the reduced denominator
can lose only one factor3 at each step. It never gains a prime.
This explains the power-of-three necessity for center membership, but
the next example disproves sufficiency.

## 4. Denominator9: two exact period3 cycles

All triples in this display have denominator9:

    (1,2,6) --0-->(1,5,3) --0-->(4,2,3) --1-->(1,2,6),
    (2,1,6) --1-->(5,1,3) --1-->(2,4,3) --0-->(2,1,6).

Every minimum is unique and all three states in each orbit differ.
The primitive periods are3 and the slot words are001 and110.
The denominator9 point(1,4,4)/9 instead strips once to the center;
(1,1,7)/9 stops at a noncentral tie.

The complete primitive denominator9 census has27 points:
three center members, six eventually tied points, and eighteen points
eventually entering the two displayed cycles. There are six periodic
states themselves. Exhaustion of smaller denominators3--8 shows that9
is the smallest denominator supporting a strict period3 cycle.
It is not the smallest nonterminating example:
(1,2,3)/6 and(2,1,3)/6 are fixed points of slots0 and1.
All these minimality statements refer to this explicitly defined map
on positive rational points.

The first cycle is also an all-length history collision. Put
p=(1,2,6)/9 and W=A_0 A_0 A_1. Then Wp=p, so W^j p=p for every j>=0.
Even knowing this noncentral terminal marker does not recover how many
blocks were applied. The marked matrices remain different:
det(W^j)=3^(-3j). The chosen center marker avoids this collision through
its strict denominator/depth law.

## 5. The period3 return carries a genuine discriminant11 operator

The chronological stripping return for001 is M=T_1 T_0 T_0:

    M=[[-7, 4, 0],
       [-4, 0, 1],
       [12,-3, 0]].

It fixes(1,2,6), preserves the sum of coordinates, and has characteristic
polynomial

    (z−1)(z²+8z+27).

On the integral sum-zero plane, use basis
u=e_0−e_2, v=e_1−e_2. Its matrix is

    L=[[-7,4],[-5,-1]],
    tr(L)=−8, det(L)=27.

The two eigenvalues are−4±sqrt(−11); the quadratic discriminant is−44.
The sum-zero plane is the transverse linear space, not a reduction
modulo9 and not the three-state orbit itself.

There is an exact cubic root

    C=−(L+3I)/2=[[2,−2],[5/2,−1]],
    C²−C+3I=0,   C³=L.

Equivalently pi=(1+sqrt(−11))/2 satisfies
pi²−pi+3=0 and pi³=−4−sqrt(−11). Its norm is3; the return norm is27.
The ring inclusion Z[pi³] inside Z[pi] has index2, directly from
pi³=−2pi−3. No inference about an elliptic conductor follows.

The half entry in C matters. Replace the original lattice by
Lambda'=Z(2u)+Zv, an index2 sublattice. It is stable, and in that basis

    C'=[[2,−1],[5,−1]],
    L'=[[-7,2],[-10,−1]]=(C')³.

This is an actual lattice change; the original basis transformation
diag(2,1) is not unimodular. Let D=−C'. Then

    D=[[-2,1],[-5,1]],  tr(D)=−1, det(D)=3,
    D³=−L',  det(I−tD)=1+t+3t².

After that index2 change there is an integral companion conjugacy:

    S=[[0,1],[1,2]], det(S)=−1,
    D*S=S*[[-1,−3],[1,0]].

The combined basis change from the original tangent lattice has
determinant−2. Thus both the sign and the lattice are retained in
identifying the local operator.

For the formal level11 eta product

    f(q)=q*prod_(n>=1)(1−q^n)²(1−q^(11n))²,

the first coefficients are1,−2,−1: after removing q, degrees at most2
come from(1−q)²(1−q²)². Hence a3=−1 directly, and the polynomial
z²−a3*z+3 is exactly D's polynomial. The interpretation as an elliptic
or modular local factor belongs to the separate level11 result; the
equality here requires only the displayed formal coefficient and matrices.
One local polynomial does not determine a global modular form,
its conductor, an orbit language, or a Collatz guard.

An additional reduction boundary is visible without approximation:
D mod3 is[[1,1],[1,1]], singular of rank1. The golden F9 clock from
[the5/8/9 transfer note](five_eight_nine_transfer_20261005.md) is invertible.
These cannot be conjugate as operators over F3. A different prime's
local operator may have that golden reduction; its prime label must
remain attached.

There is a controlled connection between finite and periodic addresses.
For W=A_0 A_0 A_1, its transverse action is L^(-1), with eigenvalues
(-4±sqrt(−11))/27, each of absolute value1/sqrt(27).
The exact rational quadratic form

    Q=[[5,−3],[−3,4]],   det(Q)=11,
    L^T Q L=27Q

gives a direct proof without approximate eigenvalues. If
W^j o−p=alpha_j*u+beta_j*v, then

    5alpha_j²−6alpha_j beta_j+4beta_j²=4/27^(j+1).

Positive definiteness therefore implies W^j o converges to
p=(1,2,6)/9. Every term is a finite centroid with exact denominator
3^(3j+1), but its limit is a nonmember. The finite-centroid language
is not closed in the real simplex. This is a geometric error measure,
not a paid rank or a Collatz convergence proof.

The contraction is specific to this block. Arbitrary infinite fan words
need not determine a unique point from their nested cells: A_2 fixes both
e_0 and e_1, so every A_2^j triangle retains their whole edge.
Finite word recovery, convergence of a selected sequence of marked
points, and shrinking of the complete nested cells are three distinct
claims.

## 6. Traces, coefficients, and controller expressions retain different data

For a two-dimensional matrix with trace t and determinant d,

    tr(D²)=t²−2d,     tr(D³)=t³−3dt.

For the displayed D, its traces begin2,−1,−5,8 and satisfy
s_(j+2)=−s_(j+1)−3s_j. Thus the geometric return has transverse trace−8.
In contrast the coefficients of1/(1+t+3t²) begin1,−1,−2,5.
They satisfy the same recurrence with different initial values; for j>=2,
s_j=h_j−3h_(j−2). In particular the period3 return trace is not the
eta product's coefficient at27. The operator, its trace, and its
reciprocal-factor coefficients are distinct observables.
The script independently truncates the formal eta product through q^27
and verifies a27=5, so this displayed coefficient comparison does not
depend on importing a Hecke theorem.

There is a simple lossless controller serialization with an explicit
boundary. Fix, for example, H->0,G->1,L->2. Every finite H/G/L expression
can be encoded as its finite fan centroid and recovered by section2.
Once recovered, the exact controller compiler can rebuild the expression's
native source guard and affine action. The supplied source, immutable
comparison source, credit state, and terminal root proof remain required
inputs to arithmetic verification. Geometry supplies no paid guard merely
by accepting a finite word.

Likewise a compact periodic address retains an infinite point orbit, not
a completed finite controller episode. The cycle001 does not assert that
the corresponding H/H/G expression is legal forever on a supplied integer.
The inherited credit proof already excludes infinite legal funded episodes;
there is no conflict, because the maps and retained predicates differ.

## 7. Reproduction and finite universe

    python -B 04-computation/experiments/centroid_membership_20261005.py
    python -B -O 04-computation/experiments/centroid_membership_20261005.py

The script independently checks all1093 words through depth6, all726
terminal-frame variants through depth4 against the inherited matrix
decoder, complete membership counts at exact denominators3,9,27,81,
and all712 primitive positive triples with denominator3 through18.
It includes the two denominator9 cycles, fixed-point and tie hostiles,
boundary points, invalid float/bool input, five finite controls of the
all-length point collision, exact return/lattice/conjugacy matrices,
and trace versus reciprocal-factor recurrences through index16.
The controlled infinite-limit proof has seven exact checks at j=0,...,6
of the denominator, quadratic error and retained-edge hostile.

Every infinite statement has a proof above; the finite counts are marked
separately. Checks use exceptions and remain active under optimization.
