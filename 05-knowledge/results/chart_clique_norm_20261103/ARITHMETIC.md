# Cubic trace fibres, an exact quartic map, and one ternary precision digit

**PROVED elementary arithmetic and algebra / CITED Eisenstein factorization /
FINITE-EXACT controls.** Local labels A1--A4 refer to this file. The field
realization of the nine q=5 traces is inherited, not newly established here.

## A1. Complete integer collision classification, for every positive q

Define F_q(t)=q^3+1-t^3+3qt for integer t and positive integer q. For distinct
u,t,

    F_q(u)=F_q(t) iff u^2+ut+t^2=3q.                       (1)

Reduction modulo3 forces u=t mod3. Write u=t+3a. Then

    q=t^2+3at+3a^2=N(t+a-a omega), omega^2+omega+1=0.       (2)

Conversely, a norm representation q=N(b+c omega) gives
(t,u)=(b+c,b-2c). If c=0, multiplication by the unit omega supplies a
representation with c!=0, since q>0. Therefore F_q is injective on Z
**if and only if q is not an Eisenstein norm**.

For completeness, the classical norm criterion says q is a norm iff every
prime congruent to2 modulo3 has even exponent in q. The count of elements
of norm q is

    r_E(q)=6 A(q),   A(q)=sum_(d|q) chi_3(d),              (3)

where chi_3 is0,1,-1 on residues0,1,2. Equivalently A=0 if an inert prime
has odd exponent; otherwise A=product_(p=1 mod3)(v_p(q)+1).
This follows from unique factorization in Z[omega]: six units; a unique
ramified prime above3; split choices e+1 above each p=1 mod3; and an inert
prime with norm p^2 above p=2 mod3. The standard prime classification and
Euclidean property are cited in Gullerud--Mbirika, *An Euler Phi Function
for the Eisenstein Integers and Some Applications*, section2 and
Proposition2.9, [primary paper](https://math.colgate.edu/~integers/u20/u20.pdf).
The application to F_q and the fibre count below are derived here.

Let f_j(q) be the number of integer target values with exactly j distinct
integer preimages under F_q. Then the complete nonsingleton-fibre census is

    f_2(q)=2 [q is a square],
    f_3(q)=A(q)-[q is a square].                          (4)

Proof: two distinct integer roots of the cubic leave the third root
-u-t integral. It coincides with one of the first two exactly at the
critical roots +/-sqrt(q). For square q these produce precisely two
two-point fibres, {s,-2s} and {-s,2s}, where s=sqrt(q)>0. Equation(2)
counts ordered pairs including the diagonal; exactly two diagonal pairs
exist when q is square. Thus

    r_E(q)-2[q square]=6 f_3(q)+2 f_2(q),

which gives(4). No fibre has more than three points. Every repeated-fibre
trace satisfies |t|<=2sqrt(q), since
u^2+ut+t^2=(u+t/2)^2+3t^2/4. All traces outside this interval are singletons.
Consequently a finite search inside the Hasse-sized interval is exhaustive
for integer collisions, even when q is not a prime power.

### The actual role of1729

All three primes7,13,19 split in the Eisenstein ring. Hence

    A(1729)=2^3=8, r_E(1729)=48, f_3(1729)=8, f_2(1729)=0.

The eight trace triples are

    {-83,37,46}, {-82,29,53}, {-79,17,62}, {-73,2,71},
    {-71,-2,73}, {-62,-17,79}, {-53,-29,82}, {-46,-37,83}.

They come from a specified norm equation, not from matching the eight
non-pair prime labels in PRIMEGAME. Here1729 is the parameter of an integer
polynomial, not a finite-field order. Its Carmichael property and its two
taxicab decompositions are separate inherited facts; neither is used in(4).

For a prime-power q=p^e the arithmetic collision criterion becomes:
collisions occur for p=3, for p=1 mod3, and for p=2 mod3 with e even.
Which colliding trace states are realized by elliptic curves must still be
checked separately. No assertion that all Hasse-allowed traces occur over
all finite fields is made.

## A2. A complete ternary precision law

For all integers q,t,

    F_q(t+3u)=F_q(t)+9[(q-t^2)u-3tu^2-3u^3].             (5)

Therefore F_q modulo9 factors through t modulo3. Its three values are
distinct, because F_q(t)=q+1-t modulo3. Thus the map on Z/9 has exactly
three fibres of size3, for every q.

When q=2 mod3, q-t^2 is a unit modulo3 for every t. For unequal integers
t,u the exact valuation is

    v_3(F_q(t)-F_q(u))=0,             if t!=u mod3;
                         v_3(t-u)+1, otherwise.           (6)

To prove the second line, write u=t+h with3|h and factor the difference
as3h[q-t^2-th-h^2/3]. The bracket is a unit. On a disk t=r+3z, equation(5)
is F_q(r)+9G_r(z), where

    G_r(z)=(q-r^2)z-3rz^2-3z^3.

Its derivative is a unit modulo3, and its reduction is a nonzero linear
map. Elementary digit-by-digit Hensel lifting therefore makes G_r a
bijection of Z_3. The entire image is exactly the union of the three
residue disks F_q(r)+9Z_3, r=0,1,2.

For every m>=1, this gives a bijection

    Z/3^m -> {y mod3^(m+1): y mod9 is F_q(r) for some r}.  (7)

For m>=2, output only modulo3^m loses one ternary source digit: each image
has three lifts. This precision loss is sharp. The positive construction
is invert.py: check the modulo9 image guard, then lift one digit at a time.
It is an arithmetic inverse, not a curve-construction algorithm.

The q restriction matters. F_1(-2)=F_1(1), and F_3(0)=F_3(3), so extending
the unit-derivative inverse to every q would be false. The complete global
integer criterion is A1; for example q=55 is injective on Z although q=1
mod3, since both inert primes5 and11 occur to odd exponent.

## A3. The supplied magic square and the older digit clock

Session143 already proves and independently checks

    #E(F125)=F_5(t)=126-t^3+15t,  t=-4,...,4,

and realizes all nine traces over F5. The supplied square

    108 148 122
    140 126 112
    130 104 144

is exactly its inherited arrangement; all rows, columns and diagonals sum
to378. It is not a new square and not a square of integer squares.

The recovered October17 Solfeggio note encodes the nine integers as

    S(a)=100a+10 d(a+6)+d(a+3), a=1,...,9,
    d(x)=1+((x-1) mod9).

Their digit rotations are the three fibres of Z/9 -> Z/3. Define the
explicit index map S(a) -> t=a-5 -> F_5(t). Digit rotation corresponds to
a shift by3 or6 modulo9, so A2 sends each digit-rotation orbit to one
residue class of the point counts modulo9:

| Solfeggio orbit | point-count orbit as a set | residue mod9 |
|---|---|---:|
|174,417,741|130,112,148|4|
|285,528,852|108,126,144|0|
|396,639,963|104,140,122|5|

Modulo27, all nine point counts are distinct and recover the nine trace
states. The difference between mod9 and mod27 is the missing ternary phase.
This map preserves the indicated quotient and index rotation. It does not
realize that rotation by an elliptic-curve automorphism, and does not assign
physical or therapeutic meaning to the integers called frequencies.
The original nine Solfeggio integers cannot themselves make a3-by3 magic
square: their sum is4995, so the required central entry555 is absent.

Incoming Session145's `portfolio/magic_frequency/PROOF.md` independently
recovers this obstruction and the trace-versus-order-map distinction, and
classifies all real magic maps in its prescribed placement. We credit that
overlap. The additional operation here is A2's all-precision arithmetic
inverse and its explicit modulo27 recovery of the digit-rotation phase.

## A4. A literal map into the inherited quartic discriminant

The repo's weighted quartic has inverse polynomial

    T(w)=w^4-w^2+2Pw-Q, P=BC, Q=AC^2.

Off C=0 its finite inverse points correspond to simple roots; its repeated
root a has gamma=0 and is excluded from that reconstruction chart. The
discriminant is parametrized by

    P=a-2a^3, Q=a^2-3a^4,
    T(w)=(w-a)^2(w^2+2aw+3a^2-1).                        (8)

These quartic mechanisms are inherited from THM-3441 and the current
Session144/prime-wait audits. The new explicit arithmetic substitution is

    a=t/sqrt(6q),
    P=(F_q(t)-(q^3+1))/(3q sqrt(6q)),
    Q=(2q t^2-t^4)/(12q^2),   (A,B,C)=(Q,P,1).           (9)

For q=5 this sends the nine actual trace-count states to the quartic's
discriminant curve over Q(sqrt30). Every one has exactly two finite
preimages: the two roots of the second factor in(8). Exceptions would
require t^2=q (a cusp, one simple root) or t^2=3q (the node, no simple
roots), and neither5 nor15 is an integer square.

At general q, the vertical P-coordinate collisions in this normalization
are exactly A1. The retained Q-coordinate separates distinct such integer
parameters except the node pair u=-t, t^2=3q. Indeed Q(u)=Q(t) implies
(u^2-t^2)(2q-u^2-t^2)=0; combine with(1). For q=1729 the eight triples of
equal P are therefore separated by their three different Q values.

This is an actual polynomial parametrization and a precise missing
coordinate, not a comparison of equal small cardinalities. It does not
identify the cubic integer fibres with the quartic's three-point boundary.
On C=0 the latter still requires the separate x=0 and gamma=0 charts:
one point on x=0, two more when B^2-A!=0, and an escaping generic sheet.
The inherited target(0,1,0) and its three rational inverse points are
independently replayed here. Incoming5e14c36ccf already extends the weighted
family to general inertia partitions; that theorem is credited, not redone.
