# Height-bounded source circuits and synchronized local lifting

Session 200. **PROVED conditional on the frozen physical C1 source/adjoint
classification:** an exact quadratic system of size O(log² n) recognizes a
C1 target over every characteristic-zero field and over any single prime
field of characteristic greater than 4n. A concrete CRT hostile shows why
solutions at separate small primes cannot be combined without synchronizing
the physical word. **OPEN:** all-target source coverage and the {7,21}
spectrum conjecture. This is a target-dependent family of finite systems,
not one fixed polynomial defining an unbounded integer predicate.

## 1. Inheritance and connection contract

The frozen dependency `dependencies/PHYSICAL_ADJOINT.md` supplies the complete
physical source grammar, the pairing K, its physical interpretation, and its
midpoint inverse. It is the Session 199 physical-adjoint packet. The nearest
hostile is the earlier profinite phantom (1,0,6,1): congruence approximants can
escape to unbounded ordinary height. The repaired near miss is bounded-height
compactness. The least-used sidecar is that an exact source word has *Boolean
choices* shared by every specialization.

The live concepts are physical source words, ordinary height, finite-field
specialization, witness synchronization, and the distinction between a field
and a product of fields. The number-field idea used here is to build an
algebraic witness whose equations themselves force its discrete provenance,
instead of assuming integrality of every auxiliary field element.

Source -> target: a pair of actual bounded source words maps to a Boolean
quadratic circuit, then to its reduction modulo a prime. The exact count is
preserved when both the height bounds and shared word choices are retained.
Reduction loses ordinary magnitude; separate reductions can additionally
lose the identity of the existential witness. Binary range sidecars restore
the first loss, while a common selector string restores the second.

## 2. Physical facts and a logarithmic word cap

The nontrivial seed is t=(1,0,2,1), in raw order (D,E,U,d). The moves are

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

Every source is such a word, and every word is physical. Its coordinates
satisfy 0<=E<=d<=D<=U. The cover U strictly increases under each letter.
For two sources let

    K(v,w)=D_v d_w+d_v D_w+E_v E_w+U_v U_w
           -E_v d_w-d_v E_w-U_v d_w-d_v U_w+d_v d_w.

The complete strong C1 scalar class is K(S,S). Moreover

    U_v U_w/2+d_v d_w <= K(v,w) <= U_v U_w,

and every represented n>=3 has a pair with both covers at most

    B(n)=floor(2 sqrt(n-1)).                            (1)

Keep the physical domain: arbitrary integral cone points are not sources.

For a self-contained word cap put G=D+2U and use positive coordinates
(a,f,e,q)=(D-E-d,U-D+E-d,E,d-E). Thus
G=3a+2f+6e+5q. The nonnegative coefficient vectors of
2G(letter(s))-3G(s), in this order, are

    O: (1,2,2,3), Z: (3,4,0,3), J: (3,0,4,3).

Consequently G grows by at least 3/2 at every actual letter. Its initial
value is 5; at cover<=B it is at most 3B. Define L=L(B) as the largest
nonnegative integer with

    5*3^L <= 3B*2^L.                                  (2)

Every actual word ending at cover<=B has length<=L. In particular L=O(log B).
An identity letter I pads shorter words to exactly L. This identity is only
padding; it adds no graph vertex. Since cover increases, every prefix of a
bounded word has all four raw coordinates in [0,B].

## 3. The height-enriched quadratic system

Fix integer n>=3, B from (1), L from (2), and ell=ceil(log2(B+1)). For each
of two traces, each time 0<=i<=L, and each raw coordinate j, introduce two
ell-bit strings b_(i,j,h), c_(i,j,h). Write

    X_(i,j)=sum_(h=0)^(ell-1) 2^h b_(i,j,h),
    Y_(i,j)=sum_(h=0)^(ell-1) 2^h c_(i,j,h).

Impose b(b-1)=c(c-1)=0, and

    X_(i,j)+Y_(i,j)=B.                               (3)

For each time i<L introduce selectors e_(i,I),e_(i,O),e_(i,Z),e_(i,J), with
Boolean equations, pairwise products zero, and sum equal to one. Impose

    X_0=t,
    X_(i+1)=sum_(a in {I,O,Z,J}) e_(i,a) A_a X_i,      (4)
    K(X_L,left,X_L,right)=n.                          (5)

Here A_a is the exact integer matrix of the corresponding displayed move,
and A_I is the identity. Every equation has total degree at most two.
One can substitute the linear expressions X,Y directly, so the only
variables needed are Boolean bits and selectors. There are
16 ell(L+1)+8L variables and O(ell(L+1)) equations. Coefficients have
O(log n) bit length; the largest range coefficient is at most B. The total
binary description length is polynomial in log n.

Call this system C_n. Its size bounds count both traces, their range slacks,
and padding, rather than treating boundedness as an informal inequality.

## 4. An explicit one-prime local-to-exact theorem

**Theorem.** For each integer n>=3 and each characteristic-zero field K,
the following are equivalent:

1. n is the Hamilton-path count of a strong tournament with exactly one
   Hamilton cycle.
2. C_n has a solution over K.
3. C_n has a solution over F_p for one prime p>4n.
4. C_n has a solution over F_p for every prime p>4n.

The field K may be a number field with nonreal embeddings. No definability
hypothesis for its ring of integers is used.

**Proof.** If n is represented, the physical midpoint theorem supplies two
sources of covers<=B. Their actual words have length<=L by (2), and identity
padding gives two traces satisfying all equations with ordinary Boolean
bits, including the binary expansion of B-X. This integer solution reduces
to every field and establishes (1)->(2),(4).

Over a characteristic-zero field, a Boolean root is literally 0 or 1.
The one-hot equations therefore select one letter at each step. Induction
from the integer seed forces all X values to be the corresponding actual
integer source coordinates. Equation (3) enforces their height, and (5)
is exact integer equality. Hence (2)->(1).

Now take p>4n. Every Boolean variable in F_p is again 0 or 1. Let X,Y be
the ordinary integers encoded by their bits. Set R=2^ell-1; then R<=2B-1
and X,Y<=R. Since B<=n and

    0<=X+Y<=2R<4B<=4n<p,

(3) forces the ordinary equality X+Y=B, hence 0<=X<=B. The selector
equations identify a unique actual letter. On inputs and outputs in [0,B],
each coordinate residual of a selected equation (4) lies between -2B and
2B. It is zero modulo p and has absolute value<p, so it is zero over Z.
Induction therefore lifts the whole trace to actual integer sources, not
just endpoints lying in a positive cone.

The physical upper bound gives 0<K<=U_left U_right<=B²<=4(n-1)<p. Also
0<n<p. Thus (5) forces exact equality K=n, proving (3)->(1). Finally (4)
implies (3) since primes exist above any integer. QED.

This is a restricted Hasse-style theorem for a deliberately bounded finite
physical language. It does not prove local solvability for an arbitrary
target. In particular translating coverage to solvability of C_n does not
settle coverage. It also does not establish a uniform Diophantine definition
of any ring of integers: n determines L, ell, and the number of equations.

The same construction works over every finite field of characteristic
p>4n, because all encoded bits and all forced state values lie in its prime
subfield. It does not apply merely from the field having cardinality>4n
when its characteristic is small.

## 5. The terminal modulus cannot ignore scalar height

Even after all bounded transitions lift, the final count can wrap.
For n=67, B=16, the actual source words

    OZO -> (5,0,14,1), OOO -> (1,0,16,1)

have K=201. Their traces meet the prescribed cap and length bound. Therefore
C_67 has a solution over F_67, since 201=67 modulo67. Yet 67 is absent from
the complete C1 pair bank at cap16. The prime 67 exceeds 4B=64, so all
range and transition equations lift exactly; the only failed implication
is congruence of the terminal count -> equality. This is a C1 gap only,
not an unrestricted Hamilton-spectrum gap. The known C3 realization is not
contradicted.

## 6. A bounded CRT hostile and the exact synchronization repair

For n=21, B=8, the complete bounded C1 pair image is

    {3,5,9,11,15,17,19,23,25,27,29,31,33,35,37,39,45,47,51}.

It omits 21. Nevertheless these literal bounded source pairs give local
solutions of the very same system C_21:

    prime 3:  words ('','')     give K=3;
    prime 5:  words ('OZ','')   give K=11;
    prime 7:  words ('OZ','ZJ') give K=35;
    prime 13: words ('OZ','OO') give K=47.

Empty words are padded with I. The modulus product is1365, larger than
B²=64 and n. Separate existence at these primes still does not imply an
exact source with K=21.

There is a literal solution over Z/1365Z obtained by applying CRT to every
variable of these four field solutions. It is not an ordinary Boolean
source word. At the first left-trace position its selectors include

    e_I=910, e_O=456, e_Z=e_J=0 (mod1365).

Both 910 and456 are nontrivial idempotents, their product is zero modulo1365,
and their sum is one modulo1365. Thus even *all* the Boolean and one-hot
polynomial equations hold. An integral domain forces a Boolean element to
be0 or1; a product of fields permits different branches on different
factors. CRT has glued different physical histories, not one physical
history. The range bits exhibit the same issue.

**Synchronization lemma.** Let m_1,...,m_r be pairwise coprime positive
integers and M their product. Suppose one fixed pair of *ordinary Boolean
selector strings* determines actual sources with U<=B, and their scalar
K is n modulo every m_i. If M>max(B²,n), then K=n.

Indeed the common word fixes one ordinary K in [1,B²]; CRT gives M|(K-n),
and |K-n|<M. No separate local state choices or field-dependent words are
allowed in this premise. The lemma also works with one shared literal
integer pair, regardless of how its certificates were obtained.

The distinction is the quantifier order:

    exists one word, for every chosen modulus: valid

is useful, whereas

    for every chosen modulus, exists a possibly different word: valid

can fail even with a fixed ordinary-height bound and a product modulus far
larger than every possible count. Uniform boundedness plus *every positive
modulus* is the older compactness equivalence; this finite-prime shortcut
needs the stronger witness synchronization.

## 7. Consequence for the number-field question

The useful transport is not to assume that defining O_K in K automatically
creates physical integer sources. Instead the finite witness compiler
forces its own discrete choices over K, so field extension creates no new
C_n solutions. Arithmetic specialization is then safe only after explicit
height control, and assembling separate specializations additionally needs
one common word. These are concrete mechanisms that can be imported from
Diophantine reasoning without claiming the external definability question
or the tournament coverage conjecture has been resolved.

A possible finishing route now has an exact contract: prove solvability of
C_n over one prime p>4n for every relevant n, or construct one bounded pair
whose same words satisfy a sufficiently large CRT modulus. Merely proving
solvability modulo each small prime independently is insufficient, as the
21 hostile demonstrates. A specialized finite-field construction could
settle coverage of a chosen class, but no such universal construction is
asserted here.
