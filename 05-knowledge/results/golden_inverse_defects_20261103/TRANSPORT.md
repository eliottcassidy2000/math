# Two realized ternary offsets transported through an infinite prefix family

**PROVED construction family.** These are actual source-word replacements,
not weighted averages or formal vectors declared realizable. The family
does not close under iterated ternary expansion and does not prove either
complete prime lane.

## T1. Three literal seeds

In the inherited65-presentation alphabet take

    U  = a a r5,               v(U) =(31,14,25),
    U2 = a a q1 a a,           v(U2)=(95,41,75),
    U4 = a a a a a q0,         v(U4)=(97,40,75).        (T1)

Their source orders are7,8,9. Every word has an explicitly stored source
relation, a legal concatenated Hamilton spine and its literal path
complement. Independent subset DP checks their counts. Set

    tau=(2,-1,0).

Then v(U2)=3v(U)+tau and v(U4)=3v(U)+2tau. The negative coordinate of tau
is a difference of actual positive states; it is never treated as a
physical source state.

The inherited matrices give the exact transport identities

    tau P_q2=tau,              tau P_b=-tau.           (T2)

For any word w in {q2,b}, let epsilon(w)=(-1)^(number of b symbols).
Multiplication proves

    h(w U2)=3h(w U)+2 epsilon(w),
    h(w U4)=3h(w U)+4 epsilon(w).                      (T3)

Thus **every prefix with an even number of triangle blocks has both
positive offsets**. The number of triangles can be arbitrarily large and
q2 blocks can occur in any order among them. Odd parity gives the two
negative offsets. The bit is necessary: for w empty the first offset is
+2, whereas w=b reverses it to-2. If the full word is stored, the bit is
computed from that word, not extra independent information.

All these targets are strong. Each source has multiple ordered components
and at least three vertices; its endpoint arc is a nonpath arc, so actual
complement closes the retained Hamilton spine into a cycle.

The operation retains the entire prefix and replaces a specified suffix.
It is a bounded symbolic source replacement, adding one or two vertices.
It does **not** retain every old edge, nor is it asserted to be a single
edge of the complement/re-marking metagraph. In particular the old r5
strong five-vertex source block cannot embed as an induced strong block
inside q1 a a, whose largest strong component has size4. This explains
how the construction escapes the insertion-only obstruction.

## T2. All five Q4 phases remain visible

The transported corrections are

| Source symbol | tau P_c |
|---|---|
|q0|(-6,-1,-2)|
|q1|(-6,0,-3)|
|q2|(2,-1,0)|
|q4|(0,4,3)|
|q5|(0,3,2)|

Only q2 among these five preserves the correction. q1 and q4 have equal
pure-block characteristic polynomials but different correction transport.
The full inverse continues to retain all five; the smaller prefix grammar
is a proved subfamily, not a quotient that identifies their phases.

## T3. The two quadratic clocks and an integral parity guard

Introduce integer coordinates

    Y=y,       X2=2x+h-y,       g=h-y.                 (T4)

Their image lattice has index2: X2=g modulo2. Reconstruction is
h=Y+g, x=(X2-g)/2, y=Y. Forgetting that congruence would invent fractional
states. On the actual seeds the bulk coordinates are

    (Y,X2,g)(U) =(25,34,6),
    (Y,X2,g)(U2)=(75,102,20),
    (Y,X2,g)(U4)=(75,102,22).

For row vectors the two actions split exactly:

    q2: [[11,18,0],[4,7,0],[0,0,1]],
    b : [[5,8,0],[2,3,0],[0,0,-1]].                  (T5)

The bulk is multiplied by3 between corresponding seeds, while g carries
the constant correction and its sign. The bulk characteristic polynomials
are X^2-18X+5 and X^2-8X-1. Their roots are

    9 +/- 2 sqrt(19),       4 +/- sqrt(17),

with norms5 and-1. This is an actual occurrence of17 and19 in the source
transfer, not a proposed identification of two prime lists. The full
three-dimensional characteristic factors include respectively X-1 and
X+1. The eigenvalue1 is precisely what permits the constant correction;
the singleton tribonacci clock has no such eigenvalue and gave the earlier
fixed-wrapper obstruction.

The clocks retain their arithmetic failures. The triangle block is
invertible over every Z/m; q2 is invertible exactly when gcd(m,5)=1.
Modulo17 the triangle bulk is4I+N with N nonzero and N^2=0, and4 has
order4. Its full order is68. Modulo19 the q2 bulk is9I+N, with N nonzero
and N^2=0;9 has order9, so the full order is171. Dropping N would lose
the prime17 or19 clock factor. The receipt independently checks these
orders on the original three-dimensional integer matrices. These are
matrix iteration clocks, not Frobenius automorphisms.

For the pure q2 prefix define h_k=h(q2^k U). It begins

    31,417,7279,128865,2283103,40451457,...

and satisfies the all-k recurrence

    h_(k+2)=18h_(k+1)-5h_k-72,
    sum h_k z^k=(31-172z+69z^2)/(1-19z+23z^2-5z^3).   (T6)

Indeed g=6 stays fixed, and the bulk satisfies its quadratic polynomial.
Equivalently h_k-6 is homogeneous with that quadratic recurrence. The
companions3h_k+2 and3h_k+4 therefore give two explicit infinite strong
families. No primality or7p coverage is inferred from their formulas.

## T4. Generating functions for the whole mixed prefix language

Let E_n,O_n count even/odd triangle-parity prefixes of vertex cost n,
where q2 costs4 and b costs3. Then

    E_n=E_(n-4)+O_(n-3),     O_n=O_(n-4)+E_(n-3),
    E_0=1, O_0=0; negative indices are zero,
    E(z)=(1-z^4)/((1-z^4)^2-z^6).                     (T7)

Retain row vectors Z^even_n,Z^odd_n for sums of (Y,X2). With
U_q=[[11,18],[4,7]], U_b=[[5,8],[2,3]], they satisfy

    Z^even_n=Z^even_(n-4) U_q+Z^odd_(n-3) U_b,
    Z^odd_n =Z^odd_(n-4) U_q+Z^even_(n-3) U_b,          (T8)

initialized by Z^even_0=(25,34), Z^odd_0=0. This is an explicit rational
matrix generating function: invert I minus the4-by-4 block matrix with
z^4 U_q on its diagonal and z^3 U_b off its diagonal. The original
Hamilton masses are F^even_n=(Z^even_n)_1+6E_n and
F^odd_n=(Z^odd_n)_1-6O_n. The two output masses are respectively

    3F^even_n+2E_n, 3F^even_n+4E_n;
    3F^odd_n-2O_n,  3F^odd_n-4O_n.                    (T9)

Here n is prefix cost; actual tournament orders are n+7, n+8 and n+9.
These formulas count marked words with multiplicity, not distinct H values
or isomorphism classes. Forgetting parity would average incompatible
signs and would no longer prove the individual constructions(T3).

## T5. The precise closure failure

The first output states do not permit a naive repeated carry update.
Exact inversion proves that the four triples

    (287,122,225), (289,121,225),
    (293,119,225), (295,118,225)

are unrealizable in the65 bank at every order. These are3v(U2)+tau,
3v(U2)+2tau,3v(U4)+tau,3v(U4)+2tau. Their H values themselves are realized
elsewhere in the bank. Thus the construction is an infinite transported
family of valid offsets, not an invariant ternary state bank.

A finite seed search with H<=10000 found98 states supporting at least one
of the required affine line lifts within that bound and four supporting
both. Three of the four arise by further transport; no universal seed
classification follows. The exact line obstruction gives a sharper next
target: change the output suffix or its full state while keeping the desired
H, and supply a closed family of such changes with a proved realization
guard. Arithmetic coverage of H alone cannot authorize that reset.

## T6. Restoring every Q4 phase enlarges the correction carrier

For any prefix over the entire65 bank, the literal seeds still give

    v(w U2)-3v(w U)=tau P_w,
    v(w U4)-3v(w U)=2tau P_w,                          (T10)

where P_w is the matrix product in the inherited reversed source-word
order. The first coordinate of this correction need
not remain2 or4. The restricted sign character in(T2) is only one part
of this full transport.

Let L be the smallest integer module containing tau and closed under all
65 matrices. Then

    L = {(2u,x,y): u,x,y in Z}.                        (T11)

Every bank matrix has even entries P[1,0],P[2,0], so the right-hand side
is invariant. Conversely the three rows tau, tau P_q4, tau P_q5 have
determinant-2 and hence generate that entire index-two lattice. Thus all
five Q4 phases force rational correction dimension3; the one sign bit
does not extend to them. For example q4 and q5 give zero current count
correction, but a subsequent singleton observes8 and6, respectively.

There is also an exact modular observation quotient. Declare two integer
corrections equivalent when every future left-bank context has the same
Hamilton correction modulo m. On L the minimum quotient is

    (Z/(m/gcd(m,2)))^3.                               (T12)

Indeed the empty, a and aa observations on(2u,x,y) are
2u,2(u+x),2(u+x+y). After dividing by2 with the correct modulus they form
a unimodular triple, so they separate exactly this quotient. All later
observations are twice integer linear forms and cannot separate more.
This is a minimum for arbitrary integer correction combinations. It is
not a proof that every combination is a realizable difference of individual
source words, nor a minimum carrier for arbitrary complement/re-marking.
The inverse from THEORY remains the realization guard.

## T7. Prime-lane test and a proved sparsity boundary

Saturating every prefix state whose positive output is at most one million
gives20 distinct prefix states,9 of even triangle parity and18 distinct
positive outputs. Their prime outputs are

    97,4729,21839,21841,81197,81199,83597,

and their sole7p output is1253=7*179. The receipt retains a literal source
word for every output. This tests an infinite construction by a value-bound
saturation; it is not a new global spectrum coverage record.

Every b or q2 prefix step multiplies h by at least5, since x+y>=h.
Consequently a prefix of length l has h>=31*5^l. With outputs bounded by X,
only O(log X) lengths occur, and at most

    2 sum_(l<=floor(log_5(X/93))) 2^l = O(X^(log_5 2))

outputs are possible. Since log_5(2)<1, this family has zero natural
density among odd integers and cannot alone complete their spectrum.
This boundary motivates a larger realizable seed bank or a different
inductive operation; it does not invalidate any of the individual maps.
The incoming [weighted-automaton obstruction](../session153_gadgets_20261103/rewriting/NOTE.md)
uses this same multiplicative-scale counting mechanism in a more general
finite affine rule system. Our inequality is the direct supplied-matrix
version and is not claimed as a new general density principle.

## T8. A literal F49 map and a complete7-divisibility clock

For the pure q2 family reduce the bulk in(T5) modulo7. Define

    F49=F7[theta]/(theta^2-4theta+5),
    (Y,X2) -> alpha=Y theta+4X2.                      (T13)

The discriminant3 is a nonsquare in F7. Multiplying alpha by theta sends
(Y,X2) to(Y,X2)[[4,4],[4,0]], precisely the actual q2 bulk modulo7.
The map is an invertible F7-linear map. It preserves the supplied q2
clock and its linear Hamilton observation; it forgets all multiples of7,
so it does not preserve primality of a lifted quotient.

Frobenius sends theta to4-theta. Thus theta^8=5, theta^24=-1 and
theta^16=4, proving theta has order48. The actual seed has

    alpha_0=3+4theta=theta^38.

Consequently alpha_k=theta^(k+38), and Frobenius on this phased orbit is

    k -> 7k+36 modulo48.                             (T14)

It is an involution, but not the unshifted k->7k map. Clock multiplication
has order48 while Frobenius has order2. The seed phase is essential.

The48 nonzero field elements are exhausted. Since H=Y+6 modulo7, its
residue6 occurs six times and every other residue seven times. The two
zero conditions for the ternary outputs are nonzero linear fibres,
explaining why each has exactly seven positions. Explicitly,

    7 | 3h_k+2 iff k mod48 in {1,4,5,19,30,40,47},
    7 | 3h_k+4 iff k mod48 in {6,16,23,25,28,29,43}.   (T15)

This is an all-k classification, including the inherited recurrence's
integer lifts, not a conjecture extrapolated from the first cycle. Each
class gives an infinite family of strong tournaments with H divisible
by7. Membership in the7p lane additionally requires the quotient to be
prime; no infinitude of such quotients is proved.
