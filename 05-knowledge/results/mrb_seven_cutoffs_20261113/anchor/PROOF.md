# Full-state loops at 7, 14 and 77, and the missing synchronization bit

**Status: PROVED algebraic classifiers and integral-word obstruction;
FINITE-EXACT complete reachable finite quotients, with exact positive-word
certificates.** The source is a rooted directed triangle. The operation
alphabet is exactly the inherited physical L/K birth language. This is not
a classification of arbitrary tournaments or a Chan-Pak proof.

The main new result is that the modulo77 quotient is **smaller than the
product of the modulo7 and modulo11 quotients by exactly12000 states**.
The lost coordinate is the parity of the final run of L births after the
last K. A complete-state arithmetic loop therefore needs synchronized
phase even before addressing the unbounded integral Hamilton count.

## A. Inheritance and the shifted forbidden pair

The inherited [binary language](../../cycle_opening_six_return_20261112/returns/BINARY_LANGUAGE.md)
and [neutral language](../../cycle_opening_six_return_20261112/returns/NEUTRAL_LANGUAGE.md)
give the physical operations and full-count observer. The
[strict-ancestry boundary](../../cycle_opening_six_return_20261112/compiler/ANCESTRY_BOUNDARY.md)
explains why these descendants have two indegree-one vertices, whereas
the successful five-core repair leaves that restricted ancestry type.
The broader low-defect classification in
[mutation/Collatz work](../../mutation_collatz_chanpak_20261112/README.md)
already shows that unit defect on the two-root stratum forces a book count.
Section F below gives the exact word preimage for this particular seed and
alphabet; it is a refinement of that inherited boundary.

The closest proved mechanism is the three-coordinate compiler. The hostile
is a positive defect congruent to1 but larger than1. The corrected near
miss is the already ruled-out scalar-to-physical Collatz map in
[central modular gauge](../../central_modular_chanpak_20261110/gauge/PROOF.md).
The least-used sidecar here is the **shared suffix parity**, beyond separate
prime-local states.

Write

    q=(H-1)/2,        r=q+1=(H+1)/2.

The path-count omissions H=7,21 become q=3,10 or r=4,11. These are two
coordinate descriptions of the same two counts, not additional missing
values. The recovered [microcosm gate note, section3](../../microcosm_gate_doubling_resolvent_20261106.md)
already supplies q and its square-sum comparison. The r coordinate makes
the desired unit-defect target particularly simple:

    H' = 3H+2    if and only if    r'=3r.

Modulo7, multiplication by3 has one fixed zero and one six-cycle. Modulo11
it has one fixed zero and two five-cycles. Modulo77 the complete cycle-length
inventory is

    1 fixed point, 2 cycles of length5, 1 of length6, 2 of length30.

Indeed the four gcd strata have sizes1,10,6,60, and the orders of3 are6
modulo7 and5 modulo11. The two exceptional integers r=4,11 coincide modulo7
but lie in different modulo77 strata: their scalar periods are30 and6.
The difference14 between H=7 and21 is the lift of a difference7 in r.

These target clocks are distinct from the actual L/K clocks below. Even if
every desired ternary step could be made ready again, six such steps change
the integer H to729H+728. Returning to a residue does not preserve H.

## B. Physical state and the last-K decomposition

For a marked indegree-one root, retain

    v=(x,b,d),  h=x+d,  d=H(T-root),  b=d-start_root(T).

Here b is the defect, not a residue or a formal free parameter. The triangle
has v0=(2,0,1). Chronological L adds a vertex beaten only by the old root
and marks the newborn; K performs the same birth and keeps the old root
marked. Every word is an actual strong rooted tournament. The inherited
exact maps are

    L(x,b,d)=(x+2b,d-b,x+d),
    K(x,b,d)=(2x+2b,0,d).

All counts in this note are total Hamiltonian-path counts. These three
coordinates are complete for future L/K count queries, not for arbitrary
selective attachments. In particular the empty, L and LL H observations
recover b,d,x by the formulas in the inherited neutral-language proof.
Here recovery is over the integers (and over odd residue moduli). The
modulo14 state below is a sufficient retained coordinate quotient; no
minimal-observer claim is made there, where division by2 is unavailable.

Let P_m be the plane b=0 modulo m. Every word can be written as a word ending
in K followed by L^j, or as L^j from the seed, which already lies in P_m.
Consequently its residue lies in

    union_(j>=0) L^j(P_m).

For m=7,11,77, the exact finite reachability certificate proves that every
point of P_m is reachable. Thus this union is exactly the reachable set.
For14 the accessible plane is b=0 with x even and d odd, and again every
point of that plane is reachable. This is a complete closure computation
over m^3 states, not a search cut off at a chosen word length.

The characteristic polynomial is t^3-t^2-t-1 and det L=1. The exact matrix
orders are48 modulo7,110 modulo11,48 modulo14,2640 modulo77. Define a phase
mask epsilon_p(v) by including bit0 if v belongs to L^(even)(P_p), and bit1
if it belongs to L^(odd)(P_p). The same tests can be made by asking when a
forward L iterate of v has b=0: the finite matrix orders are even.

## C. An explicit norm classifier modulo7

For v=(x,b,d) over F7 put

    a=6x+3b+5d,
    u=b-a,               w=5(d-b)-a,
    N=6u^2+3uw+w^2.

All these expressions are in F7. Then

    even phase is possible iff a^2-4N is in {0,3,5,6};
    odd  phase is possible iff a^2+4N is in {0,3,5,6}.        (1)

The set on the right consists of zero and the three nonsquares. Both phases
can occur. Neither phase occurs precisely on the six nonzero points of the
eigenline span(1,1,4).

**Proof.** Over F7 the characteristic polynomial factors as
(t-3)(t^2+2t+5). The expression a is the b coordinate of the eigenvalue3
projection. The complementary two-dimensional space identifies with F49.
After dividing L by3, its multiplier t on that space satisfies

    t^2-4t-1=0,   Tr(t)=4,   Norm(t)=-1,   t^8=-1.

Hence t has order16; its even powers are exactly the norm-one circle and
its odd powers exactly the norm-minus-one circle. There is a unique c in
F49 with Tr(c)=u and Tr(ct)=w. The trace-pair map is invertible, and direct
elimination gives Norm(c)=N above. The actual second coordinate obeys

    b(L^j v)/3^j = a + Tr(c t^j).

If c is nonzero, asking for b=0 at an even or odd time asks for an element
of trace -a and norm N or -N. In a quadratic extension this is equivalent
to its quadratic discriminant being a nonsquare or zero. A nonzero square
would give two distinct elements of F7, neither having the required
extension trace and norm. This proves(1). If c=0, the same test excludes
a nonzero a and includes the zero vector for both phases.

The phase counts, obtained either from(1) or the independent explicit L
orbits, are

| phase mask | neither | even only | odd only | both |
| --- | ---: | ---: | ---: | ---: |
| number of vectors | 6 | 120 | 120 | 97 |

For example each nonzero norm has eight preimages. For a=0 exactly24
nonzero vectors permit each parity, with no overlap. For each nonzero a,
four nonzero norm values permit each parity and two permit both. Including
zero gives217 states per parity,97 in their intersection,337 in their union.

## D. A Jordan-boundary classifier modulo11

For v=(x,b,d) over F11 set ell=2x+6b+d. If ell is nonzero, **both** phases
occur. On ell=0 write

    a=4x-b,       c=b-3x,
    v=a(1,3,2)+c(1,4,7).

Zero permits both phases. If exactly one of a,c is zero, neither occurs.
If both are nonzero, put z=6c/a. Then even phase occurs exactly when z is a
quadratic residue modulo11, and odd phase exactly when it is a nonresidue.

**Proof.** The characteristic polynomial is (t-7)^2(t-9). The displayed
vectors are eigenvectors for7 and9, and ell=0 is their span. In that plane,
b(L^j v)=3a7^j+4c9^j. Since7/9=2 modulo11, the zero condition is
2^j=6c/a. The element2 has order10 and is a nonsquare, so its exponent parity
is exactly the claimed Legendre sign.

Outside that plane write

    b(L^j v)=7^j(A+Bj)+9^j C,    B=9ell !=0.

For each chosen j modulo10, the equation b=0 fixes one j modulo11. The
ordinary Chinese remainder theorem makes these compatible modulo110.
Choosing an even or an odd class modulo10 gives both parities. This also
exhibits the transient-looking Jordan term as an exact finite-field clock,
not an assumption about asymptotic decay.

The phase counts are

| phase mask | neither | even only | odd only | both |
| --- | ---: | ---: | ---: | ---: |
| number of vectors | 20 | 50 | 50 | 1211 |

The two excluded eigenlines contribute20 nonzero vectors. The100 vectors
with a,c nonzero split50/50. All1210 vectors outside the plane and zero
permit both. Thus each parity set has1261 vectors and their union1311.

## E. Exact composition at14 and77

**Modulo14:** the reachable states are exactly

    x even, d odd, and phase7(v) contains parity b mod2.       (2)

The integer parity laws give x even and d odd. On their parity plane, L
toggles b and K makes b zero. Therefore the parity of a final L run is
exactly b modulo2. Each admissible F7 vector has a unique lift for each
chosen b parity, and the counts are217+217=434. Without that shared parity,
the product of separately accessible mod2 and mod7 states would contain674
states and overcount by240.

**Modulo77:** the reachable states are exactly

    phase7(v) and phase11(v) have a common parity.            (3)

Necessity follows from the same final run L^j at both primes. For sufficiency
choose local successful exponents j7 modulo48 and j11 modulo110 with that
common parity. Their congruences are compatible because gcd(48,110)=2.
They therefore have a common j modulo2640. The point L^(-j)v lies in P77,
which the finite certificate reaches; appending L^j realizes v.

The resulting exact count is

    217*1261 + 217*1261 - 97*1211 = 429807.

The naive product337*1311=441807 overcounts by

    120*50 + 120*50 = 12000.

These are precisely the incompatible even-only/odd-only pairs. This is a
failure of independently choosing local **words**, not a failure of CRT
on ordinary residue vectors. CRT assembles every vector; the common physical
word imposes the extra condition.

A small explicit hostile is v=(5,1,72) modulo77. Its modulo7 projection is
reached by KLLLKLLLL and permits only even final phase. Its modulo11
projection is reached by KKKL and permits only odd final phase. No L/K word
reaches their joint modulo77 state, at any length.

The exact closure receipt also supplies a positive word for every accessible
b=0 plane point, maximum shortest-word lengths, and farthest-plane examples.
The longest shortest word in the whole modulo77 quotient has length480.
Every certified word lifts to an actual rooted tournament by the inherited
physical birth rule. No assertion is made that arbitrary cone representatives
of the same residue are realizable.

## F. The exact integral ready states and a minimal false-ready witness

**Theorem.** From the triangle seed, an L/K word has defect exactly1 if and
only if it is K^c L for c>=0. Its state and full count are

    (2^(c+1),1,2^(c+1)+1),        H=2^(c+2)+1.              (4)

The empty word has H3 and defect0. A pure K word, including the empty word,
has (2^(c+1),0,1); it is not an exception to the theorem.

**Proof.** Throughout this language x>=2, d>=1, and s=d-b>=1. Under L,
d becomes x+d>d and s becomes x+b>=2. Under K, d stays fixed and s becomes d.
Thus s=1 occurs precisely on pure K words: once L has occurred, d>=3 and
neither later operation can restore s=1. K always produces defect0, while L
produces defect equal to the preceding s. Therefore defect1 occurs precisely
after a pure K prefix followed by L. Substitution gives(4).

This exact word statement lies inside the inherited strict two-root/book
classification. It does not exclude a return that changes old arcs or leaves
the two-root stratum, such as the successful five-core replacement.

The distinction from modular readiness already occurs at six births:

    KLLLKL :  (x,b,d)=(38,15,53), H=91.

Modulo7 its defect is1 and its H residue is the forbidden phase0. But two
more L births give H=3*91+2*15=303, while the desired unit-offset target is
3*91+2=275. The missing integral defect costs exactly28. Direct enumeration
of every word through six births verifies that this is the first length at
which b=1 mod7 and H=0 mod7 occur together; this minimality is finite-exact.

All seven H residues occur among states with b=1 modulo7, whereas actual
defect-one counts in(4) have only H residues{2,3,5}. Modulo77 the residue
quotient likewise permits all77 H residues at b=1, but(4) permits only30.
An integral lift of the defect is indispensable even after the complete
modular state is retained.

Every nonempty word strictly increases the total count: a birth changes H
by x+2b>=2. Thus L^48 and L^2640 give exact residue returns at the indicated
moduli but never exact H returns. For instance L^2640 preserves the full
state modulo77 for every source while increasing H by at least5280. Stronger
growth estimates are unnecessary for this obstruction.

## G. What this changes in the Chan-Pak search

The finite phase coordinates are now completely classified for one real
constructor language, including their composition at2,7,11. They can be used
as exact rejection rules and as positive modular witness generators. The
12000-state loss is invisible to separate local reachability tests.

The full target remains an integer equation. The unit-defect word theorem
shows that staying inside the L/K language cannot continue the generic ready
ternary induction. The prior five-core return escapes this exact ancestry
restriction. Incoming
[Session296 entrances](../../session296_multicore_metagraph_20261113/cores/ENTRANCES.md)
already supplies finite ready witnesses at405*2^k+71 for k=1..9 and several
successive targets beginning881; those searches use a richer grammar and
do not contradict this theorem. Its
[finite-atlas obstruction](../../session296_multicore_metagraph_20261113/cores/FINITE_ATLAS.md)
also retains the integral height obstruction to an infinite finite-bank rule.
No new uniform ternary replacement is asserted here.

**Connection contract:** actual rooted L/K tournaments map to full residue
states; the map preserves all future L/K H queries modulo the chosen modulus;
it loses integer height and defect. Separating primes additionally loses the
common suffix phase unless(3) is retained. The restoration sidecars are an
integer continuation state, a physical word, and the synchronized phase.
The hostile tests are the CRT phantom and H91 defect15 above.

The local concept board is: shifted H target; full three-state compiler;
prime-local norm/Jordan phase; synchronized composite modulus; integral ready
defect; five-core escape. The first three do not supply the last two, and the
new synchronized phase sharpens rather than removes that distinction.

## Reproduction

    python -B 05-knowledge/results/mrb_seven_cutoffs_20261113/anchor/check.py
    python -O -B 05-knowledge/results/mrb_seven_cutoffs_20261113/anchor/check.py

Both modes produce the saved receipt.json. The independent checks compare
algebraic classifiers, explicit local matrix orbits and complete positive-word
BFS closures. Literal tournament subset DP verifies all127 words through six
births, with permutation enumeration as a third path through three births.
The matrix laws themselves are inherited all-order proofs. No new term of
A000568 and no prime-infinitude result is claimed.
