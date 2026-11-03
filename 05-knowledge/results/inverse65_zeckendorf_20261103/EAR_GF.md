# Literal ear lifts and a state-valued generating function

Status: **FINITE-EXACT** classification of30 specified ears; **PROVED** formal
series identities. No interval theorem or isomorphism-class enumeration follows.

## E1. The two supplied integer lists have different types

Let T be the five-vertex tournament with lexicographic pair code8, bit1 meaning
i->j for i<j. It has H=9. Add x with x->i exactly for set bits in mask S.
For every nonempty proper S (30 masks), literal Hamilton-path enumeration gives

    2z^15+4z^17+2z^19+6z^23+2z^25+2z^27+2z^29
       +4z^33+4z^37+2z^41.                              (E1)

This recovers, rather than newly asserts, the support
`A={15,17,19,23,25,27,29,33,37,41}` from
[THM-4118, ear-response-lattice-and-stateful-unit-component-intervals](../../../01-canon/theorems/THM-4118-ear-response-lattice-and-stateful-unit-component-intervals.md)
and [Session152 ear work](../session152_sequences_20261103/ear/README.md).
The other supplied list
`B={2,3,5,7,11,13,17,19,23,27}` is a filter, not a prime list (27 is composite).
Keeping terms in B gives `4z^17+2z^19+6z^23+2z^27` and14 actual masks.

There is no canonical numerical bijection from the two lists. In particular,
replacing a PRIMEGAME prime coordinate29 with27 is not a harmless coordinate
renaming. The [prior history audit](../session152_sequences_20261103/history/HISTORY.md)
also distinguishes the ear support from square-sum feasibility and Ogg-prime
intersections. Those appearances do not supply a source-preserving operation.

## E2. A literal lift into65 symbols, and its exact failure

For every one of the30 children, enumerate all directed Hamilton paths p.
Relabel along p, apply the spine complement K, parse its ordered strong
components, and check whether each marked component is in the fixed bank.
When accepted, reconstruct the full target rows. They must equal the original
child relabelled by reversed p, not merely have the same H. The processing
port corresponds to the original vertex p[0].

Exactly28 children lift, with143 accepted charts in total. The two failures
are masks24 and28, both with H=41. Every one of their41 path charts has one
strong six-vertex source atom outside the bank. Their sets of repair atoms
are identical. Consequently adding any one of these41 atoms repairs both;
the least numerical marked mask is s76. Zero additions cannot work, so one
is minimal for this particular literal ear family. This is a proposed66-symbol
extension, not a change to the65-symbol computations elsewhere in this packet.

The existing65 bank attains41 at order7. That does not lift either actual
order6 ear. Source and vertex clock separate the two statements.

Every one of the14 masks selected by B lifts. A canonical path choice gives
different three-coordinate states even at equal H. The six H=23 states give
singleton responses `49,47,81,39,41,47`. The two H=19 states give31 and33.
This is a small exact hostile to using (E1) alone as a recursive carrier.

## E3. Combine recursion without discarding the state

For a row state v=(h,x,y), write `M_v=X^h Y^x Z^y`. For letter l with matrix
P_l and vertex cost c_l, define the monomial substitution S_l by replacing
variable i with `prod_j variable_j^(P_l[i,j])`. Then

    S_l(M_v)=M_(v P_l).

The formal series enumerating **presentations**, with multiplicity, is

    F(u;X,Y,Z)=X + sum_l u^c_l S_l F(u;X,Y,Z).          (E2)

Every nonempty word has a unique last letter, proving (E2). Coefficients at
fixed u-degree are finite; even after setting u=1, coefficients at fixed h
remain finite by the norm bound in INVERSE.md. Positivity of a state
coefficient is equivalent to actual reachability. Specializing Y=Z=1 loses
the data needed to iterate the substitutions, so the scalar projection is
not itself a closed recurrence.

For the14 filtered ears, start instead with the finite polynomial

    G_0(z;X,Y,Z)=sum_actual_filtered_masks z^(initial H) M_(chosen state).

Its coefficients, charts, ports, and masks are in `ear_bridge.json`. Replace
X by G_0 as the seed in (E2); the resulting series retains the initial ear
and correctly generates all its same-end continuations. If initial labels
are needed, replace z^(initial H) by a separate variable for each mask.

For a cheaper aggregate, let

    R(z)=z^17(68,32,36)+z^19(38,13,26)
         +z^23(138,83,105)+z^27(54,24,30).

Then `R(z) (I-sum_l u^c_l P_l)^(-1) e_h` sums H over all continuations by
added vertex cost. This rational series is useful for moments, but its
coefficients are sums of H values, not their support or multiplicities by H.
It cannot certify that any one particular target occurs.

## E4. Which algebra matches which operation?

The [older additive calculation](../session152_sequences_20261103/number_theory/NOTE.md)
has Apéry set modulo15
`(0,46,17,33,19,50,36,37,23,54,25,41,27,58,29)` for the semigroup generated
by A. Its conductor is44 and genus26; its value-set series is
`sum_w t^w/(1-t^15)`. B generates the additive semigroup of2 and3.
These are exact arithmetic statements about sums. Ordinal sums of tournaments
multiply H, while the seam uses a matrix action, and arbitrary ears use a
directed-cut response. Substituting one generating-function product for another
would change the operation being modelled.

Session152's [integral cut carrier](../session152_sequences_20261103/cut_induction/README.md)
and [fixed-rank moments](../session152_sequences_20261103/cut_moments/README.md)
offer a complementary exact coordinate system for arbitrary ears. Its source
is a particular parent and its image is a cut response; the cut matrix alone
has not been shown to be a closed next-parent state. The literal lift above
is the concrete bridge: it tests which actual ear children enter the closed
same-end matrix family, rather than identifying the two carriers by analogy.
