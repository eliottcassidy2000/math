# Shared primitive repairs preserving the full six-coordinate clock

**PROVED for the stated primitive C3/C5 language, every precision e>=1.**
This strengthens the pair-only legalizer
[PC2](../rooted_carry_closure_20261103/PHASE_49_50_51.md).
It does not assert preservation of all eighteen endpoint polynomials.

## PL1. Full matrix clocks lift at every precision

Let A3,A5 be the six-coordinate matrices of [JC1](JOINT_CARRIER.md).
Their exact orders modulo 7^e are

    a_e=12*7^(e-1),       b_e=48*7^(e-1).

The exact base orders and the nonzero matrices (A_j^order-I)/7 modulo7
are saved in [algebra.json](algebra.json). If B=I+7^t D with some entry
of D a unit, the binomial theorem gives B^7=I+7^(t+1)D modulo7^(t+2).
Induction gives the claimed exact lift, while reduction modulo7 forces
every possible order to be a multiple of the base order. This proves all e.

## PL2. A synchronized three-source legalizer

A primitive controller source has odd positive pair b>t>0, gcd(b,t)=1,
3 divides t, 3 does not divide b, and projective r=t/b modulo11 in

    S={0,1,3,4,5,6,7,infinity}.

The inherited sourced macros satisfy

    C3(b,t)=(61b-2t,9b),   C5(b,t)=(421b-2t,33b).

C3 is always primitive-legal, acting in cycles (4,1,5,6) and
(0,7,3,infinity). C5 is legal exactly when r!=7 and then sends r to0.
These are arithmetic primitivity guards; the actual graph programs remain
defined even outside them.

Given at most three supplied sources, choose j in {0,1,2,3} such that no
C3^j(r_i) is7. Each source forbids at most one choice. Define the word

    L_j = C3^j C5^b_e C3^(a_e-j).

It acts as the identity on all six coordinates modulo7^e. Its first C5
is legal by construction; all later C5 steps start at0 and are legal.
Its common final controller state is C3^(-j)(0), namely
0,infinity,3,7 for j=0,1,2,3. When j<=2, append the desired C5. When j=3,
first use L_1, which is legal from7 and ends at infinity, then append C5.
The final six-coordinate action is exactly A5 modulo7^e for every source.

Each neutral L has 4a_e+8b_e=432*7^(e-1) native steps, since C3 and C5
cost four and eight. Therefore the total bounds, including the desired C5,
are

    <=432*7^(e-1)+8 for one or two supplied sources;
    <=864*7^(e-1)+8 for three supplied sources.

For two sources use j in {0,1,2}; at most two choices are forbidden.
Each L adds 6a_e+8b_e=456*7^(e-1) permanent vertices and half its native
step count in temporary helpers. The root and second live port follow the
inherited replacement program: C3/C5 replace the old second live port at
their first birth. Its vertex remains in the graph, but it is no longer
a live attachment port. NR alone keeps the second port. These are constructive word
bounds, not minimum lengths. No exponentially large literal graph was
enumerated in the finite tests.

## PL3. Four supplied phases can obstruct every common primitive word

Take r=(0,7,3,infinity), realizable by four successive C3 variants of an
actual entry source. Before the first C5 in any common primitive word,
only C3 has occurred. Its action permutes this four-set, so one supplied
source is always at forbidden7. Thus no common legal word contains C5.
Nor can a power of C3 replace its required phase action: det(M3)=4 modulo7
has powers in {1,2,4}, whereas det(M5)=3. For the actual H497 entry, the
first two supplied pair residues are (1,3),(6,2), with nonzero determinant5
modulo7. Agreement on those sources would already force equality of the
two pair matrices. This is an obstruction in the
entire primitive two-letter language, not just the displayed repair scheme.
NR, other graph edits or relaxation of primitivity lie outside this claim.

[legalizer.py](legalizer.py) checks all 8+64+512=584 ordered source tuples
at each precision e=1,2,3, the full matrix identities, every intermediate
mod11 guard and the exact integer trajectory from the actual H497 k=2
entry (2828575,222057). The all-order proofs supply the unbounded extension.
