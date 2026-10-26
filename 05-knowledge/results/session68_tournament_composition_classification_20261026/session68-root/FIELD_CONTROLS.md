# Independent finite controls for the arithmetic interface

Status: FINITE-EXACT, standard-library implementation independent of the arithmetic producer. General classification and normal-basis proofs belong to the arithmetic note; this file records the precise finite test, not a general empirical proof.

## F8: rotation really is Frobenius, reflection is not a field automorphism

Use F2[a]/(a^3+a+1). The code independently checks irreducibility. A normal basis is

    a+1, (a+1)^2, (a+1)^4,

represented in polynomial-bit notation by3,5,7. Its sum is1. Mapping a three-bit decoration to its basis expansion intertwines cyclic rotation with x->x^2 on every field element. Its orbit length equals its field degree. By contrast a,a^2,a^4 (bit values2,4,6) have rank2, so using them as if they were a normal basis would lose information.

Index reflection J fixes basis index0 and swaps1,2. The test verifies J Frob = Frob^(-1) J and that x->1+J(x) has no fixed point. But J is not multiplicative: in polynomial-bit coordinates,

    J(2*2)=6,   J(2)*J(2)=4.

Thus the complemented-reflection model of converse cannot be renamed a field automorphism. The basis and the native operation must remain explicit.

## F512: field degree alone loses hierarchical symmetry

The program searches monic degree-nine binary polynomials and checks irreducibility directly, then searches for a normal element by full-rank conjugates. No unverified polynomial table is imported. In the resulting normal basis it compares two binary leaf decorations in index order:

    000000001,    001001011.

On the contiguous ternary blocks {0,1,2}, {3,4,5}, {6,7,8}, enumerate all81 permutations of W2=C3 wr C3. The nine-cycle

    (0,3,6,1,4,7,2,5,8)

belongs to this group. Read the normal-basis coefficients along that cycle, not the ordinary linear index order. Both decorations have Frobenius orbit length9 and field degree9. Their stabilizers in the full hierarchical group nevertheless have orders9 and1. The arithmetic field-of-definition statistic classifies the cyclic action, not this larger wreath action. The distinction persists when these decorations are used in the tournament construction.

## CRT: relative phase is an actual lost coordinate

For fixed binary rotation orbits of periods d,e inside a common C_m action, the code enumerates diagonal pair-orbits:

- d=e=3, m=3: three relative-phase orbits, each of size3.
- d=3,e=5, m=15: one joint orbit of size15.

These agree with gcd(d,e) relative phases and lcm(d,e) orbit size. Knowing two separate necklace classes therefore does not specify their aligned pair when the gcd is greater than one. No polynomial arithmetic or multiplication-preservation claim is inferred from these two finite cases; the general proof is stated separately in the arithmetic lane.

## Replay

    python3 check_field_bridge.py > /tmp/field.json
    python3 -O check_field_bridge.py > /tmp/field-O.json
    cmp /tmp/field.json /tmp/field-O.json
    cmp /tmp/field.json field_bridge.json

All guards are explicit RuntimeError checks and remain active under -O. `field_bridge.json` records every basis, modulus, mapped element and finite witness used above.
