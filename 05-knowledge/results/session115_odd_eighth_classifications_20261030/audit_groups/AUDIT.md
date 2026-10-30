# Independent audit of affine-cube and Sylow-D8 eighths

**Verdict: PASS for the stated family classifications.** This is an independent written-proof review and literal permutation computation, not a formal proof-assistant certification or a literature-priority assessment.

## Affine actions in arbitrary dimension

For an odd-order linear transformation h of order l over F2, the norm 1+h+...+h^(l−1) maps onto Fix(h) and acts identically there. An affine element (v,h) has odd order exactly when this norm annihilates v. Thus its conditional odd-order proportion is the reciprocal of the fixed-space size. Averaging over the actual complement elements gives the claimed uniform formula. The complement's order alone does not replace its fixed-space profile.

The dimension-three completeness proof is sound. Lagrange restricts odd subgroup orders of GL3(2) to1,3,7,21. The order-three rational canonical form and Sylow conjugacy handle3 and7. A group of order21 has a normal Sylow7 subgroup, and the Singer normalizer has order21, so every such group is exactly a Singer normalizer. This gives one conjugacy type per listed order. The even-order masses are respectively7/8,5/8,1/8,3/8.

The exact decoration lists also check: for complement orders1,3,7,21 they are respectively seven copies of1/8; 3/8 plus two copies of1/8; one1/8; and one3/8. The order-three complement therefore supplies a genuine hostile to replacing every centralizer's odd-order proportion by the reciprocal of its Sylow2 order. The Sylow-D8 proof cannot skip its extra argument.

## Universal Sylow-D8 theorem

The normal-complement input is correctly stated and applied. An independently checked author-hosted source is David Craven, *Fusion Systems*, Theorem1.12, printed page5 (PDF page10), which gives the normal p-complement conclusion under the central-normalizer hypothesis. The same source's Theorem1.1 verifies the abelian-Sylow fusion control used in the draft's transfer explanation. [Source](https://web.mat.bham.ac.uk/D.A.Craven/docs/lectures/fusionsystems.pdf).

For J=C_G(z), where z is the central involution of a Sylow D8, J/<z> has SylowV4. A normalizer element in this quotient lifts to an element normalizing the full preimage D8. It must preserve the unique nontrivial quotient element whose preimage is cyclicC4. The induced action therefore has order dividing2. It also has odd order, because the abelian SylowV4 is contained in its centralizer. The action is trivial, so the cited complement criterion applies. The quotient has odd-order proportion1/4. Each odd-order quotient element has exactly one odd-order lift among its two central lifts, and hence rho(J)=1/8.

For an involution not conjugate to z, its centralizer has SylowV4. Its normalizer action fixes that involution and is again both odd and a divisor of2. For an order-four element, its centralizer has cyclic SylowC4, with the same normalizer argument. Both centralizers have odd-order proportion1/4. The order-four elements form exactly one class because the two in D8 are already conjugate there.

If c is the number of involution classes, c is1,2 or3. One class is the central type. Summing the actual centralizer weights gives E(G)=1/8+(c−1)/4+1/4=(2c+1)/8. No classification of finite simple groups, conjectural fusion statement, or unproved reciprocal-centralizer shortcut is needed.

The equivalent independent Frobenius-criterion route also works: in a group with central involution and SylowD8, every Klein-four subgroup contains that involution, so all induced automorphism groups of2-subgroups are2-groups. It is not necessary for the published proof because the weaker Burnside criterion already suffices.

## Independent finite controls

`check_group_weights.py` imports no producer implementation. It constructs permutations, computes orders from cycle lengths, computes conjugacy orbits and centralizers literally, and tests:

- D8, S4, GL3(2), A6, A7, D8×C3 and S4×C3;
- every odd subgroup in GL3(2): one trivial,28 order-three, eight order-seven and eight order-twenty-one subgroups;
- all45 resulting affine groups;
- the exact four affine decoration spectra.

There are152 explicit exception gates. Normal and optimized execution produce byte-identical output in `results.json`. The gates encapsulate complete finite universes listed above; their number is not a claim of testing every finite group.

Run `python check_group_weights.py` and `python -O check_group_weights.py` from this directory. Python3 standard library only.

## Scope retained

The observable is even-order probability. A1/8 decoration weight and a1/8 total probability are different uses of the same rational. The results classify the stated group families and, conditionally, the contribution of any switching class whose full automorphism group belongs to the relevant family. They do not establish realization of these groups as switching-class automorphism groups, arbitrary switching-class enumeration, or a transport to Collatz content or LRC times.
