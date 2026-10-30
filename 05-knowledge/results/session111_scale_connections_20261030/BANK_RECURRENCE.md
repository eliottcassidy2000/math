# A closed recurrence from the actual shared-quotient action

**PROVED all k≥2; FINITE-EXACT independent raw-ring controls at k=2,3.**

This extends [ACTION_CARRIERS](ACTION_CARRIERS.md) beyond the accidental small-size factorization. It counts all orientation-bank states modulo a specified odd subgroup, not prime graph isomorphism classes.

## The carrier

Use R_k=(Z/5^k)[δ]/(δ²−5), H=<μ⁴>, and pair antipodal H-orbits into5^k−1 orientation coordinates. The Sylow5 subgroup of the effective normalizer is
O_k=C_(5^(k−1))×C5.
Choose generators multiplication by6 and α:δ↦(1+5^(k−1))δ. Write an element as(a,b), with a modulo n=5^(k−1), b modulo5.

An odd group acting affinely on a binary cube always has a fixed origin: sum all translates of any origin. Relative to it the action is unsigned. Let B_k(Q) be the number of Q-coloring orbits of these coordinate positions. For Q=2 this is exactly the number of bank-state orbits under O_k. For Q=2^r it counts r jointly named bank assignments in that gauge.

## Valuation layers retain the shared coordinate

For valuation d, put l=2k−d. The number of antipodal orbit coordinates in that layer is
2·5^(ceil(l/2)−1).
Scalar multiplication by6 has two orbits, each of length5^(ceil(l/2)−1). To see the length, compare6^a with an H-power after dividing a vector of valuation d: the valuation of6^a−1 is even, while a nontrivial clock displacement has odd valuation. Equality before both vanish is impossible.

The automorphism α acts trivially on clock-orbit coordinates except at d=1:
- at d=0 its top-socle change is a clock displacement of a unit;
- at d≥2 the change vanishes;
- at d=1 the action equals scalar multiplication by1+5^(k−1)=6^(5^(k−2)) modulo5^k.

Therefore the full coordinate action has:
- four cyclic quotient orbits of length5^j for every0≤j≤k−2, translated by a;
- two top orbits of length n, translated by a;
- two more top orbits of length n, translated by a+(n/5)b.

At k=2 the top coordinates are independent because a and a+b range freely over F5². At k=3, a and a+5b share their residue modulo5. Treating them as independent loses that alignment and gives the wrong count. At general k the pair is the fibre product
C_n ×_(C_(n/5)) C_n.
The lower layers retain the same common coordinate a.

## Cycle count and recurrence

With gcd(0,m)=m, the exact cycle count is
C_k(a,b)=4 Σ_(j=0)^(k−2) gcd(a,5^j)
          +2 gcd(a,n)+2 gcd(a+(n/5)b,n).

Define Ψ_k(Q)=Σ_(a,b) Q^(C_k(a,b)), so B_k(Q)=Ψ_k(Q)/5^k. The base is
Ψ_2(Q)=Q^24+8Q^16+16Q^8.

At level k+1, split a according to divisibility by5.
If a=5a', the cycle formula gives
C_(k+1)(5a',b)=4+5 C_k(a',b).
The other a are units; all gcd terms are1, hence their cycle count is4k+4. There are4·5^k=20·5^(k−1) such pairs(a,b). Consequently
Ψ_(k+1)(Q)=Q^4 Ψ_k(Q^5)+20·5^(k−1)Q^(4k+4),
and the desired short recurrence is

B_(k+1)(Q)=[Q^4 B_k(Q^5)+4Q^(4k+4)]/5,  k≥2,

with B_2(Q)=(Q^24+8Q^16+16Q^8)/25.

At k=3:
B_3(Q)=(Q^124+8Q^84+16Q^44+100Q^12)/125.
Thus
B_2(2)=692224,
B_3(2)=170141183461707171770974935790600192.

The recurrence does not start at k=1: the effective normalizer has a different base structure there. Substituting B_1(Q)=Q^4 into the displayed rule gives the wrong B_2.

## Why this is useful

The25/100/200 investigation has produced an all-order count without enumerating2^(5^k−1) bank states. The small factorization at k=2 survives at higher levels only as a common-quotient constraint. This is a precise instance of increasing alignment memory: separate top marginals are insufficient, but a shared residue is enough for this operation.

This does not yet count prime bank orbits or full graph classes. Those require retaining ideal-cut conditions and full transporters. The recurrence is a reusable action-level component toward that larger enumeration problem.

## Verification

[bank_tower.js](bank_tower.js) builds the literal rings and clock orbits at k=2,3, checks every element of O_k against the cycle formula, and verifies sign consistency. Two independent implementations agree. It also checks formal recurrence mass and integral orbit counts at Q=2,3 through k=8. Those finite checks support, but do not replace, the valuation-layer proof.
