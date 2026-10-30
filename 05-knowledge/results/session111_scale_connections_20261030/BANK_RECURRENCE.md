# A closed recurrence from the actual shared-quotient action

**PROVED all k≥2; FINITE-EXACT independent raw-ring controls at k=2,3.**

This extends [ACTION_CARRIERS](ACTION_CARRIERS.md) beyond the accidental small-size factorization. It counts all orientation-bank states modulo a specified odd subgroup, not prime graph isomorphism classes.

## The carrier

Use R_k=(Z/5^k)[δ]/(δ²−5), H=<μ⁴>, and pair antipodal H-orbits into 5^k−1 orientation coordinates. The Sylow5 subgroup of the effective normalizer is
O_k=C_(5^(k−1))×C5.
Choose generators multiplication by 6 and α:δ↦(1+5^(k−1))δ. Write an element as(a,b), with a modulo n=5^(k−1), b modulo 5.

An odd group acting affinely on a binary cube always has a fixed origin: sum all translates of any origin. Relative to it the action is unsigned. Let B_k(Q) be the number of Q-coloring orbits of these coordinate positions. For Q=2 this is exactly the number of bank-state orbits under O_k. For Q=2^r it counts r jointly named bank assignments in that gauge.

## Valuation layers retain the shared coordinate

For valuation d, put l=2k−d. The number of antipodal orbit coordinates in that layer is
2·5^(ceil(l/2)−1).
Scalar multiplication by 6 has two orbits, each of length 5^(ceil(l/2)−1). To see the length, compare6^a with an H-power after dividing a vector of valuation d: the valuation of 6^a−1 is even, while a nontrivial clock displacement has odd valuation. Equality before both vanish is impossible.

The automorphism α acts trivially on clock-orbit coordinates except at d=1:
- at d=0 its top-socle change is a clock displacement of a unit;
- at d≥2 the change vanishes;
- at d=1 the action equals scalar multiplication by 1+5^(k−1)=6^(5^(k−2)) modulo 5^k.

Therefore the full coordinate action has:
- four cyclic quotient orbits of length 5^j for every 0≤j≤k−2, translated by a;
- two top orbits of length n, translated by a;
- two more top orbits of length n, translated by a+(n/5)b.

At k=2 the top coordinates are independent because a and a+b range freely over F5². At k=3, a and a+5b share their residue modulo 5. Treating them as independent loses that alignment and gives the wrong count. At general k the pair is the fibre product
C_n ×_(C_(n/5)) C_n.
The lower layers retain the same common coordinate a.

## Cycle count and recurrence

With gcd(0,m)=m, the exact cycle count is
C_k(a,b)=4 Σ_(j=0)^(k−2) gcd(a,5^j)
          +2 gcd(a,n)+2 gcd(a+(n/5)b,n).

Define Ψ_k(Q)=Σ_(a,b) Q^(C_k(a,b)), so B_k(Q)=Ψ_k(Q)/5^k. The base is
Ψ_2(Q)=Q^24+8Q^16+16Q^8.

At level k+1, split a according to divisibility by 5.
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

## General simply ramified metallic family

The mechanism is not special to 5. Let p≥5 be prime, μ²=cμ+1, v_p(c²+4)=1, δ=2μ−c and g=μ⁴. Such metallic ramification forces p≡1 modulo4. Work in R_k=O/(p^k), k≥2, with the inherited valuation law v_δ(g^j−1)=1+2v_p(j).

Use scalar 1+p and α:δ↦(1+p^(k−1))δ. They generate the effective Sylow-p subgroup C_(p^(k−1))×Cp. This need not be the full odd-order subgroup: at p=13, scalar units also have a C3 factor. Only the specified Sylow-p action is counted. Put n=p^(k−1) and h=(p−1)/2.

The valuation-layer proof above changes only multiplicities: there are p−1 lower quotient orbits at each length p^j, and h copies of each of the two top quotient orbits. Their exact cycle count is
C_(p,k)(a,b)=(p−1) Σ_(j=0)^(k−2) gcd(a,p^j)
             +h[gcd(a,n)+gcd(a+(n/p)b,n)].

For completeness, write a nonzero vector as δ^d u. Ring units split as a scalar unit times an H-power, because their intersection is trivial and their orders multiply to the full unit count. Applying α changes the orbit of that vector by the scalar (1+p^(k−1))^d. At d=0 this is trivial; at d=1 it is the stated top translation; for d≥2 it acts trivially on the shorter scalar quotient. Comparing the even scalar-displacement valuation with the odd clock-displacement valuation proves the quotient length p^(ceil((2k−d)/2)−1). The number of antipodal coordinates is h times that length. This establishes the entire joint action.

With B_(p,k) its Q-color orbit count, the base and recurrence are
B_(p,2)(Q)=Q^(p−1)[Q^(hp)+(p−1)Q^h]^2/p²,
B_(p,k+1)(Q)=[Q^(p−1)B_(p,k)(Q^p)+(p−1)Q^((p−1)(k+1))]/p.

The proof splits a into multiples of p and units exactly as above. Multiples give p−1+p C_(p,k); units give (p−1)(k+1) cycles and occur in (p−1)p^k group elements. All counts are for the named odd action; no full graph quotient is silently substituted.

[general_bank_tower.js](general_bank_tower.js) checks literal signed actions for (p,c,k)=(5,1,2),(5,1,3),(13,3,2),(17,8,2). The last two histograms are:
- p=13: cycle counts168,96,24 with multiplicities1,24,144.
- p=17: cycle counts288,160,32 with multiplicities1,32,256.

The statement excludes non-simple ramification, the k=1 base and p=3. No valuation argument in these excluded regimes is assumed.
