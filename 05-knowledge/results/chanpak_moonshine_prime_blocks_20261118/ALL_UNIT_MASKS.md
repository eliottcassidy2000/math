# Allow every unit-response mask: an exact larger physical bank

**PROVED parametrization and all-order bank size. FINITE-EXACT enlarged fixed-continuation multiplicities.** No universal absolute placement is claimed.

Keep the physical roots L=(9,1,23,5), R=(9,1,25,5). For r=4j+2 set B=2^r, n=r+1 and s=(4B-1)/5. Write rev_n for reversal of exactly n bits. Every pair of n-bit low addresses whose coordinated edit has remote response coefficient A=1 is uniquely

    a=rev_n(x), b=rev_n(x+s), 0<=x<2B-s.

Indeed A=2^(r+1)(25-23)+5(rev_n(a)-rev_n(b))=4B-5s=1. This parametrization admits ALL signed masks, including masks differing from the earlier recursively fixed one.

An additional common bit epsilon at position n is allowed. Put a'=a+2B epsilon, b'=b+2B epsilon and

    C=23a'-25b'-5[K(a')-K(b')]-(a' mod2)+(b' mod2).

For every high word address V, the two physical output counts at addresses a'+4BV and b'+4BV differ by exactly C-2V. The nonlinear high-high terms cancel, and the high-low coefficient is exactly A=1. Thus C>=0 permits the unique collision continuation V=C/2. C is always even because both zero-tail counts are odd.

## A larger all-order family with guaranteed legal collision tails

For BOTH epsilon choices, C>0 whenever x is odd. Restricting to these explicitly selected odd addresses gives precisely

    2 * #(odd x<2B-s) = (6B-4)/5

collision address pairs at every level, compared with only 2^(2j+2) in the previous fixed-mask bank. This counts word pairs, not an unproved all-order count of distinct heights.

Here is a proof avoiding finite experimentation. For odd x, a has top bit1 and b top bit0. Let q=(s+1)/2 and u=rev_r(q)=(7*16^j-2)/5. Then

    a-b >= B-u.

To prove this, write the lower signed differences as z_i. Their reversed weighted sum is -q, so the digits -z_i form a signed binary expansion of q. Among all expansions with digits {-1,0,1}, the ordinary nonnegative binary expansion maximizes the reversed weighted digit sum. Normalize the highest negative digit by borrowing from the nearest1 above it: the pattern (1,0,...,0,-1) becomes (0,1,...,1,1), increasing that reversed sum by 3*2^i*(2^m-1)>0. The requisite1 exists because otherwise the represented number would be negative. Repeating removes all negative digits. Therefore the lower reversed digit sum is at most u, proving the bound.

For epsilon=0, use 0<=K(w)<=w and b<B to obtain

    C >=18a-25b-1 >=11B-18u+6>6.

Adding the common epsilon bit reduces C by exactly6 when the highest edit is +1. Hence both continuations remain positive.

Finally 2B-s=(6B+1)/5 is odd, proving the selected-bank count. Exact tests at r=2,6,10,14 found C<0 at every even x, but that exclusion is NOT claimed as an all-order theorem here. The complete all-mask decoder may simply evaluate and reject negative C directly.

## Target-first selection and center freedom

Compute C from each odd x and either epsilon, take V=C/2, and evaluate the exact output. This is a complete target-first inverse for the guaranteed odd-address bank, not merely for a chosen mask. Enumerating even x as well and retaining C>=0 also gives the general all-mask decoder, but the stated height bounds below are proved for the odd-address bank. The crude bounds 100B^2<H<2200B^2 hold throughout (use the preceding positive bound, C<46B+1, K<=W and W=a'+4BV). Successive B values differ by16, so the bands are disjoint. A prescribed N therefore selects at most one level, with O(sqrt(N)) candidate pairs and O(sqrt(N)log N) arithmetic work using direct checksums.

At r=2,6,10,14 the numbers of pairs are4,76,1228,19660. Independent exact enumeration finds all their heights distinct at those levels; this is not asserted globally. The maximal number of pairs sharing the SAME continuation V is respectively1,2,3,6. Thus varying masks concretely breaks the fixed-mask two-center ceiling, while preserving A=1 and the same two physical sources. This is genuine additional center freedom, not a different residual modulus.

The guaranteed odd-address bank still has only O(sqrt(X)) distinct centers through X by its proved pair count and height bounds, so it is not a universal selector. Further progress requires additional physical sources, more than one independent seam, or a different remote-response condition. The present result supplies a much larger exact bank and gives an exhaustive parametrization of the possible masks with A=1; the odd-address part has a proved quantitative classification.
