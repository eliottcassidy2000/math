# Exact cyclic divisor test and changing-core unit rays
Companion to [the five-vertex two-reservoir construction](README.md).
The independent literal replay is `verify_cyclic.py`.

Let D,E be any nonempty tournaments of profiles ((h,A)=(H(D),P_2(D))) and ((k,B)=(H(E),P_2(E))), where (P_2) counts two ordered paths allowing emptiness. Append respective transitive bags of widths m,n by ordinal sum. Use the expanded D, expanded E, and a singleton as three modules with the three cyclic arc directions. The host is strong for any m,n≥0.

With one center vertex, the two nontrivial modules can each occur either once or twice in a compressed Hamiltonian path. Both once has three cyclic skeletons; each of the two single repeats and the double repeat has one. Nonempty two-run profiles after padding are (2^mA-2h) and (2^nB-2k). Filling the compressed runs gives the following exact multiplicative shift:
 
    H = (2^m A − h)(2^n B − k) + 2hk.                 (1)

The optional-slot multiplication by (2^m) follows by independently assigning each padding vertex to an ordered slot. This law is a rotation template, not an arbitrary endpoint concatenation.

## Absolute inverse for a selected pair of cores

Given a target N and physical candidate profiles, reject when (N-2hk≤0). Otherwise factor (N-2hk) and test each oriented factor pair (u,v). There is a literal witness in this pair of core families exactly when ((u+h)/A) and ((v+k)/B) are positive integer powers of 2 (including (2^0)). These determine both widths exactly. This is an arithmetic inverse, not a claim that some candidate profile pair always passes. A finite bank contributes only (O((\log X)^2)) distinct counts through X, because the two factors in (1) are positive and each grows exponentially in its height.

## Why two deep widths cannot yield fixed +2 with a fixed bank

For any core pair, (1) implies (H≡3hk\pmod {2^{\min(m,n)}}). Compare any two families that may use different cores and different heights. If the minimum of their four widths is L, their signed gap g satisfies

    g ≡ 3(hk − h'k') mod 2^L.

For finite core banks the right-hand integer difference is bounded uniformly, so for any fixed signed gap g, once L becomes sufficiently large this congruence must be literal equality. Since all source h,k are odd, every such fixed gap is divisible by 6. In particular choosing new cores inside even a large fixed bank cannot sustain +2 when both widths on both sides tend to infinity. The same proof provides the exact numerical obstruction: if (g≠3(hk-h'k')) then (2^L≤|g-3(hk-h'k')|).

A nonzero constant gap in this regime is nevertheless physically possible. Three order-eight cores with common P₂=2176 and counts 79,85,91 have outgoing rows

    (2,20,17,7,105,15,175,63)
    (16,5,65,7,46,15,59,127)
    (8,5,17,22,67,31,175,63).

Compare the center pair (85,85) with the outer pair (79,91), putting the same width m on both modules in both hosts. The leading terms and the linear cross terms in (2^m) agree because (79+91=2·85), leaving exactly (3(85^2-79·91)=108) at every m. Both hosts are strong, of equal order (2m+17).

## Actual unit rays by holding just one width bounded

When n=0, the first-height output is (2^mA(B-k)+h(3k-B)). A matched slope with different constant terms gives exact ongoing repairs.

1. Take (D=(2,4,9,19,39,15)), (E=TT_2). Their profiles are (17,224),(1,4), giving (672·2^m-17). Taking (D'=(0,1,3,23,39,15)), (E'=TT_3) gives profiles (3,96),(1,8) and (672·2^m-15). The orders are m+9 and m+10.
2. A same-order strong pair: lower (D=TT_4,E=(4,9,26,81,163,15,55,111)) has profiles (1,16),(195,3720); upper (D'=(128,5,9,19,39,79,31,126), E'=TT_4) has profiles (241,3760),(1,16). Both give slope 56400 and respective intercepts −3135 and −3133. The two resulting counts are at order m+13. The first pair at that order begins 53,265 / 53,267.

Here one height is a physically fixed core, not an external approximation. Membership in each ray is exactly decidable by subtracting its intercept, dividing by its slope, and testing whether the result is a positive power of 2. No ordinary one-arc adjacency is asserted for these changing-core pairs. These rays do not by themselves cover every target.

