# Matching cores: physical binary carry and exact changing-core move

**Proved.** No all-target claim.

[THM-4736](../../../01-canon/theorems/THM-4736-the-delta-equals-k-zero-heavy-layer-exact-counts-missing-gaps-and-the-adjacent-matched-pair-collapse.md) finds backward-arc matchings as a structural boundary in its zero-heavy weighted-core classification. The one-edge antecedent is the [Fermat-number construction](../arithmetic_seams_20260921_synthesis.md), where a reversed extreme edge gives a binary mode.

## Exact marked-endpoint compiler

Start with a totally ordered set of marked endpoints; reverse an arbitrary set of pairs among them. Insert transitive runs of unmarked vertices in chosen gaps of the order, with every arc touching an unmarked vertex still forward. For every permutation of the marked endpoints (adjoin bottom and top sentinels):

1. A descending slot is possible only if its pair was reversed, and must contain no unmarked vertex.
2. In an ascending slot, an eligible run is one strictly between the endpoints. Any subset assigned to that slot occurs in increasing order.
3. If the ascending slot's pair was reversed, its intervening list must be nonempty. Enforce this by inclusion–exclusion over those slots.

After forbidding the slots in an inclusion–exclusion term, each vertex in gap g independently picks one of its m_g remaining eligible slots. A width w_g contributes m_g^w_g, with 0^0=1. Concatenation is the inverse of this assignment: the ordered endpoints and each vertex's slot uniquely recover the path. This is an exact signed monomial compiler independent of widths.

For r matched reversals from all r right endpoints into all r left endpoints across one central gap of width b, and no other unmarked vertices, the highest mode is r!(r+1)^b. An endpoint permutation needs r descending jumps to have r+1 eligible ascending slots for the gap: high-to-low steps and gap-eligible slots must alternate, beginning and ending with eligible slots. It traverses each of the r matched reversals once in any order; that order determines the permutation. These ascending slots are allowed even when empty: their successive right endpoint cannot be the one matched to their left endpoint since that right endpoint has already been visited; sentinels are unrestricted. The lower coefficients retain pairing: at r=3 the identity matching gives 6·4^b+20·3^b+5·2^b−2, while reverse matching gives 6·4^b+28·3^b+13·2^b−10, both for b>0. Leading-mode agreement is not contextual equivalence.

## Two matched reversals, three gaps

Let endpoints 0<1<2<3 have gaps a,b,c between consecutive endpoints. Write X for reversing (0,2),(1,3) and N for reversing (0,3),(1,2). Both are strong for a,b,c≥1. In X the overlapping backward arcs close two forward chains that meet; in N the extreme backward arc closes an increasing Hamiltonian path.

The entire signed monomial table (triples are insertion-slot counts for the three gaps):
    
    X: (1,1,1):1 (1,2,2):2 (0,1,2):-1 (2,2,1):2 (2,1,0):-1 (2,3,2):2
    N: (1,1,1):1 (1,0,1):-1 (1,2,1):1 (2,2,2):4 (2,1,2):-2 (2,3,2):2

It gives

    H_X = 2^(a+c+1) 3^b + (2^(a+1)+2^(c+1)) 2^b + 1 − [a=0] 2^c − [c=0] 2^a,
    H_N = 2^(a+c+1) 3^b + (2^(a+c+2)+1) 2^b + 1 − 2^(a+c+1) − [b=0].

At a=c=0 and b>0, H_N−H_X=2^b. Editing four endpoint arcs cancels the ternary term exactly and isolates a physical binary carry. Both endpoints are actual tournaments on the same named vertices.

## Exact 3H−2 via a width exchange

For u≥2 and b≥1, change X(u,b,u+1) into X(u−1,b+1,u+2). The outer product 2^a 2^c stays fixed, while the sum 2^a+2^c is multiplied by 3/2. Central growth multiplies the ternary mode by 3 and the binary mode by 2, so both nonconstant terms are multiplied by 3. The intercept remains +1. Hence H' = 3H−2 exactly. The edit may move an old gap vertex across a marked endpoint; no induced-extension claim is made.

The obstruction to iterating this particular exchange is exact. Keep the outer sum fixed and increase b by one. With widths sorted a≤c, let d=c−a. A redistribution that scales the sum of powers by 3/2 must increase the difference to d+2j with j≥1. The ratio is (2^(-j)+2^(d+j))/(1+2^d). At j=1 equality to 3/2 holds precisely at d=1. At j≥2 the ratio is at least (1/4+4)/2=17/8>3/2: first minimize over d≥0, then j≥2. The output d=3 cannot repeat the same rule. The next operation has to change the endpoint palette, the ternary coefficient, or the target representation.

The inherited [necklace/clock work](../collatz_necklace_20260929_fair_splits_power_clocks_basins.md) retains the growth clock and ordered carry separately. Matching cores exhibit two distinct clocks inside an actual Hamiltonian count; a relative switch annihilates the leading mode without removing it at either endpoint. This suggests a precise test for larger cores: find matching changes whose coefficient difference has prescribed support, then arrange joint scaling rather than fitting total counts. No ordinary Collatz claim follows.
