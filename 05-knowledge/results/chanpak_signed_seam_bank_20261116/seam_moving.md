# A mirrored seam can support a bank of distinct signed repairs

These are exact tournament constructions for all stated heights. They do not give a repeatable unit move or complete the odd spectrum.

## Physical geometry and high-count cancellation

Let D be any m-vertex tournament with mark V. Put the same mark on its converse D*, using the same vertex labels. Place these as two live blocks in an ordinal tape with a transitive spacer TT_c between, c≥0, and transitive outer guards TT_a and TT_b, both nonempty. Add sentinel t beaten by the whole tape and control x with t→x and x→tape except for the two marks, which beat x. The result is strong for all heights: x first reaches TT_a, and every tape vertex reaches x via t.

On (D,V), denote the usual ordered optional two-slot counts by h,q,E,L,J, where E marks the first end, L marks the second start, and J marks both. Denote by T,S the numbers of whole Hamilton paths ending and starting in V. Reversal preserves h,q,J and swaps E with L, T with S. The [two-live formula](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/chanpak_computational_tapes_20261115/HEAD.md#multiple-marks-an-exact-tape-and-a-saturated-quaternary-address) gives
```text
H = q²·2^(a+b+c) +(hE+2^c qL)(2^a−2^b)+h²−N_c(V),
N_c(V)=2hJ + 2^c L² + 1_{c=0} T².
```
The forward seam contact (a pre-x path ends in the left mark; a post-x path starts in the right mark) leaves no legal slot for a nonempty intervening spacer, hence the indicator. The reverse seam contact leaves every spacer vertex free to choose its slot, hence 2^c. All other terms follow from ordinal multiplication of h and q applied to each slot and the two single-block J corrections.

For equal outer guards, both large signed contact terms disappear in the total, even with E,L,T,S moving in both live blocks. A fixed mirror pair exposes all N_c(V) directly at a high common baseline. The two sides are physical converses, not merely sources with equal scalar h.

## A fixed seam with ±2, ±6 and ±8 repairs

Use m=7 and lex-forward code 1149404 (in lex edge order one means smaller beats larger). Both live core arc sets, all guards, and c=0 stay fixed. This is a mirror pair D,D*; it need not use an anti-automorphism of D. The mask bit index is the vertex index. The complete ordinary profile has h=59,q=940:

| Mask | E | L | J | T | S | N_0 |
|---:|---:|---:|---:|---:|---:|---:|
| 50 | 403 | 364 | 114 | 28 | 23 | 146732 |
| 11 | 356 | 365 | 110 | 23 | 24 | 146734 |
| 38 | 313 | 370 | 80 | 20 | 23 | 146740 |

At a=b=s≥1 the three strong outputs are
```text
883600·4^s+3481−{146740,146734,146732}.
```
At s=1 they are 3,391,141; 3,391,147; 3,391,149. Thus this one physical seam supports differences 2,6,8 with neither the absolute H nor the moving local observations frozen. At most eight reversals of arcs incident to x implement any pair transition. This is a bank of alternatives, not an additive semigroup of independent actions.

The low corrections contain competing terms. For example masks 50→11 give energy changes `2hΔJ=-472`, `ΔL²=729`, `ΔT²=-255`, summing to 2 (thus H decreases by 2). The 11→38 changes are −3540,+3675,−129, summing to 6. The same seam supports two genuinely different cancellations, not a duplicated primitive repair.

## A retunable separator with the same spectator budget

Use m=4, lex-forward code 4, h=5,q=32. Fix two unmarked extra transitive vertices. Either put both between the pair (c=2,a=b=s), or put one at the inner end of each outer guard (c=0,a=b=s+1). The total order and every sufficiently remote arc stay fixed, and the leading coefficient is 4096·4^s. Together with symmetric mask edits inside D,D* this gives:

| Mask | c | E | L | J | T | S | N_c |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 9 | 0 | 12 | 12 | 4 | 2 | 2 | 188 |
| 8 | 2 | 5 | 7 | 0 | 1 | 1 | 196 |
| 10 | 0 | 14 | 13 | 2 | 3 | 2 | 198 |

At every s≥1 the outputs `4096·4^s+25−{198,196,188}` support signed differences 2,8,10. A change of c switches the actual available forward seam contact rather than scaling the whole previous correction. This is the precise extra flexibility purchased by moving two pivots around the same two cores.

## Why a straight endpoint walk cannot give a long unit run here

Assume a=b and c fixed. If consecutive masks keep common h,q and make L,T arithmetic functions of a step index, with increments δ,ε, then the second difference of H is
```text
Δ²H = −2h Δ²J − 2(2^c δ² + 1_{c=0} ε²).
```
Thus three arithmetic outputs require `h | (2^c δ² + 1_{c=0} ε²)`. For c=0, an odd prime p≡3 mod4 dividing h must divide both δ and ε, because −1 is not a quadratic residue modulo p. For c>0, every odd prime p dividing h must divide δ. In particular the h=59 seam cannot have three equally spaced outputs from affine, nonzero small increments in its endpoint pair. Neither arbitrary masks nor full growing-core methods are ruled out. Either the endpoint motion itself must curve, or one must change the seam model. The physical h=59 bank illustrates varying second differences explicitly.

## How the historical random search guided the proof

A deterministic random sample of seven older theorem/result reports (seed 131316) included the October [balanced-deletion/two-cut-port packet](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/balanced_deletion_groups_ports_20261021/README.md) and September [Collatz order-versus-crossing laws](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_procgen_20260922_order_laws.md). The former distinguishes exact labelled seam corrections from ordinary projected counts; the latter distinguishes a common high-rate law from boundary crossings. Here the T² indicator is an exact seam-port correction and the second-difference obstruction formalizes why an apparently simple family cannot supply repeated arithmetic outputs. No Collatz proposition is used in the proof.

The check enumerates independent subset-DP endpoint profiles, explicitly builds the complete strong tournaments at s=1 and verifies the two exact banks. The all-height assertion comes from the formula, not extrapolating the tests.

## The twisted seam: polarization and a two-mask impossibility

The right block D* may carry W while the left carries V, rather than copying V. For equal guard lengths a=b and fixed c, high cancellation now requires equality of one signed shelf coordinate
```text
ρ_c(V)=ρ_c(W), where ρ_c(V)=hE(V)−2^c qL(V).
```
This follows from the two high terms' signed coefficient `hE(V)+2^c qL(W)−2^c qL(V)−hE(W)`. The exact remaining energy is the symmetric bilinear expression
```text
N_c(V,W)=h(J(V)+J(W))+2^c L(V)L(W)+1_{c=0} T(V)T(W).
```
Let u_c carry coordinates (2^(c/2)L, 1_{c=0}T). The twist beats the average of its diagonal endpoints by a squared distance:
```text
H(V,W) = [H(V,V)+H(W,W)+||u_c(V)−u_c(W)||²]/2.
```

**Two masks alone never yield three consecutive odds on a nontransitive D.** Indeed the only three possible distinct outputs are two diagonal ones and the symmetric cross output. If the diagonals are the ends of a three-odd interval, their gap is 4 and the cross must be their midpoint, requiring squared distance 0. Then their gap would be `2h·|ΔJ|=4`, impossible for odd h>1. If the diagonals are consecutive with gap 2, the cross cannot be below either (it is at least their average); to finish the three-odd interval above them it needs squared distance 6. At c=0 this is impossible as a sum of two integer squares: both coordinates would have to be divisible by 3. At c>0 it is impossible as a single square times 2^c. Equal diagonals cannot give three distinct counts.

**Independent controls are an orthogonality test, not just endpoint preservation.** Four states made by varying V on the left and W on the right inside one fixed high shelf have exact signed mixed defect
```text
−[2^c ΔL_left ΔL_right + 1_{c=0} ΔT_left ΔT_right].
```
Thus with no separator, commuting internal choices require orthogonal moves in two real coordinates. With a positive spacer the metric drops rank to one. One direction must then leave L unchanged for two edits to commute, and its count change is a multiple of h, or zero. In particular it cannot be a unit (+2 or −2) in this setting when h>1. This is an exact in-model obstruction to manufacturing a commuting pair of unit switches by widening an unmarked separator. It leaves asymmetric multi-mask banks and changes of porthole geometry open.
