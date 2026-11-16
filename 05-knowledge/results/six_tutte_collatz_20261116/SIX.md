# Six, prime gaps, and genuine dual pairs

## Why {3,5} is selected twice

The previous coupled degree-one classification says that b,c>1 with b+c and bc+1 powers of two are, up to order, b=2^j−1,c=2^j+1 for j≥2. Exactly one pair can consist of primes: primality of 2^j−1 forces j prime (factor the power difference), and primality of 2^j+1 forces j a power of two (factor the odd power sum). Thus j=2. Independently, 3b+1=2c holds only for j=2. The same pair is selected by primality and by one exact accelerated Collatz edge. Its swapped cubic carry is 32=3³+5 and gives 256=6³+5·8. The unswapped cube gives 1024=10³+3·8. Neither the radix-six cube nor its 3→8 exponent arrow is ordinary Collatz.

The base-six wheel is exact: the two admissible rows for primes >3 are ±1 mod 6, exchanged by negation; the odd nonunit row 3 mod 6 is fixed. Under h=2n+1 the three rows are 3,7,11 mod 12; reflection h↦2−h swaps 3 and 11 and fixes 7. The row of 21 is 9, deliberately *outside* this odd-input triple. The same wheel divides Kohl's doubling matchings by whether the lower endpoint is 1 or 2 mod 3. Six is the joint period of parity and sign mod 3, not a primality test.

## Exact prime-gap arrows, including gap six

Define A_s(n)=(3n+s)/2^v₂(3n+s) on positive odd n, s=±1. For every even gap g>0 an increasing arrow n→n+g necessarily has valuation 1 and n=2g−s, since for valuation ≥2 the image is ≤n. A decreasing arrow n→n−g has valuation k≥2 and n=g+(3g+s)/(2^k−3). Thus all downward arrows of fixed gap are read from a finite divisor list, with parity, positivity and primality checked at the end.

* Gap 2, both endpoints prime: plus 3→5; minus 5→7 and 7→5. The latter makes an exact cycle. Its uncompressed minus sheet is 5→14→7→20→10→5; in the count chart this is 11→29→15→41→21→11. The edge 21→11 is shared with the plus sheet.
* Gap 6, both endpoints prime: plus 11→17; minus 13→19 and 23→17. These are all, not heuristics about the infinitude of sexy primes.

Even if we may change sign on each step, the only 2-cycle on distinct positive odd vertices is 5↔7 in the minus branch. If the ascending step uses s and the descending step uses t, the equation is (2^k−3)(2g−s)=3g+t. For k≥3 the left side is already larger; for k=2 it gives g=−s−t, so g=2 and both signs are −1.

## Three true involutions, and where their identification fails

Tutte's regular-matroid forbidden list {F₇,F₇*,U₂,₄} is one size-two dual orbit and one fixed object. Kuratowski's K₅, K₃,₃ is a different theorem, and containment in a third graph is not duality. The matching abstraction is a set with an involution:
1. On sign labels {+,0,−}, negation has cycle type 2+1; via x↦−x it exchanges the two nonzero Collatz rules, while the neutral *label* is fixed.
2. For numbers b,c>0 the squares of geometric, arithmetic and quadratic means are an arithmetic progression centered at A²: each gap is ((b−c)/2)². Reflection of the **values** swaps the outer squares, not two physical positive-input pairs. At {3,5} the triple is precisely 15,16,17, with center the shared Collatz carry. For b,c=2^j±1 it is 2^(2j)−1,2^(2j),2^(2j)+1. The two outside values have different factor arithmetic despite reflection symmetry.
3. The four tournament classes on four vertices have Hamiltonian counts 1,3,3,5. Converse swaps source-over-triangle with triangle-over-sink and fixes *both* the transitive and the strong class. Restrict to nontransitive classes for exactly 2+1. The 3→5 rise even has a physical realization by extending a directed triangle to the strong four-vertex class. The archived Collatz construction identifies the converse pair by time reversal, not by negating the sign sheet alone.

Another finite model has the actual numbers 7,21,63. In (F₂³⊕F₂³)\{0}, the two coordinate axes have 7 nonzero vectors each, exchanged by swapping the coordinates. Among the other vectors are 21 pairs of nonzero v,w with v·w=0, the point-line flags of the Fano plane, and 28 with v·w=1. The total is 63. The 21-element stratum is fixed *setwise*, not pointwise: three flags with v=w are fixed and the rest pair off. This is a combinatorial realization of the numbers, not an isomorphism to forbidden tournament counts, and Fano polarity here is distinct from matroid duality.

Prior repository lanes: [six-row Collatz](../collatz_mod6_20260917_row_braid_typing.md), [sign and Kempe triples](../procgen_kuratowski_20260925_tait_kempe_triples.md), [four-vertex reading](../procgen_tourn_20260924_four_vertex_sheets_redei.md), and [the exact decimal/base-six packet](../decimal_1024_extensions_20261116/README.md). For Tutte's theorem see his original *A homotopy theorem for matroids, I, II*, and [Seymour's survey account](https://doi.org/10.1016/j.disc.2004.07.039).

## The two cubes share an exact six-beat clock at 7

Write A_n=(128^n−125^n)/3 and B_n=(32^n−27^n)/5. The two ratios 128/125 and 32/27 reduce to the inverse generators 5,3 of the six-element unit group mod 7. Their product is (16/15)^3, and 16³−15³=721=7·103: the relation uses the same second carry 16=3·5+1. Both sixth powers differ from 1 by valuation exactly 1 at 7 (the ratios are 12 and 3 mod 49, with sixth powers 22 and 43). Hence for both ladders the valuation at 7 is zero unless 6|n and then equals 1+v₇(n/6). At 3 the ladders differ: v₃(A_n)=v₃(n), whereas 3 never divides B_n. Thus on A_n the beat for 21 equals that for 7, and 63 divides it exactly when 18|n; on B_n the 21-beat never occurs. These are divisibility events, not a count-spectrum explanation. Compare the older [Paley/Fano lane](../collatz_mod6_20260922_paley_fano_octonion_design.md) and the [gap-six local sieve](../../hypotheses/HYP-3086-sexy-prime-pair-sieve-transfer.md).
