# Complete inversion of the ready-value images

**PROVED reduction**, with the finite heads below **FINITE-EXACT**.
Use the inherited definitions h=H(T), s=#Hamilton paths starting at r,
q=H(T-r), and F=E_r(2), the ordered two-bin covers whose first path ends
at r. A b-ready root means F-h=b. All counts are ordinary directed paths.

## O1. Weak sources do not repair a positive-H ready obstruction

Write the ordered SCC decomposition T=T1=>...=>Tt, with r in Tj. Let
P_U(2) count two ordered possibly empty spanning paths in U. Then

    F_r(T) = product_(i<j) P_Ti(2) * F_r(Tj)
             * product_(i>j) H(Ti).                         (O1)

Indeed each bin passes through the SCCs in order. Before Tj both bins
are unrestricted; inside Tj the first ends at r; afterwards only the
second can continue. Restriction and concatenation are inverse maps.
The ordinary H count is the product of all component H values.

Always P_U(2)>=2H(U) (put a Hamilton path in either bin) and F_r(U)>=H(U)
(cut a Hamilton path immediately after r). Consequently j>1 implies
F-h>=h. If0<b<h, r must lie in T1 and

    b=(F_r(T1)-H(T1))*product_(i>1) H(Ti).

For b=1 or2 every later factor is odd and hence equals1. Each later SCC
is a singleton. Thus, for h>=3, a weak b-ready tournament is exactly a
strong b-ready core followed by a transitive tail. This preserves h,b.

The h=1 boundary is different and retained: all such tournaments are
transitive, with F at position j equal to2^(j-1). There is a1-ready root
at the second vertex and no2-ready root. No statement here suppresses it.

## O2. All-order ready-value inversion has a finite, constructive head

The inherited [READY classification](../two_entry_funnels_20261103/READY.md)
has two alternatives: an indegree-one newborn of a strong parent with
exactly b starting paths, or one of the five exceptional2-ready families.
For a newborn, the parent profile satisfies

    h_child=h_parent+F_parent-b >= 2h_parent-b.        (O2)

For a requested value h, enumerate only parents with
h_parent<=(h+b)/2. Busch's **published** sharp strong-order bound
f(3q)=3*5^(q-1), f(3q+1)=5^q, f(3q+2)=9*5^(q-1) bounds the parent order.
Source: [Busch2006, Corollary1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf).
The bound and its equality examples are inherited, not new results here.

The one-start parents have the inherited rank/unrank grammar. For two
starts, delete the root: the first SCC has either one two-start entry or
two unique-start entries, and every subsequent SCC is a singleton. This
gives an explicit generator using the inherited unique-cycle core grammar
and binary tails. Canonicalize by the lesser of its two Hamilton-path
frames. Its counts agree with the proved series
`x^4(1-x)(1-x^2)/(1-3x+x^2)^2`.

`ready_inverse.py` implements this construction, not an unrestricted
graph search. Evaluating every resulting newborn and the five explicit
exception families gives a witness or a complete rejection for any
requested ready value. Work is polynomial in the numerical value h:
the order bound is O(log h), the one/two-start grammars have exponential
growth in that order, and subset counting is exponential in that order.
This is not polynomial in the bit length of h and is not an inverse for
the entire Hamilton spectrum.

## O3. H35 is outside both ready images at every size

For h35 the parent bound is18. Busch gives f(7)=25, so parent order<=6.
The complete strong one-start and two-start heads have33 and24 classes.
Filtering parent H<=18, the possible child values<=35 are

    b1: 5,9,11,13,17,19,23,29,31,33;
    b2: 9,23,29.

The five exceptional values are
`3*2^k-1,6*2^k+1,5*2^k-1,13*2^k+3,14*2^k+3`, k>=1.
None is35. This proves the all-order exclusion from the **ready** image;
O1 extends it to weak tournaments. H35 itself is realized by RETURN's
seven-vertex word `11,11,00,01`.

Both ready lanes exist at H11. The branch3*11+2=35 therefore cannot
have an H-preserving return to either ready lane, by any operation
whatsoever. This refutes closed ready-bank induction containing that
source, not the {7,21} spectrum conjecture. The earlier H17 boundary is
also useful: it has a1-ready source but no2-ready source at any order.

The complete all-order value inversion through127 uses parents only
through order8 (232 one-start and275 two-start classes). Among odd h in
[5,127], neither ready lane exists exactly at

    7,15,21,27,35,41,43,57,67,83,101,105,111,121.

The inherited [literal source bank through999](../two_port_return_induction_20261103/README.md)
realizes all the other values in this interval. Only7 and21 are spectrum
gaps in this range; the other twelve are
**interface obstructions**. The full witness bank is `inverse.json`.
The parent's admissibility bounds make this an all-order statement on
that value interval, not merely a census of small output tournaments.

## O4. What replaces the failed implication

The first failed implication is "an actual ternary output has another
ready presentation." Additional vertex names or a better clock cannot
repair the empty H35 ready fibre. The strongest survivor is a graded
return into sources with larger joint defect. RETURN realizes that
survivor at all k, while ATLAS classifies its first new neutral layer.
The new obligation is a construction acting on those graded states that
again realizes offsets2 and4. It remains OPEN.
