# Weighted equality signatures and exact CRT counts

25 October 2026. **PROVED** for all positive integers; **FINITE-EXACT** certificate ranges stated below. This is a separate counting continuation, not an alteration of the frozen `session55-crt-collision-wiring` packet. The domain is Z/nZ, square residues include zero, and objects are labeled value arrays, not choices of nine square roots.

## 1. State space and conventions

Let A_n be all 3×3 arrays with square-residue entries whose three row sums, three column sums, and two diagonal sums agree. Repetitions are allowed. Label the nine positions 0,...,8 row-major. Let Π_9 be the set of partitions of these positions. Order it by refinement: π≤σ means every π-block is contained in a σ-block. Thus the discrete partition 0̂ is minimum, the one-block partition 1̂ maximum. The meet σ∧τ consists of all nonempty intersections of a σ-block with a τ-block.

For an array A let eq(A) be its exact equality partition, and define

    M_n(π) = #{A in A_n : eq(A)=π}.
    Z_n(π) = sum over σ≥π of M_n(σ).

Z_n(π) counts arrays satisfying every equality prescribed by π; it does NOT require different blocks to have different entries. This distinction determines the direction of the zeta transform. Let N(n)=M_n(0̂), the count of nine-distinct labeled arrays.

There are Bell(9)=21,147 possible equality partitions, but most may have zero weight. Represent each partition by its exact 36-bit equality relation in lexicographic order of pairs (i,j), 0≤i<j≤8. Refinement is mask containment; meet is bitwise AND. Restricted-growth strings in the certificate give a second unambiguous encoding. Independent permutations of the nine labels are not identified.

## 2. Exact weighted CRT theorem

For gcd(a,b)=1, coordinatewise CRT gives a bijection A_a×A_b→A_ab. Line equations commute with reduction, and a residue modulo ab is a square exactly when its two reductions are squares: combine any local roots by CRT. No root multiplicity is inserted, since the objects counted are values. Two global cells agree exactly when both local pairs agree. Therefore

    M_ab(π) = sum over σ∧τ=π of M_a(σ) M_b(τ).

The formula iterates associatively over prime-power factors. For n=1, the single all-zero array gives M_1(1̂)=1, the identity for this meet convolution. In particular the total number of arrays, including repetitions, is multiplicative.

A fixed π is refined by σ∧τ if and only if π≤σ and π≤τ. Summing the convolution proves

    Z_ab(π) = Z_a(π) Z_b(π).

By Möbius inversion on the partition lattice,

    N(n) = sum over π in Π_9 of μ(0̂,π) Z_n(π),
    μ(0̂,π) = (-1)^(9−number_of_blocks(π)) product over blocks B of (|B|−1)!.

For completeness, the interval [0̂,π] is the product of partition lattices on the blocks. For a block of size k, μ(0̂,1̂)=(-1)^(k−1)(k−1)!: the partition-indexed sum of the proposed values is the sum of signs of all permutations on k letters, since each block carries (k−1)! cyclic orders. This sum is zero for k≥2 by the fixed-transposition sign-reversing bijection, and one for k=1. It is precisely the defining Möbius recurrence. Products give the displayed formula.

Consequently a reusable exact distinct-count formula is

    N(product q) = sum over π μ(0̂,π) product over prime-power q of Z_q(π).

This is an identity, not a claim of an asymptotically superior algorithm.

## 3. Lossless and lossy compression

Storing only positive M_n entries is exact: omitted entries are zero. At each multiplication form all pairwise mask intersections, multiply the weights, and ADD weights for coincident intersections. The induction above proves that this sparse support algorithm returns the entire histogram. No arbitrary partitions need be enumerated for this operation; the support is closed under the intersections that the algorithm actually produces. Enumeration of all Bell(9) states is used separately in the certificate to check the zeta and Möbius formulas.

The existence antichain compression in the frozen packet is NOT valid for counting. Discarded dominated states carry positive multiplicity and can contribute after subsequent intersections. Even retaining counts on minimal states does not recover these contributions. Likewise block-size histograms discard which shared positions collide. The frozen complementary three-block states A and B both have sizes (3,3,3), but A∧A=A whereas A∧B=0̂. Full weighted named equality relations are sufficient; minimality against all physically realizable square-array continuations is not asserted.

## 4. A factor of three: feasibility-neutral, count-doubling

Modulo 3, square residues are 0 and 1. There are exactly two admissible arrays, both constant. To see this without a search, every magic array has opposite-pair sum 2c. If c=0, each opposite pair must be (0,0). If c=1, each opposite pair must be (1,1), since the allowed sums are 0,1,2. Thus M_3(1̂)=2, with all other weights zero. For gcd(n,3)=1,

    M_3n(π)=2M_n(π),  Z_3n(π)=2Z_n(π),  N(3n)=2N(n).

The factor is neutral for existence but not for labeled counts. This is restricted to coprime CRT factors; it is not a formula for replacing 3-adic powers.

## 5. Position symmetries and scope

The dihedral group D4 acts on positions, preserves magic and squarehood, and acts freely on arrays with nine distinct entries. A nonidentity position permutation moves a cell, so invariance would force an equality between two distinct positions. Hence N(n)/8 counts D4 position-orbits. It is NOT a quotient by multiplication by square units. No independent cell permutations, square-unit scalings, or root choices enter N(n)/8. Repeated arrays need not have free D4 action.

This analysis is for quotient rings, not extension fields. It does not lift modular arrays to integer squares or solve integer magic-square existence.

## 6. Independent verification and reproducibility

`verify.py` has no nonstandard dependencies and no imports from the frozen packet. Its first enumeration uses the universal pivots (center c, top-left a, top-right b):

    a         3c−a−b     b
    c−a+b     c          c+a−b
    2c−b      a+b−c      2c−a.

The form is valid over every Z/nZ without division. Conversely, if the common sum is S, summing the four opposite pairs plus the center gives total 4S−3c, while row totals give 3S; hence S=3c. The line equations then determine the displayed array uniquely.

A second enumeration loops over the entire top row and compatible centers, reconstructs the other five entries from separate line equations, and checks all eight lines. Equality of the FULL WEIGHTED histograms, not just supports, is checked for each 1≤n≤160, and each histogram is compared with prime-power CRT meet convolution. The same independent checks cover selected products beyond this range. All 21,147 zeta coordinates and Möbius inversion are checked for products 35,55,201,36,200,77,432. Certificate examples include complete nonzero histograms for these products and both local factors. Every omitted partition has weight zero by exhaustive enumeration.

Run from this directory:

    python3 verify.py --verify
    python3 -O verify.py --verify

Runtime reports are empirical observations of this checker on the current host, not comparative performance claims. No theorem depends on Python assertions, so optimized execution retains all checks.

## 7. Actual prior all-size gluing signature and exact scope of connection

Repository source: `05-knowledge/results/session25_spectrum_wiring_20261021/session25-ports/PORT_SIGNATURES.md`, §§1–3. For arbitrary boundary size r and palette size k, that theorem fixes the SAME named boundary vertices and counts extensions f_G(π) of one named-color assignment with equality partition π. For fragments with disjoint interiors it proves

    Col_k(G union_boundary H) = sum_π (k)_{|π|} f_G(π) f_H(π).

It also proves necessity against all abstract boundary graph tests by spanning partition indicators with boundary-edge monomials. Color-name quotienting retains equality partitions; independent port-name quotienting loses attachment compatibility. This is an actual all-size counting signature, not merely a feasibility analogy.

The present source object is a square-residue magic array and target state is its weighted equality relation on the same nine cells. The operation is CRT, mapped EXACTLY to meet convolution, with distinctness the zero-collision predicate. This differs from the coloring gluing operation, which requires matching boundary assignments and has falling-factorial color weights. Our zeta transform diagonalizes the meet product; it does not identify the two products. The shared principle is that one may forget internal realizations only after retaining the full compatibility data on named ports. The destroyed information in a block-size quotient is alignment, and the cheapest hostile is A∧A versus A∧B. The needed sidecar is the labeling of all nine cells.

There is NO native tournament realization asserted. The same prior source §4 explicitly distinguishes tournament entry/exit path-run ports and says coloring counts are not identified with Hamiltonian counts. Neither the present arithmetic map nor the prior signature proof supplies such a realization, guard transfer, or metagraph theorem.
