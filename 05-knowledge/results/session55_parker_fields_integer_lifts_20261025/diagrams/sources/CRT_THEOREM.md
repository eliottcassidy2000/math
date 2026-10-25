# Exact CRT collision wiring for modular magic squares

25 October 2026. **PROVED**, with exhaustive independent checks for 2 ≤ n ≤ 160. Domain throughout: the quotient ring Z/nZ, not a field of order n. A ring is called Parker here when it has no 3×3 magic array of nine pairwise distinct square residues. Zero is an allowed square.

## 1. Universal three-pivot parametrization

For arbitrary c,a,b in Z/nZ, put

    a         3c−a−b     b
    c−a+b     c          c+a−b
    2c−b      a+b−c      2c−a.

Every line sums to 3c. Conversely every magic array has this form, without needing 2 or 3 invertible. Indeed the four opposite pairs have sum S−c; totaling all entries gives 4S−3c, while the three row sums give 3S. Subtraction yields S=3c. Taking a,b as the top corners determines the remaining entries by the line constraints. Thus enumerating (c,a,b) in Q(n)^3, where Q(n) is the square-residue set, and retaining arrays whose other six entries are squares is complete, including repeated-entry arrays.

## 2. Exact fixed-port CRT theorem

Fix row-major labels 0,…,8. List their 36 unordered pairs lexicographically. For an array M define C(M) to contain exactly the pairs of labels carrying equal entries. Encode C(M) as a 36-bit nonnegative integer. It is the equality partition as a relation, not just the multiset of block sizes.

Let n=∏q_i be the factorization into pairwise coprime prime powers. A magic array of square residues modulo n is precisely the coordinatewise CRT combination of magic square-residue arrays M_i modulo q_i. Repeats are allowed locally. Squarehood is equivalent coordinatewise: CRT-combine a chosen square root at each prime-power component. The line equations are equivalent coordinatewise as well.

Two global cells agree if and only if their local cells agree in EVERY component. Consequently

    C(CRT(M_1,…,M_r)) = C(M_1) ∩ … ∩ C(M_r),

and the global array is pairwise distinct if and only if the bitwise AND of the local masks is zero. This is an exact necessary-and-sufficient criterion and constructs a witness and its nine square roots. Cell labels must be preserved throughout; independent arbitrary permutations are not allowed. A geometric symmetry is usable only by actually applying it to an admissible local array.

## 3. Lossless dominance compression

For each local modulus retain just the inclusion-minimal realizable collision masks, with one realizing array per mask. Any discarded mask contains a retained one, and replacing it by that retained mask can only delete surviving collisions with all other components. Therefore the replacement cannot destroy a possible zero intersection. In fact the antichain of inclusion-minimal realizable global masks is computed exactly by repeated operations

    A ← Min{ X ∩ Y : X in A, Y in local_antichain },

starting at the full 36-pair relation. Retain the list of realizing local arrays at every state. The CRT witness then follows without another search. This is a positive search improvement: local tables are reusable; intermediate dominated states are deleted; no direct search over Q(n)^3 is required after the local tables exist. The worst-case antichain size is not claimed small for arbitrary n.

## 4. Orthogonal three-block partitions and an infinite sufficient family

Consider the partitions

    A = { {0,4,8}, {2,3,7}, {1,5,6} },
    B = { {2,4,6}, {0,5,7}, {1,3,8} }.

Their collision masks are respectively 0x61089888 and 0x8a054250, and A∩B is the discrete partition (bitwise AND zero). Each block of A meets each block of B in exactly one cell.

For three distinct square residues x,y,z with x+z=2y, the array

    y z x
    x y z
    z x y

is magic and realizes A. Its reflection in the vertical axis realizes B. The integer square arithmetic progression (1,25,49) remains pairwise distinct modulo any integer m that does not divide 48: the pairwise differences are 24,24,48. In particular this holds for EVERY odd m≥5.

**Corollary.** If n has coprime factors a,b≥5 that are both odd, then Z/nZ is non-Parker. More generally the same conclusion holds whenever a,b are coprime and neither divides 48. Combine the A-array modulo a with the B-array modulo b. If n has additional prime-power factors coprime to ab, append constant-zero arrays: intersecting a zero mask with any further masks remains zero. The statement for n=ab requires no extra factors. The corollary is limited to the stated coprime factorization.

Thus local Boolean Parker flags cannot classify products. The information needed is how local equalities sit on shared cell labels.

## 5. Smallest all-Parker-component counterexample

The least n whose prime-power components are all Parker but whose product is not is

    35 = 5·7.

One witness is

     0  9  1
    16 15 14
    29 21 30.

All entries are distinct modulo 35, all are squares, and every line sums to 10 modulo 35. Modulo 5 its collision partition is A; modulo 7 it is B. Both factors have exactly the minimal-mask antichain {A,B}. If instead the SAME partition A is selected at both factors, the global array still has only three distinct values, even though each local block-size histogram is the same (3,3,3). Thus even local equality block-size histograms fail; alignment is essential.

Modulo 34=2·17 every combination of the three minimal masks modulo 2 and four minimal masks modulo 17 has nonzero intersection, so 34 is Parker. Modulo 55=5·11 the complementary A,B choice works. Exhaustive full-modulus enumeration for every 2≤n<35 proves the minimality assertion, not just a search among guessed examples. Other early all-Parker-component products that become non-Parker include 45,50,55,63,65,77. The certificate contains their witnesses and exact component tables.

## 6. Verification and scope

`verify.py` computes local masks, antichains, dynamic CRT compositions and witnesses, then independently enumerates each full modulus in two ways: the universal pivot formula and a top-row/center reconstruction from individual line equations. It checks equality of their FULL mask sets and equality of their minimal masks with the CRT answer, for every n=2,…,160. Each constructed CRT witness is checked for square membership, all eight line sums, distinctness, and explicit square roots.

The theorem is exact for every n, but obtaining each local antichain is a finite computation, not a closed-form classification. The universal AP corollary is an infinite non-Parker family, not by itself a general quotient-ring classification. The additional reductions in CLASSIFICATION.md complete that classification using the separately audited prime-field result. No claim identifies F_(p^r) with Z/p^rZ, transfers p-adic lifetimes, or lifts an arbitrary modular witness to integer squares.

## 7. Exact relation to prior shared-port work

Repository source: `session51-blue-defect-wiring/sources/PATH_TRANSPORT_AND_FOREST_RULE.md`, Section 1, retains actual predecessor/successor owners on the SAME vertex, explaining why independent histograms lose transport information. Section 4 similarly retains joint based cycle holonomies rather than independently quotienting phases. The present theorem shares the rigorous information-loss mechanism: local equalities must be retained on the SAME nine cell labels. Its map is explicitly the CRT intersection map above. This is an analogy of sufficient relational data, not a map from magic squares to tournaments, a metagraph wiring theorem, or a claim that tournament guards are inherited.
