# Session68 independent proof audit

**Status: ACCEPTED final classification, literal-line boundary, arithmetic report and bounded controls.** Exact reviewed inputs are frozen in `FINAL_INPUT_MANIFEST.json`. No global prime-tournament or metagraph classification is certified.

## Evidence and inheritance

The repository mathematics was read at base `fc3768ac5fcd87022fde9e6d7e3f06c0e5bc22e8`. Concurrent `origin/main` commit `d761e4cf0ff92264b043586a6160b274dae0b29c` was then inspected read-only. Its `05-knowledge/results/tournament_composition_classification_20261026/CLASSIFICATION.md` already supplies the decorated-tree classifier, coherent split automorphism extension, and the four-type 128-realization order-twelve fibre. Its arithmetic result supplies the Sylow/digit bounds. These are inherited/concurrent mechanisms, not Session68 novelty. The new restricted all-depth polynomial and its complete exponent support are a further calculation on that mechanism. Subsequently `c9a8823436e` was inspected: `05-knowledge/results/cyclictriple_composition_20261026/README.md` already gives the complete carousel fibres, divisor/Möbius strata, automorphism and labelled counts, and converse pairing (reviewed source SHA256 `a0e6faed7d6e454699b3eaa0a5da2bbc1d613288665eb4a205d8ee69f04caa4e`). The final arithmetic note explicitly credits that overlap; the normal-basis/Frobenius and relative-phase constructions remain its additional connections. The final sample disposition also avoids claiming a historical summand/pinch error without recertifying it.

The imported realization facts were checked in the primary source: [Boussaïri–Chergui–Ille–Zaidi, 2018](https://arxiv.org/html/1805.04607v1), Theorem 10 (common strong modules), Theorem 24 (prime realizations up to converse), and Proposition 51/Corollary 57 (node choices and labelled count). These remain CITED inputs, not newly proved here. THM-3121 `path-cover-walk-content-substitution-kernel` already contains the full substitution kernel; THM-3134 `tournament-endpoint-jet-and-c3-newton-profile-transform` already contains cyclic-wreath divisibility. THM-785 `cyclic-triangle-flow-on-blue-black-metagraph` and THM-646 `line-metagraph-score-complement-skeleton` already contain the native flux law.

## 1. Decorated classification and group audit

The reduced canonical strong-module tree is lossless: ordered children at a linear node; oriented prime quotient, child ports, and complete ordinary child types at a prime node. Isomorphisms preserve the intrinsic strong modules, and the converse map assembles by induction. Adjacent linear nodes must be flattened. A nonstrong replacement at a linear node can invalidate the displayed tree.

The automorphism extension splits after choosing coherent models of equal-type children. Its quotient is the stabilizer of **ordinary child types**, not sizes or cyclic-triple fields. A coordinated section is not canonical; the proof does not establish that every abstract complement is conjugate. A noncanonical presentation computes only its partition-preserving subgroup. The documented order-nine 81/27 hostile correctly illustrates this distinction.

For terminal A/Aop colors on a depth-d C3 tree, put L=3^d and I=(L-1)/2. The exact automorphism order is `3^(L+e)` and the ordinary-class polynomial obeys

    P_0(z)=2,
    P_(d+1)(z)=(P_d(z)^3-P_d(z^3))/3+z P_d(z^3).

The diagonal correction is correct: only three identical child **types**, not merely equal exponents, admit the top rotation. The full-support induction is valid because every realized exponent has at least two types: color complement has no fixed ordinary class on an odd number of leaves. Every integer in 0..3I is obtained by a nonconstant triple, and a constant maximal triple gives 3I+1. Endpoint coefficients remain exactly two. The extremal groups, elementary abelian C3^L and the iterated cyclic wreath group W_(d+1), and their exponent/order claims are valid, including d=0. The order ratio `3^((n-4)/8)` and mass `3^I P_d(1/3)=2^L` follow exactly.

No global prime list is obtained. The word group is an iterated wreath group, not ordinary rotation of a word of length 3^d.

## 2. Full observer-fibre extension: accepted with the precise gauge

The terminal-color palette is not the whole **labelled** realization set. It does, however, meet every **ordinary isomorphism class** of the same hypergraph. To prove this, use the common strong-module tree. At each internal C3 node, a reflected permutation of its three child slots identifies either quotient orientation with the chosen standard orientation. Then normalize the three child realizations recursively. At a four-vertex terminal block, the two triangle orientations are isomorphic and the source/sink order leaves exactly A or Aop. This is an induction on actual child types; it never assumes that a quotient reflection is an automorphism of the already decorated tournament.

The hypergraph's internal three-port nodes permit S3, and each bottom hyperedge permits S3 on its three vertices, so

    |Aut(hypergraph)|=6^(L+I).

BCIZ counts L+I prime nodes and L binary linear nodes, hence

    number of labelled realizations=2^(2L+I).

A tournament type of exponent e has orbit size `2^(L+I)3^(I-e)` in that labelled fibre. Summing these orbit sizes with P_d gives the BCIZ count exactly. This is a completeness checksum, not an independent proof replacing the node normalization.

For a carousel R_m, the corresponding proof yields all binary necklaces as ordinary fibre types, hypergraph group `(S3)^m semidirect D_(2m)` where D_(2m) has order 2m, labelled realization count `2*4^m`, and period-d type orbit size `2^(m+1)d`. The outer reflection normalizes the reversed quotient; it does not identify two palette words unless an orientation-preserving rotation already does so.

## 3. Path-cover and metagraph audit

Module reversal preserves every unrestricted path-cover coefficient: split each directed path at its maximal runs in the reversed module and reverse each run. External incidence is unchanged by homogeneity. This is an involution; it need not preserve fixed vertex endpoints. Ordered spanning path covers have a free automorphism action; unordered ones need the factorial. Thus `t_c=H/|Aut|` applies to the ordinary unmerged marked-path fibre, not automatically to converse-merged degree.

The native complement partner retains path edges and reverses all other arcs. For n>=2 its exact laws are

    s_U(v)=n-1-s_T(v)+[v=first]-[v=last],
    c3(U)-c3(T)=s_T(first)-s_T(last)-1.

The separate-endpoint derivation does not apply at n=1; the script correctly excludes that case only for the second formula. The all-A regular-quotient construction has only odd scores, and the all-Aop construction only even scores. Therefore every literal partner changes the cyclic-triple count by an odd number. This proves nonclosure of each fixed-C3 family under native lines. It does not establish the converse statement that every mixed decoration exits. Zero triangle-count flux is not equality of hypergraphs.

Both endpoint histograms are exactly correct. Dividing each flux bucket by its actual automorphism order is justified because the action on complete words is free and preserves endpoint scores. The resulting 5715/17145 totals are ordinary tiling-port counts, not already degrees in a converse-merged metagraph.

## 4. Independent finite controls

`independent_controls.py` imports no producer code and reads no repository data. Explicit exceptions remain active with `python -O`. Normal and optimized outputs agree byte-for-byte.

- All 128 node-choice realizations at order twelve: an explicit 1296-element hypergraph action has four disjoint tournament orbits of sizes 16,48,48,16, with stabilizer orders 81,27,27,81. This is a BCIZ-node-choice universe, not all 2^66 tournaments.
- All 512 depth-two binary colorings: orbit BFS from the four local rotation generators, without the polynomial recurrence or an all-group stabilizer sum, produces 24 classes and exponent coefficients `[2,10,8,2,2]`.
- Direct Hamilton-word depth-first enumeration of both order-twelve witnesses, independent of the producer's endpoint Held–Karp calculation, reproduces every flux bucket and H=462915.
- All 385 unordered path covers of C3[A,1,1], generated by permutations and cuts, are bijected to the module-reversed tournament and back. The full profile is `[27,117,150,75,15,1]`.

Producer `check_classification.py` replay matches its stored receipt. Root `check_line_boundary.py` normal and optimized replays match stored receipts, including all 7887 Hamilton words over all labelled tournaments of orders 1..5. Root `check_field_bridge.py` normal and optimized replays also match stored receipts. Replay alone is not independent mathematics; the separate controls and proofs above are the independent path.

## 5. Arithmetic proof acceptance boundaries

For odd m=2k+1>=3, carousel primeness is elementary: any proper module lies in the k consecutive out-neighbors or in-neighbors of an outside vertex. In a rotated interval, two members a<b are separated by the outside vertex a+k+1. For an arc of step r in 1..k, the number of common out-neighbors is k-r, so the intrinsic maximizer is the cyclic successor. Hence the full quotient automorphism group is C_m. Canonical substitution then proves necklace classification and `Aut(T)=C3^m semidirect C_(m/d)`, with the actual block-permutation action.

For a chosen normal basis `(alpha,alpha^2,...,alpha^(2^(m-1)))` of F_(2^m), the linear coordinate map intertwines cyclic rotation with Frobenius. Least period is the minimal-polynomial degree d, and the stabilizer is Gal(F_(2^m)/F_2(x)) of order m/d. The basis is indispensable; this is not a multiplicative or tournament-composition isomorphism. Reflection J is generally only F2-linear. Complemented reflection has no fixed necklace for odd m because it changes weight w to m-w.

The depth-two wreath action is not this cyclic Galois action. The root checker explicitly uses the odometer order `(0,3,6,1,4,7,2,5,8)`, verifies its 9-cycle belongs to W2, and finds two field-degree-nine colors whose W2 stabilizers have orders 9 and 1. This avoids the false assertion that ordinary index-order rotation lies in the displayed block action.

For two fixed necklaces of periods d,e dividing m, their product has d*e ordered phase pairs, every diagonal rotation orbit has length lcm(d,e), and there are gcd(d,e) relative-phase orbits. Separate necklace classes lose this residue coordinate. The normal-basis and CRT claims are scoped precisely to these cyclic carriers.

## Acceptance / remaining actions

Core and arithmetic proofs and controls have no unresolved mathematical defect. Final producer notes include concurrent-prior credit and the labelled/ordinary full-fibre distinction. Classification Sections 4.1 and 5.1 were reread after freeze: the full-fibre normalization and deficit `delta_3=2I-e` are correct; the latter explicitly applies the inherited carry theorem. Both updated classification modes replay identically, with 29,077 checks. The final arithmetic report (SHA256 `d5883d398c362ce0e406a6055afc655e0a01e566546ee00074892c31f0f4c11f`) was reviewed in full; both modes of its checker replay identically. Its 24 original stratified draws and every selected blob hash reproduce at the pinned tree. Historical sampled claims not used as dependencies were not independently recertified. The normal-basis citation was checked against Conrad, Theorem 3.7, and the carousel definition against Belkhechine–Boudabbous, Section 2. The four-color palette really has distinct ordinary types because its singleton fixes the oriented quotient, and the finite-field tensor-product statement is the standard algebra-isomorphism consequence of the stated Frobenius orbit factorization, not a canonical equality of concrete carriers.

No computation through all order-36 permutations, all 2^27 terminal words, or all tournaments was run or required.
