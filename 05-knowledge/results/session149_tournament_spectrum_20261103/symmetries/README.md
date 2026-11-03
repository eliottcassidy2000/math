# Converse symmetry of small tournament classes

## Definitions and exact finite diagrams

A tournament has exactly one arc between every pair of vertices. Its converse reverses every arc. Converse descends to an intrinsic involution on isomorphism classes. A loop below denotes a fixed (self-converse) class; an arrow pairs two distinct converse classes. `H` counts directed Hamilton vertex-words, and `h = H/|Aut|` counts their orbits. The automorphism action on words is free.

Represent each class by the **smallest integer** mask among its relabellings. In `itertools.combinations(range(n), 2)` order, bit k equals 1 when the smaller endpoint beats the larger; bit 0 is the least significant bit. These are precisely specified codes, not an assertion of an invariant edge ordering. The numeric minimum over all bijective relabellings is invariant under isomorphism.

* n=3: `0 ↺` (H=1, h=1), `2 ↺` (H=3, h=1). These are the transitive and directed triangle classes. Thus converse alone does not pair all classes even at order 3.
* n=4: `0 ↺` (H=1, h=1), `4 ↺` (H=5, h=5), `2 ↔ 5` (H=3, h=1 for both). Mask 0 is transitive; 4 is strong; 2 is a source over a directed triangle; 5 is a directed triangle over a sink.
* n=5: `0 ↺` (H=1, h=1), `5 ↺` (3,1), `8 ↺` (9,9), `10 ↺` (9,9), `12 ↺` (11,11), `40 ↺` (15,5), `41 ↺` (13,13), `76 ↺` (15,3); `2 ↔ 11` (3,1) and `4 ↔ 9` (5,5). Distinct masks with equal displayed numbers, such as 8 and 10, need not be isomorphic; canonical masks establish their distinction.

The directed path counts are *labelled vertex-words within one fixed representative*, not numbers of all labelled tournaments. Duality reverses each word, hence H is constant on converse pairs. `class_census.json` records every class, its sorted outdegree multiset, H, automorphism order, h, and converse; the script scans all 2^(n choose 2) labelled tournaments and all permutations for n=3,4,5. Run `python check_classes.py`. The output and JSON are reproducible. These are finite exact checks, not an all-order classification of individual tournaments.

The actual converse orbit counts are respectively `(fixed, nontrivial pairs) = (2,0), (2,1), (8,2)` at n=3,4,5. The evenness theorem says the *total* number is even for all n≥3, but does not make this converse involution free.

## Prior FOLD correspondence and self-converse parameterization

This section is a self-contained restatement of the prior Session132 FOLD correspondence, not a novelty claim. [Session132, §3](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session132_fractran_tournament_recursion_20261031/README.md) and its [fold theorem](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session132_fractran_tournament_recursion_20261031/converse_folding/THEORY.md) already establish anti-involution existence and conjugacy, the class bijection and its stronger marked-groupoid version. Session132 encodes winner blocks as a pair of tournaments on one common vertex set, with a binary vertex color in odd order. The fixed-r centralizer quotient below is an equivalent presentation and its small census is a fresh computational control. This gives a bijection of unmarked isomorphism classes, not an equivalence of unmarked groupoids: folded stabilizers are centralizers and can be strictly smaller than full automorphism groups (e.g. C3).

Every automorphism group of a finite tournament has odd order: by Cauchy's theorem an even-order group contains an involution, which swaps some two vertices, impossible because it would reverse their unique arc. Let T be self-converse and g an anti-automorphism. Since g² is an ordinary automorphism, g has twice odd order; an odd power is an anti-involution. Any anti-involution has at most one fixed vertex, since an arc between two fixed vertices could not be reversed. Thus for n=2m+ε, ε in {0,1}, its cycle type is 2^m 1^ε.

Any two anti-involutions a,b of the same tournament are conjugate by an automorphism. Indeed h=ba belongs to the odd-order automorphism group, and a h a=h^-1. Choose 2k≡1 mod ord(h). Then h^k a h^-k = h^(2k)a = ha=b.

Fix one involution r on [n] with m transpositions and ε fixed points. Its centralizer is C2 wr S_m (the optional singleton is fixed). Consider tournaments on [n] reversed by r, up to this centralizer. This set is in bijection with the self-converse tournament isomorphism classes: existence of a conjugate of r proves surjectivity; for injectivity, any isomorphism between two such tournaments carries r to an anti-involution of the target, which can be conjugated back to r by a target automorphism, making the composite centralize r.

There are exactly 2^(m²+εm) labelled tournaments reversed by r: choose m internal-pair bits; between each two paired blocks choose two bits, with their reversals forced; and choose one bit for each edge orbit between the singleton and a block. Thus the above centralizer-orbit space gives an explicit exact finite model of all self-converse classes, although it does not claim an efficiently enumerated list at arbitrary order. The counts n=0..8 are `1,1,1,2,2,8,12,88,176`; they can also be obtained independently from the verified twisted-Burnside code in the source below.

## An explicitly artificial complete pairing

For each n≥3, sort all class representatives by their canonical integer masks and pair consecutive representatives, 1st with 2nd, 3rd with 4th, and so on. By the independently proved evenness theorem this is a fixed-point-free involution at every n≥3. It depends on the chosen encoding and is not an intrinsic tournament operation, a converse pairing, or an explanatory proof of evenness. Examples: n=3 `(0,2)`; n=4 `(0,2)`, `(4,5)`; n=5 `(0,2)`, `(4,5)`, `(8,9)`, `(10,11)`, `(12,40)`, `(41,76)`. `class_census.json` contains these arrays. Do not infer that this artificial pairing preserves H, scores, or strong connectivity.

Evenness is justified separately: in the local ring of rational numbers with odd denominator, the number of classes is congruent to the groupoid mass `sum_[T] 1/|Aut(T)| = 2^(n choose 2)/n!`. The valuation of this mass is at least `(n−1)(n−2)/2`, positive for n≥3. This also shows that n=0,1,2 are exceptions, each with one class.

## Sources and limitations

* [OEIS A000568](https://oeis.org/A000568): ordinary class sequence and Davis's odd-partition formula.
* John W. Moon, [Topics on Tournaments](https://mirror.csclub.uwaterloo.ca/gutenberg/4/2/8/3/42833/42833-pdf.pdf), Chapter 9, p. 116: asks the parity question. The above proof uses the odd-automorphism-group fact; do not credit this particular parity proof to Moon without further source verification.
* Repository source `session22-converse-channels/README.md` gives twisted Burnside formulas and independently checked self-converse counts. `math/01-canon/theorems/THM-4750-the-camion-function-half-speed-interpolation-and-the-camion-lift.md`, §6, supplies the groupoid parity proof.
* The code is a direct exhaustive census through n=5, not a proof of efficient general computation or an intrinsic free pairing.

Independent control: run `python verify_anti_orbits.py` for a separate adjacency-matrix enumeration of r-reversed tournaments, centralizer size and its orbit count for n=0..5. Its output is `anti_orbits.stdout`. The main class census also succeeds with `python -O`, with output identical to the normal run (`check.stdout` and `check_optimized.stdout`).
