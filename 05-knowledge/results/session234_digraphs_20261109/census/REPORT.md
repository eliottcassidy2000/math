# Exact semicomplete Hamiltonian-path scan

Convention: labelled directed Hamiltonian paths are vertex permutations satisfying each consecutive directed arc. Reverse orders count separately. Loops are zero; every unordered vertex pair has one or both arcs.

## Verified exhaustive results

| n | Labelled semicomplete digraphs | H=7 | H=21 | H=7 isomorphism classes | H=21 isomorphism classes |
|---|---:|---:|---:|---:|---:|
|1|1|0|0|0|0|
|2|3|0|0|0|0|
|3|27|0|0|0|0|
|4|729|72|0|3|0|
|5|59049|720|1760|6|16|
|6|14348907|6480|32640|9|48|

At n<=5, subset DP and independent permutation enumeration agree for **every graph**. At n=6, subset DP enumerates every graph and independent permutation enumeration agrees for **every graph with count 7 or 21**. Canonical labels minimize the row-major adjacency bitstring over every vertex permutation. This is exhaustive, not a sampling or heuristic isomorphism test.

All three n4 H7 classes are strong, asymmetric, and have two reciprocal pairs. Every n5/n6 H7 class is obtained by placing the remaining vertices as singleton components before and after one of these cores.

All sixteen n5 H21 classes are strong. Six have two reciprocal pairs; ten have three. Fourteen are asymmetric, two have automorphism group order 3 (each has 40 labelings).

There are sixteen strong n6 H21 classes, all asymmetric. Two have one reciprocal pair, eight have two, six have three. The other 32 n6 H21 classes are obtained by putting a singleton before or after one of the sixteen n5 cores.

## Small explicit counterexamples

H7, n4, code24, matrix rows:

0111
0011
1001
1000

H21, n5, code801, matrix rows:

01111
00110
00011
10001
11000

H21 with just **one** reciprocal pair, n6, code99, matrix rows:

011110
001111
000111
100011
000001
100000

The single reciprocal pair is {0,3}. Its two tournament resolutions have 17 and 19 paths. Adding the absent reverse arc to the 17-path tournament produces four additional paths; adding the other arc to the 19-path tournament produces two additional paths. Thus reciprocal arcs destroy the tournament gap at 21 without a missing pair or loop.

The second one-reciprocal-pair class is code129, rows 010100 / 001111 / 100111 / 100011 / 100001 / 100000. Its resolution spectrum is also {17,19}.

H7 requires at least four vertices; H21 requires at least five. These are certified by the smaller-order exhaustive scans. Up through order six H7 needs at least two reciprocal pairs, while H21 needs one (first possible at order six).

## Structure and all-order consequence conditional on the strong-tournament floor

Every semicomplete digraph is a uniquely ordered sum of its strong components: arcs between components all go forward. Therefore H(D)=product H(C_i). Each component must occur as one contiguous block in every Hamiltonian path, and arbitrary component Hamiltonian paths concatenate uniquely.

If one uses a proved strong-tournament lower bound excluding H<=7 for order>=5 and excluding H<=21 for order>=7, the scans above yield a complete all-order classification. A strong semicomplete digraph of order>=3 has a Hamiltonian cycle and hence a spanning strong tournament (orient unused reciprocal pairs arbitrarily). Its count is at least that of this tournament. Components of order2 have count2 and cannot appear in odd-count families.

Under that floor:
- H7: one of the three n4 cores, with arbitrary singleton chains before and after.
- H21: one of sixteen n5 or sixteen n6 strong cores with singleton chains before and after; OR one of the three n4 H7 cores and one of the two strong H3 cores (a directed triangle, or a directed triangle with one edge doubled), in either component order, with arbitrary singleton chains before, between, and after.

This gives n>=4 H7 labelled count 3(n-3)n! and isomorphism-class count 3(n-3).

For H21, sum contributions only where their core order fits:
- n5 cores: (44/3)(n-4)n! labelled; 16(n-4) classes, n>=5.
- n6 cores: 16(n-5)n! labelled; 16(n-5) classes, n>=6.
- n4+H3 cores: 4(n-5)(n-6)n! labelled; 6(n-5)(n-6) classes, n>=7.

These all-order statements depend on the cited floor, which is not proved by this scan. Finite exhaustive results above do not depend on it.

## Reproduction

Run in this directory:

python scan.py
python analyze.py
g++ -O3 scan6.cpp -o scan6
./scan6
python enrich6.py

Python n<=5 completes in a few seconds. C++ n6 including full target canonicalization and direct target verification took about 12.4 seconds on this worker. No external libraries required.

Outputs:
- results.json: complete n<=5 histogram, target representatives, SCC signatures, orbit sizes, reciprocal pairs, all witness paths and tournament-resolution spectra.
- results6.csv: complete n6 count histogram.
- classes6.csv: all n6 target classes, canonical bits, representative ternary code, strong status, reciprocal count and orbit size.
- classes6.json: same classes with explicit matrices, all paths and all tournament-resolution matrices/counts.

Ternary code convention: pairs (0,1),(0,2),...,(n-2,n-1), in that order, consume least-significant base3 digits. Digit0 means i->j only; 1 means j->i only; 2 means both. Every code is thus unique and exhaustive.

## Strong H3 check and parity warning

classify3.py independently enumerates every n<=4 semicomplete digraph. The only strong H3 cores are n3 code3 (directed triangle, orbit2, automorphism order3) and code5 (triangle with one doubled edge, orbit6, asymmetric). There are no n4 strong H3 cores. A strong tournament lower bound of at least5 for order>=4 makes this classification global. strong_h3.json records the cores.

Thus one doubled pair does not force an even path total: the smallest odd example has n3 and H3. At n3, exactly-one-reciprocal-pair labelled counts are {H2:6,H3:6}. At n4 they are {H2:36,H3:48,H5:24,H6:72,H8:12}. Adding an arc can create zero new Hamiltonian paths.
