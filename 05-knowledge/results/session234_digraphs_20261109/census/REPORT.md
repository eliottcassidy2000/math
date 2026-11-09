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

## Structure and all-order consequence using the strong-tournament floor

Every semicomplete digraph is a uniquely ordered sum of its strong components: arcs between components all go forward. Therefore H(D)=product H(C_i). Each component must occur as one contiguous block in every Hamiltonian path, and arbitrary component Hamiltonian paths concatenate uniquely.

Using the proved strong-tournament lower bound excluding H<=7 for order>=5 and excluding H<=21 for order>=7, the scans above yield a complete all-order classification. A strong semicomplete digraph of order>=3 has a Hamiltonian cycle and hence a spanning strong tournament (orient unused reciprocal pairs arbitrarily). Its count is at least that of this tournament. Components of order2 have count2 and cannot appear in odd-count families.

The cited floor and independently verified census give:
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

Python n<=5 completes in a few seconds. C++ n6 including full target canonicalization and direct target verification took about 12.4 seconds in this execution environment. No external libraries required.

Outputs:
- results.json: complete n<=5 histogram, target representatives, SCC signatures, orbit sizes, reciprocal pairs, all witness paths and tournament-resolution spectra.
- results6.csv: complete n6 count histogram.
- classes6.csv: all n6 target classes, canonical bits, representative ternary code, strong status, reciprocal count and orbit size.
- classes6.json: same classes with explicit matrices, all paths and all tournament-resolution matrices/counts.

Ternary code convention: pairs (0,1),(0,2),...,(n-2,n-1), in that order, consume least-significant base3 digits. Digit0 means i->j only; 1 means j->i only; 2 means both. Every code is thus unique and exhaustive.

## Strong H3 check and parity warning

classify3.py independently enumerates every n<=4 semicomplete digraph. The only strong H3 cores are n3 code3 (directed triangle, orbit2, automorphism order3) and code5 (triangle with one doubled edge, orbit6, asymmetric). There are no n4 strong H3 cores. A strong tournament lower bound of at least5 for order>=4 makes this classification global. strong_h3.json records the cores.

Thus one doubled pair does not force an even path total: the smallest odd example has n3 and H3. At n3, exactly-one-reciprocal-pair labelled counts are {H2:6,H3:6}. At n4 they are {H2:36,H3:48,H5:24,H6:72,H8:12}. Adding an arc can create zero new Hamiltonian paths.

## Loop decorations and complete generator

`generator.py` is genuinely standalone (all 37 strong-core matrices embedded, Python standard library only). It accepts `--target 7|21 --n N [--with-loops] [--limit K]` and emits one adjacency-matrix representative per isomorphism class as JSON Lines. It is lazy; `--limit` avoids enumerating the full requested size. `--validate` runs the core-decoration census and full loopless generation through n8 and saves `generator_validation.json`.

For every strong core we enumerated its actual automorphism group by testing all vertex permutations, partitioned all 2^n loop masks into group orbits, and checked Burnside's lemma independently by counting fixed masks. Every decorated core's Hamiltonian-path count was rechecked. Results:

| Target | Core order | Loopless cores | Loop-decorated cores |
|---|---:|---:|---:|
|3|3|2|12|
|7|4|3|48|
|21|5|16|480|
|21|6|16|1024|

The triangle contributes four decorated cores; its asymmetric doubled-edge counterpart contributes eight. The three H7 cores contribute 16 each. For H21 at n5, fourteen asymmetric cores contribute 32 each; each of two order-3-automorphism cores contributes 16. Their two nonidentity automorphisms have cycle type (1,1,3), so Burnside gives (32+8+8)/3=16. All n6 H21 cores are asymmetric and contribute 64 each.

Write x for vertex order. The ordinary generating functions for isomorphism types are:

- Loopless H7: 3 x^4 / (1-x)^2.
- Loopless H21: (16 x^5 + 16 x^6)/(1-x)^2 + 12 x^7/(1-x)^3.
- Loops allowed H7: 48 x^4 / (1-2x)^2.
- Loops allowed H21: (480 x^5 + 1024 x^6)/(1-2x)^2 + 1152 x^7/(1-2x)^3.

The exponential generating functions for labelled loopless counts are:

- H7: 3 x^4 / (1-x)^2.
- H21: ((44/3) x^5 + 16 x^6)/(1-x)^2 + 8 x^7/(1-x)^3.

Allowing loops multiplies labelled counts at order n by 2^n, equivalently substitutes 2x in those exponential generating functions. It does not simply multiply isomorphism counts by 2^n.

Each singleton SCC has two loop states, giving the decorated gap factor (1-2x)^-1. Product-family decorated coefficient is 2 orders * 48 H7 cores * 12 H3 cores = 1152. No component permutations are permitted because the SCC quotient is a directed total order.

The validator constructs every loopless generated matrix through n8, independently checks its Hamiltonian-path count with subset DP, extracts its actual SCCs, canonically labels each SCC by all vertex permutations, and confirms no two ordered canonical-SCC signatures coincide. It also checks the class count and sums n!/|Aut| against the labelled formulas. New n7/n8 values:

| Target | n | Isomorphism classes | Labelled graphs |
|---|---:|---:|---:|
|7|7|12|60480|
|7|8|15|604800|
|21|7|92|423360|
|21|8|148|5268480|

These n7/n8 are exhaustive outputs of the proved-core generator, not direct scans of all 3^(n choose 2) input graphs. Their global completeness depends on the strong-tournament floor as explained above. The validation itself is fully computational and collision-free. Default validation deliberately does not enumerate every loop-decorated n8 graph; it exhaustively checks decorated cores and loopless families, while decorated generation is available explicitly.

## All loopless digraphs through order five

`general_scan.cpp` enumerates all 4^(n choose 2) loopless digraphs for n<=5, permitting absent unordered pairs. Every graph is counted by subset DP; every H7/H21 target is independently checked by direct permutation enumeration. `general_census.csv` gives the complete target census jointly stratified by n, H, number of absent pairs, number of reciprocal pairs, and strongness. `general_verification.csv` records coverage; `general_summary.json` gives marginal totals. Reproduce with `g++ -O3 general_scan.cpp -o general_scan && ./general_scan`.

| n | All input digraphs | H7 targets | H21 targets |
|---|---:|---:|---:|
|1|1|0|0|
|2|4|0|0|
|3|64|0|0|
|4|4096|120|0|
|5|1048576|49440|7760|

All 120 order-four H7 targets are strong. Of these, 72 are semicomplete (two reciprocal pairs), while 48 have one absent pair and three reciprocal pairs.

At n5 the H7 targets have absent-pair distribution 0:720, 1:7680, 2:21240, 3:17640, 4:2160. Their reciprocal-pair distribution is 0:1440, 1:8040, 2:15960, 3:15240, 4:8280, 5:480. There are 7200 nonstrong and 42240 strong targets.

At n5 the H21 targets have absent-pair distribution 0:1760, 1:4440, 2:1320, 3:240. Their reciprocal-pair distribution is 2:720, 3:1280, 4:3960, 5:1560, 6:240. All 7760 are strong.

Allowing arbitrary diagonal entries multiplies each labelled count at n by 2^n, with the same absent/reciprocal-pair statistics and strongness.

### Global minimal orders, not restricted to semicomplete graphs

H7 requires n>=4 since n!<7 for n<=3, and an explicit n4 example is above. H21 requires n>=5: the complete bidirected n4 graph has 24 Hamiltonian paths; any other loopless n4 digraph omits at least one directed arc, excluding all (4-1)!=6 vertex permutations using that arc consecutively, so it has at most18 paths. Smaller orders have at most6. A n5 H21 example above attains the minimum. Loops change neither conclusion.

### Exact universal path-support construction

This characterizes arbitrary digraphs, including missing pairs, without claiming a finite strong-core classification outside the semicomplete class.

Fix a vertex set V. For each vertex permutation Q, let E(Q) be its set of consecutive directed arcs. Choose a proposed path family P of exactly h permutations and set U=union of E(P) over P in the family. Retain the family only if the graph U has exactly these h Hamiltonian paths (equivalently, every permutation outside P has at least one arc outside U).

On the remaining directed arcs R, form the forbidden hyperedges F_Q=E(Q) minus U, one for each Q outside P. They are nonempty. Choose any X subset of R containing no entire forbidden hyperedge. Then the graph with arcs U union X has exactly the h paths in P. Conversely, every loopless graph with h paths occurs uniquely through its actual path family P, its support U, and X=E minus U.

Proof: all paths in P use arcs of U and remain present. An outside permutation Q becomes a path exactly when its missing-from-U arc set F_Q is contained in X. Thus hyperedge avoidance is necessary and sufficient. For the converse, U is a subgraph of the original graph containing all its paths, so its path set is exactly P; any extra arcs must obey the avoidance condition. Diagonal entries can then be freely decorated.

The general scan verifies this construction for all 120 order-four H7 graphs. For each it extracts its actual path family, constructs U, verifies H(U)=7, and checks every forbidden hyperedge against the observed extra arcs. It additionally exhausts **every possible extension X above each extracted U**, independently recounts Hamiltonian paths by permutations, and confirms that hyperedge avoidance is equivalent to retaining exactly seven paths. This is stronger than checking just the original target graph.
