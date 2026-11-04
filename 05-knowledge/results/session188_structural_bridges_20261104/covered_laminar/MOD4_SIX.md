# Complete named responses modulo four distinguish every six-name tournament

**FINITE-EXACT.** The 32,768 labelled tournaments on six fixed old names have 32,768 distinct complete named two-port response vectors after coefficientwise reduction modulo four. No quotient by relabeling is taken. This is a finite modulo-four property, not an all-order reconstruction theorem. It is not presented as a new maximum-order record for exact-integer reconstruction.

## Observer and retained data

The two physical ports are x=6 and y=7, with x->y. A bit a_i=1 means i->x, and c_i=1 means i->y. F_T(a,c) is the literal Hamilton-path count of that eight-name tournament. Every Boolean coefficient of bidegree at most (2,2) is retained, including all cubics, quartics and named zeros. There are (1+6+15)^2=484 coordinates.

The old tournament code has 15 bits in lexicographic pair order (i,j), 0<=i<j<6; a bit 1 means i->j. The coordinate mask is m=a|(c<<6). Coordinates occur in increasing m among masks with popcount(a)<=2 and popcount(c)<=2.

`mod4_six_vectors.bin` contains exactly 484*32768=15,859,712 bytes. It is coefficient-major: byte offset b*32768+code is the canonical residue 0,1,2,3 of coordinate b for that physical tournament code. It has no header, padding, native integers, endian dependence or omitted zeros. Lossless compression does not change this schema. Entire 484-byte vectors are compared exactly using ordered containers; hashes are used only for file integrity, never as an equality oracle.

## Producer algorithm

`mod4_six.cpp` enumerates all 8! Hamilton words on the fixed old and port names. It rejects words using y->x. Every remaining old-old arc supplies either an old-edge bit or its complement; every old-port contact supplies exactly a_i,1-a_i,c_i or1-c_i. Expand those Boolean factors and accumulate every signed monomial modulo four. A Boolean subset zeta transform in the 15 old-edge variables evaluates the entire coefficient bank on every tournament.

There are 35,280 surviving words and 451,249 signed monomial contributions before the old-edge zeta transform. Coefficient accumulation, transform and vector comparison use exact integer/residue arithmetic. This is not interpolation from selected named moments.

## Separate validation

`validate_six.cpp` uses literal integer Hamilton-path subset DP rather than the producer's word expansion or old-edge zeta transform. It checks:

* one deterministic, arbitrary physical query for each of all 32,768 tournaments;
* all 484 lower-closed interpolation queries on each of 68 specified tournaments, then every coefficient by integer Boolean finite differences (32,912 coefficient checks);
* all 4,096 physical queries for one further fixed tournament;
* exact pairwise uniqueness by complete vectors in a separate ordered set.

The total is 69,776 literal integer DP queries. `check_six.py` independently implements the literal integer DP in Python and checks another 272 retained physical queries. Both C++ programs reproduce identical receipts with `-O2` and `-O0`; Python normal and `-O` receipts agree. The binary bank's integrity hash is pinned by the Python receipts and the packet manifest.

## A higher-degree hostile, not a truncation license

`degree_six.py` groups the same retained vectors by degree only as a finite diagnostic. Degree at most two modulo four leaves 240 collision pairs. Of these, 120 are original-module XOR pairs from fully covered laminar families, and 120 lie outside the original-module incidence span. Degree at most three separates this finite six-name universe. The complete degree-four data remain retained; no truncation claim is made at other sizes or for arbitrary ambient comparisons.

For the first pair in code order,

    T code408:   (48,49,3,7,12,28),
    U code15462: (12,13,48,52,3,19),

simultaneously reverse original modules with masks 3,12,48,63. This is a labelled instance of the inherited whole-plus-three-pairs C3[TT2,TT2,TT2] control, not a newly discovered structural family. Its modulo-four difference has 48 cubic and 48 quartic nonzero coordinates and no lower-degree terms. The exact physical query (a,c)=(6,1) gives

    F_T(6,1)=455,       F_U(6,1)=477,

so the difference is 2 modulo four. The complete nonzero modular difference and all 240 low-degree collision pairs are retained in `degree_diagnostic.json`. No exact full-response collision is claimed.

## Reproduction

    g++ -std=c++17 -O2 -Wall -Wextra -Werror mod4_six.cpp -o mod4_six
    ./mod4_six > mod4_six.json
    g++ -std=c++17 -O2 -Wall -Wextra -Werror validate_six.cpp -o validate_six
    ./validate_six > validation.json
    python check_six.py > six_python.json
    python -O check_six.py > six_python.optimized.json
    python degree_six.py > degree_diagnostic.json
    python -O degree_six.py > degree_diagnostic.optimized.json

The `-O0` C++ rebuilds have their own matching receipts. The exact all-order result in this packet is separately stated in `ORDERED_BLOCKS.md`; finite injectivity does not prove its remaining strongly connected laminar frontier.
