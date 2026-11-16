# The converse census: 33 at order 6, and the multiplicities fixed + 2·pairs

Session klein-2026-11-16-monoid. Script `scripts/converse_census.py` (every tournament class of order `n`, its count
`H`, and whether it is isomorphic to its arc reversal — brute force over permutations); receipts `converse_4..7.out`.

**33 at order 6.** Not attained at order 5 (maximum 15); at order 6 attained by exactly 5 classes: 1 self-converse,
2 converse pairs (pre-registered P-1 **confirmed**; S224 and `iso_class_fractal.out` on file say the same). The same
`1 + 2·2` shape holds for 9 and 37 at order 6, and 33 at order 7 has `1 + 2·6`.

**The multiplicity table.** Since `H(T^op) = H(T)`, each count is attained by `s + 2p` classes (`s` self-converse,
`p` pairs). Orders 4–7:

| order | classes | self-converse | counts attained only by pairs | counts attained only by self-converse classes |
|---|---|---|---|---|
| 4 | 4 | 2 | 3 | 1, 5 |
| 5 | 12 | 8 | 5 | 1, 9, 11, 13, 15 |
| 6 | 56 | 12 | 3, 11, 13, 15, 19, 23, 25, 27, 31, 43 | 1, 17, 41, 45 |
| 7 | 456 | 88 | 5, 17, 19, 23, 29, 35, 37, 39, 43, 47, 53, 55, 57, 59, 73, 89, 97, 101, 103, 109, 113, 121, 131 | 1, 61, 85, 87, 125, 127, 141, 145, 147, 155, 157, 159, 171, 175, 189 |

Order 6 in full (`H: classes, self-converse, pairs`): 1: 1,1,0; 3: 4,0,2; 5: 3,1,1; 9: 5,1,2; 11: 2,0,1; 13: 2,0,1;
15: 6,0,3; 17: 2,2,0; 19: 2,0,1; 23: 4,0,2; 25: 2,0,1; 27: 2,0,1; 29: 4,2,1; 31: 2,0,1; 33: 5,1,2; 37: 5,1,2; 41: 1,1,0;
43: 2,0,1; 45: 2,2,0. Order 7 in `converse_7.out` (e.g. 89: 14,0,7; 91: 15,1,7; 135: 8,4,2; 159: 4,4,0).

Pre-registered P-2 **confirmed** (at order 6 some counts are attained only by pairs, some only by self-converse
classes). Reading: the maximal counts at each order (15 at order 5; 41, 45 at order 6; 125–189 at order 7) are carried
by self-converse classes only, the small and middle counts mostly by pairs; the previous packet's question Q1 (is the
number of classes of count `H` congruent to the number of self-converse classes mod 2) is trivially yes, and the table
is the nontrivial content.

**The Platonic reading, corrected as the owner says.** The three flat tilings curve into the five solids as
`{3,6} ↦ {3,3}, {3,4}, {3,5}` (the triangular tiling into the three solids with triangular faces), `{4,4} ↦ {4,3}`,
`{6,3} ↦ {5,3}`; the tetrahedron is self-dual, cube and octahedron keep 12 edges and swap 8 vertices with 6 faces,
icosahedron and dodecahedron keep 30 edges and swap 12 with 20. The count `1 + 2·2` of 33 at order 6 has the same
shape: one self-converse class, two converse pairs, each pair keeping the count and swapping in- and out-degrees. The
analogy is of shape, not of mechanism (ANALOGY).
