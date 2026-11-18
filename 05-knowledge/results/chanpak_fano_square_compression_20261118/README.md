# Flags, Fano coordinates, and the two faces of a square-sum arrival

This follows the [one-spoke multiplier law](../chanpak_two_rail_pell_switches_20261118/README.md) and the [pointed deletion calculation](../chanpak_f1_pointed_rails_20261118/README.md). Three different decompositions of seven turn out to live on the same eight-state binary space, but they preserve different information.

## 1. Four charts of six, or six charts of four

A full flag of a tetrahedron is a vertex v in an edge e in a face f. It is encoded by the unique permutation (v,w,x,y) whose initial sets are those three faces. There are 24 flags. By v they form four groups of six. By e they form six groups of four (choose either endpoint and either incident face). No one of the six edges has an intrinsic rank in an unlabeled tetrahedron, because S4 acts transitively on them. A reference ordering, or a chosen flag, can impose a rank.

These flags are not the 64 labeled tournaments on four vertices: the latter make six independent direction choices. Nor does either partition of the flags by itself give 80. There is a coherent additional interpretation of 80: choose a vertex and place **balanced** binary marks on its six flags, three of each sign. Then there are 4·binom(6,3)=80 rooted records. They have relabeling orbit sizes 8,24,24,24. For a fixed root the stabilizer S3 acts regularly on the six local flags. A triple of them has stabilizer of order dividing both 3 and 6, so an orbit has size two or six. Exactly two triples (the two parity cosets) are invariant under the order-three subgroup, and the remaining 18 form three orbits of six. Marked signs on flags are not literal arcs; the 80-model is a possible extra rule, not a claim that unrestricted orientations give 80.

A basic compression of the genuine orientations is H=1+2·(number of cyclic triangular faces), since four vertices cannot pack two nontrivial odd cycles or support an odd cycle of length five. The four sorted-score classes give H=1,3,3,5 with orbit sizes 24,8,8,24. Score lists distinguish the two H=3 classes, while H erases that distinction.

## 2. One binary space, three partitions of seven

Start from a reference orientation. An edit is six bits on E(K4). Switching across a cut adds a cut vector, so there are 8 edit classes, each containing 8 vectors. A full decoder takes a spanning tree and three chord bits: switch to erase the tree bits, which determines a unique representative. The invariant four-face form records the XOR of edit bits on each triangular face. It lands in the even-weight subspace V of F2^4; its kernel is the cut space (a tree gauge or the face equations prove it), hence V is the 3-bit quotient.

The inherited [THM-4553](../../../01-canon/theorems/THM-4553-the-point-stripping-ladder-psl27-borel-torus-and-the-collatz-frobenius.md) uses {1,2,4} as quadratic residues modulo seven and an actual orbit under doubling. The earlier tetrahedral Fano construction separately uses binary reduced incidence coordinates on the six edges and adjoins an omega to fill the Fano plane. Here the **quotient** V gives a complementary flag reading. Its seven nonzero vectors are the Fano point set; every 2D subspace gives a line. Their partitions are:

- **Dyadic flag:** choose a basis (x,y,z), and read a nonzero vector by its last nonzero coordinate. The layers have sizes 1,2,4. Such a filtration is not invariant under arbitrary relabeling.
- **Tetrahedral symmetry:** permute the four faces. The all-ones vector is fixed, and the six weight-two vectors form one orbit. Thus the unmarked partition is 1+6.
- **Frobenius:** fix a face f0 and cycle the other three. The basis (e_f0+e_f1,e_f0+e_f2,e_f0+e_f3) is rotated. Identify it with (theta,theta²,theta⁴), where theta is the root of X³+X²+1 in F8. These elements are independent and sum to 1. Indeed using the polynomial basis they are theta, theta², 1+theta+theta². The rotation really is squaring, with orbits 1+3+3 on the nonzero vectors. In multiplicative exponent notation one rotating triple is precisely 1→2→4→1 mod 7. The explicit conjugacy does **not** identify that three-element orbit with three size layers.

A reference orientation matters. V describes additive **differences** of orientations. The orientation switching classes themselves are an affine space without a relabeling-invariant zero: the computed eight fibers have two with all H=3 and six with four H=1 and four H=5. In fact these groups are the S4 orbits: relabeling transitive orientations meets all six mixed fibers, and transposing two vertices of the cyclic face of an H=3 example reverses its face parity, reaching the other all-3 fiber. Thus the Hamilton count is not an invariant of the switching compression. For a fixed transitive spine the one-spoke integer ratios 1+2^b have just the anchor and one dyadic size; 7=1+2+4 uses three layers, not an extra choice of b. The quotient organizes the difference but is not a new ×7 construction.

## 3. The first 32 square-sum graphs and the shared two-port test

Write G_n for the graph on 1,…,n with x~y iff x+y is a square. The independent vertex-order replay reproduces the first 32 counts. In the convention used here a path and its reverse are identified, but a closable path is still a path with a particular cut (rather than counting its entire cycle just once). OEIS has [the cycle counts](https://oeis.org/A071984) and [the precisely matching path convention](https://oeis.org/A071983).

At n=32 the newborn has only neighbors 4 and 17. A spanning cycle exists iff G31 has a spanning path with **both** exposed endpoints 4 and 17; extend it through 4–32–17. The G31 incidence table is:

| G31 path's endpoints | Number |
|---|---:|
| neither port | 291 |
| 4 only | 34 |
| 17 only | 23 |
| both | 1 |

Thus 57 paths offer exactly one newborn endpoint extension, one offers two, and there is exactly one cycle closure. In G32, 59 paths have the newborn as an endpoint and 333 have it internal, making 392 in all. The unique spanning cycle gives 32 distinct open paths by cutting its 32 edges, only two of which have the newborn as an endpoint. This recovers the first cycle at 32 and distinguishes 32 from the full path count. A port-marked graph at birth is the correct object, not just the scalar number of paths.

For comparison, the upper two-rail tournament has 16 two-run templates with one surviving nonbase deletion face and 5 with both. Each has two width-two clone allocations. In particular its 32 one-face paths are `16×2`; the lower family has the two-run signature `(2,1)`, and a sevenfold deletion-compatible copy would have `(14,7)` rather than `(16,5)`. The two 32s do **not** label the same kind of objects: one comes from cuts of a closed cycle, the other from allocating two clones to templates with exactly one live face. The real transfer is from the square closure criterion: to get a stronger tournament repair, turn two specific one-live templates into two-live templates without changing the total, then match the actual exposed endpoint words. Merely increasing the one-face count cannot supply a missing simultaneous face.

Replay: `python tetrahedron.py | diff -u tetrahedron.out -` and `python squares.py | diff -u squares.out -`.

A second exact comparison actually deletes the last two vertices separately from each square-sum path, rather than asking about G31 endpoint exposures. Through orders 2–32 no graph has exactly 32 **oriented** paths with exactly one surviving deletion. The G32 histogram (neither, one, both) is (644,134,6). An interior deletion survives exactly when the two remaining neighbors admit a shortcut; at 31 and 32 the respective neighbor pairs (5,18) and (4,17) admit none. Thus both endpoint incidence and shortcut rules are needed when carrying this comparison elsewhere.

Finally, [a separate Paley construction](SEVEN_SHEETS.md) really does realize seven free sheets commuting with physical deletion by replacing the core. It uses the actual {1,2,4} step orbit and gives an exact unbounded 7H rail, though not one containing prime lower counts.
