# Six-vertex odd-cycle packing certificate

This is an independently checked, finite reduction, not a new theorem or a purely conceptual proof. It relies on (i) the known odd-cycle identity H(T) = Σ_k 2^k α_k for the intersection graph of directed odd cycles, (ii) Camion's theorem that every strong tournament has a directed Hamiltonian cycle, and (iii) Busch's known exact lower bound for H of a strong n-vertex tournament. Sources: Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, https://www.kurims.kyoto-u.ac.jp/EMIS/journals/EJC/Volume_13/PDF/v13i1n3.pdf; Camion theorem is also used explicitly as a classical result, not proved here.

Busch's bound equals 3·5^((n-3)/3) for n=0 mod 3, 5^((n-1)/3) for n=1 mod 3, and 9·5^((n-5)/3) for n=2 mod 3 (n≥3). It implies H≥25 for n≥7. Strong components of any tournament are linearly ordered, so their H values multiply; singleton components have H=1 and non-singleton components are strong tournaments of order at least 3, with H≥3. Hence an H=7 tournament would need a strong component with H=7; an H=21 tournament would need either a strong component with H=21 or components with counts H=3 and H=7. It suffices to exclude 7 and 21 among strong tournaments of order at most 6.

Every strong n-vertex tournament, n≥3, has a directed Hamiltonian cycle by Camion, so we may label one as 0→1→···→n−1→0. All 2^((n−1)(n−2)/2−1) orientations of the other arcs then cover every such tournament, possibly repeatedly. On n≤6 vertices, any odd cycle has length 3 or 5; the only pair of vertex-disjoint odd cycles is a pair of complementary directed triangles, and no triple of cycles is disjoint. Writing C for the number of directed odd cycles and D for the number of unordered disjoint odd-cycle pairs, the known identity becomes H=1+2C+4D. Each directed triangle has at most one disjoint triangle (its complementary triple), hence D≤C/2. Thus H=21 forces (C,D)=(10,0), (8,1), or (6,2).

The checker `cert.py` fixes the Hamiltonian cycle, enumerates all remaining orientations, counts every odd cycle once by fixing the smallest vertex as the cycle start, counts disjoint cycle pairs, and independently counts all Hamiltonian paths by permutations, checking the identity at every orientation. The complete strong cycle-containing H-value histograms are:

n=3: 3 (1 orientation)
n=4: 5 (4)
n=5: 9 (10), 11 (5), 13 (10), 15 (7)
n=6: 15 (12), 17 (12), 19 (12), 23 (36), 25 (24), 27 (12), 29 (66), 31 (36), 33 (72), 37 (126), 41 (24), 43 (60), 45 (20).

For the n=6 cycle-packing argument, the entire count restricted to C≤10 is (C,D): number of orientations:
(6,1):6; (7,0):12; (8,0):6; (9,0):12; (9,1):24; (10,2):18.
Neither target 7 nor any forced pair for 21 occurs. These numbers count labeled orientations with the distinguished 6-cycle fixed, not isomorphism classes; the table is exhaustive because of Camion's theorem.

## Optional 16-case reduction for two triangle pairs

For a six-vertex tournament, D≥2 implies C≥10 and therefore H≥29. Indeed, two complementary directed-triangle partitions can be labeled as 0→1→2→0, 3→4→5→3, 0→1→3→0, and 2→4→5→2. To see the relabeling, choose two partitions, pick the triangle of the second partition meeting the first triangle in two vertices, name those vertices 0,1 in arc order, and call its third vertex 3. Name the remaining two vertices of the first triangle 2 and of its complementary triangle 4,5 in cyclic order. The four displayed triangles follow. These determine ten arcs, leaving just 0↔4, 0↔5, 1↔4, 1↔5, and 2↔3. Put r=1[0→4]+2·1[0→5] and s=1[1→4]+2·1[1→5]. For both directions of 2↔3 the counts (triangles, five-cycles) are as follows; rows and columns are r,s=0,1,2,3:

(6,4) (6,6) (7,7) (6,8)
(7,7) (6,8) (8,8) (6,6)
(8,6) (8,8) (8,6) (7,7)
(8,6) (7,7) (8,6) (6,4)

This table is a finite hand-checkable reduction, not a computer-free proof of the full gap: the cases (C,D)=(10,0),(8,1) remain. `check_twopairs.py` verifies the table for both orientations of the fifth free arc and checks the lower bound.

## Smaller symmetry classifications for D=0 and D=1

The following reduce the cases to seven and eleven representatives respectively. The completeness of these tables is still verified computationally by the supplied checkers, so they must not be represented as fully computer-free arguments.

**D=0, strong tournaments.** Fix the Hamiltonian cycle 0→1→2→3→4→5→0. Record the remaining arcs 02,03,04,13,14,15,24,25,35 in that order, writing 1 when the lower-numbered vertex dominates the higher. Rotate the cycle and also allow reversal of its labeling combined with global arc reversal. Both transformations preserve the fixed directed cycle and (C,D). The 512 orientations have 52 orbits in all, seven of which have D=0:

bitstring | triangles | C | orbit size
000000001 | 4 | 7 | 12
000000101 | 5 | 9 | 12
000000111 | 4 | 8 | 6
000001101 | 6 | 12 | 12
000100001 | 6 | 12 | 12
000100111 | 6 | 11 | 12
000110011 | 7 | 16 | 12

No representative has C=10. `check_nopair.py` verifies orbit completeness and that (C,D) are constant within each orbit. It verifies the displayed table for all 512 fixed-cycle orientations.

**D=1, all tournaments.** Fix the unique complementary directed-triangle pair as A=(0,1,2) and B=(3,4,5), both in cycle order. Encode the nine remaining arcs by a 3×3 binary matrix M, where M[i,j]=1 iff i→(j+3). Relabel by cyclically rotating rows and/or columns; exchange A and B, which sends M to 1−Mᵀ; and permit reversal of every arc combined with reversal of the cyclic labeling of each triangle, which sends M[i,j] to 1−M[−i,−j] (indices modulo 3). These transformations preserve (C,D). They give 22 classes of all 512 matrices; the 11 classes with D=1 are listed below. Slashes separate rows of M.

matrix | triangles | C | orbit size
000/000/000 | 2 | 2 | 2
000/000/001 | 4 | 6 | 18
000/000/011 | 5 | 9 | 36
000/000/111 | 5 | 11 | 12
000/001/010 | 6 | 12 | 18
000/001/011 | 6 | 13 | 36
000/001/110 | 7 | 14 | 36
000/001/111 | 6 | 12 | 36
000/011/110 | 7 | 16 | 36
001/010/100 | 8 | 20 | 6
001/010/101 | 8 | 19 | 36

No representative has C=8. `check_onepair.py` verifies orbit completeness and (C,D) invariance for the 512 matrices and checks all table entries. This D=1 classification does not require strong connectivity.
