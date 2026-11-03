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
