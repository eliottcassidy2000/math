# Paley7 blue/black path symmetry is balance between two Fano designs

Status: PROVED structural reductions plus FINITE-EXACT exhaustive nine-orbit certificate. Specific to the Paley tournament on seven vertices; no all-order equivalence is claimed.

## Inheritance

- `01-canon/theorems/THM-211-paley-design.md` proves that the fourteen cyclic triples of Paley7 are two distinct cyclic Fano planes, together a 2-(7,3,2) design. This is inherited, not new.
- `01-canon/theorems/THM-4553-the-point-stripping-ladder-psl27-borel-torus-and-the-collatz-frobenius.md` identifies Aut(Paley7) with the order21 affine square-multiplier group.
- `05-knowledge/results/fano_tournament_move_quotients_20261021/session23b-fano/README.md` gives the fixed-Fano Q7 line-flip cube, its Hamming/syndrome quotient, and marked-path ancestry. The cube construction and regularity preservation are inherited.
- The current conference note `session35-zeta/CONFERENCE_COLLAPSE.md` recovers two distinct fixed-gauge Fano systems from source/sink switching codes. Their labeling must be transported explicitly before identifying them with a chosen set of triples here.

Targeted searches of theorem/reflection filenames and contents found no previous unused-line-balance characterization of blue/black paths. This is a bounded provenance finding, not a claim of literature priority.

## 1. Definitions and structural reduction

Let V=F7, Q={1,2,4}, and orient i→j iff j−i∈Q. Define
L+={i+{0,1,3}:i∈F7}, L−={i−{0,1,3}:i∈F7}.
These are exactly the inherited two cyclic Fano planes (a translate of Q and its negative). For a Hamiltonian path P=(p0,...,p6), let u± be the number of lines in L± containing no unordered path edge {p_i,p_{i+1}}.

Define BLUE to mean that reversal of the ordered path extends to an anti-automorphism of the tournament:
M[p_i,p_j]=−M[p_{6−i},p_{6−j}] for every i,j.
BLACK means this fails. This is the precise mirror-symmetry predicate, not merely self-converseness of the underlying unmarked tournament.

All automorphisms are x→ax+b with a∈Q; they preserve each Fano plane separately. Anti-automorphisms have nonsquare a. Path reversal is an involution, so its affine multiplier satisfies a²=1. Among nonsquares only a=−1 works. Therefore

BLUE iff p_i+p_{6−i}=b for all i,

where b=p0+p6=2p3. Thus BLUE is tested by a specific anti-affine reflection, not by searching arbitrary bijections.

The map x→−x+b swaps L+ and L− and preserves the set of path edges when BLUE holds. Hence BLUE implies u+=u− without enumeration. The converse requires the P7-specific certificate below.

## 2. Unused lines are adjacent-step turn counts

A Fano line contains at most two path edges: all three would form a triangle inside a simple path. If a_j denotes the number of lines containing j path edges, then
 a0+a1+a2=7, a1+2a2=6,
so u=a0=1+a2.

When a line contains two path edges, the edges share a vertex and are consecutive on P. Consequently a2 counts the length-two path segments whose three vertices form a line of that plane.

Write d_i=p_{i+1}−p_i∈{1,2,4}, i=0,...,5. Equal consecutive increments produce a transitive triple. Distinct increments produce a cyclic triple; their ratio determines its plane:

 d_{i+1}/d_i=2 ⇒ triple in L+,
 d_{i+1}/d_i=4 ⇒ triple in L−.

Therefore the exact formulas are
 u+=1+#{i: d_{i+1}=2d_i mod7},
 u−=1+#{i: d_{i+1}=4d_i mod7}.

Meanwhile BLUE iff the six-letter increment word is a palindrome. This identifies the sought balance as equality of the two directed turn counts on the multiplicative three-cycle Q, subject to the strong constraint that the seven additive partial sums are distinct. Arbitrary words with balanced turn counts need not be palindromic; the Hamiltonian partial-sum constraint is essential.

## 3. Complete nine-orbit certificate

Every Hamiltonian path has a unique affine-square normalization p0=0,p1=1. Indeed d0∈Q and x↦(x−p0)/d0 belongs to the automorphism group. The action is free on complete ordered paths. Thus each normalized representative has exactly21 paths in its orbit.

The nine normalized representatives are:

path       increments   (u+,u−)  BLUE  reflection b
0123456    111111       (1,1)    yes   6
0124563    112114       (2,3)    no    none
0124635    112242       (3,2)    no    none
0126345    114411       (2,2)    yes   5
0134526    121144       (2,3)    no    none
0135246    122422       (3,2)    no    none
0152346    144112       (3,2)    no    none
0152463    144224       (2,3)    no    none
0152634    144441       (2,2)    yes   4

Completeness can be checked by the120 permutations of the remaining five labels after0,1, retaining exactly those with six legal arcs, or by243 increment continuations after d0=1 retaining distinct partial sums. The supplied checker also enumerates all5040 vertex permutations and independently finds189 paths and21 automorphisms.

It follows exactly that:

BLUE iff u+=u−;
BLACK iff |u+−u−|=1.

The path counts are21 with(1,1),42 with(2,2),63 with(2,3),63 with(3,2). Equivalently the nine path orbits split1,2,3,3. The first two sets are precisely the63 BLUE paths; the latter two are the126 BLACK paths.

## 4. The path-preserving line cube

Fix ONE of the Fano planes and the ordered path P. Whole-line reversal preserves P iff the chosen line contains no path edge. Lines in one Fano plane partition all unordered vertex pairs, so different line moves have disjoint edge supports. The P-preserving whole-line choices are exactly an induced coordinate cube Q_{u±}, not merely a lower bound on available moves.

Every line is a directed triangle and its reversal preserves every vertex's score: each vertex on the line keeps one incoming and one outgoing edge on it. Hence every vertex of this path-preserving cube remains regular. This is the EVEN-projective-order case q=2. The checker verifies all1932 path/plane/cube-state instances exactly.

The two cubes have dimensions(u+,u−). Thus BLUE means equal dimensions of the two naturally available move cubes; BLACK means they differ by exactly one. This is the concrete metagraph bridge: path mirror symmetry is measured by balance between two independent design-specific mutation spaces.

Do not merge these two coordinate cubes into Q_{u++u−}: a line from each plane can share an edge, and the combined move supports are not disjoint. Nor is a whole-line reversal the same as flipping the original off-path complement. The latter reverses every off-path edge; a line move reverses only a specified three-edge block. Even the union of all unused lines can omit off-path edges lying on a line that contains a path edge.

## 5. General parity-sidecar implication, carefully scoped

For any edge-partition by blocks, a fixed path leaves a cube of whole-block choices indexed by blocks unused by its edges. This statement is combinatorial and independent of regularity.

If every block induces a regular tournament (as for odd block size in an even-order projective plane with an explicitly chosen regular local orientation), each block reversal preserves scores, so the entire unused-block cube stays regular.

Conversely, in a family for which a separate theorem proves that only the empty and full line-flip subsets preserve regularity, any nonempty path-preserving line-flip subset destroys regularity. It cannot be the full set because a nonempty path has an edge lying on some line, making that line forbidden. This last implication is elementary, but the assumed two-state regularity theorem must be supplied separately; projective-order parity by itself is not substituted for that theorem here.

## Verification and scope

Run `python3 session35-recursions/fano_paths.py`. Output `fano_paths.out` includes the nine representatives, unused-line lists, symmetry values, exact group checks, and all1932 regularity/path-preservation checks. Explicit exceptions rather than removable assertions enforce tests.

This note proves a Paley7 marked-path theorem and recovers the prior line-cube mechanism. It does not identify every self-converse tournament with a balanced two-design pair, prove blue iff balance at other orders, or transfer a theorem about whole-line flips to off-path-complement flips.

## 6. Intrinsic local wiring degree and exact departure from Paley

Because L+ and L− are disjoint and exhaust all fourteen directed triangles, the sum

 d_P=u++u−

is intrinsic: it counts directed-triangle reversals that preserve the fixed path P. Every such move also preserves regularity. In the graph of regular tournament presentations containing P, with adjacency given by one directed-triangle reversal avoiding P's edges, d_P is the degree at the Paley presentation.

The nine-orbit certificate now gives a design-independent reformulation:

At a Paley7 presentation, BLUE iff d_P is2 or4, and BLACK iff d_P=5.

This is a genuine local wiring relation with a preserved Hamiltonian path, distinct from the original off-path-complement operation. Its incidence count is
21·2+42·4+126·5=840=14·60.
Every directed triangle is avoided by exactly60 Hamiltonian paths, checked separately in the executable certificate.

Each one-triangle reversal leaves the Paley isomorphism class and enters the regular class with H=171. This transition is inherited from the session23b fixed-Fano Q7/Hamming quotient and verified here by exact path counting, not inferred solely from regularity.

Stronger: for either chosen Fano plane, the whole path-preserving cube contains Paley ONLY at its origin. Every nonempty permitted flip set has H=171. Here is the structural explanation using the inherited Paley Hamming code. A path has u≤3 by the nine-orbit certificate. Nonzero Paley-to-Paley flip words in the fixed-Fano cube have minimum weight3; the weight-three supports are precisely pencils of three concurrent original lines. Such a pencil cannot be unused by P, since it contains every edge incident with the common point, which every Hamiltonian path must enter or leave. Hence no nonzero permitted word lies in the Paley code. The inherited two-class census identifies all remaining words as H171.

The direct checker verifies378 origins (=2·189 choices of path and plane) with H189 and1554 nonorigin states with H171, for1932 path/plane/state incidences total. These cover99 distinct labelled presentations; the incidence count is not a count of distinct tournaments.

## 7. Same-order regular hostile: the criterion is not universal

Take the regular seven-vertex tournament with outneighbor bitmasks
[98,84,41,19,37,74,28].
It has H=171. Two valid Hamiltonian paths are

P=(0,1,2,3,4,5,6): BLUE, with five allowed triangle reversals;
Q=(0,1,4,2,5,6,3): BLACK, also with five allowed triangle reversals.

The executable independently checks the paths, direct mirror predicate, triangle avoidance degree and H171. Thus even at the same order and within regular tournaments, local degree5 does not generally determine color. The Paley two-design structure and its nine admissible orbit types are load-bearing. No broad degree–color conjecture is supported by the current result.

## 8. Uniform survivor: affine-reflection symmetry forces even local degree

Let p be any prime congruent3 mod4, use the Paley tournament on Fp, and let P be a Hamiltonian path reversed by the affine reflection r(x)=b−x. Let d_P count cyclic triangles containing no path edge. Then d_P is EVEN.

Proof: r is an anti-automorphism, so it acts as an involution on cyclic-triangle vertex sets and preserves the path-edge avoidance predicate. Nonfixed triangles occur in pairs. Write c=b/2. A fixed three-element subset must be{c,c+a,c−a} for a≠0, modulo a↔−a. Such a triple is cyclic iff the quadratic character χ(−2)=+1. By the elementary quadratic-character supplementary law, this holds exactly when p≡3 mod8; when p≡7 mod8 no fixed cyclic triples exist.

For p≡3 mod8 there are(p−1)/2 fixed cyclic triples. Exactly ONE contains a path edge: the triple consisting of the middle vertex c and its two adjacent path vertices. To see uniqueness, an edge between c+a and c−a would join positions interchanged by reversal; their position sum is p−1, which is even, whereas consecutive positions have odd sum. Thus no such edge occurs on P. The only path edges incident with c are the two middle edges, and they lie in the single indicated triple. Hence the number of fixed admissible triangles is(p−3)/2, which is even for p≡3 mod8. In both congruence cases, the full admissible count is even. QED.

Scope: odd d_P certifies the absence of this affine-reflection reversal. To replace that with the absence of EVERY possible reversing anti-automorphism at an arbitrary Paley order, one must separately justify the automorphism-group classification being used; no such unproved extension is implicit here. The P7 classification above already supplies that input at seven.

Even parity is not sufficient at larger Paley order. In Paley11 with quadratic-residue arcs, the Hamiltonian path
(6,4,5,3,7,0,1,10,8,9,2)
has d_P=32 but fails the direct mirror anti-automorphism predicate. The supplied checker verifies its legal arcs, all cyclic-triangle avoidance counts and failed symmetry exactly. This preserves the broad necessary parity law while decisively rejecting a general parity iff symmetry claim.

Additional replay strengthening: both complete128-word fixed-Fano cubes are independently checked in the current labeling. Each has16 H189 words and112 H171 words; the H189 words form a binary linear code with weight distribution1+7t³+7t⁴+t⁷, and every weight-three support is a concurrent pencil. This verifies the inherited code sidecar used in section6 without relying on an unrecorded labeling identification.
