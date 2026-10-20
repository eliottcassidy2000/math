# Euler tours, Hamiltonian closure, and marked torus winding

20 October 2026. Session 16 continuation. This is a staged research addendum,
not a claim that a knot invariant classifies tournament metagraphs.

## Status and inheritance

- **PROVED here:** the directed-line-graph counting identity; the edge-cocycle
  Euler homology factorization; the checkerboard even-flux obstruction; the
  cycle-reversal and paired-path flux laws; the scalar-count embedding hostile.
- **CITED, with hypotheses checked:** BEST and the directed matrix-tree theorem;
  the facial-flip/height-potential theorem of Zhang–Qian.
- **FINITE-EXACT:** the three specified order-seven models, their marked
  Hamilton-cycle slope polynomials, and the complete 2,640-orientation,
  25-flux-sector census of this particular embedded K7.
- **NOT CLAIMED:** all-order counts for these sectors, a Fibonacci bijection
  from a 273-state sector, graph-isomorphism invariance of marked knot type,
  or connectivity of the actual Hamiltonian-path-complement metagraph.

The incoming `spectrum_basins_20261020e/arithmetic.md` already proves the
Fibonacci matching interpretation, sum-of-squares identity, Farey update, and
marked torus framing `F_j F_(j+1)`. Those are inherited, not new results here.
Its `regular.md` supplies the three regular order-seven classes and their
triangle-move comparison. The independent census below does not assume that
classification. The subsequently received
`kernels_games_detours_20261020f/walks_knots.md` at shared commit
`c402d288ac84f695e0fb355444ee67bfb7917d4f` already proves the directed
Euler-word/Hamilton-word line-graph bijection and explains the information
lost by forgetting incidence ports or integer carries. That word bijection
is inherited overlap, not a novelty claim here. Section 1 makes its balanced,
cyclic-count convention explicit; the checkerboard obstruction and the marked
K7 flux census are separate additions.

The geometric companion
[geometric-forcing packet](../torus_voltage_and_halfplane_forcing_20261020/geometric_forcing_bridges.md), G2–G4, constructs
the triangular-lattice torus quotient used here. The paired-return companion
`session16-regular-flips/regular_return_theorems.md` supplies the path operation
whose edge-chain refinement appears in Section 5.

Closest mechanism: counts must retain the actual carrier. Hostile: the same
abstract directed C3 embeds along every torus knot. Corrected near miss:
Euler/Hamilton transport goes through the **directed line graph**, not the
same tournament. Least-used sidecar: integral homology of individual embedded
edges and the meridian–longitude marking. Concept board: raw H; closed versus
nonclosing paths; Euler circuits; edge-chain flux; face flips; marked slopes.

## 1. The exact Euler-to-Hamilton correspondence

Let D be a finite strongly connected directed graph with m>0 distinguishable
arcs, with equal indegree and outdegree at every vertex. Let E(D) count Euler
circuits modulo cyclic rotation of their arc words. Define the directed line
graph L(D) to have one vertex per arc of D and an arc e→f exactly when the head
of e is the tail of f. Let H count directed Hamiltonian vertex words.

**Proposition 1 (PROVED).**

    H(L(D)) = m E(D).

A Hamiltonian word of L(D) is precisely an ordered, all-edge directed trail
of D. Its initial and terminal vertices agree: for any open trail the
outgoing-minus-incoming used-edge count is +1 at its initial vertex and −1
at its terminal vertex, whereas using all arcs of balanced D gives zero
at every vertex. It is therefore an Euler circuit. Conversely, every Euler
circuit gives such a Hamiltonian word. Each cyclic circuit has exactly m
linearizations, since its arcs are distinguishable and each is used once.
This proves the formula. The Hamilton-cycle version is the direct bijection
between E(D) and Hamilton cycles of L(D), both modulo rotation.

The tournament condition is **not** preserved by L. For a tournament on at
least four vertices, take two arcs with disjoint endpoints: their line-graph
vertices are nonadjacent in both directions. Thus oddness of tournament H
must not be applied to H(L(D)).

For a regular tournament T on n=2r+1>=3 vertices, every vertex has indegree
and outdegree r and m=nr. It is strongly connected: a proper source union
S of strongly connected components would have every cross edge pointing
out, so summing outdegree minus indegree over S gives |S|(n−|S|)>0,
contradicting balance.

BEST and the directed matrix-tree theorem give, for any vertex v,

    E(T) = tau_v(T) ((r−1)!)^n,
    H(L(T)) = nr tau_v(T) ((r−1)!)^n,

where tau_v is the number of spanning in-arborescences rooted at v, equivalently
a principal cofactor of diag(outdegree)−adjacency. The convention is cyclic
Euler circuits, equivalently circuits with one specified arc first. See
[Creed–Cryan, Theorem 1, pp. 1–2](https://arxiv.org/pdf/1202.2156).

For n>=7, E(T) is even, because (r−1)! is even, while H(T) is odd by the
Rédei parity theorem. At n=3, (H,E)=(3,1); at n=5, the unique regular class
has (H,E)=(15,11). Therefore **no nontrivial regular tournament has E(T)=H(T)**
under these conventions. This is not a comparison with Euler trail length,
which is simply nr.

For completeness, regular order-five uniqueness is elementary. Fix v, with
out-neighbours x,y and in-neighbours a,b, and name x→y. Vertex y must beat
both a,b to have outdegree2. Vertex x must beat exactly one of a,b; name it a.
Then a→b is forced by a's degree. Every arc is now determined, giving the
usual regular carousel.

## 2. The sidecar that actually recovers torus winding

Fix a continuous edge map f:D→R²/Z² and a lift v~ of every vertex. For an arc
e:u→v, lift its image beginning at u~. Its endpoint is v~+z_e, for a unique
z_e in Z². For a closed walk W,

    [W] = sum_(e traversed by W) z_e in H_1(T²;Z).

**Proposition 2 (PROVED).** Every Euler circuit of a fixed embedded or mapped
D has the same homology

    omega(D) = sum_(e in E(D)) z_e.

If Euler circuits are weighted by x^p y^q for their homology (p,q), their
Laurent polynomial is exactly

    E_D(x,y) = E(D) x^(omega_1(D)) y^(omega_2(D)).

Proof: endpoint discrepancies telescope around a closed walk. An Euler
circuit uses each oriented arc once, independent of its order. Changing the
chosen vertex lifts by k_v replaces z_e by z_e+k_u−k_v; the total is unchanged
by balance. This also proves gauge invariance for every closed cycle.

Under Proposition 1 the same weight is transported as a **vertex weight**
on L(D), one weight per original edge. Consequently the homology-weighted
Hamiltonian-word polynomial of L(D) is m times the displayed monomial.

For Hamilton cycles of D itself, different subsets of edges are used, so
one must retain the full polynomial

    C_D(x,y) = sum_(directed Hamilton cycles C, modulo rotation) x^[C]_1 y^[C]_2.

This is a genuine winding distribution. Its value at (1,1) counts Hamilton
cycles, not Hamilton paths. Open paths have endpoint-relative displacement,
not an intrinsic closed-curve homology class.

An Euler circuit in an embedded graph usually revisits vertices. Its image
is then not an embedded knot. To replace it by a knot requires a local
transition/smoothing rule, and sometimes leaving the torus. A simple
essential closed curve on a torus has nonzero primitive homology, so a
nonzero nonprimitive omega is an immediate obstruction to any smoothing
that both stays simple on the torus and preserves homology.

### A sharp scalar-count hostile

For any coprime positive p,q, parametrize the corresponding standard torus
curve by t↦(pt,qt) modulo Z² and embed that torus in S³. Mark three points on
the curve and regard the intervening curved arcs as a directed C3. The
abstract tournament has H=3, Euler length 3, and E=1 for every p,q, while
the resulting torus knots have unbounded genus. Thus even the triple
(H,E,Euler length) does not recover the winding pair or knot type. Curved
edges are explicit here; a straight three-segment polygon is an unknot.
The classical positive torus-knot genus formula is
`g(T(p,q))=(p−1)(q−1)/2`, for coprime positive p,q; for a primary-paper
statement see [Jabuka–Van Cott](https://arxiv.org/abs/1907.12970).

Conversely, subdivision of a fixed embedded knot into N directed arcs gives
H=N and Euler length N on the directed N-cycle without changing knot type
(the graph is a tournament only when N=3). This distinguishes combinatorial
length from a topological invariant in the other direction too.

## 3. Checkerboard faces force even Euler flux

Let an undirected graph G be cellularly embedded in an oriented torus and
have a proper black/white face colouring: each edge separates a black and a
white face. Orient every edge according to the boundary orientation of its
black face. Call the resulting orientation D0. Its total chain is the
boundary of the sum of all black faces, so J(D0) is a cycle and

    [J(D0)] = 0.

**Proposition 3 (PROVED).** For every Eulerian orientation D of G,

    omega(D) belongs to 2 H_1(T²;Z).

Indeed J(D)−J(D0)=−2R, where R is the chain of reversed edges oriented as in
D0. Since both full chains are cycles, 2∂R=0; integral vertex chains are
torsion-free, hence ∂R=0. Thus [J(D)]=−2[R]. A nonzero even vector is not
primitive. Accordingly **no Euler tour of any such orientation can become
a single simple essential torus curve without changing its homology**.
Zero homology also excludes a simple essential torus curve; a contractible
simple curve is not excluded by homology alone.

The proof extends verbatim to any closed oriented surface. On higher-genus
surfaces zero-homologous simple curves may be noncontractible, so the last
sentence is specifically a torus conclusion. This proposition makes no
claim about a 3-dimensional smoothing that leaves the torus.

## 4. The triangular torus K7: real slopes, not scalar analogies

Let a=(1,0), b=(0,1) denote triangular-lattice basis vectors (geometric
angle 60 degrees). Put

    Lambda = {(s,t) in Z² : s+3t = 0 modulo 7},
    u=(7,0), v=(-3,1).

The quotient R²/Lambda is a torus. Its triangular lattice graph has seven
vertices, labelled by i=s+3t modulo7. Its six edge steps have residues
±1,±2,±3, hence its graph is K7. The fourteen triangular faces are the seven
translates of (0,a,b) and of (a,b,a+b); they are checkerboard-coloured by
these two kinds. This is the companion geometric packet's construction.
The basis (u,v) is the homology marking used in every number below.
The geometric companion uses the different kernel basis
`g1=(-3,1), g2=(1,2)`. The exact change is `g1=v, g2=u+2v`:
coordinates `(alpha,beta)` in that basis become `(beta,alpha+2beta)` here.
One must make this basis change before comparing slopes or marked knot types.

For an edge i→j, let (s,t) be its actual unit triangular-lattice step.
Choosing vertex lift (i,0), its edge cocycle is

    z_(i,j) = ((i+s−j+3t)/7, t).

The numerator is divisible by7. Use step a for residue1, b−a for residue2,
b for residue3, and the negatives for residues6,5,4 respectively.

Define the carousel orientation by steps {a,b,b−a}, residues {1,3,2}.
Define the Paley orientation by {a,−b,b−a}, residues {1,4,2}. Each is regular.
Their Euler fluxes are immediate, without enumerating any tours:

    omega(carousel) = 7(a+b+(b−a)) = 14b = 6u+14v,
    omega(Paley) = 7(a−b+(b−a)) = 0.

Paley is exactly the orientation of the seven black triangular boundaries.
Proposition 3 applies to every regular orientation of this fixed K7 embedding.

The independently checked counts are:

| specified model | raw H | directed Hamilton cycles C | nonclosing paths R | tau | E modulo rotation | omega |
|---|---:|---:|---:|---:|---:|---|
| carousel |175|17|56|379|48,512|(6,14)|
| carousel with triangle (0,2,4) reversed |171|15|66|386|49,408|(4,8)|
| Paley |189|24|21|392|50,176|(0,0)|

All have Euler length21. The first two already refute nondecreasing monotonicity of Euler
circuit count as a function of raw H, even at fixed order and degree. The
three together refute monotonicity in either direction.

The precise relation between H and the slope polynomial is

    H(T) = n C_T(1,1) + R(T),

where R counts Hamilton paths whose final vertex does not beat their initial
vertex. Each directed Hamilton cycle has n cuts and every closing Hamilton
path arises uniquely this way. In particular Paley's R=21 does not mean a
tournament with H=21 exists.

For this fixed marking the cycle polynomials are **FINITE-EXACT**:

    C_car(x,y) = x + 14 x²y⁵ + x²y⁷ + x³y⁷,

    C_pal(x,y) = x^(−3)y^(−7) + 7x^(−1)y^(−2) + 7y^(−1)
                 + x + 7xy³ + x²y⁷.

All cycles here are simple embedded curves in the crossing-free torus graph.
Now, and only now, choose a homeomorphism of the marked quotient torus to a
standard unknotted torus in S³, carrying u,v to its chosen meridian–longitude
basis. A cycle with primitive slope(p,q) becomes the torus knot T(p,q), up to
the declared convention and orientation. Under this fixed choice, carousel
has sixteen nontrivially knotted directed Hamilton cycles: fourteen T(2,5),
one T(2,7), and one T(3,7). Paley has two: T(2,7) and T(3,7); the other
22 have an absolute winding coordinate1 or are a primitive coordinate axis,
so are unknots. This is not a graph-isomorphism invariant without the marking.

### Paley's nonclosing paths form one orbit; knot type does not

The affine maps i↦ai+b, with a in {1,2,4} and b in F7, are21 automorphisms.
Any stabilizer of0 restricts to an automorphism of the directed C3 on its
out-neighbours {1,2,4}. This restriction is injective: the three remaining
vertices3,5,6 are distinguished by having out-neighbour sets {4},{2},{1}
inside {1,2,4}. Thus the stabilizer has at most3 elements and Aut(Paley)
is exactly this order21 affine group. Its action on ordered Hamilton paths
is free. Since R=21, all nonclosing paths form exactly one orbit.

Multiplication i↦2i is realized on the triangular torus by the lattice map

    A=[[-1,−1],[1,0]],       a↦b−a, b↦−a.

On (u,v)-homology its matrix is

    M=[[2,−1],[7,−3]],       M³=I.

It cycles the three slopes (1,0)→(2,7)→(−3,−7)→(1,0).
Consequently a single graph-automorphism orbit of directed Hamilton cycles
contains an unknot, T(2,7), and T(3,7) under the chosen S³ marking. This torus
homeomorphism cannot extend to an ambient homeomorphism preserving that
embedded marked torus in S³, because it changes knot type. Quotienting by
abstract graph automorphisms therefore loses exactly the geometric data
needed for this knot interpretation.

## 5. Local reversals, face potentials, and paired path returns

If a directed cycle C is reversed, its chain contribution changes from C
to −C. Therefore

    omega(D after C reversal) = omega(D) − 2[C].

Contractible facial flips preserve flux; winding cycle reversals can change
it without changing any vertex degree. For the carousel directed triangle
0→2→4→0, its three steps are (b−a),(b−a),b, so

    [C] = −2a+3b = u+3v,
    (6,14) − 2(1,3) = (4,8).

The raw H change175→171 is verified independently. This is a winding
triangle, not one of the fourteen faces.

For the operation K_P that reverses every nonpath arc and leaves the directed
Hamilton path P unchanged, the chain identity is

    J(K_P D) = −J(D)+2P.

After a second path Q it becomes

    J(K_Q K_P D) = J(D)−2P+2Q.

When P,Q have the same ordered endpoints, Q−P is an integral cycle and the
Eulerian return has flux change 2[Q−P]. This is the geometric sidecar to the
paired-return companion's score theorem. The paths are evaluated in the same
fixed embedded underlying graph, with their current orientations. Nothing
here identifies these derived two-step returns with a single original
Hamiltonian-path metagraph edge.

**CITED structural connection.** Zhang–Qian's Theorem3.4 applies to a cellular
surface graph with strongly connected prescribed outdegrees and with every
edge incident to two distinct faces. With no forbidden face, two orientations
are connected by counterclockwise facial flips exactly when they are
homologous. In the proof, the integral face potential z for their difference,
normalized up to a constant, gives directed flip distance
`sum_f (z(f)−min z)`. All hypotheses hold for this K7 triangulation and regular
outdegree3. Since integral torus homology is torsion-free, equality of our
full-chain flux is equivalent to their orientation-homology convention
(which uses half the full-chain difference).
[Primary source, pp.3–9 and Theorem3.4](https://arxiv.org/pdf/2011.07481).

Thus fixed flux is a structural explanation of the face-flip components,
not merely a numerical component label. This theorem is cited, not newly
proved here. Its connectivity concerns facial flips of this fixed embedding,
not arbitrary triangle flips or the actual R/K path-complement metagraph.

## 6. Complete small-order certificate and the second appearance of273

The checker generates **all labelled regular orientations of K7** by assigning
one vertex's remaining wins at a time and pruning impossible degree targets.
It does not sweep2^21 orientations and does not assume the inherited three
isomorphism classes. It finds2,640 states, with raw-H populations

    H175:720,       H171:1680,       H189:240.

There are25 realized Euler fluxes. They are listed individually, with their
raw-H distributions, in `euler_torus_controls.json`. Their fibre/component
size multiset is

    1 of size660;   6 of size273;   6 of size42;
    6 of size14;    6 of size1.

Each of the six273-state sectors has77 orientations with H175,175 with H171,
and21 with H189. The zero-flux sector contains respectively168,420,72 of
those three types. Explicit zero-flux representatives are in the JSON. Hence
**flux does not determine raw H**, even within this fixed geometry.

The earlier Fibonacci matching model also has13·21=273 states. This equality
is real but is **not** presently a structural identification. There is already
a decisive hostile to the most natural graph-isomorphism claim:

- The torus face-flip graph in flux(0,2) has273 vertices and770 edges, with
  degree multiplicities `3:14, 4:56, 5:56, 6:70, 7:42, 8:35`.
- The add/remove-one-matching-edge graph of P6 disjoint union P7 has273
  vertices and914 edges, with maximum degree11.

Those natural graphs are not isomorphic. The first failed implication is
"same cardinality implies the same local move graph". The surviving question
is whether a nontrivial, explicitly specified weighted or quotient map
exists. No Fibonacci recurrence for the torus sectors is claimed.

The complete face-flip spanning forest is stored in
`face_flip_certificate.json`: sorted orientation masks, all fourteen faces,
and each child/parent/face triple. It certifies connectivity independently
of any visual layout, and is replayed by `replay_certificate.py`.
That replayer independently proves completeness using the coefficient of
`x_0³ ... x_6³` in `product_(i<j)(x_i+x_j)`: assigning an edge's factor to
its winner counts its outdegree. The coefficient is2,640, and its largest
dynamic-programming layer has only713 degree states.

## 7. Reproduction, proof boundaries, and sources

From this directory:

    python3 verify_euler_torus.py
    python3 -O verify_euler_torus.py
    python3 replay_certificate.py

The main checker uses independent permutation and subset-DP Hamilton counts;
permutation and displacement-DP winding counts; Bareiss determinant and direct
rooted-choice arborescence counts; and direct local transition enumeration
(6^7 systems per regular order-seven model) to audit BEST. It directly counts
Hamilton paths of the line graph at orders3 and5. Gauge changes, cycle reversal,
Paley automorphisms, all degrees, all face-flip fluxes, and the273 graph hostile
are checked. All verification conditions remain active under Python `-O`.
No scaling beyond the2,640-state universe is part of this packet.

Primary imported sources:

1. [Creed–Cryan, The number of Euler tours of random directed graphs,
   arXiv:1202.2156](https://arxiv.org/pdf/1202.2156), Theorem1 and matrix-tree
   discussion. Imported only for the stated counting formula and convention.
2. [Zhang–Qian, Flips on homologous orientations of surface graphs with
   prescribed forbidden facial circuits, arXiv:2011.07481v2](https://arxiv.org/pdf/2011.07481),
   standing hypotheses in Section2 and Theorem3.4. Imported only within those
   hypotheses, which are checked above.
3. [Grinberg–Stanley, The Rédei–Berge symmetric function of a directed graph,
   arXiv:2307.05569](https://arxiv.org/abs/2307.05569), the inherited tournament
   parity theorem; see also the repository's THM-002-ocf.md and
   CORE-PAPERS.md tournament attribution guardrail.

For the inherited Fibonacci/torus-framing source and classification boundaries,
use `spectrum_basins_20261020e/arithmetic.md`, not this packet as a claim of
new discovery. A relevant but unused literature tangent is Doleshal's
[Twisted torus knots with Horadam parameters, arXiv:2501.17850v4](https://arxiv.org/html/2501.17850v4):
it concerns an additional twisting operation, not an identification of our
Euler tours or Hamilton paths with Fibonacci torus knots.
