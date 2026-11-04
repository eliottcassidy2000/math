# Planar versus spherical symmetry, Maschke, and the actual metagraph edge

Session 190, 2026-11-04. **PROVED elementary statements and FINITE-EXACT controls; independent review pending.** Classical classifications are explicitly cited. This is a connection packet, not a literature-priority claim and not a solution of the raw Hamiltonian spectrum conjecture.

## 0. Inheritance and the useful distinctions

The closest inherited mechanism is [LEM-003, `aut-acts-freely-on-ham-paths`](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md): a directed spanning path has no nonidentity automorphism. The hostile example is Paley7: its group has order21 but is noncyclic, so it cannot be a faithful finite rotation group in real three-space. The corrected near miss is equating five Platonic solids with five rotation-group families: dual solids have the same rotation group. The useful sidecar is the coefficient field in Maschke and the orientation sign in a signed action.

Relevant recovered work:

- [THM-3447, `regular-group-orbit-tournament-parity-and-involution-ties`](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-3447-regular-group-orbit-tournament-parity-and-involution-ties.md) classifies invariant orientations on a regular group orbit.
- [THM-4644, `fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion`](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md) already proves the Fermat-prime/translation-only circulant equivalence. It is inherited, not new here.
- [THM-4553, `the-point-stripping-ladder-psl27-borel-torus-and-the-collatz-frobenius`](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-4553-the-point-stripping-ladder-psl27-borel-torus-and-the-collatz-frobenius.md) already identifies the Paley group and two invariant tournament orientations.
- [Gilbreath/Fermat/Platonic synthesis](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/gilbreath_fermat_platonic_20260926.md) already separates the two fivefold counts.
- [THM-1370, `h-spectrum-omits-7-21-all-n`](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md) excludes raw counts7 and21 by strong-component floors and multiplicativity, not by group-order numerology.

The live concepts are: finite geometric actions; orientation reversal; free marked-path actions; representation characteristic; cyclic multiplier groups; orbit versus quotient information.

## 1. What the five families actually classify

The user wrote `{Z_n,Z_n,A_4,S_4,A_5}`. If the second entry was intended as the dihedral family, the standard list is

    C_n, D_n, A_4, S_4, A_5,

where `|D_n|=2n`. This is the classification of finite subgroups of **SO(3)**, up to conjugacy. It is not a classification of planar point groups. Finite subgroups of **SO(2)** are cyclic; finite subgroups of **O(2)** are cyclic or dihedral. A finite group of Euclidean plane isometries has a common fixed point, obtained by averaging any point's orbit, and hence reduces to the O(2) case after translating that point to the origin. Sources: [Etingof, Theorem10.2](https://math.mit.edu/~etingof/groups.pdf) and [Arapura, Theorem14.1](https://www.math.purdue.edu/~arapura/algebra/algebra14.pdf).

The degenerate entries need care: `C_1=1`, `C_2` is a half-turn group, and `D_1=C_2` is a reflection group in the plane. The last two are abstractly isomorphic but have different determinant characters in O(2). In SO(3), the standard embedding of a plane reflection adds a sign on the perpendicular axis, making it a half-turn. Thus counting abstract groups, marked geometric actions, and parameterized families gives different answers already at order2.

The connection to the five Platonic solids is explicit but not a five-to-five bijection. For Schlaefli symbols `{p,q}`, spherical positive curvature gives

    1/p + 1/q > 1/2,    p,q >= 3.

The solutions are `(3,3),(3,4),(4,3),(3,5),(5,3)`. Duality interchanges p and q. Its three orbits have rotation groups `A4,S4,A5`; tetrahedron is self-dual, cube/octahedron share S4, and dodecahedron/icosahedron share A5. The two infinite families C_n,D_n do not supply the remaining two Platonic solids. The precise preserved object is the rotation action; the lost coordinate is vertex/face duality.

## 2. Odd tournament symmetry forbids a literal Platonic realization

**Proposition1.** Every finite tournament automorphism group has odd order. Consequently every subgroup of a tournament automorphism group that embeds faithfully into SO(3), and in particular O(2), is cyclic of odd order.

**Proof.** A nonidentity involution of a finite vertex set has a transposed pair. It reverses that pair's arc and cannot preserve a tournament. Cauchy's theorem excludes even group order. In the SO(3) list, every noncyclic group has even order. QED.

This is a necessary geometric obstruction, not a restriction that all tournament groups be cyclic. Paley7 supplies the noncyclic odd group `C7 semidirect C3`.

If a group preserves T up to global converse, its orientation-preserving kernel has index at most2 and odd order. Any faithful SO(3) realization of that signed group must therefore be cyclic of odd or twice-odd order, or dihedral D_m with m odd. None of A4,S4,A5 qualifies. This statement concerns a faithful geometric realization of the group, not whether an abstract tournament can be drawn with crossing edges on paper.

The underlying undirected graph of an n-vertex tournament is K_n, so ordinary graph planarity also fails for n>=5. This is a separate issue from symmetry representation.

## 3. A literal role for 1, C2, 7 and21

Put `Q={1,2,4}` in F7, and define `P7` by `x -> y` iff `y-x in Q`. Set

    G = {x -> a x+b : a in Q, b in F7},          |G|=21,
    Ghat = AGL(1,7),                             |Ghat|=42.

The quadratic character of a is a homomorphism `Ghat -> {+1,-1}`. Its kernel is G; positive elements preserve P7 and negative elements carry it to its converse. Direct exhaustion over all5040 permutations confirms that these are the full ordinary and signed symmetry groups, independently reproducing the inherited Paley identification.

In the action of Ghat:

- a point has stabilizer C6 and orbit size7;
- an unordered pair has stabilizer C2 and orbit size21;
- an ordered distinct pair has trivial stabilizer and orbit size42.

Indeed an affine map is uniquely determined by the images of two distinct points. The nonidentity setwise stabilizer of `{x,y}` is `z -> x+y-z`, which exchanges the two points. Passing from ordered to unordered pairs is exactly the free reversal involution. A tournament is a choice of one orientation above each unordered pair; the Ghat action cannot preserve such a choice because this C2 swaps its two lifts. Restriction to G removes that obstruction: its21 arcs form one free orbit, and its vertex stabilizers have order3.

This is an actual orbit-stabilizer diagram linking the requested small groups and counts. It does **not** identify the point and edge counts7,21 with the forbidden raw Hamiltonian-path values7,21. In fact

    H(P7)=189=9*21,

and the full endpoint matrix has entries1 on P7 arcs,8 on reversed arcs, and0 on the diagonal. Thus every vertex starts27 paths. The9 free path orbits are the old LEM-003 mechanism, not a new special divisibility theorem.

For comparison, the eight labelled C7-invariant tournaments form the three-dimensional cube of choices from `{1,6},{2,5},{3,4}`. Six have `(H,|Aut|)=(175,7)`; the two Paley orientations have `(189,21)`. Reversing a complete translation orbit is a cube edge; reversing a single arc is not. This orbit-flip cube is not the path-complement metagraph discussed next.

## 4. Every actual metagraph edge destroys the common ordinary symmetry

For n>=3, define the **labelled path-complement metagraph** precisely. Its vertices are tournaments on a fixed n-element set. If `P=(v1,...,vn)` is a directed Hamiltonian path of T, let `T_P` keep the n-1 arcs of P and reverse every other arc. Join T and T_P. This operation is involutive with the same P. Different P give different neighbours: the common arc set recovers P and its directed order uniquely.

**Theorem2, edge rigidity.** For every such edge,

    Aut(T) intersection Aut(T_P) = 1.

More strongly, its unordered edge `{T,T_P}` has trivial stabilizer under the simultaneous relabelling action of S_n. Consequently Aut(T) acts freely on the actual set of distinct labelled neighbours of T, and for any nontrivial permutation group G the G-fixed tournament locus is an independent set in this metagraph.

**Proof.** The common arc set of T and T_P is exactly the directed spanning path P. Every permutation preserving the two endpoints either separately or as an unordered pair preserves that intersection. A directed path has a unique source and a uniquely determined successor at every step, so every vertex is fixed. The path-to-neighbour map is equivariant and injective, giving freeness and the fixed-locus conclusion. QED.

For P7 this yields189 distinct neighbours in nine G-orbits of size21. Every neighbour breaks every original nonidentity automorphism; it may have new automorphisms, which this theorem does not forbid. A quotient neighbour is not a labelled neighbour, and two G-orbits may become identified after forgetting the source and relabelling the target arbitrarily.

**Theorem3, the residual C2.** Enlarge the relabelling action to `S_n x C2`, where the second factor is global converse. The stabilizer of any unordered labelled metagraph edge is either1 orC2. Its only possible nonidentity element is

    global converse followed by r_P,
    r_P(v_i)=v_(n+1-i).

It exists precisely when r_P is an anti-automorphism of T, and then fixes T and T_P separately. No signed symmetry interchanges the two endpoints.

**Proof.** A signed symmetry preserving the edge must preserve its directed common path up to the applied global converse. In the positive case Theorem2 applies. In the negative case the underlying permutation is forced to be r_P, so there is at most one candidate, of order2. It cannot interchange the endpoints: the pair `{v1,vn}` is off P for n>=3, and converse followed by exchanging its endpoints leaves its original arc orientation unchanged, whereas T_P reverses it. Thus the candidate must fix each endpoint. QED.

This gives the requested pair `{1,C2}` directly at the level of a metagraph edge: ordinary edge symmetry is always1; signed residual symmetry is1 orC2. Forgetting path direction by itself would also give the usual undirected-path C2, but that is a different action.

### 4.1 All-order orbit counts for this edge carrier

Normalize P to `0->1->...->n-1` and put

    q=(n-1)(n-2)/2,
    f=floor((n-1)/2),
    c=(q+f)/2.

There are q off-path bits. Interchanging the endpoints of the *metagraph edge* identifies the bit word x with its complete complement; this is not a relabelling of the tournament. Hence the number of S_n-orbits of unordered labelled metagraph edges is

    E_n = 2^(q-1).                                      (1)

Global converse followed by renormalizing the reversed path acts on these q coordinates by `(i,j) -> (n-1-j,n-1-i)`. It has f fixed coordinates and c cycles. Its fixed bit words number2^c. No word is sent to its complement, since the endpoint pair `{0,n-1}` is a fixed coordinate. Thus precisely2^(c-1) ordinary edge-orbits have C2 signed stabilizer. Burnside gives

    E_n^signed = (2^(q-1)+2^(c-1))/2.                   (2)

These formulas hold for n>=3; they do not assert anything about black-edge support, vertex-class reconstruction, connectivity, or a cut-switching aggregate. The first rows `(n,E_n,E_n^signed,C2 classes)` are

    (3,1,1,1), (4,4,3,2), (5,32,20,8),
    (6,512,272,32), (7,16384,8320,256),
    (8,1048576,525312,2048).

Finite verification uses two representations: all labelled tournaments and their actual edges for n=3,4,5; and the normalized off-path cube for n=3,...,7. The n=8 row is a formula evaluation, not a cube census.

### 4.2 The Paley signed control

Its seven anti-automorphism involutions are `x -> b-x`. Each fixes exactly9 Hamilton paths under path reversal, by exact enumeration. The189 paths split under Ghat's signed action into three orbits of size21 and three of size42. Equivalently, the nine ordinary G-orbits become six signed orbits: three remain individually fixed and six pair up. This is a finite control, not an all-Paley-prime formula.

## 5. Maschke says which field permits splitting, not which shape exists

Maschke's theorem states that every finite-dimensional representation of G over a field k is completely reducible exactly when `char(k)` does not divide `|G|`. Characteristic0 always qualifies. The usual proof averages a projection by `|G|^-1 sum_g gPg^-1`. It does not classify geometric groups or produce a preferred orientation. See [Etingof et al., Theorem4.1.1 and Proposition4.1.2](https://math.mit.edu/~etingof/representationtheorybook.pdf).

In characteristic2, the trivial group and every odd tournament group are in the semisimple regime. C2 is not: its swap matrix on F2^2 is nonidentity with `(g-I)^2=0`; the constant line has no invariant complementary line. Hence the field characteristic separates the ordinary tournament symmetry from the residual reversal precisely where averaging by2 would fail. Over characteristic not2, the familiar projections `(I+r)/2` and `(I-r)/2` split the involution's two eigenspaces; over F2 they are unavailable.

For F21, characteristics3 and7 are the bad ones. In characteristic7, even the natural seven-point permutation representation does not split off its constant line: a translation-invariant linear functional has all equal coefficients and sends the all-ones vector to7a=0. A projection onto constants would have to send that vector to itself. In characteristic3, the nonsplit regular representation of the quotient C3 can be inflated to F21. These are concrete failures, not a claim that every representation in bad characteristic is nonsplit.

### 5.1 The sharp real dimension is6, with an SU(3) model

Let S be the integer skew adjacency matrix of P7, with entries+1 on arcs and-1 in the opposite direction. Direct multiplication gives

    S 1 = 0,       S^2 = J - 7I,

where J is the all-ones matrix. On the real zero-sum space `W=1^perp`, of dimension6, the operator `I_c=S/sqrt(7)` squares to-I and is orthogonal. It therefore makes W a complex three-dimensional Hermitian space. Ordinary automorphisms commute with I_c and are complex-linear/unitary; anti-automorphisms anticommute with it and are conjugate-linear/antiunitary.

An explicit faithful complex model is

    D=diag(zeta,zeta^2,zeta^4),  zeta=exp(2*pi*i/7),
    P e_0=e_2, P e_1=e_0, P e_2=e_1.

Then `D^7=P^3=1`, `PDP^-1=D^2`, and both determinants are1. The21 matrices `D^a P^b` are distinct, giving an SU(3) embedding and hence a real SO(6) embedding.

**Proposition4.** The minimum dimension of a faithful real linear representation of F21 is6.

**Proof.** A faithful representation is nontrivial on its normal C7. After complexification, some nontrivial seventh-root eigencharacter k occurs. Conjugation by the order3 element forces the full character orbit `{k,2k,4k}`. Reality forces its conjugate orbit `{-k,-2k,-4k}`. These six nonzero residues are disjoint and exhaust F7*, so the complexified dimension, equal to the original real dimension, is at least6. The model above attains6. QED.

Thus the obstruction is sharper than failure to draw a symmetric heptagon: the complete Paley group cannot act faithfully by real matrices of dimension2,3,4, or5 at all. Its C7 subgroup does act by planar rotations. This is not a faithful planar action of the whole group.

### 5.2 The same group needs only three binary dimensions

Let t act by cyclically permuting the seven coordinates of F2^7. Since

    x^7-1=(x+1)(x^3+x+1)(x^3+x^2+1)

over F2, the natural permutation module is the direct sum of the kernels of these three coprime factors evaluated at t, of dimensions1,3,3. Multiplication of labels by2 preserves each cubic summand, so this is an F21 decomposition. Both cubic summands are irreducible: C7 is transitive on their seven nonzero vectors. Label negation exchanges the two summands.

Equivalently, write F8=F2[a]/(a^3+a+1) and identify label i in Z7 with a^i in F8*. The affine maps

    i -> 2^j i+b,       j=0,1,2,

become the F2-linear maps `x -> a^b x^(2^j)`. This realizes

    F21 = GammaL(1,8) < GL(3,2),     [GL(3,2):F21]=8.

The exact check constructs all168 invertible binary3x3 matrices and verifies that F21 is the normalizer of the Singer C7. Dimension3 is minimal over F2, since GL(1,2) and GL(2,2), of orders1 and6, cannot contain an element of order7.

This is a genuine characteristic-sensitive bridge: minimum faithful real dimension6, minimum faithful binary dimension3. It is not a contradiction, and the larger group GL(3,2) cannot itself preserve the tournament because it has involutions. Maschke applies to F21 over F2, not to GL(3,2) over F2.

**Lost-information warning.** Reducing S itself modulo2 erases the orientation: every off-diagonal entry becomes1. The binary decomposition uses the retained permutation action and its cyclic operator, not the sign-forgetting reduction of S. Negation swaps the two cubic pieces; Maschke does not prefer one or supply the missing tournament orientation.

## 6. Which other fives really connect?

The five known Fermat primes are `3,5,17,257,65537`; this is not a theorem that exactly five exist. The exact old tournament connection is

    p-1 is a power of2
      iff every circulant tournament of prime order p has Aut=C_p.

The reason is an involution obstruction: a nontrivial subgroup of the cyclic2-group F_p* contains-1, which cannot preserve a tournament connection half-set S. If p-1 has an odd prime factor r, choose S as one member of every paired orbit `{aH,-aH}` for an order-r subgroup H; then translations together with H preserve the tournament. The full-group conclusion uses the inherited prime-degree affine theorem, not merely this multiplier calculation. At p7 the odd multiplier C3 survives and produces21; at a Fermat prime no such multiplier survives. This gives a structural link between Fermat primes, the trivial multiplier group, and the Paley group rather than a five-element matching.

Classically Gauss-Wantzel also connects these primes to constructible polygons via `n=2^a product distinct Fermat primes`. The five Platonic solids instead satisfy the curvature inequality in Section1. The previously studied square-sum sizes `{15,16,17,23,25}` come from a Hamilton-path existence window and have no group-action map supplied here. Equal cardinality is not evidence for an invariant-preserving correspondence. The exact transfers retained here are duality, multiplier exclusion, stabilizer indices, and representation-field changes.

## 7. Scope and next useful tests

- The forbidden raw H values7,21 remain an inherited all-order result. This packet gives no new proof that those are the only missing positive odd values.
- Theorem2 applies to the labelled path-complement metagraph, not arbitrary arc-flip, switching, aggregate, black-support, or unmarked-isomorphism graphs.
- Equations(1)-(2) count edge orbits under specified actions. They do not reconstruct adjacency among unlabelled tournament vertices or identify any black/blue coloring.
- The signed edge test is executable from a single normalized Hamilton path: reflect its off-path bit word. This can be added to a metagraph atlas without enumerating the full signed group.
- The free-neighbour action gives an immediate audit gate for future wiring data: every Aut(T)-orbit of labelled neighbours must have size|Aut(T)|. A smaller orbit is a concrete decoding or labelling error.
- The SO(3) and field-dimension obstructions supply cheap hostile probes for any proposed geometric interpretation. A picture preserving only C7 must not be labelled a faithful realization of F21.

Run `python3 verify.py --check` and `python3 -O verify.py --check`. All finite gates remain active in optimized Python. The stored JSON records the enumerated universes and distinguishes the n8 formula row from exhaustive rows.
