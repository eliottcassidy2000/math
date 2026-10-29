# Eleven squares, a golden clock, and the metric boundary at 7:3

Session91, 2026-10-29. **PROVED:** the algebraic, flat-torus and finite-metric statements below. **CITED / INHERITED:** the packing facts and golden contact-family dictionary. **FINITE-EXACT:** the accompanying controls. No external-priority claim or major open-problem closure is made.

## 1. Recover the packing problem before transferring its representation

The physical problem packs eleven congruent closed unit squares inside a square container. Independent rotations are allowed; boundary contact is legal, while open interiors must be disjoint. The [primary formalization repository](https://github.com/Queuingtheorydotcom/11SquaresFormalized), inspected on 2026-10-29, reports a completed proof using Lean's kernel and native numerical certificates, rather than kernel-only checking. It gives

    T = (6u+4)/(1+2u-u²) = 3.8770835900228141773...,
    5u⁸−10u⁷−2u⁶+14u⁵+12u⁴−6u³+2u²+2u−1=0,
    9/25 < u < 37/100.

This session does **not** rebuild that formalization. Its global optimality is cited with the source's trust model, not established by our small arithmetic tests. [Stromquist's 2003 paper](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v10i1r8/0) is the older rigorous source for a lower bound and for the need to leave the 0/45-degree orientation restriction. Its historical frontier must not be mistaken for the present source report.

Two distinct prior repository records matter:

- `05-knowledge/results/square11_octic_langlands_collatz_20261006.md` derives the octic field and its S8 Galois group, including nonconstructibility. Those results are inherited, not recomputed here.
- `05-knowledge/results/eleven_squares_routing_bridges_20261007.md` treats Pingyou LTD's supplied *The Fibonacci Geometry of Eleven Squares* as an accepted premise. It describes a two-parameter Trump **contact family**, not every physical packing. Its marked-torus representation preserves two signed clearances and the enclosing side. For fixed shape parameter, forgetting the mark loses the side/clearance information. It explicitly does not assign one individual planar square to each torus cell. The original PDF was not available in this checkout and was not independently reread here.

The source snapshot is `90e704a8f35614d2085778d1b518f03ff0d750ca`. Relevant arithmetic predecessors are `THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md`, with its Session74 all-period correction, and the Session90 exact carrier/quotient distinctions. The live board is contact inequalities / finite quotient / actual metric / odd affine symmetry / shortest-vector choices. The underused coordinate is the metric or chosen displacement lift, not another cardinality match.

## 2. The golden eleven-cell quotient and its nonperiodic lift

Use

    Q = [[0,1],[1,1]],       Q²=Q+I,
    C = 4I−Q = [[4,−1],[−1,3]],       det C=11,
    M_w=I+wQ,               B_w=M_w C.

For any real w with 1+w−w² nonzero, define the flat torus R²/B_w Z² and its eleven-point subgroup M_w Z²/B_w Z². The arithmetic statements below make no assertion that an arbitrary such w belongs to the physical packing window.

The quotient Z²/C Z² is identified with F11 by

    (a,b) -> a+4b mod11.

Both columns of C map to zero and its index is11, proving the identification. Q acts by multiplication by4, which has order5. This is precisely the arithmetic period-five cell action recovered in THM-4592; the common language with cyclic Zeckendorf words is inherited.

But the lift is not periodic:

    Q⁵ = [[3,5],[5,8]],
    Q⁵−I = Q³ C,            Q⁵ != I.

The first equation explains why the fifth iterate fixes all eleven cells; it does not make the ambient torus map or a physical packing return to the same metric configuration. The integer quotient has discarded the displacement Q³Cz.

The shape formula is equally explicit:

    B_w = (4−w)I+(3w−1)Q,
    det B_w = 11(1+w−w²).

At w=1/3, B_w=(11/3)I. The inherited contact-family source says that this square-torus parameter lies outside its physical parameter window. We use it only as an algebraic control, not a feasible packing. The eleven-point index stays11 while the shape and metric change with w.

## 3. A sharp metric obstruction, and what the clock really preserves

**Proposition 1.** No nonzero real symmetric form G satisfies QᵗGQ=G. In particular the golden step is not an isometry for any positive-definite flat metric.

Write G=[[a,b],[b,c]]. Then QᵗGQ=[[c,b+c],[b+c,a+2b+c]]. Equality forces c=a, c=0, and2b=0, hence G=0.

The correct arithmetic norm has a different signature:

    N = [[2,1],[1,−2]],       Qᵗ N Q = −N.

Thus Q² preserves the split quadratic form 2(a²+ab−b²), while Q reverses its sign. This is a Lorentz-type arithmetic norm, not a Euclidean squared distance or square-packing clearance.

**Proposition 2.** For every k>=1 and every fixed positive-definite metric, the ratio of the largest to smallest singular values of Q^k is at least phi^(2k). This lower bound is attained in the ordinary Euclidean metric.

The eigenvalues are phi^k and (−phi^−1)^k. The largest singular value is at least the largest eigenvalue magnitude, and applying that statement to the inverse bounds the smallest singular value above by phi^−k. Their ratio is at least phi^(2k). Since Q is real symmetric, the ordinary metric attains equality.

In particular a five-step identity on the eleven cells can hide distortion at least phi^10 in the ambient lift. This is an exact model of a general warning: periodic symbolic coordinates do not give periodic geometry when the lift has nontrivial integer displacement or expanding directions.

Neither proposition contradicts the contact paper's quarter-turn statement. That operation moves the frame, marked lift and clearance functionals together. It is a different operation from holding the Euclidean metric fixed and applying Q. The result identifies the required sidecars; it does not refute the supplied packing dictionary.

## 4. An exact crystallographic distinction between 7:3 and 11:5

For a prime p≡3 mod4, let

    G_p = {x -> ax+b on F_p : a is a nonzero square, b in F_p}.

It has odd order p(p−1)/2 and preserves the Paley tournament. The arithmetic eleven-cell action contains G_11=11:5, whereas the familiar Paley7 group is G_7=7:3 of order21.

**Theorem 3.** G_p admits a faithful action by isometries of some flat two-torus if and only if p=3 or p=7.

**Necessity.** An isometry of a flat torus has a translation part and a lattice-preserving orthogonal linear part. For a finite odd subgroup, the linear image is a cyclic rotation group. Its generator has integer trace2cos(theta), so its odd order is1 or3. The translation kernel is abelian and normal.

For the nonabelian Frobenius group C_p:C_q with q=(p−1)/2>1, every abelian normal subgroup lies in C_p. Indeed, if such a subgroup contains an affine map with multiplier a!=1, commutators with translations generate C_p inside that subgroup; this contradicts commutativity with the multiplier. Therefore its quotient by the translation kernel has order at least q. The linear image has order at most3, so q<=3, giving p=3 or7. The p=3 group is already cyclic and acts by translations.

**Construction at7.** On R²/Z² use metric

    G = [[2,−1],[−1,2]],
    R = [[0,−1],[1,−1]],       v=(3,1)/7.

Then R³=I, RᵗGR=G, and Rv≡2v modulo Z². Translation by v and R generate the faithful group C7:C3. On the seven points xv, this is exactly x->2^j x+b, whose multiplier subgroup is {1,2,4}. It preserves the Paley orientation x->y iff y−x is a quadratic residue. QED.

This is an exact geometric reason to distinguish the pairs7/21 and11/55. It is **not** the theorem excluding7 and21 as Hamiltonian-path counts. Those exclusions use a different statistic and proof.

A simpler corollary already rules out the eleven-cell clock as a Euclidean symmetry: multiplication by4 fixes the zero cell and has order5 on the others. A torus isometry fixing a point has linear order1,2,3,4 or6, so its restriction cannot have order5. Finite arithmetic permutations need not be realizable by isometries of the geometric carrier in which their labels were discovered.

## 5. A toroidal unit-distance K7, and the exact lost displacement

The construction above has an additional strong metric property. Its seven points are pairwise equidistant in the **geodesic torus metric**.

For a nonzero residue d modulo7, a displacement lift has coordinates (a,b)/7 with (a,b)≡(3d,d) modulo7. Its squared length is

    2(a²−ab+b²)/49.

The numerator norm is positive and is0 modulo7, hence at least7. For d=1,...,6, the representatives

    (3,1), (−1,2), (2,3), (−2,−3), (1,−2), (−3,−1)

attain7. Every nonzero pair difference consequently has squared torus distance2/7. Scaling lengths by sqrt(7/2) produces a unit-distance K7 on a flat torus.

This does **not** produce a planar unit-distance K7 or a new lower bound on the chromatic number of the plane. At most three distinct points in the Euclidean plane are pairwise equidistant: place two at (0,0),(1,0); a third must lie at (1/2,±sqrt(3)/2), and the two alternatives are distance sqrt(3) apart. No fourth can be added.

The quotient lost an operation-relevant coordinate. Each pair is permitted to choose its own shortest lattice displacement. Those choices cannot be integrated into a single consistent set of planar vertex lifts with all pair distances retained. A concrete triangle already witnesses the problem: shortest displacements for residues1,1,−2 can be chosen as (3,1),(3,1),(1,−2), divided by7; they sum to (1,0), not zero. Their nonzero lattice sum is invisible modulo the torus lattice.

This is the geometric analogue of retained port alignment in a graph codec, but now there is an explicit map and an explicit cycle defect. A necessary condition for edge displacements to come from planar vertex positions is zero sum around every cycle; on a connected graph this condition is also sufficient, by integrating along a spanning tree. That is exactly the same tree-versus-cycle mechanism as the previous session's trijections, with values in R² instead of F2. Preserving individual edge lengths does not pay this global integration condition.

## 6. Fourteen faces, two Fano planes, and two Hamming codes

This is the standard triangular torus mechanism, not a new realization of K7; see Richard Schwartz, [Vertex-Minimal Paper Tori, section 2](https://www.math.brown.edu/reschwar/Papers/PaperTori.pdf), for the classical seven-vertex,21-edge,14-triangle pattern and the distinction between intrinsic flat and extrinsic polyhedral realizations. Set

    e=(3,1)/7, f=(-1,2)/7, L'=Ze+Zf.

Both vectors and e+f have squared metric length2/7. The lattice L' contains Z² because (1,0)=2e−f and (0,1)=e+3f. Its covolume is1/7, so L'/Z² has7points. The six shortest lattice directions ±e,±f,±(e+f) have quotient labels ±1,±2,±3. Thus quotienting the ordinary triangular tiling gives every edge of K7 exactly once.

The14 oriented faces, with all arithmetic modulo7, are

    (r,r+1,r+3), (r,r+3,r+2),   r=0,...,6.

The first triangle lifts to (0,e,e+f), the second to (0,e+f,f), translated by a representative of r. All lifted faces have the same orientation. Every edge has two incident faces with opposite edge orientations, and every vertex link is a6-cycle. Hence this is a closed orientable surface; its Euler characteristic is7−21+14=0, so it is a torus. This occurrence of7and21 is therefore a genus-one cell count, not a Hamiltonian-path statistic. It is an intrinsic flat-torus construction, not a planar drawing or a claim about an isometric polyhedral embedding in R³.

Forget face orientation and divide faces into the two translation orbits

    A_r={r,r+1,r+3}, B_r={r,r+2,r+3}.

Each orbit is a Steiner triple system on seven labels: the six ordered nonzero differences in either base block exhaust F7*. Consequently every unordered vertex pair occurs exactly once in each family. Each is a Fano plane. In the Paley7 orientation, the first oriented faces have differences1,2,4, and the second have differences3,6,5; both families are directed cyclic triples, with opposite geometric orientations. There are exactly14 cyclic triples in this regular tournament: among its35 triples, the21 transitive ones are counted uniquely by their source, giving7*binomial(3,2)=21. Thus these faces are exactly its cyclic triples, not merely selected examples.

The F2 row spans of the two7-by7 face-incidence matrices are [7,4,3] Hamming codes. This finite claim is checked completely by enumerating each span: each contains16words, with weight enumerator1+7z³+7z⁴+z⁷. Their intersection is precisely the zero and all-one words. Accordingly this particular geometry produces two *distinct embedded* Hamming codes rather than one canonical code after labels are forgotten. A change of labels identifies their abstract code types but loses which Fano system a face belongs to.

This supplies an explicit bridge to the repository's square-class/Hamming construction: both use seven-point binary incidence with the same abstract code type. It does not identify the arithmetic primes with torus labels without an additional incidence-preserving bijection, nor turn Fano matroids into ordinary graph cycle matroids. Indeed, graph incidence, triangular face incidence, and arithmetic square-class incidence are different linear operators. The useful common object is the specified binary linear map, not an unlabeled numerical seven.

## 7. Verification and limitations

`check_golden_metric.py` uses only exact integers and fractions. It verifies the quotient labels, Q⁵−I factorization, lattice/Gram formulas on a rational grid, nonisometry, the explicit21-element isometry action, its Paley orientation, the seven torus distances, oriented face gluing, six-cycle links, Fano pair incidence, Hamming spans and their exact intersection. There are10,296 active checks. All checks remain active under Python optimization.

The all-order group classification and distortion bound are proofs above, not extrapolations of matrix searches. The script does not rebuild the packing proof, classify all eleven-square contact configurations, or solve LRC(14), Hadwiger–Nelson, Collatz, or an isomorphism problem. The concrete gain is a metrically typed bridge: golden arithmetic clocks, flat-torus symmetries and planar geometric lifts have different exact symmetry and cycle constraints.
