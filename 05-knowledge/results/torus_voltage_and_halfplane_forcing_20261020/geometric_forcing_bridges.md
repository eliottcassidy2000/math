# Half-planes, toroidal quotients, and the missing winding coordinate

2026-10-20. **PROVED elementary / classical deduction:** G1--G4 below. **FINITE-EXACT:** 160 checks, including independent Hamiltonian-path enumeration and subset DP. **SOURCE CLAIM / NOT BUILD-VERIFIED:** the supplied September 2026 six-colour lower bound. No new Euclidean chromatic bound, six-chromatic point set, or novelty claim is made.

The productive bridge is concrete: the triangular lattice has a seven-state toroidal quotient equal to K7. Orient its same edges in two different ways and the quotients become the regular carousel tournament and the Paley tournament. Their lifts have radically different directed closed walks. The missing coordinate is the displacement of a lifted walk, not its degree, endpoint residue, or ordinary spectral trace.

## 1. Inheritance and concept board

- Closest proved mechanisms: THM-418, `01-canon/theorems/THM-418-lattice-chromatic-dichotomy-of-the-plane.md`, gives the three-colour ceiling for a single triangular lattice at every norm scale. The complete boundary-palette relation in `05-knowledge/results/boundary_palette_compression_20261011.md` already specifies the correct gluing object.
- Canonical hostile: a quotient has extra identifications. In session 9, the finite norm-plane PSD/cut hostile additionally showed that radial positive-definiteness does not recover an actual Boolean colouring law (`session9-plane-coloring/proofs/finite_norm_triangle_bridge.md`). Neither phenomenon becomes a Euclidean certificate by reinterpreting it.
- Corrected near miss: inherited `finite_algebraic_witness_20261007.md` and `plane_algebraic_witness_20261011.md` keep finite six-colour existence conditional on the supplied primary paper. Large three-colourable exteriors are not six-colour forcing: THM-4614 proves a marked-port extension theorem; THM-4661 keeps its stronger exterior bounds and equality-tester scopes separate.
- Least-used sidecar here: a lattice-valued voltage on each quotient edge. Its sum identifies exactly which quotient returns really close in the plane.

Live board: **half-plane exhaustion; amenable action; exact unit incidence; labelled palette relation; torus quotient; winding-sensitive spectrum**. Tournament orientation is intrinsic only after the three forward unit directions are declared; it is not an invented comparison among colours.

## 2. G1: every nonempty open or closed half-plane has the plane's chromatic number

Write U(X) for the graph on X with edges at Euclidean distance exactly one. For a half-plane H and every finite integer k,

    U(H) is k-colourable if and only if U(R²) is k-colourable.

Restriction proves the reverse direction. For the forward direction, any finite set F of plane points can be translated entirely into H. Restrict a proper colouring of H to that translate, then translate its labels back. Thus every finite subgraph of U(R²) is k-colourable. Graph compactness gives a colouring of the entire plane.

The same proof works for every X containing translates of every finite plane configuration, in particular every X containing disks of arbitrarily large radius. Scaling is not used, since it would change the unit-distance predicate. A fixed-width infinite strip does not satisfy this disk hypothesis.

This is a standard deduction from de Bruijn--Erdős compactness, not a new lower bound. It concerns unrestricted colourings. The argument does not preserve measurability through compactness and makes no assertion equating measurable chromatic numbers.

**What half-planes offer computationally.** Cut a finite candidate by x=0, assigning boundary vertices to one chosen side. Every unit edge crossing the cut has both endpoints in the width-one seam |x|≤1. It is therefore exact to retain each side's full possible labelled colourings on its seam vertices, then test whether some pair satisfies every cross-seam inequality. Colour names can be quotiented only if one retains the relative palette permutation. Counting seam colours, or knowing that both sides separately admit k colours, is insufficient. The two endpoints (−1/2,0),(1/2,0) already witness a missing cross-seam constraint.

This is a specialization of the repo's complete boundary-relation method. It is not a bounded-size separator theorem: the seam may contain arbitrarily many vertices.

## 3. G2: a three-colourable lattice has a seven-chromatic toroidal quotient

Let ω=exp(iπ/3) and identify the triangular lattice with Z² via (a,b)↦a+bω. Its norm form is

    Q(a,b)=a²+ab+b²,

and its six unit steps are

    S={±(1,0), ±(0,1), ±(1,−1)}.

The colouring c(a,b)=a−b modulo 3 changes across every step, while a unit triangle requires three colours. Thus its chromatic number is exactly three.

Now define

    q(a,b)=a+3b modulo 7.

The six steps have the six nonzero images modulo seven. Consequently the quotient graph on Z/7Z is K7, whose chromatic number is seven. The kernel Λ has basis (−3,1),(1,2), with determinant −7 and each basis vector of squared length seven. Thus this is a genuine rank-two periodic quotient, not a graph obtained by deleting vertices of the lattice.

The quotient map is a graph homomorphism, so its automatic chromatic inequality is χ(lattice)≤χ(quotient), not the reverse. A colouring of the quotient is exactly a Λ-periodic colouring of the lattice. The lattice's usual three-colouring is not Λ-periodic: its value changes along both displayed kernel generators. Therefore the quotient's non-five-colourability excludes one periodicity constraint, not arbitrary plane colourings. There is no section of the quotient that realizes all K7 edges as simultaneous unit distances among seven distinct plane points. Already K4 cannot be so realized: the two common unit neighbours of one unit edge are at distance √3, not one.

### Spectral form of the same error

The quotient adjacency matrix is J−I, with spectrum 6 once and −1 six times. The infinite lattice's Fourier symbol is

    λ(s,t)=2cos(s)+2cos(t)+2cos(s−t)
          =|1+exp(is)+exp(it)|²−3.

Its minimum is −3, attained at s=2π/3,t=−2π/3. The seven-state quotient samples only (s,t)=(2πj/7,6πj/7). It excludes those minimizing phases, leaving the quotient's nontrivial eigenvalues at −1. A large chromatic lower bound computed on the quotient is therefore compatible with a three-colourable lift. Its operator acts on a different space: periodic functions rather than all lattice modes.

This is the exact obstruction to using a torus closure or a finite residue spectrum as a lower bound for the original plane graph without a lift theorem.

## 4. G3: zero winding is the correct return predicate

For each directed quotient edge r→s choose its unique unit step d satisfying q(d)=s−r. A walk lifts uniquely once its initial lattice vertex is fixed. It returns to that vertex exactly when

    sum of its step vectors = (0,0).

A quotient return requires only that this sum lie in Λ. Forgetting the sum therefore confuses genuine closed walks with paths that end at another translate of their start.

Equivalently, choose representatives t(r)=(r,0). The edge voltage is

    ν(r,s)=t(r)+d−t(s) in Λ.

For a closed quotient walk, the sum of its voltages equals the sum of its actual unit steps. This expression is independent of representatives: changing them adds a coboundary, which telescopes around the closed walk. A formal adjacency matrix with Laurent-monomial voltage weights restores the information. The coefficient of the identity monomial in its nth trace counts precisely the zero-winding closed walks. Evaluating every voltage monomial at one gives the ordinary quotient trace and counts all windings instead.

### First exact hostile

The quotient triangle 0→1→2→0 uses steps (1,0),(1,0),(1,−1), with total (3,−1)≠0. It is a triangle only after a nonzero period has been identified with zero.

Of K7's 35 unordered triangles, precisely 14 have zero winding and 21 do not. Thus

    tr(A_K7³)=210,
    zero-winding part of that trace=84.

Per starting lattice vertex, the true triangular lattice has 12 closed walks of length three, while the quotient has 30. The exact script verifies this independently by displacement dynamic programming and matrix/spectral counts. Through length eight, the return pairs (lift,quotient) are

    n=0: (1,1); n=1: (0,0); n=2: (6,6); n=3: (12,30);
    n=4: (90,186); n=5: (360,1110); n=6: (2040,6666);
    n=7: (10080,39990); n=8: (54810,239946).

The first two moments conceal the error; the third exposes it. The repaired observable is a group-ring coefficient, not the unmarked trace.

## 5. G4: the same construction produces carousel and Paley tournaments

Keep the quotient q and orient one edge in each opposite pair of triangular-lattice directions.

### Carousel orientation

Take forward steps

    D_C={(1,0),(0,1),(−1,1)}.

Their residues are {1,3,2}, so the quotient is the regular carousel C7: r→s iff s−r∈{1,2,3} modulo seven. Every vertex has outdegree three and indegree three.

The lifted orientation is acyclic: the real-valued integer height h(a,b)=a+2b increases by 1,2,1 along its three steps. No positive-length directed walk closes in the lattice, even though the quotient has directed cycles. Its six length-three returns per vertex are all nonzero-winding artifacts.

### Paley orientation

Reverse only the (0,1) direction:

    D_P={(1,0),(0,−1),(−1,1)}.

Their residues are {1,4,2}, the nonzero quadratic residues modulo seven, so the quotient is the Paley tournament P7. Degrees are again three in and three out.

Here the three steps sum to zero. If a length-n walk uses them x,y,z times, zero displacement says x=z and y=z. Hence n=3m and the exact number of closed directed walks based at one lattice vertex is

    (3m)!/(m!)³ if n=3m, and zero otherwise.

The quotient agrees at length three (six returns), but not at length four (12 quotient returns, none in the lift). Changing one direction therefore changes not just a label: it changes whether the local triangular relation survives before quotienting.

### Checked finite tournament data

Full permutation enumeration and an independent subset dynamic program both give:

- Carousel C7: 175 Hamiltonian directed paths; 17 Hamiltonian directed cycles modulo starting vertex.
- Paley P7: 189 paths; 24 cycles modulo starting vertex.

Their adjacency characteristic polynomials, calculated exactly from traces by Newton identities, are respectively

    t⁷−14t⁴−28t³−28t²−14t−3,
    t⁷−14t⁴−21t³−42t²−28t−24.

The Paley polynomial factors as (t−3)(t²+t+2)³. If α,β are the roots of t²+t+2 and u_n=αⁿ+βⁿ, then

    u_0=2, u_1=−1, u_n=−u_(n−1)−2u_(n−2),
    quotient returns per vertex at length n=(3ⁿ+3u_n)/7.

The division by seven is valid because the Cayley tournament is vertex-transitive, so every diagonal entry of Aⁿ equals its trace divided by seven. The recurrence and spectral factorization are independently checked against displacement DP through length 12. The initial return counts are 1,0,0,6,12,30,108. This finite recurrence forgets the lattice winding; the zero-winding lift instead has the multinomial formula proved above.

These are finite-quotient statements. In particular a Hamiltonian cycle of the carousel quotient is not a directed cycle of its acyclic lattice lift. The underlying quotient graph is K7 in both cases, so its undirected spectrum cannot distinguish these two orientation choices. The declared direction gauge is part of the object.

The accompanying [Euler-torus packet](../euler_torus_hamilton_transport_20261020/EULER_TORUS_BRIDGES.md) supplies the separate Euler-tour homology calculation. Euler tours, Hamiltonian paths, and embedded knots are not identified here.

## 6. Banach--Tarski supplies a boundary, not a forcing mechanism

The orientation-preserving plane isometry group is R² semidirect SO(2), an extension of abelian groups; adjoining reflections is a finite extension. As an abstract discrete group it is amenable. This supports invariant means and, for countable isometry subgroups acting on compact colouring spaces, invariant probability laws. By contrast the three-dimensional rotation group contains the nonamenable subgroup used in Banach--Tarski constructions. The plane does not admit the corresponding bounded-area duplication by finitely many rigid motions. Nonmeasurable subsets or the axiom of choice alone do not supply that missing action.

A self-contained [classic lattice nonduplication proof](CLASSIC_LATTICE_NONDUPLICATION.md) makes this boundary concrete for square-lattice translations and triangular-lattice isometries. It allows arbitrary pieces and uses quadratic box interiors versus linear boundary defects; it does not infer a universal plane or three-dimensional theorem from the lattice cases.

Amenability itself implies no useful chromatic ceiling. The Cayley graph on Z with steps ±1,...,±r has chromatic number r+1: {0,...,r} is a clique, and reduction modulo r+1 gives a colouring. This is an abstract graph example, not a unit-distance realization in the plane. Conversely the usual Cayley graph of a rank-two free group is a bipartite tree. Neither high chromatic number nor a large quotient spectrum diagnoses nonamenability.

The supplied six-colour paper uses amenability at a specific, legitimate interface: create a stationary law on proper labelings of the algebraic plane. It then needs substantially more than amenability to turn that law into ordinary measurable fields. Its central Haar-rigidity theorem and its geometric measurable obstruction are the load-bearing additional assertions; see `source_audit.md`. Averaging alone does not discharge them.

## 7. The exact bridge and its loss ledger

| Source | Map / target | Preserved | Lost | Required sidecar / decisive test |
|---|---|---|---|---|
| Plane unit graph | Restrict to a half-plane | Finite colourability by translation plus compactness | No constructive global colouring and no measurable regularity | Exact finite unit-distance witness for a lower bound |
| Triangular lattice | q to seven-state torus | Each local step; degree six | Distinct translates become one vertex | Λ-valued edge voltage; reject nonzero-winding returns |
| Three forward directions | Tournament quotient | Pairwise orientation after a declared gauge | Acyclicity and Euclidean closure | Height or zero-displacement test |
| Voltage adjacency | Set all monomials to one | Total quotient walk counts | Winding distribution | Extract identity coefficient before interpreting plane returns |
| Two cut fragments | Seam palette gluing | All proper colourings if complete labelled relations retained | Interior labels deliberately forgotten | Complete relation, exact cross-edges, compatible relative palette permutation |
| Arbitrary colouring | Amenable stationary law | All finite pointwise unit constraints in the law | Usual-topology strong continuity | Actual analytic transfer, not PSD or symmetry alone |

## 8. Productive next target and stopping boundary

For six-colour forcing beyond a triangular lattice, first require an exact finite Euclidean realization that escapes the inherited three-colour carrier. At palette five the inherited circle-saturation criterion is the right finite target: a five-colourable fragment whose every colouring uses all five labels on the neighbours of a prospective centre. Adding that centre then forces a sixth colour. The complete boundary relation and its omitted-state proof are essential; pairwise palette counts are not enough.

Any proposed torus/residue candidate should first pass the cheap voltage test developed here. If its apparent contradiction disappears on the actual lift, do not enlarge the same periodic search and report the extra quotient constraints as progress toward a Euclidean lower bound. The K7 example is an exact hostile control for that test.

This packet stops with completed bridge proofs and checked finite controls. It does not start an unbounded search, assert a new six-colour graph, or report the supplied Lean theorem as independently verified.

## Reproduction and sources

    python3 check_geometric_forcing.py
    python3 -O check_geometric_forcing.py

Both runs agree; `exact_checks.json` contains the full triangle census, quotient spectra, directed moment controls through length 12, and the tested finite patches. Checks use explicit exceptions rather than assertions, so optimized Python does not disable them.

- N. G. de Bruijn and P. Erdős, *A colour problem for infinite graphs and a problem in the theory of relations* (1951), Theorem 1: https://research.tue.nl/en/publications/a-colour-problem-for-infinite-graphs-and-a-problem-in-the-theory-/
- Stan Wagon, *The Banach--Tarski Paradox*, Chapter 10, invariant measures and the plane-isometry obstruction: https://www.cambridge.org/core/books/abs/banachtarski-paradox/measures-in-groups/6F9AE5EC8448E16EEFBD9B9796F4FADF
- Primary paper and exact source audit: `source_audit.md`, with pinned URLs and source blob hashes in `source_manifest.json`.
