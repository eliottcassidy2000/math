# Hard guards, signed incidence, and the exact leakage rank

**Session91, 2026-10-29. PROVED below; finite checks separately scoped.**
No external-priority claim is made. This is a new application and synthesis
relative to the repository sources inspected, not a claim to have classified
all constrained quantum models.

## Inheritance and operation contract

The Session90 decorated Walsh theorem and the incoming
`odd_sign_spectra_trijections_20261028/THEORY.md` diagonalize **free** coordinate
flips using signed support orbits. Their count correspondences do not remove a
hard feasibility condition. THM-4078
`even-graph-triangle-quotient-spectrum-and-boolean-noncommutation.md` already
warns that changing an operation to its Boolean support changes its algebra;
its graphs are degree-even, not Royle-even. The signed-incidence mechanisms in
Session85 are a second ancestor of the proof below. These are inherited ideas,
not new names for the same counts.

Source pin: `90e704a8f35614d2085778d1b518f03ff0d750ca`.
Live concepts: free flip spectrum / hard feasibility / independent-set exchanges /
source weights / bipartite sign kernel / geometric constraint graphs.

Let H be a finite simple graph on n named coordinates. A bit string is feasible
if its support is an independent set of H. Write I(H) for the number of
independent sets, including the empty set. Let P be the diagonal projector
onto these states in the real vector space on all 2^n bit strings. Put

    B = sum_i X_i,       A = PBP restricted to the feasible states,
    L = PB(1-P),         D = LL^t.

X_i flips bit i. A is the actual legal-flip adjacency. L records every named
one-step passage from a feasible to an infeasible state. All entries initially
have unit weight. The operators are not transition probabilities.

The elementary block identity is

    PB²P − A² = PB(1-P)BP = D.                       (1)

Thus D counts precisely the two-step excursions that leave feasibility and
return. It is positive semidefinite over the reals. Replacing free evolution
by legal evolution loses this second-order term, even when a coding preserves
the number of states. This is a standard projection identity, not a new
identity of matrix algebra.

## 1. Exact rank theorem for every constraint graph

Let b(H) be the number of **bipartite connected components** of H, counting
isolated vertices as bipartite components. For the graph with no vertices,
b=0 and I(H)=1.

**Theorem.** Over Q or R,

    nullity D = nullity L^t = 2^b(H),
    rank D = rank L = I(H) − 2^b(H).                 (2)

This includes disconnected graphs, isolates and the empty graph.

### Step 1: every relevant forbidden state gives an edge or a pin

A forbidden support Y has a feasible predecessor Y\{v} exactly when deleting
v removes every edge induced by Y. If it has two such vertices u,v, every
induced edge must contain both, so its conflict graph is exactly the single
edge uv. There are then exactly two predecessors. If it has one predecessor,
all its conflict edges form a star with at least two edges and the deleted
vertex is its center. Any other forbidden Y has no feasible predecessor.

Construct an exchange graph J on feasible independent sets: connect S,T when
some forbidden Y has exactly those two predecessors. These edges exchange
one occupied vertex for an adjacent unoccupied vertex of H. Let p(S) count
forbidden Y having S as their unique predecessor. Then

    D = Q(J) + diag(p),                             (3)

where Q is the signless Laplacian. Equivalently,

    x^t D x = sum_{ST in E(J)} (x_S+x_T)²
                + sum_S p(S) x_S².                 (4)

The boundary incidence has exactly one column per such edge or pin. Hence its
kernel on an exchange component is one-dimensional precisely when that
component is bipartite and has no pin, and is zero otherwise. This is the
ordinary signed-incidence mechanism with positive edge signs.

### Step 2: a connected H has no higher-particle unpinned component

Assume H is connected. The empty independent set is an isolated unpinned
exchange state, contributing one kernel dimension. Singleton states induce
exactly H in J and have no pins. They contribute one further dimension iff
H is bipartite.

For an independent set S with at least two vertices, choose a closest pair
u,v in graph distance, with distance d>=2. If d=2, their common intermediate
vertex sees two occupied neighbors, and adding it gives a forbidden star
with unique predecessor S: S is pinned. If d>2, move u one step along a
shortest u-v path. The new vertex cannot be adjacent to any other occupied
vertex: otherwise its distance from the old u would give an occupied pair
at distance at most two, contradicting minimality of d. Thus the exchange is
legal and decreases the minimum occupied-pair distance. Repeating reaches a
pinned state. Every exchange component of particle number at least two
therefore contains a pin and has zero kernel. This proves (2) for connected H.

For connected K1, the empty and singleton states both survive, giving nullity2;
there are no forbidden states and no contradiction.

### Step 3: disconnected graphs factor exactly

For H=H1 disjoint-union ... disjoint-union Hc, feasible spaces tensor-factor.
A forbidden state reached by one flip has precisely one infeasible component.
A second flip can restore feasibility only in that same component. Therefore

    D_H = sum_j I tensor ... tensor D_Hj tensor ... tensor I.

For positive semidefinite matrices, the kernel of this sum is the tensor
product of the individual kernels. Its dimension is the product of the
connected nullities, which is 2^b(H). QED.

## 2. Nonzero source weights preserve rank, not necessarily the integer lattice

Replace B by B_w=sum_i w_i X_i, where every named weight w_i is a nonzero
real number. The projection identity still gives D_w=L_w L_w^t. In a
boundary equation with forbidden Y, substitute

    x_S = (product_{i in S} w_i) y_S.

Every predecessor coefficient becomes the same nonzero factor
product_{i in Y} w_i. The equation reduces to y_S+y_T=0 or y_S=0, exactly as
before. Thus (2) holds for arbitrary nonzero real weights, including signed
ones. The empty product is one.

This rescaling need not be an integral unimodular change of basis. It preserves
rank and the real kernel dimension, not Smith invariants or a modular model
when weights cease to be units. Zero weights remove operations and are outside
this theorem's uniform contract.

## 3. The minimum second-order interface is almost the whole state space

Any real factorization D=YY^t needs at least rank(D) columns; a spectral
factorization attains this number. Thus an exact linear interface retaining
this two-step return Gram needs at least I(H)−2^b(H) scalar channels.
For connected nonbipartite H this is I(H)−1, despite the guard being described
by only the edges of H. This is a linear-algebraic rank bound, not a bound on
bit entropy, succinct formulas, decision certificates, or arbitrary nonlinear
encodings.

For a path with n>=1 vertices, I(H)=F_(n+2) and H is bipartite, so the rank is
F_(n+2)−2. For a cycle with n>=3 vertices, I(H)=L_n and the rank is

    L_n − 1 − [n is even].                          (5)

The Fibonacci/Lucas counts grow like powers of the golden ratio. A two-state
transfer matrix counts these feasible words, but cannot thereby serve as this
almost-full-rank return interface. The two tasks ask different questions.

Free Walsh diagonalization is still an exact change of basis on the ambient
cube. What fails is dropping P. Under Fourier transform, a hard guard becomes
an additional convolution operator on support characters. It is not in general
diagonal in the free-flip basis. Count equality, signed-support positivity,
and the edge budget e=v−1 do not discharge this extra operator.

## 4. A genuine quantum-adjacent interpretation, with geometry kept separate

The legal operator A is the standard independent-set/PXP flip Hamiltonian:

    A = sum_i (product_{j adjacent to i} P_j^0) X_i

on the feasible subspace, where P_j^0 projects a neighbor to bit0. The graph
H is a **blockade/constraint graph**, not the state graph of A.
Li, Giudici and Pichler, [Phys. Rev. Research6,023146(2024)](https://doi.org/10.1103/PhysRevResearch.6.023146),
SectionII, describes the constrained Hilbert space and projection viewpoint.
Stastny, Büchler and Lang, [Phys. Rev. B108,085138(2023)](https://doi.org/10.1103/PhysRevB.108.085138),
studies realization of local constraints in planar Rydberg structures. These
primary papers establish the physical modeling connection; they are not cited
as sources of the specific rank theorem (2).

A geometric realization of H can require distance inequalities and auxiliary
sites. It is not automatically a unit-distance graph, a square-contact graph,
or the physical eleven-square packing. Nor is this PXP model a QCD model.
The exact map here is from a constraint graph to its feasible-bit Hilbert
space and guarded flip operator; geometry is a separate input obligation.

## 5. Verification

`check_guard_leakage.py` constructs both the literal boundary matrix and the
independent full-cube two-step moment formula. It tests every labelled simple
graph through five vertices,72 deterministic graphs at orders6–8, and cycles
through10. Rational elimination verifies the ranks. The source-weight gauge
is tested on every graph through4, including negative weights. The next note
records exact golden and modular controls at cycle length5.

The final checker performs 251,825 active checks. Normal and optimized Python outputs match. All checks remain active
under `-O`; only the standard library is used. The all-order claim follows
from the proof above, not extrapolation from the finite census.
