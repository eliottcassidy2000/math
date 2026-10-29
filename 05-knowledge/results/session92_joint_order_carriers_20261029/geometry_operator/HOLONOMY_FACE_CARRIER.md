# A hidden sign can turn a product kernel into a face enumerator

**Session92, 2026-10-29. PROVED elementary construction; finite checks separated below.**
No external-priority claim. Source pin `b07ec51b204e`.

## Result

For every hypergraph \(\mathcal F\) on \([m]\), with distinct hyperedges of size
at least two, we construct a real symmetric signed single-bit-flip operator
\(B_{\mathcal F}\) on \(2m\) bits. Every matrix entry has exactly the same
absolute value as the ordinary free-flip operator \(\sum_i X_i\). The hard
constraint graph is always the same matching of \(m\) disjoint edges.
Nevertheless the leakage Gram satisfies

\[
\dim\ker D_{\mathcal F}
 =\#\{S\subseteq[m]: \text{no }F\in\mathcal F\text{ is contained in }S\},
\qquad \operatorname{rank}D_{\mathcal F}=3^m-\dim\ker D_{\mathcal F}.
\]

Thus every finite simplicial complex containing its named vertices occurs as
the support indexing the surviving kernel lines. For graph hyperedges this
is the independent-set complex, and the nullity is its independent-set count.
This is an explicit construction, not a tournament-count identity or a proof
of efficient enumeration. The exponential matrix has not made counting cheap.

One full hyperedge \([m]\), for \(m\ge2\), reverses just one undirected
feasible/forbidden coupling, changes the rank from \(3^m-2^m\) to
\(3^m-2^m+1\), and agrees with the untwisted model on **every proper
occupied-module sector**. Arbitrarily high-order hostile controls therefore
have a two-entry Hermitian matrix modification.

## Inheritance, corrections, and concept board

Anchor: exact assembly of local operators. Niche: signed boundary incidence.
Wildcard: arbitrary simplicial complexes as surviving dark sectors.
Board: guarded excursions / signed cycles / spectator ports / simplicial faces /
local-versus-joint observations / arithmetic square classes.

* Closest proved mechanism: Session91
  `operators/GUARDED_FLIP_LEAKAGE.md` identifies the ordinary leakage with
  signless incidence and proves nullity \(2^{b(H)}\), including arbitrary
  **named onsite weights**. We do not reproclaim that theorem: our amplitudes
  depend on the other occupied ports, and therefore are not onsite weights.
* Classical ancestor: Session85 `SIGNED_INCIDENCE_RELATIONS.md` and its
  Petersen suppression example distinguish balanced signed cycles, free
  parameters, and integral parity obstructions. We use the real kernel, not
  its integer Smith form.
* Older underused sidecar: THM-2353,
  `target-plaquette-holonomy-and-cross-axis-slice.md`, retains multiplicative
  plaquettes separately from additive interaction; THM-2814,
  `projective-allocation-square-holonomy-and-idempotent-provenance-no-go.md`,
  distinguishes normalized coefficients from projective line data. The
  present operator gives an actual realization of a nonflat square in a
  different model, **not** the physical LRC allocation those papers leave open.
* THM-2708, `c3-hermitian-gain-holonomy-discriminant-gate.md`, uses an F4
  holonomy to distinguish a triangle lift from a nine-cycle. That is an
  analogy in the mechanism, not a map from our real operator to that
  discriminant module.
* Corrected near miss: keeping all absolute amplitudes does not keep the
  return operator. Conversely, keeping the full named contextual boundary
  amplitudes does keep it. We claim neither identical full local responses
  under all spectator assignments nor a physical realization as ordinary
  independently driven PXP dynamics.

Primary background checked: Germina, Hameed and Zaslavsky,
[On Products and Line Graphs of Signed Graphs](https://people.math.binghamton.edu/zaslav/Tpapers/cpl.pdf),
Lemma2.2, gives signed-incidence/Laplacian rank as vertex count minus balanced
component count. Their Cartesian products retain factor-dependent signs;
ours deliberately permits spectator dependence and is not such a signed
Cartesian product. Reff,
[Spectral Properties of Complex Unit Gain Graphs](https://arxiv.org/abs/1110.4554),
provides the broader unit-gain incidence framework. These classical facts
are background, not new results claimed here; a direct proof follows.

## 1. Exact operator and ports

Each constraint edge is one module with two physical bits. Its three allowed
states are \(\emptyset,0,1\): no bit occupied, first port occupied, second
port occupied. Both occupied is forbidden. Let \(P\) project onto the \(3^m\)
feasible states, and put

\[
A=PBP,\quad R=(1-P)BP,\quad D=R^tR=PB^2P-A^2
\]

on the feasible space. Start all hypercube edge amplitudes at +1.
A forbidden state adjacent to feasibility has exactly one double-occupied
module i; every other module is empty or occupies one port. It has precisely
two feasible predecessors, called \(s_0,s_1\) according to that module's port.
Set its boundary row to

\[
(Rx)_Y=x_{s_0}+(-1)^{p_i(s)}x_{s_1},\qquad
p_i(s)=\sum_{F\in\mathcal F:\max F=i}
       \prod_{j\in F\setminus\{i\}}[s_j=1]\pmod2.
\]

These assignments never conflict: a physical hypercube edge has a unique
forbidden endpoint. Assign the same real amplitude to its reverse direction.
Leave all other hypercube edges +1. Hence B is real symmetric and has the
same absolute entries as the free flip operator. In particular A is unchanged:
all altered edges leave feasibility.

## 2. Face-enumeration theorem and proof

A row exchanges port0 and port1 without changing which modules are occupied.
Thus the boundary exchange graph splits into cubes \(Q_{|S|}\), one per
occupied-module subset S. Isolated vacuum is the zero-dimensional cube.
Each row imposes a nonzero ratio between its endpoints. On a connected cube
there is therefore at most one kernel dimension: choose a value at one
vertex and propagate it. A negative cycle product forces that value to zero.
If all signs are untwisted, \(x_s=(-1)^{\#\{j:s_j=1\}}\) gives a nonzero
kernel vector.

If S contains no hyperedge, every monomial vanishes identically on its cube,
so precisely this one kernel line survives.

Conversely suppose S contains a hyperedge. For some i, the polynomial p_i
restricted to S is nonzero. Its surviving monomials are distinct squarefree
nonconstant monomials; the Boolean multilinear monomials are linearly
independent as functions, so p_i is nonconstant. Some derivative
\(\partial_jp_i=p_i(x_j=0)+p_i(x_j=1)\), necessarily with j<i, is a
nonzero function. Choose a vertex where this derivative is1. The mixed
(i,j)-square has gain product

\[
(-1)^{\partial_jp_i+\partial_ip_j}=-1,
\]

because p_j depends only on coordinates less than j and thus never on i.
The four extra minus signs from converting sum equations to propagation
ratios cancel. The square kills the kernel on this entire connected cube.
Summing the surviving one-dimensional blocks proves the formula. Positivity
\(D=R^tR\) over R, or the same real matrix over C with adjoint, identifies
its kernel with that of R. This step is **not** exported to finite fields.

Empty and singleton hyperedges are excluded: the construction always retains
the vacuum and all one-module kernel lines. Empty hypergraphs and m=0 are
allowed; the formula gives nullity2^m and nullity1 respectively.

The argument does not require \(\mathcal F\) to be an antichain. Distinct
hyperedges cannot cancel as Boolean monomials. Redundant hyperedges leave the
kernel count unchanged but can change actual coefficients and the Gram.

## 3. Small hostile and all-order family

For m=2, \(\mathcal F=\{\{1,2\}\}\), both modules alone have the same
three-state operator as before. In the fully occupied sector order states
as (00,01,10,11). The untwisted block is

```
2 1 1 0
1 2 0 1
1 0 2 1
0 1 1 2
```

of rank3. Reversing the edge (10,11) gives

```
 2  1  1  0
 1  2  0  1
 1  0  2 -1
 0  1 -1  2
```

of rank4. The full feasible-space ranks are5 and6. The rank **gap** is one;
the matrix difference itself has rank two, not one. The diagonal and
entrywise absolute values of these Grams agree. A port-labelled gauge cannot
remove the negative square product.

For one hyperedge of size r=m, only the boundary row with module m doubled
and every earlier module at port1 changes. Thus every sector omitting even
one module is identical, but the full cube loses its kernel line. For m=3
the ranks are19 and20; all one- and two-module occupied restrictions agree.
This is the cheap arbitrary-order hostile generator requested by the broader
composition investigation.

For graph-valued twists, each pair square detects exactly its graph edge.
All such squares are flat iff the twist graph is empty. In the hypergraph
case the same square test needs its spectator assignment: omitting that
assignment loses the higher-arity monomial. It is not enough to store only
one square value per pair.

## 4. Typed bridges and limits

* **To tournament run lists:** shared shape is independent child responses
  versus a joint assembly obstruction. Our source is a signed boundary
  operator, target a cube gain graph, map the two-predecessor rows to gain
  edges; preserved predicate is zero second-order leakage. No actual map to
  the tournament Hamilton-path counts is known. Their missing joint order
  is not proved to equal our flux bit.
* **To graceful/square-sum incidence:** the literal shared object is a
  signed linear equation on graph endpoints. Sum versus difference changes
  cycle balance. Arithmetic realization guards, distinct labels, and
  nonnegative squares do not transfer under this linear correspondence.
* **To LCA/tree carriers:** a spanning tree fixes propagation from one
  amplitude, but every omitted chord needs its gain consistency test.
  The example is an exact warning for ports assembled along a tree: tree
  data alone do not certify later cycle closure. No LRC conclusion follows.
* **To quantum language:** D is a positive return Gram; its kernel comprises
  feasible amplitudes with no one-step escape under this specified B.
  Ordinary onsite drives factor differently. We claim no experimentally
  implemented interaction, QCD identification, or mutually unbiased basis
  construction.
* **To enumeration:** kernels enumerate faces because each admissible
  occupied subset contributes exactly one line. This is a useful formal
  translation for hostile testing, not an algorithmic speedup for independent
  sets or A000568. The classical signed-incidence law supplies the linear
  mechanism; the new local construction supplies the chosen face complex.

## 5. Reproduction

Run `python3 check_holonomy_carrier.py` and
`python3 -O check_holonomy_carrier.py`. The checker uses exact rational
elimination, exhausts all hypergraphs on at most three named modules, all64
simple graphs on four modules, and full-hyperedge controls on four and five
modules. It constructs physical hypercube operators independently for m2,3,
checks the projection identity entrywise, verifies unchanged magnitudes and
one changed undirected edge, and checks each sector rank. Normal and -O
checks use explicit exceptions, not removable assertions.

Final local receipt: **4,577 checks in each of normal and `-O` modes**,86
primary families, plus the displayed repeated small controls. Both exit0.
The physical controls also verify equality of the entire legal operator A
and rank2 of the Gram difference. Independent auditor scripts are separately
owned in `audit/check_flux_independent.py` and
`audit/check_hyperflux_independent.py`; their evidence is not included in this
producer's count.
