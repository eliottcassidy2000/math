# Independent Boolean-sector audit and the exact two-state boundary repair

**PASS: analytic proof accepted and exact controls independently rebuilt.**
The [producer theorem](overnight10_20260906_boolean_tokens.md)
correctly identifies the native induced Boolean triangle sector, its vertex
budget, full-rectangle adjacency spectrum, and extra resonant zero modes.
No repair to those statements was required. A separate proof read by the
producer also accepts the root referee's pendant-pair cap extension below,
including its integral basis change and kernel lift. The result concerns induced
adjacency, not the full ambient operator, its Laplacian, or directed toggle
multiplicities. The separate extension below is PROVED ELEMENTARY, with
finite-exact controls; it preserves the integral cokernel at one resource cap.

## Proof audit

The four possible numbers of old edges in a toggled triangle exhaust the
maximum-degree-two source and target condition. At zero or three old edges
the cycle number changes; at one or two it inserts or deletes a vertex on
one cycle. In particular a toggle meeting two distinct cycles cannot survive
the degree condition. Sorting the cycle lengths, including ties, yields
exactly the path-token moves after adding the coordinate index. The isolated
vertex budget becomes the displayed sum cap and is necessary as well as
sufficient. The coordinatewise-minimum path proves the distance formula
inside that cap.

On an ordered path, a legal one-particle move does not cross another token.
Its exterior sign is positive. Thus the additive compound, not the
multiplicative exterior matrix, is the Boolean adjacency. Wedges of the
ordinary sine eigenbasis give the claimed sum spectrum and determinant
eigenvectors with their multiplicities. The c=N and M=3 masks are sound.
For even N and odd c, the signed subset enumerator proves equal parity
classes. The three-cosine identity supplies distinct orthogonal zero modes
in the claimed infinite family. Equality of nullity is claimed only in the
explicit 56-state case, where two independent exact paths verify it.

The representative multiplicities count actual triangle subsets rather
than ordered vertex triples: growth has m_s*s*I choices and shortening has
m_s*s choices. Their detailed balance agrees with the orbit volumes.
The Laplacian loses twice the number of occupied path edges from its diagonal;
further restriction loses the external boundary degree. Those changes prevent
the proposed adjacency spectrum from being silently transferred to either
other operator. The longer-cycle and ambient-boundary hostiles are valid.

The primary sources were read directly during this audit. The definition
of token graphs agrees with [Fabila-Monroy et al., Token Graphs](https://arxiv.org/abs/0910.4774).
The additive exterior construction and positive path signs are classical;
see [Mallory--Raz--Tamon--Zaslavsky, Section 3 and Fact 13](https://people.math.binghamton.edu/zaslav/Tpapers/xbal.pdf),
printed pp. 4-6 and 9. Neither source is cited for the repository's native
cycle-sector identification or for the cap extension below.

## A cap preserving the complete integral cokernel

Use the producer notation c>=1, M>=4, d=M-3, with cd>=2. Let A_full be the
adjacency of the full sector S(cM;c,M), and A_cap that of S(cM-2;c,M).
The condition cd>=2 ensures cM-2>=3c and a nonempty capped sector.

**Theorem.** There is an integer unimodular congruence

```
A_full ~ [[0,1],[1,0]] direct-sum A_cap.             (1)
```

Consequently A_full has exactly two more nonzero Smith factors, both units,
than A_cap; their integer cokernels and adjacency nullities are isomorphic.
Their real inertia differs by one positive and one negative eigenvalue.
This is not equality of their nonzero eigenvalues.

**Proof.** The only removed states are

```
u=(M,...,M),       v=(M-1,M,...,M).
```

The state u is a leaf and v its unique neighbor. Ordering u,v before the
retained states puts the Boolean adjacency in the exact form

```
A_full = [[0,1,0], [1,0,b^T], [0,b,A_cap]].          (2)
```

Replace every retained basis vector e_i by e_i-b_i e_u, leaving e_u,e_v
unchanged. This change is integral and unitriangular, so its determinant is
one. Direct multiplication turns (2) into (1): the cross term at v cancels,
while e_u has no retained adjacency or diagonal term and hence induces no
correction on A_cap. Smith equivalence and real inertia follow immediately.
In particular every kernel vector z of A_cap lifts uniquely to
`(-b^T z,0,z)` in the full ordering. The converse is obtained from the u row,
which forces the v coordinate to zero. This also proves nullity equality
without any spectral theorem or numerical rank inference.

This is an application of elementary pendant-pair elimination, whose full
matrix calculation is included; no priority is claimed for that algebraic
operation. The recovered coordinate is its exact location at the two-state
resource boundary of this native induced sector.

**Corollary.** For odd h>=3, c=3, M=3h-1, the capped native sector at
`n=3M-2=9h-5` has balanced parity classes and nullity at least h-1. The two
removed states have opposite edge parity, so balance survives, while (1)
preserves the producer's nullity. At h=3 this gives n22 with 54 vertices,
parity 27+27, and nullity exactly two. Full rectangles at n>=3M remain the
producer's separate theorem; the intervening one-state cap need not preserve
nullity. At n23 in this example it is one.

The cap map preserves the integral cokernel and kernel via the displayed
lift but discards the two nonzero spectral directions. It neither extends
vectors to the ambient Eulerian graph nor gives a new Hasse-jet Smith law.
The resemblance to integral-observer calculations is a method connection,
not a transfer between the matrices.

## Independent exact universe

The verifier imports no producer or repository code. Its literal graph path
constructs edge sets and triangle symmetric differences, then checks degree
and connected components. Fresh complete sectors/caps (n,c,M) are
(9,2,6),(11,3,5),(15,3,5),(15,4,6), plus the critical (24,3,8).
In total 122,243 actual triangle toggles are tested. Every retained ordered
pair is checked for Boolean adjacency, exact directed multiplicity, and
orbit-volume detailed balance. Breadth-first distances are compared with
the sorted-coordinate formula for all vertex pairs.

A separate exact cyclotomic engine checks 41 literal Slater vectors, including
every mode for (N,c)=(3,2),(4,2),(5,2),(6,3) and the two zero modes at (8,3).
It verifies the eigenvector equation in each coordinate, rather than only
matching characteristic polynomials. Integer matrix rank independently
confirms the critical nullity. The complete finite cap bank is

| n | vertices | parity index | nullity |
|---:|---:|---:|---:|
|19|45|3|3|
|20|49|1|1|
|21|52|2|2|
|22|54|0|2|
|23|55|1|1|
|24|56|0|2|

This bank is FINITE-EXACT and is not a formula for every cap. Fifteen fresh
rectangles 1<=c,d<=4 with cd>=2 verify (1) by full integer matrix multiplication
and determinant-one basis changes, including a 70-state control.

```
python -B 04-computation/overnight10_20260906_boolean_independent_audit.py
python -B -O 04-computation/overnight10_20260906_boolean_independent_audit.py
```

Both runs pass **17,463 always-active gates** with identical LF transcripts.

```
source f588933bec981dca1ac912192dfd0986f53f9953ed47827e22df8434b593ca70
output c5cfb1ba4da87f1dc69267e2fbffa540aede40fc18af45ec05f5e73c90ce5c23
```

All artifacts are outside the worktree until root integration. The full
ambient spectrum and the general capped-sector spectrum remain OPEN.

**Filing:** root integrated these audited artifacts in the tenth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
