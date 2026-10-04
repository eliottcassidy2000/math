# Five, six, seven: surfaces, retained memory, and certified Collatz observations

2026-10-04. **PROVED** elementary constructions and obstructions in the four
linked proof notes; **FINITE-EXACT** for their stated exhaustive universes;
**CITED** for the external graph-flow results. Collatz, Tutte's 5-flow
conjecture, and LRC(14) remain **OPEN**. No historical priority claim is made.

The strongest connection recovered in this session is that a small extra
coordinate can change which theorem is available. A local rotation turns a
graph into a surface map. A circulation completes a triangle-to-star flow
quotient. A marked axis completes a Collatz multiplier. The coordinates are
different, but in each case we can specify the exact forgotten information,
exhibit a hostile pair, and restore a faithful representation.

## 1. Inheritance, portfolio, and the six live objects

**Anchor:** Collatz source-aware route certificates and their observation
boundary. **Niche:** Paley/Fano maps with their rotations and flags.
**Wildcard:** Tutte flows and the surface-coloring obstruction. The live board
throughout is **guarded word / marked axis / pair orientation / rotation /
circulation and homology / finite observation**.

The closest proved mechanism is the exact affine source displacement
in [the mediant-tree proof](collatz_mediant_tree_K_20260930.md), joined to
the incoming [all-length boundary compiler](collatz_boundary_compiler_20261004.md)
and [completed inverse-route codec](inverse_ray_ternary_addresses_20261004.md).
The canonical hostile is loss of source or carry: a contracting formal map
can still increase its least legal positive source, as in the compiler's
`165 ->167` example. Its earlier `165 ->31` descent is retained; the example
does not refute the still-open first-stopping coincidence conjecture.

The corrected geometric near miss is in
[the September affine audit](collatz_blueprint_20260921_affine.md): the
unrestricted Collatz branch group fixes infinity and is not a discrete
triangle group. The least-used sidecars recovered here are local star
rotations, the projective sign sheet, seven triangle circulations, and the
word's finite fixed point. Relevant prior routes include:

* [THM-4529, Petersen/Heawood in Paley coordinates](../../01-canon/theorems/THM-4529-petersen-and-heawood-families-in-paley-coordinates-anti-circulant-all-odd-tournaments.md): actual triangle/star constructions, with the failed Hamiltonian-path parity transfer already removed;
* [THM-2996, prime-modular affine-defect trichotomy](../../01-canon/theorems/THM-2996-prime-modular-affine-defect-trichotomy-and-spherical-quartic-uniqueness.md): marked triangle-group quotients, not identification of a finite quotient with its universal geometry;
* [THM-4464, projective tournament and output sheet](../../01-canon/theorems/THM-4464-checksum-projective-tournament-and-output-sheet-cocycle.md): determinant orientations and a missing common lift sign;
* [noble subdivision with flags](noble_subdivision_flags_20261004.md): vertex/edge/face incidence and the distinction between cell and flag transitivity;
* [the Paley bridge audit](collatz_paley_bridge_20261001.md): a period-three coding is not a global symmetry of the integer Collatz function.

These sources were recovered before deriving the new interfaces. None of
their established mechanisms is claimed as newly discovered here.

## 2. What five, six, and seven actually measure

For a triangular tessellation with q triangles at each vertex, or the dual
cubic map with q-sided faces, the triangle excess is

\[
 \epsilon(q)=\frac12+\frac13+\frac1q-1.
\]

| q | excess | intrinsic geometry |
|---:|---:|---|
|5|1/30|spherical|
|6|0|Euclidean|
|7|-1/42|hyperbolic|

Thus the ordering is **5 spherical, 6 flat, 7 hyperbolic**. This is a theorem
about the specified incidence type, rather than a classification of every
object carrying one of those numbers.

The owner's three root-ray exponents `5+6b,3+6b,7+6b` retain their earlier
meaning for `H(n)=(3n+1)/2`, in source order `3,5,1 mod6`. The corresponding
fully accelerated valuations are `6+6b,4+6b,8+6b`; every one of those
single-letter forward maps contracts. Root deletion explains the offset in
the third ray. The first geometric number5 occurs in the first ray, but
its next exponent11 has negative triangle excess, while remaining on the
same arithmetic ray. No curvature classification follows from the residue
row. See [the sixth-clock convention proof](sixth_clock_branches_20261004.md)
and the [drift note](geometry_collatz_drift_carriers_20261004.md).

## 3. A canonical Paley construction changes surface after one bit flip

Use seven point labels P_p and seven line labels L_t, with incidence
`p-t in {1,2,4} mod7`. This is the Heawood graph. Its local reference orders
come from the cyclic list `(1,2,4)`. Each of the point and line parts has one
common rotation sign in a Paley-equivariant map.

**PROVED:** equal signs give seven hexagonal faces on a torus; opposite
signs give three 14-gonal faces on a genus-three surface. The mechanism is
the explicit two-step face permutation

\[
 \psi^2(p,d)=(p+(a-1)d,ab\,d),\qquad
 a=2^{\varepsilon_L},\ b=2^{\varepsilon_P}\pmod7.
\]

When ab has order three, faces have length six. When ab=1, a nonzero
translation has order seven and faces have length fourteen. This explains
the result without reading genus off a group order.

Following the tournament's actual directed neighbor cycles selects
`(epsilon_P,epsilon_L)=(-1,+1)`, hence genus three. Reversing just the
in-neighborhood convention selects the torus. Both maps have vertex-,
edge-, and face-transitive full map groups of order42, with two flag
orbits. This is a concrete continuation of the noble/flag work; it is
about combinatorial maps, not Euclidean realization of a noble polyhedron.

| Retained data | equal signs | opposite signs |
|---|---:|---:|
| vertices / edges |14 /21|14 /21|
| faces |7 hexagons|3 fourteen-gons|
| genus |1|3|
| dual |K7|K3 with seven parallel edges per pair|
| cycle-space dimension |8|8|
| face-boundary dimension |6|2|
| first-homology dimension |2|6|

The genus-three Heawood map is **not** the Klein heptagonal map: its faces
have length14 and it has14 vertices. Independent `(2,3,5)` and `(2,3,7)`
matrix controls give the dodecahedral and Klein carriers with respective
`(V,E,F)=(20,30,12)` and `(56,84,24)`. All16,384 labelled Heawood rotations
were enumerated; exactly four preserve the fixed Paley subgroup. Proof,
complete universe, and primary Klein-map context are in
[the rotation note](geometry_567_flag_carriers_20261004.md).

## 4. Two famous graph problems enter through an exact flow quotient

Replacing each of the seven Paley triangles by a star gives an exact
sequence for any finite abelian group A:

\[
 0\longrightarrow A^7\longrightarrow
 \operatorname{Flow}_A(K_7)\longrightarrow
 \operatorname{Flow}_A(\mathrm{Heawood})\longrightarrow0.
\]

The seven lost coordinates are constant triangle circulations. A nowhere-zero
star has exactly `|A|-3` nowhere-zero triangle lifts. Thus the fibre count
is `(|A|-3)^7` when all projected edges remain nonzero. A constant nonzero
Paley circulation projects to zero everywhere, so the guard cannot be
omitted.

The surface version exposes the exact limit of planar coloring/flow
duality. There are48 nowhere-zero `F2^2` flows on Heawood. In the torus map,
none is a face boundary, because its dual K7 cannot be colored with four
colors. In the genus-three map, six are boundaries. The total flow space
is unchanged; its boundary subspace depends on the rotation. For `Z/5`
the corresponding counts are1692 total, zero torus boundaries, and12
genus-three boundaries. This puts the four-color connection on explicit
objects, with homology as the missing coordinate.

Tutte's 5-flow conjecture supplies a second precise target. Seymour's
general 6-flow theorem is established; a counterexample minimal within a
fixed-surface-genus class has no face boundary shorter than seven, by
[Möller--Carstens--Brinkmann](https://onlinelibrary.wiley.com/doi/abs/10.1002/jgt.3190120208).
For a cubic orientable cellular map with all face lengths at least seven,
Euler gives `V<=28(g-1)`. The cited short-face reduction and this elementary
conditional inequality are separate steps.

Five here is a flow bound, six is a proved bound and separately a flat
face size, and seven is a minimal-counterexample face threshold. Their
roles are not interchangeable. Indeed Heawood admits an explicit integer
3-flow even in its negative-curvature embedding. The current
[flow-reconfiguration preprint](https://arxiv.org/html/2512.17342v4) refutes
the analogous connectivity claim for5-flows, not Tutte's existence
conjecture. Full proofs, finite counts and literature boundaries are in
[the flow note](geometry_flow_obstructions_20261004.md).

## 5. The Collatz geometry that preserves the actual endpoint

A legal valuation word w with j letters and exponent sum A has map

\[
 F_w(n)=\frac{3^j n+B}{2^A},\quad
 r_w=\frac{3^j}{2^A},\quad a_w=\frac{B}{2^A-3^j}.
\]

As a real upper-half-plane map it is hyperbolic for every nonempty word:
`3^j!=2^A`. Its axis is the vertical line through a_w. The identity

\[
 F_w(n)-n=(r_w-1)(n-a_w)
\]

recovers descent when the integer source and exact valuation guards are
retained. The axis and multiplier together determine the supplied word;
the multiplier alone forgets the ordered carry B. Two word maps have the
same axis exactly when their marked words are positive powers of the same
primitive word. Distinct patterns require ordered composition data.

No fixed additive rational weighting of valuation letters1 and2 can
classify all multiplier drift signs. Its threshold is rational; the true
threshold is `log_2(3)-1`. The inherited -17 word also supplies actual
arbitrarily long growing positive prefixes with a negative formal local
score. The faithful drift coordinate is `j log3-A log2`, or its exact
integer comparison, rather than a rational curvature surrogate. See
[the axis and drift proof](geometry_collatz_drift_carriers_20261004.md).

The projective tournament offers another exact observation. Modulo a prime
`p>3` with `p=3 mod4`, a word multiplies the Paley pair orientation by
`chi_p(3)^j chi_p(2)^A`. Across all these primes, these multipliers retain
exactly the two parities `(j mod2,A mod2)`. Words1 and3 have identical
signatures but opposite drift. This loss concerns the orientation
multipliers; it does not limit every possible finite-field observable.

At p7 the normalized eight-point determinant tournament has21 strict
symmetries, whereas its switching class has168. Keeping a signed lift
gives a faithful336-element SL2 action on16 vertices with eight antipodal
ties. Those ties are essential. Actual Collatz words still fix infinity;
the other projectivities require an observer change that is not a Collatz
operation. This recovers the physical-frame boundary of
[THM-2626, Paley--Borel frame](../../01-canon/theorems/THM-2626-paley-borel-projective-frame-torsor-and-physical-c13-boundary.md).

## 6. A theorem about finite observations, with both routes already home

**PROVED:** for every modulus M and observation depth D, two positive odd
integers can have identical residues for their first D accelerated edges,
the same first D valuations, and a common later endpoint u, while one has
grown to u and the other has fallen to u. Both then reach1 with explicitly
constructed first-hit certificates and equal odd first-hit rank.

Write `M=2^h 3^k m`, with `gcd(m,6)=1`, and put

\[
 j=D+h+1,\quad u=\frac{2^{3^j+1}-1}{3},\quad
 L=2j\varphi(3^{j+k}m).
\]

The sources are

\[
 n=\frac{2^j(u+1)}{3^j}-1,
 \qquad n'=\frac{2^j(2^L u+1)}{3^j}-1.
\]

They satisfy `n<u<n'`, with exact words `1^j` and `1^(j-1)(1+L)` to u.
The first D edges all have valuation1, and all D+1 observed vertices
agree modulo M. The smallest control is `3 ->5 ->1` and `53 ->5 ->1`.
The full statement covers arbitrary M,D; the finite audit has240 parameter
pairs,122 expanded independent replays, and one symbolic pair whose sources
each exceed `10^35` binary digits.

This isolates a genuine limit on fixed finite observations of later drift.
It does not exclude adaptive precision, source height, retained full words,
or a convergence proof using such data. Both sources already converge.
The construction and signed/projective controls are in
[the observation proof](paley_geometry_observations_20261004.md).

## 7. Consequences for the next research step

The recommended Collatz target is now **adaptive certificate selection for
the exceptional least positive cylinder members**. The incoming sharp
bound `B/[Q(Q-P)]<=1121/3328` already handles higher positive lifts of a
contracting coarse cylinder. A selector should retain the original source,
the exact guard, and the endpoint relative to its axis. Test any proposed
finite state summary first against the alias construction above and the
known165 boundary, before trying to prove completeness. Universal selection
remains open; the present session makes its required information clearer.
There are two global obligations: reach a usable contracting word from
every source, and settle the exceptional least members when encountered.
The sharp carry bound alone supplies neither universal selection step.

For the graph lane, the next bounded experiment is to transport **flows
together with the seven circulation coordinates and their homology class**
through the inherited Petersen/Heawood triangle moves. The decisive test is
whether a proposed nonzero-flow repair can be made without creating a zero
edge. The `q-3` fibre count gives an exact local budget, while the constant
Paley circulation is the immediate hostile. A general graph-existence
claim needs more than the two embeddings computed here.

For the recursive storage lane, use typed records: a map stores its dart
permutations; a Collatz macro stores its source guard and ordered affine
data; a flow quotient stores its kernel coordinate. This is a design
principle supported by three separate mechanisms, not an asserted
equivalence between topology and Collatz. LRC(14) receives a boundary
lesson about observer changes, but no new runner certificate or progress
claim from these graph examples.

## 8. Reproduction and independent audit

Each proof note links its standalone script and exact output. Run the four
scripts with normal Python and again with `-O`; all mathematical checks
remain active and the two outputs agree. The geometry and flow notes were
independently cross-audited, as were the Collatz-axis proof and the
all-modulus observation construction. No new scarce theorem identifiers
were reserved, and no open result was promoted on finite evidence.
