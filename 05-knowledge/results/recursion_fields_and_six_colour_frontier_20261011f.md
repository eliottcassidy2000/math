# Recursion, divisor records and the six-colour equality frontier

2026-10-11. **PROVED scoped constructions and obstructions; FINITE-EXACT
candidate graphs and certificate checks; CITED sequence identification.**
No six-colour equality fragment has been found. All three new point sets
below are exactly four-chromatic. The plane six/seven target and universal
Collatz coverage remain **OPEN** in this packet.

## 1. What a successful configuration must contain

The inherited target is a finite unit-distance graph F, distinct ports u,v
at unit distance with their edge omitted, one supplied proper six-colouring,
and a refutation of every proper six-colouring with different port colours.
Adding uv then gives a seven-chromatic unit graph. Thus this is a full
existence problem, not a final edge-edit after a purely arithmetic escape.

Closest proved mechanism: complete marked boundary support and critical-edge
extraction in [existence/counting alternation](existence_count_equality_and_runner_20261011e.md).
Canonical hostile: a full point set with a proper four-colouring separates
every omitted unit pair, even if its field is large. Corrected near miss:
a norm-one rotation does not guarantee a realizable endpoint displacement.
Least-used relevant sidecar: the norm valuation at the obstructing prime.

The session board was **unit incidence / marked boundary / defect support /
prime decomposition / attainable displacement / actual colour witness**.
The anchor was six-colour equality; the niche was historical recursion;
the wildcard was the user's reciprocal equation. Incoming commit d6246592
supplied the [golden and taxicab boundary gadgets](golden_taxicab_boundary_20261011.md):
useful exact list-colouring interfaces, whose prescribed leaf colours still
need enforcement when used inside an unprecoloured construction.

## 2. Field escape can be made systematic, but does not force more colours

[The defect/field packet](unit_defect_frontier_20261011f.md) proves that for
any finite prime bank S, setting

\[
 P=\prod_{p\in S,\ p\ne2}p,\qquad
 t=2P,\qquad D=8P-1,\qquad
 \omega=\frac{2t-1+\sqrt{-D}}{2t}
\]

gives a unit rotation whose imaginary quadratic field splits at every prime
in S. Any conjugation-stable extension containing that field has no
conjugation-fixed prime over S. This removes all corresponding instances of
the inherited residue screen, including those with higher residue degree.
It does not remove other primes or other colouring methods. The inherited
screen is [THM-4558, residue colourings of field planes](../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md),
with the MISTAKE-576/579 corrections retained.

A second guard is essential. In the original Moser field
M=Q(sqrt(-3),sqrt(-11)), a rational norm z conjugate(z) has even 2-adic
valuation. The first four even parameters passing our displayed split-subfield
test at primes through19, namely22,32,38,46, fail that displacement guard.
The next passing parameter60 has an actual endpoint:

\[
 z=\sqrt{33}+3\sqrt{-3},\quad |z|^2=60,\qquad
 \omega_{60}=\frac{119+\sqrt{-239}}{120},\quad
 |z-\omega_{60}z|=1.
\]

Using nine bent diamond modules reaches z. Joining a rotated copy and then
thickening by a unit triangle produces these exact full unit graphs:

| Point set | Vertices | Unit edges | Certified chromatic number |
|---|---:|---:|---:|
| Bent diamond chain |28|46|4|
| Chain plus rotated copy |55|93|4|
| Triangular thickening |153|405|4|

All coincidences and all unit pairs are computed in an exact eight-dimensional
rational basis. Literal four-colourings are exported in
[the complete coordinate output](unit_prime_escape_candidate_20261011f.out).
The chain already has an extra unit edge obstructing three colours; it must
not be described as a three-colourable equality seed. These are actual
arithmetic escapes, and failed six-equality candidates with known witnesses.

The [all-length diamond-chain construction](unit_field_escape_recursion_20261011f.md)
goes further: for every m>=1 it gives a four-critical unit graph with6m+1
vertices and10m+1 edges. Taking m=2P escapes every prime in S while retaining
chromatic number four. Deleting any edge gives a genuine three-colour equality
fragment. This shows, at every finite scale, why more arithmetic freedom alone
does not supply the missing six-colour relation.

## 3. Recursion must preserve the relevant interface

The [historical atlas](recursion_escape_atlas_20261011f.md) recovers ten
operation families, with theorem IDs, exact source files, missing data and
hostile controls. Its useful distinctions include:

| Historical operation | Preserved quantity | Consequence for this target |
|---|---|---|
| Tournament endpoint padding / ordinal join | Hamiltonian count stays H / becomes H squared | Neither statement preserves a six-colour boundary relation |
| Tournament substitution | Vertex count becomes n squared | It differs from ordinal join even on C3; operation order matters |
| SCC/Lyndon type codec | Full addressed type can be recovered | Reversible coding can change chromatic number and supplies no metric placement |
| Barycentric / edgewise refinement | Incidence flags / multiplicative lattice scale | Both have inherited low-colour structures; equal face counts lose adjacency |
| Guarded affine Collatz words | Source cylinder, common future and valuation fuel | Repetition still needs entry and terminal certificates |
| Elliptic duplication | Curve equation and rationality | Positivity and geometric forcing require separate tests |

An additional geometric obstruction rules out a tempting shortcut:
Mycielski lifting any graph with a degree-three vertex creates two vertices
with three common neighbours. That K2,3 cannot be injectively realized by
plane unit distances: two distinct unit circles have at most two intersections.
An abstract chromatic-raising recursion therefore needs a different local
incidence construction before it can help here.

There is also a useful exact relaxation. A six-equality seed at **any positive
port distance d** suffices. Chain m translated copies until md>=1/2, rotate
the chain so its two far tips are unit separated, and add that edge. For an
n-vertex seed this uses at most2m(n-1)+1 points. Copy maps retain every actual
coincidence; no intermediate colourability is assumed. Edge-critical extraction
then recovers a unit-port equality fragment. This broadens geometric search
parameters but is existentially equivalent to the original problem.

## 4. The two number sequences record an object and its capacity

The first list contains the strict record arguments of the divisor function
tau; the second contains their record values. In particular

\[
 \tau(24)=8,\quad\tau(36)=9,\quad\tau(240)=20,
 \quad\tau(360)=24.
\]

The final24 in the second supplied list corresponds to the next first-list
entry360. These are [OEIS A002182](https://oeis.org/A002182) and
[A002183](https://oeis.org/A002183). The exact recursion in
[the divisor/equation packet](divisor_equation_recursions_20261011f.md)
enumerates nonincreasing prime exponents and retains strict divisor records.

The exponent box is the richer object behind tau=product(a_i+1). Its cover
graph is bipartite; its full comparability graph contains a plane-unit-impossible
K4 once the exponent sum reaches three. Thus retaining the box is necessary
to choose an incidence rule, and neither of these immediate choices provides
the desired colouring constraint. Counts are useful capacities, not a
replacement for the geometry or marked boundary relation.

## 5. The reciprocal equation supplies an exact guarded recursion

For A/(B+C)+B/(A+C)=C/(A+B), put s=A+B. The maps

\[
 x=\frac{4(s+2C)}{2s+C},\qquad
 y=\frac{12(A-B)}{2s+C},\qquad
 A:B:C=(8-x+y):(8-x-y):(4x-8)
\]

identify its valid rational solutions with a guarded chart on
E: y squared=x(x squared+5x-32). The original pair-sum denominators must
remain nonzero. Positive representatives are exactly the real points with
2<x<4. Starting at P=(4,4), elliptic addition reaches a positive primitive
integer triple at8P:

\[
 (A,B,C)=(301361533449900458837600,
 49105016933436320224063,
 316629033253501281102807).
\]

All original fractions are checked exactly. The resulting nondegenerate
triangle supplies a genuine unit rotation outside M, rather than an unrelated
circle parametrization. However, its joined field retains the inherited
prime11 five-colour screen. This is a successful equation construction and
an eliminated six-equality carrier, not a reason to claim general coverage.
No leastness, full elliptic rank or infinite positive family is asserted.

## 6. A recursive state that keeps the decisive information

For a graph G and any k-colour assignment c, retain the actual monochromatic
edge set M_G(c), then keep the inclusion-minimal such sets A_k(G). For every
deleted set D,

\[
 G-D\text{ is k-colourable}\iff
       \exists A\in A_k(G)\text{ with }A\subseteq D.
\]

If a minimal A is nonempty, every proper k-colouring of G-A identifies both
ends of **each** edge in A. If G is already k-colourable, the minimal family
is only the empty set, ruling out forced equality on edges removed from E(G).
When G is the full induced unit graph, this rules out all omitted-unit-edge
equality fragments on that point set. A two-edge unit path with unit-separated
ports shows why that full-incidence qualification matters. This is an exact
obstruction representation, not a new
existence proof or an efficient general colouring algorithm.

For recursive gluing, retain a separate antichain for each equality partition
of the labelled boundary. Compose by unions of compatible defects and remove
dominated sets. The two-edge path and one edge at palette two demonstrate why
forgetting the port pattern before gluing loses the triangle obstruction.
The formal compiler authenticates the source graph, exact labels and palette.

The next concrete requirement is a geometric overlap that changes this full
boundary support: preserve one equal-port witness while hitting **every**
separating assignment. Private-leaf list gadgets supply candidate constraints;
their own colour patterns still need to be forced or jointly eliminated.
Ordinary independent series/parallel duplication cannot supply the missing
zero by multiplying positive support. The new arbitrary-distance compiler
allows the ports to be chosen for rigidity first and their final unit distance
handled afterward.

## 7. Evidence and scope

Five experiment/output pairs accompany the four detailed proof packets:

| Packet | Exact checks |
|---|---:|
| Divisor records and reciprocal equation |1,897|
| Historical atlas and any-distance compiler |237|
| All-length field-escape recursion |2,239|
| Defect profiles and prime/norm guards |2,901|
| Three exact geometric candidates |109|
| Total |7,383|

Normal and optimized Python outputs agree for each packet. Independent peer
audits cover the mathematics and controls; the profile audit found and repaired
Python equality aliases, with nine retained hostile tests. The detailed notes
give reproduction commands, finite universes and saved-output hashes. Counts
of checks are validation evidence, not a measure of progress toward a global
colouring or Collatz proof.
