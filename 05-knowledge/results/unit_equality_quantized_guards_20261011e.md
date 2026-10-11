# Fixed-universe deletion, colour-orbit quanta, and geometric search guards

**PROVED:** the elementary orbit divisibility and fixed-universe implications,
and the extreme-point refinement of the inherited cherry bound.
**CITED:** the two external edge-bound inputs to the 27-vertex necessary guard.
**FINITE-EXACT:** the declared enumeration and authenticated small controls.
**OPEN:** a six-colourable unit-distance fragment with unit-separated,
nonadjacent ports forced to the same colour. No such fragment is constructed.

## Inheritance and working board

The closest mechanism is the named boundary-count table in
[boundary_context_weights_20261011d.md](boundary_context_weights_20261011d.md).
The hostile is the four-vertex unit rhombus: one blocked Kempe swap in each
alternative colour does not exclude other separating colourings. The corrected
near miss is normalized rarity mistaken for a zero absolute count. The
least-used relevant sidecar is the critical core, together with the old
unit-distance edge bounds:
[THM-421, cherry floor](../../01-canon/theorems/THM-421-unit-distance-3n-floor.md)
and [THM-431, small-point edge table](../../01-canon/theorems/THM-431-u21-equals-57-and-the-3N-floor-tightens-to-25-28.md).
The current board is **fixed universe / named ports / colour orbit / critical
core / geometric boundary / conditional survivor**.

The anchor is the six-colour equality fragment; the niche is an exact integer
threshold for deleting its separating assignments; the wildcard is the
extreme-point obstruction behind an apparently feasible small graph.

## 1. The exact target and its chromatic premise

Let F be a finite simple graph with distinct nonadjacent ordered ports u,v.
Fix k>=2. Let e count extensions of port colours (0,0), and d extensions of
(0,1). Thus

\[
 P_F(k)=ke+k(k-1)d.
\]

The target is e>0 and d=0. If the vertices have distinct plane positions,
every edge has length one and the omitted port edge also has length one,
then F+uv is a unit-distance graph with chromatic number k+1. Conversely,
deleting an edge from a subgraph-minimal (k+1)-chromatic unit graph gives
exactly such an interface. This is an equivalence of finite witnesses, not
an easier conjecture that can be assumed as a construction step.

Every equality tester needs all k colours: if F admitted a colouring using
at most k-1, either the ports already differ or one can recolour u with an
unused colour. The latter keeps all edges proper and separates the ports.
The fragment is a unit-distance **subgraph**, not the induced graph of its
point set, since the omitted pair is itself at unit distance.

## 2. Orbit quantization supplies a larger zero threshold

Let a^E_j and a^D_j count partitions of V(F) into j nonempty independent
sets with the ports in the same or different blocks. They count unlabelled
colour classes, not just arbitrary set partitions. Write (a)_b for a falling
factorial. Assigning the remaining colour names gives

\[
 e=\sum_{j=1}^{k}(k-1)_{j-1}a^E_j,\qquad
 d=\sum_{j=2}^{k}(k-2)_{j-2}a^D_j.                 \tag{1}
\]

Suppose independently that chi(F)>=s, where 2<=s<=k. No partition with
j<s occurs. Every remaining coefficient in each sum is an integer multiple
of its first possible coefficient. Therefore

\[
 (k-1)_{s-1}\mid e,\qquad (k-2)_{s-2}\mid d.       \tag{2}
\]

At k=6 the equal coefficients are 1,5,20,60,120,120, and the different
coefficients are 1,4,12,24,24. In particular

\[
 \boxed{\chi(F)\ge5,\quad e>0,\quad d<24
              \quad\Longrightarrow\quad d=0.}    \tag{3}
\]

The positive witness, the lower chromatic bound, and the absolute upper
count must all concern the same actual fragment. Inequality (3) does not
provide those premises. The version with chi(F)>=4 instead has threshold12;
with chi(F)>=3 it has threshold4. Without a chromatic premise the universal
integer threshold remains1. Equality at the threshold does not suffice.

This is where a specific existence certificate and a general counting
certificate genuinely combine. The small verifier authenticates both
numerical premises by complete enumeration; it never accepts a caller's
claimed bound on trust. For large graphs, one needs independently checked
colouring-unsatisfiability and count-bound certificates in their place.

## 3. The right set on which exponential deletion can finish

Keep a finite vertex set fixed. Let F_0 subset F_1 subset ... subset F_r
be graphs on that same set, with the same omitted ports, and retain one
proper equal-port k-colouring throughout. Write D_i for the actual set of
fixed-(0,1) proper colourings of F_i. Then D_i subset D_(i-1), and, while
the denominator is positive,

\[
 d_r=d_0\prod_{i=1}^r\rho_i,
 \qquad \rho_i=|D_i|/|D_{i-1}|.                    \tag{4}
\]

Suppose an authenticated adaptive rule keeps rho_i<=sigma<1, an independent
bound gives d_0<=M, and chi(F_r)>=s. Then

\[
 M\sigma^r<(k-2)_{s-2}\quad\Longrightarrow\quad d_r=0. \tag{5}
\]

In particular sigma=5/6 or6/7 is meaningful **if it bounds the actual
conditional factors**. It is not automatic from the palette or from the
fraction deleted under an unconditioned assignment. The finite ambient
universe is crucial: adding new vertices introduces new assignments and
changes the absolute-count budget. The accompanying
[amplifier packet](existence_count_amplifiers_20261011e.md) proves that
private-interior serial/parallel repetition cannot attain d=0 from d>0 at
k>=3, despite arbitrarily small normalized ratios.

An explicit finite positive control starts with seven isolated vertices at
the exact inherited Moser positions and the ports of one omitted spindle
edge. Add the remaining ten edges, always preserving an equal-port witness,
and greedily choose the edge deleting the most surviving separating
assignments. At k=3, their counts are

\[
243,162,81,54,36,18,12,6,4,2,0.
\]

At k=6 the **same source and edge-addition rule** gives

\[
7776,6480,5184,4320,3600,2880,2400,1920,1600,1328,1056.
\]

The exact latter factors are
5/6,4/5,5/6,5/6,4/5,5/6,4/5,5/6,83/100,66/83.
This illustrates both the real role of5/6 and the surviving obstruction:
every available edge has been used, but separating assignments remain.
The rule has a verified positive small-palette control and an explicit
failure at the target palette. It is not a supplier of (5) for arbitrary
larger configurations.

## 4. A boundary refinement of the inherited counting floor

Any two distinct points have at most two common unit neighbours. Hence for
a unit-distance graph on n points,

\[
 \sum_v\binom{\deg v}{2}\le2\binom n2.
\]

If its minimum degree is delta>=2, this implies
n>=1+delta(delta-1)/2. Equality would force every degree to equal delta
and every vertex pair to have exactly two common neighbours. Consequently
every vertex neighbourhood would induce a 2-regular graph.

Choose a vertex strictly extreme for a generic linear functional. Its unit
neighbours lie in an open semicircle around it. Two such points are adjacent
only when their angular coordinates differ by pi/3. A finite set of angles
in an interval of length pi with these adjacencies has no cycle: take the
smallest angle in a purported cycle. Its only possible larger neighbour is
the angle plus pi/3, so it cannot have two cycle neighbours. Thus this
neighbourhood cannot be 2-regular. Equality is impossible, proving

\[
                    n\ge2+\delta(\delta-1)/2.       \tag{6}
\]

The Shrikhande graph has16 vertices, degree6, every codegree2, and cyclic
neighbourhoods. It saturates the counting relaxation and fails the geometric
boundary condition. THM-421 already uses this hostile; (6) explains the
equality failure directly. This elementary guard is not a new best-known
lower bound on high-chromatic unit graphs.

## 5. Recovered edge bounds give a stronger 27-vertex guard

Inside F+uv choose a subgraph-minimal7-chromatic graph H. It must contain
uv since F is six-colourable. Its minimum degree is at least6, so (6)
already gives n=|V(H)|>=17.

**CITED input 1.** Kostochka--Yancey's critical-edge theorem, specialized to
k=7, gives |E(H)|>=ceil((10n-7)/3). See
[Ore's conjecture on color-critical graphs is almost true](https://arxiv.org/abs/1209.1050).
**CITED input 2.** The unit-distance upper bounds in Alexeev--Mixon--Parshall,
[Theorem1](https://arxiv.org/html/2412.11914v2), already routed through
THM-431 above, exclude this edge count for every17<=n<=26. The last rows are:

| n | Critical graph requires at least | Unit graph has at most |
|---:|---:|---:|
|25|81|78|
|26|85|84|
|27|88|90|

Thus every proposed six-colour unit-separated equality fragment has at
least27 vertices. These particular bounds stop at27; they do not assert
that a27-point candidate exists or that27 is the best available bound.
The input publications are cited, not re-proved by the accompanying code.

## 6. What is retained and what is still missing

| Map | Preserved predicate | Lost coordinate / required sidecar | Decisive test |
|---|---|---|---|
| Colourings to independent-set partitions | Exact equality/separation count after orbit factors | Actual colour names; retain a witness assignment | Independent fixed-name enumeration |
| Assignments to a survival fraction | Relative surviving count | Absolute ambient count and lower orbit quantum | Check (5), not merely sigma<1 |
| Unit fragment to critical core | A smaller seven-chromatic unit subgraph containing the added edge | Noncritical vertices, unnecessary edges | Retain inherited injective coordinates |
| Geometry to codegree bounds | Necessary counting inequalities | Position of an extreme point | Shrikhande hostile |

The next construction obligation is a fixed finite point universe outside
the four-colourable Moser-field bank, with a retained equal-port colouring,
and enough actual unit constraints to kill every separated assignment.
The next proof obligation is an authenticated source-specific count bound
below the relevant quantum. Neither an unconstrained graph selector nor
an asymptotic probability estimate supplies the missing geometric object.

## Reproduction and audit

Run `python 04-computation/experiments/unit_equality_quantized_guards_20261011e.py`.
The universe is all1098 labelled simple graphs on2..5 vertices and all5325
nonedge port markings, palettes2..7. Partition-based formulas are compared
with independent fixed-name backtracking for every graph on2..4 vertices.
The familiar count1099 is recovered by adjoining the one-vertex graph;
there is no extra palette-threshold claim in that total.
Positive controls are the abstract diamond/K7-minus-edge receipts and the
actual Moser deletion at palette3. Hostiles include false chromatic premises,
forged count bounds, invalid witnesses, Boolean/integer aliases, the same
Moser deletion at palette6, and Shrikhande as a non-geometric counting
extremum. The saved `.out` has48,734 checks; ordinary and optimized Python
must agree. Algebraic and geometric proof steps above are separate from the
bounded enumeration. Independent agent audit covered orbit divisibility,
the extreme-point proof, and the cited-bound consequence.
