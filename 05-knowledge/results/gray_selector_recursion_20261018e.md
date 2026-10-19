# Four charts, ten coordinates, and the recursion carried by the missing bits

**Status: PROVED elementary chart, Boolean-ring and marked-codec statements;
FINITE-EXACT reproduction and exhaustive small-code controls; CITED classical
quaternary-code construction. Geometric selector realization and the major
conjectures remain OPEN.** Self-contained proofs and independent agent audits;
no historical-priority claim or new theorem ID.

2026-10-18, continuation of the marked boundary-fibre work. This packet answers
the owner's request to recover the old “path of least resistance” mechanisms
behind the ten-coordinate relation's four charts.

The [completed equality-case classification](hamming_matching_atlases_20261018e.md)
continues this result: 28 maximum charts yield exactly 105 optimal
four-chart partitions, with an explicit dual-matching dictionary and a
retained central phase.

## 1. Inheritance and the objects that must stay separate

The anchor is the exact relation in
[ten-label Hamming replacements](ten_label_matching_hamming_20261018n.md),
consumed as component-flip charts by
[the pinned-fibre atlas](boundary_fibre_atlas_20261018d.md).
The niche is the earlier tournament mirror/blue recursion; the wildcard is
the arithmetic carry in a quaternary representation.

- Closest proved mechanisms: bipartite components have one free binary phase;
  the [marked four-tournament section](partial_four_tournament_20261011c.md)
  has two independently reversible named arcs.
- Canonical hostile: gluing mutually alternative charts on the same ports
  imposes intersection, not union. The old odd-cycle and pinned-path guards
  remain necessary. A tournament path count also merges the two diamond
  states, so it is not a decoder for the four choices.
- Corrected near miss: [THM-452, doubling tiling coset law](../../01-canon/theorems/THM-452-doubling-tiling-coset-law.md)
  and MISTAKE-065 retain the actual affine defect; MISTAKE-642 retains the
  Hadamard normalization. The new carry hostile below passes every unmarked
  count and still corrupts the source.
- Least-used relevant sidecar: [THM-830, deletion deck and mirror-current
  calculus](../../01-canon/theorems/THM-830-b3-deletion-deck-mirror-current-calculus.md),
  Section 6, keeps both mirror halves, the boundary trace and their defect.
  This gives an actual older carrier for our present relation.

Board: disjoint component phases; marked four-state selector; quaternary
carry; mirror defect; Boolean multiplication; arithmetic modulus versus
coordinate count. The old Fibonacci resistance ladder supplies a comparison
of operations, not a theorem identifying these six objects.

All codes below are binary linear subspaces on **named coordinates**.
Decimal numbers such as 170 are bit masks with coordinate 0 least significant.
Let I2={00,11}. Our relation is

~~~
C = RM(1,3) direct-sum I2
  = span(15,240,768,170,204),       dim C=5, |C|=32.
D = span(15,240,768),              dim D=3, |D|=8.
~~~

The code RM(1,3), also the extended length-eight Hamming code H8, consists of
all affine binary functions on the eight addresses in F2^3. The distinguished
pair is coordinates 8,9. One chart is an affine set with independent flips
on disjoint nonempty blocks; coordinates outside those blocks are fixed.

## 2. The exact general cost of component charts

For a binary linear code C of dimension k, define

\[
\nu(C)=\max\{r:\ d_1,\ldots,d_r\in C\setminus\{0\},
                   \ \operatorname{supp}(d_i)\cap
                   \operatorname{supp}(d_j)=\varnothing\ (i\ne j)\}.
\]

**Proposition 1 (PROVED).** The minimum number of component-flip affine
charts whose union is exactly C is

\[
\boxed{\operatorname{charts}(C)=2^{\kappa(C)},\qquad
       \kappa(C)=\dim C-\nu(C).}
\]

“Exactly” excludes illegal extra words. Empty charts can be discarded. If
r+span(d1,...,ds) is a chart contained in C, then r and all di belong to C.
Disjoint nonempty supports make the di independent, so s<=nu(C) and each
chart has at most 2^nu(C) words. Covering 2^k words needs at least
2^(k-nu(C)) charts, even with overlapping charts.

Conversely choose a maximum disjoint packing and let D be its span. Its
2^(k-nu(C)) cosets partition C and are eligible charts. This proves equality,
including C={0}, whose one chart fixes every coordinate.

Every such chart has an abstract pinned bipartite-graph realization. In each
block choose a root. Join another coordinate to it by an odd path when their
prescribed offset bits differ, and an even path when they agree; length one
or two suffices. Pin every uncovered coordinate. Hidden length-two vertices
have a unique shade once the boundary is chosen. This is a binary inequality
graph, not a claimed unit-distance embedding or an unrestricted six-colour
realization. Fixed ternary pair labels are still part of the fibre method.

If 1 belongs to C, a maximum packing covers every coordinate: otherwise
1+sum(di) is a further nonzero disjoint codeword. Thus the optimal charts of
our Hamming/repetition family need no boundary pins. Arbitrary codes do.

**Corollary 1 (PROVED, algebraic interpretation).** One chart suffices
exactly when C is closed under coordinatewise multiplication, i.e. bitwise
AND. More generally, nu(C) is the largest dimension of a Boolean subring
contained in C. Subrings need not contain the ambient all-one vector.

Disjoint-support generators span a Boolean ring. Conversely the minimal
nonempty supports of a finite Boolean subring are pairwise disjoint:
intersecting distinct atoms would have a smaller nonzero product. They span
the ring, by successively removing atoms from any word. Its identity is the
union of its active coordinates. This proves both statements, including zero
and proper-support codes.

This is a precise interpretation of “least resistance” **for this interface**:
retain as many independent component bits as possible; the remaining k-nu
bits select among charts. It is not a claim of universal optimal compression
under arbitrary encodings. C itself is one unrestricted linear subspace;
the four-chart obstruction concerns the more restrictive component geometry.

## 3. Which recursions change that cost?

**Proposition 2 (PROVED).**

\[
\nu(C\oplus E)=\nu(C)+\nu(E),\qquad
\kappa(C\oplus E)=\kappa(C)+\kappa(E).
\]

A disjoint packing projects to disjoint nonzero words in each summand after
zero projections are removed. Every packed word uses at least one summand,
so its number is at most the sum of their two packing maxima. Pure packings
from the separate summands attain equality. Diagonal duplication
C -> {(c,c):c in C} preserves dimension and disjointness, hence preserves
kappa. Appending I2 also preserves kappa: it adds one dimension and one
independent component.

For m>=1, every proper nonzero affine function on F2^m has weight 2^(m-1);
the other nonzero function is constant one. Hence nu(RM(1,m))=2, achieved by
the two complementary half-cubes. Consequently, for every q>=0,

\[
\boxed{\operatorname{charts}(RM(1,m)\oplus I_2^q)=2^{m-1}.}
\]

| Relation | Length | Dimension | Maximum independent disjoint flips | Charts |
|---|---:|---:|---:|---:|
| RM(1,2) | 4 | 3 | 2 | 2 |
| H8=RM(1,3) | 8 | 4 | 2 | 4 |
| H8 plus one I2 | 10 | 5 | 3 | 4 |
| H8 plus two I2 | 12 | 6 | 4 | 4 |
| RM(1,4) | 16 | 5 | 2 | 8 |
| Five independent I2 pairs | 10 | 5 | 5 | 1 |

The Plotkin step

\[
RM(1,m+1)=\{(f,f+\epsilon\mathbf1):f\in RM(1,m),
                                             \epsilon\in\mathbf F_2\}
\]

therefore doubles the chart count. Literal duplication does not. Both double
length; only the former introduces the additional interaction/selector bit.
This separates operations that earlier informal recursion pictures could
otherwise merge.

There is also an exact graphic specialization. For the binary **Eulerian
edge-subset** code Z1(G;F2), nu is the maximum number of edge-disjoint simple
cycles. Every nonzero even subgraph contains a cycle; shrinking each member
of a disjoint packing to one cycle proves the upper bound, and cycles give
the reverse bound. Thus

\[
\operatorname{charts}(Z_1(G;\mathbf F_2))
 =2^{\,|E|-|V|+c(G)-\text{maximum edge-disjoint cycle count}}.
\]

K4 has dimension 3 and packing number 1: four charts on six coordinates.
K5 has dimension 6 and packing number 3: eight charts on ten coordinates.
The latter packing is triangles 012 and 034 and the cycle 1-3-2-4-1.
These are not the “automorphism-even graphs” in the tournament
equinumerosity theorem. Even that vocabulary distinction carries information.

## 4. Literal arithmetic modulo four, with an explicit inverse

The modulo-four connection is stronger than the number of charts. Consider

\[
Q=\{(a,a+b,a+2b,a+3b,2t):
                a,b\in\mathbf Z/4,\ t\in\mathbf F_2\}.
\]

It has 32 words and additive group Z4^2 x Z2. Apply the Gray map

\[
0\mapsto00,\quad1\mapsto01,\quad2\mapsto11,\quad3\mapsto10.
\]

This construction is classical: the first four coordinates are generated
by Krotov's displayed matrix A^(1,0), with rows 1111 and 0123.
The Gray map preserves Lee distance as Hamming distance.
See [Krotov, Z4-linear Hadamard and extended perfect codes,
Sections 2–3](https://arxiv.org/pdf/0710.0199).
Our contribution here is the explicit marked transport to the inherited
chart/partial-tournament interface, not a new quaternary code.

Write a=a0+2a1 and b=b0+2b1. List the two Gray bits of each coordinate
consecutively, and move source positions 0,...,9 to

~~~
(0,1,4,7,2,3,6,5,8,9).
~~~

Then its binary image is exactly

\[
c=170a_0\oplus204b_0\oplus15u\oplus240v\oplus768t,\qquad
u=a_1,\quad v=a_1\oplus b_1\oplus a_0b_0.                 \tag{4.1}
\]

Here the masks are multiplied by bits; sums are XOR. Conversely,

\[
\begin{aligned}
a_0&=c_0\oplus c_1,& b_0&=c_0\oplus c_2,\\
u&=c_0,&v&=c_4,&t&=c_8,\\
a&=a_0+2u,&b&=b_0+2(v\oplus u\oplus a_0b_0).
\end{aligned}
\]

For a direct proof, write j=j0+2j1 and let eta choose the first/second
Gray bit. The low bit of a+jb is a0+b0*j0; its high bit is
a1+b0*j1+b1*j0+a0*b0*j0, all in F2. The Gray bit is high+eta*low.
The retained address change

\[
(\eta,j_0,j_1)\mapsto
(X_0,X_1,X_2)=(\eta,j_1+\eta j_0,j_0)
\]

turns it into a1+a0*X0+b0*X1+(b1+a0*b0)*X2. This proves (4.1)
at each named coordinate, not just equality of the two word sets.

Fixing the two parity bits (a0,b0) gives precisely the four cosets of D.
These selectors form F2^2, not the cyclic group Z4. Actual mod-four addition
is transported with its carry:

\[
\gamma(x+y)=\gamma(x)\oplus\gamma(y)
 \oplus ((x\bmod2)(y\bmod2))(1,1).
\]

**Sharp hostile.** At a=b=1,t=0, the correct binary word is 150; omitting
a0*b0 gives 102. Nevertheless omitting that term preserves the entire code
and each of its four parity-selected charts, since the free b1 absorbs it.
All counts and unmarked coverage tests pass the wrong map. Only retaining
the original source and its inverse detects this failure.

This is the exact methodological connection to the earlier
[golden digit/carry work](golden_digit_carry_20261003.md): a normalized value
or set of outputs can stay correct while the route from a particular source
is lost. The arithmetic operations in the two constructions differ.

## 5. The old blue tiling recursion literally contains this relation

Order-six fixed-path tournament tilings have ten coordinates. In
[THM-477, the blue code](../../01-canon/theorems/THM-477-the-blue-code.md),
the mirror involution has four two-element orbits and two fixed tiles.
Choose and retain the coordinate identification
sigma=(01)(23)(45)(67), with 8,9 fixed. Its fixed space is
B=I2^4 direct-sum F2^2.

Every word of our code is uniquely

\[
(u,u+a,u+b,u+a+b,\ v,v+a,v+b,v+a+b,\ w,w).                \tag{5.1}
\]

Thus (1+sigma)c=a*255. Set R=C intersect B. Then

\[
C=R\sqcup(170+R),\qquad
R=D\sqcup(204+D).
\]

This is a genuine two-stage lossless fold: first retain defect a, then
selector b. On R, collapsing the four mirror pairs gives
(u,u+b,v,v+b), exactly the even-parity length-four code RM(1,2).
Keep w for the remaining repeated pair. R has dimension 4 and nu=3, so the
inner two-chart cost is also minimal. The two selectors yield the exact four.

More intrinsically, THM-830 Section 6 represents a full old tiling as
(h,h+delta,z), where h is the four off-boundary coordinates of a chosen
mirror half, delta its four-bit defect and z the two fixed bits. The full
half includes those two trace coordinates as well. Impose three conditions:

~~~
h has even parity;       delta = a*1111;       z=(w,w).
~~~

Writing h=(u,u+b,v,v+b) and interleaving mirror halves gives (5.1).
This is a restriction of the old full witness carrier, with an inverse.
It is not a claim that C lies entirely in the blue fixed locus.

Indeed B-perp=I2^4 lies inside every self-dual subcode C' of B, since
C' subset B implies B-perp subset C'-perp=C'. On the two remaining fixed
coordinates self-orthogonality forces the line {00,11}. Hence the only
length-ten binary self-dual code contained in B is I2^5. The Hamming code
must escape four simultaneous mirror equalities; the defect bit records how.

The related old [THM-291, Mode-B recursion](../../01-canon/theorems/THM-291-mode-b-multilinear-recursion.md)
adds a source and sink, increases tournament order by two and preserves H.
It does not send H to H+2. The corrected boundary increment is 2n-1 free
tiles when order n becomes n+2. The old n−2/n−1 arm error must not return.

## 6. The four-vertex tournament is the actual marked controller

Retain arcs 0->1->2->3 and 0->3. Orient 02 backwards iff a=1 and
13 backwards iff b=1. The unique marked frame has partial outdegrees
(2,1,1,0). Its four completions, proved in
[the partial-four packet](partial_four_tournament_20261011c.md), are:

| (a,b) | Tournament state | Hamiltonian paths |
|---|---|---:|
| (0,0) | T, transitive | 1 |
| (1,0) | +, triangle above sink | 3 |
| (0,1) | −, source above triangle | 3 |
| (1,1) | S, strong | 5 |

Use (a,b)=(c0+c1,c0+c2) from (5.1), and retain the three phase bits
(u,v,w)=(c0,c4,c8). This is a bijection between the 32 codewords and a
marked tournament state together with three bits. The two named arc flips
are exactly XOR170 and XOR204, preserving those phase bits. Converse plus
path reversal swaps a,b. Neither arbitrary tournament relabellings nor a
third arc flip descends to this four-state interface without more data.

This is the strongest connection to the early tournament microcosm:
the small object is a reusable **controller with specified ports**. Replacing
it by its scalar H merges + and −. Replacing it by an unmarked graph gadget
does not supply a physical choice between our four boundary relations.

There is a separate authentic mod-four spectral law in
[THM-584, antipodal spectrum](../../01-canon/theorems/THM-584-complement-is-antipodal-map-level-parity-spectrum.md).
On the full tournament arc cube with d coordinates, eigenvalues are d-2k;
full converse acts by (-1)^k on Fourier degree k. At n=5, d=10, the converse-even
spectrum is 2 mod4 and the converse-odd spectrum is 0 mod4. This connects
ten and four through a specific operator and grading.

The inherited proof's relabelling action was repaired in this session:
it is x -> P_g x XOR epsilon_g, a signed affine cube action. Swapping
vertices 0 and 1 sends the all-forward word 0 to word 1, so Hamming weight
is not preserved. Fourier degree is preserved by the signed character
action. The theorem's spectral conclusion survives. The reversal mask
is the same kind of indispensable gauge information as the carries above,
not the same numerical invariant.

## 7. What ten and 105 do, and do not, identify

The current ten is 8+2: four quaternary evaluation positions plus a
marked even fifth position. The unique weight-two word is 768, supported
on {8,9}. Every coordinate automorphism preserves this pair. Therefore no
10-cycle on coordinates preserves C; in particular cyclic rotation sends
768 to 513, which is absent. Cyclic translation on the ten coordinate
positions cannot preserve this code. This does not exclude arbitrary
abstract actions of Z/10 on its word set.

Other old tens have genuine but different carriers:

| Old construction | Exact connection | Information that must remain |
|---|---|---|
| Mod-10 multiplication wedge | Fifteen addresses split into five loops, four fixed path arcs and six free arcs | Fold signs, addresses and orientation choices |
| Dotted diagonal on four labels | Ten positions = four repeated pairs plus six distinct pairs | Endpoint labels; its natural S4 orbits have sizes 4 and 6, with no invariant pair |
| Five-vertex tournaments | Ten arc coordinates, plus the mod-four spectral grading above | Affine relabelling mask; C is not an S5-invariant family |
| Six-vertex fixed-path tilings | Ten tiles = four mirror pairs and two fixed tiles | Mirror defect, common boundary, and chosen coordinate identification |
| Present Hamming relation | Eight affine-function samples plus one repeated pair | Distinguished pair, two selector bits, three component phases |

For the first row the actual loop/path/free-arc bijection is in
[modular atoms](tiling_modular_atoms_20261003.md), Sections 1–2. Two free
arcs both receive residue 6 mod10, so residues alone do not recover their
addresses. The dotted-diagonal comparison is in
[medial pair tiles](medial_pair_tiles_20261004.md). Neither cardinality match
automatically transports graph incidence or the desired colouring predicate.
Even an edge-coordinate permutation cannot turn C into the K5 Eulerian
code: the latter has dimension 6 and no weight-two word.

There is a useful exact old/new 105 bridge. Fix the distinguished pair.
The remaining eight labels have 7!!=105 perfect matchings. Each matching
determines an even-pair subcode and exactly two H8 extensions. Thus it
becomes an edge of the 30-vertex, degree-seven Hamming replacement graph,
with 105 edges. The inherited
[ten-label packet](ten_label_matching_hamming_20261018n.md), Section 3,
constructs its point-plane incidence identification in PG(3,2).
This is a matching-to-edge map, not merely two lists of size 105.
The older golden/binary clocks modulo 105 are separate; no map transporting
their time dynamics to these 105 incidence edges is established here.

## 8. A sharper compression target, and the open geometric obligation

Our code is not AND-closed: 170 AND 204 = 136 is not a codeword. Its smallest
Boolean-ring completion has nine independent atoms: eight singleton
coordinates and the distinguished pair. It is F2^8 direct-sum I2, with 512
words. Removing the selector constraints to obtain one component chart
therefore admits **480 illegal extra boundary words**.

This gives three exact alternatives for an interface:

1. Preserve the 32-word relation using at least four component charts.
2. Retain a more expressive constraint language that enforces the missing
   interactions; one arbitrary linear-code object already does so.
3. Enlarge the relation, but explicitly account for all newly admitted words.

An incoming checkpoint, acc88b273c, supplies a sharper selector obstruction
in [marked geometries and selectors](../../04-computation/experiments/marked_geometries_selectors_20261018o.py).
Its finite-geometry maps retain their separate labels and are not
unit-distance constructions. The relevant elementary obstruction is
majority(0,170,204)=136 outside C; our checker independently replays it.

**Corollary 2 (PROVED).** A binary linear relation C is the exact boundary
support of a finite graph with fixed colour lists of size at most two,
allowing hidden vertices and pins, iff kappa(C)=0. At each boundary port,
the bits name a fixed ordered pair of distinct actual colours. Lists, edges
and pins must be fixed independently of the requested boundary word.

For necessity, take three proper list-colourings. At every vertex a majority
colour exists because its list has size at most two. If adjacent vertices
had the same majority colour, each would have that colour in at least two
of the three colourings. Those two index sets intersect, contradicting
properness of one original colouring. Pins are automatically retained.
Thus majority of full extensions is a full extension, proving closure of
the projected boundary support. Since 0 is in C, majority(0,x,y)=x AND y
forces Boolean-ring closure, equivalent to kappa=0. Sufficiency follows
from the component-graph construction with a common two-colour palette
and the required pins. Equivalently, the fixed-list edge constraints are
2-CNF clauses, and majority preserves each clause.

The proof also works for any odd number of colourings. It recovers the same
intersection-of-majorities mechanism as the previous
[odd-fibre preservation result](odd_prime_function_cycle_20261018d.md),
now applied to the boundary support of a fixed-list network. It preserves
existence, not the number of extensions.

Therefore neither binary controller vertices nor arbitrarily large networks
of overlapping two-colour lists can realize our four-chart union. Somewhere
a genuinely larger colour domain or a more expressive constraint is needed.
The incoming checkpoint provides a useful local candidate: a triangle whose
vertices each have the fixed list {0,1,2}, with three pendant binary ports,
implements not-all-equal on those ports. But with unrestricted six colours
the same graph accepts all eight binary port tuples. The checker records all
extension counts. Thus forcing the list and composing at the actual ports
are essential unpaid obligations, not a completed geometric selector.

The all-order formula and concrete carry distinguish faithful recursion
from apparently cheaper compression that weakens the target predicate.
The geometric next obligation is now a finite, authenticated selector:
realize the union of the four named eight-word cells while retaining
the marked ports, and exclude every additional separating six-colouring.
Four separate bipartite witnesses prove individual feasibility, not that
their union is realizable by a coupled unit-distance graph. Identifying
their ports enforces intersection; unrestricted ternary labels introduce
additional repairs. Both issues remain OPEN, as do Hadwiger–Nelson, Collatz,
planar Jacobian and LRC(14). No inference to any of them is made from these
finite encodings.

The same incoming checkpoint's actual ROOT-clock work retains integer
arrival clocks before comparing Fourier channels. We use only its transfer
lesson: a channel average, a chart union and a source-specific response are
different consumers. Its finite clock bounds do not supply the missing
selector, and our chart invariant supplies no Collatz energy bound.

For the “least resistance” phrase's other recovered meaning, see
[strong tournament / golden rigidity](strong_tournament_balance_golden_rigidity_20261014.md),
Section 3: series-parallel unit-resistor ladders give Fibonacci ratios.
Modulo 5 the Fibonacci matrix Q has full state/Pisano order 20; the ladder
rung Q^2 has full state order 10 and projective order 5 (MISTAKE-615).
Those clocks and quotients must be distinguished. This is a useful parallel
to our 32 source states versus four selectors, but the recurrences and cost
functions are different. The valid shared research move is to
identify the operation, retained state, quotient and decoder before
transferring an invariant.

## 9. Exact reproduction and audit boundary

~~~
python -B -X utf8 04-computation/experiments/gray_selector_recursion_20261018e.py
python -B -O -X utf8 04-computation/experiments/gray_selector_recursion_20261018e.py --verify
~~~

The standard-library checker exports the
[JSON receipts](gray_selector_recursion_20261018e.json).
The finite universe is every binary linear code on 0 through 6 named
coordinates: counts 1,2,5,16,67,374,2825, totaling 3290. For all 91 codes
through length 4, an independent literal minimum-set-cover solver enumerates
all eligible charts and confirms Proposition 1. Larger cases authenticate
the attained coset partition, packing and Boolean closure criteria; the
universal result follows from the proof, not extrapolation.

It also checks all 32 source parameter triples, all 1024 quaternary-code
addition pairs with carry, full graph colour assignments for all four charts,
the odd-cycle and missing-pin hostiles, the blue fixed part and its unique
self-dual subcode, K3/K4/K5 Eulerian codes, all pairs of length-three codes
under direct sum, and RM(1,m)+I2^q for 1<=m<=6, 0<=q<=3. It retains both
successful maps and the missing-carry/rotation/AND counterexamples.
Every binary relation on two bits is checked for majority preservation;
the listed and unrestricted triangle gadgets are enumerated separately.

The completed packet has 31,317 generation checks; optimized replay
adds one artifact-equality check. Separate agents independently audited the
algebra, all 32 Gray parameters, the mirror construction, tournament codec
and all small codes through length five. No test uses optimization-disabled
assertions. Named selectors, gauges, pins and the full boundary relation
remain part of every claimed transport.
