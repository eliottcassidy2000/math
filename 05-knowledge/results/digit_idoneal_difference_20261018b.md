# Ordered digit fibres, idoneal class labels, and exact colour differences

2026-10-18. **PROVED:** the elementary carry, reconstruction, Cartesian-box
and marked action statements below. **FINITE-EXACT:** the declared small
graphs, nine-vertex tournament and original eight-port counts. **CITED:**
quadratic-order genus theory. No new idoneal completeness theorem, geometric
equality fragment, or six-colouring of the plane is claimed.

The useful compression is an ordered difference code with its cycle equations.
Replacing that code by independent sets of differences is usually not exact.
The nine digit words provide a concrete example of the phase and orientation
that must be retained; an idoneal class group provides a different example of
a quotient that really is complete, under its stated hypothesis.

## 1. Targeted recovery and the source types

The exact list

$$396,417,528,639,741,852,963,174,285$$

is recovered in §4 of
[the repaired Solfeggio atlas](preprint_connections_solfeggio_digits_20261017.md).
It consists of the three cyclic presentations of each fibre of
$\mathbb Z/9\to\mathbb Z/3$, using the printed digit9 for residue0. The
atlas's §5 matches to similarly numbered theorem IDs are explicitly labelled
analogy or numerology. They do not establish tournament origins.

There is nevertheless an authentic tournament source:
[THM-871, Fermat-rung rigidity](../../01-canon/theorems/THM-871-fermat-rung-rigidity.md)
identifies $C_9\{1,3,4,7\}$ with $C_3[C_3]$, with these exact three fibres.
Its corrected torsor statement distinguishes a presentation orbit from its
quotient set. We use the labelled nine-vertex construction, not the coincident
cardinalities in that theorem.

For idoneal numbers the current source is
[THM-4566, Hilbert-class-field planes](../../01-canon/theorems/THM-4566-hilbert-class-field-planes-are-the-idoneal-planes.md),
with the order/field distinction in
[golden tail and idoneal carriers, §6](golden_tail_idoneal_20261017c.md) and the
conductor calculation in
[escape and renewal, §§6–7](escape_renewal_idoneal_20261017d.md).
MISTAKE-626 preserves the two status layers: the classical literature's
qualified completeness statements, and this repository's adopted completeness
modulo its accepted openai/math #003 input. This packet neither re-audits that
input nor uses completeness of the65-element list. MISTAKE-628 also forbids
transporting the atlas's classical odd-preimage interpretation to a different
seed transport.

The concept board is **ordered lift / integer carry / named palette / cycle
holonomy / class quotient / marked automorphism**. The canonical hostile is
three individually possible differences that do not sum to zero around a
cycle. The least-used sidecar is the full joint relation between differences,
rather than their separate sets of possible values.

## 2. The digit list has a reversible all-depth carry operation

For depth $d\ge1$, put $M=3^{d-1}$ and write a canonical residue modulo$3M$ as

$$x=r+Mt,\qquad0\le r<M,\quad t\in\{0,1,2\}.$$

The affine map $B(x)=4x+1$ has the exact skew-coordinate rule

$$r'=(4r+1)\bmod M,\qquad
 t'=t+\left\lfloor\frac{4r+1}{M}\right\rfloor\pmod3. \tag{1}$$

This follows from $4Mt\equiv Mt\pmod{3M}$. It is invertible: recover
$r=4^{-1}(r'-1)\bmod M$, recompute the carry, and subtract it from$t'$.
At$M=1$ the lower coordinate is0. The API validates canonical residues, rather
than silently reducing a supplied wrong coordinate.

Moreover

$$B^{3^{d-1}}(x)\equiv x+3^{d-1}\pmod{3^d}. \tag{2}$$

Indeed $B^k(x)=4^kx+(4^k-1)/3$, and induction by cubing gives
$4^{3^{d-1}}\equiv1+3^d\pmod{3^{d+1}}$. This proof retains the extra
factor3 of precision before division. It is not valid to divide a residue
known only modulo$3^d$ by3.

At$d=2$, the word beginning at$a$ is

$$[a,a-3,a-6]\pmod9. \tag{3}$$

Printing0 as9 gives exactly the nine words above. They follow the **reverse**
of the$B^3$ generator in(2). Forgetting the starting phase gives only one of
three fibres; forgetting orientation identifies the chosen generator with its
inverse. Reading the word as an ordinary decimal integer gives yet another
object, with residue$3a\pmod9$. It does not preserve the ordered lift.

### The intrinsic tournament map

Orient $u\to v$ on$\mathbb Z/9$ iff $v-u\in\{1,3,4,7\}$. In coordinates
$u=r+3t,v=s+3z$, its rule is

$$u\to v\iff
\begin{cases}s-r=1\pmod3,&r\ne s,\\z-t=1\pmod3,&r=s.\end{cases} \tag{4}$$

Thus it is exactly three directed triangles arranged in a directed triangle.
The map$B$ preserves the connection set under multiplication by4, so it is an
actual automorphism, not a permutation action imposed on an unrelated census.
It is a nine-cycle, and$B^3$ rotates each fibre.

The fibre partition is intrinsic: a pair within one fibre belongs to exactly
one directed triangle, while a pair in different fibres belongs to exactly
three. Every automorphism therefore permutes the fibres cyclically and rotates
each inner triangle independently. The automorphism group has order
$3^3\cdot3=81$. The checker independently verifies this among all$9!$
label permutations. Formula(4), the direction convention, and the marked
starting vertex make this a lossless representation; the decimal magnitudes
have no role in its proof.

## 3. Difference codes retain all colourings when cycles are retained

Let$G=(V,E)$ be a finite connected labelled simple graph, choose an orientation
of every edge, a root$o$, and a named palette$\mathbb Z/q$, $q\ge2$.
For a colouring$c$, define

$$\delta_{uv}=c(v)-c(u)\pmod q. \tag{5}$$

**Exact reconstruction.** An edge assignment$\delta$ comes from a colouring
iff its signed sum around every cycle is zero. For a fixed root colour the
colouring is unique. It is proper iff every edge difference is nonzero.

Proof: sum(5) around cycles for necessity. For sufficiency, integrate along
the unique paths of a fixed spanning tree. Every non-tree edge is satisfied
exactly when its fundamental-cycle equation holds. The same argument proves
uniqueness. This is the elementary potential/gain description, not a new
general graph-colouring theorem.

Consequently only$n-1$ tree differences and one root colour are independent;
the$m-n+1$ non-tree edges are checked by an integer cycle matrix$C$. The
script's `tree_code` retains the literal graph and orientation, constructs$C$,
and `decode` returns the exact colouring or `None` for inconsistent gains.
`None` rejects that supplied assignment; it does not mean the graph is
uncolourable.

At$q=6$ the cycle equations split exactly modulo2 and3. With the canonical
component residues, the reconstructed colour is
$c=3c_2+4c_3\pmod6$. Properness is the coupled condition that the two
component differences are not both zero **on the same edge**. Separately
finding two component assignments is insufficient unless their edge supports
are compatible. Relabelling a named palette arbitrarily also need not preserve
arithmetic differences: swapping names1 and2 sends the equal differences of
$(0,1)$ and$(1,2)$ to2 and5. Root translation is harmless; an arbitrary colour
permutation requires its own retained gauge.

## 4. Sharp obstruction to replacing joint support by a Cartesian box

Let each edge$e$ have a nonempty allowed list$S_e\subseteq\mathbb Z/q$.
The lists may include0 for this theorem.

**Cartesian-box criterion.** Every selection from$\prod_e S_e$ integrates
to a colouring iff both conditions hold:

1. every edge lying on a cycle has a singleton list;
2. the resulting fixed gains have zero sums on cycles.

Bridge lists may vary arbitrarily. For necessity, change just one cyclic
edge between two distinct list values. A simple cycle through that edge
changes its sum by a nonzero signed amount, so both selections cannot
integrate. For sufficiency all cycles use only the fixed edges and the
reconstruction theorem applies. This is an iff for **every Cartesian
selection**; it is not a criterion for whether at least one selection works.

In particular, for$q\ge3$ the full nonzero edge box is integrable for every
choice exactly when the connected source is a tree. At$q=2$ the nonzero list
is the singleton$\{1\}$ instead, and the condition is precisely bipartiteness.
This palette boundary is why the binary and ternary parts cannot be treated
as interchangeable independent choices.

The minimal proper-colouring hostile is a triangle. Every oriented edge
individually permits all five nonzero$\mathbb Z/6$ differences. Choosing$+1$
around the directed cycle has sum3 and is impossible. With sorted edge order
$(01,02,12)$ the rejected tuple is$(1,5,1)$, whereas$(1,2,1)$ reconstructs
the proper colouring$(0,1,2)$.

### Repeated direction labels: a useful exact generalization

If several edges must share a signed variable$x_j$, substitute those
variables into the fundamental cycles to obtain an integer matrix$C$.
For nonempty lists$S_j$ choose any$a_j\in S_j$. Then

$$\prod_j S_j\subseteq\ker(C\bmod q)
\iff
Ca=0\pmod q\quad\hbox{and}\quad
C_{ij}(s-a_j)=0\pmod q\quad(\forall i,j,s\in S_j). \tag{6}$$

Necessity follows by varying a single coordinate; sufficiency follows by
adding the coordinate differences. Unlike an ordinary edge variable, a
repeated variable can have a nontrivial annihilator: at$q=6$, the constraint
$3x=0$ allows the whole list$\{2,4\}$. This example is an abstract labelled
gain system, not three equal Euclidean displacement vectors forming a
triangle. Equation(6) is a cheap complete test for a proposed *universal box*.
It does not solve arbitrary list-support existence without joint information.

### The retained eight-port source

For the actual old source$C_5\vee C_4$, with edges

$$01,04,05,07,12,23,34,56,67,$$

the independent cycle-polynomial count and the exact tree enumeration agree:

$$\#\{c:c(0)=0,\ c\text{ proper}\}
=\frac{(5^5-5)(5^4+5)}{6^2}=54\,600. \tag{7}$$

There are327,600 labelled six-colourings and570 equality patterns. All nine
individual edge-difference sets are$\{1,2,3,4,5\}$, whose Cartesian product
has1,953,125 elements. The tree code examines78,125 assignments and retains
exactly54,600. It recovers the same570 patterns as the original boundary
table, without identifying a marginal list with full support.

This gives a concrete consumer for shared exceptional phases: first retain
which constraints belong to the **same phase and source**, then integrate
their joint differences. It neither constructs13 new geometric pairs nor
turns their number into a sufficient equality certificate. The parent
exceptional-phase search still owes complete incidence and an actual
boundary-colouring obstruction.

**Incoming exact complementary compiler.** Section7, G1 of
[source orientation and boundary contexts](source_orientation_boundary_contexts_20261018k.md)
(read from `origin/main` before its integration) fixes an additive colouring
on Cartesian physical factors and takes the exact sets of differences of
independently selectable factor pairs. There a contact family is proper iff
zero is absent from the sum of those sets. This is a true iff because every
combination of factor pairs is an actual contact in that declared family.
Our box obstruction concerns the different operation of independently
selecting gains on edges of one graph and asking for a common colouring.
Those gains share vertices and owe cycle consistency. Neither result
contradicts the other. The distinction is independence of the **actual
source choices**, not the use or avoidance of a Cartesian data structure.

## 5. Idoneal genus labels: a complete quotient with a different operation

**CITED classical input.** For the proper class group$A=\operatorname{Cl}(-4n)$
of the indicated quadratic **order**, the genus map is$A\to A/A^2$.
It is injective iff$A$ has exponent at most2, the idoneal condition.
[Kani, *Idoneal Numbers and some Generalizations*, Theorem6](https://www.labmath.uqam.ca/~annales/volumes/35-2/PDF/197-227.pdf)
gives this equivalence. [Watkins's two tables](https://magma.maths.usyd.edu.au/~watkins/papers/IDONcorr.pdf)
distinguish the order discriminant from the underlying field and its
fundamental discriminant. These sources do not identify arbitrary65-element
sets by a canonical map. The repository's adopted additional completeness
input was already separated in§1; no step below uses it.

A small exact instance is$n=30$, discriminant$-120$. The complete reduced
primitive forms are

$$ (1,0,30),\ (2,0,15),\ (3,0,10),\ (5,0,6). $$

All are self-inverse, so this four-element class group is$\mathbb F_2^2$.
The two marked genus characters$\chi_8,\chi_{-3}$ label these forms by
$00,01,10,11$: represented units1,17,13,11 give the labels respectively.
The script independently enumerates the reduced forms and checks the
characters on represented120-units with$|x|,|y|\le8$. Their class invariance
is the cited genus theorem, not an inference from that finite box.

There is an exact, but **marked**, three-cycle on its nonprincipal classes:

$$Q(a,b)=(b,a+b),\qquad Q^3=I,\qquad
10\mapsto01\mapsto11\mapsto10. \tag{8}$$

Choosing one initial nonprincipal class conjugates the ordered three-lift
cycle to this automorphism. This preserves equality and the chosen cyclic
successor. It does not preserve addition from$\mathbb Z/3$, because every
homomorphism$\mathbb Z/3\to\mathbb F_2^r$ is zero. Nor can translation by
one binary class have order3. The operation in(8) is an **automorphism**, not
translation and not an asserted automorphism of a number field. Its inverse
reverses the associated directed triangle. An orientation and a character
basis are genuine sidecars, rather than canonical consequences of the number30.

The missing-coordinate hostile is already in the recovered route:
$\operatorname{Cl}(-44)$ and$\operatorname{Cl}(-76)$ each have three classes,
but trivial genus quotient. Their full reduced forms are respectively

$$ (1,0,11),(3,-2,4),(3,2,4),\qquad
(1,0,19),(4,-2,5),(4,2,5). $$

Binary genus data erase their order-three class direction completely.
Inverting the conductor also changes the order to its maximal-order
localization; an order-level obstruction cannot be imposed on unrestricted
field coordinates. This guards the intended geometric interpretation.

### A six-colour action from the four class labels

Identify the four marked labels of$\mathbb F_2^2$ with the vertices of$K_4$.
Its six unordered edges are an exact six-colour palette. An edge$\{x,y\}$
has a nonzero direction$d=x+y$; the two edges of each direction are opposite.
Give an edge shade0 when it contains the marked origin0, and shade1 otherwise.
Thus colour$2(d-1)+s$, using the binary labels$d=1,2,3$, recovers exactly a
three-pair palette with one binary shade in each pair.

**Marked palette lemma.** The24 permutations of the four vertices act on
these six colours as exactly the pair-preserving permutations with an even
number of shade flips. The full pair-preserving group has48 elements.
Its missing central involution is

$$\tau:\{x,y\}\longmapsto\mathbb F_2^2\setminus\{x,y\},
\qquad \tau(c)=c\mathbin{\mathrm{xor}}1. \tag{9}$$

Proof: every permutation of the four vectors is affine: after translating
the image of0 to0, its arbitrary permutation of the three nonzero vectors
is linear. A linear map permutes the three directions and preserves shade0.
Translation by$t\ne0$ preserves shade in direction$t$ and flips the other
two directions, since$\{t,t+d\}$ contains0 exactly when$d=t$. These are
all four even flip vectors, followed by all six direction permutations.
The opposite-edge involution flips all three shades, commutes with every
pair-preserving permutation, and doubles the subgroup to48. It cannot come
from a vertex permutation: it sends a triangle's three edges to a star's
three edges. Equivalently its flip parity is odd. The checker verifies all
24 vertex actions, all48 pair actions, and this centrality directly.

This is the actual operation used by the [112-point support hostile](shared_phase_threshold_hostile_20261018b.md):
each copied source may apply$\tau$ to its colours to repair its single spoke.
It is a palette action and can preserve every boundary colouring even when
no source-fixed graph retraction exists. The identification with
$\operatorname{Cl}(-120)$ uses the explicitly marked character basis above;
it imports neither a field automorphism nor a new geometric edge. Conversely,
restricting palette operations to vertex-induced actions would incorrectly
discard this valid odd-flip repair.

## 6. Reproduction and current boundary

The [script](../../04-computation/experiments/digit_idoneal_difference_20261018b.py)
and [JSON receipt](digit_idoneal_difference_20261018b.json) retain literal
graphs, ordered residues, cycle matrices, forms and character labels.
Public APIs reject Boolean/float aliases, noncanonical residues, malformed
or disconnected source graphs, and empty lists. No assertion is needed for
production validation or controls.

The finite universe is all residues through ternary depth7; every permutation
of the nine tournament vertices; all44 connected labelled graphs on at most4
vertices with root-fixed palettes2,3,4; all$7^4$ nonempty-list boxes for a
triangle with one bridge at palette3; all216 triangle gain assignments at
palette6; the complete78,125 tree assignments of the original eight-port
source; the stated reduced-form/represented-unit controls; and all24/48
marked$K_4$ vertex/pair palette actions.

```
python -B -X utf8 04-computation/experiments/digit_idoneal_difference_20261018b.py
python -B -O -X utf8 04-computation/experiments/digit_idoneal_difference_20261018b.py
```

Normal and optimized runs agree at **23,190 exact checks**. The LF-normalized
JSON SHA256 is
`2e6abe1c2203ac72bd3a834f6ec8e76cf06eb129a6ecab5870eb48a286ed90c7`.

Independent proof/API review and optimized replay agree with the saved
receipt, including the all-depth carry precision, joint cycle constraints,
marked class-group scope, and the 24/48-action distinction.

The next obligation is to use a *joint* difference relation whose support
actually excludes a separating boundary colouring at one authenticated
geometric phase. The digit automorphism, idoneal genus labels and full
per-direction difference marginals alone do not discharge it.
