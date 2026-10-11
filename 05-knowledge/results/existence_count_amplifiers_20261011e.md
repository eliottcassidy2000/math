# Existence versus count amplification at a marked two-port interface

**Status: PROVED elementary composition, support and selector statements;
FINITE-EXACT bounded controls. No six-colour unit-distance equality tester
is constructed.** The ports, palette and source graph are retained throughout.

The desired six-colour interface has a positive equal-endpoint extension
count and a zero different-endpoint extension count. This packet separates
that exact support target from making different endpoints unlikely. It gives
both an obstruction for a fixed composition grammar and a concrete operation
outside that grammar which creates an exact equality interface.

## 1. Inheritance and the retained observable

The closest inherited mechanism is the full boundary-weight table in
[boundary_context_weights_20261011d.md](boundary_context_weights_20261011d.md),
especially its conditional counts, named-colour orbit factors and singleton
context selectors. Its hidden-leaf/isolate example preserves Boolean
extension support while changing counts. The existing
[boundary_palette_compression_20261011.md](boundary_palette_compression_20261011.md)
retains actual source graphs, ports and colour witnesses. The current
[context_prime_root_20261011.md](context_prime_root_20261011.md), A3--A4,
similarly distinguishes an actual quantitative inequality from compatible
finite shadows; its arithmetic statements are not dependencies here.

The live objects are: two-port conditional counts, Boolean support,
serial/parallel composition, an interior constraint graph, and unit-distance
realizability. The canonical hostile is a positive integer count whose
normalized probability tends to zero. The corrected near miss is treating
the same ratio, such as $5/6$, as the same observable in different palettes
and contexts. The useful sidecar is the marked source graph: changing only
which two vertices are ports can reverse an equality interface into an
inequality interface.

Fix an integer $k\ge2$. A gadget is a finite simple graph $F$ with distinct
ordered ports $u,v$. Let $e$ count proper $k$-colourings extending the fixed
assignment $(c(u),c(v))=(0,0)$, and let $d$ count those extending $(0,1)$.
Colour permutations show that these counts apply to every fixed equal or
different ordered pair respectively. They are counts per fixed assignment,
not already multiplied by its orbit size. Thus

\[
 P_F(k)=ke+k(k-1)d,\qquad
 p_D=\frac{(k-1)d}{e+(k-1)d}.                         \tag{1}
\]

The probability is defined only when $F$ is $k$-colourable. The equality
target is exactly $e>0,d=0$; $(0,0)$ is an uncolourable graph, not a tester.

| Source and target | Exact map and preserved predicate | Lost information / required sidecar | Cheapest test |
|---|---|---|---|
| Marked graph to two counts | Fixed equal/different assignments; exact extension support | Counts do not recover its edges, coordinates or literal colouring | Retain graph and ports; independently enumerate a small source |
| Two gadgets to a composite | Private interiors, specified port identifications | A graph operation need not have a simultaneous unit embedding | Authenticate vertex maps; check geometry separately |
| Counts to probability | Equation (1) | Absolute count scale is discarded | $K_{2,r}$ has $p_D\to0$ while $d>0$ |
| Interior graph to selector | Join each port to every interior vertex | Unit-separated ports cannot host an arbitrary common-neighbour graph | One interior edge gives the abstract three-colour diamond |

## 2. The exact two-coordinate algebra

Parallel composition identifies both pairs of ports and otherwise takes
disjoint interiors, with no extra edges. Series composition identifies the
second port of the first gadget with the first port of the second and
retains the two outer ports. For count pairs $(e,d)$ and $(f,g)$:

\[
\begin{aligned}
 (e,d)\mathbin{\|}(f,g)&=(ef,dg),\\
 (e,d)\mathbin{\circ}(f,g)&=
 \bigl(ef+(k-1)dg,\;eg+df+(k-2)dg\bigr).             \tag{2}
\end{aligned}
\]

For equal outer colours, the shared series vertex either has that colour
or one of the other $k-1$ colours. For different outer colours it may have
the first, the second, or one of the $k-2$ remaining colours. This proves
(2) directly. Parallel counts multiply because the interiors are independent
after fixing the boundary assignment. Duplicate boundary edges do not cause
a problem: the same zero/one boundary restriction is idempotent.

The full named-colour kernel is $M=eI+d(J-I)$. On functions constant away
from a specified colour its two-state matrix is

\[
 A(e,d)=\begin{pmatrix}e&(k-1)d\\d&e+(k-2)d\end{pmatrix}.
\]

Series composition is ordinary matrix multiplication. Parallel composition
is entrywise multiplication of the **full** $k$-state kernels, or componentwise
multiplication of $(e,d)$; it is not entrywise multiplication of these
two-state matrices. The two eigenvalues are $e+(k-1)d$ and $e-d$. Consequently

\[
 \theta(e,d)=\frac{e-d}{e+(k-1)d},\qquad
 \theta(a\circ b)=\theta(a)\theta(b).                \tag{3}
\]

For $d,g>0$, write $R_1=e/d,R_2=f/g$. Then parallel composition multiplies
$R_1R_2$, while series gives

\[
 R_{\rm series}=\frac{R_1R_2+k-1}{R_1+R_2+k-2}.
\]

If $R_1,R_2\ge1$, this lies between $1$ and
$\min(R_1,R_2)$. Subtracting $R_1$ has numerator
$-(R_1-1)(R_1+k-1)$, and symmetrically for $R_2$.
Thus series composition does not improve either already equality-favouring
ratio; parallel composition amplifies it, subject to its unchanged support.

**Finite positive-support obstruction.** For $k\ge3$, every finite
serial/parallel expression whose leaves all have $d>0$ still has $d>0$.
Parallel multiplication is positive; in series the term $(k-2)dg$ is
positive. Induction proves the claim for every finite expression, regardless
of its size or intermediate count ratios. In particular, a grammar starting
only from edges cannot construct an equality-only tester at these palettes.
This is a marked two-terminal grammar restriction, not a theorem forbidding
all graph gadgets or all overlapping gluing operations.

The boundary $k=2$ is real: an edge has $(e,d)=(0,1)$; two edges in series
have $(1,0)$. In Boolean terms the series different-support condition has
the third disjunct $D_1\wedge D_2$ precisely when $k\ge3$.

## 3. Amplification, its geometric cap, and the absolute-count obligation

For $r\ge0$, give the two nonadjacent ports $r$ independent common neighbours.
The resulting graph $K_{2,r}$ has

\[
 e=(k-1)^r,\qquad d=(k-2)^r,\qquad
 \frac de=\left(\frac{k-2}{k-1}\right)^r.             \tag{4}
\]

For each $k\ge3$ this makes $p_D\to0$, but every finite $r$ has $d>0$.
At $k=3$, $d=1$; at $k\ge4$, $d$ itself grows. Thus convergence of
normalized probabilities is not a finite existence certificate.
The single-neighbour ratio is $4/5$ at palette six, $5/6$ at palette seven,
and $6/7$ at palette eight. In contrast, attaching a leaf rather than an
isolated vertex to a fixed existing graph changes **both** channels by the
common ratio $(k-1)/k$. That $5/6$ at palette six cancels from $p_D$.

If the two ports are at unit separation in the Euclidean plane, their common
unit neighbours are in the intersection of two unit circles: there are at
most two. When both exist, their separation is $\sqrt3$, so they cannot
carry an additional unit edge between them. This caps the particular
common-neighbour-only bank at $r\le2$. At palette six its probabilities are

\[
 p_D(r=0)=5/6,\quad p_D(r=1)=4/5,\quad p_D(r=2)=16/21.
\]

The last value is a lower bound only for this bank, not for general
unit-distance gadgets with further vertices.

There is a valid conditional way to turn a probability estimate into an
exact support conclusion. Put $Z=e+(k-1)d$. If separately certified bounds
$Z\le U$ and $p_D\le\eta$ satisfy

\[
 U\eta<k-1,                                           \tag{5}
\]

then $(k-1)d=Zp_D<k-1$, hence the nonnegative integer $d$ is zero.
Both premises must concern the same source and palette. Equality in (5)
is insufficient: the three-colour $K_{2,2}$ has $Z=6,p_D=1/3,d=1$.
Using the truthful mass and probability from (4) never satisfies the
strict threshold. The program's `integer_zero_test` checks only implication
(5); it does not authenticate caller-supplied analytic bounds. A separately
proved larger divisibility quantum can strengthen the integer threshold;
that is an additional premise, not a consequence of a small probability.

## 4. A positive selector: connect the interiors

Let $J$ be any finite simple graph. Form $F_J$ by adjoining two nonadjacent
ports, each adjacent to every vertex of $J$. The interior has $k-1$ available
colours under an equal assignment and $k-2$ under a different assignment.
Therefore

\[
 e=P_J(k-1),\qquad d=P_J(k-2).                         \tag{6}
\]

It follows exactly that $F_J$ is an equality tester at palette $k$ iff
$\chi(J)=k-1$. Here an empty interior has chromatic number zero and gives
$(1,1)$, not an equality tester. This is a constructive sufficient system
with a decidable finite obligation on a supplied $J$.

The minimal instructive comparison at palette three is:

* two independent interior vertices give $K_{2,2}$ and $(e,d)=(4,1)$;
* adding their one mutual edge gives the diamond and $(e,d)=(2,0)$.

The extra edge couples the previously independent interiors and escapes
the positive-product grammar. The marked terminal role matters: the same
diamond with its two adjacent degree-three vertices as ports has $(0,1)$.
Indeed, it is obtainable as an edge in parallel with two length-two paths
between these adjacent ports. The obstruction in §2 does not concern the
unmarked graph alone.

Taking $J=K_{k-1}$ gives the abstract selector $((k-1)!,0)$, including
$(120,0)$ at palette six. It does not supply a unit-separated plane gadget:
its $k-1$ common neighbours violate the two-circle bound when $k\ge4$;
even the two adjacent neighbours of the diamond cannot both be unit from
unit-separated ports. The diamond instead admits a unit realization whose
nonadjacent ports are separated by $\sqrt3$.

There are actual unit-separated positive controls at palette three. The
inherited exact Moser spindle constructor produces an eleven-edge unit
graph. Deleting any one edge and retaining its endpoints as ports gives
$(e,d)=(4,0)$; the program checks every case and exports an actual equal-port
colouring. The deleted edge is absent from the source graph while its
endpoints remain at unit distance. Thus these are unit-distance subgraphs,
not induced unit-distance graphs. This known small positive control does
not extrapolate to a six-colour construction.

## 5. Every genuine equality tester has a high-degree core

Suppose $F$ has $e>0,d=0$ at palette $k$, with distinct ports $u,v$.
There is no edge $uv$ in $F$, since an equal-port colouring exists. The graph
$G=F+uv$ is not $k$-colourable; giving $u$ a fresh colour in any proper
$k$-colouring of $F$ shows that it is $(k+1)$-colourable. Choose a
vertex-critical $(k+1)$-chromatic subgraph $H$ of $G$. It must include the
added edge: otherwise it is a subgraph of the $k$-colourable graph $F$.

Every vertex of $H$ has degree at least $k$, since a $k$-colouring of its
vertex deletion would extend over a vertex of degree below $k$. Hence the
two ports have degree at least $k-1$ already in $F$, and at least $k-1$
other vertices have degree at least $k$ in $F$. At palette six this means
ports of degree at least five and at least five interior vertices of degree
at least six. It is a necessary filter, not sufficient, and does not claim
that irrelevant pendant vertices elsewhere in $F$ have high degree.

This identifies a concrete next target: keep a unit-distance realization,
construct interior dependence surviving these degree requirements, and
authenticate an equal witness together with zero different support (or a
truthful strict quantitative bound such as (5)). Increasing the repetition
count of the private common-neighbour bank cannot meet that target.

## 6. Implementation and reproduction

The owned experiment is
[existence_count_amplifiers_20261011e.py](../../04-computation/experiments/existence_count_amplifiers_20261011e.py).
`Counts` retains the palette and exact nonnegative integer counts;
`TwoPort` retains the actual graph and named ports. `compile_tree` validates
and literally counts each supplied leaf, fresh-copies interiors, and returns
the resulting graph and exact algebraic pair. `witness` performs finite
source-colouring enumeration, with no efficiency claim. No count pair alone
is accepted as a source-graph certificate. Public inputs are checked before
equality can exploit Boolean/integer or float/integer aliases.

The declared universe is: 1,792 pair-algebra cases with palettes 2--8 and
each input count 0--3; all 51 edge expressions with one through four leaves
at palettes 2--6 (255 cases); all 76 interior graphs on zero through four
labelled vertices at palettes 2--6 (380 cases); exact fan and leaf controls;
all eleven Moser edge deletions; the necessary-degree filter on all graphs
through order four and every nonadjacent port pair at palettes 2--4;
fourteen malformed public inputs. The direct named-colour sum and literal
source enumeration are independent of the transfer-matrix calculation.

```powershell
python -B -X utf8 04-computation/experiments/existence_count_amplifiers_20261011e.py
python -B -O -X utf8 04-computation/experiments/existence_count_amplifiers_20261011e.py
```

Both modes produce 11,443 always-active checks. The saved output is
[existence_count_amplifiers_20261011e.out](existence_count_amplifiers_20261011e.out).
Its LF-normalized SHA256 is
`4efe4828ec4dc426d98b1e2355a0573ed26fa5896d164217350850691a2c1f91`.
These finite controls accompany the all-finite-expression proofs; they are
not a census of possible unit-distance equality gadgets.
