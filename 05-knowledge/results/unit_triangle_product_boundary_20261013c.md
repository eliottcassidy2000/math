# A colour-preserving geometric repair of the ternary frame

2026-10-13. **PROVED:** the all-prime construction, its complete induced
incidence, the labelled three-colouring equivalence, and the sharp
three-factor boundary below. **FINITE-EXACT:** the 81-point algebraic model
and declared support controls. **OPEN:** a comparable repair preserving the
needed six-colour boundary relation. This packet constructs three-chromatic
graphs, not a new four-, six- or seven-chromatic unit-distance graph. No
literature priority claim is made.

## 1. Inheritance and the changed operation

The closest mechanism is the
[finite abelian Cayley realization criterion](unit_cayley_realization_20261013a.md).
Its obstruction preserves every commutative square and every closing group
relation. The ternary frame is a useful hostile: its low chromatic number
does not make that complete labelled graph geometrically realizable. Here
we deliberately delete most of one direction family's edges, preserving
the **entire labelled three-colouring relation**, and then realize the
remaining graph with all its actual unit incidences.

The current
[THM-4628, every conjugation-fixed prime](../../01-canon/theorems/THM-4628-residue-colourings-of-the-unit-distance-graph-of-q-zeta55-are-kloosterman-forced-at-every-conjugation-fixed-prime.md)
and [composite residue extension](tree_resonance_composite_pell_20261013b.md)
close the fixed-prime residue-colouring route described there. We do not
replace their corrected quantifiers with a one-prime exception, or reverse
a residue graph homomorphism into a Euclidean embedding. Our construction
is a different, explicit graph over $\mathbb Q(\zeta_{3p})$.

The [marked fish cover](fish_marked_elliptic_cover_20261013.md) retains
simultaneous square conditions, not just individual directions. That lesson
is used here concretely: retain original vertex addresses, exact distances,
and complete boundary-colouring support together. Neither its unresolved
full-field fish problem nor its genus claims are premises of our proof.

The concept board is: commutative squares, geometric incidence, colouring
phase, marked boundary, and palette size. The anchor is a valid replacement
interface; the niche is an all-prime trace test; the wildcard is the least
number of triangle factors admitting an extra unit incidence. No identity
with the Conway or Hofstadter recurrence is asserted.

## 2. An all-prime exact incidence theorem

Fix a prime $p\ge5$, put $r=p-1$, and let
$\omega=e^{2\pi i/3}$ and $\xi=e^{2\pi i/p}$. Use the three unit-triangle
vertices

$$b_0=0,\qquad b_1=1,\qquad b_2=1+\omega.$$

Their six ordered nonzero differences are exactly $\mu_6$, the sixth roots
of unity. Retain the source address set $V=\mathbb F_3^r$ and place it by

$$z(j)=\sum_{i=0}^{r-1}\xi^i b_{j_i}\in\mathbb Q(\zeta_{3p}). \tag{1}$$

Let $H_r=C_3\mathbin\square\cdots\mathbin\square C_3$ be the graph on these
addresses whose edges differ in exactly one coordinate. Let $B_r$ add to
$H_r$ the three edges on the constant addresses $\mathbf0,\mathbf1,\mathbf2$.

**Incidence theorem.** Map (1) is injective and its full induced unit graph
is exactly $B_r$. Consequently it has $3^r$ vertices and $r3^r+3$ edges.

**Proof.** The relative degree of $\mathbb Q(\omega,\xi)$ over
$\mathbb Q(\omega)$ is $p-1$. For example, the Eisenstein proof for
$\Phi_p(X+1)$ works at any prime over $p$ in $\mathbb Q(\omega)$, since
$p\ne3$ is unramified there. Thus
$\operatorname{Tr}(1)=p-1$ and
$\operatorname{Tr}(\xi^a)=-1$ for $a\not\equiv0\pmod p$.

Write a difference of two placed addresses as
$d=\sum_i\delta_i\xi^i$, with
$\delta_i\in\{0\}\cup\mu_6$, and let $m$ be the number of nonzero
coefficients. The relative trace identity is

$$\operatorname{Tr}(d\bar d)
  =pm-\left|\sum_i\delta_i\right|^2.                       \tag{2}$$

If $d=0$, equation (2) and the triangle inequality give
$pm\le m^2$. Since $m\le p-1$, this forces $m=0$, proving injectivity.
If $|d|=1$, they give

$$pm-(p-1)\le m^2,
  \qquad (m-1)(m-(p-1))\ge0.                              \tag{3}$$

For $1\le m\le p-1$, only $m=1$ and $m=p-1$ remain. In the latter case
equality in the triangle inequality forces all $\delta_i$ to be identical.
Conversely, support one has unit length, and equal full support gives
$d=-\delta\xi^{p-1}$, also of unit length. A fixed ordered triangle
difference uniquely specifies its ordered pair of endpoints. Thus equal
full support is precisely a pair of constant addresses. This proves the
complete incidence statement, including exclusion of coincidences and
of every other possible cross edge. It is not just a lower-edge check.

## 3. Complete labelled three-colouring support survives

Regard the labelled palette as $\mathbb F_3$. Every proper three-colouring
of $H_r$ has the unique form

$$c(j)=c_0+\sum_{i=0}^{r-1}\varepsilon_i j_i,
  \qquad c_0\in\mathbb F_3,\quad\varepsilon_i\in\{1,-1\}. \tag{4}$$

Indeed, each coordinate line is a triangle, hence a permutation of the
three colours and therefore an affine function with slope $1$ or $-1$.
On a two-coordinate face, the three row permutations are pairwise disjoint
in each column. After normalizing the first row, the two others are the two
derangements of three symbols, which are its cyclic translates. Thus the
slope in one coordinate cannot change with any other coordinate. Connecting
such faces proves (4) globally.

Let $F_r$ be the full Cayley frame obtained by adding all edges of direction
$\pm\mathbf1$ to $H_r$. There are $3^r$ such additional edges. Every one
has the same colour condition in (4):

$$\sum_i\varepsilon_i\not\equiv0\pmod3.                  \tag{5}$$

Even the one chord $\{\mathbf0,\mathbf1\}$ enforces (5). Therefore, on
the **same labelled vertex set**,

$$\operatorname{Col}_3(F_r)
 =\operatorname{Col}_3(B_r)
 =\operatorname{Col}_3(H_r+\{\mathbf0\mathbf1\}).        \tag{6}$$

This is equality of colourings, not an inference from equal counts. Since
$r$ is even, the number of sign vectors with zero sum modulo three is

$$\frac{2^r+2}{3}.$$

For completeness, a root-of-unity filter gives
$[2^r+(\omega+\omega^{-1})^r+(\omega^2+\omega^{-2})^r]/3$,
which is the displayed expression because both parenthesized sums are
$-1$. Thus all three graphs in (6) have exactly

$$3\left(2^r-\frac{2^r+2}{3}\right)=2^{r+1}-2             \tag{7}$$

labelled proper three-colourings. They have triangles and at least one
such colouring, so their chromatic number is exactly three.

The full frame $F_r$ has $r+1$ cyclic direction subgroups in
$\mathbb F_3^r$, so the inherited Cayley theorem forbids its faithful planar
unit realization. Deleting $3^r-3$ of its edges gives the explicit induced
realization $B_r$ without changing (6). No minimality among arbitrary
edge-deletion repairs is claimed. For $r\ge4$, the three constant vertices
have degree $2r+2$ and every other vertex degree $2r$, so $B_r$ is not a
Cayley graph: it has escaped the original symmetric category.

### Exact boundary replacement and its palette limit

For any chosen subset of original vertices, equation (6) preserves the
boundary assignments, their equality partitions, their exact extension
counts, and their colouring lifts. If an abstract external graph is glued
along any of these original labels, the two complete colouring sets are
intersected with the same extra constraints. Hence the replacement still
preserves the entire labelled three-colouring relation in that context.
This is a concrete instance of the source-bound interfaces in
[boundary palette compression](boundary_palette_compression_20261011.md).

An abstract interface replacement does not provide coordinates for an
external geometric context: its retained ports must actually coincide and
all cross distances must be checked. A support packet alone authenticates
an abstract graph; the point tuple and complete incidence audit authenticate
the geometric source separately.

The palette restriction is sharp. Choose any deleted frame edge $\{u,v\}$.
It is a nonedge of $B_r$. In a proper three-colouring of $B_r$, recolour both
$u,v$ with one fresh fourth colour. The result remains proper on $B_r$ but
violates the restored edge in $F_r$. Thus (6) fails already with four
available colours, and in particular supplies no six-colour equivalence.
By the current
[THM-4614, three-colourable exterior contexts](../../01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md),
the known six-colour target also requires information this three-colour
construction does not supply. We do not transfer that theorem to different
fields or port configurations; the explicit fourth-colour hostile already
proves the limitation here.

## 4. Four factors are the first possible extra-incidence depth

**Sharp factor boundary.** Every faithful planar unit realization of
$H_r$ is induced for $r\le3$. The construction at $p=5$, $r=4$ shows this
fails at four factors.

First, each commutative coordinate square has four distinct unit-edge
vertices, so it is a rhombus. The opposite increments agree. Propagating
these identities separates any faithful realization additively into one
equilateral triangle per coordinate, up to an overall translation. Each
factor's six nonzero differences therefore forms a rotated copy of $\mu_6$.

A putative extra unit edge differing in two coordinates has unit increments
$a,b$ with $|a+b|=1$. Their angle is $120$ degrees, so the two factor
direction sets coincide. Choosing opposite differences then makes two
different source addresses coincide, contradicting faithfulness.

For three coordinates, suppose $a+b+c=d$ with all four vectors unit.
If $a+b=0$, a cancelling pair already contradicts faithfulness. Otherwise,
the unordered pair of unit-circle points having the nonzero midpoint
$(a+b)/2$ is unique. From $a+b=d-c$ it follows that
$\{a,b\}=\{d,-c\}$. Again two increments from distinct factors cancel.
This excludes every extra edge for $r\le3$.

Thus retaining a few squares does not imply that all higher cross incidences
are harmless. At four factors the trace identity gives a deliberate extra
triangle, and its colouring effect is evaluated with the complete support,
not inferred from pairwise direction independence.

## 5. Exact finite object, APIs and reproduction

The concrete implementation uses the integral power basis of
$\mathbb Q(\zeta_{15})$, with
$\Phi_{15}(X)=X^8-X^7+X^5-X^4+X^3-X+1$. The original labels are all
81 tuples in lexicographic $\{0,1,2\}^4$ order.

| Object | Exact count |
|---|---:|
| Triangle-product vertices | 81 |
| Product edges | 324 |
| Complete frame edges | 405 |
| Actual induced unit edges | 327 |
| Deleted frame edges | 78 |
| Product labelled three-colourings | 48 |
| Repaired and full-frame labelled three-colourings | 30 |
| Difference words in $\{0\}\cup\mu_6$ | 2401 |
| Unit difference words | 30 |

The three diagonal ports are original indices $(0,40,80)$. Before the added
triangle they have 18 all-equal and 30 all-distinct labelled extensions;
afterward exactly the latter 30 survive. The implementation retains every
valid sign vector and palette origin. `compile_support` constructs the
complete relation, `audit_support` recomputes it, `lift` checks the literal
source graph, and `boundary_counts` retains extension multiplicities.
`colour_safe_edges` means true in **every** retained colouring. Its empty
support case is explicitly vacuous, not a certificate of colourability.

The geometric checker enumerates every original pair and validates exact
coefficient and vertex types before comparison. Its cyclotomic multiplication
and conjugation are independently checked against the inherited larger
cyclotomic representation. Complete backtracking three-colour enumeration
on the one- and two-factor cases is independent of formula (4). Additional
controls cover the fourth-colour hostile, missing phases, wrong types,
duplicate coordinates, selected boundary marginals, the inherited frame
collision, and the all-prime phase-count formula at $p=5,7,11,13$.
Large-prime graphs are covered by the proof, not by an unstated census.

[Script](../../04-computation/experiments/unit_triangle_product_boundary_20261013c.py)
and [saved output](unit_triangle_product_boundary_20261013c.out):

```text
python -B -X utf8 04-computation/experiments/unit_triangle_product_boundary_20261013c.py
python -B -O -X utf8 04-computation/experiments/unit_triangle_product_boundary_20261013c.py
```

All checks are active under optimization. Normal and optimized runs match
the saved output at **6,023 exact checks**. LF SHA256:
`9a7560a17fd7efc6fc278a13ab8815314e4c43e6ee610a812a7304231829d3ed`.

The next specific obligation is to find a geometrically realizable edge
repair preserving the required higher-palette port relation, or to certify
the exact loss and compensate it with a realizable external context. The
source graph, metric ports and full colour relation are now separate checked
inputs to that task; no current six- or seven-colour plane conclusion follows.
