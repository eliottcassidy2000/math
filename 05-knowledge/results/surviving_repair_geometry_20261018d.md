# A genuine cubic insertion that retracts onto the old graph

2026-10-18. **PROVED:** the odd-degree one-neighbour lemma, a complete
cubic-field contact reader, and the retraction consequence of the contact
certificate below. **FINITE-EXACT:** four bounded geometric recipes, their
complete induced edges, and all twenty-one field-contact decisions for the
specified fractional phase. **OPEN:** a six-colour equality fragment.

The main candidate has **1,482 points and 6,118 unit edges**, retaining the
old 714 points literally. All 768 new points genuinely use the cubic
extension. With the eight ports fixed, 470 new vertices and 334 old-to-new
edges survive degree-five deletion. Nevertheless, an explicit graph
retraction fixes every old vertex. Consequently **every old colouring
extends**, not merely the 570 saved source patterns.

The decisive obstruction is a complete field-contact calculation: the six
fractional branches have no possible unit contacts with one another, even
if their finite old-field base sets are enlarged. Each branch contacts the
old field only by a unit translation that survives forgetting its shift.

## 1. Inheritance and the new test

The source is the complete 714-point receipt in
[coupled fractional geometry](coupled_fractional_geometry_20261018c.md).
It retained 94 extra factor edges in its relative six-core, but
[coupled support branches](coupled_support_branches_20261018c.md) retained
all 570 source patterns. The older 2,646-point Cartesian-core example is
another reason not to equate a large core with lost boundary support.

The corrected near miss is separate: the fixed-anchor spindle triangle
already defeats anchor-only two-colour-pair charts. Its failure is not
caused by the 94 mixed edges, and does not exclude full six-colourings.
The present test therefore fixes each complete old colour word, rather
than testing only another restricted palette ansatz.

The concept board is **literal old source / true cubic coordinates / exact
cross contacts / relative core / retraction / full colouring extension**.
The anchor is a genuine cubic insertion; the niche is a complete
old-field neighbour reader; the cheap hostile is a cubic unit matching
which breaks one projection while preserving all colourings.

We use the current field conventions of
[THM-4700, unit-path ring](../../01-canon/theorems/THM-4700-the-unit-path-ring-of-the-six-colour-field-short-reciprocal-walks-and-the-heron-census-of-coefficient-triangles.md).
The triangle-free exterior bound in
[THM-4661, odd-degree carriers](../../01-canon/theorems/THM-4661-odd-degree-carriers-a-failed-six-colour-exterior-has-at-least-twenty-vertices-and-equality-testers-are-seven-colour-certificates.md)
does **not** apply here: our old field contains $\sqrt{-3}$. Its
one-old-neighbour mechanism does apply and is proved directly below.

## 2. Fixed field, phase, and literal source

Let

$$K=\mathbb Q(\zeta_{15},j),\qquad j^2=-11,
\qquad F=K(w),\quad w^3+2w^2-4w-4=0.$$

Here $[K:\mathbb Q]=16$, $[F:K]=3$, and conjugation fixes $w$ while
sending $j$ to $-j$ and $\zeta_{15}$ to its inverse. The cubic polynomial
has no rational root. It remains irreducible over $K$: a factorization of
a cubic would give a root of degree three inside a degree-sixteen field.

Use the exact unit

$$u=\frac{w+3-j}{w+3+j},\qquad
w=j\frac{1+u}{1-u}-3. \tag{1}$$

Thus $K(u)=F$; its cubic coefficient is not merely an unused ambient
coordinate. Let $C$ be the eight literal source points, $M$ the spindle,
$D=\{0,1,\eta\}$, and $S$ the six old fractional shifts, all exactly as
in the inherited receipt. Put

$$R=C+(M+D)\subset K,\qquad A=R+S.$$

The sets have $|R|=128$ and $|A|=714$. The JSON preserves the complete
coordinates of $A$ as vertices $0,\ldots,713$ and its 2,914 edges.
The first eight labels remain the original ports. Every sum is a literal
set of points; all address coincidences are identified.

Write $\delta=\zeta_{15}^5$, $\eta=1+\delta$ and

$$q=\delta\frac{-1+\sqrt{-15}}4,\qquad z=uq,
\qquad s_i=\eta^i\quad(0\leq i<6).$$

All these phases and directions have norm one. Define

$$T_s(z)=\frac{2z-s}{2-\bar s z},\qquad
t_s=2(T_s(z)-z),\qquad
G=A\cup\bigcup_{i=0}^{5}(R+t_{s_i}). \tag{2}$$

The denominator is nonzero because $|\bar s z|=1<2$; the code also checks
it exactly. These are six different exterior cosets modulo $K$, each
containing 128 new points. Equation (2) is the main 1,482-point candidate.
The common phase $z=uq$ is fixed throughout; no unrelated phase is
substituted at an individual contact.

## 3. An exact field-contact reader

### Odd-degree one-neighbour lemma

Let $K\subset\mathbb C$ be conjugation-stable and $L/K$ a finite extension
of odd degree. An exterior point $x\in L\setminus K$ has at most one
unit-distance neighbour in $K$.

Indeed, two distinct neighbours $a,b\in K$ give, by subtracting their
norm equations, $\bar x=\lambda x+\mu$ with coefficients in $K$ and
$\lambda\ne0$. Substitution into $(x-a)(\bar x-\bar a)=1$ gives a
nonzero quadratic over $K$. But the minimal degree of $x$ divides
$[L:K]$ and is odd and greater than one. This is impossible. Conjugation
stability of $L$ is not needed for this conclusion; it is available for
our explicit cubic algorithm.

In particular, between any two distinct additive $K$-cosets in $F$,
unit edges are a translation matching: there is at most one old-field
difference that can make an edge. This statement retains the coset labels
and the actual translation. It is not a bound on the number of points
in a coset.

### Complete cubic decision

For $x\in F\setminus K$, expand in the retained basis:

$$x=x_0+x_1w+x_2w^2,\quad
\bar x=y_0+y_1w+y_2w^2,\quad
x\bar x=N_0+N_1w+N_2w^2.$$

A neighbour $a\in K$ must satisfy

$$y_1a+x_1\bar a=N_1,\qquad
y_2a+x_2\bar a=N_2. \tag{3}$$

If $\Delta=y_1x_2-y_2x_1\ne0$, solve this two-by-two system in $K$.
The candidate is valid exactly when the second solution equals the
conjugate of the first and $N(x-a)=1$. Both checks are necessary and
sufficient; no finite-place hit is accepted as a geometric edge.

If $\Delta=0$, there is no neighbour. The two nonconstant coefficient
vectors of $x$ and $\bar x$ are nonzero and proportional, so
$\bar x=\lambda x+\mu$ with $\lambda\ne0$. A neighbour again gives
a nonzero quadratic over $K$, contradicting degree three. Thus returning
`None` in this branch is a complete decision, not a failed linear solve.

The public `k_neighbour` reader validates exact field coefficients and
requires an exterior input. Positive controls use $x=a+u$ and other
marked unit translates; the real cubic generator $w$ tests the singular
determinant branch. Old-field and floating-point aliases are rejected.

## 4. The complete six-branch contact certificate

For every unit $z,s$, direct clearing of denominators gives the retained
multiplicative phase identity

$$2(T_s(z)-z)+s=\bar s\,z\,T_s(z),\qquad
N\bigl(2(T_s(z)-z)+s\bigr)=1. \tag{4}$$

For the specific phase in (2), the complete cubic reader gives:

* each $t_{s_i}$ lies outside $K$ and has unique $K$-neighbour $-s_i$;
* each of the fifteen differences $t_{s_j}-t_{s_i}$, $i<j$, lies
  outside $K$ and has **no** $K$-unit neighbour.

The JSON records all twenty-one decisions. The first six also have
the independent identity (4); the fifteen negative decisions use the
complete algorithm (3), not a bounded search over elements of $K$.

This finite exact coefficient calculation proves a statement for
arbitrarily enlarged base sets. For **any** $A'\subset K$ and any
$R_i\subset A'$, form the complete unit graph on

$$A'\cup\bigcup_i(R_i+t_{s_i}). \tag{5}$$

Its branches are disjoint from $A'$ and from one another. A point
$r+t_{s_i}$ has an old neighbour precisely when $r-s_i\in A'$;
that is its only old neighbour. Distinct new branches have no edges.
Within a branch, translation preserves exactly the unit graph of $R_i$.

Consequently

$$\pi(a)=a\ (a\in A'),\qquad
\pi(r+t_{s_i})=r \tag{6}$$

is a graph retraction onto the old graph. Internal edges retain their
unit differences. An old contact maps to $r\leftrightarrow r-s_i$,
also a unit edge. There are no other cases. Thus restriction and pullback
preserve the **entire** colouring support on every old subset, for any
number of colours. Enlarging the $R_i$ inside this same old field cannot
turn these six branches into a forcing gadget.

For the actual candidate, (6) is stored for every vertex and checked on
all 6,118 complete edges, including all aliases and all old vertices.
This source-fixed certificate is stronger than testing a particular
collection of successful colourings.

## 5. Four bounded geometric recipes

The search tested exactly these recipes before stopping. All rows contain
the old graph $A$ literally; no old edge or point is removed.

| Added points | Total points / edges | New points | New points in eight-port core | Old/new edges | Old/new edges in that core |
|---|---:|---:|---:|---:|---:|
| $C+u(M+D)+S$ | 1482 / 6034 | 768 | 648 | 192 | 192 |
| $R+\{2(uq-q),\,2(T_{s_i}(uq)-q)\}$ | 1610 / 6452 | 896 | 549 | 0 | 0 |
| $R+\{0,t_{s_0},\ldots,t_{s_5}\}$ | 1482 / 6118 | 768 | 470 | 480 | 334 |
| union of the first and third rows | 2250 / 9238 | 1536 | 1118 | 672 | 526 |

Every new point in every row is outside $K$. The final row has no
edge between the two exterior families and no alias joining them; it merely
combines their separate obligations. It was not sent to a larger SAT
search.

The rotated-seed baseline has its own exact explanation. Since
$[K(u):K]=3$, the elements $1,u,u^{-1}$ are linearly independent over
$K$. A unit difference $b+ua$, with $a,b\in K$ both nonzero, would
contradict that independence after expanding its norm. Thus the new
block has only Cartesian factor edges. The map
$c+ua+s\mapsto c+a+s\in A$, fixing $A$, preserves every actual edge.
Old-to-new contacts can only be the unit rotated-seed steps from their
old centre. This proves the baseline retraction independently of its
finite edge census.

The second row also has a checked retraction to $R\subset A$, but its
complete finite edge list has no old-to-new edge at all. The union of
the first and third rows has the combined checked retraction. All four
maps are retained and replayed during verification; only the main
candidate's complete coordinate/edge/map payload is exported separately.

For the main candidate, fixing the eight original ports leaves 942 core
vertices, of which 470 are new. Fixing **all 714 old vertices** instead
leaves 1,205 vertices and 4,944 edges, of which 491 vertices are new.
These are different relative-deletion problems. Neither core size
overrides the retraction certificate.

## 6. Colour consumers, hostile, and next obligation

`audit_retraction` checks exact vertex types, the literal old induced
subgraph, fixation of every old vertex, and every mapped edge without
collapse. `extend_colouring` then consumes a supplied proper old
six-colouring and pulls it back. It verifies the full output, including
all 714 prescribed old colours. This graph consumer does not infer a
metric embedding from an abstract edge list; the geometry is provided
by the separately authenticated recipe.

The checker applies it to all 570 saved full old words. The sibling
[variable pair fibres](variable_pair_fibres_20261018d.md) separately
solves all 570 extension problems with every old colour pinned. This is
independent evidence for the extension claim; the universal theorem is
the retraction, not extrapolation from that bank. The terminal pair
$(2,4)$ retains squared distance $(3+\sqrt5)/2$, so remains nonunit.

A cheap positive hostile marks the next gate precisely. The four points
$0,1,u,1+u$ form a unit rhombus, with graph $C_4$. Forgetting the $u$-layer collapses its
two matching edges, so that proposed retraction fails. Yet the proper
colouring $(0,1,1,2)$ exists. More generally, $A'\cup(A'+u)$ is the
Cartesian two-layer graph and every old six-colouring extends by adding
one modulo six on the new layer. Breaking one projection is therefore
not equality forcing.

For a future pair of fractional branches, the complete reader can first
ask whether their difference has a $K$-neighbour $a$. If none, those
branches cannot couple at any old-field scale. If $N(a)=1$, an actual
edge still maps to an old unit edge under the forgetful map. A candidate
which defeats this map must realize a contact with $N(a)\ne1$, including
$a=0$ for a projected loop, and retain the appropriate old-field
difference in its finite sets. It must then defeat all surviving full
colourings while keeping an equal-colouring witness. The rhombus hostile
shows why even that contact gate is only necessary for this route.

## 7. Reproduction and audit

[Script](../../04-computation/experiments/surviving_repair_geometry_20261018d.py)
and [JSON](surviving_repair_geometry_20261018d.json) retain the field basis,
phase, all contact decisions, literal old inclusion, full main candidate,
peel orders, retraction, four finite-recipe summaries, and the supplied
colour-word replay hash. Finite-place rejection uses only primes at which
every coordinate denominator is a unit and retains the actual conjugate
embedding. Surviving norms are tested exactly. A separate one-place
rejection path independently recovers the full main edge list.

```
python -B -X utf8 04-computation/experiments/surviving_repair_geometry_20261018d.py
python -B -O -X utf8 04-computation/experiments/surviving_repair_geometry_20261018d.py
```

Normal and optimized runs agree with the saved JSON at **28,486 exact
checks**. Its LF-normalized SHA256 is
`6cd2fbbf6cbd03836cf4b43e947e2bb2aa521af8a30cfa544faf65aea2961177`.
The main coordinate/edge payload hash is
`512de33d3227a2a0cc29bd8e7918808cf620adc8b9193d26469476ebb8894ca8`.

**Independent proof/API and retraction-consumer audit PASS.** The reviewer
checked all 1,482 typed map entries, 714 fixed old vertices, 6,118 mapped
edges and 570 complete pullbacks against every edge and pin, with
5,157,698 independent gates. The complete cubic-reader proof, including
the singular determinant case, and the arbitrary-base star retraction
also passed the proof skim. The draft's geometric word “square” was
corrected to “rhombus”; the union statement explicitly concerns edges
between the two exterior families, not their internal edges.

No conclusion here assumes a successful full-support solver or a new
unit-terminal equality fragment. Only these three owned packet files
were edited; publication and shared navigation are handled separately.
