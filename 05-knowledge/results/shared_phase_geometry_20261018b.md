# Shared phases: a uniform rotation obstruction and a fractional escape test

2026-10-18. **PROVED:** the guarded reduction from simultaneous unit-distance
events to one finite residue root, and the resulting all-phase exclusion for
the stated finite layer banks. **FINITE-EXACT:** their residue certificates
and the complete 56-point fractional configuration. **OPEN:** a configuration
which retains an equal six-colouring and excludes every separating one.

The largest tested bank has 128 layers and 1,024 points. At every unit phase,
at most ten new pairs of original anchor labels can be forced by its edges.
This is below the necessary thirteen-pair threshold. Increasing this bank
from sixteen to sixty-four powers therefore does not produce an equality
fragment. A separate radius-two fractional construction escapes the loopless
projection premise, but an explicit relative degeneracy certificate still
extends every colouring of the original source.

## 1. Inheritance and the literal source

We use the exact field, embedding, anchors, and Möbius interfaces of
[cubic boundary search](cubic_boundary_search_20261018a.md), without changing
that frozen packet. Write

$$F=\mathbb Q(\zeta_{15},j,w),\qquad j^2=-11,\quad
w^3+2w^2-4w-4=0,\quad 8/5<w<2.$$

Complex conjugation sends $\zeta_{15}$ to its inverse, $j$ to $-j$, and fixes
$w$. This is the inherited degree-48 field. In particular it **contains**
$\sqrt{-3}$. The older exterior lower bounds for odd extensions of fields
without $\sqrt{-3}$ do not apply here.

Let $C=(c_0,\ldots,c_7)$ be the eight original points, in their original
order. Their complete unit graph $B=C_5\vee C_4$ has edges

$$01,04,05,07,12,23,34,56,67.$$

It has nineteen nonedges and 570 proper equality patterns using at most six
colours. The exact coordinates, field basis, and complete edges are copied
into this packet's JSON; no replacement source is used.

The inherited scalars are denoted $b,c$, and

$$u=\frac{w+3-j}{w+3+j},\qquad u\bar u=1.$$

Their full definitions remain in the cited source and in the JSON. The
current concept board is **source labels / difference provenance / shared
phase / common residue root / projected loops / full boundary support**.
The relevant corrected near miss is restriction to a finite norm-one circle:
a split residue place need not carry complex conjugation to finite-field
Frobenius. The canonical hostile below uses an actual event root outside that
finite circle. The least-used retained sidecar is the *common* phase of all
edge events, rather than separate existence of their roots.

## 2. One shared phase, one root bucket

For distinct layer coefficients $\lambda_0,\ldots,\lambda_{s-1}\in F$, set

$$x_{i,a}(q)=c_a+(\lambda_i-\lambda_0)q,\qquad |q|=1. \tag{1}$$

Layer zero is literally $C$. For $i<j$, put $v=\lambda_j-\lambda_i$ and
$d=c_b-c_a$. Expanding the squared distance and using $\bar q=q^{-1}$ gives

$$|d+vq|^2=1
\quad\Longleftrightarrow\quad
P_{v,d}(q)=d\bar v+(d\bar d+v\bar v-1)q+v\bar d\,q^2=0. \tag{2}$$

For each original nonedge we retain both orientations. In this source all
38 oriented differences are distinct, each with its original ordered pair.
Thus a root event has a definite bit among the nineteen nonedge labels.
Repeated witnesses of the same pair are combined by Boolean OR, not counted
as new constraints.

Choose a prime $p$ and a ring homomorphism from a localization of the source
coordinate ring to $\mathbb F_p$. The script validates all three defining
polynomials, the primitive fifteenth root, and every rational denominator.
It evaluates $z$ and $\bar z$ **separately at this same place**. Write the
reduction of (2) as $A+Bq+Cq^2$.

**Common-root lemma.** Suppose every relevant constant coefficient $A$ and
leading coefficient $C$ is a unit at this place. Every simultaneous set of
characteristic-zero events at an algebraic phase $q$ is contained in a single
bucket indexed by a root $\beta\in\mathbb F_{p^2}$ of their reduced
quadratics.

To prove it, extend the place to a number field containing $q$. Dividing one
event equation by its unit leading coefficient makes it monic, so $q$ is
integral at that extension. The reciprocal equation, divided by the unit
constant coefficient, similarly makes $q^{-1}$ integral. Reduction is
therefore defined and nonzero. All event equations reduce at the same place
and the same $\beta$. Since $\beta$ satisfies a nonzero quadratic over
$\mathbb F_p$, it lies in its quadratic extension, even if the ambient
residue field is larger. If $q$ is transcendental over $F$, none of the
nonzero quadratics can vanish, so there are no new events.

Enumerating **all** quadratic roots and OR-ing their pair masks can add false
positive constraints after reduction. It cannot delete an actual simultaneous
constraint. In particular it gives an upper bound on the number of actual
new projected pairs. Neither a residue hit nor a maximum bucket is asserted
to lift to a geometric configuration.

Two further guards make the projection $x_{i,a}\mapsto a$ valid uniformly:

* A coincidence of two different layer points would require
  $v\bar v=d\bar d$. We exclude these equalities at the residue place for
  every anchor difference, and exclude $v=0$ as well.
* A unit edge whose two endpoints project to the same anchor would require
  $v\bar v=1$. We exclude this at the place too.

Failure of any guard makes the program refuse the all-phase certificate. It
does not silently omit the offending layer or edge. Within a layer the exact
source graph is unchanged. With these guards, all source and copied points
remain distinct for every unit phase and the projection has no loops.

### Why twelve projected pairs are still insufficient

[The exceptional-factor support theorem](exceptional_factor_support_20261018a.md)
proves that adding at most twelve pairs to this eight-vertex source still
permits a proper six-colouring separating each original nonedge. Equivalently,
adding one specified original nonedge as an edge cannot yet create a
seven-chromatic obstruction. The sharp necessary threshold is thirteen
distinct pairs; it is not a sufficient test for a geometric equality fragment.

Apply a separating colouring of the actual quotient to every copied point.
It gives a separating colouring of the complete unit graph in (1). Thus the
common-root bound of at most twelve proves exclusion for every phase. This
conclusion concerns each target pair separately; the bound alone does not
claim that all 570 source patterns survive in the rotation banks.

## 3. Exact finite certificates with an all-phase consequence

The banks use $\lambda_{k,i}=ku^i$, with $k=b$, $k=c$, or both, and
$0\le i<r$. The subtraction of $\lambda_0q$ in (1) is essential for literal
source retention. No value of $q$ is sampled in this certificate.

The three places are specified by $(p,\zeta_{15},j,w)$:

$$ (10321,3678,1363,394),\quad(10531,846,3442,9367),\quad
(100621,2204,9259,47934). $$

All leading, constant, collision, and loop guards pass in the following
complete finite universes. “Maximum” means the largest residue pair-mask
size, hence an upper bound on actual simultaneous new pairs.

| Prime | Radii | Powers per radius | Layers / points | Residue roots | Maximum |
|---:|:---|---:|---:|---:|---:|
| 10321 | $b$ | 16 | 16 / 128 | 8302 | 4 |
| 10321 | $c$ | 16 | 16 / 128 | 7588 | 4 |
| 10321 | $b,c$ | 16 | 32 / 256 | 26560 | 7 |
| 10531 | $b$ | 16 | 16 / 128 | 8178 | 3 |
| 10531 | $c$ | 16 | 16 / 128 | 8274 | 4 |
| 10531 | $b,c$ | 16 | 32 / 256 | 27864 | 7 |
| 100621 | $b,c$ | 32 | 64 / 512 | 128280 | 6 |
| 100621 | $b,c$ | 64 | 128 / 1024 | 392584 | 10 |

The last row checks all $38\binom{128}{2}=308864$ event polynomials. Its
root-bucket histogram, indexed by distinct projected pair count, is

$$1:311898, 2:22432, 3:24304, 4:18064, 5:9520,\quad
6:4348, 7:1396, 8:504, 9:112, 10:6.$$

The JSON stores all source differences and their provenance, layer labels,
maximum buckets, guards, and a hash of the complete sorted bucket table.
It retains a concrete actual-source event root outside the finite norm-one
circle. The polynomial controls independently enumerate every quadratic root
over $\mathbb F_{p^2}$ for $p=3,5,7$. A planted family
$(q-1)(q-a_j)$ has thirteen distinct pair labels at $q=1$, ensuring that the
implementation can retain the threshold; this is an algebraic control, not
a geometric construction. The leading-unit hostile $1-q+3q^2$ demonstrates
why a root escaping to infinity at a bad place cannot be ignored.

## 4. A genuinely fractional radius-two escape

We now deliberately change the radius to $k=2$ and use fractional maps. This
is outside the preceding $b,c$ rotation universe. Put

$$\delta=\zeta_3,\quad \eta=1+\delta,\quad
t=\frac{-1+\sqrt{-15}}4,\quad q=\delta t,\qquad
T(z)=\frac{z-1/2}{1-z/2}. \tag{3}$$

All these elements belong to the inherited field $F$: for example
$\sqrt{-15}=\sqrt5(1+2\delta)$. Since $2t^2+t+2=0$ and
$\delta^2+\delta+1=0$, we get

$$|q|=1,\quad 2q^2+\delta q-2(1+\delta)=0,
\qquad 2(T(q)-q)=\delta. \tag{4}$$

Here $T(q)-q=(q^2-1)/(2-q)$ and $q\ne2$, so the denominator is legitimate.
Equation (4) gives eight actual same-anchor unit edges between the original
layer and its $T$-layer. Thus the loopless projection premise has really
failed; these are not residue-only events.

For $s\in\{1,\eta\}$ write
$T_s(z)=(z-s/2)/(1-\bar s z/2)$, and let $R_u(z)=uz$.
Use the eight maps

$$\mathrm{id},T_1,T_\eta,T_1\circ T_\eta,T_\eta\circ T_1,
R_u,T_1\circ R_u,T_\eta\circ R_u$$

at the single exact phase (3), and the literal coordinates
$c_a+2(g_i(q)-q)$. All maps and poles are checked. There are **64 marked
addresses, 56 distinct points, and 74 complete unit edges**. The map
$T_\eta\circ T_1$ fixes this phase although it is not the identity rational
map: its eight addresses alias the original eight points. The JSON retains
the entire raw-to-point map, every coordinate, every edge, and every map.
The first two layers alone have sixteen points and 27 edges: their two
nine-edge copies, eight matching edges, and one additional cross edge.
Their incidence also passes an independent unsieved exact-norm comparison.

**Positive extension certificate.** Holding the eight original points fixed,
the program removes each of the other 48 vertices while its degree in the
remaining graph is at most five. Reinsert them in reverse order; at most
five colours are forbidden at each insertion. Therefore every labelled
proper six-colouring of the original source extends. The receipt includes
the removal order and explicit lifts of all 570 canonical source patterns;
palette permutations give all their labelled versions.

This is a concrete fractional escape from the projection obstruction, but
it does not remove any boundary support. It also explains why the aliases
cannot be thrown away as “duplicate map words”: they are equality of actual
points at this phase and must share a colour.

## 5. Remaining target and exact interface

The companion
[shared-phase threshold hostile](shared_phase_threshold_hostile_20261018b.md)
supplies a different actual metric control: thirteen new projected pairs,
112 points, and no source-fixing retraction, yet all 570 source patterns
still extend by palette changes. Together the controls separate three
conditions: enough projected pairs; failure of a retraction; elimination of
all separating colourings. Only the last, together with an equal-colouring
witness, certifies the desired equality fragment.

The next productive target is a coupled fractional family whose exterior
survives the relative degree-five deletion test, followed by its complete
boundary CSP. For general Möbius layers the cleared event equation can have
degree four. A common-residue sieve must then include irreducible factors of
degrees three and four, or equivalent projective residue roots, as well as
poles, leading degeneration, coincidences, and projected loops. The
quadratic-extension sieve proved here does **not** certify that larger
universe. No uniform fractional-family exclusion or successful equality
fragment is claimed.

The map (1) preserves source labels and exact unit events; the residue map
preserves necessary event equations but loses liftability. The pair-mask
quotient loses full colouring compatibility. The retained sidecars are the
literal source, conjugate-place evaluation, common root, event provenance,
and full graph/colour lifts for the fractional packet. In particular a
residue maximum, a contact count, and a single separating witness are not
interchangeable certificates.

## 6. Reproduction and audit

[Script](../../04-computation/experiments/shared_phase_geometry_20261018b.py)
and [JSON receipt](shared_phase_geometry_20261018b.json) use exact integers,
rational field coordinates, and always-active validation. The public
compiler retains a named source and radius bank, rejects Boolean/float
parameter aliases, authenticates the chosen place, and fails closed on
missing unit guards.

```
python -B -X utf8 04-computation/experiments/shared_phase_geometry_20261018b.py
python -B -O -X utf8 04-computation/experiments/shared_phase_geometry_20261018b.py
```

Normal and independent optimized runs match at **43,238 exact checks**. The LF-normalized JSON
SHA256 is `134f845b169e829d04a42e79b065edca2bdcb1bbdcebed3eb6b1ae7efe0393af`.
**Independent audit PASS:** the complete note and production APIs, common-root
integrality, all quadratic-extension roots, conjugate-place/collision/loop
guards, fractional identity, actual aliases, and all 570 source-fixed lifts.
The independent optimized run took 96 seconds and exactly reproduced the
saved JSON. All local links resolve.

The prepublication review corrected one wording error: the separating test
adds the marked edge; it does not identify its endpoints. The implemented
pair-mask criterion and all computations already used the separating test.
