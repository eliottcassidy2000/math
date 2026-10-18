# Cubic phase layers: a generic retraction and exact exceptional candidates

2026-10-18. **PROVED** rational-identity classification and finite exceptional-phase reduction. **FINITE-EXACT** complete induced graphs and explicit colouring lifts for the ten configurations below. **OPEN**: a six-colourable unit-distance fragment with equal endpoints and no separating six-colouring. No such fragment is produced here.

The positive search reduction is this: for a fixed finite base and a fixed finite list of rational circle maps, every edge that can defeat the natural layer projection occurs at a root of an explicit polynomial of degree at most four over the marked coefficient field. Generic phases give no new boundary restriction. Exact low-degree phases do create new edges in the experiments, but another retraction preserves all the old boundary support.

## 1. Inheritance and the field boundary

The current sources are:

- [THM4661, odd-degree carriers and the six-colour exterior boundary](../../01-canon/theorems/THM-4661-odd-degree-carriers-a-failed-six-colour-exterior-has-at-least-twenty-vertices-and-equality-testers-are-seven-colour-certificates.md). The proved failed-exterior floor is **20**, superseding the earlier 17-point result. The further floor 22 is explicitly **FINITE-EXACT**, not a proof-checked theorem.
- [THM4700, the unit-path ring and coefficient-triangle census](../../01-canon/theorems/THM-4700-the-unit-path-ring-of-the-six-colour-field-short-reciprocal-walks-and-the-heron-census-of-coefficient-triangles.md). In particular, neither of its coefficients $b,c$ admits a unit translation chord. Its reflected links and mixed quadrilaterals remain available.
- [THM4701, finite unit-circle graphs and Petersen interfaces](../../01-canon/theorems/THM-4701-kuratowski-graphs-the-petersen-graph-and-the-finite-unit-circle-graphs.md): a finite residue graph, its planar realization, and a colouring lift are separate objects.
- [The golden contact circuit](golden_contact_circuit_20261017i.md), whose selected incidence systems preserve all 570 original eight-port patterns; and the source-bound exact field/complete-incidence implementation [golden forcing search](../../04-computation/experiments/golden_forcing_search_20261017.py).

Use the marked field

$$F=\mathbb Q(\zeta_{15},j,w),\qquad j^2=-11,\qquad
w^3+2w^2-4w-4=0,$$

with the real root $8/5<w<2$, $\bar j=-j$, and $\bar\zeta_{15}=\zeta_{15}^{-1}$. Its inherited degree is 48. In particular **$F$ contains $\sqrt{-3}$**. The triangle-free odd-carrier hypotheses in THM4661 do not apply to this field, so neither its 20-point nor its 22-point floor is imported into these examples.

Write $N(z)=z\bar z$, $d_0=1+\zeta_5$, and

$$u=\frac{w+3-j}{w+3+j},\qquad N(u)=1.$$

The exact definitions of $b,c$ are inherited from the circuit: with $a=\bar d_0$,

$$m_2=-\frac{20+11j}{19},\quad m_1=\frac{20-11j}{19},\quad
c=\frac{am_2+1-a^2m_1+a^3}{N(d_0)-1},\quad b=cd_0+am_1-a^2.$$

The program verifies $k\ne0$ and $N(k)\ne1$ for both $k=b,c$. THM4700 supplies the stronger property

$$\nexists\omega\in F:\quad N(\omega)=1,\quad N(k(\omega-1))=1.\tag{1}$$

Indeed a pair of distinct unit points $q,q'$ with $N(k(q'-q))=1$ would give $\omega=q'/q$. This use of the inherited no-chord theorem is explicit; a finite failed search is not substituted for it.

Concept board: literal source coordinates; reflected circle links; exceptional low-degree phases; complete port support; alternative retractions. The hostile is a new geometric edge that still removes no boundary colouring.

## 2. Exact rational identity classification

Let $F\subset\mathbb C$ be conjugation-stable, let $k\in F^*$ with $N(k)\ne1$, and let $h\in\operatorname{PGL}_2(F)$ be a nonconstant rational map preserving the unit circle. Conjugation of a rational function on that circle means $\bar z=1/z$.

For $d\in F$,

$$N\bigl(d+k(h(z)-z)\bigr)\equiv1\tag{2}$$

holds exactly in the following cases:

1. $h(z)=z$ and $N(d)=1$;
2. $d=0$, $h(z)=\omega z$, $N(\omega)=1$, and $N(k(\omega-1))=1$;
3. $N(d)=1$ and

$$h(z)=T_{d/k}(z):=\frac{z-d/k}{1-(\bar d/\bar k)z}.\tag{3}$$

These are rational-function identities, not classifications inferred from sampled phases.

**Proof.** If $h$ is affine, preservation of the circle forces $h(z)=\omega z$. The norm of $d+k(\omega-1)z$ is identically one only if its constant or linear coefficient vanishes, giving cases 1 and 2. If $h$ has its pole at zero, circle preservation forces $h(z)=\omega/z$; the norm in (2) has a nonzero second harmonic and cannot be constant.

Otherwise write

$$h(z)=\frac{az+b}{cz+1},\qquad c\ne0.$$

Circle preservation and the nonzero determinant imply $N(a)=1$, $c=a\bar b$, and $N(b)\ne1$. Set $W=d+k(h-z)$. It has genuine poles at infinity and $-1/c$. The identity $W(z)\overline{W(1/\bar z)}=1$ therefore forces zeros at zero and $-\bar c$. Consequently

$$W(z)=\lambda\frac{z(z+\bar c)}{1+cz},\qquad N(\lambda)=1.$$

Comparing the leading coefficient gives $\lambda=-kc$, hence $N(k)N(c)=1$. Its value at zero gives $b=-d/k$. As $N(c)=N(b)$, this yields $N(d)=1$. Comparing the linear coefficient gives $(N(k)-1)(a-1)=0$. Thus $a=1$, which is (3). Conversely direct substitution proves each listed identity. All cancellations used in this argument are excluded by the nonzero Möbius determinant and $k\ne0$. $\square$

The reflection formula (3) is the same usable link retained by THM4700. The new point is the complete identity classification, including the rotation exception that must be excluded before asserting a projection theorem.

## 3. Generic phases preserve the entire source relation

Fix distinct $c\in C\subset F$, pairwise distinct nonconstant unit-circle Möbius maps $g_0,\ldots,g_{r-1}$ over $F$, and a scale $k$ satisfying (1). For a unit phase $q$ set

$$x_{i,c}=c+k\bigl(g_i(q)-g_0(q)\bigr).\tag{4}$$

This normalization retains $x_{0,c}=c$ literally. Suppose $q$ is algebraic of degree greater than four over $F$. Then all coordinates in (4) are distinct, and their **complete induced unit graph** retracts onto the original induced graph on $C$ by

$$\pi(x_{i,c})=c.\tag{5}$$

In particular every labelled $s$-colouring of $C$ extends by composition with $\pi$, for every palette size $s$. Restriction gives the converse. Thus the complete colouring relation on any marked subset of $C$ is exactly unchanged, including equal/separating endpoint status. This is not only an equality of counts.

**Proof and exact edge criterion.** For a pair of maps, write

$$g_j-g_i=R(z)/D(z),\qquad \deg R,\deg D\le2.$$

For $d=c'-c$, put $A=dD+kR$. If $P^{\#}(z)=z^2\overline{P(1/\bar z)}$ uses fixed grade two, the unit event has numerator

$$P_{i,j,d}(z)=A(z)A^{\#}(z)-D(z)D^{\#}(z),\qquad\deg P_{i,j,d}\le4.\tag{6}$$

At a nonzero unit phase away from poles, $P_{i,j,d}(q)=0$ is **equivalent** to the actual unit-distance condition. A degree-greater-than-four phase can satisfy it only if (6) is the zero polynomial. Apply §2 to $h=g_j\circ g_i^{-1}$. Property (1) excludes case 2. Hence every edge has $N(d)=1$; more precisely it is either an old same-map edge or $g_j=T_{d/k}\circ g_i$.

A coincidence has the degree-at-most-two numerator $A$. If it vanishes at $q$, it is a rational identity. Two nonconstant circle maps differing by a constant must be equal: their images are the entire unit circle, whose nonzero translate is a different circle. Thus only identical map/base labels could coincide, and these were excluded. The original layer is induced because its coordinates are unchanged. This proves (5). $\square$

The program checks six maps for each coefficient:

$$\mathrm{id},\ R_u,\ T_{1/k},\ T_{\eta/k},\
T_{1/k}\circ R_u,\ T_{\eta/k}\circ R_u,
\qquad R_u(z)=uz,\quad\eta=1+\zeta_{15}^{5}.$$

They give 48 distinct formal coordinates over the original eight ports, 40 beyond the source, with all 570 old six-colour patterns preserved by the theorem. An actual phase is $q=\zeta_{11}$: $\mathbb Q(\zeta_{11})$ has degree five over $\mathbb Q(j)$, while $[F:\mathbb Q(j)]=24$, so their intersection over $\mathbb Q(j)$ is trivial and $[F(\zeta_{11}):F]=5$. The program also verifies the two conjugate quintic factors multiply to $\Phi_{11}$; it does not approximate these coordinates numerically.

The no-chord hypothesis is material. In the same triangle-containing field, $k=1$ and $\omega=\eta$ give $N(\omega-1)=1$. Two layers of the same anchor then have a unit edge which projects to a loop, identically in the phase. This is not an exceptional-phase phenomenon.

## 4. A finite positive search reduction

For any fixed finite $C$ and map list satisfying the hypotheses, every pair that projects to a nonedge has a **nonzero** polynomial (6). All possible phase specializations at which the natural projection first fails therefore lie among the roots of a finite list of polynomials of degree at most four over $F$.

This is an effective finite algebraic candidate list. Its consumers must discard zero phases, poles, phases with $N(q)\ne1$, and coincident vertices. They must then reconstruct **all** unit edges, including edges not used to propose that phase, before testing the full boundary relation. Roots introduced by denominator clearing are not geometric certificates. A minimal polynomial factor is retained with its source coefficient field and labels.

The parent [exceptional-factor support sieve](exceptional_factor_support_20261018a.md) provides a further necessary test: collect all projected conflicts belonging to the same algebraic phase factor and test whether their quotient graph can still separate the marked endpoints. It is a preliminary obstruction, not a substitute for the full geometry; equal-endpoint survival must also be checked. This packet supplies exact field coefficients and event polynomials for that consumer.

## 5. Ten complete exceptional configurations

The original eight coordinates are retained literally from the cited source. Their induced edges are

$$01,04,05,07,12,23,34,56,67,$$

a pentagon and a four-cycle meeting at vertex 0. These have 570 canonical labelled-colour equality patterns with at most six colours.

For each $k\in\{b,c\}$ and $r=2,4$, use $g_i(z)=u^iz$ and add the anchors

$$d_j=u^j-k(u^j-1),\quad 1\le j<r,\qquad
C=(P_0,\ldots,P_7,d_1,\ldots,d_{r-1}).\tag{7}$$

At phase $q=1$, (4) has $r(7+r)$ distinct vertices. The exact complete incidence census gives:

| phases | source anchors | total vertices | exterior to full source | beyond eight marked ports | complete unit edges |
|---:|---:|---:|---:|---:|---:|
| 2 | 9 | 18 | 9 | 10 | 19 |
| 4 | 11 | 44 | 33 | 36 | 39 |

The complete graph is precisely $r$ disjoint copies of the old nine-edge port graph, with all the extra copied anchors isolated except for the $r-1$ spokes

$$0\ \longleftrightarrow\ (j,d_j),\quad 1\le j<r.$$

Here $(j,d_j)$ is the actual point $u^j$, of unit norm, whereas $N(d_j)\ne1$. Thus these spokes genuinely defeat the canonical projection (5). This is a positive geometric escape from the generic obstruction, not a claim that new boundary support has been eliminated.

In fact a different retraction survives. Fix every original source anchor; send each copied port to its old port and each copied additional anchor to old port 1. Every spoke maps to the existing edge $01$. This is a retraction onto the full nine- or eleven-anchor source, including its isolated anchors. It proves that **every** source colouring extends, not merely the 570 selected port patterns. The JSON retains explicit six-colour lifts for all 570 patterns as finite independent controls. Each complete graph has chromatic number exactly three: it contains the original $C_5$, and an explicit three-colouring is supplied.

For the $j$th designed event put $v=k(u^j-1)$ and $d=d_j$. After removing the extraneous zero-phase factor from fixed-grade (6), its exact quadratic is

$$A+(N(d)+N(v)-1)q+Cq^2,
\quad A=d\bar v,\quad C=v\bar d.
\tag{8}$$

Since $d+v=u^j$ is unit, (8) factors as $(q-1)(Cq-A)$. The second root

$$q_*=A/C$$

is unit and distinct from 1 in all six four-phase tests. For each scale and each $j=1,2,3$, the full 44-point graph at $q_*$ has exactly 37 edges: four old nine-edge copies and only the one selected spoke. There are no coincidences. The same repaired retraction applies, so all 570 patterns survive here too. These are the six additional complete packets, not an enumeration of all possible exceptional factors for the map family.

The next productive target is therefore simultaneous exceptional factors producing genuinely coupled constraints, followed by a complete six-colour boundary audit. More cubic phases or more isolated spokes cannot by themselves defeat the exhibited retraction. Neither a high vertex count nor an edge missing from the canonical projection proves a six-colour equality fragment.

## 6. Exact receipt and validity boundaries

The [script](../../04-computation/experiments/cubic_boundary_search_20261018a.py) exposes `Context.event_polynomial(d,k,left,right)`: maps are ordered $(a,b,c,d)$ tuples, acting by $(az+b)/(cz+d)$; polynomials are ascending tuples of exact 48-coordinate field elements. `Context.specialize` retains the literal source, verifies the unit phase and map denominators, rejects coincidences, and computes the complete induced graph. JSON coefficients use sparse triples `(basis_index,numerator,denominator)` in the explicitly recorded ordered tensor basis.

Every absent edge in the inherited complete census has either an exact norm rejection or a good-prime norm rejection; surviving modular candidates are checked in characteristic zero. The $q=1$ packets additionally have an unsieved, independent all-pairs exact norm census. The second-root packets use the complete good-prime census and recheck every retained edge exactly. No sampled subgraph is substituted for a full graph.

`retraction` checks the exact induced source edges, typed labels, source-fixing property and every actual edge. The resulting lift is an authenticated sufficient certificate for preservation of all source colourings. It is not a theorem that all support-preserving extensions admit a retraction.

During the unpublished draft audit, a cache could accept float or Boolean aliases after a valid exact field element had warmed the cache. Public field validation now occurs before private memoization. Four warm-cache regressions retain this correction mechanism. Invalid field dimensions, non-circle/constant maps, repeated layers and Boolean vertex labels are also rejected.

Reproduce normal and optimized output with:

```
python -B -X utf8 04-computation/experiments/cubic_boundary_search_20261018a.py
python -B -O -X utf8 04-computation/experiments/cubic_boundary_search_20261018a.py
```

The [JSON](cubic_boundary_search_20261018a.json) records all ten complete coordinate/edge packets, full labels, repaired retractions and explicit lifts. The mathematical scope is an exact obstruction and finite candidate reduction; global six/seven colourability and a six-equality fragment remain open here.

**Final audit.** Independent full proof/API review and optimized replay passed after the cache repair: **71,957 checks**, normal/optimized/saved JSON equal after LF normalization, SHA256 `a900e9670c26f29df7e5dbba84b3305920f916927fdfccd46f6d74c97759a649`. All local source links resolve. No inherited paper theorem or global colouring conjecture was re-audited or promoted by this verification.
