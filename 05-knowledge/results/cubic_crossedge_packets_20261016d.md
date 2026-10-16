# Cubic unit circuits and compatible cross edges

2026-10-16. **PROVED:** the new direction has degree three over the old
cyclotomic field, supplies two exact 32-step unit walks with a common
endpoint, and has a trace projection with an explicit positive variance
defect. The enlarged degree-24 base field still admits five colours.
**FINITE-EXACT:** five complete induced geometric assemblies, each exactly
three-chromatic, with the original eight ports fixed.
**OPEN:** a six-colour forcing context and the plane's six/seven alternative.
Self-audited; no independent referee or historical priority claim is made.

## 1. Inheritance and the live interfaces

The closest proved construction is G7 of
[elliptic cross-edges](elliptic_crossedges_20261016c.md): an odd cubic
extension supplies a norm-one direction of exact quintic-basis width three.
G2 of that note proves that the previous two rational elliptic contact
families, with both signs and all rotation indices, form an infinite tree.
G4--G5 close the two-seam, shared-parameter linkage route. We do not
repeat larger finite versions of those closed tests.

The incoming
[seventeen-vertex extension theorem](hn_elliptic_bridge_20261016c.md),
including its odd-degree corollary, preserves every old six-colouring
when at most seventeen exterior vertices are adjoined. It already applies
to the new cubic carrier. A field-escaping ear or a new edge alone is
therefore an insufficient target.

The hostile is genuine algebraic field escape with an unchanged full
boundary relation. The corrected near miss is counting abstract contact
families without retaining which cross-edges coexist in one geometry.
The underused sidecars are the integer cubic generator, the exact
coordinate identifications of a unit-word circuit, and the prime at 19.

Anchor: compatible cubic cross-edges and native Collatz cancellation.
Niche: the local colouring of the enlarged base field. Wildcard: the
trace/variance comparison with refinement energy. The board is fixed
ports / actual common endpoint / cubic marking / old-field residue /
full induced graph / source-clock interval coverage. The
[Collatz companion](collatz_interval_matching_20261016d.md) gives the
source-preserving interval construction and its independent open target.

## 2. The two field towers and a retained five-colour ceiling

Keep the previous literal notation

\[
 K=\mathbb Q(\zeta_5,\sqrt{-11}),\quad
 L=\mathbb Q(\zeta_{55}),\quad \xi=\zeta_{11}.
\]

Let theta be the root in `(4/5,1)` of
`2theta^3+2theta^2-2theta-1=0`, and put

\[
 w=2\theta,\quad f(w)=w^3+2w^2-4w-4=0,
 \quad K'=K(\theta),\quad L'=L(\theta).
 \tag{1}
\]

The inherited degree statements are
`[K':Q]=24`, `[L':Q]=120`, `[L':K']=5`, and `[L':K]=15`.
The discriminant of f is `592=16*37`. Its three roots are real: there
is one in each of `(-4,-2),(-1,0),(1,2)`, by sign changes and degree.
The cubic rational-root test and degree coprimality establish the stated
extensions; no additional square root is being adjoined.

**C1 (PROVED).** The full unit-distance graph on K' satisfies

\[
                         3\le\chi(K')\le5.
 \tag{2}
\]

**Proof.** The recovered K colouring in
[rhombus overlap, G2](rhombus_overlap_elliptic_20261011.md)
uses a conjugation-fixed prime over 19, with residue field `F_(19^2)`.
Here `f(4)=76=0 mod19` and `f'(4)=60=3 mod19`. Hensel lifting therefore
gives a root in `Q_19`. Embedding w as that root extends the selected
valuation of K to K' without extending its residue field. Complex
conjugation fixes the chosen root, so the extended valuation remains
conjugation-fixed.

For every geometric unit u in K', `u*bar(u)=1` gives `2v(u)=0`.
Every unit direction is therefore a local unit. In each connected graph
component, subtract a chosen vertex; all remaining vertices are sums of
unit directions and lie in the valuation ring. Reduction maps every unit
edge to the norm-one graph on `F_(19^2)`.

Write that field as `F19[i]`, `i^2=-1`. The old verified five-colouring is
`(a,b) -> c(a+b mod19)`, where

```text
c = (4,3,2,4,3,2,1,0,4,1,0,4,3,2,1,0,2,1,0).
```

The new checker directly verifies all twenty norm-one directions and all
nineteen translates again. Pull back the colouring componentwise. The
retained unit pentagon proves the lower bound. QED.

Thus changing elliptic generators while staying in K' cannot produce a
six-chromatic graph. This does not assert that all prescribed six-colour
boundary assignments extend across every K' context. The colouring also
does not extend by this argument to L': its eleventh root requires a
larger residue field. The prime, its fixedness and the actual residue
field, rather than odd degree alone, carry this upper bound.

## 3. A genuine cubic circuit made of unit steps

The inherited direction is

\[
 v=\frac{\theta-\theta\xi+\xi^2}
          {\theta-\theta\xi^{-1}+\xi^{-2}},\qquad v\bar v=1.
 \tag{3}
\]

Its exact width-three expression over K' is already proved in the
preceding note. The new tensor implementation independently verifies
the division and norm in the 120-dimensional basis `zeta55^i w^j`,
`0<=i<40, 0<=j<3`.

Define

\[
 P=1-\xi-v+v\xi^{-1},\qquad Q=v\xi^{-2}-\xi^2.
 \tag{4}
\]

Equation (3) gives `theta P=Q`. Also P is nonzero: P=Q=0 would force
`v=xi^4` and then `P=(1-xi)(1+xi^3)`, impossible for an eleventh root.
Consequently `L(v)=L(theta)`, and v has degree **exactly three over L**.
This is a different field comparison from its width over K'.

**C2 (PROVED).** There is the exact unit-word relation

\[
       \xi^6\bigl(2Q^3+2Q^2P-2QP^2-P^3\bigr)=0.
 \tag{5}
\]

Every monomial after expansion is a unit direction `v^i xi^j`.
Reducing j modulo eleven leaves 28 nonzero integer coefficients, whose
positive and negative masses are each 32. Thus (5) gives two literal
32-step unit walks from zero to the same point. The full coefficient
list and ordered walks are retained in the receipts.

The mechanism is exact substitution of `theta=Q/P` into its cubic.
The displayed expression is cubic in v over L and is nonzero as a
polynomial; its leading coefficient has xi-degree six and cannot vanish
modulo `Phi11`. It is a multiple of the minimal polynomial of v.
No minimality of the 64-term unit relation or of either walk length is
claimed: the carrier already has cyclotomic relations.

The control theta=1 still supplies a unit, namely `v=xi^2`, but (5)
becomes `xi^6(1-xi^2)^3`, which is nonzero. Thus this closing seam really
uses the cubic parameter; it is not a closure that holds for every
Hilbert90 direction. A common endpoint is nevertheless only a graph
incidence. It has not by itself increased a chromatic requirement.

## 4. Complete assemblies and where the coupling actually occurs

Every test keeps the original pentagon and rhombus coordinates as its
first eight ordered vertices. Let alpha be zeta5, eta the original rhombus
unit, and gamma=(5+sqrt(-11))/6. All sums below are literal coordinate
sums, followed by deduplication and a test of **every** pairwise distance.

| Declared assembly | Vertices | Edges | Degeneracy | Chromatic number |
|---|---:|---:|---:|---:|
| Sums of at most two signed `v^k xi^h`, k=-1,0,1 | 2251 | 4364 | 2 | 3 |
| Two-layer xi/v bank plus alpha, eta and gamma directions | 2381 | 4684 | 5 | 3 |
| Two-layer bank using all five K-fixing conjugates of v | 9944 | 19634 | 5 | 3 |
| Eight port translates of the cubic walk bank | 1828 | 3949 | 4 | 3 |
| Five K-fixing Galois copies of that walk bank, plus ports | 1139 | 1200 | 2 | 3 |

The script fixes the exact direction sets and ordering. For the first
three rows, h runs from zero through ten and both signs are included;
the full lists distinguish the mixed and conjugate banks. Each lower
chromatic bound is the original pentagon. Each upper bound is an explicit
colour vector checked on every authenticated edge, not a timeout or a
claim based on degeneracy.

The cubic walk bank uses four stated orderings of the 32 positive steps
and the 32 negative steps. Its 229 distinct points retain all eight
length-32 walks and their common endpoint. Translating it by all eight
ports creates four coordinate identifications. **Every** actual edge
has an address pair sharing either a port-translate index or a walk-bank
vertex index. No edge changes both coordinates without such an alternative
address. Thus this particular translated bank still fails to supply an
additional mixed constraint. It is a finite address test, not a theorem
about arbitrary translations of arbitrary cubic circuits.

For the last row, the five automorphisms fix K and theta and send
`xi -> xi^h`, `h in {1,3,4,5,9}`. They preserve norm-one equations because
they commute with conjugation, but are not asserted to preserve every
real Euclidean distance. The completed graph has fifteen cross-copy
edges whose endpoints share no copy membership. **Both endpoints of every
one of those edges lie in L**: the cubic coordinates vanish. There are no
cross-copy edges incident to a point requiring the cubic extension in
this finite assembly. The three-colouring alone would have hidden that
more informative failure mechanism.

All five assemblies preserve the entire original six-colour port
relation, consisting of 570 proper equality patterns. This follows from
their three-colourable exteriors and
[THM-4614, three-colourable contexts](../../01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md),
with its port-census hypotheses transported as in G7 of the preceding
note: `[L':K]=15` is odd, so a common unit neighbour of two K ports still
lies in K. The same nine double-port points and absence of triple-port
points remain. This uses that precise theorem, not the generally false
inference that a three-colourable whole graph preserves every boundary
assignment.

The live geometric question is now narrower: build a compatible induced
assembly whose **essential cross-copy edges involve the cubic coordinates**,
and then test the full boundary relation. The current circuit closes a
new seam but the tested operations do not turn it into that coupling.

## 5. Trace fixes the ports but pays a variance defect

Let `Tr=Tr_(L'/L)` and `m(u)=Tr(u)/3`. Since
`Tr(w)=-2` and `Tr(w^2)=12`, the mean is particularly cheap in the retained
three-block basis. It fixes every old port and commutes with addition and
conjugation.

**C3 (PROVED).** For every u in L' of norm one,

\[
 1-|m(u)|^2
   =\frac13\sum_{r=1}^3|u_r-m(u)|^2\ge0,
 \tag{6}
\]

where u_r are the three embeddings over L. Equality holds exactly when
u belongs to L. Indeed all three cubic roots are real, so each embedding
commutes with complex conjugation and each u_r has norm one. Expand the
squares. Vanishing means that the three conjugates coincide; separability
then places u in the base field.

For the new direction v, the defect is strictly positive. The checker
stores its exact L-coefficients and verifies
`Tr((v-m(v))*conj(v-m(v)))=3*(1-m(v)*conj(m(v)))` in rational arithmetic.
Thus projection to the old field preserves every source port but strictly
shortens every unit edge whose difference leaves L. It cannot transfer
the completed unit graph as a unit graph without the variance data.

This is a precise comparison with the inherited Collatz Gram remainder:
a linear projection preserves a source interface while losing a positive
quadratic quantity. In geometry, three conjugate coordinates restore
the metric. In Collatz, physical detail rows restore the quadratic form,
while original source atoms and integer clocks restore the native paths.
Neither variance identity supplies the other's open forcing or
cancellation bound. The incoming marked Lyness work likewise retains
its differential; a modulus alone does not certify a polynomial Keller
map. No planar Jacobian conclusion is added here.

## 6. Reproduction and audit scope

[Program](../../04-computation/experiments/cubic_crossedge_packets_20261016d.py).
Exact receipts:
[powers](cubic_crossedge_packets_20261016d_powers.json),
[mixed](cubic_crossedge_packets_20261016d_mixed.json),
[conjugates](cubic_crossedge_packets_20261016d_conjugates.json),
[circuit](cubic_crossedge_packets_20261016d_circuit.json), and
[conjugate circuits](cubic_crossedge_packets_20261016d_circuit_conjugates.json).

```text
python3 -B 04-computation/experiments/cubic_crossedge_packets_20261016d.py --mode circuit
python3 -B -O 04-computation/experiments/cubic_crossedge_packets_20261016d.py --mode circuit
```

Replace the mode by `powers`, `mixed`, `conjugates`, or
`circuit_conjugates` for the other universes. Normal and optimized runs
agree byte for byte in all five modes. Coordinates use the integer
relations `Phi55(z)=0` and (1). Evaluation at primes 331, 661 and 991,
with authenticated roots of both polynomials and denominators prime to
the filters, can only reject nonedges. Every survivor is tested by its
full exact norm. All old port pairs also receive a direct unfiltered test.

Receipts retain the generator coordinates, fixed port coordinates, full
edge sets, checked colour witnesses, peeling orders, circuit words,
coordinate hashes and the exact reconstruction recipe. These are finite
universes, not censuses of all unit directions in L'. The two independent
algebra presentations, the theta=1 failed closure, the residue table and
the trace defect provide positive and hostile controls beyond the colour
solver. No unchecked solver nonexistence claim is used.
