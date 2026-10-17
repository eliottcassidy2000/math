# Graded anchor deficits and valuation-aware cubic hub coupling

2026-10-17, clock_duality_oct15. **PROVED, scoped:** exact hole-moment
bounds, graded source-preserving controller gluing, an abstract obstruction
to one global threshold, complete elimination of the middle cubic
coordinate, and the four-branch prime-19 gate. **FINITE-EXACT:** the declared
native-controller and induced-graph experiments below. **OPEN:** uniform
Collatz anchor bounds, independently positive ROOT scores, universal entry,
and six-colour forcing. Self-audited; no independent referee or priority
claim. No numerical test is promoted to a universal statement.

## 1. Inheritance, incoming work, and the concept board

The anchor is uniform Collatz residual coverage; the niche is coupled
edges from the even-extension hub; the wildcard is graded arithmetic
compatibility. The board is **original atom capacity / integer clock /
physical anchor holes / coefficient field / individual prime branch /
fixed colour boundary**. The closest proved mechanisms are A1--A3 of
[the interval atlas](collatz_interval_atlas_20261016e.md) and B5--B8 of
[cubic adelic coverage](cubic_adelic_20261016e.md). Canonical hostiles are
the halving controller's cubic energy and the hub product whose
coefficient field omits j. The least-used sidecars are boundary-tail
positions and the separate primes over 19.

The inheritance pass incorporated origin/main through `fe0cea0d3cc4`.
The [independent S38 audit](../../04-computation/experiments/s38_audit_six_colour_20261016/REPORT.md)
replaces an unsaturated-lattice inference by an explicit twenty-step
relation, confirms the descent with the corrected square oracle, and
excludes eighteen- and nineteen-vertex five-cores. The twenty-two-vertex
failed-exterior floor remains **FINITE-EXACT**, with its stated carrier
and port hypotheses. We do not apply its odd-degree attachment lemma to
our even-extension hub.

The incoming [dense Collatz census](collatz_explicit_formula_dense_census_20261016.md)
resolves a sampling alias with irrational-density samples and records the
failure of `e_X <= 1` at `2^41`. Its exact source parity energies and our
independent finite-controller upper bounds are different objects. The
census does not supply the missing anchor estimate. This session tests
the latter directly, without a ROOT oracle.

One incoming claim requires correction: the audit's assertion that -11
is a square at every finite completion of F4 is false at 19. Section 7
repairs the justification while retaining the existing everywhere-local
solubility theorem. This is a new correction to the audit, not a reopening
of the already repaired descent.

## 2. Exact hole moments replace a maximum-gap requirement

Keep the interval atlas's notation: `T` is the shortcut map stopped at 1,
`t(n)` is the ROOT-supported full-clock sign (zero off ROOT),
`w(n)=3[n=2 mod 3]-1`, `H(M)=sum_{n<=M} w(n)t(n)` and
`a_N=H(4N)-H(N)`. A point certificate at a has residual B_a and pairs
oppositely signed original unit atoms at the same literal odd terminal.
Every path carries its full native integer clock. A1 transports this
certificate to row N at additional residual cost at most `4|N-a|`.
Controllers at different anchors may use different depths and terminals.

Let I be the integer interval `[L,R)`, with nonempty selected anchors
`a_1<...<a_m` inside I and `B_a <= T_I`. Put

\[
 \ell=a_1-L,\qquad r=R-1-a_m,\qquad g_i=a_{i+1}-a_i.
\]

**D1 (PROVED).** For `delta(N)=min_a |N-a|`, its exact first and second
moments are

\[
 R_1=\frac{\ell(\ell+1)+r(r+1)}2+\sum_i\left\lfloor\frac{g_i^2}4\right\rfloor,
\]
\[
 R_2=\frac{\ell(\ell+1)(2\ell+1)+r(r+1)(2r+1)}6
 +\sum_i\begin{cases}(g_i^3+2g_i)/12,&g_i\text{ even},\\
                         (g_i^3-g_i)/12,&g_i\text{ odd}.
         \end{cases}
\]

Consequently

\[
 \sum_{N\in I}a_N^2\le |I|T_I^2+8T_I R_1+16R_2.       \tag{D1}
\]

**Proof.** On a tail the distances are consecutive integers. In an
interior gap of length g the distances at its g-1 interior rows are
`min(t,g-t)`. Summing their first two powers gives the formulas; anchors
contribute zero. A1 supplies an actual transported certificate at each
row, bounded by `T_I+4 delta(N)`. Expand its square and sum. QED.

For a single octave with `T_I <= beta sqrt(X)`, the condition `R_2=O(X^2)`
suffices for energy `O(X^2)`. The exact arithmetic target is therefore a
**cubic moment of holes, including both tails**. It does not require
every gap to be `O(sqrt(X))`: an isolated gap of order `X^(2/3)` is allowed
by this bound. Merely counting successful anchors, or reporting their
interior maximum gap, loses the tails. The checker gives two 32-anchor
sets with identical interior gaps 16 but unequal coverage costs.

## 3. Grading is a genuine relaxation of a global threshold

**D2 (PROVED).** Partition `[X,2X)` into disjoint physical integer intervals
I_j. Supply separate selected anchors, finite native controllers and
thresholds T_j within each interval. Then

\[
 D_X:=\sum_{N=X}^{2X-1}a_N^2
 \le \sum_j\bigl(|I_j|T_j^2+8T_jR_{1,j}+16R_{2,j}\bigr)
 \le \left(\sqrt{\sum_j|I_j|T_j^2}
              +4\sqrt{\sum_jR_{2,j}}\right)^2.          \tag{D2}
\]

The first inequality is D1 on disjoint intervals. The second is Minkowski
on the full row vector, or Cauchy--Schwarz on its cross term. No common
terminal controller is assumed. Each physical row has its own original
source capacities and authenticated paths. Splitting those physical
certificate records retains the same rowwise bound; it need not leave an
anchor inside every newly split subinterval. Reoptimizing with that extra
constraint is a different task.

**Sufficient uniform target, OPEN.** For every epsilon>0 and all dyadic X,
construct such certificates with both

\[
 \sum_j |I_j|T_j^2\le C_\epsilon X^{2+\epsilon},\qquad
 \sum_j R_{2,j}\le C_\epsilon X^{2+\epsilon}.            \tag{U}
\]

D2 then gives `D_X <= 25 C_epsilon X^(2+epsilon)`. Through A3 of the
interval atlas and
[THM-4658, two-octave energy descent](../../01-canon/theorems/THM-4658-collatz-two-octave-memory-energy-descent-critical-lines.md),
this suffices for (W). Constants 7 and 49 are convenient finite thresholds,
not requirements for (W). Even (W) alone does not imply ROOT positivity or
universal convergence. The arithmetic construction in (U) is unproved.

The choice of interval and controller may be adaptive: the proof is a
finite pointwise certificate argument, not a probabilistic union bound.
What cannot be omitted is the cost of every uncovered physical row.

**D3 (PROVED abstract hostile).** A single global threshold is strictly
stronger than an energy bound, even for nonnegative four-Lipschitz data
with an anchor at every row. Let `X=2^(8m)`, `H=2^(5m)` and `L=2^(6m)`,
with m>=1. In the middle of an octave put a plateau of L rows of height H,
linear ramps of slope four and length H/4 on either side, and zero
elsewhere. Set the anchor bound B_a equal to this profile. Its true
square sum is at most

\[
 H^2L+H^3/6\le 2X^2.
\]

For any global threshold T>=H, D1 pays at least `XT^2>=X^(9/4)`.
For T<H, the unselected plateau leaves an interior gap of at least L+1,
whose `16R_2` alone exceeds `(4/3)L^3 >= X^(9/4)`. Thus every global
threshold has a superquadratic invoice. Three physical slabs (zero,
plateau plus ramps, zero) pay at most `H^2(L+H/2)<=2X^2`, using threshold
H only in the middle slab and every row as an anchor. This is an abstract
scalar hostile, **not a claimed Collatz profile**. It explains precisely
why heterogeneous thresholds can avoid a loss that additional samples
alone cannot repair.

**D4 (PROVED adaptive compression).** Suppose a supplied controller gives
residual B_N at every row of the octave and `E=sum B_N^2`. Partition into
bins of length b dividing X. Select an actual minimum-residual row in
each bin, with its original controller and full clocks. Then

\[
 Q^2:=\sum_I b\min_{N\in I}B_N^2\le E,\qquad
 R_2\le\frac{X(b-1)(2b-1)}6,
\]
\[
 D_X\le\left(\sqrt E+4\sqrt{X(b-1)(2b-1)/6}\right)^2.  \tag{D4}
\]

The first bound is the minimum-versus-mean inequality. For any selected
row in a bin, its distance-square sum is maximized at an endpoint, where
it is `sum_{t=0}^{b-1}t^2`. Apply the same vector triangle inequality as in
D2. This compresses an authenticated all-row certificate to X/b anchors;
it does not establish the fine energy E independently. The selector is
finite and deterministic (choose the first minimizer), and has no ROOT
query. If only E is saved, D4 proves existence and the bound, but does not
supply the selected row indices.

## 4. Native finite controls and the first aggregation failure

[Program](../../04-computation/experiments/collatz_anchor_holes_20261017a.py).
The inherited schedule is fixed before the new scale tests:
`b=2^floor(k/2)`, lower-middle anchors in bins of length b, and stage
`floor(k/2)+2+(anchor index mod 3)`. Original-source replay and each new
terminal's literal extension are inherited from
[the audited repair](../../04-computation/experiments/collatz_interval_atlas_repair_20261016e.py).
The new runs do not inspect t(n) or use known convergence to select pairs.

For each slab, enumerate every distinct authenticated B_a as its threshold,
retain anchors with B_a<=T, and minimize the exact integer D1 score. This
is complete only over that declared finite threshold menu. The 2/4/8-slab
menus were introduced after the single-threshold failure at k21; they
are exploratory repairs, not preregistered universal predictions.

| k | Native sharp cone / X^2 | One threshold | Two slabs | Four slabs | Eight slabs |
|---:|---:|---:|---:|---:|---:|
| 18 | 16.66345834 | 28.45040894 | 25.42461646 | 24.59427893 | 24.12415087 |
| 20 | 20.19465803 | 36.24531700 | 32.24232826 | 30.39055576 | 29.37109992 |
| 21 | 29.87986488 | 53.12561712 | 42.28884541 | 39.88954658 | 38.89970609 |
| 22 | 63.13877658 | 111.39989863 | 93.13487271 | 85.04206655 | 82.57042788 |

At k21 the best single threshold is 8754 with 1,585 retained anchors;
its right tail is 2,560 and largest interior gap 6,144. Its energy bound
exceeds 49 X^2, while two slabs succeed without any extra trajectory
steps. The sharp cone already succeeds; the gain repairs a coarse
coverage statistic, not a failure of the native controller itself.
At k22 the fixed mixed-depth schedule itself fails the 7X threshold.
All four threshold menus also fail; they are upper bounds, so this is not
a lower bound on the actual current. The final common stage 15, however,
has exact all-row residual energy `646432378890360`, or
`36.74542647845 X^2`. D4 at b=512 proves existence of 8,192 anchors with
score at most `44.06709420 X^2`; b=256 gives `40.31779944 X^2`. This
repair uses no extra trajectory steps. The selected row indices were not
stored, and are not represented as a delivered list. The 7X comparison is
checked without rounding: for `R=X(b-1)(2b-1)/6`, test
`S=49X^2-E-16R>=0` and `64ER<=S^2`. The native run authenticates every
source through `8X=33554432`, with integer clocks at most 404 and more
than 1.6 billion literal extension-step checks. It is one exhaustive
compositional replay, not two independent large simulations.

Small exact tests enumerate all one-, two- and three-anchor sets on
intervals of lengths 1 through 12 and compare both formulas with direct
row-distance sums. The plateau hostile is independently summed at m=1,2.

## 5. Solve the cubic metric before assembling the graph

Put `B=Q(i,sqrt7,sqrt11)`, `j=i sqrt11`, and retain the real root
`w^3+2w^2-4w-4=0`, `1<w<2`. Set

\[
 u=\frac{w+3-j}{w+3+j},\qquad
 a=\frac{20+11j}{19},\qquad
 u^3+a u^2-\bar a u-1=0.
\]

The cubic remains irreducible over B: the real cubic field has degree
three, B has degree eight, and u recovers w by the prior trace identity.
Use capital A,B_1,C for displacement coefficients to avoid confusing the
middle coefficient with the base field. B5 of the previous note gives

\[
 |A+B_1u+Cu^2|=1
\]

if and only if

\[
 \bar A B_1+C\bar B_1=K:=aC\bar A-\bar C A,             \tag{G1}
\]
\[
 |A|^2+|B_1|^2+|C|^2+\bar a C\bar A+a\bar C A=1.       \tag{G2}
\]

**G1 (PROVED).** If A,C are nonzero and such a unit exists, then
`Delta=|A|^2-|C|^2` is nonzero and its middle coefficient is uniquely

\[
 \boxed{B_1=\frac{AK-C\bar K}{\Delta}.}                \tag{G3}
\]

Conversely, for nonzero Delta, (G3) satisfies (G1), so checking (G2)
completely decides the unit condition. There is no omitted singular
family. Indeed if Delta=0, put `r=C/A`, so `r bar r=1`. The necessary
compatibility `AK=C bar K` becomes
`r^3+a r^2-bar(a)r-1=0`, a root in the base field, contradicting cubic
irreducibility. QED. The same argument works over any conjugation-stable
base with this irreducible cubic, not just our degree-eight field.

Thus an alphabet search need only enumerate the two outer differences;
the middle one is forced. The exact positive control is
`u^3=1+bar(a)u-a u^2`; replacing its middle coefficient by `bar(a)+1`
is rejected. An independent degree-24 quotient algebra checks the metric
on all 64 triples from `{0,1,-1,j}^3` and every accepted candidate.

**G2 (PROVED prime gate).** Suppose A,B_1,C are integral at a
conjugation-fixed prime P over 19. The numerator `20+11j` is a P-unit:
its norm is 1731, nonzero modulo 19. Multiply (G1) by 19 and reduce to get

\[
 C\bar A=0\pmod P.
\]

Hence at least one outer difference vanishes at each such P. In our base
there are four fixed primes: `sqrt7=+/-8`, `sqrt11=+/-7` in F19, while i
lies in F361 and conjugation acts nontrivially. The four vanishing bits
must be retained individually. Having j in the coefficient field removes
the previous field gate but does not pay these valuation conditions.
This necessary gate is not sufficient for (G1)--(G2).

## 6. The hub alphabets, circuits, and boundary tests

Let G be the literal nine-point Gamma9 hub from B7 of the previous note:

\[
 0,1,b,c,d,1+b,1+c,b+d,c+d,\qquad
 b=(\sqrt7-3i)/4,\ c=-i,\ d=(3+i\sqrt7)/4.
\]

The [difference-pair program](../../04-computation/experiments/cubic_hub_coupling_20261017a.py)
tests the three base alphabets `V=G+t jG`, t=1,8,19. Each has 81 points
and exactly 2,401 differences. Every denominator divides four, so G2
applies without qualification. All 5,760,000 ordered nonzero outer pairs
are covered by the gate; every surviving pair then meets three rejecting
finite-field screens (primes 617,641,653), with singular residues retained
whenever compatibility does not reject them. An exact norm calculation
is required for any final acceptance.

| t | After four-prime gate | After modular screens | Three-support unit directions |
|---:|---:|---:|---:|
| 1 | 0 | 0 | 0 |
| 8 | 4,608 | 0 | 0 |
| 19 | 228,096 | 0 | 0 |

For t=1 every nonzero difference has zero vanishing bits. For t=8 there
are 48 differences on each of the complementary masks 0101 and 1010;
this value was chosen after inspecting the branch masks. For t=19 there
are 48 nonzero differences vanishing at all four primes. Neither partial
repair supplies the remaining norm equation. Exhaustive base-graph checks
give 126 edges and chromatic number three in all cases. By the one-or-three
support theorem and cubic independence, each 531,441-vertex cubic product
`V+uV+u^2V` is therefore Cartesian and three-colourable; this conclusion
uses the complete difference search, not enumeration of all product edges.

The [short-circuit assembly](../../04-computation/experiments/hub_short_circuit_20261017a.py)
uses the incoming exact twenty-unit relation, independently checked in our
w-basis. Four seeded orderings (seed 20261017) give 77 partial-sum vertices.
Keep all eight original ports, then add translated copies of G, or of its
four rotations by `1,v_1,v_3,v_4`. The full carrier is
`M=Q(zeta55,w,sqrt7,sqrt11)`, of degree 480. Every candidate edge surviving
three finite-field homomorphisms receives an exact norm test. All old
port pairs are independently retested without a filter.

| Assembly | Vertices | Full induced edges | Mixed-address cubic edges | Proper colouring | Old six-patterns removed |
|---|---:|---:|---:|---:|---:|
| Translated hub | 699 | 1,815 | 0 | 3 | 0 of 570 |
| Rotated hub | 2,523 | 6,932 | 3 | 3 | 0 of 570 |

A mixed-address edge means that no representation of its endpoints
shares a shape coordinate or a partial-sum coordinate. All three rotated
mixed edges have cubic displacement, but none has a displacement outside
L'. Two touch vertices outside L'; their new radical coordinates cancel.
They therefore do not yet consume the hub's new geometry. Each of the
570 restricted-growth six-colour source patterns has its own SAT extension,
replayed against every edge and the literal fixed port colours. Exactly-one
colour clauses are included; a proper graph colouring alone would not
verify a prescribed port assignment.

### Paying the obstruction constructively

The failure of the integral alphabets suggests allowing the very
coefficients forced by G3. With b the hub's unit direction above, set

\[
 V_* = G\cup\{b\bar a,-ba\},\qquad
 \mathcal V=V_*+uV_*+u^2V_*.
\]

**G3 (PROVED construction).** Both 0 and `bu^3` lie in this grid, because

\[
 bu^3=b+(b\bar a)u+(-ba)u^2.
\]

Their difference is unit, uses all three cubic coordinates, and retains
the even hub extension: b is outside L' while u is a nonzero element of
L'. This is a genuine hub-dependent cubic cross-edge. The denominator-19
coefficients deliberately leave G2's integral-coefficient hypothesis;
G3's exact metric, rather than a weakened residue criterion, authenticates
the resulting edge. The old ports can be adjoined unchanged.

The [valuation-grid program](../../04-computation/experiments/hub_valuation_grid_20261017a.py)
checks the entire induced graph after adjoining those ports: **1,337
vertices, 5,096 edges, six three-coordinate cross-edges, all six retaining
the new radical displacement**. The whole graph has chromatic number
three and every one of the 570 old six-colour patterns extends. Coupled
edges are replayed in both the independent u-quotient and w-basis
arithmetic. This crosses the field/valuation barrier but does not yet
produce a forcing context. The six mixed edges have only two unoriented
directions, `bu^3` and `b/u`; the second follows from the reciprocal
identity `u^-1=u^2+a u-bar(a)`. A relation-supported chord is not by
itself a colour obstruction.

The next declared extension takes the four hub directions `1,b,c,d`
and adds `h bar(a), -ha` for each h. Its 17-point base gives a graph with
**4,919 vertices, 12,170 edges and 24 coupled edges in eight unoriented
directions**; eighteen of those edges retain the new radical displacement.
It too has chromatic number three and removes zero of the 570 old patterns.
There is an explicit explanation: colour the nine original base points
`(0,1,1,1,1,0,0,2,2)` and every added coefficient point zero. Colour
`A+B_1u+Cu^2` by the sum of its three base colours modulo three. The
program checks this rule against every induced grid edge, then extends it
to the six old ports outside the grid by a 729-case search. This supplies
an independent constructive check of the SAT three-colouring. Every
subset of these four direction enrichments is a subgraph, so repeating
that finite recipe cannot produce a four-chromatic obstruction.

The next port-sensitive target is consequently a coupled edge or gadget
that frustrates this additive base-colour rule and constrains the new
coefficient colours through actual adjacency. A new algebraic unit alone
is insufficient: here the denominator-bearing coefficients remain free
to take colour zero. Forcing their colours requires compatible attachment
contexts with the original port assignments retained.

## 7. A local-solubility correction with the global theorem preserved

The incoming independent audit's claim that -11 is square at *every*
finite F4-place is **REFUTED**. For
`F4=Q(tau)`, `tau^4-55tau^2+605=0`, modulo 19 the four roots are
`1,4,15,18`, all simple. Hence the completions there are Q19, where
`-11=8 mod 19` has Legendre symbol -1. This is an exact minimal hostile
to that universal statement, not a failure of local solubility of C3.
At D=10, `(z0,z1,z2)=(4,6,6)` solves all three defining C3 equations
modulo 19, with nonzero derivatives in the z coordinates. These finite
points Hensel-lift.

The valid full proof was already given in B6 of
[cubic adelic coverage](cubic_adelic_20261016e.md): infinity points at 2,3,7;
at all other finite places use good reduction and the even-valuation
Kummer image of `(33,33)`. Thus everywhere local solubility, the nonzero
Sha[2] class, index two, and odd-extension nonexistence remain proved
with the repaired
[THM-4662, fish-curve descent](../../01-canon/theorems/THM-4662-the-fish-curve-has-no-point-over-the-real-quartic-field-so-the-fish-and-ei21-are-not-in-q-zeta55-or-any-odd-carrier.md).
The audit report and its program's overbroad comment are corrected in
place, with this lineage. Checking only the bad primes never justified
an assertion about a chosen chart at every good prime.

## 8. Connections that retain the information needed by the consumer

| Source -> target and map | Preserved predicate | Lost information / required sidecar | Cheapest decisive test |
|---|---|---|---|
| Native anchor certificates -> physical hole moments, D1 | Original-row signed cancellation and integer clock | Anchor count loses tails and residual height | k21 single threshold versus two slabs |
| Threshold families -> graded partition, D2 | Pointwise coverage of every row | A scalar global cutoff loses location of expensive rows | Plateau has O(X^2) energy but global invoice >=X^(9/4) |
| Cubic edge -> outer coefficient pair, G3 | Exact unit norm, with forced middle coefficient | Equal outer norms appear singular; irreducibility excludes that branch | u^3 positive control and perturbed middle negative control |
| Integral differences -> four residue masks, G2 | Necessary vanishing on every fixed prime | Rational prime label loses which branch is paid; masks lose norm | t=8 complementary branches still fail the remaining screens |
| Exact unit circuit -> rotated hub assembly | Literal coordinates, ports, full induced edges | Word closure alone loses boundary restrictions and new-radical participation | 3 cubic cross-edges, all 570 old patterns still extend |
| Torsor -> local points | Local solubility | A uniform chart need not work, and Sha retains global obstruction | 19 has finite points but no stated infinity points |

These are explicit maps with different predicates, not an identification
of Collatz paths with unit-distance walks. The common useful move is to
retain the coordinate on which the next certificate actually depends,
then eliminate only coordinates whose recovery is proved. The planar
Jacobian analogue is the existing inverse-clock warning: an invariant
curve or a rational target reencoding does not recover the marked source
Jacobian. This session makes no new planar-Jacobian claim.

The remaining research obligation is arithmetic: prove (U)'s graded
weighted coverage costs uniformly, and construct a fixed-port graph whose
new-radical cubic edges restrict its six-colour boundary. Neither the
prime masks nor the sample-energy spectrum settles those obligations.


## 9. Reproduction and limits of the receipts

[Receipt manifest](../../04-computation/experiments/anchor_hub_20261017a/manifest.json)
records SHA256 hashes of the programs, their local dependencies and all
saved outputs, together with the environment and incoming Git revision.
Native k21/k22 outputs retain the original wrapper's format; the final
wrapper adds the proved postprocessing fields. Compare their native
anchors, stage energies, endpoint/clock hashes and extension audits when
reproducing, rather than expecting identical whole-file bytes.

```sh
python3 -B 04-computation/experiments/collatz_anchor_holes_20261017a.py --k 21
python3 -B 04-computation/experiments/collatz_anchor_holes_20261017a.py --k 22
python3 -B 04-computation/experiments/collatz_anchor_holes_20261017a.py --receipt 04-computation/experiments/anchor_hub_20261017a/native_k22.json
python3 -B 04-computation/experiments/cubic_hub_coupling_20261017a.py --mode jhub
python3 -B 04-computation/experiments/cubic_hub_coupling_20261017a.py --mode jhub8
python3 -B 04-computation/experiments/cubic_hub_coupling_20261017a.py --mode jhub19
python3 -B 04-computation/experiments/hub_short_circuit_20261017a.py --mode translated
python3 -B 04-computation/experiments/hub_short_circuit_20261017a.py --mode rotated
python3 -B 04-computation/experiments/hub_valuation_grid_20261017a.py --mode one
python3 -B 04-computation/experiments/hub_valuation_grid_20261017a.py --mode four
```

Normal and `python3 -O` outputs are compared for the finite geometry
searches and the new k21/k22 postprocessing; runtime checks use explicit
exceptions, not removable assertions. This tests optimization robustness,
not mathematical independence. Independence comes from the alternate
quotient arithmetic, direct hole sums, the explicit product colouring,
and literal source-clock replay. No second large native replay is claimed.
The receipt retains failed and successful policies, all full induced edge
lists, graph colourings, old port coordinates where applicable, native
source/clock digests, selected threshold-anchor sets, and the distinction
between an existential adaptive selection and an exported list of rows.
