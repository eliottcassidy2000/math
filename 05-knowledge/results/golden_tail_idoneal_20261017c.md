# Coverage-tail debt, the norm-19 carrier, and quadratic orders

2026-10-17, clock_duality_oct15. **PROVED, scoped and self-audited:** the
tail-debt implication and peak bound; the field/incidence arguments below.
**FINITE-EXACT:** explicit colouring witnesses prove preservation of all570
old boundary patterns in the enlarged degree192 carrier. **CITED:** the
idoneal/GRH results, with the order-versus-field distinction retained.
**OPEN:** the uniform native tail debt, (W), independent ROOT positivity,
universal Collatz convergence, and six-colour forcing. No independent
referee or priority claim. A closed proposed route is not a solution of
the ambient problem.

## 1. Inheritance, incoming work, and the concept board

Anchor: the uniform Collatz coverage-tail estimate. Niche: mixed golden/cubic
attachments at19. Wildcard: quadratic-order genus quotients and the possible
extra idoneal number. The board is **literal source / integer clock /
unpaid tail / marked port / residue fibre / quadratic conductor**.

The closest proved mechanisms are C1--C4 and G1--G7 of
[the graded-tail packet](graded_tail_carrier_20261017b.md), A1--A3 of
[the interval atlas](collatz_interval_atlas_20261016e.md), and Lemma R of
[THM-4558, residue colourings](../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md).
Hostiles: a large constant terminal coordinate; the fixed-stage k22 failure;
and a residue quotient identifying a golden neighbour with the marked3/4
port. Corrected near misses: replacing a field by its larger advertised
carrier, or a quadratic order by its maximal order. The least-used sidecars
are the actual height of a coverage defect and the ability to recolour
one literal point inside a residue fibre.

The live remote was fetched at startup and during the experiment; there
were no incoming commits beyond5285ba6006. The recent
[S38 explicit-formula note](collatz_explicit_formula_dense_census_20261016.md)
was read, including its failed energy-constant test at2^41 and its conditional
THM-4660. Nothing here promotes that conditional amplitude to an unconditional
one. The recovered
[THM-4566, idoneal planes](../../01-canon/theorems/THM-4566-hilbert-class-field-planes-are-the-idoneal-planes.md)
already supplies the genus-field connection; its global completeness clause
is explicitly relative to the repo's separate #003 zero-free assumption.
The present arguments need neither that assumption nor completeness.

Method cards used: search the statement before the method; retain local
incidence through counting; expose the obstruction before choosing a scale;
and correct the object before sharpening the technique. No new general
method card is claimed.

## 2. Pay the tail that survives, with its physical rows retained

Keep the original weights w(n)=3[n=2 mod3]-1, the ROOT-supported sign t,
and a_N=H(4N)-H(N). Each native anchor has an authenticated literal matching,
actual source labels, and full integer clocks. At scale X let

\[
 C_j(N)=\min_{(a,B,\text{controller})\in\mathcal A_j}
           (B+4|N-a|),\qquad X\le N<2X,
\]

where the finite anchor banks are nested. Thus old certificates survive
every subsequent refinement even when newly selected anchors or terminal
controllers change. The interval atlas transports the selected matching to
the original annulus at N. This is a certificate envelope, not an evaluation
of eventual ROOT membership.

Write U_j(t)=#{N:C_j(N)>=t} and K_j=max_(t>=1)t²U_j(t)/X².
Define the **unpaid half-tail debt**

\[
 A_j=\max_{t\ge1}\frac{t^2}{X^2}
          \left(U_{j+1}(t)-\tfrac12 U_j(t)\right)_+.       \tag{1}
\]

**T1 (PROVED conditional implication).** For every such finite chain,

\[
 K_{j+1}\le\tfrac12 K_j+A_j,\qquad
 K_J\le2^{-J}K_0+\sum_{j<J}2^{-(J-1-j)}A_j.             \tag{2}
\]

Proof: multiply the defining inequality for each height by t²/X² and
take the maximum; then iterate. This standard distribution-function
calculation is not claimed as a new analytic inequality.

Append the trivial anchor(X,4X), so C_0<8X and K_0<=64X. If a native
construction supplies J=ceil(log2 X) finite refinement rounds with
A_j<=A_X, then K_J<=64+2A_X. Consequently an independently proved
A_X=X^{o(1)} would establish the previous packet's critical tail target,
and hence the dyadic discrepancy estimate(W), using its logarithmic
layer-cake bound and
[THM-4658, two-octave descent](../../01-canon/theorems/THM-4658-collatz-two-octave-memory-energy-descent-critical-lines.md).
Neither(W) nor(2) alone establishes independent ROOT positivity or convergence.

This formulation pays rare high defects rather than demanding that every
row fall below a fixed multiple of sqrt(X). It does not supply the missing
arithmetic estimate. With unrestricted choice of a successful finite
controller, such a debt hypothesis can simply restate the tail target.
A meaningful next strengthening is a specified native policy, for example
choosing among2,4,8 further reset stages, whose debt is bounded independently
of the already desired conclusion. This policy claim remains OPEN.

### Exact finite continuation, including the hostile tests

The computation reuses the authenticated native arrays and individual
anchor receipts from the preceding packet. It retains every earlier anchor
when adding the next stage. No extra trajectory steps or assumed ROOT signs
are generated. The maximum in(1) is evaluated at every histogram jump;
between jumps the positive count difference is constant and t² increases,
so the upper endpoint is sufficient. The receipt saves the maximizing
height and both literal row counts, not only a decimal score.

| k | Native stage change | Exact-debt decimal A_j | Maximizing height |
|---:|---:|---:|---:|
|18|11 to13|1.73334074|1150|
|18|13 to17|0.84616459|844|
|20|12 to14|1.44252589|2638|
|20|14 to18|0.49092239|1246|
|22|13 to15|7.27550859|12910|
|22|15 to19|2.61188698|6406|

These six comparisons are FINITE-EXACT, not evidence of a settled uniform
constant. In particular the k22 debt is substantially larger. Keeping the
earlier anchors slightly improves some saved cone energies, without changing
any source or clock. The large k22 computation is not rerun or described as
new orbit verification; its input hashes and inherited audit scope are retained.

The physical source sidecar is now explicit. At the maximizing height6406
for k22, stages15 to19, the old expensive set has582 interval components.
Only one has a positive survivor-minus-exit count: the half-open row interval
[4,873,234, 8,388,608), with2,981,841 survivors and533,533 exits. Its positive
debt alone is2.85555228; the other components reduce the signed global debt
to2.61188698. At stages13 to15 the sole positive component is
[5,678,932, 8,388,608), with2,194,535 survivors and515,141 exits. These are
intervals of physical annular rows N, not intervals of terminal labels.
The remaining four tested transitions have2,4,2,2 positive components,
respectively. This supplies a concrete right-boundary renewal target instead
of an unlocated scalar discrepancy. For any height, summing positive excess
separately over old components is a stronger sufficient debt bound; no
component estimate is inferred just from its length.

A harder-looking good-lambda test, U_(j+1)(2t)<=U_j(t)/8 above a cutoff,
was also evaluated. Most passing high-cutoff ranges are vacuous because2t
already exceeds the new maximum. The receipt records those thresholds;
they are not presented as cancellation evidence. A constant profile C=2X
is four-Lipschitz and survives the identity contraction, but forces a linear
cutoff and has debt2X. Thus local Lipschitz control and norm contraction
alone cannot pay(1). Exhaustive small monotone-profile controls check the
recurrence independently of native data.

## 3. Why the tail target must allow isolated peaks

**T2 (PROVED).** Any nonnegative integer four-Lipschitz profile C on X
consecutive rows with H=max C<8X satisfies

\[
                         H^3\le32K_X X^2.               \tag{3}
\]

At a maximizing row, at least floor(H/8)+1 rows of the clipped interval
are within floor(H/8) steps. Their values are at least ceil(H/2).
Substitute this population and height into K_X. This also handles a maximum
at an endpoint; H=0 is immediate.

The exponent is sharp for this information. At X=m³ take
C_j=max(m²-4j,0), 0<=j<X, with m even. Then max C=X^(2/3), while K_X is
bounded and sum C²=O(X²). These are abstract profiles, not claimed Collatz
controllers. Thus the critical tail target permits peaks much larger than
sqrt(X); replacing it by a uniform square-root cutoff is a genuine
strengthening. The legitimate consequence of the tail target is
|a_N|<=X^(2/3+o(1)), in addition to its stronger average-square information.
It is not a certificate that each individual orbit reaches1.

## 4. The mixed field and the norm-19 square root are still closed routes

Use the fixed ports from the preceding packet, in their original order,

\[
0,1,1+\alpha,1+\alpha+\alpha^2,-\alpha^4,
(3-j\sqrt5)/8,3/4,(3+j\sqrt5)/8,
\]

where alpha=zeta5, j=i sqrt11, phi=(1+sqrt5)/2. Put

\[
B=\mathbb Q(i,\sqrt7,\sqrt{11}),\quad
f(w)=w^3+2w^2-4w-4=0,\ 1<w<2,\quad
K=B(w,\phi),\quad r=\sqrt{4\phi-3}>0.
\]

The mixed unit v=(w+3+phi-j)/(w+3+phi+j) lies in K and already recovers
w and phi by G7. The previous recovered noble-polyhedron chart is
(x+1/2)²=phi-3/4, with r=2x+1. The field degrees are48 for K and96 for K(r).

**G1 (PROVED incidence and place).** In K, and also in K(r), the only
neighbours of the external port2 are1 and phi; the only neighbours of
external port4 are0 and1-phi; external port3 has no neighbours in the field.
There is a conjugation-fixed prime over19 with residue field F19[i], where
conjugation is x+iy -> x-iy, and the five internal ports0,1,5,6,7 reduce to

\[
                    (0,0),(1,0),(17,15),(15,0),(17,4).    \tag{4}
\]

Proof of incidence: alpha remains quadratic, with its nontrivial automorphism
alpha -> alpha^-1 fixing the field and commuting with conjugation. A common
unit neighbour of the two port2 conjugates must be real. The circle equation
has roots1 and phi. The port4 equation gives0 and1-phi. The port3 conjugates
have real part1/2 and imaginary-part square phi+3/4>1, so neither has a
real unit neighbour. The exact ten-point check in G7 of the preceding packet
also proves that phi and1-phi have no unit edges to the five internal ports
or to each other.

Here alpha is absent from K because a cyclic quartic cyclotomic field is
absent from the multiquadratic base, and the subsequent cubic extension
cannot acquire a missing quadratic extension. It is still absent from K(r):
at a prime over19 with phi=15, the element4phi-3 has valuation one, so K(r)/K
ramifies; K(alpha)/K is unramified at19. They cannot be the same quadratic
extension. The place is explicit: modulo19 choose sqrt7=8, sqrt11=-7,
w=4, phi=5 and i²=-1. All the real roots lie in F19, so the prime is fixed
by conjugation. At this branch4phi-3=17=6², so r can be adjoined without
enlarging the residue field. At the other golden branch phi=15 it is instead
ramified, again with unchanged residue field and a unique prime. Hensel
lifting applies to the unramified roots; f'(4)=60 is nonzero modulo19.

**G2 (PROVED with finite explicit witnesses).** Every proper six-colour
assignment of the old eight ports extends to the entire unit graph on
K(r) together with those ports. This includes K and the mixed generator v.

Proof: use the inherited five-colouring c(x+y) of the19-residue unit graph,
where

    c = (4,3,2,4,3,2,1,0,4,1,0,4,3,2,1,0,2,1,0).

Its3610 undirected edges are all replayed. Rotate the linear projection to
ax+by+t with a²+b²=2, permute its five colours, and reserve a sixth colour.
There are exactly12 proper equality patterns on the five internal ports.
The saved12 templates realize each by recolouring only the literal marked
points assigned the sixth colour. Those points are independent because the
prescribed pattern is proper; no other point uses that colour.

Each golden contact can now keep its base colour or, if it conflicts with
its external port, take the reserved colour. These two point repairs cannot
conflict with any repaired internal port, by G1's exact incidence. Port3
has no new neighbours. The script directly replays this construction for
all570 old patterns, not merely for an unlabelled graph.

For completeness of the infinite lift, every unit difference is integral at
the fixed prime: v_P(u)+v_P(bar u)=2v_P(u)=0. Colour the valuation ring by
its residue, with the stated finite point exceptions. Distinct additive
cosets of that valuation ring have no unit edges, so colour the other cosets
by translated residue maps. This does not assume that finite rhombus closure
has already connected the whole valuation ring. QED.

The old residue-only hostile is repaired by information *inside* the fibre.
Ramifying at19 does not itself destroy the residue colouring. This is the
precise reason the proposed norm-19 escape is insufficient.

## 5. Even adding the pentagon direction preserves all570 patterns

Let E=B(w,alpha,r), a degree192 field. It contains K(r), all eight ports,
the mixed generator, and the recovered noble-polyhedron chart.

**G3 (PROVED with finite explicit witnesses).** Every proper six-colour
assignment of the old ports extends to G(E). Hence no finite graph confined
to E restricts their six-colour boundary relation.

The eight usable19-residue branches have alpha=(2,±4) or(7,±3), with
sqrt11=±7. In each, conjugation is Frobenius19 and the residue field is
F361. At phi=5 the r-extension splits locally; at phi=15 it ramifies
without changing that residue field. All old ports are integral.

The five-colour-plus-one-point-repair method extends565 patterns. Its five
misses, in the original port order, are

    01021202, 01202121, 01212121, 01232121, 01232323.

These are failures of that template family, not graph obstructions. A
complete enumeration gives exactly19 five-colour partitions of the projected
19-vertex circulant, all translates of c. Enlarging its palette to six,
and keeping the eight marked points separate from their residue classes,
repairs every miss. For each pattern the script saves an explicit19-entry
colour table, branch, and projection. It verifies the lift on all3610
residue edges, all20 residue neighbours of each literal port, and the9 old
port edges. The five patterns have respectively8,16,24,32,32 positive
branch/projection templates among160 tested each.

The construction colours generic points by the19-entry table and assigns
the prescribed colours at the eight literal ports. The table forbids each
port's colour on all residue neighbours of that port. Thus all actual edges
to every unmarked field point are proper. Edges between marked ports are
checked in their actual graph. Translate the construction on the other
valuation-ring cosets as in G2. This proves the infinite statement from
finite positive witnesses, without inferring planar edges from residue edges.

An earlier bounded SAT search on the369-vertex split-residue graph returned
UNKNOWN on all40 attempted pattern/branch cases. The saved exploratory
receipt retains this result. The successful19-state search is a constructive
reroute around that search obstruction; no UNKNOWN result was called UNSAT.

### A genuinely different next carrier

Adjoining a real sqrt3 destroys each inherited fixed19 prime:3 is nonsquare
in F19, while its roots in F361 are exchanged by Frobenius19; complex
conjugation must fix the real sqrt3. Thus no extension of that prime can
still have conjugation equal to its old residue involution. More precisely
the polynomial splits over F361 and conjugation exchanges the two primes.
This removes this specific five-colour certificate; other primes must still
be audited, and no chromatic lower bound follows.

Alternatively, adjoining zeta11 over the existing j has degree5. Since
ord_11(19)=10, the residue field above19 enlarges from F19² to F19¹⁰;
conjugation becomes Frobenius19^5. The old361-state certificate no longer
applies. This gives an exact reason to return to eleventh-root directions,
rather than add more points generated by v inside E. The required next
object is still an authenticated finite boundary gadget in an escaped
carrier, with its actual unit cross-edges and original port labels.

A first exact place audit of E(sqrt3) checked every unramified odd rational
prime below5000 for a conjugation-fixed place with residue F_(p²). The
qualifying primes are1319,1559,1979,4259,4919. This search explicitly
excludes the bad primes2,5,7,11,19,37 and does not classify higher residue
degrees. Its criterion is p=3 mod4, p=-1 mod5, quadratic residues3,7,11,
a root of f in F_p, and a golden root with4phi-3 square in F_p. The first
branch has w=1112 and phi=400 modulo1319. This is a finite place inventory,
not a proof of a chromatic lower bound. The actual unit triangle
0,v,eta v, eta=(1+i sqrt3)/2, gives a minimal positive geometric control
for the carrier change; unit triangles are standard and do not themselves
force six colours.

## 6. Idoneal numbers: the conductor explains11,7,8,15,21

An idoneal n is characterized by exponent at most2 of the proper class group
of the **order** with discriminant-4n. The genus map is
Cl(-4n) -> Cl(-4n)/Cl(-4n)². It loses no class information exactly in this
exponent-two case. This is a precise quotient-completeness property, not
a statement that multiplication or the number line itself has exponent two.

The independent reduced-form enumeration gives:

|n|Order discriminant|Field discriminant|Conductor|Order class group|Idoneal?|
|---:|---:|---:|---:|---|---|
|7|-28|-7|2|trivial|yes|
|8|-32|-8|2|C2|yes|
|11|-44|-11|2|C3|no|
|15|-60|-15|2|C2|yes|
|19|-76|-19|2|C3|no|
|21|-84|-84|1|C2 x C2|yes|
|37|-148|-148|1|C2|yes|

For11 the three forms are(1,0,11),(3,-2,4),(3,2,4); for19 they are
(1,0,19),(4,-2,5),(4,2,5). Their genus quotient is trivial but their class
group is not. Both underlying fields have class number one. The standard
quadratic order class-number formula at conductor2 multiplies the field
class number by2-(D/2), when the unit index is one: the multiplier is1 for
D=-7,-15,2 for D=-8, and3 for D=-11,-19. Thus11 and19 expose exactly the
extra coordinate discarded by replacing an order with its field.

This is also why the two occurrences of65 must not be identified: the
known idoneal orders and the known fundamental discriminants with one
class per genus are different lists. See the explicit correspondence in
[Watkins's note](https://magma.maths.usyd.edu.au/~watkins/papers/IDONcorr.pdf).
The table above is reproduced here from exact reduced forms, not inferred
from those list lengths.

Without the repo's adopted #003 input, the classical possible66th example
has a precise caveat. [Kani, Corollary23 and
Remark24](https://www.labmath.uqam.ca/~annales/volumes/35-2/PDF/197-227.pdf)
allows at most one new squarefree n>1848: either a single n=1 mod4, or a
pair n,4n with n=2 mod4. GRH excludes both cases. The relevant functions
are quadratic Dirichlet L-functions; classical RH for zeta alone is not
the stated hypothesis. [Elsenhans--Kluners--Nicolae, Theorems2 and4](https://arxiv.org/pdf/1803.02056)
give the corresponding field-exponent completeness under absence of Siegel
zeros, with at most one exceptional field unconditionally.

**Do not regress the inherited status:** under the #003 result already
adopted in THM-4566, the65-member list is complete and there is no66th
example. The preceding paragraph explains the classical dependency with
that input removed; it does not reopen the repo's completed implication.
The current [primary formalization scope](https://raw.githubusercontent.com/openai/math/main/lean/docs/003.md)
states the7/8 half-plane for every Dirichlet L-function. This session did
not rebuild that external Lean development or audit its entire analytic
proof. Its status here is CITED/adopted input, not newly verified work.
Likewise the current [Erdos508 statement](https://raw.githubusercontent.com/google-deepmind/formal-conjectures/main/FormalConjectures/ErdosProblems/508.lean)
records the externally proved lower bound6 and upper bound7; the search in
this packet concerns the remaining six-colour boundary obstruction.

THM-4566 already maps the squarefree1 mod4 slice to Hilbert-class-field
planes. A hypothetical new odd member would add a field to that particular
family; it would not by itself specify its chromatic number. A hypothetical
even pair belongs to a different slice. No implication between existence
of an extra idoneal number and chi(R²)=6 or7, or Collatz convergence, is
proved here or by the cited theorems. The shared methodological question is
whether a quotient retains the exceptional object needed by the consumer.

## 7. Arithmetic operations and what their changes of representation retain

The [squares/doubles foundry](collatz_sqdbl_20260925_squares_doubles_foundry.md)
already supplies the strongest relevant exact dictionary: n -> q^n in an
infinite cyclic group sends addition to multiplication, evenness to being
a square *inside that group*, and3n+1 to q(q^n)³. Its off-diagonal
multi-prime construction even diverges; a lossless change of notation is
not a new descent proof. This work is recovered rather than reproved.

Although addition and multiplication are individually commutative, their
actions do not commute. For A_b(x)=x+b and M_a(x)=ax,
M_a A_b=A_(ab) M_a. Their order retains the translation(a-1)b. In Collatz
words the corresponding affine numerator is load-bearing, in addition to
the multiplier and the parity gate. This is the same reason a quadratic
field name cannot replace its conductor, or a residue cannot replace a
marked vertex. These are comparisons of information loss, not conjugacies
between the three problems.

The golden ring offers a second exact, differently typed table. The order
of phi in(Z/mZ)[phi], phi²=phi+1, for m=7,8,11,15,19,21,37 is respectively
16,12,10,40,18,16,76. The composite cases retain prime-power digits and
combine by CRT. At19 the scalar roots5 and15 have orders9 and18, and only
phi=15 makes w+phi=0 at w=4. This is an actual link between the norm-19
branch and a doubled residue clock. It neither identifies these unit groups
with the quadratic ideal class groups nor transfers a clock bound to
Collatz. The factor3 in the conductor11 order is also not our cubic field:
x³-x²+x+1 has discriminant-44, whereas f has polynomial discriminant592
and its integral field discriminant is148.

## 8. Connection contracts and remaining obligations

|Map|Preserved predicate|Lost data / repair|Decisive result|
|---|---|---|---|
|Native matching -> physical tail|Original weighted annulus and full clocks|A mean alone loses high-level populations|Exact debt(1); uniform native bound OPEN|
|Class group -> genus quotient|All genus characters|Cl² and conductor|The C3 examples at11 and19 are invisible|
|Field point -> residue|Unit edges at a fixed prime|Literal port inside its fibre|Point exceptions repair all570 patterns in E|
|Square-root chart -> golden displacement|The exact algebraic displacement|Sheet and unit-distance metric|Ramification at19 preserves the colouring obstruction|
|Adjoin sqrt3 or zeta11|Original field points and unit edges|Old fixed prime or old residue size|An explicit escape condition, no forcing yet|

The next Collatz test should identify arithmetic structure at the saved
maximizing debt heights, with original row intervals retained, rather than
fit another uniform constant from six values. The next geometry search
should operate beyond E and first replay the full boundary relation. An
extra idoneal number is neither an available oracle for those tasks nor a
reason to erase their separate quantifiers.

## 9. Reproduction and audit boundary

Scripts and receipts are under
`04-computation/experiments/golden_tail_20261017c/manifest.json`.
The normal/optimized replay, independent residue-witness verification,
input hashes and exact integer tail scores are recorded there. Native
array authentication is inherited from the preceding packet's C1 audit.
All finite colour conclusions use positive, explicitly replayed colourings;
no SAT-unsatisfiability claim is needed. Literature claims above are CITED;
their hypotheses are not adopted as unconditional new canon.
