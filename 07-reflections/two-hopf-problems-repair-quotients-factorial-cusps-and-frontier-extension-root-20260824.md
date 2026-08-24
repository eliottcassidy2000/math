# Two Hopf problems, repair quotients, factorial cusps, and the next mathematical frontiers

> **SYNTHESIS / RESEARCH REFLECTION, 2026-08-24.** Canon links below control
> truth.  The two source manuscripts remain **PREPRINT CLAIM / UNDER AUDIT**.
> No analogy in this document promotes either headline theorem, proves
> `LRC(14)`, or proves `JC(2)`.

## Executive result

The productive common theme is not “Hopf” and not a theorem-level bridge
between Riemannian and complex geometry.  It is a proof grammar:

```text
symmetric degenerate scaffold
  -> complete the residual zero stratum
  -> minimize or normalize away transverse/gauge variables
  -> identify the first residual class modulo allowed repairs
  -> retain branch/address/orientation data before taking a quotient.       (1)
```

Six unconditional results came out of that grammar.

1. [THM-3990](../01-canon/theorems/THM-3990-componentwise-harmonic-obstruction-and-repair-quotient.md)
   proves the exact componentwise Laplacian repair criterion and its finite
   graph/Farkas analogue.  The real repair quotient can be torsion-free while
   the integral graph cokernel retains torsion; the triangle already has
   Smith factors `(1,3,0)`.
2. [THM-3991](../01-canon/theorems/THM-3991-periodic-unimodular-toric-cusp-factorial-euler-obstruction.md)
   proves the periodic unimodular fan law
   `chi(W)=d*n!`.  A unique homology-sphere fibre in this grammar forces
   `d*n!=2`; irreducibility forces complex torus rank `n=2`.
3. [THM-3992](../01-canon/theorems/THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual.md)
   collapses the first reduced `2:3` planar cusp-pole cell to an exact mixed
   residual.  Polynomiality forces `h=gamma*s`, and the deleted source line
   maps to a centered node, not a cusp.  The two normalization addresses are
   now an oriented conductor-incidence problem.
4. [THM-3994](../01-canon/theorems/THM-3994-double-resultant-collision-separates-two-address-and-length-two-seams.md)
   proves that the two double-resultant seams excluded by THM-3972 have
   different completed geometry: `cr=3` is two reduced transverse addresses,
   whereas `4cr=3` is one curvilinear length-two centre whose Rees graph has
   an `A1` singularity and local class group `Z/2`.
5. [THM-3995](../01-canon/theorems/THM-3995-scale-two-parity-hole-support-and-integer-variance-tariff.md)
   turns the four oriented endpoint holes in the scale-two LRC row into an
   exact support cap, its sharp integer variance tariff, and a new sufficient
   gate.  It is a real narrowing of the row, not a proof of `LRC(14)`.
6. [THM-3996](../01-canon/theorems/THM-3996-etale-node-address-balance-cycle-and-nonproperness-dichotomy.md)
   proves the complete node-address conservation law.  Away from the
   nonproperness locus, each normalized component has equal incoming and
   outgoing address degree, every address lies on a directed cycle, and the
   two THM-3992 clutches force either an additional address or a Jelonek
   value.  Distinct companion owners put the extra address in their packet;
   a two-edge common-owner cycle can only be a proper subpacket.

The source ledgers remain deliberately stricter:

- [Brendle--Hung `S2 x S2` ledger](../05-knowledge/reference/CORE-PAPERS-BRENDLE-HUNG-S2XS2-2026-08-24.md):
  the omitted `z(h_a)=0` branch is now independently **FINITE-EXACT**, and a
  sensitive generic point gives exact cancellation in `V_bc`; the global
  `V_bc` identity and the full curvature argument remain open referee gates.
- [Hopf/`S6` ledger](../05-knowledge/reference/CORE-PAPERS-HOPF-S6-2026-08-24.md):
  the displayed matrices, Smith forms, conductor pushout, and twist
  presentation have exact audits.  The infinite-fan analytic quotient,
  attaching maps, nearby cycles, and global recognition remain open.

## 1. Inheritance pass and portfolio

| lane | closest proved mechanism | canonical hostile | corrected near miss | least-used sidecar |
|---|---|---|---|---|
| Anchor: `LRC(14)` | THM-3910 tariff equality and sheet variance | AP13/Goddyn--Wong owner masking | endpoint residue without orientation | signed exposed-event chain |
| Anchor: source audits | compact residual-stratum minimization; exact lattice Smith | mutable notebook state; reversed conductor branch | raw gluing-coordinate noninjectivity | completed stratum / presentation gauge |
| Niche: planar `2:3` cell | THM-3989 Laurent conductor and scalar moment | same nodal boundary with bracket zero | treating early integration constants as obstructions | two labelled node clutches |
| Wildcard: cusp fan rank | toric orbit Euler additivity | cell counts without stabilizers/orientation | extrapolating the `A2` cusp to arbitrary rank | quotient index and simplex orbit |
| Boundary collision | THM-3972 graph blowup normalization | one scalar double resultant | multiplicity as place count | scheme length distributed by address |

The wildcard overtook the anchor long enough to become THM-3991; the niche
overtook it long enough to become THM-3992 and THM-3994.  `LRC(14)` remains
the principal open anchor.

## 2. The two Hopf problems are mathematically distinct

Brendle--Hung's [arXiv:2608.19068](https://arxiv.org/abs/2608.19068) concerns
positive sectional curvature on `S2 x S2`.  The manuscript at
[alpo.ge/s6.pdf](https://alpo.ge/s6.pdf) concerns an integrable complex
structure on `S6`.  Neither problem supplies a theorem about the other.

The Riemannian construction starts from a nonnegative-curvature metric,
minimizes the perturbed curvature over nearby two-planes, proves a positive
quadratic coefficient away from a completed torus `Sigma`, and tries to make
the first cubic residual positive on `Sigma` by a conformal Laplacian repair.
The analytic-complex construction starts from a rank-four triangle-group
representation, fills two elliptic fibres by logarithmic transforms, fills a
unipotent cusp by an infinite toric fan quotient, and tries to recover `S6`
from conductor gluing, twists, nearby cycles, and attaching matrices.

Their shared grammar is narrower and useful:

```text
Brendle--Hung:  quadratic zero locus -> completed Sigma -> cubic / im Delta;
S6 manuscript: toric central fibre -> conductor branches -> clutch / gauge;
JC boundary:   index zero locus -> normalization addresses -> divisor / principal;
LRC:           tariff equality locus -> exposed owner events -> defect / response.
```

The preserved object is a residual class modulo a lawful repair operator.
Positivity, complex integrability, polynomial invertibility, and loneliness
are all destroyed by the transfer and must be reintroduced independently.

## 3. The concept board

The live board had seven objects.

1. **Completed residual zero stratum.** Open charts are not the invariant;
   completion determines connected components, boundary flux, and interface
   conditions.
2. **Repair quotient or adjoint cokernel.** The obstruction is not a raw
   coefficient but its class modulo allowed changes.  For an operator `A`,
   the dual obstruction lives in `ker(A^*)`.
3. **Conductor/address incidence.** Normalization forgets which branches were
   glued.  A labelled dual graph is required before a class-group quotient.
4. **Periodic fan Euler budget.** A fundamental simplex count becomes a
   factorial obstruction after orbit additivity.
5. **The `2:3` cusp-log cell.** Early Laurent rows are target gauges; the first
   real residual occurs later and forces a node.
6. **Signed LRC owner events.** Unoriented endpoint coverage is too coarse;
   the chain boundary remembers which owner wall creates each jump.
7. **Rational elliptic/Picard--Fuchs spine.** The `3,4,infinity` local
   monodromies suggest a much smaller search space than arbitrary integral
   matrices.

The useful board updates were:

- completing `Sigma` turned two apparent open pieces into one torus and one
  mean obstruction;
- quotienting the `S6` twist presentation turned a rank-two raw kernel into
  standard Seifert gauge, retracting a false moduli claim;
- retaining address distribution split one double resultant into two smooth
  blowups versus one `A1` graph singularity;
- retaining the scalar moment in the `2:3` Laurent system converted a free
  integration constant into `h|s`, hence `h=gamma*s`;
- retaining LRC orientation sharpened `u congruent +/- vt` to
  `u congruent -vt (mod 14)`, but owner masking showed that this still is not
  the full event sidecar.

## 4. Repair quotient: the exact transferable theorem

Let a compact residual locus be a disjoint union of closed connected
components `Sigma_j`.  THM-3990 proves that a smooth correction `chi` with

```text
epsilon*(w-Delta chi)>0
```

exists exactly when every component average of `w` has sign `epsilon`.
Indeed one can solve the Poisson equation separately and flatten `w` exactly
to those averages.  This cleanly identifies the obstruction:

```text
C^infinity(Sigma) / im Delta  ->  R^(number of components).               (2)
```

Three boundaries matter.

1. If components meet, the correct operator includes conductor/interface
   conditions; branchwise averages alone are incomplete.
2. For weighted, equivariant, or singular operators, constants need not span
   `ker(A^*)`.
3. Over integers, the real image loses Smith torsion.  A triangle graph has a
   real Laplacian solution for `(1,-1,0)` while the integral cokernel retains
   order three.

This is the precise reason the same proof move can feed the three frontiers
without identifying their conclusions:

```text
curvature: smooth Laplacian cokernel;
JC:        principal-divisor valuation cokernel;
LRC:       integer response-matrix cokernel;
S6:        conductor/attaching lattice cokernel.                         (3)
```

## 5. The factorial toric obstruction and its escape hatches

THM-3991 considers a periodic rank-`n` unimodular simplicial fan with
translation quotient of index `d`.  Every top-dimensional fundamental simplex
contributes one torus fixed point, and a unimodular fundamental cell has `n!`
oriented simplex chambers.  Orbit additivity gives

```text
number of irreducible components = d,
chi(central fibre)             = d*n!.                                  (4)
```

If this is the unique exceptional fibre in a fibration whose total space is
an integral homology sphere, the Euler invoice forces `d*n!=2`.  An
irreducible central fibre has `d=1`, hence `n=2`.

This explains why the rank-two `A2` cusp is exceptional inside that exact
grammar.  A higher-rank continuation must explicitly leave at least one
hypothesis:

- use nonunimodular cones and account for quotient singularities;
- allow stabilizers and use a Burnside/orbifold correction;
- introduce additional singular fibres whose Euler defects cancel;
- use a non-simplicial fan and retain incidence multiplicities;
- abandon the unique-fibre homology-sphere invoice.

“Try rank three” without naming the escaped hypothesis is now dead.

## 6. The planar quintic route and the forced nodal residual

The older common-zero quintic route is already sharply closed in its stated
normal cubic-cover scope.  [THM-3887](../01-canon/theorems/THM-3887-binary-cubic-common-zero-quintic-one-tangent-obstruction.md)
rules out the seductive single unibranch quintic cusp.  [THM-3890](../01-canon/theorems/THM-3890-universal-quintic-common-zero-resolvent-class-group-dichotomy.md)
then treats every quartic tangent cone, including repeated-root pencils:

```text
Cl(quadratic resolvent) = Z       if the quartic is nonsquare,
                       = Z/5     if the quartic is a square.              (5)
```

Neither group has `3`-torsion, and a would-be tangent escape makes the
quintic factor.  Thus a normal common-zero cubic cover cannot survive at
degree five.  Honest positive targets begin at degree six, at unit-ideal
nonmonogenic forms, or outside that common-zero grammar; THM-3906 and
THM-3907 are the current typed controls.

The different cusp-log route has now contracted just as sharply.  In the
first reduced depths `(2,3)`, THM-3992 proves

```text
h=gamma*s,
C^2-A^3+(3a^2/4)A+a^3/4
  =gamma*u+(3a/(2gamma))*p+R(p,y),       R in (p^2,y).                    (6)
```

On the deleted line,

```text
A=X^2+a,              C=X^3+(3a/2)X,              a!=0.                 (7)
```

The `a=0` cusp would make both `A_x` and `C_x` vanish at one source point,
contradicting the Keller equation.  Therefore `(7)` is a node.  Its two
addresses center the translation, force `b in s^2 k[s]`, and give the
companion factor

```text
G(A,C)=t Q(x,t),
Q(x,0)=gamma*(x^2+3a/(2gamma^2)).                       (8)
```

The two roots in `(8)` are simple transverse clutches.  The remaining
question is discrete and global: are their companion branches on one
component, producing a dual-graph cycle, or on two components, producing a
forest?  A raw coefficient expansion cannot answer that question because it
has discarded component ownership.

The next computation should therefore factor the complete companion divisor,
label both target-node branches, build the exceptional/pole/ramification
incidence matrix, and compute both its integral Smith form and
`H^1(dual graph,mu_3)`.  The quadratic `Q(x,0)` sees only clutches on the known
line, not every source address over the node.  The interior pullback curves
are not the completion-boundary primes of THM-3951, so a two-owner path is not
itself an obstruction.  Outside the nonproper-value locus, finite-etale
address balance forces every edge of the complete graph onto a cycle;
distinct companion owners must therefore pay at least one extra node address.
A complete forest instead detects nonproperness at the node.  This exact
dichotomy is
[THM-3996](../01-canon/theorems/THM-3996-etale-node-address-balance-cycle-and-nonproperness-dichotomy.md).
The generic-degree gap strengthens it: outside the Jelonek locus the complete
fibre has `d>=3` points, so some third address exists even when the two known
clutches already form a common-owner cycle.
On any completed graph, a forest has no clutch holonomy and a cycle is the
only place a new three-primary class can live; existence of that class would
still not make it a Keller boundary basis.  This repairs the overreach
recorded in MISTAKE-469.

## 7. Repeated resultants: a new exact warning theorem

THM-3994 completes the two excluded seams in THM-3972's first-height
constant-`(c,r)` family.  With

```text
D=1-tv^3,             N=3r+3v+cv^2,
Xi=Res_v(D,N)=c^3+27r^3t^2+27(1-cr)t,                                  (9)
```

one has

```text
disc_t(Xi)=-27(cr-3)^2(4cr-3).                                         (10)
```

At `cr=3`, the double root is the sum of two transverse length-one
basepoints.  At `4cr=3`, exact local coordinates give the single base ideal
`(L,V^2)`.  Its Rees chart is

```text
V^2=LZ,                                                                (11)
```

an `A1` singularity with local class group `Z/2`.  Both rows have resultant
order two and total intersection length two.  The scalar sees neither the
support distribution nor the graph singularity.

This is the cleanest transfer of the “complete the zero locus” theme.  It is
also a warning for the confluent degree-eight experiment: one projective point
and a repeated eliminant never imply one normalization place.

## 8. `LRC(14)`: what the repair-operator lens actually bought

On THM-3910's tariff equality locus with `k=floor(m/|B|)>=1`, the sheet count
has only the values `0,k,k+1`.  Its distributional boundary satisfies

```text
t_* partial(1_G)=partial N_t
                =k partial(1_B)+partial(1_E),
E={N_t=k+1} intersect B.                                (12)
```

Thus every genuine boundary of `B` must be the signed image of an actual
body safe-wall event.  For a reduced target endpoint `w=c/d`, jump `j`, and
owner `u`, the individual incidence test is

```text
14 d gcd(t,u) divides 14u c-t j d.                      (13)
```

At a scale-one wall `w=(14b+tau)/(14v)` with `j=-tau`, this becomes

```text
14 divides u+vt,
v gcd(t,u) divides ub+tau*(u+vt)/14.                    (14)
```

The oriented shadow is `u=-vt mod 14`, not an unsigned choice.  Exact
enumeration nevertheless shows that all sixteen scale-one survivor types are
address-feasible: choosing `t=13+14pq` and owners `p,q` realizes every
correctly oriented wall at this coarse level.  Hence endpoint arithmetic
alone cannot close a scale-one type.  The missing data are exposure and
masking in the actual sets `G_u`.

The scale-two survivor `(2,1,9)` behaves differently.  Odd `t` leaves four
interior parity holes of total length `4/(63t)`, so a failure can have support
at most

```text
s_t=4(t-1)/(63t).                                       (15)
```

Writing `m=t*mu(G)` and `theta_t={m/s_t}`, the exact integer tariff strengthens
to

```text
Var(N_t)>=m^2(1-s_t)/s_t+s_t theta_t(1-theta_t).         (16)
```

Together with the variance upper bound this yields the sufficient gate

```text
t*mu(G)/r > sqrt(4(t-1)/(3(59t+4))).                    (17)
```

The individual oriented endpoint equations `(12)--(14)` remain a
**FINITE-EXACT SCOUT RESULT** pending separate canon promotion.  The parity
support and tariff `(15)--(17)` are now **PROVED + VERIFIED-EXACT +
INDEPENDENTLY HOSTILE-AUDITED** in
[THM-3995](../01-canon/theorems/THM-3995-scale-two-parity-hole-support-and-integer-variance-tariff.md).
They do not close the scale-two row or `LRC(14)`.  The faithful next object is
the signed exposed-event matrix, not the individual address matrix.  On
equality rows that survive, the next scalar should be the cubic
falling-factorial response

```text
integral N_t(N_t-1)(N_t-2),                             (18)
```

whose Fourier expansion retains triple phases discarded by variance.

## 9. The rational elliptic spine behind the `S6` source

The displayed local monodromy orders `3,4,infinity` should be searched first
through rational elliptic surfaces, not through arbitrary `4 by 4` integral
matrices.  The Euler-compatible singular-fibre packet

```text
IV* + III + I1,               8+3+1=12,                                (19)
```

is the primitive-width candidate; `IV+III+I5` is a width-five control.  The
most valuable independent reconstruction would:

1. build a Weierstrass model with packet `(19)`;
2. derive its Picard--Fuchs equation and local exponents;
3. compare the resulting integral period lattice with the manuscript's
   rank-four representation;
4. recover, or fail to recover, the half-weight square root and `O(-1)`
   torsor class without using the manuscript's analytic construction.

Agreement would independently certify the finite monodromy spine.  Failure
would localize the discrepancy before the infinite fan and topology enter.
Neither outcome alone decides whether the analytic quotient is compact or
diffeomorphic to `S6`.

## 10. Ranked next frontiers

### A. Global `V_bc` certificate

Rebuild `V_bc` in an immutable exact state, reduce it under the orthonormal
frame relations, and archive a nonzero hostile.  A symmetry proof is cheaper
only if the action on every tensor, plane coordinate, Hessian, and sign is
written explicitly.  This is the shortest remaining load-bearing audit of
the Brendle--Hung quadratic positivity step.

### B. `S6` analytic gauge quotient and attaching maps

The displayed twist scalar

```text
p=12 ell0-4 ell1-3 ell2
```

has kernel generated by `(1,3,0)` and `(1,0,4)`, exactly the standard Seifert
presentation moves; adjoining `(0,-1,1)` gives a unimodular basis on which
`p=1`.  Thus signed `p` is complete for oriented labelled Seifert data, while
raw coordinate noninjectivity does not exhibit analytic moduli.  The exact
centralizer acts by integer shifts on the punctured period family, but the
order-four affine filling permits only even shifts.  The first filled marked
address is `[c0-c2] in C/(2Z)`.  Compute the Cech overlap/cusp class of the
remaining generator `b=2`; only its vanishing can promote the local
conjugacies to a completed analytic equivalence.  Independently reconstruct
the oriented attaching maps and nearby-cycle extensions; Smith forms of
displayed matrices cannot certify their geometric provenance.

### C. Nodal companion completion

For `(8)`, normalize every companion component and record its map to the
normalization of the target nodal cubic.  Proper components should be tested
against the fact that affine line has no nontrivial connected finite etale
cover in characteristic zero.  By
[THM-3996](../01-canon/theorems/THM-3996-etale-node-address-balance-cycle-and-nonproperness-dichotomy.md),
the complete address graph is balanced away from the nonproper locus:
distinct companion owners force extra addresses, whereas a complete
two-address connected packet forces a common owner and a two-edge cycle.  A
complete forest forces the node into the nonproper locus.  Since a
noninvertible Keller map has
generic degree at least three, the two-edge cycle can only be a connected
subpacket: a finite full node fibre has at least one further address.  This is
the precise replacement for the retracted boundary-forest shortcut.

There is a sharper **CONDITIONAL** completion test using
[THM-3968](../01-canon/theorems/THM-3968-canonical-vector-different-affine-plane-boundary-obstruction.md).
If this lane supplies an actual finite flat cubic normal completion whose
unique boundary prime has
tame inertia two, then a boundary point over the node consumes fibre length
at least two.  Together with the two known reduced affine addresses this
would exceed rank three.  In that completed cubic model the node is therefore
not Jelonek and has exactly one additional affine address.  The missing
hypothesis is precisely that the cusp chart under study is the actual finite
normal completion, not merely a rational/local model.

### D. Oriented conductor Smith compiler

Given a normal cubic cover with labelled boundary primes, compile:

```text
vertices = normalized components,
edges    = conductor branches with orientation and multiplicity,
columns  = principal-divisor responses,
outputs  = Smith form, H1(dual,mu_3), canonical/ramification vectors.       (20)
```

Apply it first to THM-3992's two clutches and THM-3994's two seams.  The
triangle from THM-3990 is the hostile showing why real rank is insufficient.

### E. Iterated repair towers

If the repaired coefficient vanishes on a smaller stratum, restrict the next
operator and recompute its adjoint kernel.  Develop this for weighted smooth
strata, normal crossings, and finite graphs.  The target theorem is an
iterated Fredholm/Farkas certificate with explicit interface conditions, not
an informal “go to fourth order.”

### F. Classification of factorial-law escapes

For rank `n>=3`, classify finite-index periodic fan quotients with stabilizers,
nonunimodular cones, or extra fibres that can reduce the Euler invoice to two.
Use Burnside counts and local intersection cohomology.  A single explicit
escape would create a new complex-geometric frontier; a no-go under mild
singularity bounds would be a substantial theorem independent of the `S6`
claim.

### G. Signed LRC response complex

On the remaining seventeen types, build the exact matrix of *exposed* events
of every `G_u`, with owner, orientation, period, and masking.  Test its strict
Farkas dual first; only surviving rows merit the cubic response `(18)`.  AP13,
Goddyn--Wong, and a `k=0` row are mandatory hostile controls.

### H. Degree-six and unit-ideal positive search

The degree-five common-zero grammar is closed by THM-3890.  Search the
smallest degree-six multibranch infinity packet and the unit-ideal
nonmonogenic packet by the same order:

```text
equisingular deformation space
  -> linearized Keller equations
  -> completed address scheme
  -> resolvent units/class group
  -> oriented conductor and canonical vector.                           (21)
```

This is narrower than another raw coefficient census and retains exactly the
sidecars that killed the quintic mirage.

## 11. Connection contracts

| source | target | map | preserved | destroyed | required sidecar | cheapest decisive test |
|---|---|---|---|---|---|---|
| Brendle residual torus | general repair problem | conformal tensor to Laplacian response | quotient class and sign | curvature formula | completed components, uniform remainder | exact `ker A*` / Farkas test |
| `S6` conductor quotient | JC boundary | normalization and branch gluing to incidence matrix | integral cokernel | complex integrability, Keller equation | branch labels, orientations, units | Smith plus hostile reversed edge |
| THM-3991 fan | higher-rank cusp | simplex orbit count | Euler budget | attaching maps, analytic compactness | stabilizers, cone indices, extra fibres | Burnside Euler census |
| THM-3992 node | cubic-cover boundary | source addresses to dual edges | local clutch type | full fibre and ownership | normalized factors, Jelonek flag | census fibre or prove node nonproper |
| THM-3994 resultant | any eliminant collision | scalar root to base scheme | total length | support distribution and graph singularity | completed ideal/Rees algebra | primary decomposition plus Jacobian |
| LRC tariff equality | repair matrix | lawful owner event to response column | signed left-kernel obstruction | chronology if owners are dropped | exposure, masking, period | exact event matrix on 17 types |
| `3,4,infinity` monodromy | rational elliptic model | local exponents to Kodaira fibres | integral monodromy conjugacy | global toric quotient | lattice marking, cusp width | Weierstrass/Picard--Fuchs rebuild |

## 12. Rejected bridges and stopping reasons

- There is no direct Riemannian-curvature theorem about integrable complex
  structures, and no reverse implication.
- Positivity order has no analogue in the complex Jacobian coefficient ring;
  only the quotient/cokernel mechanism transfers.
- One projective point is not one normalization place.
- A repeated resultant is not a nonreduced point until the common-zero scheme
  is computed.
- A raw kernel in gluing coordinates is not a moduli space before presentation
  gauge is removed.
- A local `Z/2` or `Z/3` class is not a global boundary obstruction before
  the labelled incidence map and its consumer are identified.
- Endpoint-address feasibility is not an LRC body because owner masking and
  chronology have been discarded.
- A finite jet solution is not algebraization or termination.
- Tournament analysis is not used here: no intrinsic complete pairwise
  orientation is present.  Components, branches, events, and proof obligations
  are the faithful vertices.

The session's strongest change is therefore not a headline claim.  It is a
reduction in the dimension of the unknowns.  The quintic common-zero route is
closed in its normal scope; the `2:3` cusp-log cell is forced to a centered
node; the double-collision seams are scheme-theoretically separated; and the
remaining work is now concentrated in three typed objects: a global
quadratic identity, an analytic gluing quotient, and an oriented conductor or
owner-event incidence matrix.
