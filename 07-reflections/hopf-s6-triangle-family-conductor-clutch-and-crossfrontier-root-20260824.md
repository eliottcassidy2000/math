# Hopf/S6 triangle family: conductor, clutch response, and cross-frontier tasks

**Long-session synthesis -- 2026-08-24.**  This reflection treats the
108-page manuscript at [alpo.ge/s6.pdf](https://alpo.ge/s6.pdf) as a
**PREPRINT CLAIM / UNDER AUDIT**.  The durable truth sources produced by the
session are the
[source and referee ledger](../05-knowledge/reference/CORE-PAPERS-HOPF-S6-2026-08-24.md),
[exact matrix audit](../04-computation/hopf_s6_triangle_monodromy_snf_audit_20260824.py),
its [frozen output](../05-knowledge/results/hopf_s6_triangle_monodromy_snf_audit_20260824.out),
and the unconditional local
[THM-3955](../01-canon/theorems/THM-3955-node-cotangent-normalization-kernel-and-conductor-torsion.md)
and [THM-3957](../01-canon/theorems/THM-3957-triple-normal-crossing-cotangent-conductor-kernel-and-normalization-cokernel.md).

No sentence below promotes the manuscript's claimed compact complex
threefold, its topology, or the Hopf problem to proved status.

The session's compact outcome is:

```text
native gluing operation
  -> labelled integral presentation
  -> Smith response and ordinal
  -> kernel/fibre sidecar,

nonnormal object
  -> normalization plus conductor incidence
  -> restriction/extension ledger
  -> consumer-specific obstruction.                         (1)
```

These are two compatible compilers.  The first tells us which integer a gluing
produces; the second tells us what normalization erased before that integer was
computed.

## 1. Inheritance pass

### Closest proved mechanisms

1. [THM-3756](../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md)
   is the clean ordinal precedent: the odd-square shell can be called `n`, but
   a second ordinal is needed to recover a Pythagorean node inside the shell.
2. [THM-773](../01-canon/theorems/THM-773-prime-seven-sheet-monodromy-and-tournament-fibre.md),
   [THM-778](../01-canon/theorems/THM-778-centered-christoffel-endpoint-skew-product.md),
   and [THM-794](../01-canon/theorems/THM-794-unbounded-full-active-prime-seven-packets.md)
   already separate a coarse tournament state from sheet translations,
   ownership, and legal continuation.
3. [THM-2047](../01-canon/theorems/THM-2047-phase-height-toric-arrangement-for-lrc.md)
   is the strongest existing LRC base carrier because it retains cell, phase,
   height, sign, and wall labels.
4. [THM-3944](../01-canon/theorems/THM-3944-repeated-factor-double-torus-one-place-square-conductor-collapse.md)
   and MISTAKE-466 distinguish the original nonnormal regular locus from the
   conductor complement and the full normalization.
5. [THM-3950](../01-canon/theorems/THM-3950-a1-internal-split-denominator-debt-and-equianharmonic-shadow.md)
   and [THM-3951](../01-canon/theorems/THM-3951-affine-plane-boundary-incidence-forest-and-equianharmonic-survivor-nonentry.md)
   show why a normalized surface still needs a boundary-incidence graph before
   one can recognize an affine-plane chart.
6. [THM-3166](../01-canon/theorems/THM-3166-falling-factorial-order-join-path-colour-transform.md)
   supplies the operation-response discipline: compile a native operation
   before comparing its numerical shadow.

### Canonical hostile examples

```text
ordinal collision:  Pythagorean triples (5,12,13) and (15,8,17)
                    share the odd-square shell B+C=25;

clutch collision:   for (a,b)=(3,4), ell and ell+(1,3,0)
                    have the same p;

normalization loss: B=k[x,y]/(xy), where x dy is ambient and nonconormal
                    but maps to conductor torsion and dies after torsion-free
                    normalization;

LRC state collision: THM-773 has rows with the same coarse tournament data
                    but different legal next events;

boundary hostile:   THM-3951 forbids two distinct boundary primes of an A2
                    chart from meeting twice, although each normalized curve
                    may separately be smooth.                              (2)
```

These hostiles all have the same grammatical shape—one scalar or normalized
object survives while a future operation changes—but the missing coordinates
are different.  They cannot be merged without typing the consumer.

### Corrected near misses

- MISTAKE-233 blocks any move from the words "triangle group" or "cusp" to an
  LRC modular argument.  A map must preserve loneliness, not vocabulary.
- MISTAKE-466 blocks replacing the original regular locus by the full
  normalization.  THM-3955 shows the same operation can erase a node
  cotangent class; THM-3957 adds conductor-axis torsion and a normalization
  cokernel at a triple crossing; THM-3944 involves a Kummer class.
- The CDP conflict cannot be decided by the local node alone.  The node proves
  a sheaf-level loss mechanism; globalization and all hypotheses remain open.
- Two homology computations in one manuscript are not automatically
  independent when both consume the same cusp quotient and attaching maps.

### Least-used relevant sidecar

The conductor **incidence multigraph**, not just the conductor ideal, is the
least-used sidecar.  THM-3951 uses an incidence forest to close a Jacobian
survivor.  The S6 manuscript's cusp fibre instead identifies opposite sides of
a hexagon, leaving three double curves and two triple points.  Normalization
alone erases both the self-identification and which branches meet at a triple
point.  Any cup-product, nearby-cycle, or differential calculation which uses
only the normalized `dP6` is therefore under-typed.

## 2. Anchor / Niche / Wildcard portfolio

| lane | question | operation | result |
|---|---|---|---|
| **Anchor -- Hopf/S6 audit** | do the manuscript's displayed algebra and topology maps pass hostile checks? | exact transcribe, exterior powers, Smith forms, two presentation controls | **FINITE-EXACT PASS** for displayed matrices; global analytic/topological gates remain open |
| **Niche -- conductor algebra** | can the claimed failure of a normalization reduction be isolated without trusting the manuscript? | compute `Omega` at `xy=0` and `xyz=0` before and after torsion-free normalization | **THM-3955/3957 PROVED**, characteristic-free |
| **Wildcard -- cross-frontier compiler** | can monodromy/clutch ideas generate lawful LRC, JC, and ordinal tasks? | write source/target/map/predicate/loss/sidecar/test contracts | four precise task programs; false Mahler/ABC/IUT bridges rejected |

The niche produced the unconditional theorem.  The anchor produced a sharp
audit boundary.  The wildcard turned the incoming construction into new work
without making the construction a theorem dependency.

## 3. Live concept board

The board stabilized at seven concepts.

| concept | representation | native operation | invariant/predicate | missing coordinate |
|---|---|---|---|---|
| triangle lattice family | `Delta(a,b,infinity)->SL_4(Z)` | monodromy composition | finite orders plus square-zero cusp | period positivity and analytic quotient |
| clutch word | labelled translations at three boundary components | regluing | presentation matrix | analytic isotopy and branch choice |
| Smith response | invariant factors of an integer matrix | cokernel | rank and torsion | geometry, labels, metric semantics |
| ordinal evaluator | `p=ab ell0-b ell1-a ell2` | unit increment in a gluing coordinate | selected group order | kernel representative and sign |
| conductor ledger | original / complement / normalization | restrict and extend | exact obstruction locus | incidence, torsion, consumer |
| affine local system | chamber graph plus sheet translations | wall crossing and holonomy | normalized return | owner, metric lift, legal next event |
| boundary incidence | labelled multigraph/dual complex | resolution and contraction | cycles/forest conditions | ambient surface type |

Every meaningful pull was compared to every board coordinate:

- the clutch evaluator is a Smith response, not the analytic family;
- the conductor ledger can change which generators legitimately enter a
  presentation;
- an affine local system is useful only on a base which retains the predicate;
- a boundary multigraph may obstruct an `A2` chart but is allowed to have
  cycles on a different ambient surface;
- the ordinal is lawful for scheduling even when its fibre is large;
- exact monodromy data says nothing about period positivity or properness;
- a tournament is usable only where an intrinsic pairwise orientation already
  exists.

## 4. Exact local result: where normalization loses the cotangent class

Let

```text
A=k[x,y],        B=A/(xy),        Btilde=k[x] direct_sum k[y],
E=Omega_(A/k) tensor_A B=B dx direct_sum B dy.             (3)
```

Let `T` be torsion in `Omega_(B/k)`, meaning elements annihilated by a
nonzerodivisor, and put `K=ker(E -> Omega_(B/k)/T)`.  THM-3955 proves

```text
T=k*tau,             tau=x dy=-y dx,
Omega_(B/k)/T ~= Omega_(Btilde/k),
K=(y)dx direct_sum (x)dy
 =B*d(xy)+B*(x dy),
0 -> B*d(xy) -> K -> T -> 0,
K/B*d(xy) ~= B/(x,y) ~= k.                                (4)
```

The proof survives characteristic two: then `x dy=y dx`, and the conormal map
`B -> E`, `1 |-> y dx+x dy`, remains injective.

The subtle point is the location of the loss.  The torsion-free normalization
map in (4) is an isomorphism.  The extra one-dimensional kernel occurs before
it, in the ambient map `E -> Omega_B/T`.  Thus the invalid shortcut is

```text
ambient restriction
  = conormal contribution + torsion-free fibre contribution,               (5)
```

because the preimage of conductor torsion is a third term.

This was the first unconditional local advance of the session.  It proves the
node mechanism asserted in the S6 manuscript's response to CDP20; the triple
sequel below strengthens the loss ledger.  Neither result globalizes a local
form or verifies the claimed threefold.

## 4A. The triple point changes the normalization map itself

The node theorem cannot simply be copied to the manuscript's two claimed
triple points.  For

```text
B=k[x,y,z]/(xyz),
r=yz dx+xz dy+xy dz,
K=(yz)dx direct_sum (xz)dy direct_sum (xy)dz,             (4a)
```

THM-3957 proves

```text
0 -> B*r -> K -> T -> 0,
T ~= Btilde/B,                    Ann(T)=(xy,xz,yz),

0 -> T -> Omega_B -> Omega_Btilde
  -> k[x]dx direct_sum k[y]dy direct_sum k[z]dz -> 0.     (4b)
```

Thus there are two distinct losses:

1. the ambient-to-torsion-free kernel `K`, generated by the three branch
   forms modulo one conormal diagonal;
2. the failure of torsion-free differentials to fill all normalized
   differentials, measured by coefficient mismatches along the three axes.

The torsion has rank one along each double axis and fibre dimension two at the
triple point; the normalization cokernel has triple-point fibre dimension
three.  A normalized `dx` placed only on the plane `y=0` is the minimal
hostile: it disagrees with the zero branch on their common `x`-axis and cannot
descend.

This is a genuine extension of the session, not a manuscript dependency.  It
holds over every field and in the completed or convergent local models.  It
also changes the global audit task: a conductor calculation on `W` must retain
both the torsion sheaf and the branch-mismatch cokernel.

## 5. Exact displayed monodromy and clutch response

The finite audit verifies the manuscript's displayed matrices exactly:

```text
T1^3=I,        T2^4=I,        T2^2!=I,
T0=(T1 T2)^(-1)=I+N,          N^2=0, rank N=2,
A1=(T1^(-1))^t, A2=(T2^(-1))^t, M0=(T0^(-1))^t,
A1 A2 M0=I.                                             (6)
```

The joint one-dimensional lattice coinvariant is exact:

```text
SNF([A1-I | A2-I])=(1,1,1,0).                          (7)
```

On exterior degrees `q=0,1,2,3,4`, the cusp invariant ranks are

```text
(1,2,4,2,1),                                            (8)
```

with saturated images for the displayed matrices.  The degree-two joint
coinvariant has Smith diagonal `(1,1,1,1,1,0)`, and the manuscript's quadratic
class has coinvariant value 12.  The cusp conductor/pushout matrix and the
special-surface map have respective Smith data

```text
(1,1,1,1,0,0),       (1,1,1,0).                        (9)
```

These reproduce the integral arithmetic the later topology wants.  They do
not prove that the analytic attaching maps induce those matrices.

There is now a second, independent conductor calculation conditional only on
the stated oriented quotient.  With `F=dP6`, hexagonal conductor `partial`,
and quotient double locus `D`,

```text
chi(F)=6,            chi(partial)=6,          chi(D)=2,
chi(W)=6+2-6=2,
H_0..H_4(W)=(Z,Z^2,Z^4,Z^2,Z),
pi1(W)=Z^2.                                               (9a)
```

The degree-two pushout map has rank four and Smith diagonal
`(1,1,1,1,0,0)`.  The hexagon loop maps to a commutator in the free group on
the two double-locus loops, which gives the nonabelian statement in (9a).
If just one paired branch is orientation-reversed, the matrix instead has
rank five and Smith diagonal `(1,1,1,1,2,0)`.  Euler characteristic stays two,
but `Z/2` appears.  Thus "opposite sides" and incidence counts are not enough;
the equal holomorphic branch degrees are load-bearing.

This verifies the finite topology of the stated quotient, not its realization
by the infinite fan, the cusp-neighborhood retraction, or any global `S6`
claim.

For general coprime boundary orders `(a,b)` the displayed presentation grammar
becomes

```text
R=[[a,-ell1],[b,ell2-b ell0]],
p=ab ell0-b ell1-a ell2,
det R=-p.                                                (10)
```

Under the admissibility conditions `gcd(ell1,a)=gcd(ell2,b)=1`, the gcd of all
entries of `R` is one.  Hence for `p!=0`

```text
coker R ~= Z/|p|.                                       (11)
```

This general algebraic statement is independent of whether any corresponding
complex torus family exists.

## 6. The ordinal lesson: call the response n, retain its fibre

The odd-square Pythagorean convention says that after selecting the `n`th odd
square, it is often cleaner to call that square simply `n`.  The same move is
lawful here: after selecting a positive clutch response `|p|`, call it task
`n` and work inside its fibre.

What is not lawful is replacing the gluing object by `n`.  For `(a,b)=(3,4)`,

```text
phi(ell0,ell1,ell2)=12ell0-4ell1-3ell2,
(1,3,0) in ker phi.                                     (12)
```

Thus a lossless scheduling address can be

```text
(rank of |p|, shortlex rank of a chosen kernel representative,
 admissibility witness, ordered gluing word).            (13)
```

The evaluator and its native unit operations remain useful:

```text
ell0 += 1  changes p by ab,
ell1 += 1  changes p by -b,
ell2 += 1  changes p by -a.                             (14)
```

This is a procedural task generator.  It creates neighboring clutch problems
with a known response delta while retaining the fibre data needed to ask
whether their analytic realizations, homology, or moduli differ.

The triangular-number lesson also survives in the background: a finite
difference such as `T(Z+2)-T(Z-2)=2+4Z` is valuable because it is the response
of a named native shift.  Equation (14) is the clutch analogue.  Neither
formula alone reconstructs the state on which the shift acted.

## 7. Connection contract A: conductor three-site ledger for planar JC

**Source.**  The exact node computation (4), together with the manuscript's
claimed nonnormal cusp fibre.

**Target.**  MISTAKE-466, THM-3944, and the current nonnormal boundary orders
in the planar Jacobian program.

**Map.**  For a nonnormal surface `Y`, replace the binary regular/singular flag
by

```text
Y_reg(original)
  -> Ytilde minus nu^(-1)(conductor)
  -> Ytilde,                                             (15)
```

and attach the conductor incidence multigraph.

**Preserved predicate.**  A class or section on the exact site where it is
defined, together with its restriction and extension behavior.

**Destroyed information.**  Conductor-supported torsion, branch labels,
self-identification, and finite boundary characters when one jumps directly
to the full normalization.

**Needed sidecar.**  Conductor ideal, preimage branches, incidence graph,
valuation/torsion vector, and the target consumer.

**Cheapest decisive tests.**

1. the node `xy=0`, where the missing cotangent kernel is one-dimensional;
2. THM-3944's nonzero `(2,1)` Kummer class on the original regular locus which
   does not extend over `A2`;
3. for each current JC boundary order, compute all three sites and ask which
   Cardano/Kummer characters survive each arrow.

**Survivor.**  Nonnormality is not itself an exclusion criterion.  A candidate
is excluded only when the actual predicate fails on the actual site consumed
by a Keller chart.

## 8. Incoming signal: THM-3950--3958 changes the sidecar

During this session, incoming work promoted THM-3950 and THM-3951 to
**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED**.  THM-3950 shows
that a normal quadratic surface and two Cardano classes can still carry an
equianharmonic residual.  THM-3951 closes that survivor only after proving
that boundary incidence for an `A2` open surface is a forest and that two
distinct boundary primes cannot meet at two distinct smooth points.

This changes the conductor board in two ways:

1. normalization plus local characters is still insufficient;
2. the dual/incidence multigraph is an operation-bearing sidecar, not a
   decorative picture.

The direct theorem transfer to the Hopf fibre fails.  THM-3951's forest comes
from resolving a completion of `A2`; the S6 manuscript explicitly wants a
quotient with opposite-side identifications and triple points.  Its incidence
multigraph is allowed to cycle because its open surface is not asserted to be
`A2`.  The surprising connection is therefore a validity rule:

```text
normalize, retain boundary incidence, then apply the ambient-specific
recognition theorem.                                     (16)
```

This rule is strong enough to guide both frontiers and weak enough not to
manufacture a false equivalence.

The next incoming promotions sharpen the same rule in three orthogonal ways:

- [THM-3952](../01-canon/theorems/THM-3952-minimal-mobius-internal-split-carriers-are-four-critical-colors-and-nonentry.md)
  classifies exactly four critical infinity colors and then uses the boundary
  incidence obstruction rather than a normalization scalar;
- [THM-3953](../01-canon/theorems/THM-3953-rationally-split-hidden-cubic-ramification-triangle-nonentry.md)
  makes three distinct split ramification primes form a forbidden boundary
  triangle, including incidences at singular points on an already-normal
  surface;
- [THM-3954](../01-canon/theorems/THM-3954-extra-common-debt-creates-A-singularities-and-nonunibranch-residual-boundary.md)
  shows that extra common debt gives a normal `A_(3m-1)` local surface while
  one residual prime still has two normalization addresses at the same point;
- [THM-3956](../01-canon/theorems/THM-3956-split-hidden-cubic-integrality-and-repeated-root-trichotomy.md)
  closes the `k(t)`-split hidden-cubic lane: monicity removes rational
  denominators, and the repeated-root branches are eliminated by reducibility,
  ramification-unit, or class-group gates;
- [THM-3958](../01-canon/theorems/THM-3958-one-hidden-root-principal-different-and-pure-power-boundary.md)
  closes the exactly-one-hidden-root lane by showing that the reduced relative
  different is the sum of two required boundary primes, hence becomes a
  forbidden nonconstant unit on any same-field affine-plane open.

All five are **PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED** in
their stated scopes.  Together they refute the coarse heuristic that normality
itself is the decisive coordinate.  The decisive coordinate is the labelled
incidence pattern and ramification predicate consumed by the ambient
recognition theorem. THM-3958/3960 close the one-root and natural
one-parameter monogenic lanes; proved THM-3959 closes the centered degree-five
color rows; and proved THM-3961 classifies normality for irreducible
arbitrary-`q` monogenic cubics. Proved THM-3962 closes coefficient-constant
cylinders and proved THM-3963 closes the moving scalar `P^2` debt. Proved
THM-3964 closes the displayed scalar-coefficient graph-root family; THM-3965
rules out a one-place discriminant in its exact constant-`(g,h)` unit-ideal
family; THM-3966/3968 give class/Euler and canonical/different boundary
invoices. THM-3967 closes every irreducible `deg_P(q)<=2` row, and THM-3969
closes the collision-free `Xi in k^*` affine-`P` graph packet. THM-3970/3971
remain RESERVED only. None supplies a Hopf dependency here.

## 9. Connection contract B: affine local systems over the LRC carrier

**Source.**  The manuscript's sequence "monodromy, translation twist,
coinvariants" and the exact one-dimensional displayed coinvariant (7).

**Target.**  The prime-seven sheet data of THM-773/778/794 over THM-2047's
labelled phase-height cell complex.

**Map.**  Over an endpoint chamber use the affine sheet state

```text
k in F_7^r,
k |-> k-sum_(a in crossed block A) w_a^(-1)e_a,          (17)
```

then quotient only the diagonal translation `<1>`.

**Preserved predicate.**  Exact sheet coverage and normalized holonomy;
diagonal translation permutes all seven sheets.

**Destroyed information.**  Marked sheet after quotient, and—if the base is
also collapsed—owner, metric chronology, global carry, and legal next event.

**Needed sidecar.**  Owner-to-sheet assignment, inverse steps, endpoint order,
phase/height, global carry, and metric lift.

**Cheapest decisive tests.**

1. replay the two THM-773 states sharing coarse tournament data but having
   different next events;
2. compute normalized holonomy on the unbounded full-active packets of
   THM-794;
3. compare a cellular boundary computation with the affine-local-system
   coinvariant computation, mirroring the manuscript's two topology routes.

**Tournament boundary.**  The tournament is the coarse base because LRC has
an intrinsic pairwise endpoint relation.  The owner-labelled sheet assignment
and continuation holonomy form its stalk.  Hopf lattice vectors themselves
have no intrinsic pairwise orientation, so there is no Hopf tournament.

No statement here improves the `1/14` loneliness threshold.  It produces a
better-typed finite carrier on which an LRC proof obligation can be computed.

## 10. Connection contract C: the Smith response compiler

**Source.**  Presentation (10) and the exact agreement of its determinant and
cokernel.

**Target.**  Open Smith programs HYP-1787, HYP-6835, and the operation-response
discipline behind THM-3166.

**Map.**  Compile each native gluing, continuation, or extension operation into
a labelled integral relation/chain matrix.  Take Smith form only afterward.

**Preserved predicate.**  Integral rank, cokernel, and torsion of exactly the
modeled operation.

**Destroyed information.**  Metric loneliness, positivity, ownership,
orientation labels after forgetting the basis, and response to an operation
not present in the matrix.

**Needed sidecar.**  Row/column generator labels, orientation convention,
operation word, and target evaluator.

**Cheapest decisive test.**  Construct the LRC continuation boundary twice:
cellularly and from local-system coinvariants.  Require Smith forms **and the
labelled comparison map** to agree.  Equal rational rank or equal mod-two rank
is an intentionally hostile insufficient control.

## 11. Connection contract D: a general triangle-family task factory

The pair `(3,4)` should not be fetishized.  Equation (10) gives a purely
arithmetic task family for coprime `(a,b)`.

### Objects

```text
coprime orders (a,b),
finite-order T1,T2 in SL_4(Z),
square-zero unipotent T0=(T1T2)^(-1),
admissible local translations,
period/fan/log-filling passports,
clutch response p.                                       (18)
```

### Representations

- integral matrices and exterior powers;
- affine local systems and coinvariants;
- fans, conductor dual complexes, and boundary labels;
- relation matrices and Smith forms;
- ordinal plus kernel representative.

### Operations

- change one clutch coordinate using (14);
- dualize or take exterior powers;
- vary `(a,b)` through coprime pairs;
- normalize while retaining the conductor incidence;
- change scale from local node to global fibre;
- compare two independent chain models.

### Cheapest staged search

1. enumerate small finite-order pairs in `SL_4(Z)` with square-zero cusp;
2. reject pairs whose joint coinvariant has the wrong rank or torsion;
3. compute exterior invariants and candidate Euler contributions;
4. enumerate admissible twists by shortlex ordinal in `|p|` fibres;
5. only then attempt period maps, positivity, fan quotients, and log fillings;
6. require collapse and nearby-cycle computations to agree through a labelled
   comparison map.

This factory is useful even if the S6 manuscript fails: it separates finite
searchable algebra from the genuinely analytic realization problem.

## 12. Rejected transfers and stopping reasons

### Mahler's `3/2` problem

The numbers `3` and `4` and a torus orbit do not map to the rational-base
orbit `{(3/2)^n x}`.  The repo's Mahler work consumes digit carries,
half-circle tails, and stabilization in ordinary integers.  Triangle
monodromy preserves none of these.  **STOP: no trajectory/predicate map.**

### ABC conjecture

The clutch equation for `p` is an integer linear evaluator, not a coprime
additive triple with radical and height.  Rewriting three signed summands as
`A+B=C` would be cosmetic because admissibility and the target are different.
**STOP: no radical/height consumer.**

### Inter-universal Teichmuller theory

Typed charts, transition words, and monodromy are a general geometric grammar,
not an instantiation of IUT's alien copies or disputed height comparison.  The
manuscript supplies no Frobenioid/Hodge-theater functor and no weighted-height
membership step.  **STOP: vocabulary only.**

### Direct planar Jacobian transfer

A compact complex threefold fibration is not a polynomial Keller map of the
affine plane.  Only the normalization/conductor operation and its incidence
sidecar transfer.  **STOP: ambient and target predicate differ.**

### Unstable homotopy

Knowing homotopy groups of spheres does not test integrability of an almost
complex structure.  The manuscript's final topology consumes ordinary
fundamental group, integral homology, smoothness, and recognition, not a novel
unstable-homotopy invariant.  **STOP: no integrability detector.**

## 13. Generated next tasks

The session leaves eight precise tasks, ordered by expected information gain.

1. **Cusp analytic rebuild.**  Reprove free/proper action and compactness for
   the infinite fan quotient with explicit fundamental domains and estimates.
2. **Conductor realization bridge.**  The oriented finite chain model is now
   exact; prove that the infinite-fan analytic quotient realizes its branch
   degrees, deformation retraction, specialization maps, and attaching word.
3. **Global cotangent audit.**  Starting from THM-3955 and THM-3957, compute
   both the torsion sheaf and branch-mismatch cokernel on the claimed `W`, then
   determine whether any local generator globalizes with the required twist.
   Keep local nonzero, global nonzero, descent, and ambient extension separate.
4. **Half-weight period audit.**  Check automorphy, multiplier, square-root
   branch, and the exact `O(-1)` torsor class independently.
5. **LRC affine-local-system prototype.**  Put THM-773 transitions over a
   small labelled THM-2047 chamber complex; retain owner and metric lift;
   compute holonomy and strongly connected components.
6. **Two-path Smith comparison.**  Compile that LRC prototype both cellularly
   and by coinvariants, then compare labelled Smith presentations.
7. **JC three-site census.**  For each surviving boundary grammar after
   THM-3942, 3946, 3949, incoming 3950--3954, 3956, and 3958, record original regular locus,
   conductor complement, normalization, incidence graph, and every surviving
   Kummer/Cardano character.
8. **Coprime triangle scout.**  Enumerate small `(a,b)` matrix pairs satisfying
   the finite algebraic gates, organize twists by `(rank |p|, kernel rank)`,
   and record the first obstruction before any analytic construction.

Every task has a positive control and a hostile control.  None begins by
assuming the manuscript's main theorem.

## 14. Final truth boundary

### PROVED

- THM-3955's node and THM-3957's triple-crossing cotangent sequences;
- the general two-by-two Smith conclusion (11) under its stated coprimality
  and admissibility assumptions;
- the repo theorems cited with proved status, including planar-Jacobian
  THM-3950--3954, THM-3956, and THM-3958--3963.

### FINITE-EXACT

- all displayed integer-matrix, exterior-power, and presentation checks in the
  frozen audit;
- conditional on the stated oriented conductor quotient, `chi(W)=2`,
  `pi1(W)=Z^2`, and free homology ranks `(1,2,4,2,1)`, with a one-reversed-
  branch `Z/2` hostile;
- the 54,320-tuple hostile sweep in its stated finite universe.

### READ / NO IMMEDIATE CONTRADICTION

- the period, filling, and gluing narrative on a first full pass;
- the manuscript's local explanation of the CDP20 normalization loss.

### OPEN / UNDER AUDIT

- analytic period and positivity details;
- infinite fan quotient and cusp degeneration;
- global conductor section and CDP20 conflict;
- integral nearby cycles, attaching maps, and both homology routes;
- compact smooth threefold, diffeomorphism to `S6`, and a complex structure on
  `S6`;
- every LRC(14), planar JC, Mahler `3/2`, ABC, or IUT frontier consequence.

The main advance is therefore not a claimed solution by proxy.  It is an exact
local theorem, a reproducible finite certificate, and a set of better-typed
research programs generated from the manuscript's strongest mechanisms.
