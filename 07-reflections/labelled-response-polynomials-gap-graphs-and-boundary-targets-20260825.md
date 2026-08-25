# Labelled response polynomials, gap graphs, and boundary targets

**Research synthesis, 2026-08-25.**  The proved advances are
[THM-4097](../01-canon/theorems/THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension.md),
[THM-4099](../01-canon/theorems/THM-4099-squarefree-gap-transfer-and-mixed-insertion-boundary.md),
[THM-4102](../01-canon/theorems/THM-4102-selected-order-ten-strong-ear-solid-interval.md),
[THM-4103](../01-canon/theorems/THM-4103-jc23-theta-boundary-ramification-and-degree-response.md),
and
[THM-4109](../01-canon/theorems/THM-4109-ap7-weight-seven-gap-atlas-and-sharp-pair-overlap.md).
The multiaffine response object and unequal-arc trapezoid below are elementary
exact deductions; the proposed transfers and next attacks are **OPEN**.
LRC(14), JC(2), the global tournament `H`-spectrum conjecture, and the
four-dimensional Euclidean Kakeya conjecture remain **OPEN**.

## Inheritance pass

- **Closest proved mechanisms.**  THM-4092 supplies the parity-weighted comb
  compiler; THM-4094 supplies exact insertion incidence; THM-4053 supplies
  the last live `(2,3)` Jacobian boundary.
- **Canonical hostiles.**  AP7 bank `{8,9,11,13}` covers the whole safe
  interval; the transitive triple and directed triangle have identical
  proper insertion faces but different mixed coefficients; deleting the
  mandatory Jacobian vertex `(4,2)` preserves genus and total ramification
  while changing the ramification packet.
- **Corrected near miss.**  MISTAKE-499 separates support of tournament
  correction terms from failure of an equality: nonzero terms can cancel.
  The analogous warning in the other lanes is that a zero one-pair margin or
  a valid scalar degree does not decide the full object.
- **Least-used sidecars.**  For LRC it is labelled pair-overlap geometry; for
  tournaments it is the full gap/block assignment; for JC it is the target
  image and Galois orbit of each boundary point.

## Live concept board

| lane | object retained | operation | invariant | decisive loss under scalarization |
|---|---|---|---|---|
| LRC anchor | intersection polynomial of danger combs | add one overlap coefficient, then all pairs | literal intersection measures | which speeds overlap on which component |
| tournament niche | squarefree gap polynomial | insert labelled blocks into base-word gaps | every face Hamiltonian count | mixed placement and repaired cuts |
| JC wildcard | labelled infinity-puncture packet | map punctures to finite target or infinity | ramification index and Galois block | target value, collision, leading coefficient |
| arithmetic firewall | rational/algebraic/transcendental approximation packet | order or valuation comparison | chosen sign/order | height, determinant, decay, product formula |
| Kakeya/sixty firewall | phase-coloured line packet | shift phase, place and rescale tubes | address and direction | basepoint, shading, transversality, scale ancestry |

The board changed the session's main question.  Instead of asking whether a
scalar statistic is large enough, ask for the smallest labelled response
algebra in which the next operation is still compositional.

## Anchor: AP7 becomes a weighted gap graph

Let `V` be a finite outlier bank and retain THM-4109's interval `J`.  Define
the multiaffine intersection polynomial

```text
I_V(Y)=int_J product_(v in V)(1+Y_v 1_(U_v)(theta)) dtheta.       (1)
```

Literal expansion gives

```text
[product_(v in S)Y_v] I_V
   =mu(J intersect intersection_(v in S)U_v)                    (2)
```

for every `S subseteq V`.  Thus the exact union measure is the alternating
evaluation of all nonempty coefficients, while the first-moment compiler
keeps only degree one.  THM-4109 keeps one labelled degree-two coefficient

```text
M_g(u)=[Y_u Y_(u+g)]I_V.                                      (3)
```

The full second-factorial-moment certificate is the still stronger exact
inequality

```text
mu(J intersect union_v U_v)
 <=sum_v [Y_v]I_V-(1/2)sum_(v<w)[Y_vY_w]I_V.                   (4)
```

Equation `(4)` follows pointwise from
`1_(F>0)<=F-(1/2)binom(F,2)` for `0<=F<=4`.  It explains why `197` and `232`
are sharp for the selected-pair certificate but need not be sharp after the
next lawful operation.  This is the LRC analogue of THM-4099's squarefree
gap polynomial: singleton faces do not determine a mixed face.

The natural carrier is not a tournament.  Its vertices are speeds, its edge
weight is

```text
c(v,w)=mu(J intersect U_v intersect U_w),                       (5)
```

and ties are genuine.  For two odd vertices of nearby speeds, THM-4109 proves
that the limiting edge weight is zero exactly at gaps `2,6,10`.  Those three
gaps form a sparse zero-edge graph.  A hostile three-odd configuration tries
to place its two unused edges in that graph while the designated gap is `4`
or `8`.  This is why the configurations

```text
gap 4:  {u-2,u,u+4}, {u,u+2,u+4}, {u,u+4,u+6}, ...
gap 8:  {u-2,u,u+8}, {u,u+2,u+8}, {u,u+6,u+8}, ...             (6)
```

are the correct first hostiles.  Orienting these edges would discard the
quantity that `(4)` needs.

### The even--odd edge is an unequal-arc trapezoid

For even `e` and odd `o`, put

```text
x=e theta,       s=(2o-e)theta,
A=(-1/14,1/14),  B=(-1/7,1/7) in R/Z.                          (7)
```

The simultaneous danger condition is `x in A` and `x+s in B`.  Its frozen
phase length is therefore

```text
K(s)=mu(A intersect (B-s))
    = 1/7                     if ||s||<=1/14,
      3/14-||s||              if 1/14<=||s||<=3/14,
      0                        if 3/14<=||s||<=1/2.             (8)
```

This is an exact trapezoid, the unequal-width companion to THM-4109's
triangle `H`.  The next uniform proof can therefore be organized around the
three-star response

```text
K((2o_1-e)theta)+K((2o_2-e)theta)+K((2o_3-e)theta),            (9)
```

together with endpoint residue laws.  The source is the labelled bank, the
target is the degree-two truncation of `(1)`, the map is coefficient
extraction, and the preserved predicate is the Bonferroni survivor margin.
Scalar total overlap destroys which zero-gap pattern caused the loss.  The
needed sidecar is the three offsets `(2o_i-e)` and the clipped endpoint cell.

The cheapest decisive test is now precise: classify the zero-edge odd
triples in `(6)`, derive an exact residue-ray law for the three trapezoids in
`(9)`, and determine the true full-pair uniform floors.  A bounded scan is
only reconnaissance; positivity outside its declared box cannot be promoted
without the residue-ray tail proof.

## Niche: tournament insertion should stay squarefree for one more step

THM-4099 proves that for a base tournament `B` and labelled inserted set `I`,

```text
[X_S]Z_(W/B)=H(W[B union S]).                                  (10)
```

For two inserted vertices this is a four-state transfer.  The transitive
triple/directed-triangle hostile proves that all proper-face counts and the
first-step count profile can agree while the mixed coefficient differs.
This is exactly the failure mode repaired by `(1)`--`(4)` in the LRC lane.

THM-4097 then certifies every odd `H` from `85` through `2881` at order nine,
and its independently reconstructed atlas shows that neither order-eight
`H=613` class reaches `623` by a nonconstant ear; a parent with `H=99` does.
So scalar parent size is not an ancestry coordinate.  THM-4102 uses one
labelled parent for every known order-nine value and already forces the
global allowed prefix through `14655` at order ten.

The next operation should not be another scalar parent histogram.  Take
three inserted vertices, retain the eight squarefree faces, and use the
defect filtration: a degree-`r` coefficient sees only base words with at most
`r` bad adjacencies.  The concrete target is an eight-state rank-three
transfer that explains which order-nine parents fill the two solid
order-ten intervals and which mixed faces control the first gap `14657`.

Connection contract:

```text
source:       labelled base words and inserted blocks;
target:       eight squarefree face coefficients;
map:          allocate disjoint inserted blocks to exposed gaps;
preserved:    Hamiltonian-path count on every induced inserted face;
destroyed:    chronological ancestry after coefficient summation;
sidecar:      actual gap allocation and defect layer;
cheap test:   compare parents with equal H but different rank-three image. (11)
```

This remains tournament mathematics because orientations are intrinsic arcs.
It supplies no arithmetic hierarchy of rational, irrational and
transcendental numbers: THM-4088/4095 prove that order tournaments are blind
to those types without quantitative height and error sidecars.

## Wildcard: the Jacobian survivor is a finite Nielsen-class problem

THM-4103 replaces an unbounded degree question by the exact source-infinity
packet

```text
{1,2,2,3,3,3,7},                                             (12)
```

which spends all fourteen Riemann--Hurwitz units.  The index-one point is
forced to target infinity, and the two index-two points are one quadratic
Galois orbit.  The only remaining degree/profile rows are

```text
n=7:   (a,b,c,d)=(0,2,1,1),
n=12:  (a,b,c,d)=(0,0,3,0),
n=21:  (a,b,c,d)=(0,0,0,0).                                (13)
```

The scalar generating product

```text
(1+z)(1+z^2)^2(1+z^3)^3(1+z^7)                             (14)
```

is the same kind of shadow as a tournament face histogram or an LRC first
moment.  The Galois-aware product replaces the two `z^2` atoms by one locked
`z^4` block; edge labels split the three index-three atoms.  What is still
missing is not another degree sum but the target value of each labelled
puncture and its leading `(A,C)` coefficient.

There are two orthogonal next attacks.

1. **Local target-response ideal.**  Expand `A,C` on the index-seven and
   three index-three edges, eliminate the permitted Laurent coefficients,
   and test each row in `(13)` for a compatible finite target point.  The
   hostile is the deleted `(4,2)` support: it preserves genus and total
   ramification but changes the packet, so the full Newton support must stay
   in the elimination.
2. **Genus-one Nielsen classes.**  For `n in {7,12,21}`, enumerate transitive
   permutation data with local cycle lengths `(12)` and the elliptic relation
   `[alpha,beta] product_i sigma_i=1`, respecting the Galois-locked pair and
   the finite/infinite partition `(13)`.  Nonexistence would eliminate a
   degree before coefficient algebra; existence supplies a hostile showing
   that monodromy alone is insufficient.

The second proposal is a genuine finite group problem, not Tournament
Analysis: the vertices are sheets and the native operations are permutations,
commutators and branch cycles, not pairwise orientations.

## Why the older themes fit only through sidecars

### The sixty clock and four-dimensional Kakeya

THM-4035's AP tail, Fibonacci modulo ten, and triangular numbers modulo
thirty share a pointed `C_60` address, but they have different evaluators.
The response-polynomial lesson sharpens that firewall: a common phase can
index labelled coefficients; it cannot supply the coefficient values,
basepoints, or compatibility law.  In four-dimensional Kakeya, direction is
only the label.  A lawful analogue of `(1)` would also retain tube placement,
shading, three- and four-way transversality, and ancestry across scales.
Without that multiscale sidecar, a 60-phase direction atlas is a finite
laboratory, not a Euclidean Kakeya estimate.

### Sun's 2-4-6-8 hole

The exact hole at `896315812331399` has complete residue support modulo every
modulus.  This is the additive-basis version of proper-face blindness:
nonempty local projections do not determine a common bounded integral fibre.
The missing coefficient is height-labelled compatibility, not a new modulus.
The role-labelled representation hypergraph has ties and multiplicities, so
it is not naturally a tournament either.

### p-adic zeta and matching logic

The supplied p-adic-zeta manuscript and matching-logic preprint remain under
their recorded source-audit statuses.  Their usable contribution here is a
quantifier discipline: local signs, reachability, or valuation order can be
correct while the global carrier is absent.  THM-4089's next-case negative
margins show that retuning an already scalarized certificate cannot recover a
lost geometric input.

## Ranked mathematical obligations, not ranked analogies

1. Prove the three-trapezoid AP7 residue law and determine full-pair floors
   for gaps `4,8`; test the exact zero-gap triples first.
2. Build the rank-three tournament squarefree transfer and measure the first
   order-ten interval boundary against equal-`H` parent hostiles.
3. Enumerate the three JC Nielsen-class universes, then feed survivors to a
   labelled Laurent target-response elimination.
4. Replace the Kakeya direction clock by a multiscale placement/shading
   tensor and attack it with broad/narrow hostile pairs having the same phase.
5. For Sun and p-adic lanes, attach bounded height or adelic magnitude before
   asking a local packet to imply a global arithmetic conclusion.

The reusable research move is now concrete:

```text
form labelled face/intersection responses
-> compose the next operation there
-> identify the first scalar quotient collision
-> restore only the missing block/orbit/endpoint label
-> run a hostile and a positive control
-> scalarize only after the target predicate is decided.                    (15)
```

No cross-problem implication is claimed.  The shared advance is a smaller
state design for each frontier and a precise test for whether that state is
still compositional.
