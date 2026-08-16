---
id: THM-3479
title: "Literal half-twist relation-current two-transplant certificate"
status: >
  PROVED STRUCTURAL + VERIFIED-EXACT + FINITE-EXACT + INDEPENDENTLY AUDITED.
  One explicit all-91-unit relation is realized by two primitive q=11 owner
  tuples.  U_full has literal masks, positive delayed mass, all 169 coarse
  unrestricted THM-2334 target aggregates nonzero, and an exact guard-packet-
  row-refined 13^3 endpoint bank whose five role values give nonzero bridge,
  two K4 factors, and product in all 72 labelled charts.  U_clock has literal
  masks at the fixed common centre 1/22, the same-clock delayed word, positive
  delayed mass, a two-embedding nonconstant endpoint bank, and one explicit
  THM-2331 atomic address term.  Exact q=27 and q=51 CRT lifts preserve only
  their stated mod-7 and 13-adic decorations, including q=51's affine k mod 3
  character.  The natural C13-equivariant 13-character-fibre-to-edge map is
  obstructed; the 72 role charts are labelled, non-equivariant, and require an
  owner-order orientation gauge.  No grouped relation coefficient, all-91-
  unit aggregate, ancestry/bispectrum, physical current, scalar-row exclusion,
  or LRC(14) conclusion follows.
source: codex/relation-current-bridge/2026-08-15
audit: >
  standalone deterministic exact companion; independent periodic-subtraction
  and boundary-refinement endpoint engines; exact finite-field nonvanishing
  under certified primitive embeddings; normal/-O/stored LF-identical replay;
  dependency and sidecar hash pins; AST no-assert/security gate; immutable
  core and carrier semantic digests; exact 13^3 refined endpoint bank;
  independent graph-factor reconstruction; independent relation, CRT,
  automorphism, quotient/back-map, and common-ancestry scope audit; promotion
  audit accepted 2026-08-16
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2331-two-sided-septimal-address-embedding-in-marked-current
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-3398-general-finite-mode-sheet-cover-cochain
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
  - THM-3429-prime-fibre-activity-descent-for-mixed-order-half-twist-seven-covers
  - THM-3453-global-literal-half-twist-cap-seven-support-classification
  - THM-3461-literal-half-twist-common-centre-lifts-and-q83-rank-nine-boundary
  - THM-3473-three-times-p-eight-owner-private-sheet-partition-and-irredundancy
related:
  - THM-2512-lawful-interaction-cut-bundle-transplant-and-replica-dichotomy
  - HYP-9032-the-transplant-trichotomy-rehoming-the-91-stalk-laws
script: 04-computation/lrc_half_twist_relation_current_bridge_thm3479.py
output: 05-knowledge/results/lrc_half_twist_relation_current_bridge_thm3479.out
script_sha256: ad2a620cdc238f28e3384698b2c612f38cdf2566bd56b76d1cbabcc03107ec0b
output_sha256: e1e4355290c493d8f9518a4bba3825173bb3137c3d50df421816946e6b1f1207
semantic_sha256: 1bf53086d3da347dde443175462be716da6e9dac54c96582718d19ec8fddff21
carrier_semantic_sha256: 3f17c2206feec73da48a989ab2150ceb1c7d1bc275c77291df476d882957581a
hash_basis: LF-normalized bytes
---

# THM-3479 -- literal half-twist relation-current two-transplant certificate

**PROVED STRUCTURAL + VERIFIED-EXACT + FINITE-EXACT + INDEPENDENTLY
AUDITED.**

This theorem records a concrete positive relation-current transplant and its
exact stopping boundary.  The two owner tuples below realize different parts
of the desired conjunction and must not be silently identified.  Promotion
incorporates the exact U_full guard-packet-row-refined endpoint bank and its
independent graph audit, while retaining every physical-current and ancestry
boundary found by the independent end-to-end audit.

## 1. The common relation and target quotient

Use relation-coordinate order

```text
(c1,c2,c3,H,q1,q2,q3,q4,q5).                         (1)
```

The explicit relation is

```text
a=(-27,-27,-27,20110798,-41,-27,-27,-27,38).         (2)
```

Every coordinate of `a` is a unit modulo `91`, and

```text
a mod 91=(64,64,64,71,50,64,64,64,38),
a mod 13=(12,12,12,6,11,12,12,12,12).                (3)
```

For each of the two tuples `U` below, the exact integer equation is

```text
a.U=0,                                                (4)
```

and `gcd(U)=1`.  The six THM-2309 owner-packet rows have rank six modulo
thirteen.  Adding the selected target axis leaves rank six, while adding the
relation row raises the rank to seven:

```text
(rank(packet),rank(packet+target),rank(packet+a))
  =(6,6,7).                                           (5)
```

Thus `(2)` is a live exact relation whose THM-2334 quotient retains the
two-dimensional target plane.  Equation `(5)` is a relation-lattice statement;
it is not endpoint nonvanishing or ancestry.

## 2. The U_full transplant

Put

```text
U_full=(13,2197,742586,1,183,27,131,53,313)           (6)
```

in relation order.  In THM-2334 current order this is

```text
(1,183,27,131,53,313,13,2197,742586).                 (7)
```

Modulo `22`, the raw and sign-normalized q=11 residues are

```text
raw       =(13,19,20,1,7,5,21,9,5),
normalized=(9,3,2,1,7,5,1,9,5).                      (8)
```

Their literal masks, in the order `(1)`, are

```text
({2,8},{3,7},{5},{0,10},{1,9},{4,6},{0,10},{2,8},{4,6}). (9)
```

The six-mask witness `(1,2,3,5,7,9)` partitions all eleven sheets exactly
once.  The delayed marked intersection has exact positive mass

```text
411318338170045 / 524041104621345129.                 (10)
```

At the exact scales

```text
T_DEN=483730250419703196,
NN=81750412320929840124,                              (11)
```

the coordinate-twist bank has

```text
gamma(0) =248447851579556601771,
gamma(v2)=510954897124935772821,
gamma(v2)-gamma(0)=262507045545379171050               (12)
```

in the first certified finite-field image.  All `168` nontrivial twists
differ from `gamma(0)`.  Exact inverse finite Fourier transformation gives

```text
#{q in F_13^2:A(q)!=0}=169.                           (13)
```

Nonzero reduction in the certified primitive image proves exact algebraic
nonvanishing.  The complete gamma and target banks have respective digests

```text
afcdb043eb1bf8095c313473a3d3bdcf4ce027f86b01b5f3cecbc7c87e6484b3,
726423df1b9e1c93b356966e5c3c386669e2f6b19da0bf8818606204eb2e9ee5. (14)
```

An independent exact-boundary engine obtains interval counts

```text
(|E_0|,|E_v2|,|Q_a|)=(33810,34560,28730)              (15)
```

and matches the fast engine's endpoint values under two certified embeddings.
The values in `(13)` are the unrestricted THM-2334 aggregates `A(q)`.  They
are not the coupled all-`91`-unit projectors `B(q)`.

## 3. The U_clock transplant

Put

```text
U_clock=(65,2197,742586,5,661549,655231,658533,661445,291). (16)
```

Its `13`-adic valuation profile is

```text
(1,3,5,0,0,0,0,0,0).                                 (17)
```

The raw and normalized q=11 residues are

```text
raw       =(21,19,20,5,9,5,7,15,5),
normalized=(1,3,2,5,9,5,7,7,5),                      (18)
```

and the literal masks are

```text
({0,10},{3,7},{5},{4,6},{2,8},{4,6},{1,9},{1,9},{4,6}). (19)
```

Again `(1,2,3,5,7,9)` partitions every q=11 sheet exactly once.  More
strongly, the literal source mask and the delayed owner word occur at the
same fixed clock

```text
c=1/22.                                                (20)
```

The current-order base distances are

```text
(5,9,5,7,7,5,1,3,2)/22,                              (21)
```

and after the delayed dilation `R=169` they are

```text
(9,3,9,5,5,9,7,1,8)/22.                              (22)
```

The four strict margins used by the two Boolean patterns are

```text
(2/77,13/154,2/77,41/154).                            (23)
```

The complete THM-3398 cochain is zero at `(20)`, and every projective wedge
is zero.  This is a common-centre realization, not a nonzero projective
current.  Its exact delayed mass is

```text
1397606991636352080199080533 /
1692517471993352536064760510465.                      (24)
```

The independent boundary engine has interval counts

```text
(147372,147404,136158).                               (25)
```

At the zero twist and `v2`, two certified finite-field embeddings give

```text
gamma(0)=
 (56767723330345680038743661041266194,
  65234233976034532625816110096140982),

gamma(v2)=
 (34870555972766792317398130208739733,
  74671298727704698408794173004769050),               (26)
```

with both differences nonzero.  Hence the U_clock endpoint bank is exactly
nonconstant.  The full `169`-twist inverse transform was not computed, so in
particular the named aggregate `A(e_c2)` remains open.

Finally, the explicit endpoint harmonics

```text
u=(-3,-3,1,2840374,-3,3,2,-3,-48974),
v=(24,24,29,-17270424,38,30,29,24,-49012)             (27)
```

give exact address `(2)` at

```text
(X,m,Y)=(13,1,742599),                                (28)
```

with maximum absolute endpoint heights `2840374` and `17270424`; every
endpoint coordinate is nonzero modulo seven.  This is one THM-2331 atomic
term.  It does not prove the grouped orbit coefficient `C(a;X,m)` nonzero.

## 4. Exact q=27 and q=51 CRT decorations

The q=27 tuple is

```text
(28405,7599423,18279868269,3459,2016,2757,1041,3693,
 11163142875).                                        (29)
```

Its residues modulo `54` are

```text
(1,3,3,3,18,3,15,21,3),                              (30)
```

and `(3,15,18,21)` partitions all `27` sheets exactly once.  The q=51 tuple
is

```text
(70993,7199569,30105550319,5825,4200,5214,7684,4421,
 18313194875),                                        (31)
```

with residues modulo `102`

```text
(1,1,11,11,18,12,34,35,23).                          (32)
```

The seven-mask witness `(1,11,12,18,23,34,35)` has multiplicity profile

```text
1^42 2^4 3^3 4^2.                                    (33)
```

Both lifts preserve

```text
mod 7=(6,6,5,1,0,6,5,4,3),
v_13=(1,3,5,0,0,0,0,0,0).                            (34)
```

For q=51, reduction modulo `34` and the affine lift characters are

```text
bar r=(1,1,11,11,18,12,0,1,23),
k mod 3=(0,0,0,0,0,0,1,1,0).                        (35)
```

On base fibres `y=0,2`, owners `1` and `35` have disjoint active sections;
on `y=1` they agree.  Thus `(35)` retains exactly the THM-3429 affine-lift
hostile which quotient order alone destroys.

## 5. Why this is a two-transplant theorem

The common mechanism is now explicit:

```text
one all-91-unit relation a
  -> primitive owner tuples satisfying a.U=0
  -> literal q=11 masks and a positive delayed intersection
  -> THM-2334 endpoint character tests.               (36)
```

The positive signal is that relation compatibility, literal support, delayed
mass, and endpoint nonconstancy coexist.  The obstruction is the quantifier
"on the same tuple":

```text
U_full:  complete unrestricted 169-target nonvanishing,
         plus the complete guard-packet-row-refined 2197-twist bank and
         all-72 labelled graph-factor nonvanishing,
         but no same-clock common-centre delayed realization proved;

U_clock: same-clock common-centre delayed realization and exact endpoint
         nonconstancy,
         but no complete 169-target or refined-bank calculation. (37)
```

Therefore `(36)` is not yet one physical bridge row.  The split in `(37)` is
the strongest exact boundary, not a cosmetic distinction.

## 6. The guard-packet-row-refined quotient and its labelled back-map

This section is only about `U_full`.  Write `w` in THM-2334 current order and
put

```text
K_w=ker(w:F_13^9 -> F_13),
r_i=e_i-(w_i/w_q1)e_q1 in K_w.                       (R1)
```

Let `L_w` be the span of the six published owner-packet rows.  In the coarse
quotient `G_w=K_w/L_w`, with dual basis `(v1,v2)`, both tuples have

```text
[r_H]=[r_q5]=(1,0).                                  (R2)
```

Thus the lawful coordinate response is `P_i=A_w([r_i])`, not the ill-typed
expression `gamma_w(e_i)`, and every coarse role chart has zero `H--q5`
bridge.  In particular, the coarse 169-target bank by itself cannot power the
carrier determinant.

Delete only the `H`-labelled **row of the six-row owner-packet span**, and
call the resulting rank-five span `L_w^-`.  This algebraic operation does not
delete the Boolean guard condition: both endpoint patterns `PATTERN_E` and
`PATTERN_QA` still contain their literal `guard_safe` factor.  Exact row
reduction gives

```text
dim(K_w/L_w^-)=3,
(L_w^-)^perp/<w> = <v1,v2,e_H>.                      (R3)
```

For the fixed certified primitive embedding with

```text
p=572252886246508880869,  zeta^13=1,  zeta!=1,
```

the refined transform and normalized inverse are

```text
gamma^-(alpha,beta,tau)
  = Gamma_w(alpha v1+beta v2+tau e_H),

A_w^-(x,y,t)
  =13^(-3) sum_{alpha,beta,tau in F_13}
     gamma^-(alpha,beta,tau) zeta^(-alpha x-beta y-tau t). (R4)
```

The quotient projection has the exact back-map

```text
A_w(x,y)=sum_{t in F_13} A_w^-(x,y,t).               (R5)
```

The `tau=0` character slice reproduces both frozen coarse digests, and `(R5)`
is checked at all four coarse bases occurring among the role classes.  The
eight retained role labels have five distinct refined classes:

| role labels | class in `F_13^3` | `A_w^-` modulo `p` |
|---|---:|---:|
| `c1` | `(0,0,0)` | `405336876493642499425` |
| `c3` | `(0,1,0)` | `518539850465495448196` |
| `c2,q3,q4,q5` | `(1,0,0)` | `503604956476841920373` |
| `H` | `(1,0,1)` | `320618948602619577408` |
| `q2` | `(1,12,0)` | `15703541686881447885` |

All five displayed values are distinct.  The surviving coarse collision is
therefore refined just enough to give

```text
A_w^-(1,0,1)-A_w^-(1,0,0)
  =389266878372286537904 mod p !=0.                  (R6)
```

The remaining role collision `c2=q3=q4=q5` is real; individual graph edge
weights need not all be nonzero.  It does not kill the two complete `K4`
tree polynomials.

Fix the owner-order orientation gauge `u_i -> u_j` when `i<j`.  A labelled
role chart sends `H` to the unique hub `u5`, `q5` to the unique leaf `u7`,
bijects `(c1,c2,c3)` with one three-vertex wing and `(q2,q3,q4)` with the
other, and may exchange the two wings.  Hence there are exactly

```text
2 * 3! * 3! = 72                                    (R7)
```

charts, forming one target-graph automorphism orbit.  In every chart the
zero counts for

```text
(H--q5 bridge, left K4 sum, right K4 sum, product)
```

are exactly `(0,0,0,0)`.  A separately implemented graph route reproduces
all 72 rows and the chart digest

```text
b7d8c2c9860e4f1aa542b1c85fdb7b65cf4985aba5a81a84ff3a324834d51c51. (R8)
```

Nonzero reduction proves nonvanishing of these unrestricted cyclotomic
endpoint aggregates and their displayed graph factors.  It does not supply
order, sign, positivity, chronology, or a common ancestry coupling.

The fixed integer relation `a` has refined residue class `(1,0,6)`, not one
of the five role classes above.  Even its residue aggregate is a sum over an
entire congruence class, not the grouped exact-address coefficient
`C(a;X,m)`.  Likewise no all-unit `B(q)` is present: the refined endpoint
transform has no mod-seven mask/character coordinate.  These are type
boundaries, not unperformed simplifications.

## 7. The 7 by 13 carrier: exact equivariance obstruction

This section uses proved THM-3473 and the frozen edge set from its
FINITE-EXACT incidence sidecar.  Its private-owner packets are

```text
{u1,u4,u5,u6}, {u2,u3,u5,u8}, {u5,u7},               (38)
```

whose two-section has thirteen edges.  The proposed connection contract is:

| field | exact content |
|---|---|
| source | the nine-coordinate U_full/U_clock relation tuples and their `F_13^2` coordinate-twist bank |
| target | the thirteen labelled edges of the private-support two-section in `(38)` |
| candidate map | choose one translation-stable thirteen-character fibre and biject it with the thirteen edges |
| predicate sought | nonzero edge weights with nonzero bridge weight and nonzero weighted tree sums on both `K4` blocks |
| information destroyed | exact relation address, q=11 masks, clock, endpoint phase, private-sheet counts, and the `k mod 3` state |
| required sidecar | for a labelled determinant, the role chart and owner-order orientation; for absolute `H^1` or physical realization, phase/holonomy and common ancestry |
| cheapest hostiles | the unique edge `(u5,u7)` and the two sixteen-term `K4` tree polynomials |

There is no native label dictionary: the source labels

```text
(c1,c2,c3,H,q1,q2,q3,q4,q5)                          (39)
```

and the target labels `(u1,...,u8)` are disjoint.  More decisively, a
thirteen-character fibre has a regular translation of order thirteen.  Exact
enumeration of all `8!` vertex permutations gives for the target graph

```text
|Aut(G)|=72,
element orders={1,2,3,4,6},
edge-orbit sizes=(1,6,6).                             (40)
```

An equivariant bijection would conjugate the regular source translation to a
target graph automorphism of order thirteen, contradicting `(40)`.  Hence

```text
there is no C13-equivariant bijection from a native 13-character fibre
to the thirteen private-support edges.                (41)
```

Equation `(41)` does not forbid the explicit 72-chart construction in Section
6.  The two maps have different types.  The impossible map is a
`C13`-equivariant bijection from thirteen source characters to thirteen
edges.  A role chart instead maps eight named relation roles to eight
vertices, then takes their thirteen oriented edge differences; it is labelled,
non-equivariant, and noncanonical.  Thus the positive all-72 calculation does
not weaken `(41)`, and `(41)` does not invalidate it.

For each role chart, the edge-difference current is a coboundary in
`B^1(G;K)`.  All six absolute cycle pairings vanish, so its absolute `H^1`
class is zero.  The signed weighted matrix-tree determinant is not an
absolute-`H^1` invariant and can be nonzero on such a coboundary, once the
owner-order orientation gauge is fixed.  Phase or holonomy remains necessary
for a nonzero absolute `H^1` class or a physical realization, but not for the
finite weighted-determinant statement proved here.

## 8. Exact scope boundary

The package does **not** prove any of the following:

```text
C(a;X,m)!=0 for either tuple;
A(e_c2)!=0 for U_clock;
the 13^3 refined endpoint bank for U_clock;
B(q)!=0 for a coupled all-91-unit target fibre;
a source-native, canonical, or C13-equivariant endpoint-to-edge map;
a nonzero private-support absolute H1 flux or physical phase/holonomy;
a lawful Boolean THM-2512 physical current;
common ancestry or root-character bispectrum nonvanishing;
exclusion of one of the 165 scalar rows;
LRC(14).                                               (42)
```

The projective wedge of U_clock is zero, not a hidden current.  The atomic
term `(27)` is support, not grouped noncancellation.  U_full's `(13)` and
`(R4)` are unrestricted mod-thirteen aggregates, not all-unit `B(q)`.  The
refined endpoint worker first marginalizes its left and right atoms and only
then forms `phase*ax*by`; it emits no shared ancestry key.  Exact hostiles show
that equal left/right marginals and equal products can have different joint
intersection mass, while finite-field residues also forget characteristic-
zero order and sign.  Thus neither the five values nor their graph factors
upgrade to a Boolean/common-ancestry current.  THM-2512's signed ANOVA bundle
and the private-support carrier remain separately typed.  LRC(14) is OPEN.

## 9. Deterministic exact package

Run from the repository root:

```bash
PYTHONHASHSEED=0 python -B 04-computation/lrc_half_twist_relation_current_bridge_thm3479.py
PYTHONHASHSEED=1 python -B -O 04-computation/lrc_half_twist_relation_current_bridge_thm3479.py
```

Both LF-normalized transcripts equal

```text
05-knowledge/results/lrc_half_twist_relation_current_bridge_thm3479.out
```

byte for byte.  The companion uses only the standard library, contains no
`assert`, forbids dynamic execution, and pins every proved dependency.  Its
fast endpoint engine uses periodic subtraction and a scaled sweep; the
independent engine refines every exact boundary, checks exact midpoints, and
directly intersects all word preimages.  The carrier audit independently
enumerates all target graph automorphisms.

The promoted refined package is replayed by

```bash
PYTHONHASHSEED=0 python -B 04-computation/lrc14_guard_deleted_refined_endpoint_role_probe_20260816.py
PYTHONHASHSEED=1 python -B -O 04-computation/lrc14_guard_deleted_refined_endpoint_role_probe_20260816.py
python -B 04-computation/lrc14_guard_deleted_refined_endpoint_role_graph_audit_20260816.py
python -B 04-computation/lrc_endpoint_ufull_frozen_five_common_ancestry_gate_20260816.py
```

The first two commands replay the complete 2,197-twist bank.  The graph audit
reconstructs the 72 factor rows independently of the endpoint engine, and the
common-ancestry gate supplies both a positive scalar Boolean realization and
a hostile showing why the physical-current upgrade is unavailable.

The immutable LF hashes are

```text
script:  ad2a620cdc238f28e3384698b2c612f38cdf2566bd56b76d1cbabcc03107ec0b
output:  e1e4355290c493d8f9518a4bba3825173bb3137c3d50df421816946e6b1f1207
core semantic:
         1bf53086d3da347dde443175462be716da6e9dac54c96582718d19ec8fddff21
carrier semantic:
         3f17c2206feec73da48a989ab2150ceb1c7d1bc275c77291df476d882957581a
```

The exact core companion and pinned finite sidecars prove the scoped
statements above.  The independent end-to-end audit accepted promotion with
the quotient, map, cohomology, guard, and physical-current repairs now written
explicitly into this theorem.
