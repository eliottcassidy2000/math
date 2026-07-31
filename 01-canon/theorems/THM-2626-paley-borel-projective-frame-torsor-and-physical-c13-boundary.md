---
id: THM-2626
title: "Paley--Borel projective frame torsor and physical C13 boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Reduction
  of the modular C2*C3 generators into
  PSL2(F13) produces three order-seven Hurwitz trace classes.  The ordered
  norm of seven conjugates of the parabolic U^a is central exactly when
  t-a, respectively t+a, lies in {plus or minus 3, plus or minus 5,
  plus or minus 6}.  Six oriented charts cover every nonzero exponent, and
  the sharp minimum is three charts, with exactly two minimum-cardinality
  covers.  The P13-coset space is an 84-state projective frame bundle over
  P1(F13) with six-valued tangent fibre; its 78-state affine locus is a
  regular torsor for the Frobenius group C13 semidirect C6 and the arc set
  of Paley(13), whose local frame graph is C6.  This endpoint model is only
  Borel-equivariant.  Over every odd finite field, an ordered endpoint basis
  is recovered losslessly from its difference, Gram matrix, and oriented
  determinant; over F13 this recovers the 84-state projective frame, while
  over F2 the identity/endpoint-swap pair is an exact collision.  The LRC
  application still lacks a physical metric intertwiner.  Every closing chart has one all-affine
  seven-cycle lifting to six frame cycles.  THM-2605 realizes the underlying
  root word positively with one coherent affine phase, but physically supplies
  only the normal C13 translation layer: its positive Paley-pair kernel is
  flat in all six frame directions.  The same root word can have tangent
  holonomy one or four; a physical future-owner/frame edge
  and endpoint residues remain absent.  This is a finite nonabelian design theorem, not an LRC row
  exclusion or a proof of LRC(14).
source: codex-2026-07-27-paley-borel-frame
depends_on:
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2602-commutative-vertex-insertion-and-ordered-transition-curvature-no-go
  - THM-2605-inverse-root-dipole-connection-and-mixed-square-invoice
related:
  - THM-2603-hurwitz-projective-root-owner-atlas-and-nonabelian-seven-edge-trivialization
  - THM-2597-six-vertex-bicycle-modular-abelianization-cycle
  - THM-2599-rootwise-opposite-shift-paired-slice-law
  - THM-2593-charged-target-section-atlas-and-minimal-c91-holonomy-trivialization
  - THM-2600-constant-six-middle-rail-common-x-atlas-and-uniform-bockstein-section
  - THM-2601-linear-bockstein-sheet-coordinate-and-nonlinear-target-monodromy
  - THM-2604-unshifted-future-root-accessibility-and-selector-cross-mixing-boundary
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2620-endpoint-pair-parabolic-transvection-and-translation-gauge-boundary
  - THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary
script: 04-computation/psl2f13_frame_norm_thm2626.py
output: 05-knowledge/results/psl2f13_frame_norm_thm2626.out
script_sha256: 5e3e6b20d48690abcede98b2142695be5dbf144ddfbeb7cb9327c17a9f39d4da
output_sha256: 05eaf0a62a7beba59d39b156103bbdc35ebe6e4ad06396ad5be6b02ae7ac8d54
hash_basis: LF-normalized bytes
---

# THM-2626 -- the modular frame sees ordered seven-edge curvature

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2603 already owns the Hurwitz quotient, ordered norm atlas, and abstract
84-frame carrier.  This theorem retains that inheritance for provenance and
adds the nonoverlapping affine Borel/Paley model, the exact physical
`C13`-only and frame-flat boundary, and the oriented Gram reconstruction
with its characteristic-two failure.

THM-2602 proves that commutative root-vertex insertion cannot recover an
ordered transition.  The modular free factors provide an exact finite model
of the missing operation.  In `PSL_2(F_13)`, a parabolic order-thirteen root
shift and an order-seven projective frame motion can have a trivial ordered
norm in one orientation and a nontrivial norm in the reverse orientation.
The phenomenon is trace-theoretic, covers every nonzero parabolic exponent,
and has a sharp three-chart atlas.

The same group also explains the recurring number `84`.  Its quotient by a
Sylow-thirteen translation subgroup is a projective root together with one
of six tangent frames:

```text
84 = 14 projective roots * 6 frame colours
   = 78 affine root/frame states + 6 states over infinity.       (1)
```

This is the strongest valid bridge and also the stopping boundary.  The
current LRC carrier supplies neither the six missing infinity profiles nor a
physical tangent-frame transition.  More sharply, THM-2605 already realizes
every finite root word positively with one coherent affine paired phase,
while a root word alone does not determine the frame holonomy.  The exact new sidecar must therefore be an ordered
root-plus-owner-frame edge, not another root accessibility lemma.

## 1. The modular `C2*C3` quotient is literal

Work over `F_13` and use projective matrix equality.  Put

```text
S   =[[0, 1],[-1,0]],
R   =[[0, 1],[-1,1]],
U^a =[[1, a],[ 0,1]],
g_t =[[0, 1],[-1,t]].                                    (2)
```

Then

```text
S^2=R^3=1,                 U=R^(-1)S,
g_t=S(SR)^t                in PSL_2(F_13).                (3)
```

The matrices `S,U` generate `PSL_2(F_13)`, so (3) is the reduction of the
standard surjection

```text
PSL_2(Z)=C2*C3 -> PSL_2(F_13).                             (4)
```

For

```text
t in {3,5,6},
```

direct Cayley--Hamilton reduction gives

```text
g_t^7=-I.                                                 (5)
```

Thus each `g_t` has projective order seven.  The three trace-square classes
have one further exact arithmetic symmetry.  Squaring an order-seven
generator induces

```text
y |-> (y-2)^2,
9 -> 10 -> 12 -> 9,                                      (6)
```

which is the cyclic group
`(Z/7Z)^*/{plus or minus 1}`.  Equation (6) is an arithmetic `C3` on trace
classes.  It is not an assertion that the three norm charts form a physical
LRC `C3` orbit, nor an identification with the order-three modular generator
on a Boolean carrier.

## 2. The order-seven trace polynomial

Let `M` be a noncentral determinant-one matrix with trace `x`.  Define

```text
f_0=0, f_1=1, f_(n+1)=x f_n-f_(n-1).
```

Cayley--Hamilton gives

```text
M^7=f_7(x)M-f_6(x)I.                                     (7)
```

In `F_13`,

```text
f_7(x)
 =x^6-5x^4+6x^2-1
 =(x-3)(x-5)(x-6)(x+3)(x+5)(x+6).                       (8)
```

Consequently

```text
M^7 is scalar
 iff tr(M) lies in T_7={plus or minus 3, plus or minus 5,
                         plus or minus 6}.                (9)
```

No scalar exception occurs in the application below: both relevant matrices
have a fixed off-diagonal entry `plus or minus 1`.

## 3. Forward and reverse ordered norms

For `a in F_13`, define seven conjugate parabolics

```text
A_k=g_t^k U^a g_t^(-k),                 0<=k<=6.          (10)
```

Order is load-bearing.  Telescoping before any entrywise calculation gives

```text
N_t^+(a)=A_0 A_1 ... A_6
        =(U^a g_t)^7 g_t^(-7),

N_t^-(a)=A_6 A_5 ... A_0
        =g_t^7(g_t^(-1)U^a)^7.                           (11)
```

The two moving traces are

```text
tr(U^a g_t)=t-a,                 tr(g_t^(-1)U^a)=t+a.    (12)
```

Equations (9)--(12) prove the exact criterion

```text
N_t^+(a) is central iff t-a in T_7,
N_t^-(a) is central iff t+a in T_7.                       (13)
```

Restricting to nonzero parabolic exponents gives the six oriented charts:

| trace `t` | forward exponents | reverse exponents |
|---:|---|---|
| `3` | `{6,8,9,10,11}` | `{2,3,4,5,7}` |
| `5` | `{2,8,10,11,12}` | `{1,2,3,5,11}` |
| `6` | `{1,3,9,11,12}` | `{1,2,4,10,12}` |

Their union is all of `F_13^*`.  The statement concerns every nonzero
**parabolic exponent or marker step** `a`; it must not be paraphrased as a
claim about every physical root point.

The polynomial factorization reported by the original exact probe is also
explained.  At `t=3`, the lower-left and upper-right entries of the forward
norm are respectively

```text
c(a)=a(a-6)(a-8)(a-9)(a-10)(a-11),

b(a)=10a(a-4)(a-6)(a-8)(a-9)(a-10)(a-11).               (14)
```

The common nonzero roots in (14) are precisely the forward row at `t=3`.

### 3.1 Orientation is genuinely noncommutative

At `(t,a)=(3,6)`,

```text
N_3^+(6)=-I,
N_3^-(6)=[[3,4],[1,6]],                                  (15)
```

while at `a=7`,

```text
N_3^+(7)=[[6,9],[12,3]],
N_3^-(7)=-I.                                              (16)
```

By contrast, for every homomorphism `phi` to an abelian group, conjugation
and order disappear:

```text
phi(N_t^+(a))=phi(N_t^-(a))=phi(U^a)^7.                  (17)
```

This is an exact finite realization of THM-2602's boundary.  The global
abelianization `PSL_2(Z)^ab=C6` sees neither (15) nor (16).

## 4. The sharp minimum-cardinality oriented atlas has three charts

Each of the six oriented charts in Section 3 contains five of the twelve
nonzero exponents.  Two charts therefore cannot cover `F_13^*`.  Exhausting
the `64` subfamilies gives exactly two covers of size three:

```text
(3,+),(3,-),(6,+),
(3,+),(3,-),(6,-).                                       (18)
```

Thus three is sharp, and the two rows in (18) are the only
minimum-cardinality covers.  They are not the only inclusion-minimal
covers: for example all four charts at `t=3,5` form a genuine
inclusion-minimal cover.  This distinction is why the exact claim is an
optimization theorem, not a complete blocker classification.
This triple is an oriented set-cover theorem.  It is not canonically the
order-three free factor and does not select its own chart on an LRC packet.

## 5. The exact 84-state projective frame bundle

Put

```text
G=PSL_2(F_13),                     |G|=1092,
P_13=<U>,                          |P_13|=13,
B=N_G(P_13),                       |B|=78.                (19)
```

Here `B` is the upper-triangular Borel.  The natural projection

```text
G/P_13 -> G/B=P^1(F_13)                                    (20)
```

has

```text
|G/P_13|=84,                   |P^1(F_13)|=14,
B/P_13 = F_13^*/{plus or minus 1}
       = QR_13
       isomorphic to C6.                                  (21)
```

Every fibre in (20) therefore has six states.  After choosing the affine
chart, its preimage has

```text
13*6=78                                                    (22)
```

states, and the six remaining states lie over infinity.  One explicit,
noncanonical comparison with the six labelled septimal owner colours is

```text
chi:QR_13 -> F_7^*,                  4^j |-> 3^j,

1->1, 4->3, 3->2, 12->6, 9->4, 10->5.                   (23)
```

The word **noncanonical** is essential.  The theorem constructs no physical
identification between the frame fibre (23) and an LRC owner character.

### 5.1 The affine chart is a regular Frobenius-group torsor

The product count in (22) has a stronger intrinsic law.  Write an element
of the upper Borel as

```text
h=[[r,b],[0,r^(-1)]].                                    (23a)
```

Represent a frame over an affine point by `(x,eta)` with
`eta in QR_13`.  Direct projective and tangent transport gives

```text
h.(x,eta)=(r^2 x+rb,r^2 eta).                            (23b)
```

The map `h |-> h.(0,1)=(rb,r^2)` is a bijection

```text
B -> F_13 x QR_13.                                       (23c)
```

Hence the `78` affine states are a regular left torsor for

```text
B=C_13 semidirect QR_13 = C_13 semidirect C_6,
(c,kappa)(d,lambda)=(c+kappa d,kappa lambda).             (23d)
```

The normal `P_13` has six thirteen-state orbits, one at each frame value;
the quotient `B/P_13=C_6` permutes those layers.  The fibre over infinity
is one six-state Borel orbit rather than part of the regular affine orbit.

This supplies the exact co-occurrence law which a root/chronology sidecar
and an owner-frame sidecar would have to satisfy.  Around a cyclic path of
affine transitions `(c_i,kappa_i)`, frame holonomy is
`product_i kappa_i`, while root holonomy is the **twisted** sum

```text
c_1+kappa_1 c_2+kappa_1 kappa_2 c_3+... .                (23e)
```

Separate vanishing tests for an untwisted `C_13` sum and a `C_6` product do
not imply (23e).  This is an abstract semidirect torsor model, not a
physical old/future-`C_13` bibundle: identifying a chronological deck with
`P_13` and choosing a lawful frame layer remain missing.

### 5.2 The affine torsor is the Paley arc set, but only affinely

There is a canonical co-occurrence graph inside the abstract affine model.
Let `Paley(13)` have vertex set `F_13` and join `x,y` when
`y-x in QR_13`.  Since `-1 in QR_13`, this is an undirected graph, not a
tournament.  The map

```text
(x,eta) |-> the oriented edge x -> x+eta                  (23g)
```

identifies the `78` affine frame states with its `78` arcs, and (23b) is
exactly the induced Borel action on those arcs.  Hence `B` is regular on the
arc set.  The graph has parameters

```text
SRG(13,6,2,3),             39 edges,                     (23h)
```

and the six neighbours of zero induce the literal cycle

```text
1-4-3-12-9-10-1.                                         (23i)
```

Thus the six frame colours acquire a canonical local `C_6` co-occurrence
relation.  It is closely related to THM-2597's partial-cube `C_6`, but the
whole Paley graph is not a partial cube: (23h) gives triangles through every
edge.  A nonsquare multiplier swaps the Paley graph with its complement;
the full affine group is regular on all `156` non-diagonal ordered pairs.

The endpoint reading (23g) stops exactly at the Borel boundary.  For

```text
M=[[a,b],[c,d]],
```

the natural second endpoint and the tangent frame differ by

```text
M(x+eta)-M(x)=eta/[(cx+d)(c(x+eta)+d)],
tangent increment=eta/(cx+d)^2.                           (23j)
```

Away from poles and for `eta!=0`, equality forces `c=0`.  In the closing
example `M=[[7,6],[12,3]]`, the Paley arc `0->1` has transported base `2`;
its natural second endpoint is `0`, whereas tangent transport labels the
virtual edge `2->5`.  Therefore the Paley model does not turn the Hurwitz
tangent cycle into a genuine two-root path or close the LRC owner edge.

There is now a proved but still physically conditional chronological
comparison.  THM-2611 proves that faithfully restoring its forgotten
chronological orbit has the sharp form of a `C_13` torsor.  Abstractly
combining that torsor with the frame fibre `B/P_13` has

```text
13*6=78                                                    (23f)
```

states, exactly the size of the affine chart in (22), while projective
compactification adds the six frames over infinity.  Equation (23f) is now
explained by the extension (23d), but no constructed physical
intertwiner identifies THM-2611's chronological `C_13` with the subgroup
`P_13=<U>` used here.

There are also two different order-six objects here.  The fibre `B/P_13` is
a local subgroup quotient.  The modular abelianization in Section 3.1 is a
global quotient which destroys order.  They are abstractly cyclic of the
same size but are not the same construction.  The exact companion verifies
that the normal closure of `[S,R]` is all `1092` elements of `G`, so `G`
itself has no nontrivial abelian quotient through which (23) could descend.

## 6. The exact `6+7` affine boundary

On the fourteen projective points, the three Hurwitz matrices have cycles

```text
g_3: (infinity,0,9,2,1,7,3)   (4,12,10,11,8,5,6),

g_5: (infinity,0,8,4,1,10,5)  (2,9,3,7,6,12,11),

g_6: (infinity,0,11,5,1,8,6)  (2,10,3,9,4,7,12).        (24)
```

Thus deleting infinity leaves a directed six-vertex path and an all-affine
seven-cycle.  This is an intrinsic partial permutation `P6 disjoint union
C7`, not a tournament.

The same split holds for every closing norm chart.  Put

```text
M_(t,a,+)=U^a g_t,
M_(t,a,-)=g_t^(-1)U^a.                                  (25)
```

For each of the `30` labelled `(t,orientation,a)` pairs in Section 3,
equation (13) says that (25) has projective order seven.  Its fixed-point
discriminant is one of

```text
5,6,8,                                                       (26)
```

all nonsquares in `F_13`.  It therefore acts fixed-point-freely on
`P^1(F_13)`, as two seven-cycles.  Exactly one contains infinity, so the
other consists of seven affine points.

The preimage of that affine cycle under (20) has `42` states.  Since an
order-seven element cannot fix a `P_13` coset, those states split into six
seven-cycles.  Across the thirty labelled chart/exponent choices this gives

```text
30*6=180 labelled chart/exponent frame cycles.            (27)
```

Equation (27) is a labelled count.  Distinct chart labels are not asserted
to give distinct underlying subsets.

## 7. The tangent-frame cocycle

The six-frame lift is explicit.  For

```text
M=[[alpha,beta],[gamma,delta]] in SL_2(F_13)
```

and an affine point away from the pole, put

```text
M.(x,eta)
 =((alpha x+beta)/(gamma x+delta),
   eta(gamma x+delta)^(-2)),             eta in QR_13.    (28)
```

The square makes (28) independent of the sign of the `SL_2` representative.
On a finite seven-cycle of a closing matrix (25), the seven tangent factors
multiply to one because `M^7` is projectively the identity.  Consequently
each of the six initial frames returns.

For the sharp forward chart `(t,a)=(3,6)`,

```text
M=[[7,6],[12,3]],

finite cycle:       (0,2,7,9,8,11,1),
QR_13 factors:      (3,1,9,4,12,12,10),
chi-owner factors:  (2,1,4,3,6,6,5).                    (29)
```

Both products in (29) are one in their respective groups.  The companion
checks the lift twice: once from (28), and independently by enumerating the
exact action on all `84` cosets.

The decomposition in Sections 5--7 has suggestive sizes, but it licenses no
silent semantic naming.  In particular:

```text
six affine points in the broken cycle  -/-> six direct LRC roles;
seven points in the closed cycle       -/-> seven owner clocks;
six tangent frames                     -/-> physical owner characters.   (30)
```

Every arrow in (30) needs a constructed map.

## 8. Positive physical root words already exist

Fix any nonzero exponent and choose one closing chart from (18).  Let

```text
(h_0,h_1,...,h_6)
```

be its all-affine cycle, ordered by (25), and append `h_7=h_0`.  THM-2605
proves that one invariant affine phase realizes every ordered physical root
word on a positive same-ancestry opposite-shift paired carrier.  Therefore,
for its common cylinder depth `ell`,

```text
mu(E(h_0,...,h_6,h_0;N))=13^(-8ell)>0                  (31)
```

for every permitted block delay `N>=ell`.  This applies after THM-2605's
lawful all-root role selection on every hypothetical scalar-cover row and,
unlike the earlier rootwise THM-2599 choice, retains one coherent affine
shift phase along the entire word.

Equation (31) is a real positive Boolean root itinerary.  It is not a
projective frame path.  Its delay blocks are not identified with the seven
THM-2542 chart edges; its root digit is not identified with a target section
or relation residue.  Its fixed source role/paired blocker is source-side
provenance; no owner/repair endpoint or future owner is retained at the
eight stops.

Thus finite root-word accessibility is saturated.  Another tree language,
finite root atlas, or forbidden-word search cannot supply the coordinate
which (31) deliberately forgets.

### 8.1 The proved physical layer is `C_13`, not `C_13 semidirect C_6`

The fixed-phase typing in THM-2605 makes the stopping boundary exact.  For a
fixed physical role `k` and affine phase `v`, its root states are

```text
P_v={(r,q): q=v-kr}.                                    (31a)
```

The lawful opposite shift is the normal translation action

```text
(r,q) |-> (r+Delta,q-k Delta).                           (31b)
```

A general abstract Borel transition `r'=s r+t`, with `s in QR_13`, would
instead require

```text
q'=s q+(1-s)v-kt.                                       (31c)
```

The companion checks (31c) on all `158,184` choices of
`(k,v,s,t,r)`.  The root displacement is uniform in `r` exactly when `s=1`.
Thus (31b) realizes the normal `C_13`; no proved operation realizes the five
nontrivial `C_6` scalings.

There is nevertheless a positive static Paley survivor.  If `ell` is the
common root-cylinder depth in THM-2605, then for every sufficiently separated
two-time block (`N>=ell`) and every Paley arc `(r,eta)`,

```text
K_(v,N)(r,eta)
 =mu(D_r intersect T^(-N)D_(r+eta))
 =13^(-2ell)>0.                                         (31d)
```

This is genuine simultaneous root co-occurrence, but it is constant in the
six square directions `eta`.  Consequently all five nontrivial `C_6`
Fourier modes vanish.  The repo has therefore constructed a positive Paley
**edge kernel**, not a positive tangent-frame observable.  The next decisive
object is a positive carrier-valued quadratic-residue observable that obeys
the scaling law in (23b); relabelling a fixed role by `k^2` does not do so.

## 9. Same root path, different frame holonomy

The loss in Section 8 is sharp within `PSL_2(F_13)`.  For any prescribed
affine edge `x->y`, the translation

```text
T_(y-x)
```

sends `x` to `y` and has tangent multiplier one.  The determinant-one matrix

```text
H_(x,y)=[[2,7(y-4x)],[0,7]]                              (32)
```

also sends `x` to `y`, but its tangent multiplier is four.

Use the seven-root cycle in (29).  One lift takes translations on all seven
designated edges.  A second lift replaces only the edge `0->2` by

```text
[[2,1],[0,7]].                                            (33)
```

Both paths visit exactly

```text
(0,2,7,9,8,11,1,0),                                     (34)
```

but their tangent returns are respectively

```text
1 and 4 in QR_13,                                        (35)
```

or, after (23),

```text
1 and 3 in F_7^*.                                        (36)
```

Hence even the complete ordered root word (34) does not determine the
six-valued frame holonomy.  This is stronger than an endpoint-marginal
hostile: every successive root pair agrees.

## 10. The cheapest lawful LRC sidecar

The proved finite-group implication is

```text
chosen oriented Hurwitz chart + nonzero parabolic exponent
 -> central ordered norm in PSL_2(F_13)
 -> six formal all-affine frame cycles.                   (37)
```

To turn (37) into the ordered object demanded by THM-2602, a physical
construction must supply seven positive joint kernels on one common ancestry
carrier,

```text
K_j((q,kappa),(q',kappa'))>=0,             0<=j<=6,       (38)
```

with all of the following retained before marginalization:

```text
adjacency: outgoing (q',kappa') at j is incoming at j+1;

root law: q'=M_j q;

frame law:
 kappa'=kappa chi((gamma_j q+delta_j)^(-2));

positivity: one of the six lifted all-affine cycles is one positive
            sevenfold fibre product.                     (39)
```

For the root-plus-frame transition debt alone, (38)--(39) are the cheapest
decisive test.  A practical probe is to refine one THM-2605 affine cylinder word by
the future owner at all seven stops and inspect the six lifted cycles.

Inside the affine Borel there is an equivalent two-root formulation.  Use
(23g) and write the frame state as a pair `(q,z)` with
`eta=z-q in QR_13`.  Two such pairs determine a unique Borel transition:

```text
kappa=(z'-q')/(z-q) in QR_13,
c=q'-kappa q,
(q',z')=(kappa q+c,kappa z+c).                            (39a)
```

Thus the cheapest purely root-valued producer is one positive joint kernel
for **two simultaneous physical roots** satisfying (39a) on common
ancestry.  THM-2605 realizes either prescribed finite root word with one
coherent phase, but it does not prove positivity for their intersection;
applying it twice and multiplying the masses would be an illegal
independence step.  Equation (39a) turns the owner-frame debt into a concrete
paired-word experiment without confusing it with the non-affine Hurwitz
endpoint transport ruled out by (23j).

For a THM-2334 relation current, one further field is still mandatory:
both endpoint relation residues must survive on the same path, followed by
a nonzero difference character.  The owner-frame lift does not manufacture
that endpoint coordinate.

There is also a fixed-gauge warning.  The central product in (11) is the
holonomy of an **interleaved projective connection**.  It includes the
`g_t` frame motion.  Until that motion is physically realized and compared
with THM-2542's affine deck, centrality of `N_t^+` or `N_t^-` is not the
identity

```text
K=S_(-7a)                                                  (40)
```

in THM-2602's existing thirteen-root gauge.  Silently deleting the
interleaved frame changes the theorem.

## 11. An oriented Gram sidecar recovers the frame in odd characteristic

The missing finite frame is not intrinsically unrecoverable.  It is lost by
a specific quotient.  Let `k=F_Q` have odd characteristic, let

```text
U=[L R] in GL_2(k),        v=L-R,
G=U^T U=[[A,B],[B,C]],     Delta=det(U).                  (41)
```

Then

```text
U |-> (v,G,Delta)                                           (42)
```

is injective.  If another basis has the same triple, its ratio with `U` lies
in `SO_2(k)` and fixes the nonzero vector `v`.  It is therefore the identity.
For anisotropic `v` this follows from the orthogonal splitting.  For
isotropic `v`, choose `h` with `det(v,h)=1`; writing `Oh=h+t v`, preservation
of the bilinear form forces `t=0`.

There is an explicit inverse.  Put

```text
n=v.v=A-2B+C,             mu=v.R=B-C,
v=(x,y).                                                    (43)
```

If `n!=0`, then

```text
R=n^-1 (x mu-y Delta, y mu+x Delta),      L=R+v.           (44)
```

If `n=0`, then `x!=0`.  With

```text
h=(0,x^-1),       kappa=v.h=y/x,
mu=Delta kappa,
t=(C-Delta^2 h.h)/(2 Delta kappa),                        (45)
R=Delta h+t v,             L=R+v.                          (46)
```

The image satisfies `v.v=A-2B+C` and `det(G)=Delta^2`; in the isotropic
branch the oriented compatibility in (45) is additional.

Equivalently,

```text
C_U=[v,Delta^-1 R] in SL_2(k),
U <-> (Delta,C_U) in k^* x SL_2(k).                        (47)
```

At `Q=13`, all `26,208` endpoint bases have distinct full keys.  Forgetting
`Delta` leaves exactly `3,744` singleton isotropic keys and `11,232` double
keys.  Projection of `C_U` to `PSL_2(F_13)/P_13` gives the 84-state frame
bundle, uniformly with `312=26*12` endpoint bases per frame.

The prime-two boundary is sharp.  Over `F_3`, all 48 full keys are distinct;
over `F_2`, the identity basis and endpoint-swap basis have the same

```text
v=(1,1),       G=I,       Delta=1.                         (48)
```

Thus characteristic two destroys orientation because `-1=1`; neither the
bare `V4` state nor a partial-cube/graceful shadow can restore it.

This is a static algebraic reconstruction theorem, not yet a physical LRC
lift.  THM-2620 supplies an endpoint allocation plane but no canonical dot
product agreeing with THM-2596's lattice metric modulo 13.  A lawful
intertwiner must first retain that Gram form on one positive chronological
carrier.  Even then (42) supplies neither the affine translation `c` required
by THM-2622, nor owner/word semantics, endpoint phase, or a relation current.

## 12. Scope ledger

The following stronger readings are false or unproved.

1. The `84` cosets and the regular affine Borel torsor are not a constructed
   physical LRC packet or old/future-`C_13` bibundle.
2. The six states over infinity are not supplied by the present target
   atlas.  Restriction to the `42`-state all-affine branch is formal until a
   positive owner-frame lift is built.
3. The comparison `chi` in (23), the oriented chart in (18), and the affine
   root/target identification are choices, not canonical LRC selectors.
4. A scalar norm `plus or minus I` is an algebraic identity, not positive
   mass, a nonzero Fourier coefficient, or a relation current.
5. THM-2600's uniform `q=0` section is a vertex and does not move around a
   projective cycle.  THM-2601's nonlinear sheet coordinate does not
   implement the Mobius transformation (28).  THM-2593 explicitly assumes
   the root/target affine identification which is still missing physically.
6. THM-2605 gives the root word and one coherent paired affine phase, and
   (31d) gives a positive static Paley edge kernel.  It supplies only the
   normal `C_13` action: the kernel is frame-flat, so it does not give the
   future owner-frame sequence or terminal endpoint transport.
7. Neither the arithmetic `C3`, the `6+7` split, nor the six tangent states
   define a tournament or a physical free-factor action.
8. No scalar row is removed.  The LRC ledger remains unchanged and LRC(14)
   remains open.

## 13. Exact companion

Run

```text
python 04-computation/psl2f13_frame_norm_thm2626.py
python -O 04-computation/psl2f13_frame_norm_thm2626.py
```

Both modes must byte-match

```text
05-knowledge/results/psl2f13_frame_norm_thm2626.out
```

after LF normalization.  The dependency-free script uses only integer
arithmetic and explicit exceptions; it contains no floating point, random
fixture, external algebra package, or truth-bearing `assert`.  It verifies:

- all `2,184` determinant-one matrices and `1,092` projective classes;
- the modular identities (3), all three order-seven powers, and the full
  `1,092`-element normal closure of `[S,R]`;
- subgroup sizes `13,78`, coset counts `84,14`, and all fourteen fibre-six
  invoices, together with the regular affine Borel action and its exact
  semidirect-product law;
- the `SRG(13,6,2,3)` Paley arc model, its local frame `C_6`, the affine
  complement swap, and the explicit non-affine endpoint/tangent hostile;
- the `158,184`-cell fixed-phase action (31c) and exact vanishing of all five
  charged modes of the positive static Paley kernel;
- all oriented endpoint Gram keys over `F_3,F_5,F_7,F_13`, both reconstruction
  branches, the determinant-forgetting fibre formulas, the uniform
  `312`-endpoint invoice over all 84 frames, covariance under modular
  coordinate changes, and the characteristic-two identity/swap hostile;
- the exact trace factorization, the arithmetic trace-class cycle, every
  direct/telescoped forward and reverse norm, and (14)--(16);
- all `64` oriented chart subfamilies and the two sharp covers (18);
- the projective cycles (24), all thirty labelled closing pairs, and their
  fixed-point discriminants;
- independently, twelve frame orbits per closing pair and six over each
  all-affine base cycle;
- all `180` labelled chart/exponent tangent-frame cycles from (28), including
  the frozen sequence (29); and
- the same-root/different-frame hostile (32)--(36), together with explicit
  zero sentinels for physical-carrier, relation-current, and row-exclusion
  claims.

The finite proof, group action, frame typing, and LRC scope passed an
independent hostile audit.  The companion was replayed in normal and
optimized modes against the stored transcript after the Borel-torsor
strengthening.
