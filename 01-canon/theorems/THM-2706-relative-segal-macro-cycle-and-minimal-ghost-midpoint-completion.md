---
id: THM-2706
title: "Relative Segal macro cycle and minimal ghost-midpoint completion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For a
  deterministic map T and carrier W, the
  fixed-grading composition map from two degree-one nerve arrows to a
  degree-two endpoint arrow is injective with exact defect
  {x,T^2x in W, Tx notin W}.  A central commuting involution can make this
  defect survive after even holonomy cancels.  On the LRC delayed carrier the
  source defect has exact certified mass floor 1/2635437714, and a new fully
  decorated ordinary-lattice endpoint-sampled D^2 macro cycle runs between
  4/17 and 13/17 through forced forbidden midpoints 1/17 and 16/17.  Every integral
  degree-one affine factorization has the same forbidden midpoint.  The
  minimal ghost-midpoint completion is exact, but no semantic/current cospan,
  principal C2 bibundle, row exclusion, or LRC conclusion is proved.
source: root/relative-segal-macro-cycle-2026-07-28
depends_on:
  - THM-2616-cross-time-target-future-diagonal-and-principal-action-no-go
  - THM-2623-guard-safe-danger-cospan-and-residual-unit-wall
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2693-odometer-skew-product-three-event-escape-and-uniform-delayed-depth-four-nilpotence
  - THM-2698-central-half-odometer-full-local-cycle-and-semantic-sidecar-boundary
related:
  - THM-2292-common-catalytic-section-and-helly-calibration-nerve
  - THM-2315-marked-target-gain-corolla-and-pairwise-composition-boundary
  - THM-2518-perron-inverse-branch-owner-word-cospan-recovery
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2644-odd-torsor-purity-return-gate-and-nonlinear-fixed-branch-decoder
  - THM-2658-balanced-lift-helly-circular-arc-gain-nerve-and-wrap-boundary
  - THM-2680-dilation-reversed-two-edge-clock-fibre-products-and-source-drift-boundary
  - THM-2697-filtered-affine-handoff-germ-category-and-base-signature-holotopy-boundary
  - THM-2701-literal-singleton-word-one-step-dilation-nilpotence
script: 04-computation/lrc14_relative_segal_macro_cycle_thm2706.py
output: 05-knowledge/results/lrc14_relative_segal_macro_cycle_thm2706.out
script_sha256: 3083256d4735357a8bf0d9216fa5c8e318ac737fc8a9ebe7d58e038a9042ed6a
output_sha256: 4e8c054f182aece20e55e244c2de64ba0c5bd418e0637d05d6e8781c2adbcbb3
secondary_script: 04-computation/lrc14_relative_segal_ghost_preword_transit_thm2706.py
secondary_output: 05-knowledge/results/lrc14_relative_segal_ghost_preword_transit_thm2706.out
secondary_script_sha256: 4c81a475bb7380f441b2ea1ad435decca67e8d0523961db98a1e8a62169f2ae4
secondary_output_sha256: afde267be8baf80c25d7f419c8a75b05baa066cb4bd6427f684dbc8644b0a725
hash_basis: LF-normalized bytes
---

# THM-2706 -- relative Segal macro cycle and minimal ghost-midpoint completion

## Status and scope

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The proof below
builds an abstract graded action category and an exact LRC carrier inside it.
It does not claim that the corresponding semantic/current category has been
built.

## 1. Endpoint arrows and unit-step factorizations

Let `T:X->X` be a map and `W subset X`.  Let `C_W(T)` be the full subcategory
on objects `W` of the action category `X semidirect N`: a degree-`n` arrow is

```text
(x,n): x -> T^n x,                 x,T^n x in W.          (1)
```

Composition adds degrees.  Write

```text
E_n^W(T)=Mor_n(C_W(T))
        ={(x,T^n x):x,T^n x in W},                        (2)

N_n^W(T)=E_1^W(T) x_W ... x_W E_1^W(T)
        ={(x,Tx,...,T^n x):T^j x in W for every 0<=j<=n}. (3)
```

Thus `N_n^W(T)` is the `(1,...,1)` multidegree stratum of the ordinary nerve
of `C_W(T)`.  Composition gives the forget-intermediate map

```text
c_n:N_n^W(T)->E_n^W(T),
(x,Tx,...,T^n x) |-> (x,T^n x).                           (4)
```

Because `T` is deterministic, `c_n` is injective.  At degree two its image
and exact defect are

```text
im(c_2)={(x,T^2x):x,Tx,T^2x in W},

Def_2(W,T)=E_2^W(T)\im(c_2)
 ={(x,T^2x):x,T^2x in W, Tx notin W}.                    (5)
```

Its source projection is

```text
SrcDef_2(W,T)={x:x,T^2x in W, Tx notin W}.                (5a)
```

Consequently `c_2` is surjective if and only if

```text
W intersection T^(-2)W subset T^(-1)W.                  (6)
```

This is not a failure of the Segal axiom: `C_W(T)` is an honest category and
its nerve is Segal.  It is strict nonsurjectivity of the **fixed-grading
factorization map** from `(1,1)` nerve simplices to degree-two endpoint
arrows.  Equivalently, `C_W(T)` need not be generated in degree one.

There is an equivalent relation-theoretic formulation which keeps that
distinction explicit.  Put

```text
R_n=E_n^W(T) subset W x W,              R_0=Delta_W.       (6a)
```

With ordinary composition of relations,

```text
R_m after R_n subset R_(m+n).                              (6b)
```

Thus `n |-> R_n` is a normalized lax action of `N` in the
inclusion-ordered monoidal poset `Rel(W)`.  Its first strictness is exactly

```text
R_1 after R_1=im(c_2) subsetneq R_2                       (6c)
```

when `Def_2(W,T)` is nonempty.  This is a lax-action/fixed-grading
factorization defect, not the failure of an unbuilt category.

## 2. Central-involution subdivision lemma

Let `H:X->X` satisfy

```text
H^2=id,                       HT=TH,
S=HT.                                                       (7)
```

Then

```text
S^2=T^2.                                                    (8)
```

If

```text
x,Sx,S^2x in W,                 Tx notin W,                (9)
```

the two unit `S` arrows form a path in the enlarged commuting
`<H,T>`-action category (equivalently in `C_W(S)`).  If that central
extension is explicitly parity-graded, with `H` representing the nontrivial
element of `C_2` and each `S` edge assigned parity one, their total parity is
zero.  Since `S^2=T^2`, their composite has the
same underlying endpoint transformation/relation as the degree-two `T`
arrow `(x,T^2x)`.  Yet by `(5)` that `T` endpoint arrow has no subdivision
into two unit `T` arrows in `W`.  Thus trivial even holonomy does not make
one-step descent effective.

The hypotheses are sharp.  If `Tx in W`, the deterministic factorization is
present and unique.  If `HT!=TH`, then `(HT)^2` need not equal `T^2`.  If the
middle object is forgotten rather than retained, endpoint data alone cannot
distinguish the two cases.

## 3. Quantitative raw delayed-word defect

Put

```text
p=13,      T(y)={13y},      H(y)={y+1/2},
S(y)={13y+1/2}.                                           (10)
```

Let `W` be THM-2693's raw delayed carrier

```text
D_(13^3) intersection D_14^c intersection D_27^c
 intersection D_40^c intersection D_53^c intersection D_66^c
 intersection D_(2*13^5)^c.                              (11)
```

Both

```text
y_-=11/24,                         y_+=13/24               (12)
```

are fixed by `S`.  Their displayed factor distances are

```text
(1/24,5/12,3/8,1/3,7/24,1/4,1/12),                      (13)
```

so they lie strictly in `W`.  Their `T` midpoints are `23/24` and `1/24`;
the target distance there is `11/24`, so both midpoints lie strictly outside
`W`.  Since `S^2=T^2`, each fixed `S` loop gives a diagonal arrow in
`Def_2(W,T)`.

The defect is open.  The first endpoint inequality to bind under `T^2` is
the high-speed safety factor.  Its phase margin is

```text
1/12-1/14=1/84,
```

its speed is `2*13^5=742586`, and the endpoint derivative is `13^2=169`.
Hence for either sign in `(12)`, every

```text
|y-y_+/-| < r,                  r=1/10541750856            (14)
```

still satisfies `y,T^2y in W` and `Ty notin W`.  The two intervals are
disjoint, so the source-coordinate defect has

```text
mu(SrcDef_2(W,T)) >= 4r = 1/2635437714.                   (15)
```

This is a raw delayed-coordinate statement.  THM-2698 separately provides a
fully decorated positive half-odometer cycle and its much smaller physical
packet cylinder.

## 4. A nonconstant ordinary-lattice endpoint-sampled macro cycle

There is a stronger positive control which does not descend from the
half-odometer fixed loop.  Put `R=13^6`, let `D(x)={13x}`, and use `A_a`
only for pure translation:

```text
A_a(x)={x+a/R},                 M_K=A_K after D^2.         (16)
```

The two strict source-one THM-2640 packet points are

```text
x0=39123022/82055753,       pi(x0)=4/17,
(j,shallow,owner,carry,h,kappa,edge,root)
  =(8,1,4,3,3,0,left,7),

x1=41305372/82055753,       pi(x1)=13/17,
(j,shallow,owner,carry,h,kappa,edge,root)
  =(3,4,1,1,9,1,left,4).                                  (17)
```

Here `pi(x)={Rx}`.  Every geometric rail, dynamically typed present factor,
delayed word, predecessor-carry cell, future-half-digit cell, private
half-tooth, and root condition is strict; the retained private row is
primitive/unit.  Owner-to-next-shallow gluing is literal.

For

```text
K0=4472391,                    K1=1956127                  (18)
```

exact rational arithmetic gives

```text
M_(K0)(x0)=x1,                 M_(K1)(x1)=x0.              (19)
```

The bare `D^2` roots are respectively `2` and `12`; the root increments are
`2K0=2` and `2K1=8` modulo thirteen, landing at roots `4` and `7`.
Both packet slacks are

```text
11/853068347561612.                                      (20)
```

Therefore the endpoint-sampled three-state macro cylinder about `x0` has
exact radius and length

```text
rho=11/24364485074707200332,
2rho=11/12182242537353600166.                             (21)
```

Only the center is an exact two-cycle.  The open cylinder certifies the
three-event finite horizon `x,M_(K0)x,M_(K1)M_(K0)x` in the same packet
cells; perturbed points need not return to themselves.

On the delayed base,

```text
4/17 --T^2--> 13/17 --T^2--> 4/17,                       (22)
```

and both endpoints lie strictly in `W`: the seven factor distances are

```text
(1/17,5/17,6/17,7/17,8/17,8/17,2/17).                   (23)
```

The forced unit-step midpoints are

```text
T(4/17)=1/17,                    T(13/17)=16/17.           (24)
```

At both, the target distance is `4/17>1/14`, so they lie strictly outside
`W`.  Thus `(19)` is a positive, nonconstant ordinary-lattice
endpoint-sampled affine macro cycle: both endpoint charts are strict and
their owner/next-shallow labels agree.  The full affine edges use two
different maps `M_(K0),M_(K1)`, so they are not themselves arrows of the
single-map category `C_W(T)`.  Rather, they lift the two strict delayed-base
relation defects

```text
(4/17,T^2(4/17)),                 (13/17,T^2(13/17))       (24a)
```

in `Def_2(W,T)`.  This is not yet a physical semantic/current cospan: the
factors between the displayed endpoint charts have not been assembled into
one lawfully co-shifted table.

## 5. Every integral affine subdivision fails at the same midpoint

For `a in Z/RZ`, put

```text
F_a=A_a after D,                 F_a(x)={13x+a/R}.         (25)
```

Every integral degree-one factorization of `M_K` has the form

```text
M_K=F_b after F_a,                   13a+b=K mod R.        (26)
```

But

```text
pi(F_a x)={13 pi(x)}=T(pi(x)),                             (27)
```

independently of `a`.  Therefore no choice of fibre lifts in `(26)` can move
the middle phases in `(24)` back into `W`.  The obstruction is universal
over all integral affine factorizations of the displayed macro edges, not a
failure of the selected `K0,K1` presentation.

## 6. Minimal midpoint completion and the live bridge

For arbitrary `(W,T)`, define the required midpoint carrier

```text
Mid_2(W,T)=T(W intersection T^(-2)W),
Ghost_2(W,T)=Mid_2(W,T)\W.                                (28)
```

After enlarging the object set from `W` to `W union Mid_2(W,T)`, every old
endpoint macro in `E_2^W(T)` factors through `Mid_2(W,T)`.  This set is
minimal: any proposed intermediate object set through which every old
degree-two endpoint arrow factors by the two fixed degree-one `T` legs in the
ambient action category must contain the unique points `Tx` in `(28)`.  This
is not a minimality assertion among arbitrary cospans or alternate leg maps.
The factorization does not occur inside
`C_W(T)` when `Tx notin W`.  In the exact controls above, `Ghost_2` contains

```text
23/24, 1/24, 1/17, 16/17.                                (29)
```

The physical ghost audit is more informative than the set-theoretic
completion.  At delayed phases `1/17` and `16/17`, respectively,

```text
(h,kappa)=(0,1),                         (12,0),            (30)
```

both inherited delayed sectors are empty and every inherited pair-prefix
clock is empty in both sectors.  Hence neither midpoint carries an old
THM-2640 aggregate unit row.  Reusing that row after deleting the failed word
would be an illegal change of universe, and guard-danger rows realized at
other points cannot be attached pointwise.  This does not rule out rebuilding
a genuinely new private transit row with a new normalization; such a rebuild
is a separate changed-universe theorem and is not asserted here.

If one honestly omits the old delayed word and the old aggregate unit-row
predicate, the remaining pointwise source-one atlas is nevertheless large.
For the forward midpoint, shallow/owner clock `4` and rails `8,9` give the
exact cascade

```text
371293 -> 53042 -> 8160 -> 2925 -> 2925 -> 5850,           (31)
```

and all `5850` final cells obey both integral affine carry/root split laws.
The first is

```text
x=41513423/82055753, N=2441966, j=9, carry=7,
(h,kappa,edge,root)=(0,1,left,2),
(a,b)=(1485215,4471832), source/middle roots=(7),(2).      (32)
```

For the reverse midpoint, shallow/owner clock `1` and rails `2,3` give

```text
371293 -> 53042 -> 8160 -> 2479 -> 2479 -> 4958,           (33)
```

and again all `4958` final cells obey both split laws.  Its first cell has

```text
x=38392136/82055753, N=2258360, j=2, carry=0,
(h,kappa,edge,root)=(12,0,left,2),
(a,b)=(4459563,1903516), source/middle roots=(7),(12).     (34)
```

There is also a sharp positive hostile to simply declaring these new transit
objects to be ordinary same-row events.  On the fixed canonical nine-speed
control

```text
(1,14,27,40,53,66,13,13^3,2*13^5),                       (35)
```

all `2925` forward rail points (both edges, hence `5850` final cells) and
`2193` of the `2479` reverse rail points (both edges, hence `4386` final
cells) are strictly safe for every displayed speed.  The first safe forward
cell is `(32)`, with x-radius slack `45/2297561084`; the first safe reverse
cell is

```text
x=38400313/82055753, N=2258841, j=2, carry=0,
(h,kappa,edge,root)=(12,0,left,2),
(a,b)=(4460044,1897263), source/middle roots=(7),(12),     (36)
```

with x-radius slack `227/2297561084`.  Thus both macro edges have strict
all-safe canonical-row midpoints.  A new transit grammar must add event
semantics, not merely factor safety.  This is a fixed-control safe-exit
signal, not an identification with every live profile and not a ledger
decrement.

There is a compact integral holotopy shadow.  Along the completed delayed
four-cycle

```text
4/17 -> 1/17 -> 13/17 -> 16/17 -> 4/17,                  (37)
```

the strict danger signatures in the canonical nine-coordinate control are

```text
{c1,c2}, {guard}, {c1,c2}, {guard}.                       (38)
```

Write `E` for an endpoint-type vertex and `G` for a ghost-type vertex.  The
two endpoint invoices land on the two ghost phases after one and three base
steps.  Hence on the free abelian type quotient their incidence is

```text
(B_*+B_*^3)[E]=2[G].                                     (39)
```

Before quotienting by the two phase symmetries this is the all-ones `2 x 2`
incidence matrix.  On the invariant rank-one quotient its cokernel is
`Z/2`.  Reducing modulo two therefore erases the load-bearing positive
multiplicity; an integral/nonnegative completion needs two labelled arms or
a lawful cancellation.  Equation `(39)` is an incidence invariant, not a
semantic current.

These facts isolate four honest next architectures:

1. physicalize the ghost layer as a labelled transit object and construct
   a **new** transit grammar with its private unit rows recomputed, rather
   than importing the now-empty old word;
2. construct a lawful direct `D^2` endpoint/current cospan which may
   marginalize the forbidden midpoint; or
3. on the central half branch, construct the principal semantic `C_2`
   bibundle and actual odd-edge intertwiners, retaining the `S` subdivision
   before taking the even composite; or
4. prove that the fixed-control safe exits in `(35)--(36)` transfer to the
   actual live-row typing, which would turn the ghost obstruction itself into
   a closure mechanism.

THM-2315 warns that a bare one-skeleton or 2-Segal summary cannot reconstruct
the intermediate pullback.  THM-2292 and THM-2658 give the exact
counterindication: a common calibration or complete balanced component-gain
section would produce an actual middle witness, but endpoint pairs alone are
an incomplete graph and the old `W` middle is empty.  THM-2680 exhibits the
dual failure direction (a quotient one-skeleton has walks while the physical
two-step path is empty).  THM-2644 cannot consume the macro loop: it still
needs one nonnegative transition, a common gauge, and a lawful common middle
for both quadratic compositions.

THM-2701 is nonoverlapping.  It proves that the literal singleton-word
one-step language on its canonical row dies at length six and that only an
enlarged edge-debt grammar recurs.  Here the old middle word is already empty,
while a degree-two endpoint relation survives after skipping that forbidden
middle.  No chain map from this endpoint relation to THM-2701's terminal-word
or edge-debt grammar is supplied.

It also gives a no-terminalization check on `(37)`.  With
`Q=Q_a union Q_b`, THM-2701 says the five-arrow literal nerve stratum
`N_5^Q(B)` is empty.  A degree- and chronology-preserving functor from the
recurrent completed four-edge endpoint/ghost path to `C_Q(B)` would iterate
to such a five-arrow path, a contradiction.  Thus a physicalized ghost must
remain a genuinely relative/coloured object, retain the debt, or change the
chronology; it cannot be an old terminal-word object with the debt forgotten.

No semantic endpoint current, direct macro cospan, principal `C_2` sidecar,
row exclusion, or LRC(14) conclusion is proved.

## 7. Compact exact referee

Run normally and with `-O`:

```text
python 04-computation/lrc14_relative_segal_macro_cycle_thm2706.py
python -O 04-computation/lrc14_relative_segal_macro_cycle_thm2706.py
```

The two modes byte-match.  The primary referee checks `(7)--(15)`, both physical
packet rows in `(17)`, `(18)--(24)`, primitive units, roots, all packet
slacks, and the exact macro-cylinder length.  The universal assertion
`(26)--(27)` is the displayed one-line affine calculation, not an exhaustive
search over `R^2` presentations.

The secondary exact referee

```text
python 04-computation/lrc14_relative_segal_ghost_preword_transit_thm2706.py
python -O 04-computation/lrc14_relative_segal_ghost_preword_transit_thm2706.py
```

checks the empty old profiles, every cell in `(31)` and `(33)`, every affine
carry/root split, and every strict fixed-control safety statement in
`(35)--(36)`.  It deliberately does not manufacture a replacement transit
word or aggregate unit row.
