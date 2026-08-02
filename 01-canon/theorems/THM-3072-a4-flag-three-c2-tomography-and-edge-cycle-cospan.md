---
id: THM-3072
title: "A4 flag three-C2 tomography and edge-cycle cospan"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  The regular twelve-flag A4 torsor has three conjugate C2
  quotients.  One is the six edges of K4, one is the six oriented
  Hamiltonian cycles, and the third is their missing conjugate chart.  The
  edge/cycle joint label is pointwise injective, with image three disjoint
  K2,2 components, but separate quotient tables have only rank nine and
  lose an exact three-dimensional V4-character sector.  All three quotient
  averages reconstruct by P0+P1+P2=I+2P_V4.  This is finite permutation-
  module tomography, not a physical quartic, Keller, Farey, or LRC
  intertwiner.
source: root-2026-08-01-a4-three-quotient-tomography
depends_on:
  - THM-3067-tetrahedral-modular-two-three-flag-quotient-and-origin-loss
related:
  - THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame
  - THM-2968-quartic-edge-and-oriented-cycle-s4-complements
  - THM-2971-discriminant-cover-edge-orientation-sextic-algebra-intertwiner
  - THM-2975-modular-six-sheet-schreier-graphs-and-farey-partial-cube-boundary
  - THM-3049-k4-matching-monomial-tropical-root-extraction-clutch
script: 04-computation/modular_a4_flag_three_quotient_tomography_thm3072.py
output: 05-knowledge/results/modular_a4_flag_three_quotient_tomography_thm3072.out
script_sha256: 856cab6338ef68560ef0b1ef5cef164a771a0a2b671a44abb0f8dec7bc75fb76
output_sha256: e1605bb1c1ddec2d18209ab3ea003738b1b2282317d3728fd8711c735d73b18e
hash_basis: LF-normalized bytes
---

# THM-3072 -- A4 flag three-C2 tomography and edge-cycle cospan

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

## 1. Three binary quotients of one tetrahedral flag torsor

Let

```text
V=F_2^2,                 D=V\{0},                       (1)
```

and choose a linear map `rho` of order three on `V`, so `rho` cyclically
orders `D`.  Use the twelve flags from
[THM-3067](THM-3067-tetrahedral-modular-two-three-flag-quotient-and-origin-loss.md):

```text
Omega=V x D.                                               (2)
```

For `j in Z/3` define

```text
tau_j(x,d)=(x+rho^j d,d).                                  (3)
```

Since over `F_2`

```text
d+rho d=rho^2 d,                                           (4)
```

the three `tau_j` are the three nonidentity elements of a right `V_4`:

```text
tau_j^2=1,             tau_0 tau_1=tau_2.                  (5)
```

They commute with the left affine `A_4=V_4 semidirect C_3` action.  Put

```text
X_j=Omega/<tau_j>.                                         (6)
```

Every `X_j` is therefore a six-sheet transitive `A_4`-set of type
`A_4/C_2`.  The right ternary move

```text
R(x,d)=(x,rho d)                                           (7)
```

cyclically conjugates the three deck subgroups and therefore transports the
three quotient charts.  It does not act on any one fixed chart.  This is the
exact repaired meaning of the no-descent hostile in THM-3067.

## 2. Edges, oriented cycles, and the third chart

The first quotient is the ordinary edge set of the affine tetrahedron:

```text
q_0(x,d)={x,x+d}.                                          (8)
```

Its deck involution is `tau_0`.

Let an oriented Hamiltonian cycle be taken modulo cyclic rotation but not
reversal.  Define

```text
q_2(x,d)
 =[x, x+d, x+d+rho d, x+rho d].                           (9)
```

Translation by `rho^2 d=d+rho d` rotates this cycle by two places.  Exact
comparison of fibres gives

```text
q_2(omega)=q_2(omega')
 iff omega'=omega or tau_2 omega.                         (10)
```

Thus `X_2` is exactly the six oriented Hamiltonian four-cycles used by
[THM-2968](THM-2968-quartic-edge-and-oriented-cycle-s4-complements.md).

The remaining quotient is another edge chart:

```text
q_1(x,d)={x,x+rho d},                                     (11)
```

with deck involution `tau_1`.  It is isomorphic to the edge set, but it is
not the same labelled quotient of `Omega`.

The odd reflection `J` of THM-3067 satisfies

```text
J tau_0 J=tau_0,            J tau_1 J=tau_2.              (12)
```

So it fixes the ordinary edge chart and swaps the other edge chart with the
oriented-cycle chart.  This binary operation is separate from the ternary
cycle `(7)`.

## 3. The common-atom cospan

The pointwise map

```text
(q_0,q_2): Omega -> X_0 x X_2                             (13)
```

is injective.  Indeed, equality in the first coordinate puts two flags in
one `tau_0` orbit; equality in the second puts them in one `tau_2` orbit;
the two order-two deck groups intersect trivially.

Its image has a particularly rigid form.  For each `d in D`, let

```text
M_d={{x,x+d}:x in V}                                      (14)
```

be the two-edge perfect matching, and let `C_d` be the two orientations in
the image of `(9)` at fixed `d`.  Then

```text
image(q_0,q_2)=disjoint_union_(d in D) M_d x C_d.         (15)
```

Thus the bipartite incidence graph of the image is exactly three disjoint
copies of `K_(2,2)=C_4`.  These three components are indexed by the three
perfect matchings which form the positive `K_4` carrier in
[THM-3049](THM-3049-k4-matching-monomial-tropical-root-extraction-clutch.md).

Equation `(13)` is a **joint-label** result.  If an actual atom carries both
labels, it recovers its flag.  It says nothing yet about recovering a signal
from two separately aggregated quotient tables.  That distinction is the
load-bearing point of the next section.

## 4. Exact three-view tomography

Let `k` be a field of characteristic different from two and put

```text
W=k^Omega.                                                (16)
```

Let `tau_j` also denote pullback on `W`, and define the quotient averages

```text
P_j=(I+tau_j)/2,                                          (17)
P_V=(I+tau_0+tau_1+tau_2)/4.                              (18)
```

Then `im(P_j)` is precisely the six-dimensional space of functions pulled
back from `X_j`, while `im(P_V)` is the three-dimensional direction
quotient.  The commuting `V_4` algebra gives, for distinct `i,j`,

```text
P_i P_j=P_V,                                              (19)
P_0+P_1+P_2=I+2P_V.                                      (20)
```

Consequently every `f in W` is reconstructed from the three quotient
tables by

```text
f=P_0 f+P_1 f+P_2 f-2P_V f,                              (21)
```

and `P_V f` is already computable as `P_i P_j f` from any two displayed
tables.  Thus no fourth table is hidden in `(21)`.

The dimension ledger is sharp:

```text
dim im(P_j)=6,
dim im(P_V)=3,
dim(im(P_i)+im(P_j))=9,
rank[f |-> (P_i f,P_j f)]=9,
rank[f |-> (P_0 f,P_1 f,P_2 f)]=12.                      (22)
```

Equivalently, the regular module splits into the four right-`V_4` character
sectors

```text
W=W^V direct_sum W_0 direct_sum W_1 direct_sum W_2,
dim W^V=dim W_j=3,                                       (23)
```

where `W_j` is the nontrivial character sector whose kernel is
`<tau_j>`.  Hence

```text
im(P_j)=W^V direct_sum W_j.                              (24)
```

Any two quotient averages lose exactly the third nontrivial sector.

For the physically suggestive pair `X_0` (edges) and `X_2` (oriented
cycles), this blind sector is `W_1`.  An explicit integral basis is obtained
one direction at a time.  Write uniquely

```text
x=a d+b rho d,              a,b in F_2,                  (25)
```

and, for fixed `d_0`, set

```text
f_(d_0)(x,d)=0                         if d!=d_0,
f_(d_0)(x,d)=(-1)^a                    if d=d_0.           (26)
```

Then the three `f_(d_0)` are independent and

```text
P_0 f_(d_0)=P_2 f_(d_0)=0,
P_1 f_(d_0)=f_(d_0).                                    (27)
```

So separate edge and oriented-cycle marginals do **not** recover a general
flag signal, even though the joint point-label map `(13)` is injective.  The
missing coordinate is exactly a cross-view coupling, not another point
census.

## 5. Relation to the quartic and modular frontiers

The three-dimensional common sector `W^V` is the same information type as
the cubic-resolvent `S_3` quotient in
[THM-2950](THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame.md):
it remembers the three matching directions and forgets affine position.  The
three additional sectors in `(23)` are the three conjugate
`V_4`-character channels above that quotient after a cyclic orientation is
chosen.  They are not the two `S_4` descent classes classified by THM-2975.

Over the quartic discriminant cover,
[THM-2971](THM-2971-discriminant-cover-edge-orientation-sextic-algebra-intertwiner.md)
proves that the edge and oriented-cycle sextic **algebras** become
isomorphic.  Equations `(8)--(15)` give a different object: a common regular
set-theoretic double cover of their six-sheet `A_4` actions.  No polynomial
map from a quartic source to `Omega`, no multiplication-preserving rank-
twelve algebra, and no common physical coefficient is constructed here.
The exact algebraic sidecar remains THM-2971's discriminant orientation and
primitive-coordinate map.

[THM-2975](THM-2975-modular-six-sheet-schreier-graphs-and-farey-partial-cube-boundary.md)
shows that the marked six-sheet modular actions have short finite kernels
and lose the Farey word and Gram owner.  The present theorem does not repair
those losses.  It instead explains why the ternary move cannot live on a
single six-sheet chart: it rotates the **three charts** in `(6)`.  A Farey or
LRC transplant would still need a common source atom carrying the three
quotient labels together with the load-bearing Gram/owner/current data.

No SFC(4), degree-four point-cap Keller, `JC(2)`, `DC(2)`, GMC, LRC, N-body,
or supergroup conclusion follows.

## 6. Exact evidence

Run

```text
python 04-computation/modular_a4_flag_three_quotient_tomography_thm3072.py
python -O 04-computation/modular_a4_flag_three_quotient_tomography_thm3072.py
```

Both modes byte-match the LF-stored transcript

```text
05-knowledge/results/modular_a4_flag_three_quotient_tomography_thm3072.out.
```

The dependency-free companion checks all twelve flags; all three deck
relations and six-sheet quotients; the edge, skew-edge, and oriented-cycle
fibre laws; left-`A_4` equivariance; the three `K_(2,2)` joint components;
the rational projector identities; every rank in `(22)`; the explicit blind
basis `(26)--(27)`; and the odd-reflection action `(12)`.  Every truth-bearing
gate uses an explicit exception rather than Python `assert`, so optimized
mode verifies the same claims.

LF-normalized SHA256:

```text
script  856cab6338ef68560ef0b1ef5cef164a771a0a2b671a44abb0f8dec7bc75fb76
output  e1605bb1c1ddec2d18209ab3ea003738b1b2282317d3728fd8711c735d73b18e
```

**QED, conditional only on independent hostile audit for promotion.**
