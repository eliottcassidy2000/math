---
id: THM-2889
title: "Dicyclic reverse-action joint carrier and skew-lift separation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Selector
  provenance and carry-edge reversal are the two independent semantic
  characters chi_s=det(QB,-) and chi_c=det(QA,-).  No single
  Q8-to-C2 character realizes both cyclic inversion actions, and the
  literal seam event pattern (0,1,0) is not a character.  Each separately
  gives a differently labelled Dic_338=C169 semidirect Q8 shadow of
  order 1352.
  The carry shadow has the exact via/direct clutch C^13*(-1), of order
  26.  The stable-diagonal coordinatewise-inversion envelope realizing
  both as full C169 actions is B=(C169_s x C169_c) semidirect Q8 of
  order 228488, with center C2,
  derived subgroup C169^2 x C2, abelianization V4, and minimum faithful
  splitting-field degree four.  Its balanced torus element acts by
  diag(E,1,1,E^-1), exactly the raw-Prony Hermitian weights.  This
  identifies only the torus diagonal: Hermitianization kills the
  quaternionic center/global phase, the original Prony two-plane is not
  Q8-stable, and no physical four-channel operator, unmarked current
  descent, row exclusion, or LRC(14) proof follows.
audit: >
  Two independent referees reconstructed the selector/carry character
  split, both labelled Dic_338 presentations, the order-26 clutch, the
  order and conjugacy censuses, the bi-dicyclic center/derived/quotient
  anatomy, the four-weight minimal-degree argument, the typed right horn
  action, and the q7 filler boundary.  They independently specialized the
  two- and four-dimensional matrices in both endpoint fields and checked
  the balanced Hermitian diagonal, the central-sign loss, and the
  noninvariant Prony plane.  Normal, optimized, and stored tracked outputs
  agree; the companion contains no Python assert.  The referees accepted
  the abstract-versus-physical, class-function, scalar-character,
  Hermitianization, row-exclusion, and LRC boundaries.  Verdict: ACCEPT.
source: root/lrc-dicyclic-reverse-action-2026-07-29
depends_on:
  - THM-2868-two-origin-endpoint-projector-and-projective-kummer-carry-descent
  - THM-2882-event-twisted-all-q-coefficient-carry-lift
  - THM-2884-macro-semantic-diagonal-horn-carrier-and-origin-even-boundary
  - THM-2886-stepped-origin-provenance-transport-on-the-v4-horn
  - THM-2887-quaternionic-arf-lift-of-the-semantic-v4-and-global-carry-no-go
related:
  - THM-2851-q3-q11-q7-bockstein-holonomy-and-realization-sidecar
  - THM-2874-endpoint-kummer-galois-clutch-and-bockstein-seam-transgression
  - THM-2878-endpoint-factor-exit-carry-transducer-and-flat-vertex-boundary
script: 04-computation/lrc14_dicyclic_reverse_action_joint_carrier_thm2889.py
output: 05-knowledge/results/lrc14_dicyclic_reverse_action_joint_carrier_thm2889.out
script_sha256: 1a2a664253c4a38591b859c9391896084f129a43b6b274b03d706bdb6c4a77c5
output_sha256: a6674dbc4db9c80cc0151a6f014962fa8b2d09930b64918ce6a1c8b005217942
lean_module: 04-computation/lean/TournamentH7/TournamentH7/LRCC169CarrySelector.lean
lean_sha256: 4027012d46f70fce92d1210f70204f4b55e27e0204affee3b8ba6d6622132d05
hash_basis: LF-normalized bytes
---

# THM-2889 -- dicyclic reverse-action joint carrier and skew-lift separation

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2887 ruled out the semidirect direction

```text
C169 -> Aut(Q8).                                        (1)
```

It did not rule out the reverse direction.  The reverse action exists,
but a label audit exposes a sharper fact: selector provenance and
carry-edge reversal require two different characters of the semantic
`V4`.  One cyclic inversion action can use either character, not both.
A single dicyclic group still retains the whole `Q8`, so it can evaluate
the other character as an observable; what it cannot do is realize both
as the same cyclic inversion character.

The minimal stable envelope realizing both as full `C169`
coordinatewise-inversion actions on the chosen diagonal therefore has
two cyclic coordinates and four linear channels.  Its balanced torus
action is exactly the Hermitianized Prony transport, but that
Hermitianization forgets the quaternionic center.  This theorem
identifies that scoped abstract carrier and the precise physical debt
that remains.

## 1. Two characters, not one

Use the THM-2887 coordinates

```text
QA=(1,0),             QB=(0,1),             QAB=(1,1)
```

on `V4=Q8/{+/-1}`.  THM-2886's selector-pair XOR parity on the target
path `(Q0,QA,QAB)` is

```text
(0,1,1).                                                (2)
```

There is a unique linear character with these values:

```text
chi_s(h)=det(QB,h).
```

Its values and kernel are

```text
chi_s(QA,QB,QAB)=(1,0,1),
ker(chi_s in Q8)=<QB>={1,QB,-1,-QB}.                    (3)
```

In particular the successive `q11 -> q7` semantic direction is `QB`,
which `chi_s` fixes.  A character which instead keeps `QA` flat and
reverses `QB` is

```text
chi_c(h)=det(QA,h),
chi_c(QA,QB,QAB)=(0,1,1),
ker(chi_c in Q8)=<QA>.                                  (4)
```

On `(Q0,QA,QAB)`, however, `chi_c` has values `(0,0,1)`, not `(2)`.
Therefore:

```text
no Q8 -> C2 character both
  recovers the selector parity (0,1,1)
  and reverses the QB carry edge.                       (5)
```

There is a second, independent failure.  The literal reduced carry-event
pattern on the two successive edges and their composite is

```text
(q3->q11, q11->q7, q3->q7)=(0,1,0).                    (6)
```

Because `QA+QB=QAB`, no character has the values `(6)`.  The event is
Bockstein curvature, not a vertex action.  These two finite no-go
statements are also sorry-free theorems

```text
no_single_character_types_selector_and_carry
carry_event_pattern_not_character
```

in the root-imported Lean module `LRCC169CarrySelector`.

## 2. The two labelled dicyclic shadows

Let `C=<c> ~= C169`.  Each character defines a reverse semidirect
product

```text
G_s=C semidirect_(chi_s) Q8,
G_c=C semidirect_(chi_c) Q8.                            (7)
```

The groups are abstractly isomorphic but differently labelled.  Both
are

```text
Dic_338 =
  <a,x | a^676=1, x^2=a^338, x a x^-1=a^-1>,           (8)
```

using the convention `|Dic_n|=4n`.

For the selector shadow one may take

```text
c=a^344,      <c>=<a^4>,
QA=x,         QB=a^169,       QAB_section=QA*QB.        (9)
```

Then `QA,QAB` invert `c` while `QB` fixes it.  Swapping the labels
`QA,QB` gives the carry shadow, in which `QB,QAB` invert `c` and `QA`
fixes it.

The section convention matters.  The computational normal form in `(9)`
uses `QAB_section=QA*QB`; the canonical THM-2887 section uses

```text
QA*QB=-QAB_canonical.                                  (10)
```

The difference is precisely the central sign.

For either shadow,

```text
|C|=169,       |Q8|=8,       C intersect Q8={1},
|G_s|=|G_c|=1352.                                      (11)
```

The exact element-order census is

```text
1^1, 2^1, 4^678, 13^12, 26^12, 52^24,
169^156, 338^156, 676^312.                             (12)
```

There are `341` conjugacy classes, with `(order,class size)^count`

```text
(1,1)^1, (2,1)^1, (4,2)^1, (4,338)^2,
(13,2)^6, (26,2)^6, (52,2)^12,
(169,2)^78, (338,2)^78, (676,2)^156.                  (13)
```

Moreover

```text
Z(G)={+/-1},
[G,G]=<a^2> ~= C338,
G_ab ~= V4.                                            (14)
```

Hence every one-dimensional character kills the `C169` subgroup.

## 3. The carry-shadow order-26 clutch

In the carry labelling, `QA` fixes `c`, `QB` inverts it, and the canonical
section satisfies `(10)`.  Exact multiplication gives

```text
(c^8 QA)(c^9 QB)
  =(c^13*(-1))(c^4 QAB_canonical).                     (15)
```

Thus the via/direct defect is the single element

```text
K=c^13*(-1),
|K|=26,
QB K QB^-1=K^-1,
K^2=c^26.                                              (16)
```

The vertical base-thirteen carry and the quaternionic central sign have
therefore combined into one exact order-`26` clutch.  This is the clean
positive result in the carry shadow.  It does not recover selector
parity, since its target values are `(0,0,1)`.

## 4. The coordinatewise-inversion joint carrier

To realize both characters as inversion actions on full `C169`
coordinates, retain them on separate cyclic coordinates:

```text
B=(C169_s x C169_c) semidirect Q8,                     (17)

QA:(s,c)->(-s,c),
QB:(s,c)->(s,-c),
QAB:(s,c)->(-s,-c).                                    (18)
```

Equivalently,

```text
B ~= G_s x_(Q8) G_c,
|B|=1352^2/8=228488.                                   (19)
```

Its exact anatomy is

```text
Z(B)={+/-1},
[B,B]=C169^2 x {+/-1},             order 57122,
B_ab=V4.                                               (20)
```

To avoid dihedral notation ambiguity, the central quotient is

```text
B/<-1> ~=
  (C169 semidirect C2) x (C169 semidirect C2),          (21)
```

where each factor has order `338`.

The full order census is

```text
1^1, 2^1, 4^57798, 13^168, 26^168, 52^8112,
169^28392, 338^28392, 676^105456.                      (22)
```

The `14789` conjugacy classes have census

```text
(1,1)^1, (2,1)^1,
(4,338)^2, (4,57122)^1,
(13,2)^12, (13,4)^36,
(26,2)^12, (26,4)^36,
(52,338)^24,
(169,2)^156, (169,4)^7020,
(338,2)^156, (338,4)^7020,
(676,338)^312.                                        (23)
```

### Minimal full-C169 coordinate rank and linear degree

Within direct sums of full `C169` coordinates with coordinatewise sign
actions, one coordinate can carry only one nonzero `V4 -> C2` character.
Requiring `chi_s,chi_c` as separate inversion actions therefore makes
the full-`C169` coordinate rank two minimal.  This does not exclude a
smaller auxiliary odd cyclic coordinate outside the full-address model,
nor claim that `B` is minimal among all groups retaining `Q8` and both
characters as observables.

On the dual of `C169^2`, the sign action has orbit census

```text
size 1: 1,             size 2: 168,             size 4: 7056.  (24)
```

One axis orbit has size two but spans only one cyclic direction.  A
generic orbit, such as the orbit of `(1,1)`, has size four and spans both
directions because `2` is a unit modulo `169`.  Hence a faithful
semisimple splitting-field representation has degree at least four, and
degree four is sharp.

The inherited physical address initially supplies only the diagonal
`{(n,n)}`.  Its `Q8` orbit contains `(1,1),(-1,1),(1,-1),(-1,-1)` and
generates all `C169^2`.  Thus `B` is the minimal stable abelian envelope
of that diagonal **inside this abstract two-character model**.  The
off-diagonal states have not been physically constructed.

## 5. Typed right horn action and filler separation

The group law must retain its action side.  With right multiplication in
the section gauge `QA*QB=QAB_section`, define

```text
e8=( 8, 8,QA),
e9=(-9, 9,QB).                                         (25)
```

The minus sign in the selector coordinate is forced because the
intermediate `QA` state inverts that coordinate.  For each ancestry
`a in C13`, the diagonal `q3` state follows

```text
(13a+3,13a+3,Q0)
 -> (13a+11,13a+11,QA)
 -> (13a+20,13a+20,QAB_section).                       (26)
```

The canonical direct edge is `(4,4,QAB_canonical)`, and

```text
e8*e9=(13,13,-1)*(4,4,QAB_canonical).                  (27)
```

All thirteen ancestry starts satisfy `(26)--(27)`.  Replacing `-9` by
`+9`, or forgetting whether the action is on the left or right, is false.

The two q7 fillers lift to

```text
(13a+7,13a+7,QAB),
(13(a+1)+7,13(a+1)+7,QAB).                             (28)
```

They are distinct, square to the central sign, and are conjugate by
`(91,91)`, since `2*91=13 mod169`.  This is pointwise separation, not
invariant separation: the fillers are conjugate and every class function
still identifies them.  A physical consumer must retain the marked lift.

## 6. Four-channel Pauli carrier

Over either endpoint field, choose a primitive `169`-th root `s` and
`i^2=-1`.  Put

```text
D(s)=diag(s,s^-1),
X=[[0,1],[1,0]],
Z=diag(1,-1),

c_s=D(s) tensor I,            c_c=I tensor D(s),
QA =i X tensor I,             QB =i Z tensor X.        (29)
```

Then

```text
QA^2=QB^2=-I,          QA*QB=-QB*QA,
QA inverts c_s and fixes c_c,
QB fixes c_s and inverts c_c.                          (30)
```

The `169^2` torus matrices are distinct, do not contain `-I`, and the
four semantic directions have distinct monomial permutation supports.
Thus `(29)` is faithful.  The weight proof above makes its degree four
minimal.

The tensor factors suggest

```text
selector/origin channel tensor Prony/carry channel,
```

but this is an interpretation of the abstract representation, not a
constructed physical operator.

## 7. Prony and Hermitian boundary

Let `xi` be the inherited primitive `2366`-th root and set

```text
omega=xi^182,              rho=xi^42.                  (31)
```

Then `rho` has order `169` and `rho^13=omega^3`.  The normalized THM-2868
projective coordinate is

```text
t_r=U_r/V_r=xi^(955+546r),
t_(r+1)=omega^3 t_r.                                   (32)
```

Cyclotomic inversion gives `bar(t_r)=t_r^-1`.  In the selector shadow,

```text
c:t->rho t,        QA:t->-bar(t),
QB:t->-t,          QAB:t->bar(t).                      (33)
```

The carry shadow swaps the `QA,QB` roles.  The THM-2882 lifted line also
admits the exact state-dependent relative-channel gauge

```text
gamma(a,q)=rho^(13a+q) omega^(-3(a+q-3)),              (34)
```

for which all `13^3` transitions in each field satisfy

```text
gamma(L_h(a,q))*omega^(3(h+kappa(q,h)))
  =rho^h*gamma(a,q).                                   (35)
```

This is not a common scalar gauge on `U+V`.

### Raw frame and scalar obstruction

On the wrapped edge, with `E=omega^4`,

```text
diag(E,1)
  =omega^2 diag(omega^2,omega^-2).                     (36)
```

The determinant-one factor is the dicyclic rotation `c^130`.  In the
gauge-normalized honest `C169` frame, however, displacement `+9` is
`c^9`; the two frames must not be identified.

The common scalar `omega^2` cannot be supplied by a one-dimensional
character.  Both cyclic circles lie in `[B,B]`, so every such character
kills them.  This is the group-theoretic form of the THM-2886 failure of
scalar descent through `U+V`.

### Balanced Hermitianization

The exact four-channel torus element is

```text
c_s^130 c_c^130 = diag(E,1,1,E^-1).                    (37)
```

The same diagonal is the conjugation action of the raw Prony matrix
`diag(E,1)` on the Hermitian rank-one coordinates

```text
(U*bar(V), |U|^2, |V|^2, V*bar(U)).                    (38)
```

This is an equality of the **balanced torus diagonals only**.  The
faithful representation `(29)` sends the quaternionic center to `-I4`,
whereas Hermitianization kills every global phase and sends that center
to the identity.  Moreover `QA` sends the literal Prony plane spanned by
the first two tensor coordinates to the complementary plane; the
original two-channel space is not `Q8`-stable.

Thus Hermitianization is the correct four-channel linear completion of
the relative phase, but it is not a faithful realization of the
quaternionic sign and does not itself produce the physical current.

## 8. Verification and exact scope

The companion pins THM-2868, THM-2882, THM-2884, THM-2886, and THM-2887.
It checks:

- all four semantic characters and both no-go patterns;
- the complete `1352`-element dicyclic and `228488`-element joint order
  censuses;
- all `341` and `14789` conjugacy classes;
- center, derived, abelianization, quotient, and dual-weight anatomy;
- all thirteen right-action and filler rows;
- faithful two- and four-dimensional specializations in both fields;
- `2197` gauge transitions per field; and
- the balanced Hermitian diagonal, central loss, and Prony-plane exit.

Normal, optimized, and stored outputs agree byte for byte.  Two
independent hostile audits reconstructed the group and coefficient
arguments, and the Lean character theorems build without `sorryAx`.

This theorem proves an exact abstract carrier `B` and a minimal faithful
splitting-field representation of `B`.  It does **not** prove:

- a physical `Q8` action on a horn packet;
- physical realization of the off-diagonal `C169^2` states;
- descent of the THM-2886 marked subcopy to the canonical unmarked
  current;
- a lawful simultaneous Hermitian polarization retaining the central
  sign;
- positivity, a row exclusion, or LRC(14).

The cheapest decisive physical test is now precise: construct the
right-action edge `e9=(-9,+9,QB)` on one common marked packet and test
whether the via endpoint differs from the direct endpoint by
`(13,13,-1)`.  This simultaneously tests the second character, ancestry
carry, and quaternionic sign before scalar recombination.
