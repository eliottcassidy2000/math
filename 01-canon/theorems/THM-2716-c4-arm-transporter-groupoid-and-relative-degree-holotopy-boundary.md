---
id: THM-2716
title: "C4 arm-transporter groupoid and relative-degree holotopy boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Relative to the
  fixed C4 extension and named endpoint/ghost objects, the two odd arms of
  the THM-2706 macro cycle generate the minimal two-object/eight-arrow
  transporter groupoid.  Each cross-Hom is a principal C2 bitorsor; its two
  sections and the two macro/debt alignments form unbased C2 torsors.  The
  integral arm norm acts by 2 on the invariant line and by 0 on the sign
  line; THM-2710's independent half sheet changes this to character
  boundaries +2 and -2, not cancellation.  A pointed-cone addendum proves
  that quarter-turn cancellation is intrinsically signed.  The conditional
  Phi7 middle-factor cancellation kills a nonzero global debt through a
  four-dimensional kernel and is not a physical projector or edge gain.
  This is a relative Z2-graded support/incidence groupoid, not a physical
  one-step action, semantic word, endpoint current, row exclusion, or LRC(14)
  conclusion.
source: root/c4-arm-transporter-holotopy-2026-07-28
depends_on:
  - THM-2041-frobenius-stability-of-exact-period-projectors
  - THM-2701-literal-singleton-word-one-step-dilation-nilpotence
  - THM-2706-relative-segal-macro-cycle-and-minimal-ghost-midpoint-completion
  - THM-2707-full-physical-lift-fibre-common-simplex-and-packet-scc
  - THM-2710-central-half-phase-literal-word-nilpotence-and-prescribed-clock-invisibility
related:
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2697-filtered-affine-handoff-germ-category-and-base-signature-holotopy-boundary
script: 04-computation/lrc14_c4_arm_transporter_thm2716.py
output: 05-knowledge/results/lrc14_c4_arm_transporter_thm2716.out
script_sha256: 9091b4f9ac921d3dcfb3602858f74fb62ff427b023c75cc500af09e6bdd3b4e2
output_sha256: 71100d3716ce86879cceb5bccbedcff9dfcf45f86e2bc7feebc255f10550af90
secondary_scripts:
  - 04-computation/lrc14_paired_quarter_turn_debt_local_system_20260728.py
secondary_outputs:
  - 05-knowledge/results/lrc14_paired_quarter_turn_debt_local_system_20260728.out
secondary_script_sha256s:
  04-computation/lrc14_paired_quarter_turn_debt_local_system_20260728.py: abf625bd6ba7c874e4f8b0fcddc70c8ee269277116fcf5c1365aa0a4b3b21fb6
secondary_output_sha256s:
  05-knowledge/results/lrc14_paired_quarter_turn_debt_local_system_20260728.out: c8a8245272bc943dce12dcbb6483432d46d538e695f3f7e5e286a03171696bd1
hash_basis: LF-normalized bytes
---

# THM-2716 -- C4 arm-transporter groupoid and relative-degree holotopy boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2706 leaves a literal physical four-cycle after adjoining its forced ghost
midpoints, while THM-2701/2710 show that the old semantic word does not travel
around that cycle.  The right categorical object is neither the parity
quotient alone nor a claimed physical `C2` action.  It is the transporter of
the fixed `C4` extension.  This theorem computes that transporter, its
coefficient-two integral boundary, and the first exact obstruction to turning
it into a semantic/current bridge.

## 1. The relative transporter

Let

```text
C4=<B | B^4=1>,             J=B^2,
H=<J>={1,J},                P=C4/H={E,G}.                (1)
```

Here `E` and `G` are the named endpoint and ghost types.  Define

```text
G = C4 action P,
[p;r]: p -> p+r mod 2,       p in Z/2, r in Z/4.         (2)
```

Composition adds the `r` labels modulo four.  Forgetting `r` gives

```text
pi:G -> Pair(P),                                             (3)
```

whose fibres are

```text
Hom(E,E)=Hom(G,G)={0,2},
Hom(E,G)=Hom(G,E)={1,3}.                                (4)
```

Precomposition by `H_E` and postcomposition by `H_G` make `Hom(E,G)` a
principal `H`-`H` bitorsor; the analogous statement holds for every fibre of
`pi`.  On the incoming fibre, postcomposition with the nontrivial target deck
arrow is

```text
J_G o 1 = 3,                  J_G o 3 = 1.               (5)
```

This is the deck action over the `H` quotient.  It is not categorical
inversion, which reverses the source and target, and it is not the physical
reflection functor, which moves both endpoints.

All assertions follow directly from

```text
[p;r]^(-1)=[p+r;-r].                                      (6)
```

In particular, the smallest subgroupoid of `G` containing the two identities
and the exact labelled arms `u_1,u_3:E->G` contains

```text
1_E,J_E, 1_G,J_G, u_1,u_3,u_1^(-1),u_3^(-1),             (7)
```

and hence all eight arrows.  Identifying the arms kills `J`; retaining the
fixed `C4` incidence forces `(7)`.  This minimality is relative to the fixed
extension and the named objects, not a canonical physical realization.

## 2. The unbased section and alignment torsors

A functorial section of `(3)` chooses `u_1` or `u_3`; the opposite cross arrow
is its inverse.  There are exactly two sections, target-deck translation swaps
them, and no section is `H`-invariant.  Thus the arm fibre is canonical as an
unbased torsor, but it has no canonical `+1` arm.

The same torsor occurs between the two proved four-cycles.  The THM-2706
macro phases and the THM-2701 debt phases are

```text
macro:  4/17 -> 1/17 -> 13/17 -> 16/17 -> 4/17,
types:     E       G        E         G       E;

debt:   7/170 -> 91/170 -> 163/170 -> 79/170 -> 7/170,
types:      B         A          B         A        B.    (8)
```

A `B`-equivariant type map taking `E` to `A` and `G` to `B` has offset `1`
or `3`, and these two alignments are exchanged by `H`.  Consequently the
abstract incidence groupoids are isomorphic, but canon supplies neither an
`H`-invariant alignment nor a physical chain map between their predicates.

## 3. Integral linearization and the half sheet

The full lift of the unique cross arrow of `Pair(P)` is the norm

```text
u=B+B^3=B(1+J).                                           (9)
```

Its two-by-two incidence matrix is all ones.  It acts by multiplication by
`2` on the `J`-invariant line and by zero on the sign line.  This is exactly
the coefficient-two invoice in THM-2701/2706.  Reduction modulo two erases
the positive multiplicity; it does not solve the integral problem.

For THM-2710, adjoin the independent half-sheet coordinate `epsilon in C2`
and let

```text
r:(p,epsilon) -> (p+r,epsilon+r) mod 2.                  (10)
```

The difference `p-epsilon` is invariant, so this action groupoid is the
disjoint union of two copies of `G`.  Both odd arms toggle `epsilon`.
Fourier transformation in that sheet therefore gives the exact character
boundaries

```text
+2, -2,                                                       (11)
```

not cancellation.  THM-2710 does not supply an enriched `H`-stable semantic
fibre or current intertwiner; it closes the unchanged literal endpoints and
the displayed half/foreign-base cycles.

### Pointed-cone guardrail

The signedness is forced, not cosmetic.  Let `C` be a pointed cone in a real
vector space, let `T` be invertible, and let `v!=0`.  Then either of the
following hypotheses forbids symmetric cancellation:

```text
Tv,T^(-1)v in C;                                  or
v in C and T(C) subset C.                                (11a)
```

In the first case, `(T+T^(-1))v=0` would put `Tv=-T^(-1)v` in
`C intersection (-C)`, forcing `Tv=0`.  In the second, multiplying the same
equality by `T` gives `T^2v=-v`; forward invariance puts `T^2v` in `C`, so
pointedness again forces `v=0`.  In particular `T^2=-I` preserves no nonzero
pointed cone, even under the one-sided condition `T(C) subset C`; when
`T(C)=C`, both formulations apply.  A quarter-turn cancellation can therefore
be a signed local coefficient differential, but not by itself THM-2644's one
nonnegative physical transition.

## 4. Relative degree is not physical chronology

The parity functor grades both cross arms by relative degree one in `Z/2`.
Their chronological colours remain `1` and `3`.  The second arm is a
three-step forward arrow.  It can be read as an inverse lift only after
restricting to the exact `B^4=1` orbit and localizing modulo four.  It is not a
global one-step physical map on a neighbourhood.  Thus

```text
principal arm torsor
  != physical principal C2 bibundle,
  != conversion of +3 into a global +1 step,
  != semantic word or endpoint current.                 (12)
```

The displayed THM-2706 physical banks give an immediate hostile:

```text
raw forward/reverse cells:   5850, 4958,
safe forward/reverse cells:  5850, 4386.                 (13)
```

Unequal cardinalities exclude an involution exchanging either pair of full
selected banks.  This does not exclude an enlarged, newly normalized transit
grammar.

## 5. The THM-2707 SCC does not furnish the deck action

THM-2707 gives `3346` fixed-packet addresses with one common open interval and
a complete directed eleven-partite support SCC.  Its physical lift
translations lie in

```text
Z/R,                         R=13^6=4826809.              (14)
```

This cyclic group has odd order.  Hence `2k=0` or `4k=0` implies `k=0`; no
nontrivial homomorphic lift translation realizes `J` or `C4`.  Its closed
support words use tailored increments whose sum is zero, not powers of one
order-four generator.  A chosen four-cycle in the SCC can host the objects of
`(7)` only after an external noncanonical selector.  Common support is not a
deck involution, semantic following atom, or current covariance.

## 6. Frobenius-stable labels and the exact remaining bridge

There is one rigorous cross-frontier survivor.  In

```text
A=F_13[C7]=F_13[u]/(u^7-1),
```

absolute Frobenius sends `u^j` to `u^(13j)=u^(-j)`.  It is the inversion
automorphism, so every already nonzero or unit aggregate label remains
nonzero or a unit under every `13^r` iteration.  Over `F_169`, coefficient
conjugation and exponent inversion together fix each primitive character
idempotent.  Over `F_13`, the six nontrivial characters descend only as the
three inverse pairs `(1,6),(2,5),(3,4)`.

There is a precise quotient-loss warning for the conditional four-row `BABA`
comparison.  If those rows are first identified in one coefficient algebra,
their displayed product `U` lives in

```text
A_prim=F_13[z]/(Phi_7)
  ~= K_3 x K_5 x K_6,
K_a=F_13[z]/(z^2+a z+1),                 dim_F13 K_a=2.  (14a)
```

The symmetric order-four holonomy has exact coordinates

```text
S=U^42+U^(-42)
 =(0,0,8,9,9,8) in the polynomial basis
 =(-2,0,-2) in K_3 x K_5 x K_6.                         (14b)
```

Thus projection to `K_5` has a four-dimensional kernel `K_3 x K_6` and kills
this nonzero global debt.  It does not prove cancellation in `A_prim`.
Moreover
`U` is conditional on a still-unproved common transport of the four rows, and
`U^42` is their cycle holonomy, not an edgewise coefficient gain.  No
target-preserving physical projector or edge transport is supplied here.

This is preservation, not creation: it does not prevent the original
summands from cancelling before an aggregate label is formed.  Therefore the
next bridge is precisely a functor or physical cospan from a semantic-current
groupoid into this unit-labelled transporter.  It must provide a common
selector, phase/amplitude transport, and a lawful exit; more row arithmetic or
another unlabelled support cycle cannot supply those data.

## 7. Exact verification and audit

Run

```text
python 04-computation/lrc14_c4_arm_transporter_thm2716.py
python -O 04-computation/lrc14_c4_arm_transporter_thm2716.py
python 04-computation/lrc14_paired_quarter_turn_debt_local_system_20260728.py
python -O 04-computation/lrc14_paired_quarter_turn_debt_local_system_20260728.py
```

The dependency-free companion enumerates the eight-arrow closure, identities,
inverses and associativity; the hom fibres, deck swap and two sections; the
actual macro/debt phase and type arrays; the norm eigenlines; both half-sheet
components and character boundaries; the odd-deck obstruction; the
`F_13[C7]` Frobenius permutation; and the supplied canonical bank counts in
`(13)`.  Normal and optimized runs equal the stored output exactly.

The secondary companion verifies the three irreducible quadratic factors,
the `2+2+2` CRT dimensions, the conditional product and its three component
orders, both global powers in `(14b)`, the nonzero polynomial symmetric debt,
its `(-2,0,-2)` CRT image, and the middle projection's rank-two/four-kernel
ledger.  Its normal and optimized runs also byte-match the stored output.

An independent hostile audit rederived the groupoid, minimality, bitorsor,
section and alignment torsors, norm-two/half-sheet boundary, relative-degree
typing, THM-2707 odd-deck obstruction, and physical-scope no-go.  It required
the source/target pre/postcomposition convention, the qualified THM-2710
wording, the pointed-cone hypotheses, the conditional/holonomy distinction,
the full CRT loss ledger, and executable norm/type checks now present above.
No theorem defect remained after those repairs.

`QED.`
