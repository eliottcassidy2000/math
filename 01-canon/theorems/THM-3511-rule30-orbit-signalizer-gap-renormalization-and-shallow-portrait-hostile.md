---
id: THM-3511
title: "Rule 30 orbit-signalizer gap renormalization and shallow-portrait hostile"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The packed Rule 30 map is a
  three-state binary-tree automorphism.  Its first-return section at each
  innovation cylinder is an exact orbit-signalizer whose square, restricted
  below the first new moved bit, is the next signalizer.  This gives a
  three-letter, one-control-bit finite presentation of every innovation gap.
  It is not a proved finite graph or a bounded-gap theorem.  Exact physical
  signalizers at scales 7 and 17 have the same depth-three portrait but gaps
  8 and 5, so shallow portrait quotients require an overflow/owner sidecar.
  No Rule 30 prize consequence is claimed.
source: root/rule30-normalized-displacement/orbit-signalizer/2026-08-16
audit: >
  PASS (2026-08-16), independent proof, scope, hostile, and replay audit.
  The auditor rederived the right-action section order, orbit-signalizer
  recurrence, marked and all-phase owners, max-plus implication, eventual
  dimension and Ahlfors-regularity consequences, and the bounded-activity
  boundary.  Independent packed controls reproduced the physical shallow
  hostile.  Ordinary, optimized, and stored transcripts are byte-identical.
depends_on:
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
  - THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary
  - THM-3493-rule30-dyadic-wrap-atlas
  - THM-3500-rule30-dyadic-section-cut-defect-and-cross-depth-valuation-carrier
  - THM-3507-rule30-normalized-dyadic-displacement-sibling-trace-and-assouad-spectrum
related:
  - THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum
  - THM-3503-rule30-odometer-ultrametric-regrading-and-orbit-closure-dimensions
external:
  - I. V. Bondarenko, N. V. Bondarenko, S. N. Sidki, and F. R. Zapata, "On the conjugacy problem for finite-state automorphisms of regular rooted trees," Groups, Geometry, and Dynamics 7 (2013), 323--355, arXiv:1011.2227 (CITED for orbit-signalizer terminology and finite order-graph method only)
  - Eric S. Rowland, "Local Nested Structure in Rule 30," Complex Systems 16 (2006), 239--258, DOI 10.25088/ComplexSystems.16.3.239 (CITED for the identical nestedness sequence and values through scale 40)
script: 04-computation/rule30_orbit_signalizer_gap_thm3511.py
output: 05-knowledge/results/rule30_orbit_signalizer_gap_thm3511.out
script_sha256: 2ce110f0b8e9c71c3d298aaf07e8e6c02b70d33e5671bc763f3f3b490caa5445
output_sha256: 8c7599690b4eb013f94c3928c7fd6906979bdd4f5e1a1e866cfbfd2753ae8a51
hash_basis: raw bytes after audit promotion
---

# THM-3511 -- Rule 30 orbit-signalizer gap renormalization and shallow-portrait hostile

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-3500 makes the dyadic cut defect physical, while THM-3507 identifies the
next innovation gap with uniform signed sibling cancellation.  The remaining
question is what object actually renormalizes from one dyadic time scale to
the next.  It is not a fixed list of Hasse moments.  It is a first-return
**section** of the finite Rule 30 tree automorphism.

## 1. Inheritance, conventions, and the live board

Read an infinite binary word from least to most significant bit.  The packed
Rule 30 map and its two noninitial Mealy states have wreath recursions

```text
A=(A,B),             B=(C,B)sigma,             C=(A,B)sigma,       (1)
```

where `sigma` interchanges `0` and `1`.  Products are right actions: `gh`
means apply `g` first and then `h`, and

```text
(gh)|_u = g|_u h|_(u^g).                                      (2)
```

The boundary ray

```text
xi=10^infinity                                                   (3)
```

is the integer seed `1`, and `xi^(A^t)` is the packed state `R_t`.

The written factor lists in THM-3463 use the opposite display convention:
later displayed factors act first, so their section input is strict suffix
activity.  Reversing such a list gives the right-action words used here.  For
example, this theorem writes

```text
B^2|_0=CB,
```

where THM-3463's factor-list convention displays the same automorphism as
`BC`.  This is a notation conversion, not a different section order.

The inheritance pass is:

1. closest proved mechanism: THM-3500's lossless dyadic section cut and
   THM-3507's sibling trace;
2. canonical hostile: `s_7` and `s_17` below have the same shallow portrait
   but different gap costs;
3. corrected near miss: keep the first-return section itself, not a bounded
   portrait or marginal Hasse list; and
4. least-used sidecars: the marked child, the ordered product of its two
   sibling sections, and ordinary subtraction borrow.

The live board is the rooted-tree section, orbit signalizer, first-return
induction, anti-sibling common prefix, ordered owner, borrow, overflow, and
max-plus cycle potential.

## 2. The exact first-return signalizer

Retain THM-3507's notation

```text
q_m=2^m,       v_m=nu_2(R_(q_m)-1),
d_m=v_(m+1)-v_m.                                      (4)
```

Let

```text
u_m=10^(v_m-1),       s_m=(A^(q_m))|_(u_m).            (5)
```

The word `u_m` has length `v_m`.  Since the exact seed period at that width
is `P_(v_m)=q_m`, its `A`-orbit has length `q_m`.  Thus `s_m` is literally an
element of the orbit-signalizer

```text
OS(A)={ (A^|Orb_A(u)|)|_u : u in {0,1}^* }.            (6)
```

The term is imported from Bondarenko--Bondarenko--Sidki--Zapata; every
identity below is derived directly from (1)--(5).

### Theorem 2.1 (gap renormalization)

Every `s_m` is active at its root, `s_0=B`, and

```text
boxed:
d_m=nu_2( [ (0^infinity)^(s_m^2) ]_2 ),               (7)

s_(m+1)=(s_m^2)|_(0^(d_m)).                            (8)
```

Here `[z]_2` is the 2-adic integer whose low-to-high digit ray is `z`.
Equivalently, `s_m^2` fixes the marked zero ray through its first `d_m`
digits and its section after those digits is active.

### Proof

The return `A^(q_m)` fixes `u_m`, but it changes the next seed bit because
`v_m` is the exact first changed coordinate.  Its section `s_m` therefore
interchanges the two children.  At `m=0`, `u_0=1` and (1) gives `s_0=A|_1=B`.

Since `A^(q_m)` fixes `u_m`, the product rule (2) gives

```text
(A^(2q_m))|_(u_m)=s_m^2.                               (9)
```

The seed tail after `u_m` is `0^infinity`.  By definition, `A^(2q_m)` first
changes the seed at coordinate `v_(m+1)=v_m+d_m`.  Hence `s_m^2` first changes
the zero tail at its coordinate `d_m`, proving (7).  The longer fixed prefix
is

```text
u_(m+1)=u_m 0^(d_m).                                   (10)
```

Taking the section of (9) at the appended zeros proves (8).

This is a first-return or Kakutani induction on nested seed cylinders.  It is
also the exact anti-sibling common-prefix picture: `d_m` is the Gromov-product
length by which the two odd sibling directions cancel.  Neither analogy adds
a finiteness claim.

### 2.1 Exact identification with Rowland's nestedness sequence

Rowland's shifted and reflected Rule 86 coordinate defines `lambda_I(t)` as
the initial black-run length in row `t-1`, and puts `a(m)=lambda_I(2^m)`.
Under the coordinate equivalence in THM-3458, this is the initial low-bit run
of ones in `R_(2^m-1)`.  Write `Phi` for the packed update represented by
`A`.  If a packed word `x` has exactly `L` initial one bits, then its low-bit
rule gives

```text
nu_2(Phi(x)-1)=L.                                     (10a)
```

Indeed bit zero remains one, the intervening bits cancel, and bit `L` is the
first disagreement with the seed.  Therefore

```text
boxed: a(m)=lambda_I(2^m)=nu_2(R_(2^m)-1)=v_m.        (10b)
```

Thus Rowland's published values through `m=40` are finite historical data for
the same sequence, not a merely similar statistic.  His paper proves the
nestedness/strict-increase mechanism and reports the continuation through
`v_40=93`; the cross-depth signalizer recurrence (7)--(15) is the additional
object proved here.  The cited continuation is not used as a verifier gate,
and no literature novelty claim is made.

There is nevertheless an exact **CITED-FINITE** corollary.  The published
list has `v_m<2^m` for every `3<=m<=40`.  THM-3493's dyadic wrap atlas,
together with `v_0=1,v_1=3,v_2=4`, therefore gives

```text
W intersect [1,2^41-1]={1,2,3,4},
[5,2^41-1] subset H.                                (10c)
```

The endpoint `2^41-1=2,199,023,255,551` is inherited from Rowland's cited
finite computation; it is not independently replayed by this companion and
has no asymptotic consequence.

### 2.2 Why Rowland's fixed-diagonal doubling does not bound `d_m`

Rowland's period-doubling proposition in his Section 5 uses the opposite,
fixed-diagonal chart.  In centered seed coordinates `A_t(j)`, write

```text
b_k(t)=A_t(t-k),              ell_m(t)=A_t(m-t).
```

Then `ell_m(t)=b_(2t-m)(t)`: a fixed `ell`-column cuts diagonally across the
right-edge sheet.  The two exact recurrences are correspondingly different:

```text
b_k(t+1)=b_k(t) xor (b_(k-1)(t) or b_(k-2)(t)),

ell_m(t+1)=ell_(m-2)(t) xor ell_(m-1)(t) xor ell_m(t)
             xor ell_(m-1)(t)ell_m(t).               (10d)
```

When Rowland's survival gate `ell_(m-1)=0` holds, the second line reduces to
`ell_m(t+1) xor ell_m(t)=ell_(m-2)(t)`.  Its doubling test and THM-3463's
`epsilon` are both odd cycle-holonomy obstructions.  They are not the same current: a
Rowland doubling makes `ell_m` nonwhite and blocks an adjacent doubling,
whereas the physical right-edge lift has `epsilon_3=epsilon_4=1`.  Thus the
shared cohomological template transfers, but the tempting stronger
no-adjacent-gap constraint does not.  Outside the white gate, Rowland's
one-period fiber can synchronize and erase its owner; the packed lift remains
an invertible translation extension.  Finally, the center is the moving
observer `c_t=ell_t(t)`, not any fixed `ell_m` column, so fixed-column eventual
periods do not determine a center value.

## 3. A three-letter, one-control-bit finite presentation

Represent a product of states by a word `w in {A,B,C}^*`, applied from left
to right.  Put

```text
alpha(w)=number of letters B or C mod 2.               (11)
```

Thus `alpha(w)` is its root activity.  For `x in {0,1}`, let `tau_x(w)` be
the section word at the root input `x`.  It is produced by one left-to-right
two-state transducer.  Immediately before reading a letter `g`, with control
bit `x`, output

```text
control x=0:      A -> A,      B -> C,      C -> A,
control x=1:      A -> B,      B -> B,      C -> B,      (12)
```

and then replace

```text
x <- x xor 1_[g in {B,C}].                             (13)
```

The final control is `x xor alpha(w)`.  These are exactly the sections and
root permutations in (1), so no automaton quotient has been taken.

Let `w_m` represent `s_m`.  Since `alpha(w_m)=1`, the marked zero section of
the square is the **ordered** word

```text
h_m=tau_0(w_m) tau_1(w_m).                             (14)
```

Consequently (7)--(8) are the explicit recursion

```text
boxed:
w_0=B,
d_m=min{r>=1: alpha(tau_0^(r-1)(h_m))=1},
w_(m+1)=tau_0^(d_m-1)(h_m).                            (15)
```

In particular `|w_m|=2^m`.  Formula (15) is a finite **presentation** of the
gap path: three output letters plus one control bit.  Its represented group
element may still carry unbounded information.  Calling (15) a finite graph
would erase precisely the question that remains open.

The order in (14) is load-bearing.  It records that the marked child is `0`,
then the active first copy sends it to child `1`.  The other child has section
word `tau_1(w_m)tau_0(w_m)` and cannot be substituted silently.

## 4. Marked signed units and the borrow sidecar

Put

```text
Z(s)=[(0^infinity)^s]_2.                               (16)
```

The seed and its first two dyadic translates share the prefix `u_m`, and
their successive tails are

```text
0^infinity,       (0^infinity)^s_m,
                  (0^infinity)^(s_m^2).               (17)
```

Therefore THM-3507's marked normalized displacements become

```text
boxed:
U_m(0)=Z(s_m),
U_m(q_m)=Z(s_m^2)-Z(s_m),                              (18)

U_m(0)+U_m(q_m)=Z(s_m^2)
                =2^(d_m)Z(s_(m+1)).                   (19)
```

This is the sibling trace as a literal first-return section identity.  The
first unit in (18) needs no borrow because its ordered base tail is zero.  To
decode the second signed difference, put `x_j=bit_j Z(s_m)`,
`y_j=bit_j Z(s_m^2)`, `h_j=x_j xor y_j`, and start `c_0=0`.  THM-3507's
ordinary subtraction transducer remains exact:

```text
u_j=h_j xor c_j,
c_(j+1)=((1-y_j)x_j) or (c_j(1-h_j)).                 (20)
```

The section word preserves the ordered endpoints from which (20) is derived;
an XOR portrait alone does not.

## 5. Full-phase owners and the principal sibling ratio

The marked section `s_m` is enough for the gap and the units in (18), but it
is not the full all-phase carrier.  For a phase `t`, let `u_(m,t)` and
`eta_(m,t)` be the length-`v_m` prefix and remaining tail of `iota(t)`, and
put

```text
s_(m,t)=(A^(q_m))|_(u_(m,t)).                          (21)
```

The three phase tails are then

```text
eta,          eta^s,          eta^(s^2),
               where s=s_(m,t).                       (22)
```

Thus, exactly,

```text
U_m(t)      =[eta^s]_2-[eta]_2,
U_m(t+q_m) =[eta^(s^2)]_2-[eta^s]_2.                  (23)
```

If the principal anti-sibling ratio is normalized as

```text
G_m(t)=-U_m(t+q_m)/U_m(t),                            (24)
```

then the owner pair `(s_(m,t),eta_(m,t))` recovers it literally:

```text
boxed:
G_m(t)=-([eta^(s^2)]_2-[eta^s]_2)
          /([eta^s]_2-[eta]_2),                       (25)

nu_2(G_m(t)-1)=d_m.                                   (26)
```

Both denominator and numerator in (25) are odd units.  The sign in (24) is
what makes `G_m(t)=1 mod 2^(d_m)`.

For integer `0<=t<q_m`, let

```text
r_(m,t)=(A^t)|_(u_m).                                 (27)
```

Commutativity of `A^t` and `A^(q_m)`, followed by (2), gives

```text
u_(m,t)=u_m^(A^t),
eta_(m,t)=(0^infinity)^(r_(m,t)),
s_(m,t)=r_(m,t)^(-1) s_m r_(m,t).                     (28)
```

An arbitrary rooted-tree conjugacy preserves common-prefix lengths, so (28)
transports the valuation in (26).  It is not additive in the 2-adic integer
coordinate and need not preserve the signed ratio itself.  Therefore feeding
a shifted phase tail into the marked `s_m` alone recovers neither (23) nor
(25); one must retain the phase owner `s_(m,t)`, equivalently `r_(m,t)`.

More intrinsically, if `a` is the first `d_m` digits of `eta`, the complete
phase renormalization is

```text
(s,eta) -> ((s^2)|_a, shift^(d_m)(eta)).               (29)
```

THM-3507 proves that the cost in (29) is phase-independent.  The states and
borrow paths need not be.

## 6. Exact shallow hostile and finite lower bound

Recursion (15) gives the physical signalizer words.  On the depth-three
binary tree, both `s_7` and `s_17` induce

```text
(0,1,2,3,4,5,6,7) -> (7,6,5,4,3,2,1,0),             (30)
```

but

```text
d_7=8,          d_17=5.                               (31)
```

Their depth-four portraits, listed by the images of `0,...,15` in
hexadecimal, are respectively

```text
s_7:  7e5cb290f6d43a18,
s_17: f6543a987edcb210.                               (32)
```

This is an actual Rule 30 hostile, not an ambient automorphism pair.  A
depth-three state sends both transitions to overflow and cannot determine
the edge cost.  The first `23` physical signalizers

```text
s_0,s_1,...,s_22                                      (33)
```

are pairwise distinct already on the depth-four tree.  Hence any exact finite
orbit-signalizer graph using the actual automorphisms as vertices has at
least `23` vertices.  Syntactic word length alone is not used for this claim.

The companion freezes the universe `0<=m<=22`.  It checks all `23`
renormalization costs, direct sections of `B^(2^m)` through `m=16`, marked
units and independent packed states, every phase owner through `m=7`, the
borrow decoder, (30)--(33), and ordinary/optimized determinism.  Its canonical
state transcript digest is

```text
d06c6d7291bf043d028760f0dc1af50a225aa103a50a14084715a0ba9672e5ca. (34)
```

This is finite exact evidence for the universal identities proved above and
a finite exact hostile.  It is not evidence that only finitely many
signalizers occur.

## 7. The exact finite-graph/max-plus target

Let `R(s)` and `d(s)` denote the operations in (7)--(8) on an active marked
signalizer.  A successful certificate consists of a finite set `S` of
finite-state automorphisms such that

```text
B in S,       R(S) subset S,                          (35)
```

together with exact Mealy-equivalence checks for every declared transition.
If a potential satisfies

```text
d(s)+psi(R(s))-psi(s)<=C                 for s in S,  (36)
```

then telescoping gives

```text
v_m<=Cm+O(1),       dim_H X>=1/C.                     (37)
```

For a finite directed graph, (36) is feasible exactly when every directed
cycle has mean cost at most `C`; this is the max-plus cycle-mean dual of
THM-3502's Perron subeigenvector certificate.  In the deterministic closed
case the actual gap word is eventually periodic.  If its eventual cycle has
`L` edges and total cost `D`, then the innovation word is eventually periodic
in state width with density `L/D`.  The Hausdorff, packing, lower/upper box,
lower-Assouad, and Assouad dimensions in THM-3503 and THM-3507 all equal
`L/D`; the orbit closure is also Ahlfors `L/D`-regular under the pushed-forward
Haar measure.

No such finite closed `S` is proved here.  The Rule 30 automorphism is not in
the bounded-activity class for which the cited orbit-signalizer paper proves
automatic finiteness: all `2^n` level-`n` sections of `A` lie in the
nontrivial set `{A,B,C}`, so its activity is `2^n`.  Equations
(30)--(32) are the mandatory overflow test for any proposed quotient.

## 8. Preservation, loss, and scope

| map | preserved | destroyed | required sidecar / cheapest test |
|---|---|---|---|
| `A^(q_m)` at the seed cylinder -> `s_m` | exact first-return action below the fixed prefix | absolute prefix `u_m` | retain `v_m`; check (5) directly |
| `s_m -> (d_m,s_(m+1))` | exact gap and next marked return section | no information if the full section is retained | ordered marked child `0`; (7)--(8) |
| `w_m -> depth-K portrait` | action through `K` tail bits | every deeper section and overflow length | overflow state; `s_7/s_17` at `K=3` |
| `(s_(m,t),eta) -> G_m(t)` | complete signed sibling ratio | none | ordered tail plus borrow |
| `s_m -> shifted phase` without `r_(m,t)` | marked valuation analogy only | phase owner and signed unit | transport section in (28) |
| finite closed signalizer graph -> max-plus bound | all accumulated gap costs | internal word presentations after equality audit | Mealy equivalence and every cycle mean |

There is no intrinsic tournament: the two children are ordered by the marked
orbit and squaring.  This theorem proves no finite signalizer graph, bounded
gap, positive lower or Hausdorff dimension, center nonperiodicity, balance,
prediction lower bound, or Rule 30 prize.  No literature novelty or priority
claim is made.

Run

```bash
python3 04-computation/rule30_orbit_signalizer_gap_thm3511.py
python3 -O 04-computation/rule30_orbit_signalizer_gap_thm3511.py
```

and compare both byte-for-byte with the stored output.  The companion uses no
`assert` gates.
