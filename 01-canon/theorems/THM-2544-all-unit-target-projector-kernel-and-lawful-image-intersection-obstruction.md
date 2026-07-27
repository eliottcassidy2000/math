---
id: THM-2544
title: "All-unit target projector kernel and lawful-image intersection obstruction"
status: >
  PROVED + VERIFIED-EXACT.  For every primitive nine-coordinate scalar
  word of THM-2309 type, the unrestricted mod-thirteen target pushforward U
  on the full mod-91 relation-residue space has rank 169 and kernel dimension
  91^8-169.  The all-coordinate-91-unit pushforward J has rank zero at the
  sharp septimal-support-one boundary and rank 169 whenever at least two
  scalar coordinates are septimal units.  In the latter case the joint map
  (U,J) has rank 338, J restricted to ker(U) is onto all 169 target
  coordinates, and ker(U) contains 169 independent target-coordinate
  switches, including the designated 168 nonzero-target switches.  A
  169-term current supported on the zero-septimal CRT section has full
  unrestricted target support and identically zero all-unit target
  aggregate.  Consequently even the exact value of A does not determine B
  on the arbitrary current space.  For actual THM-2334 currents, Uc=A and
  Jc=B exactly, so all remaining force must come from the constrained lawful
  Abel-current image: the strong uniform condition is avoidance of the
  nonzero-target kernel, while the rowwise existential condition has a
  different quantifier.  No hostile vector is proved lawful, no covering row
  is exhibited or excluded, no same-ancestry arrival field is constructed,
  and LRC(14) remains open.
source: codex-2026-07-27-all-unit-target-kernel
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2325-prescribed-target-gain-full-lattice-91-unit-needle-bank
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2349-first-depth-one-delayed-shallow-restart
related:
  - THM-2536-deep-comb-selector-flow-target-drift-and-centered-toothpick-boundary
  - THM-2538-anchored-transverse-gain-and-common-ancestry-arrival-boundary
  - THM-2540-weighted-live-event-kakeya-flux-and-transverse-gain-boundary-refinement
  - THM-2541-canonical-typed-row-full-target-plane-support
  - THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile
script: 04-computation/lrc14_all_unit_target_projector_kernel_thm2544.py
output: 05-knowledge/results/lrc14_all_unit_target_projector_kernel_thm2544.out
script_sha256: 36bb8ff1a2901add0cea3e2e46d63c66615eed223869704fad38688f563012fb
output_sha256: 9b0254e3ac8e9c7303445f6d3801d9f515ed5536cd8122b46a0daeaf8021724b
hash_basis: working-tree bytes (LF)
---

# THM-2544 -- the all-unit projector is transverse to the unrestricted one

**PROVED + VERIFIED-EXACT.**

THM-2541 proves that all `169` unrestricted target aggregates `A(q)` are
nonzero on one canonical typed non-cover row.  THM-2334 carefully warns that
the needed object is instead the all-`91`-unit aggregate `B(q)`.  This theorem
computes the exact information loss between those two objects before any
analytic structure is used.

The answer is maximally hostile:

```text
arbitrary residue currents:
  (A,B) are independent coordinates;

lawful Abel-boundary currents:
  a highly constrained nonlinear image inside that ambient space;

remaining problem:
  determine how that lawful image meets the B-nonzero-target kernel.       (1)
```

Thus the gap is not a missing finite-field support count.  It is a theorem
about the analytic image which has not yet been proved.

## 1. The two exact projectors

Use the scalar setup of THM-2309
[`THM-2309-owner-aligned-pivot-packets-and-visible-height-separation.md`]
with `d=9`.  Thus `w in Z^9` is primitive, its six guard/unit entries are
nonzero modulo `13`, its owner and two target blockers vanish modulo `13`,
and

```text
K_N={y in (Z/NZ)^9:y.w=0 mod N}.                         (2)
```

Fix an owner-aligned packet `L_13 subset K_13`.  THM-2309 gives

```text
dim K_13=8,                 dim L_13=6,

K_13=L_13 direct_sum span(e_a,e_b),

G:=K_13/L_13 isomorphic to F_13^2.                       (3)
```

Let `pi:K_13->G` be the quotient map and put

```text
Omega=K_91,
Q:Omega->G,                  Q(y)=pi(y mod 13),

S={y in Omega: every coordinate y_i is a unit mod 91}.  (4)
```

For any field `F`, let `V=F^Omega`, with coordinate vectors `delta_y`.
Define two linear maps

```text
U:V->F^G,        (Uc)(q)=sum_(Q(y)=q)c(y),

J:V->F^G,        (Jc)(q)=sum_(Q(y)=q, y in S)c(y).       (5)
```

`U` forgets the mod-seven residue and the mod-thirteen packet coordinate.
`J` first imposes the all-coordinate-unit mask.  When `c` is the exact
THM-2334 residue current, these are precisely `A` and `B`; Section 5 proves
that identification without an informal change of modulus.

## 2. CRT fibres and the nonunit section

Primitivity makes the scalar functional onto modulo both primes.  CRT and
THM-2334 (25) give

```text
K_91 isomorphic to K_7 x K_13,

|Omega|=91^8,

|Q^(-1)(q)|=|K_7| |L_13|=7^8 13^6                 (6)
```

for every `q in G`.  In particular `Q` is onto.

There is an especially useful section which is completely disjoint from
`S`.  Write `q=(q_a,q_b)` under (3), put

```text
x_q=q_a e_a+q_b e_b in K_13,

v_q=CRT(0,x_q) in K_7 x K_13=K_91.                 (7)
```

Every coordinate of `v_q` is divisible by `7`, so `v_q notin S`, while
`Q(v_q)=q`.  Therefore the `169` disjoint rows of `U` have independent
coordinate sections:

```text
rank U=169,

dim ker U=91^8-169.                                (8)
```

This part uses no septimal-support hypothesis.

## 3. The sharp all-unit rank dichotomy

Let

```text
s=|{i:7 does not divide w_i}|.                      (9)
```

Primitivity gives `s>=1`.  At the sharp boundary `s=1`, the unique
septimally supported coordinate of every vector in `K_7` must vanish.
Hence `S` is empty and

```text
s=1  ->  J=0.                                      (10)
```

Now assume `s>=2`.  The number of all-coordinate-nonzero points of `K_7`
is the THM-2325 character count

```text
N_7(s)
 =6^(9-s)[6^s+6(-1)^s]/7
 >=5*6^7
 =1,399,680.                                       (11)
```

At `13`, every affine target fibre is

```text
X_q=q+L_13.                                        (12)
```

The nine forbidden equations `x_i=0` are proper affine hyperplanes in the
six packet coordinates because THM-2309's packet is bright in all nine
columns.  Their normals span the six-dimensional dual because the packet
has row rank six.  The essential-arrangement estimate used in THM-2325
therefore gives

```text
|X_q intersection (F_13^*)^9|
 >=(13-9+6-1)12^5
 =9*12^5
 =2,239,488                                        (13)
```

for **every** `q in G`.  THM-2325 stated (13) for nonzero target vectors
because it was counting projective needles; its proof does not use
`q!=0`, and (13) supplies the previously unneeded zero-fibre boundary too.

Pair (11) and (13) by CRT.  Every target fibre contains at least

```text
N_7(s) 9*12^5
 >=3,134,566,563,840                               (14)
```

all-`91`-unit points.  Hence every row of `J` is nonempty and the rows have
disjoint supports.  It follows that

```text
s>=2  ->  rank J=169,
           dim ker J=91^8-169.                     (15)
```

Equations (10) and (15) are the exact all-unit rank dichotomy for this
scalar type.

## 4. Split surjectivity and the 169 kernel switches

Continue under `s>=2`.  Choose one

```text
u_q in S intersection Q^(-1)(q)                   (16)
```

for each target, and retain the nonunit point `v_q` from (7).  Define

```text
h_q=delta_(u_q)-delta_(v_q).                        (17)
```

Then

```text
Uh_q=0,                    Jh_q=delta_q.            (18)
```

Thus the `h_q` are linearly independent, `J|ker U` is onto, and

```text
rank(U,J)=338,

rank(J|ker U)=169,

dim(ker U intersection ker J)=91^8-338.            (19)
```

The subbank

```text
{h_q:q!=0}                                         (20)
```

has exactly `168` independent members and changes the `168` designated
nonzero-target `B` coordinates one at a time without changing `A` at all.

More strongly, for arbitrary vectors `a,b in F^G`, the current

```text
c_(a,b)
 =sum_q a(q)delta_(v_q)+sum_q b(q)h_q              (21)
```

satisfies

```text
Uc_(a,b)=a,                    Jc_(a,b)=b.           (22)
```

So `(U,J):V->F^G x F^G` is a split surjection.  Taking `a(q)=1` and
`b(q)=0` gives the requested hostile

```text
c_host=sum_q delta_(v_q),

support(Uc_host)=G,             Jc_host=0.          (23)
```

This is stronger than failure of a support implication.  Even after the
**exact value** of `A=Uc` is fixed, `B=Jc` ranges over all of `F^G` on the
ambient affine fibre.  In particular no function of `A` alone, linear or
nonlinear, determines `B` on arbitrary residue currents.

## 5. Identification with THM-2334's analytic current

For `rho<1`, let

```text
c_rho(y)=C_(rho,91)(y),             y in K_91,       (24)
```

be THM-2334's absolutely convergent residue pushforward
[`THM-2334-relation-residue-current-and-character-twist-pushforward.md`].
Then finite regrouping gives

```text
(Uc_rho)(q)
 =sum_(pi(y mod 13)=q) sum_(r in Lambda; r=y mod 91) C_rho(r)

 =sum_(r in Lambda; pi(r mod 13)=q) C_rho(r)

 =A_rho(q).                                         (25)
```

The mask defining `J` is exactly THM-2334 (46)--(47), so

```text
(Jc_rho)(q)=B_rho(q).                               (26)
```

There are finitely many residue coordinates.  Passing to the Abel boundary
in (25)--(26) therefore yields, with `c(y)=C_91(y;1-)`,

```text
Uc=A,                         Jc=B.                 (27)
```

The hostile (23) and switches (17) are elements of the same finite vector
space as (27), but they have not been shown to satisfy the analytic
factorization, word positivity, phase covariance, or incidence constraints
which define a THM-2334 current.  They are ambient obstructions, not
counterexamples to a lawful-current theorem.

## 6. The lawful-image intersection is the live object

For a hypothetical covering row `w` among THM-2349's `165` rows, let

```text
L_law(w) subset C^(K_91)                            (28)
```

be the set of all boundary vectors (24) obtained from THM-2349-licensed
choices of delayed owner, literal word, marked unit triangle, and endpoint
data.  Let

```text
L_cov=union_(hypothetical covering rows w)L_law(w),

J_* = pr_(G minus {0}) o J.                         (29)
```

For one licensed datum, THM-2334 (49) is exactly

```text
c notin ker J_*.                                    (30)
```

Consequently the strong statement that **every** licensed covering-row
current survives is exactly the image-avoidance condition

```text
L_cov intersection ker J_* = empty.                (31)
```

If one is free to choose favourable delayed data separately for each row,
the correct weaker quantifier is instead

```text
for every hypothetical covering row w,
    L_law(w) is not a subset of ker J_*.             (32)
```

For a proposed canonical selection `s(w) in L_law(w)`, uniform survival is
the corresponding intersection statement

```text
s(Cov) intersection ker J_* = empty.                (33)
```

These distinctions matter.  One lawful kernel hit would refute (31), but
would not refute the existential route (32).  Conversely, exhibiting the
ambient hostile (23) says nothing about any of (31)--(33), because membership
in `L_cov` is the unproved step.  If no covering row exists, (31) is vacuous;
it cannot by itself serve as an LRC proof.

This reframes the current frontier as a restriction problem:

```text
describe equations, inequalities, recurrences, or equivariant sidecars
cutting out L_law inside C^(K_91), then test their transversality to ker J_*.
                                                                  (34)
```

Another unrestricted target-support computation cannot settle (34).

## 7. What the two positive controls do not transfer

### THM-2541

THM-2541
[`THM-2541-canonical-typed-row-full-target-plane-support.md`] proves
`A(q)!=0` for all `169` targets on its canonical typed non-cover datum.
Equation (21) shows that even those exact `169` complex values are compatible,
in the ambient current space, with `B=0` and with every other prescribed
`B`.  THM-2541 can constrain `B` only through additional equations satisfied
by its particular analytic current, not through the value or support of `A`.

### The 42-cut ancestry pairing

THM-2545
[`THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile.md`]
has now given the complementary event-level diagnosis: the same cut data can
be tensored with aligned and swapped two-root couplings while preserving every
recorded marginal and changing the physical diagonal hit from full mass to
zero.  The projector calculation here is the residue-current analogue of that
transportation hostile.

The exact artifact

```text
04-computation/lrc14_cut_bundle_ancestry_pairing_opus_20260727.py
```

finds `432/432` nonzero equivariant cut/stalk pairings and retains target
charge and owner type.  Its own scope also proves that both the kappa-torsor
and cut-torsor invariant averages vanish; the surviving object retains those
labels as sidecars.  It is not the all-unit address projector (5), and no
typed map from its equivariant coefficient space to `V=C^(K_91)` has been
constructed.

Therefore it presently imposes no restriction on the exact ambient fact

```text
J(ker U)=C^G.                                       (35)
```

To affect (31)--(34), it needs a branch-transplant/intertwining map `T` whose
image contains the relevant lawful currents and for which the cut/stalk
nonvanishing controls `J_*T`.  Averaging away its two torsor sidecars cannot
be that map, because those invariant readouts are exactly zero.  THM-2545's
exact tensor hostile confirms that retaining the current recorded data does
not repair this loss.  The 42-cut result is promising input to a lawful-image
description, but it does not repair `A -> B` by itself.

The same semantic guard applies to the later-arrival program: a nonzero
all-unit relation-address aggregate is still not a genuine later Boolean
arrival field on one common ancestry base.  That extra landing map remains
separate.

## 8. Exact referee

The companion uses the canonical typed word

```text
w=(1,14,27,40,53,66,13,2197,742586).                (36)
```

It reconstructs THM-2309's exact owner packet, checks its rank-six pivot and
nine bright columns, builds a common all-unit `K_7` witness, and constructs
both sections `u_q,v_q` for all `169` targets.  It then verifies the hostile
(23), all `169` switches (17), and the `168` nonzero-target subbank without
enumerating `K_91`.  For this word `s=8`, and the exact lower bound furnished
by (11)--(13) is

```text
3,224,137,125,888 all-unit points per target fibre.  (37)
```

The script performs `2,756` optimization-safe exact checks.  Reproduce with

```bash
python3 04-computation/lrc14_all_unit_target_projector_kernel_thm2544.py
python3 -O 04-computation/lrc14_all_unit_target_projector_kernel_thm2544.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_all_unit_target_projector_kernel_thm2544.out
```

byte-for-byte after LF normalization.

## 9. Scope

This theorem proves an exact finite-dimensional obstruction and identifies
the lawful-image intersection which a future analytic theorem must control.
It does **not** prove that the hostile or any switch is a lawful Abel current,
does not transfer the 42-cut sidecars into the all-unit projector, does not
construct a same-ancestry later arrival, does not remove a scalar row, and
does not prove LRC(14).
