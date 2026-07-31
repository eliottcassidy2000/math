---
id: THM-2851
title: "Q3/Q11/Q7 Bockstein holonomy and realization sidecar"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The cyclic
  carry-curvature lemma, q3/q11/q7 unit-holonomy triangle, flat
  affine-address comparison, sharp 13-state sidecar invoice, and oriented
  carry derivative are exact.  The 449-sheet and E3-horn inputs are pinned
  PROVED dependencies.  No physical E3 contraction, row exclusion, or
  LRC(14) conclusion is asserted.
source: root/lrc-q3-q11-q7-bockstein-holonomy-2026-07-28
depends_on:
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
  - THM-2829-q11-semantic-reselection-and-fine-ancestry-phase-obstruction
  - THM-2835-q11-semantic-word-horn-and-bockstein-blind-support-no-go
  - THM-2839-prime-power-unit-mass-full-spectrum-and-q11-response-provenance
  - THM-2847-q3-q11-transverse-endpoint-edge-and-e3-realization-horn
related:
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
  - THM-2460-idempotent-semantic-word-copy-and-word-block-cosupport-boundary
  - THM-2744-relative-present-unit-repair-and-root-zero-overlap-clutch
  - THM-2820-boolean-idempotent-rigidity-and-norm-top-cotangent-jet-no-go
  - THM-2842-ordered-positive-cone-vandermonde-multiplier-observability
script: 04-computation/lrc14_q3_q11_q7_bockstein_holonomy_thm2851.py
output: 05-knowledge/results/lrc14_q3_q11_q7_bockstein_holonomy_thm2851.out
script_sha256: 2227f59c717095da0f2042096ada145de4e3661530c9aa2cc9020f42c8237a8b
output_sha256: 424fd2e83a618f862a5ee1b5f073a93fe236d10cdc5412eab1b54dec5e537eac
hash_basis: LF-normalized bytes
---

# THM-2851 -- q3/q11/q7 Bockstein holonomy and realization sidecar

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The theorem identifies the first genuinely higher attachment behind the
THM-2847 rank-one horn.  The coefficient sequence splits over its cyclotomic
field, but the natural residue lifts on its comparison nerve and the affine
address transport do not have the same curvature:

```text
comparison-nerve lifts:  L_9 L_8 = T L_4,
affine address lifts:    A_5 A_3 =   A_8.                  (1)
```

Here `T` is one positive ancestry carry.  A flat equivariant contraction
therefore kills `T`.  Retaining even its first residue requires an additional
13-state fibre; retaining target residue and carry together requires at least
169 states.  This lower bound is sharp abstractly, but no such physical
sidecar is constructed here.  The exact signed carry derivative does exist:
it has rank `12` on the first carry quotient and rank `13^5-1` on the full
ancestry orbit.

## 1. Cyclic carry curvature

The following lemma is independent of LRC.

Let `p` be prime and `d>=1`.  Write a point of `C_(p^(d+1))` uniquely as

```text
n=p a+q,       a in C_(p^d),       q in {0,...,p-1}.       (2)
```

For `h in {0,...,p-1}`, define the natural lifted residue translation

```text
L_h(a,q)
 =(a+kappa(q,h), q+h mod p),

kappa(q,h)=floor((q+h)/p),                                (3)
```

and let

```text
T(a,q)=(a+1,q).                                            (4)
```

For all `h,k`,

```text
L_k L_h
 =T^epsilon(h,k) L_(h+k mod p),

epsilon(h,k)=floor((h+k)/p).                               (5)
```

Indeed, the ancestry increment on the left is

```text
floor((q+h)/p)
 +floor(((q+h mod p)+k)/p)
 =floor((q+h+k)/p),
```

while the right side contributes

```text
floor((q+(h+k mod p))/p)+floor((h+k)/p).
```

These are equal.  Associativity gives, equivalently,

```text
epsilon(h,k)+epsilon(h+k mod p,l)
 =epsilon(k,l)+epsilon(h,k+l mod p),                       (6)
```

so `epsilon` is the standard normalized carry 2-cocycle.

It is not a coboundary.  Iterating the lifted generator gives

```text
L_1^p=T.                                                    (7)
```

Already after reducing `a mod p`, a homomorphic section
`C_p -> C_(p^2)` would have to send `1` to an element congruent to `1 mod p`
and of additive order dividing `p`.  Every such element is prime to `p` and
has order `p^2`, a contradiction.  Thus the first carry quotient is the
nonsplit extension

```text
0 -> C_p -> C_(p^2) -> C_p -> 0.                           (8)
```

The full extension

```text
0 -> C_(p^d) -> C_(p^(d+1)) -> C_p -> 0                   (9)
```

is likewise nonsplit.  The class in `(8)` is the Bockstein generator.

## 2. Flat realization kills the carry

Let `X` carry a flat `C_p` action `B_h`, so

```text
B_k B_h=B_(h+k mod p).                                    (10)
```

Suppose a map `F` from the lifted state space to `X` intertwines every
natural lift with this flat action:

```text
F L_h=B_h F.                                               (11)
```

Then `(7)` and `(10)` give

```text
F T=F L_1^p=B_1^p F=F.                                    (12)
```

Equivalently, any one unit-carry triangle from `(5)` already forces `(12)`
after cancelling its direct lifted edge.  Therefore exactly one of the
following occurs:

```text
carry-blind contraction:  F T=F;

carry-faithful realization:
  the target must be enlarged by a nontrivial T-fibre.     (13)
```

At the first carry quotient, a faithful deterministic `T`-fibre has at
least `p` states because it contains a free `C_p` orbit.  If the target
residue `q` is also retained, there are `p` fibres and hence

```text
|joint state space|>=p^2.                                  (14)
```

Here the joint realization is required to transport the central `T` action
along the residue action.  This is automatic for an intertwiner of the
natural lift, because `T L_h=L_h T`; hence faithfulness in one residue fibre
forces the same `p`-orbit in every fibre.  Without this joint-equivariance
hypothesis, one isolated faithful fibre would not imply `(14)`.

Equality is attained by `C_(p^2)` itself: every residue fibre is one
`C_p` torsor.  There is still no preferred section of those fibres, by
`(8)` and THM-2611.  Thus the sharp sidecar is a torsor, not a canonical
basepoint.  Full fidelity to `a in C_(p^d)` costs `p^d` states per residue.

This state invoice concerns finite deterministic faithful encodings.  It is
not a lower bound for untyped probabilistic, real, or quantum encodings.

## 3. The exact q3/q11/q7 triangle

Now set `p=13`, `d=5`.  Start at `(a,q)=(a,3)`.  The three relevant lifted
residue arrows are

```text
q3  --h=8--> q11,
q11 --h=9--> q7,
q3  --h=4--> q7.                                          (15)
```

Since `8+9=13+4`, equations `(3)--(5)` give

```text
L_8(a,3)       =(a,11),
L_9 L_8(a,3)   =(a+1,7),
L_4(a,3)       =(a,7),

L_9 L_8=T L_4.                                            (16)
```

Every nontrivial character of the first carry quotient detects this loop:

```text
chi_r(T)=zeta_13^r !=1,       r=1,...,12.                  (17)
```

The affine-address side is flat.  With

```text
S=13^5,       M=13^6,

A_t(y)=y+t S(y-7 mod13) mod M,
beta(q)=2q mod13,
n(q)=6716+beta(q)S mod M,                                 (18)
```

one has, for every `q,h`,

```text
A_(2h)n(q)=n(q+h).                                        (19)
```

In particular,

```text
n(3)=2234474,
n(11)=3348353,
n(7)=378009,

A_3 n(3)=n(11),
A_5 n(11)=n(7),
A_8 n(3)=n(7),
A_5 A_3=A_8.                                              (20)
```

Thus `(16)` has central holonomy one while `(20)` has holonomy zero.
The affine address remembers the target-residue edge and forgets the
semantic ancestry carry.

The `L_8` q3-to-q11 leg in `(15)` is the natural-extension comparison
attached to the physical THM-2847 address edge, and `L_4` is the direct
comparison arrow closing that nerve triangle.  Neither is asserted to be a
449-sheet physical semantic correspondence.  Only the `L_9` leg has the
proved 449-sheet coefficient-support incidence below, and even `L_9` is not
a physical current action on the complete packet.

## 4. The 449-sheet horn attachment is on the carry leg

THM-2835 proves, on the cell `(s,t,clock)=(0,5,1)`, a set `Lambda` of
exactly `449` whole-cylinder labels such that

```text
QB_source(a)=QA_target(q11,a)=1,       a in Lambda.         (21)
```

The beta-selected natural lift is exactly `h=9`, and it sends every label to

```text
(a,11) -> (a+1,7),

QAB_target(q7,a+1)=1,
QA_target(q7,a+1)=QB_target(q7,a+1)=0.                      (22)
```

The cell `(0,5,1)` is the first cell of THM-2847's exact 20-cell E3-only
horn.  On that horn, before the E3 gate, rows `(q0,q3,q11,q7)` and columns
`(origin00,origin12,QA,QAB)` give

```text
A =
[[1,1,  0,  0],
 [1,0,  0,  0],
 [1,1,449,  0],
 [0,0,  0,449]],             det(A)=-449^2.                (23)
```

The E3 projection has rank three and kernel generated by the pure QAB
column.  The complementary macro block sends this generator to
`449 delta_q7`.  Equation `(22)` identifies the support-level attaching leg:
it is the `T`-shifted endpoint of the q11 semantic carry arrow.

The ancestry shift is not forced by the number `449`.  THM-2835 proves that
the corresponding support fillers are carry-blind: after forgetting the
orientation `a -> a+1`, the same cardinality can remain.  Thus `(22)` gives
the canonical natural-lift provenance of the attaching leg, not a recovery
of that provenance from support alone.

This resolves an apparent paradox.  The vector-space sequence in THM-2847
splits over `K_0=Q(zeta_1183)`, while the state-level cyclic extension
`(8)` does not split equivariantly.  There is no ordinary linear `Ext`
obstruction; the nonzero class lives in the realization groupoid, where a
flat address action has forgotten `T`.

## 5. The oriented carry derivative

THM-2839 proves that the horn mask

```text
h_L(Z)=sum_(a in Lambda) Z^a in Z[C_(13^5)]              (23a)
```

has augmentation `449=7 mod13` and is a unit over both
`F_13[C_(13^5)]` and `Q[C_(13^5)]`.  Orient the natural carry edge in
`(22)` and form

```text
partial_T h_L=(Z-1)h_L
 =1_(Lambda+1)-1_Lambda.                                  (23b)
```

The two supports are disjoint: every source label is `156 mod169`, while
every successor is `157 mod169`.  Hence `(23b)` has coefficients in
`{-1,0,1}` and exact signed `l1` norm `898`.

In characteristic thirteen, put `e=Z-1`.  Since `h_L` is a unit,

```text
partial_T h_L=e*(unit),

ord_e(partial_T h_L)=1,

partial_T h_L mod e^2=7e.                                 (23c)
```

Thus the carry derivative generates the augmentation ideal and its image in
the cotangent line `(e)/(e^2)` is nonzero.  Over the complex numbers,

```text
(partial_T h_L)(xi)=(xi-1)h_L(xi).                        (23d)
```

THM-2839 gives `h_L(xi)!=0` for every `13^5`-th root `xi`.  Consequently
the derivative vanishes only at the trivial character and its cyclic
convolution rank is exactly

```text
13^5-1=371292.                                             (23e)
```

Pushing `(23b)` to the first carry residue uses
`156=0 mod13` and `157=1 mod13`:

```text
partial_1 h_L=449(delta_1-delta_0) in Q[C_13].             (23f)
```

It has rank `12`, and all twelve nontrivial carry characters are nonzero.
This is a signed algebraic model on the sharp 13-state carrier from
`(13)--(14)`, not by itself a deterministic sidecar.  Its orientation pairs
two separately positive blocks.
Canon does not supply lawful positive/current access to that pairing while
retaining the fixed source and E3/complement macro factors.  Hence `(23b)`
is a concrete charged reference candidate, not the missing physical
realization.

The mechanism is the same abstract move as the multiplier jets in THM-2842:
an ordinary scalar/augmentation forgets a first-order direction, while one
cotangent multiplier exposes it.  This is a structural bridge, not a
GMC-to-LRC implication.

## 6. Sharp consequence and sharp boundary

Let a proposed carry-to-response map identify the natural semantic lifts
with the flat affine address lifts:

```text
F L_h=A_(2h) F.                                           (24)
```

Then `(12)` or the specific triangle `(16)--(20)` forces

```text
F T=F.                                                     (25)
```

Therefore such a map cannot simultaneously be:

```text
flat-address equivariant,
and
faithful to the ancestry orientation a -> a+1 selecting
the natural-lift provenance in (22).                        (26)
```

To retain that carry at its first residue, one must add a 13-state fibre.
Together with the 13 target/address states, the sharp joint deterministic
state count is `169`; equality is the nonsplit `C_169` torsor.  Choosing an
actual identification with a response orbit additionally needs one
basepoint, or a zero-holonomy family of local basepoints, exactly as in
THM-2611.

This conclusion is conditional on using the natural-lift strategy.  It does
not prove that every conceivable physical horn contraction must expose this
carry.  A non-equivariant construction or a different physical transporter
could in principle bypass `(24)`.

The sharp carry-blind hostile is the projection

```text
F_0(a,q)=q.                                                 (27)
```

It satisfies `F_0 L_h=B_h F_0`, identifies both sides of `(16)`, and retains
all three target labels while forgetting `T`.  One may also retain the
printed QA/QAB support names and the number `449`; those data still do not
recover whether the q7 label arose from `a` or `a+1`.  Hence the theorem
locates a provenance/basepoint obstruction, not a support obstruction.

Nor does the abstract `C_169` lift provide the missing LRC current.  The
natural semantic arrow in `(22)` does not lawfully co-shift the fixed source
and macro E3 factor; THM-2835's literal physical transport obstruction
remains.  The theorem proves the missing coordinate and its minimal state
cost, not its physical realization.

## 7. Connection contract and next decisive test

```text
source:
  THM-2847 q3/q11 transverse address edge and 20-cell E3 horn;
  THM-2835 449-sheet QA(q11,a)->QAB(q7,a+1) carry leg;

target:
  the first ancestry-carry quotient C_169 -> C_13;

map:
  compare L_9 L_8 with L_4, and compare A_5 A_3 with A_8;

preserved:
  q3/q11/q7 labels, the exact (0,5,1) horn cell, 449-sheet
  QA/QAB multiplicity, affine addresses, and E3-kernel attachment;

destroyed by the flat address contraction:
  the central ancestry translation T and all twelve of its charged
  C_13 characters;

needed sidecar:
  lawful physical/current access to the signed derivative (23b), or an
  equivalent C_13 carry fibre over the address bank, plus a selected
  carry-to-response basepoint/zero-holonomy section and a target-relative
  macro transporter;

cheapest decisive test:
  refine the (0,5,1) endpoint current by a mod-13 ancestry character,
  apply the a->a+1 carry leg, and test whether one nontrivial carry
  character survives while the fixed source and E3/complement macro
  factors are retained.  Total cancellation records the sharp hostile;
  a survivor supplies the missing polarized macro-truth reference.
```

The q/t graph of THM-2847 has only thirteen states.  Equations `(14)` and
`(17)` explain why another Fourier transform of that graph cannot discover
the carry: it needs a second fibre coordinate, not another character name
on the same graph.

No scalar row is excluded.  The ledger remains `165`, and LRC(14) remains
open.

## 8. Exact evidence

The companion pins the promoted THM-2835, THM-2839, and THM-2847
script/output hashes and checks their exact 449-sheet, full-spectrum,
20-cell, and rank-one mapping-cone statements.  Independently it checks:

- `6,591` natural-lift composition instances;
- all `2,197` carry-cocycle triples;
- the generator winding `L_1^13=T`;
- all thirteen candidate sections of `C_169 -> C_13`;
- all `169` target/address semiconjugacies;
- the q3/q11/q7 semantic and affine triangles;
- all twelve nontrivial first-carry characters;
- the rank-12 reduced derivative, full rank `371292`, first
  augmentation-adic order, and signed norm `898`; and
- the determinant, ranks, kernel, and complementary attachment in `(23)`.

Run

```text
python 04-computation/lrc14_q3_q11_q7_bockstein_holonomy_thm2851.py
python -O 04-computation/lrc14_q3_q11_q7_bockstein_holonomy_thm2851.py
```

Both modes must byte-match the stored output after LF normalization.

The independent hostile audit rederived the general carry cocycle,
`L_1^13=T`, the nonsplit `C_169 -> C_13` extension, the flat-intertwiner
identity `FT=F`, and the sharp fibrewise `13/169` state invoice.  It checked
that `L_8,L_4` are comparison-nerve arrows, while only `L_9` has the pinned
449-sheet coefficient-support incidence and even that incidence is not a
physical current action.  It independently reconstructed the oriented
derivative, its disjoint `156/157 mod169` supports, signed norm `898`,
cotangent class `7e`, ranks `12` and `371292`, and all scope boundaries.
Normal, optimized, and stored transcripts agreed exactly.  No remaining
mathematical or evidence defect was found.
