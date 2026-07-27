---
id: THM-2567
title: "Deep-coloured duty-replica cycle and augmentation cancellation"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  Every nonzero rational nonnegative C_13^3 table with the lawful paired
  diagonal zero has a nonzero physical target-zero anchor in every nonzero
  deep colour.  The canonical duty commutator therefore forces six equal
  nonzero target replicas in each colour.  Exact target-null circulation
  cancels the entire family after summing deep colour, so the construction
  is a coloured cycle rather than an uncoloured mixed face.  The sharp
  nonnegative hostile 1_(s!=0)1_(r!=t) realizes target-zero value -12/169
  and all first-target nonzero values 1/169 in every nonzero deep colour,
  with exact augmented cancellation.  THM-2563 supplies the hypotheses, but
  its partial-bare target and THM-2562's duty carrier are not yet one lawful
  deep-refined physical carrier.  No row exclusion or LRC(14) conclusion.
source: root-holotopy-2026-07-27-coloured-duty-cycle
depends_on:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2562-canonical-duty-commutator-line-rank-and-anchor-rigidity
  - THM-2563-paired-dipole-deep-target-corner-and-partial-bare-boundary
related:
  - THM-2556-reynolds-duty-curvature-and-fibre-covariance-mixed-cell
  - THM-2564-six-tooth-doubly-centered-tomography-and-section-holonomy-boundary
script: 04-computation/lrc14_deep_coloured_duty_cycle_thm2567.py
output: 05-knowledge/results/lrc14_deep_coloured_duty_cycle_thm2567.out
script_sha256: 5c0cb3eb5dd6608dfa7fbb3e25c52dd244d179155d9ea8af44efc394c75a2e20
output_sha256: 9098db9b53bc65195e7d4f5cf0617982440b86153952d060fd75feffbe5c6f97
hash_basis: working-tree bytes (LF)
---

# THM-2567 -- the duty replicas form a deep-coloured cycle

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT
PENDING.**

THM-2563 supplies exactly the paired three-coordinate geometry that the
THM-2562 duty calculation did not yet see.  The composition has a clean but
surprising outcome.  Every deep colour separately carries a rigid nonzero
six-replica duty face, while the complete family cancels exactly after the
deep colour is forgotten.

The result is therefore a finite holotopy statement: all coloured faces are
nonzero, but their augmentation is a boundary.  A deep-refined common carrier
is load-bearing; an uncoloured marginal cannot see the obstruction.

## 1. Lawful paired table and physical target coordinates

Put `p=13`.  Let

```text
H:F_p^3 -> Q_(>=0),                  H not identically zero,  (1)
```

with coordinates `(r,s,t)`, and assume the lawful paired diagonal zero

```text
H(t,s,t)=0                    for every s,t.                  (2)
```

Use the normalized transform

```text
B(m,b,h)
 =p^(-3) sum_(r,s,t) H(r,s,t) zeta^(mr+bs+ht).               (3)
```

THM-2365's physical target is not the raw pair `(b,h)`.  Its coordinates are

```text
q=(q_1,q_2)=(b,m+h).                                        (4)
```

For a fixed deep colour `m`, define the physical target density

```text
A_m(q_1,q_2)=B(m,q_1,q_2-m).                                (5)
```

In particular, the target-zero coefficient in deep colour `m` is the
inverse-character slice

```text
A_m(0,0)=B(m,0,-m),                                         (6)
```

not `B(m,0,0)`.

## 2. Every nonzero deep colour has a target-zero anchor

For every `m in F_p^*`,

```text
A_m(0,0)!=0.                                                (7)
```

### Proof

Set `u=r-t` and aggregate

```text
G(u)=sum_(s,t) H(u+t,s,t).                                  (8)
```

The profile `G` is rational, nonnegative, and nonzero.  Equation (2) gives

```text
G(0)=0.                                                     (9)
```

Changing variables in (6) gives

```text
p^3 A_m(0,0)=sum_u G(u) zeta^(mu).                          (10)
```

If (10) vanished for one `m!=0`, the rational coefficient polynomial of
`G` would be divisible by `Phi_13`.  Its degree is at most twelve, so all
thirteen coefficients would be equal.  The zero coefficient (9) would then
make them all zero, contradicting `H!=0`.  Galois conjugacy proves (7) for
all twelve nonzero colours. QED.

Nonnegativity is used only to infer `G!=0` from `H!=0`; a signed table could
cancel in (8).  Rationality is essential to the cyclotomic step.

## 3. Six canonical duty replicas in every colour

Let `n(q)` be the unnormalized canonical duty of THM-2562 and, for a gain
`g!=0`, let

```text
C_g=sum_(j=0)^6 tau_g^j,

K_g=[diag(n),C_g].                                          (11)
```

THM-2562 equation (14) applies to every complex target density:

```text
(K_g a)(jg)=-d_g a(0),               j=1,...,6,             (12)
```

where

```text
d_g
 =229692,       g lies on a coordinate axis,
 =440232,       g_x g_y!=0 and g_x+g_y=0,
 =440244,       otherwise.                                  (13)
```

Applying (12) to (5) and using (7) gives, for every `m!=0`, every `g!=0`,
and every `j=1,...,6`,

```text
(K_g A_m)(jg)=-d_g A_m(0,0)!=0,                             (14)

||K_g A_m||_2>=sqrt(6)d_g |A_m(0,0)|.                       (15)
```

Thus the paired table has a rigid six-replica duty face in **every**
nonzero deep colour.  No selection over `m` is paid.

## 4. Deep augmentation cancels all replicas

The same diagonal law (2) gives the exact target-null circulation

```text
sum_(m in F_p) A_m(q)=0             for every q in F_p^2.    (16)
```

Indeed, summing (5) over `m` inserts

```text
sum_m zeta^(m(r-t))=p 1_(r=t),                               (17)
```

and the surviving diagonal is zero.  Since `K_g` is linear,

```text
sum_m K_g A_m=0                    for every g!=0.            (18)
```

Equations (14) and (18) are compatible only because the deep-zero colour
supplies the exact opposite augmentation.  The uncoloured target table sees
zero even though all twelve nonzero deep-coloured faces are rigidly nonzero.

This is the coloured duty-replica cycle.  The deep character is not optional
decoration; it is the coordinate carrying the nontrivial class.

## 5. Sharp nonnegative hostile

The exact table

```text
H(r,s,t)=1_(s!=0) 1_(r!=t)                                (19)
```

has mass `1872`, is rational and nonnegative, and has both the `s=0` plane
and `r=t` diagonal zero.  Direct character summation gives, for every
`m!=0`,

```text
A_m(0,0)=-12/169,

A_m(b,0)=1/169                 for every b!=0,

A_m(b,q_2)=0                  for q_2!=0.                    (20)
```

The deep-zero colour has

```text
A_0(0,0)=144/169,

A_0(b,0)=-12/169              for b!=0,                     (21)
```

so (16) is coefficientwise exact.  Consequently the twelve nonzero colours
each have replicas

```text
(K_g A_m)(jg)=12d_g/169,
```

while the zero colour has `-144d_g/169`; their sum is zero.

The additional `r=0` plane in THM-2563 does not remove the mechanism.  The
stronger table

```text
1_(r!=0)1_(s!=0)1_(r!=t)                                  (22)
```

has all three zero loci, mass `1728`, and all twelve inverse-character
anchors nonzero.  The companion verifies this independently of (19).

## 6. Exact application and the carrier boundary

THM-2563's table `R_(r,s,t)` is rational, nonnegative, nonzero, and obeys

```text
R_(0,s,t)=0,
R_(r,0,t)=0,
R_(t,s,t)=0.                                                (23)
```

Therefore Sections 2--4 apply formally: every nonzero deepest colour has a
nonzero inverse-character target-zero anchor and a six-replica canonical
duty face; the complete deep augmentation cancels.

This does **not** yet identify a physical THM-2556 mixed curvature.  There
are two distinct typing debts:

1. THM-2563 retains the moving-endpoint residue but not the fixed old-head
   residue, so its `q` is still partial-bare rather than the completed
   left-minus-right THM-2334 target;
2. THM-2562's `K_g` lives on the canonical quotient-duty/Reynolds carrier.
   No theorem yet refines that carrier to the same deep colour and endpoint
   atom as `R`.

Thus (14) is an exact algebraic face on the partial-bare target coordinate,
and (18) is the exact obstruction to forgetting deep colour.  Calling their
composition a physical mixed face before constructing the common carrier
would repeat the cross-base error recorded in MISTAKE-281.

## 7. Holotopy interpretation and comparison with tomography

The family `(A_m)_m` lies in the kernel of deep augmentation, while every
nonzero colour has nonzero duty boundary.  In chain language,

```text
deep-coloured faces --augmentation--> uncoloured face
      nonzero                              zero.              (24)
```

This is a cycle created by a missing common coloured 2-cell, not by lack of
Fourier support.

THM-2564 cannot remove it.  Any six tooth profiles may be prescribed
independently in its signed codomain, including the same six-replica vector
in every chart.  Its tooth slope also connects owner phase to intervention
co-shift, not deep colour to the physical target (4).  The correct next map
must refine the Reynolds/duty square on one common deep-coloured endpoint
carrier, not add another signed tomography bank.

## 8. Exact companion and stopping boundary

Run

```bash
python3 04-computation/lrc14_deep_coloured_duty_cycle_thm2567.py
python3 -O 04-computation/lrc14_deep_coloured_duty_cycle_thm2567.py
```

Both executions must reproduce

```text
05-knowledge/results/lrc14_deep_coloured_duty_cycle_thm2567.out
```

byte-for-byte.  The dependency-free referee computes the complete transform
of (19) under the physical reindexing (4), checks all `169` circulation
identities, all `168*12*6` nonzero-colour replica rows, and all
`168*169` augmented duty cancellations.  It also checks the three-zero-locus
hostile (22) and all `2,028` diagonal-free singleton controls.  There are
`67,233` explicit assertions.

The theorem proves that the paired-dipole corner has a deep-coloured duty
cycle and isolates exact augmentation cancellation.  It does not build the
deep-refined Reynolds carrier, recover the fixed left residue, choose an
endpoint colour or future root, produce a Hall edge, exclude a scalar row,
or prove LRC(14). **QED.**
