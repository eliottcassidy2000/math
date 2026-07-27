---
id: THM-2611
title: "Principal C13 bibundle lift torsor and holonomy section obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  A deterministic faithful lift of a quotient which has forgotten a free
  C13 orbit needs at least thirteen sidecar states, and equality makes the
  sidecar a C13 torsor.  Equivariant identifications of two fixed-action
  C13 torsors form another C13 torsor: there are exactly thirteen, with no
  canonical member.  Around a cyclic atlas, the accumulated kernel
  translation is the complete obstruction to a parallel section; when it
  vanishes there are exactly thirteen sections.  For C91 -> C7 the thirteen
  lifts of one positive clock edge are delta(c)=78+14c.  This recovers
  THM-2607's gauge count, explains THM-2601/2608's thirteen root-to-section
  references, and proves that THM-2610's chronological quotient needs a
  genuine thirteen-sheet ancestry bibundle before its old and future C13
  actions can be compared.  No such physical bibundle or selected section
  is constructed, and no scalar row or LRC(14) conclusion follows.
source: deep-energy-audit-2026-07-28-principal-deck-holotopy
depends_on:
  - THM-2601-linear-bockstein-sheet-coordinate-and-nonlinear-target-monodromy
  - THM-2607-constant-six-rail-boundary-holonomy-invoice
  - THM-2608-alternative-rail-clock-collapse-and-missing-transition-index
  - THM-2610-chronological-paired-slice-marked-triangle-graft-and-action-axis-boundary
related:
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
script: 04-computation/lrc14_principal_c13_bibundle_holotopy_thm2611.py
output: 05-knowledge/results/lrc14_principal_c13_bibundle_holotopy_thm2611.out
script_sha256: 866b1eecff0e220ba1bdaf695d5dee32d382a857f261988a34d8122948aec209
output_sha256: 4ef1c37448ee4ad5e819a04200ecb30603783d7acfa0d770160aab3f63bea829
hash_basis: LF-normalized bytes
---

# THM-2611 -- the missing chronological object is a principal `C_13` bibundle

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The number thirteen occurs three times at the current root-holotopy
frontier.  THM-2607 has thirteen gauges after its rail correction kills the
seven-clock holonomy.  THM-2601, as resolved explicitly in THM-2608, admits
thirteen equivariant deep-root-to-target references.  THM-2610 preserves a
marked triangle beside a later full shift spectrum, but positive Koopman
time has forgotten the old root deck.

These are one finite-extension mechanism.  A quotient has collapsed a
principal `C_13` fibre.  Recovering the fibre is a thirteen-state problem;
identifying two recovered fibres is itself a thirteen-choice torsor; and a
cycle of local identifications has a global section precisely when its
kernel holonomy is zero.  This theorem proves that classification.  It does
not assert that any of the thirteen abstract choices is a lawful physical
ancestry map.

## 1. Faithful recovery of a collapsed orbit costs thirteen states

Let `K=C_p` act freely and transitively on a finite set `P`; thus `P` is a
`K`-torsor and `|P|=p`.  Suppose a quotient observable

```text
F:P -> B
```

is `K`-invariant.  On this orbit it is constant, so it cannot distinguish
any two source states.

A deterministic sidecar recovery is a finite `K`-set `A` and a
`K`-equivariant map

```text
j:P -> A                                                   (1)
```

which is injective on the collapsed fibre.  Then

```text
|A|>=p.                                                    (2)
```

Indeed, if `k` fixes `j(x)`, equivariance gives
`j(kx)=k j(x)=j(x)`.  Injectivity and freeness of `P` force `k=0`.
Therefore the orbit of `j(x)` is free and has `p` elements.  Equality in
(2) says that all of `A` is this orbit, so `A` is itself a `K`-torsor.
Conversely `A=P` and `j=id` attain equality.

Thus thirteen is the sharp finite deterministic state cost when `p=13`.
The statement is deliberately about faithful finite `K`-sets.  It is not a
lower bound for untyped real, probabilistic, or quantum encodings.

## 2. Equivariant fibre identifications form a torsor, not a point

Let `P,Q` be `K`-torsors with the action generator fixed.  Every
`K`-equivariant map

```text
phi:P -> Q                                                 (3)
```

is a bijection.  Choosing one `x_0 in P`, the value `phi(x_0)` may be any
of the `p` points of `Q`, and then equivariance determines every other
value.  Hence

```text
|Iso_K(P,Q)|=p.                                            (4)
```

Postcomposition by translations of `Q` is free and transitive on (4).
Precomposition by inverse translations of `P` gives the commuting action
from the other side.  Thus `Iso_K(P,Q)` is the principal `(K,K)` bibundle
of fixed-action identifications.  It has no preferred point without one
additional basepoint or transport rule.

This fixed-action clause is load-bearing.  If `K=F_p` and one also permits
an arbitrary automorphism of the action generator, the maps are

```text
phi_(kappa,c)(r)=kappa r+c,
kappa in F_p^*, c in F_p,                                 (5)
```

so there are `p(p-1)`, not `p`.  Neither family becomes canonical merely
because both endpoint label sets have been printed as `F_p`.

## 3. A cyclic local system has zero holonomy or no section

Let `P_0,...,P_(n-1)` be `K`-torsors and let

```text
tau_i:P_(i-1) -> P_i,                 i mod n,             (6)
```

be fixed-action equivariant bijections.  Their cyclic composite is an
equivariant automorphism of `P_0`, hence translation by a unique

```text
Hol(tau) in K.                                             (7)
```

A parallel cyclic section is a tuple `(x_i)` satisfying

```text
x_i=tau_i(x_(i-1))                         for every i.    (8)
```

Such a section exists if and only if

```text
Hol(tau)=0.                                                (9)
```

If (9) holds, every one of the `p` choices of `x_0` extends uniquely, so
there are exactly `p` sections.  If it fails, translation by a nonzero
element of the free torsor has no fixed point, so there are none.

Choose temporary origins in the `P_i` and write the transition coordinates
as `g_i in K`.  Then

```text
Hol(tau)=sum_i g_i.                                        (10)
```

Changing origin `i` by `h_i` changes `g_i` by the coboundary
`h_(i-1)-h_i`; the sum (10) is unchanged.  Equations (9)--(10) are the
finite torsor form of the Cech calculation, with its equality case included.

## 4. Exact `C_91 -> C_7` normal form

For distinct primes `p,q`, put `E=Z/(pq)`.  The CRT idempotents

```text
e_p=q(q^(-1) mod p),              e_q=p(p^(-1) mod q)     (11)
```

identify the kernel fibre `C_p` and quotient `C_q`.  In the live case
`(p,q)=(13,7)`,

```text
iota(r)=14r,
pi(z)=z mod 7,
e_p=14,                            e_q=78.                 (12)
```

Every fibre of `pi` is a principal `C_13` set.  The complete fibre of
positive-clock edge lifts is

```text
pi^(-1)(1)={delta(c):c in F_13},
delta(c)=78+14c mod 91.                                  (13)
```

There are exactly thirteen choices.  Multiplication by thirteen has kernel
`iota(F_13)` and sends every member of (13) to the same quotient edge `13`:

```text
13 delta(c)=13 mod 91,
13 iota(r)=0 mod 91.                                     (14)
```

So positive Koopman time forgets exactly the fibre which distinguishes the
thirteen lifts.

For a seven-edge loop with kernel corrections `c_ell`, (13) gives

```text
sum_ell delta(c_ell)=iota(sum_ell c_ell).                 (15)
```

Adding the constant marker correction `a` on each edge changes the
right-hand side to

```text
iota(7a+sum_ell c_ell).                                  (16)
```

By Section 3 the lifted root local system has a section exactly when (16)
is zero; then it has exactly thirteen.

## 5. The three frontier appearances

### 5.1 THM-2607: the gauges are precisely fibre origins

THM-2607's constant-six rails have

```text
c_ell in {0,7},                 n=#{ell:c_ell=7}.          (17)
```

Equation (16) is `iota(7(a+n))`, so it vanishes exactly for

```text
a=-n mod 13.                                               (18)
```

After the quotient clock vertex is fixed and the theorem's explicitly
conditional rail/deck identification is granted, the thirteen sections in
Section 3 are exactly THM-2607's thirteen kernel-fibre origins/additive
vertex gauges.  The theorem does not strengthen that conditional grant.

### 5.2 THM-2601/2608: thirteen references are the splitting torsor

THM-2601's target successor `S` is one thirteen-cycle.  Give the physical
root set its `+1` action and the target-sheet set its `S` action.  The
fixed-action equivariance equation

```text
phi(r+1)=S(phi(r))                                         (19)
```

has exactly the thirteen solutions

```text
phi_c(r)=S^r(c),                    c in F_13.             (20)
```

This is Section 2 with `P` the root torsor and `Q` the target-sheet torsor.
The exact companion independently checks the explicit successor from
THM-2601 and all thirteen maps.  Proved THM-2608 finds the same family and
correctly treats it as an unresolved reference torsor.  Summing over all its
points is gauge-invariant algebra; it does not select a physical
root-to-target map.

### 5.3 THM-2610: chronology supplies the quotient, not the bibundle

THM-2610 has an old root/target action and a separately installed future
paired-shift action.  On the primal state spaces both are regular
`C_13` actions.  But for every `L>=1`,

```text
T^L(x+r/13)=T^Lx.                                         (21)
```

Thus the chronological map factors through the quotient which is constant
on the old `C_13` orbit.  Section 1 proves that any finite deterministic
faithful lift of that orbit must add at least thirteen ancestry states.
At the minimum it is a `C_13` torsor, and Section 2 says that identifying it
with the future shift torsor still has thirteen fixed-action choices.

THM-2610's nonzero future characters prove coexistence on a physical
two-time packet, but a character is not a chosen point of this bibundle.
Numerically setting a future colour equal to an old residue likewise does
not choose an equivariant identification.  A genuine positive result must
construct a common physical ancestry fibre, its two commuting actions, and
one lawful section or holonomy-zero family of local sections.

## 6. Equality controls and scope

All hypotheses above are sharp.

1. With only twelve sidecar states no free `C_13` orbit exists, so faithful
   equivariant recovery is impossible.  Thirteen regular states attain the
   bound.
2. Forgetting the fixed generator enlarges thirteen translations to 156
   affine semilinear identifications; it does not select one.
3. On a cycle, zero holonomy gives all thirteen parallel sections, while a
   single nonzero edge translation gives none.
4. A coefficient sum over all thirteen candidate identifications can be
   gauge invariant without being the coefficient of any one nonnegative
   common-ancestry carrier.  The theorem is a classification and lower
   bound, not a positivity transfer.

No physical `C_13` bibundle is produced here.  In particular there is no
same-event next target-section index, adjacent-clock carrier, semantic
owner arrow, completed relation current, scalar-row exclusion, or proof of
LRC(14).  The live ledger remains `165`.

## 7. Exact companion

Run

```text
python 04-computation/lrc14_principal_c13_bibundle_holotopy_thm2611.py
python -O 04-computation/lrc14_principal_c13_bibundle_holotopy_thm2611.py
```

The dependency-free exact companion verifies the `C_91` exact sequence,
all seven thirteen-point fibres, all thirteen positive-edge lifts and their
common Koopman image, the thirteen fixed-action and 156 semilinear affine
maps, the complete 1,664 binary seven-cycle/marker pairs and their section
counts, and the thirteen explicit THM-2601 successor intertwiners.  Every
check uses exact integer arithmetic and an explicit optimized-mode guard.

Normal and optimized runs byte-match the stored transcript after LF
normalization.

QED.

The independent hostile audit rederived the finite deterministic
equivariant minimum, including its equality case; checked that
`Iso_K(P,Q)` is a bitorsor rather than a free `K x K`-set; separated the
thirteen fixed-action maps from the 156 semilinear maps; and independently
verified the cyclic holonomy criterion, the `C_91` formulas, every
application boundary, both optimized and normal transcripts, and the
declared LF-normalized hashes.  It found no physical-ancestry promotion
hidden in the proof.
