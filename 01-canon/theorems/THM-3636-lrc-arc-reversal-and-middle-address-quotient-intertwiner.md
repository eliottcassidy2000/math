---
id: THM-3636
title: "LRC arc reversal and middle-address quotient intertwiner"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On the pinned six-point coefficient
  module, adjoining native directed-arc reversal to THM-3625's split address
  algebra resolves both doubled character multiplicities.  The restricted
  address projector onto the middle two points identifies
  R4sharp/K2sharp with a two-plane and LambdaSharp factors through it by an
  invertible matrix.  This is static finite-field algebra, not chronology,
  current, characteristic-zero transport, row exclusion, or LRC(14).
source: kps-s189 + agents Dalton/Descartes / THM-3625 multiplicity-debt continuation, 2026-08-21
audit: >
  PASS -- an independent reconstruction recovered S, Pi_W, Pi_K, the complete
  character-corner table, the eight-dimensional crossed algebra and its
  F_p^2 + M_2(F_p) + F_p^2 type.  It independently refuted the endpoint map,
  recovered the middle quotient, exact LambdaSharp factor and determinant,
  and matched every maintained operator/space digest.  The refutation is an
  explicit hostile in the corrected theorem, not a surviving claim.
depends_on:
  - THM-3625-lrc-point-by-branch-split-four-character-address-algebra
  - THM-3615-lrc-pointed-root-difference-lift-flag
related:
  - THM-3602-lrc-centered-a4-flag-inside-pointed-six-carrier
script: 04-computation/lrc_arc_reversal_middle_address_quotient_thm3636.py
output: 05-knowledge/results/lrc_arc_reversal_middle_address_quotient_thm3636.out
script_sha256: 796a80daae5401c37f280ed97d7c80b40e47c06c23c3d03a844f38f797595001
output_sha256: c96d5b7b107c6b67a0500020727f349478908baffcce87ce336642c316ede82a
semantic_sha256: 6d230ba65767bdb7080f86d391213c9ef0e6d0fcb7edf3b7e2a7ea6629bd5c73
hash_basis: raw LF bytes for files; canonical JSON for semantic ledger
---

# THM-3636 -- arc reversal resolves the address multiplicity debt

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**  This restores the native reversal of
the six directed point/root carriers.  It gives a precise static quotient
and repairs one tempting but false physical interpretation of THM-3615's
rank-two augmentation.

## 1. Exact coefficient flag

Work over the pinned field

```text
F=F_p,                  p=755373809845391722745761.      (1)
```

Use the six-point order

```text
(0,0),(1,0),(1,6),(3,6),(3,12),(2,12).                 (2)
```

In the raw/centered coefficient frames reconstructed from THM-3625 and
THM-3615, respectively, the inherited spaces have the exact forms

```text
K=K2sharp
 =span{e_0+e_1,e_4+e_5},

R=R4sharp
 =span{e_0+e_1,e_2,e_3,e_4+e_5},

W=P6sharp intersect Q6sharp
 =span{e_0,e_1,e_4,e_5}.                               (3)
```

The equality of the two displayed coordinate descriptions of `K` is checked
after independently lifting the inherited marginal flag; it is not assumed
from point labels.

## 2. Arc reversal and the crossed algebra

Let

```text
S=(0 1)(2 3)(4 5)                                     (4)
```

act on coefficient rows.  Thus `S^2=I`.  If `A` is THM-3625's commutative
four-dimensional address algebra, exact multiplication gives

```text
dim_F <A,S>=8,
S A S=A,
<A,S> is product-closed and noncommutative.             (5)
```

Write the four address-character spaces as

```text
P6sharp=E_0 direct_sum E_1 direct_sum E_2 direct_sum E_3,
(dim E_0,dim E_1,dim E_2,dim E_3)=(2,1,1,2).           (6)
```

The involution preserves `E_0,E_3`, splitting each into one `+` and one `-`
line, and exchanges `E_1,E_2`.  The restriction spans of the crossed algebra
on

```text
E_0,                  E_1 direct_sum E_2,             E_3
```

have dimensions

```text
2,                              4,                     2. (7)
```

The address spectral projectors separate the three blocks, so `(5)--(7)`
identify the algebra abstractly as

```text
<A,S> isomorphic to F^2 direct_product M_2(F) direct_product F^2. (8)
```

This is a split statement over the one pinned field, not a
characteristic-zero representation theorem.

The endpoint projector and rigidity projector are

```text
Pi_W=diag(1,1,0,0,1,1) in A,
Pi_K=Pi_W (I+S)/2.                                     (9)
```

Then `Pi_K` is a central rank-two idempotent of the crossed algebra and

```text
im Pi_K=K,                  W intersect Fix(S)=K.       (10)
```

Thus the added involution supplies exactly the two signs that the commutative
address algebra could not see.

## 3. The middle projector gives the exact quotient

Put

```text
Pi_M=I-Pi_W=diag(0,0,1,1,0,0).                         (11)
```

This projector already belongs to `A`.  In the power basis
`(I,A_3,A_3^2,A_3^3)` its coordinates are

```text
(231956577553530055552921,
 258341843836556363678220,
 526705381605750340113858,
 221092595255266945805223).                            (12)
```

On the whole six-space, `ker Pi_M=W`, so THM-3625's global no-go for selecting
`K` by an address-algebra kernel is untouched.  On the restricted four-space
`R`, however,

```text
ker(Pi_M|R)=K,
im(Pi_M|R)=span{e_2,e_3}.                              (13)
```

Equivalently, the middle-coordinate map

```text
delta_mid:R -> F^2,              delta_mid(c)=(c_2,c_3) (14)
```

induces the exact isomorphism

```text
R/K isomorphic to F^2.                                  (15)
```

This resolves a *restricted quotient*; it does not contradict the inability
of `A` to isolate `K` inside all of `P6sharp`.

## 4. LambdaSharp factors through the middle quotient

Let `A2sharp` be the canonical rank-two augmentation image from THM-3615,
with the canonical row basis used by the exact companion.  If row vectors act
on the right, there is an exact factorization

```text
LambdaSharp|R = delta_mid T_mid,                        (16)
```

where

```text
T_mid=
 [126498113787680818370196  646561219255993169342961]
 [384727693242765231857657  452865069702498262026030], (17)

det T_mid=275730381649850587765623 !=0 mod p.           (18)
```

The equality is checked on an ordered basis of `R` before any rowspace
canonicalization.  This ordering point is load-bearing: canonicalizing the
domain without applying the same row operations to the images creates a
false failed factorization.  Equations `(13),(18)` also recover

```text
ker(LambdaSharp|R)=K                                   (19)
```

by an independent coordinate route.

The two augmentation rows in physical difference/relation coordinates have
digests

```text
profile:  4af0bf56a78f56ef0af7d3cd6ff6d11b046d43421dc2e6006158b70616e3ede6,
expanded: 54e969e70b3b42d3ad5ecd4601ee63451bdfbcec18f368ef7f975b744d223ea3. (20)
```

## 5. The endpoint-imbalance interpretation is false

The geometrically tempting endpoint map is

```text
delta_end(c)=(c_0-c_1,c_4-c_5).                        (21)
```

But `(3)` gives the stronger hostile

```text
delta_end|R=0,                     rank(delta_end|R)=0. (22)
```

Since `LambdaSharp|R` has rank two, it cannot factor through `(21)`.  The
successful map `(14)` reads the middle directed-arc pair, not endpoint
imbalance.  This correction agrees with the older third-digit response atlas,
where the full-parent arc-reversal defect is localized on the same middle
pair while both endpoint pairs have zero defect.

For completeness, endpoint differences do have a lawful domain:

```text
W=K direct_sum span{e_0-e_1,e_4-e_5},
W intersect R=K.                                       (23)
```

Thus `(21)` identifies `W/K`, while `(14)` identifies `R/K`; confusing the
two transverse complements is exactly the failed implication caught by the
hostile audit.

## 6. Verification and strict boundary

The deterministic companion pins the promoted THM-3625 and THM-3615 theorem,
script, and output triples, rebuilds the complete branch-resolved tensor, and
checks `(3)--(22)`.  It also verifies the exact operator and augmentation
digests, contains no Python `assert` node, and executes the same gates under
normal and optimized Python.

Reproduce with

```bash
python3 04-computation/lrc_arc_reversal_middle_address_quotient_thm3636.py
python3 -O 04-computation/lrc_arc_reversal_middle_address_quotient_thm3636.py
```

Every conclusion is static over one finite field.  In particular, no child
space in the existing third-digit bank is arc-invariant, and the old
`68 -> 78` quotient growth still forbids calling this projector a digit
transition.  The theorem supplies no chronology, source/arrival typing,
physical current, characteristic-zero lift, row exclusion, or proof of
LRC(14).  The independent hostile audit closes the promotion gate without
changing that boundary.  **QED.**
