---
id: THM-3654
title: "LRC fixed-branch rigidity eigencriterion"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.
  On THM-3636's four-dimensional relation space, literal restriction of the
  centered six-point tensor to the fixed branch r=6 has alpha eigenspace
  exactly the rank-two rigidity plane and beta eigenspace exactly the middle
  coordinate plane.  Its alpha-defect realizes LambdaSharp, up to the pinned
  invertible middle-coordinate change.  This is static finite-field algebra,
  not chronology, current, characteristic-zero transport, row exclusion, or
  LRC(14).
source: kps-s191 / THM-3647 literal-branch continuation, 2026-08-21
depends_on:
  - THM-3647-lrc-single-reversal-paired-branch-spectral-projector
related:
  - THM-3636-lrc-arc-reversal-and-middle-address-quotient-intertwiner
  - THM-3625-lrc-point-by-branch-split-four-character-address-algebra
script: 04-computation/lrc_fixed_branch_rigidity_eigencriterion_thm3654.py
output: 05-knowledge/results/lrc_fixed_branch_rigidity_eigencriterion_thm3654.out
script_sha256: 544736158db5e05ff595a50a2ed67741c2f3030cdc9eed33307a1515efe5515b
output_sha256: 8405640656b510b43325999ab74845bd933c2c6047ca1667a1d97f920a0a53c4
semantic_sha256: 4cab3b1ef03eebf75afcd7f968bf2e1de7729711de3c519261b645fc60000369
hash_basis: raw LF bytes
---

# THM-3654 -- a literal fixed branch detects the LRC rigidity plane

**PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE
AUDIT.**  Unlike THM-3647's abstract spectral reconstruction, this theorem
returns to the branch-resolved tensor and identifies the corresponding
literal restriction map.  The restriction remains a static observable.

## 1. The inherited flag and fixed address

Work over

```text
F=F_p,                  p=755373809845391722745761.     (1)
```

In the six-point coefficient coordinates of THM-3636, put

```text
K=span{e_0+e_1,e_4+e_5},
R=span{e_0+e_1,e_2,e_3,e_4+e_5},
M=span{e_2,e_3}.                                      (2)
```

Thus `dim K=2`, `dim R=4`, and `R=K direct_sum M`.  Let `A_6` be the
singleton fixed-branch address operator of THM-3647 and let

```text
alpha=485148529749842658283160,
beta =533745727702934015947346.                       (3)
```

THM-3647 gives, with row vectors acting on the right,

```text
A_6=alpha Pi_W+beta Pi_M,
Pi_W=diag(1,1,0,0,1,1),
Pi_M=I-Pi_W.                                          (4)
```

Consequently

```text
A_6-alpha I=(beta-alpha)Pi_M.                         (5)
```

The two scalars in `(3)` are distinct modulo `p`.

## 2. Literal branch restriction intertwines the address

Let `G_0,...,G_5` be the six centered point generators reconstructed from
the parent branch tensor, and let `G_{i,6}` be the centered profile obtained
by restricting point `i` to the literal branch digit `r=6`.  Exact inversion
of the branch-resolved point cores gives

```text
G_{i,6}=sum_j (A_6)_{ij} G_j             for i=0,...,5. (6)
```

For a coefficient row `c`, define the branch-summed and literal-branch
profiles

```text
T(c)=sum_i c_i G_i,
L_6(c)=sum_i c_i G_{i,6}.                             (7)
```

Equations `(6)--(7)` give the literal intertwining identity

```text
L_6(c)=T(c A_6).                                      (8)
```

This is checked in the full `13 by 13` difference/relation profile, before
any eigenspace or quotient conclusion is used.

## 3. The rigidity criterion

The centered point generators are independent, so `(5)` and `(8)` imply,
for every `c in R`,

```text
L_6(c)-alpha T(c)
  =T(c(A_6-alpha I))
  =(beta-alpha)T(c Pi_M).                             (9)
```

On `R`, the kernel and image of `Pi_M` are `K` and `M`, respectively.
Therefore

```text
c in K     iff     L_6(c)=alpha T(c).                 (10)
```

The complementary statement is equally exact:

```text
{c in R:L_6(c)=beta T(c)}=M.                          (11)
```

Thus the four-dimensional relation space splits into the two literal
branch-response eigenspaces

```text
R=K_alpha direct_sum M_beta.                          (12)
```

## 4. The same defect realizes LambdaSharp

THM-3636 supplies an invertible matrix

```text
T_mid=
 [126498113787680818370196  646561219255993169342961]
 [384727693242765231857657  452865069702498262026030] (13)
```

and the factorization `LambdaSharp(c)=delta_mid(c)T_mid` on `R`.  Since
`(5)` gives

```text
delta_mid(c)=
  delta_mid(c(A_6-alpha I))/(beta-alpha),             (14)
```

one obtains the literal-defect formula

```text
LambdaSharp(c)=
  [delta_mid(c(A_6-alpha I))/(beta-alpha)] T_mid.     (15)
```

The companion checks `(15)` on an ordered basis of `R` before row-space
canonicalization.  Because `T_mid` and `beta-alpha` are invertible, the
literal branch defect, the middle-coordinate quotient, and `LambdaSharp`
have exactly the same kernel `K`; they differ only by invertible coordinates
on `R/K`.

## 5. Reproduction and boundary

Reproduce with

```bash
python3 04-computation/lrc_fixed_branch_rigidity_eigencriterion_thm3654.py
python3 -O 04-computation/lrc_fixed_branch_rigidity_eigencriterion_thm3654.py
```

The assertion-free companion raw-pins the audited THM-3647 theorem, script,
and output; reconstructs the centered six-point generators and the literal
branch-6 tensor; verifies `(2)--(15)`, both eigenspaces, all maintained
digests, source LF, and the canonical semantic digest.  Normal and optimized
streams must match the stored transcript.

The phrase *literal branch* means literal restriction in the pinned static
coefficient tensor.  It does not make `A_6` a temporal transition, identify a
lawful current operation, lift the statement to characteristic zero, exclude
a physical row, improve a lonely-runner bound, or prove LRC(14).  The next
typing problem is to identify whether any source/current label operation
implements `(8)` without choosing a noncanonical section.  **QED.**
