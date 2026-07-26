---
id: THM-2409
title: "Unfiltered septimal source completion and word-phase boundary"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING. In the
  source-deletion branch of THM-2407 with terminal word Q=1, translating
  the missing source danger through C_7 partitions the source-deleted
  packet exactly. For any fixed nonzero C_13 target/deep coefficient,
  the seven inserted-source coefficients lie in Q(zeta_13), sum to the
  original coefficient, and are either flat or have every six nonzero
  septimal Fourier colours. This follows because Phi_7 is irreducible
  over Q(zeta_13). In the nonflat branch each joint colour supplies a
  source residue nonzero mod 7; in the flat branch only source residue
  zero survives. Flatness is sharp even for a nonflat target current and
  strictly positive even base factors. With a genuine delayed word,
  R=13^k is a unit mod 7, so the word's source factor must be shifted
  too; then the danger partition no longer reconstructs the fixed-word
  deletion current. No terminal-word repair, all-91-unit address, row
  exclusion, or LRC(14) conclusion is claimed.
source: codex-2026-07-26-septimal-source-completion
depends_on:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
  - THM-2407-owner-or-source-deletion-target-current-dichotomy
related:
  - THM-2400-clean-parent-root-gauge-quotient-and-target-slope-boundary
  - THM-2408-endpoint-prony-resultant-clock-separation-and-shared-node-boundary
  - THM-2410-full-coordinate-projector-local-gram-and-integrated-phase-boundary
script: 04-computation/lrc14_unfiltered_septimal_source_completion_thm2409.py
output: 05-knowledge/results/lrc14_unfiltered_septimal_source_completion_thm2409.out
script_sha256: 6af42ee995bf31d5cd3d71ea85eb0ec7fa252eccd4081840485c406844f6c380
output_sha256: 93f7c06f0bc30dbbf3be25cfda86f90a2639c711f6a3f15aae2b23d34b50215f
hash_basis: working-tree bytes (LF)
---

# THM-2409 -- the unfiltered deletion current is septimally all or flat

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-2407 leaves one sharply typed alternative. Either the genuine
positive source owner carries all target colours, or the packet with the
stationary source factor deleted does. In the second branch one may try
to restore the missing source coordinate by translating its danger
factor through seven phases.

With no terminal word this operation is exact. The seven translated
dangers partition one, so their joint currents sum to the
source-deletion current. Cyclotomic disjointness then gives a complete
alternative:

```text
source phase nonflat
  -> all six nonzero septimal source colours survive;

source phase flat
  -> only septimal source residue zero survives.             (1)
```

The flat branch is real. A delayed word does not remove it for free:
the word is target-neutral modulo thirteen but active modulo seven, so
a lawful full source twist must move the word factor too and loses the
partition identity.

## 1. The exact seven-shift partition

Use the source-deletion branch of THM-2407 with the unfiltered terminal
word

```text
Q=1.
```

Let `U_(s,t)` be its source-deleted present packet and retain the same
deep probe

```text
Delta_r(x)=d(c_3x-r/13).                             (2)
```

THM-2407 supplies, for each prescribed first-target colour `b!=0`,
some `alpha!=0` and `tau` such that

```text
B_U(alpha,b,tau)!=0.                                (3)
```

For the omitted source label `j`, put

```text
d_(j,ell)(x)=d(c_jx-ell/7),       ell in F_7.       (4)
```

The seven strict-open translates of the interval of length `1/7`
partition the circle away from their null endpoints:

```text
sum_(ell in F_7)d_(j,ell)(x)=1                     (5)
```

almost everywhere.

At the indicator boundary, insert the same `d_(j,ell)` into the left
present and bare right copies of the missing source factor. Their
product is again `d_(j,ell)`. Define

```text
H_ell(r,s,t)
 =integral_T U_(s,t)(x)d_(j,ell)(x)Delta_r(x) dx,  (6)

C_ell
 =1/13^3 sum_(r,s,t)
    H_ell(r,s,t)zeta_13^(alpha r+b s+tau t).        (7)
```

Equations (5)--(7) give the exact coefficient identity

```text
sum_ell C_ell=B_U(alpha,b,tau)!=0.                  (8)
```

All interval breakpoints in (6) are rational. Therefore every table
entry is rational and

```text
C_ell in K:=Q(zeta_13).                             (9)
```

No Poisson idempotence is used: (5)--(8) are indicator-boundary
identities, after which the existing Abel boundary theorem applies.

## 2. Why `Phi_7` remains irreducible

The cyclotomic fields

```text
Q(zeta_7) and Q(zeta_13)
```

have coprime odd conductors. Their intersection is `Q`. One elementary
proof is by ramification: every finite prime ramified in the
intersection would have to ramify in both fields, hence divide both
`7` and `13`; a nontrivial number field cannot have discriminant one.
Consequently

```text
[K(zeta_7):K]=6,                                   (10)
```

and `Phi_7=1+X+...+X^6` is irreducible over `K`.

Define the normalized septimal transform

```text
Ctilde(e)=(1/7)sum_(ell in F_7)C_ell zeta_7^(e ell).
                                                               (11)
```

If `Ctilde(e)=0` for one `e!=0`, then the degree-at-most-six
polynomial

```text
P(X)=sum_(ell=0)^6 C_ell X^ell
```

is divisible by `Phi_7` over `K`. Hence all seven coefficients are
equal. The converse is immediate. Thus:

```text
(C_ell) nonflat
  iff Ctilde(e)!=0 for every e!=0;                 (12)

(C_ell) flat
  iff Ctilde(e)=0 for every e!=0.                  (13)
```

Equation (8) also gives

```text
Ctilde(0)=B_U(alpha,b,tau)/7!=0.                   (14)
```

This proves the alternative (1) globally, not colour by colour.

## 3. Exact relation-address meaning

The variable `ell` is the lawful coordinate translation of the source
present and bare factors in THM-2334. Therefore a nonzero joint
coefficient

```text
Ctilde(e),            e!=0 mod 7,                  (15)
```

selects source relation residue `e` modulo seven. Expanding the fixed
target/deep coefficient as in THM-2365 supplies an exact frequency
triangle with:

```text
first-target colour b!=0 mod 13,

source residue e!=0 mod 7,

deep multiplier gcd(m,91)=1.                       (16)
```

This does not say that the source coordinate is nonzero modulo
thirteen, nor that every other relation coordinate is a unit modulo
seven. It is a septimal source completion of one unfiltered aggregate,
not an all-`91`-unit relation address.

In the flat branch every inserted-source coefficient equals

```text
C_ell=B_U(alpha,b,tau)/7,                           (17)
```

and the whole nonzero source spectrum vanishes. The unfiltered target
triangle remains alive only on source residue zero.

## 4. Two sharp flat-source hostiles

Flatness cannot be removed from (1) by nonnegativity.

### 4.1 Exact rational finite table

Let the target marginal be

```text
h=(1,2,2,...,2) in Q^13.
```

It has all twelve nonzero target colours and exact charged energy

```text
12/169.
```

Put

```text
C_ell(s)=h_s/7              for every ell.          (18)
```

Then the source profile is perfectly flat for every target cell, while
the target current is nonflat. This is an exact nonnegative rational
finite-projector hostile.

### 4.2 Positive even one-circle control

The same phenomenon occurs before any finite abstraction. Choose

```text
epsilon=delta=1/2,

A_s(x)=1+epsilon cos(2pi x-2pi s/13),

G(x)=1+delta cos(2pi x),                             (19)
```

and use the source danger `d(Nx-ell/7)` with any integer `N>=3`.
The base functions in (19) are strictly positive, and their unshifted
forms are real and even.

The Fourier support of `A_sG` lies in

```text
{-2,-1,0,1,2},
```

whereas the source danger has support in `N Z`. Only frequency zero
meets. Exact integration gives

```text
integral_T A_s(x)G(x)d(Nx-ell/7)dx
 =1/7 (1+1/8 cos(2pi s/13)),                        (20)
```

independently of `ell`. Its target modes `b=+1,-1` have amplitude

```text
1/112.                                              (21)
```

Thus even a strictly positive one-circle target current can be
septimally flat.

## 5. Why a real terminal word breaks the partition

Return to a genuine THM-2305 word transported by

```text
R=13^k.
```

Modulo thirteen, `R` is zero and the word is target-neutral. Modulo
seven,

```text
R=(-1)^k mod 7,                                     (22)
```

so its source factor is active.

A full source-coordinate twist of the THM-2334 current must translate:

```text
the left present source factor,
the bare right source factor,
and the transported-word source factor              (23)
```

with the latter moving by `R ell/7`. If the word is held fixed, (5)
still partitions the deletion packet, but the Fourier character omits
the word charge `R beta_j` and is not the full relation-address
coordinate. If the word is shifted lawfully, its factor becomes
`Q_ell` and only

```text
sum_ell d_(j,ell)Q_ell                              (24)
```

appears. There is no identity equating (24) with the fixed-word
source-deletion packet.

The failure is already exact on `C_7`. Let

```text
D_ell(r)=1_(r=ell).
```

For a fixed word `W_0=1-D_0`,

```text
sum_ell D_ell W_0=W_0.
```

For the lawfully co-shifted words `W_ell=1-D_ell`,

```text
sum_ell D_ell W_ell=0.                              (25)
```

Thus the active word can destroy the partition completely.

## 6. Consequence and exact boundary

On the THM-2407 deletion branch:

```text
Q=1:
  every nonzero target/deep triangle is either
  source-flat or carries all six nonzero C_7 source colours;

Q a real delayed word:
  the complete source twist is a coupled present/word convolution,
  not the unfiltered partition.                    (26)
```

The next useful sidecar must therefore prove nonflatness from physical
endpoint geometry, or control the mod-seven coupled word convolution in
one common phase gauge. A lone source danger probe is insufficient.

No terminal word is restored, no all-`91`-unit address is produced, no
scalar row is excluded, the ledger remains `165`, and LRC(14) remains
open.

## 7. Exact companion

The dependency-free companion uses integer and `Fraction` arithmetic
only. It:

- checks the seven-shift danger partition;
- exhausts all `2^7=128` Boolean source profiles;
- verifies the all-or-flat reduction over a symbolic
  `Q(zeta_13)` coefficient basis;
- checks the rational finite-table and smooth low-frequency flat
  hostiles; and
- verifies the fixed-word partition and its complete failure under the
  coupled active-word shift.

Run:

```bash
python3 04-computation/lrc14_unfiltered_septimal_source_completion_thm2409.py
python3 -O 04-computation/lrc14_unfiltered_septimal_source_completion_thm2409.py
```

Both commands must reproduce

```text
05-knowledge/results/lrc14_unfiltered_septimal_source_completion_thm2409.out
```

with the LF hashes in the frontmatter.
