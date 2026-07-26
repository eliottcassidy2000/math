---
id: THM-2409
title: "Unfiltered septimal source completion and word-phase boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the
  source-deletion branch of THM-2407 with terminal word Q=1, translating
  the missing source danger through C_7 partitions the source-deleted
  packet exactly. For every nonzero target colour b, the seven
  inserted-source target coefficients lie in Q(zeta_13), sum to the
  nonzero deletion coefficient, and have anchored value zero at the
  original owner phase. Since Phi_7 is irreducible over Q(zeta_13),
  every six nonzero septimal source colours therefore survives jointly
  with every b. Source-colour-retained diagonal cancellation then
  produces a nonzero deep colour and an exact source-residue-refined
  fixed-frequency triangle with gcd(m,91)=1. The joint mixed energy is
  at least 9 rho^2/2798978. Without the anchored owner zero, flatness is
  sharp even for a nonflat target current and strictly positive even
  base factors. With a genuine delayed word,
  R=13^k is a unit mod 7, so the word's source factor must be shifted
  too; then the danger partition no longer reconstructs the fixed-word
  deletion current. More sharply, the THM-2334 exact-address gauge
  permits only the skew-diagonal source/word character omega=R ell;
  every relative word offset is gauge. No terminal-word repair,
  all-91-unit address, row exclusion, or LRC(14) conclusion is claimed.
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
script_sha256: d556d4873c6be2f704618160dfcb0f10f3d72022df8eb45b6e5f2f10b8ade67a
output_sha256: 3e3d4a98512bfb29e87cb1f8d30b260ea42fa62a00c5b7d3412dc976a076d2ff
hash_basis: working-tree bytes (LF)
---

# THM-2409 -- the unfiltered deletion current is septimally all or flat

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2407 leaves one sharply typed alternative. Either the genuine
positive source owner carries all target colours, or the packet with the
stationary source factor deleted does. In the second branch, translating
the missing danger through seven phases restores a genuine nonzero
septimal source residue.

With no terminal word this operation is exact. The seven translated
dangers partition one, so their joint currents sum to the
source-deletion current. Cyclotomic disjointness then gives a complete
alternative:

```text
every b!=0 and every kappa!=0
  -> a joint source/target coefficient survives
  -> some nonzero deep colour and exact 91-unit deep leg.     (1)
```

The anchored zero of the original owner phase is load-bearing. Without
it, a flat source profile is real. A delayed word does not preserve the
positive argument for free:
the word is target-neutral modulo thirteen but active modulo seven, so
a lawful full source twist must move the word factor too and loses the
partition identity.

## 1. The exact seven-shift partition

Use the source-deletion branch of THM-2407 with the unfiltered terminal
word

```text
Q=1.
```

Put

```text
rho=delta=mu(S)>0,
```

where `S` is the clean set used by THM-2403. For `Q=1`, this is exactly
the overlap parameter in THM-2407.

Let `U_(s,t)` be its source-deleted present packet and retain the same
deep probe

```text
Delta_r(x)=d(c_3x-r/13).                             (2)
```

Write the target marginal of `U` as `u_s`. In the deletion branch,
THM-2407 proves, for every `b!=0`,

```text
uhat(b)=(1/13)sum_s u_s zeta_13^(bs)!=0.            (3)
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
 =integral_T U_(s,t)(x)d_(j,ell)(x)Delta_r(x) dx,   (6)

z_ell(b)
 =1/13 sum_s sum_r H_ell(r,s,0)zeta_13^(bs).        (7)
```

Equations (5)--(7) give

```text
sum_ell z_ell(b)=uhat(b)!=0.                        (8)
```

At `ell=0`, the inserted factor is the original positive source-owner
danger. The defining branch of THM-2407 says that this owner's target
marginal is flat. Hence

```text
z_0(b)=0                    for every b!=0.         (9)
```

All interval breakpoints in (6) are rational. Therefore every table
entry is rational and, for fixed `b`,

```text
z_ell(b) in K:=Q(zeta_13).                         (10)
```

No Poisson idempotence is used: (5)--(10) are indicator-boundary
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
[K(zeta_7):K]=6,                                   (11)
```

and `Phi_7=1+X+...+X^6` is irreducible over `K`.

For each fixed `b!=0`, define the normalized joint source/target
transform

```text
Z(kappa,b)
 =(1/7)sum_(ell in F_7)z_ell(b)zeta_7^(kappa ell). (12)
```

If `Z(kappa,b)=0` for one `kappa!=0`, then the degree-at-most-six
polynomial

```text
P_b(X)=sum_(ell=0)^6 z_ell(b)X^ell
```

is divisible by `Phi_7` over `K`. Hence all seven coefficients are
equal. Equation (9) would make them all zero, contradicting (8).
Therefore

```text
Z(kappa,b)!=0
  for every kappa!=0 and every b!=0.               (13)
```

The trivial source colour also survives:

```text
Z(0,b)=uhat(b)/7!=0.                               (14)
```

This proves (1). The conclusion holds separately for every prescribed
target colour, while one global anchored-zero argument handles all six
source colours.

There is also an exact energy floor. Parseval gives

```text
sum_kappa |Z(kappa,b)|^2
 =(1/7)sum_ell |z_ell(b)|^2.
```

Under `z_0=0` and `sum_ell z_ell=uhat(b)`, Cauchy is sharp at
`z_1=...=z_6=uhat(b)/6`. Removing the trivial source colour yields

```text
sum_(kappa!=0)|Z(kappa,b)|^2
 >=|uhat(b)|^2/294.                                (15)
```

This equality boundary is realized by a nonnegative owner table.  Give
the original source phase the flat target row

```text
o_(0,s)=1/14,
```

and, for each `ell!=0`, put

```text
o_(ell,0)=0,             o_(ell,s)=1/84 for s!=0.
```

Then `z_0(b)=0` for every `b!=0`, while the other six `z_ell(b)` are
equal and their sum is `uhat(b)`.  Tensoring with `r=t+1` also preserves
the empty deep diagonal.  Thus the constant `1/294` cannot be improved
from the anchored-zero and sum constraints alone.

Let

```text
D_mix=sum_(b!=0,kappa!=0)|Z(kappa,b)|^2.
```

THM-2407's deletion-branch target-energy floor now gives

```text
D_mix>=9rho^2/2798978,                             (16)

max_(b!=0,kappa!=0)|Z(kappa,b)|>=rho/4732.         (17)
```

## 3. Exact relation-address meaning

Retain the source character in the full table:

```text
B(kappa,alpha,b,tau)
 =1/(7*13^3) sum_(ell,r,s,t)
    H_ell(r,s,t)
    zeta_7^(kappa ell)
    zeta_13^(alpha r+b s+tau t),                   (18)

J(kappa,alpha,b)=sum_tau B(kappa,alpha,b,tau).     (19)
```

Equations (7), (12), and (19) give

```text
J(kappa,0,b)=Z(kappa,b)/13!=0.                     (20)
```

Every `H_ell` retains the moving deepest-safe factor, so
`H_ell(t,s,t)=0`. At `t=r=0`, source-Fourier-retained diagonal
cancellation gives

```text
sum_alpha J(kappa,alpha,b)=0.                      (21)
```

Thus, for every `kappa!=0,b!=0`, some `alpha!=0` and `tau` satisfy

```text
B(kappa,alpha,b,tau)!=0.                           (22)
```

The variable `ell` is the lawful coordinate translation of the source
present and bare factors in THM-2334. Hence (22) selects source relation
residue `kappa` modulo seven. The absolutely convergent expansion gives
a **source-residue-refined fixed-`(N,m)` coefficient** with:

```text
first-target colour b!=0 mod 13,

source residue kappa!=0 mod 7,

deep multiplier gcd(m,91)=1.                       (23)
```

Here `m=alpha mod 13`, and the centered deep coefficient kills every
nonzero `7|m`. The corresponding unrefined fixed-`(N,m)` coefficient may
cancel after source residues are summed; (22) retains the residue.
Moreover, `alpha` and `tau` are selected only after the source transform
has selected `(kappa,b)`.  The resulting triangle need not be the
preselected THM-2407 triangle, and it cannot be intersected with a
THM-2410 survivor without an additional common selector.

Quantitatively, Cauchy first distributes each `J(0)` cancellation over
the twelve nonzero `alpha`, and then each `J` over thirteen `tau`.
Therefore the eligible full-table energy satisfies

```text
sum_(kappa,b,alpha!=0,tau)|B|^2
 >=D_mix/(169*12*13).                              (24)
```

This does not say that the source coordinate is nonzero modulo
thirteen, nor that every other relation coordinate is a unit modulo
seven. It is a septimal source completion of one unfiltered aggregate,
not an all-`91`-unit relation address.

## 4. The anchored zero is load-bearing

The positive conclusion uses both parts of the deletion branch:
`z_0(b)=0` and `sum_ell z_ell(b)!=0`. If the anchored owner zero is
discarded, flatness remains possible despite nonnegativity, target
charge, and even diagonal-zero structure.

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
C_ell(s)=h_s/7              for every ell.          (25)
```

Then the source profile is perfectly flat for every target cell, while
the target current is nonflat. Here `C_0` has the same nonzero target
charge, so the anchored-zero hypothesis fails exactly.

One may also retain the deep diagonal. Put

```text
M_(ell,s)=(6/7)1_(s!=0),

H_ell(r,s,t)=M_(ell,s)1_(r=t+1).
```

This nonnegative rational table is source-flat, diagonal-zero, and
fires every `b!=0` together with every nonzero deep `alpha` at
`tau=-alpha`. Thus diagonal cancellation alone does not replace (9).

### 4.2 Positive even one-circle control

The same phenomenon occurs before any finite abstraction. Choose

```text
epsilon=delta=1/2,

A_s(x)=1+epsilon cos(2pi x-2pi s/13),

G(x)=1+delta cos(2pi x),                            (26)
```

and use the source danger `d(Nx-ell/7)` with any integer `N>=3`.
The base functions in (26) are strictly positive, and their unshifted
forms are real and even.

The Fourier support of `A_sG` lies in

```text
{-2,-1,0,1,2},
```

whereas the source danger has support in `N Z`. Only frequency zero
meets. Exact integration gives

```text
integral_T A_s(x)G(x)d(Nx-ell/7)dx
 =1/7 (1+1/8 cos(2pi s/13)),                       (27)
```

independently of `ell`. Its target modes `b=+1,-1` have amplitude

```text
1/112.                                             (28)
```

Thus even a strictly positive one-circle target current can be
septimally flat when its original source phase is not pinned to zero.

## 5. Why a real terminal word breaks the partition

Return to a genuine THM-2305 word transported by

```text
R=13^k.
```

Modulo thirteen, `R` is zero and the word is target-neutral. Modulo
seven,

```text
R=(-1)^k mod 7,                                    (29)
```

so its source factor is active.

A full source-coordinate twist of the THM-2334 current must translate:

```text
the left present source factor,
the bare right source factor,
and the transported-word source factor             (30)
```

with the latter moving by `R ell/7`. If the word is held fixed, (5)
still partitions the deletion packet, but the Fourier character omits
the word charge `R beta_j` and is not the full relation-address
coordinate. If the word is shifted lawfully, its factor becomes
`Q_ell` and only

```text
sum_ell d_(j,ell)Q_ell                             (31)
```

appears. There is no identity equating (31) with the fixed-word
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
sum_ell D_ell W_ell=0.                             (32)
```

Thus the active word can destroy the partition completely.

### 5.1 Skew-diagonal descent obstruction

The tempting repair is to introduce an independent relative word offset.
Algebraically that exposes real information, but it is not a lawful
exact-address coordinate.

Use THM-2334's decomposition variables `(u,beta,v)` at clock `R`. A putative
septimal character which shifts present/bare occurrences by `ell` and the
transported-word occurrence by `omega` is

```text
chi_(ell,omega)(u,beta,v)
 =e_7(ell.u+omega.beta-ell.v).                    (32a)
```

Under the exact-address gauge action

```text
(u,beta,v) -> (u+p,beta+q,v+p+Rq),
```

the multiplier of (32a) is

```text
e_7((omega-Rell).q).                               (32b)
```

For every `q`, the gauge contains `(p,q)=(-Rq,q)`. Therefore

```text
chi_(ell,omega) descends to the exact-address quotient
iff omega=Rell mod 7.                              (32c)
```

The lawful source/word characters are precisely this skew diagonal. A
relative offset `omega_j-Rell_j` for one source occurrence is gauge-invariant
only when it is zero.

The off-diagonal algebra is nevertheless exact. Put

```text
epsilon=R mod 7 in {+1,-1},

P_(ell,r)=D_ell W_(epsilon ell+r).
```

Its two-dimensional normalized transform factors as

```text
Phat(k,h)=Dhat(k-epsilon h)What(h).                 (32d)
```

For the singleton hostile `D_ell(x)=1_(x=ell)` and
`W_a=1-D_a`, the lawful diagonal at `epsilon=1` vanishes while every
off-diagonal averaged offset is `1/7`. Thus polarization sees exactly the
information that (32c) declares gauge-dependent.

Changing the clock changes only `epsilon=(-1)^k`; it does not create a
transverse offset. Translating a THM-2305 word changes the literal terminal
stratum rather than producing another character of the same current.
Candidate THM-2414 points to the correctly typed replacement: retain the
parent-dependent affine digit `floor(Ry) mod 7` and its torsor map on one
same-parent cell. It supplies one physical offset per parent, not seven
independent copies.

Any lawful repair must therefore retain the parent cell, source owner,
terminal word, clock, both septimal phase torsors, their affine transporter,
the common root sheet/right-endpoint gauge, and target/deep colours in one
signed joint coefficient. An off-diagonal Fourier energy without that
sidecar cannot repair (31).

## 6. Consequence and exact boundary

On the THM-2407 deletion branch:

```text
Q=1:
  every b!=0 and kappa!=0 has a joint source/target mode,
  and diagonal cancellation lands a nonzero deep colour;

Q a real delayed word:
  the complete source twist is a coupled present/word convolution,
  not the unfiltered partition.                    (33)
```

The next useful sidecar must control the mod-seven coupled word
convolution in one common phase gauge. A lone source danger probe is
decisive for `Q=1` but is not a full-current source twist once the word
is active.

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
- verifies the sharp anchored-source energy `1/294`, the mixed floor
  `9/2798978`, the `4732` max-mode denominator, and a positive equality
  table;
- checks the rational finite-table, diagonal-zero, and smooth
  low-frequency flat hostiles; and
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

## 8. Independent audit

Two independent hostile audits rederived:

- the lawful `Q=1` rerun of THM-2407 and the exact identities
  `z_0(b)=0`, `sum_ell z_ell(b)=uhat(b)`;
- irreducibility of `Phi_7` over `Q(zeta_13)` and simultaneous survival
  of all six nontrivial source colours for every `b!=0`;
- the sharp source-energy constant `1/294`, its positive equality
  table, the mixed floor `9rho^2/2798978`, and the max-mode denominator
  `4732`;
- the normalizations `J(kappa,0,b)=Z(kappa,b)/13` and
  `D_mix/(169*12*13)` for the eligible full-table energy;
- the source-residue-refined absolute `(N,m)` extraction and
  `gcd(m,91)=1`, including the fact that it need not retain the
  preselected THM-2407 triangle; and
- the exact failure of the seven-shift partition for a lawfully
  co-shifted active word.

Both audits replayed normal and optimized execution against the stored
transcript and checked the LF hashes.
