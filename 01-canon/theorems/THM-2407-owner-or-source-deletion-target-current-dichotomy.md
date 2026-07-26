---
id: THM-2407
title: "Owner-or-source-deletion target-current dichotomy"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.
  In the repaired THM-2403 owner pivot, the shallow source factor is
  stationary under both target shifts. The packet U with that one
  factor omitted decomposes target-equivariantly as O+H_+, where O is
  the genuine positive source-owner packet and H_+ is the fully-all-safe
  complement. Scalar cover gives U_00=O_00 pointwise, while all three
  tables retain nonnegativity and the deepest diagonal zero. The
  rational C_13 source-owner marginal is globally either nonflat or
  flat. In the first branch O fires every nonzero target colour; in the
  second its charged spectrum vanishes identically and U inherits every
  THM-2403 colour from H_+. The same diagonal argument then gives an
  exact 91-unit deep triangle in the selected bank for each first-target
  colour. In the second branch U is a same-base source-factor-deletion
  current, not a canonical nine-factor owner current. Common-base
  equality does not identify twisted coefficients, and THM-2370's clone
  hostile forbids recovering that missing orientation from squared
  deletion data alone. No row exclusion or LRC(14) conclusion is
  claimed.
source: codex-2026-07-26-owner-or-source-deletion-dichotomy
depends_on:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
  - THM-2403-clean-toothpick-unequal-slope-target-axis-imbalance
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2401-common-filter-endpoint-or-first-death-certificate
script: 04-computation/lrc14_owner_or_source_deletion_dichotomy_thm2407.py
output: 05-knowledge/results/lrc14_owner_or_source_deletion_dichotomy_thm2407.out
script_sha256: 6241cf55b45e2b1785047bda91623340e87550174d6a4c2d33a3e4f6561703a8
output_sha256: 4c78863c470c450077624f4b202c9f40d01451a0b96d8b6d88b1c4655ff38ca4
hash_basis: working-tree bytes (LF)
---

# THM-2407 -- every target colour lands either in the owner or in its source deletion

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-2403 supplies a lawful fully-all-safe present current with every
nonzero first-target colour. That current is not the selected positive
shallow-owner current. The missing source factor, however, is stationary
under the two target dipoles. Restoring it therefore gives an exact
target-equivariant Boolean partition, rather than the generally
noncovariant deletion routing obstructed in THM-2370.

Rational prime-cyclic rigidity turns this partition into one global
alternative:

```text
genuine positive source owner has all twelve target colours

or

the packet with only that source factor deleted has all twelve.
                                                               (1)
```

The second branch shares its untwisted physical slice with the owner,
but it does not share the owner's twisted coefficients. This distinction
is the theorem's sharp stopping boundary.

## 1. The three lawful packets

Use the repaired THM-2403 owner pivot. Thus `j` is the positive shallow
source label from THM-2349, `a` is the other low blocker,

```text
c=c_3,
```

and the omitted pivot unit is `q_*`. Choose the two target dipoles as in
THM-2403:

```text
eta=e_a-e_(q_i),          ell=e_c-e_(k_c),           (2)
```

where `q_i,k_c,q_*` are distinct ordinary unit labels. In particular,
the source blocker is target-neutral:

```text
eta_j=ell_j=0.                                      (3)
```

Write

```text
d_j(x)=d(c_j x),             g_j(x)=1-d_j(x).       (4)
```

Let `U_(s,t)` be the product of all shifted present factors in the
THM-2403 packet except the stationary source-blocker factor. Define

```text
O_(s,t)=U_(s,t)d_j,          P^+_(s,t)=U_(s,t)g_j.  (5)
```

Thus:

- `O` is the genuine exclusive source-owner present packet;
- `P^+` is the fully-all-safe packet of THM-2403; and
- `U` has the source factor deleted and retains every other present
  factor.

Strict-open endpoints are null, so the Boolean partition is exact
almost everywhere:

```text
U_(s,t)=O_(s,t)+P^+_(s,t)       for every s,t.      (6)
```

Because of (3), (6) commutes with the actual target action. No clean-set
projector or frozen moving factor has been inserted.

Choose a positive rational THM-2305 terminal word belonging to the
source `j`, transported by a sufficiently large clock

```text
Q(Rx),                    R=13^k.                   (7)
```

THM-2349 makes the untwisted owner overlap positive. Increasing the
clock if necessary, the THM-2397 two-BV estimate also gives the positive
clean overlap `rho_R` used by THM-2403. Since `13|R`, the word is
target-neutral.

For `X` equal to `O,U`, or `P^+`, define

```text
H_X(r,s,t)
 =integral_T X_(s,t)(x)Q(Rx)d(c x-r/13) dx.         (8)
```

All three tables are nonnegative rational tables. Equation (6) gives

```text
H_U=H_O+H_+                                         (9)
```

entry by entry.

Every packet retains the moving `c`-safe factor. Hence

```text
H_X(t,s,t)=0                for every s,t,           (10)
```

for all three choices of `X`.

At `(s,t)=(0,0)`, scalar cover says that the fully-all-safe present
packet vanishes pointwise:

```text
P^+_(0,0)=0.
```

Consequently

```text
U_(0,0)=O_(0,0),

H_U(r,0,0)=H_O(r,0,0)       for every r.            (11)
```

This is an exact common bare slice. It is not an equality of the
twisted tables.

## 2. The global all-or-flat alternative

Put

```text
o_s=sum_r H_O(r,s,0),

h_s=sum_r H_+(r,s,0),

u_s=sum_r H_U(r,s,0).                               (12)
```

Then

```text
u_s=o_s+h_s,                u_0=o_0,                (13)
```

because `h_0=0`. Normalize the `C_13` transform by

```text
xhat(b)=(1/13)sum_(s in F_13)x_s zeta^(bs).         (14)
```

THM-2403 proves

```text
hhat(b)!=0                  for every b!=0.         (15)
```

The owner marginal `(o_s)` is rational. By the prime-cyclic
all-or-flat theorem of THM-2398, exactly one of the following occurs.

### Owner branch

If `(o_s)` is nonflat, then

```text
ohat(b)!=0                  for every b!=0.         (16)
```

Thus one fixed bank, the genuine positive owner `O`, carries all twelve
first-target colours.

If

```text
o_s=a_s/D,        a_s in Z_(>=0),

S=sum_s a_s,      g=gcd_s a_s,                     (17)
```

THM-2398 also gives the denominator-sensitive floor

```text
|ohat(b)|>=g^6/(13 D S^5)     for every b!=0.       (18)
```

There is no denominator-free positive lower bound.

### Source-deletion branch

If `(o_s)` is flat, then

```text
ohat(b)=0,

uhat(b)=hhat(b)!=0            for every b!=0.       (19)
```

Now the one fixed source-deleted bank `U` carries all twelve colours.
It inherits the exact THM-2403 aggregate bounds

```text
sum_(b!=0)|uhat(b)|^2>=27rho_R^2/28561,

max_(b!=0)|uhat(b)|>=3rho_R/338.                    (20)
```

The branch decision is global. There is no colour-by-colour mixture.

## 3. Every selected colour gives a deep triangle

For the bank `X` selected by Section 2, define

```text
B_X(alpha,b,tau)
 =1/13^3 sum_(r,s,t)
    H_X(r,s,t)zeta^(alpha r+b s+tau t),             (21)

J_X(alpha,b)
 =sum_tau B_X(alpha,b,tau).                         (22)
```

At `t=0`, (12)--(14) give

```text
J_X(0,b)=xhat(b)/13,                                (23)
```

where `x=o` in the owner branch and `x=u` in the deletion branch.
Equation (10) gives the exact target-line cancellation

```text
sum_alpha J_X(alpha,b)=0.                           (24)
```

Fix `b!=0`. The selected transform in (23) is nonzero. Therefore
some `alpha!=0` has `J_X(alpha,b)!=0`, and then some `tau` has

```text
B_X(alpha,b,tau)!=0.                                (25)
```

The THM-2365 target typing and absolutely convergent exact-frequency
expansion apply unchanged:

```text
q=(b,alpha+tau),

B_X(alpha,b,tau)
 =sum_(m=alpha mod 13) sum_N A^X_(N,m)(q).          (26)
```

Thus (25) supplies exact `N,m` with

```text
A^X_(N,m)(q)!=0,

13 does not divide m.                               (27)
```

The centered deep-danger coefficient vanishes when `7|m`, so every
live term also has `7` not dividing `m`. Hence

```text
gcd(m,91)=1.                                        (28)
```

This holds separately for each prescribed `b in F_13^*`.

In the owner branch, (25)--(28) belong to the genuine positive
nine-present-factor source-owner current. In the deletion branch, the
source factor is the constant `1`; its Fourier harmonic is forced to
zero. The latter is a lawful exact source-factor-deletion current with
the same untwisted slice (11), but it is not an all-coordinate owner
address.

## 4. Sharp branch controls

The dichotomy cannot be strengthened by always choosing `O` or always
choosing `U`. At the marginal level, indexed by `s=0,...,12`, take

```text
h=(0,1,1,...,1).                                    (29)
```

For the deletion branch let

```text
o=(1,1,...,1),

u=o+h=(1,2,2,...,2).                                (30)
```

Then `u_0=o_0`, the owner is flat, and

```text
energy(o)=0,

energy(h)=energy(u)=12/169.                         (31)
```

For the owner branch let

```text
o=(2,1,1,...,1),

u=o+h=(2,2,...,2).                                  (32)
```

Now the source-deleted bank is flat, while

```text
energy(o)=energy(h)=12/169,

energy(u)=0.                                        (33)
```

Both controls are nonnegative rational tables and satisfy the common
base identity. They show that cancellation between `O` and `H_+` can
erase either candidate spectrum.

There are also nonnegative `13^3` diagonal-zero realizations of both
controls: put all mass at `(r,t)=(1,0)`. Then (10), (23), and (24) hold
exactly, and every nonzero marginal colour extracts a nonzero
`alpha!=0`. The companion checks these realizations.

## 5. What common-base equality does not restore

Equation (11) is stronger than separate positive masses: it identifies
the entire untwisted physical slice. It still does not imply

```text
B_U(alpha,b,tau)=B_O(alpha,b,tau)                   (34)
```

away from `s=t=0`. The source-deletion branch therefore cannot be
renamed a canonical owner current.

THM-2370 gives the sharp physical no-reference boundary. Its clone cube
has:

```text
one common bare endpoint,
equal and collinear signed singleton deletions,
zero mixed Boolean ANOVA in every degree >=2,
and the complete squared deletion spectrum,         (35)
```

yet Boolean complementation reverses which endpoint is terminal without
changing those squared data. Thus common-base equality plus unsigned
deletion tomography cannot restore the missing source orientation or
terminal phase. A signed charged reference, or a same-clock
target-covariant owner reconstruction, is still necessary.

Two tempting follow-ons have exact stopping rules.

1. Shifting the source danger through its `C_7` phase is a genuine full
   coordinate twist only when every occurrence of that factor is
   shifted. For `Q=1` this gives a useful unfiltered positive control.
   For a real THM-2305 word, however, `R mod 7` is `+1` or `-1`; the
   transported source factor is septimally active. Shifting only the
   present danger is not a canonical relation-address twist.

2. Endpoint-Prony recurrences can force a clock survivor when the owner
   and deletion exponential sums have coprime node polynomials. A shared
   node permits exact cancellation `O_h=-H_h` at every later clock,
   despite nonnegative underlying packets. Resultants therefore require
   endpoint-node separation in one common gauge before they improve
   (1).

## 6. Consequence and remaining bridge

The theorem removes per-colour owner routing from the ledger:

```text
one globally selected lawful bank
  -> all twelve first-target colours
  -> an exact 91-unit deep triangle for each colour. (36)
```

What remains is a single typed alternative:

```text
positive canonical owner current

or

same-base source-factor-deletion current.            (37)
```

Closing (37) requires a signed source reference, a lawful septimal
source-phase action including the transported word, or an exact
same-clock reconstruction of the missing source harmonic. None is
proved here. No scalar row is excluded, the ledger remains `165`, and
LRC(14) remains open.

## 7. Exact companion

The dependency-free companion uses integer and `Fraction` arithmetic
only. It:

- exhausts all `2^13=8192` Boolean `C_13` marginals and verifies that
  precisely the two flat words lose a nonzero character;
- checks both sharp nonnegative branch controls and their exact charged
  energies;
- realizes each selected marginal as a sparse nonnegative
  diagonal-zero `13^3` table; and
- verifies `J(0,b)` and a nonzero `alpha` for all twelve colours.

Run:

```bash
python3 04-computation/lrc14_owner_or_source_deletion_dichotomy_thm2407.py
python3 -O 04-computation/lrc14_owner_or_source_deletion_dichotomy_thm2407.py
```

Both commands must reproduce

```text
05-knowledge/results/lrc14_owner_or_source_deletion_dichotomy_thm2407.out
```

with the LF hashes in the frontmatter.
