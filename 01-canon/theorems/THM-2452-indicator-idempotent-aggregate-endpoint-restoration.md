---
id: THM-2452
title: "Indicator-idempotent aggregate endpoint restoration"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. At the
  indicator boundary, a complete Boolean mask already
  present on the moving endpoint can be copied to the partial bare
  endpoint without changing the nonnegative THM-2365 table:
  orthogonal masks vanish and a matched mask squares to itself. The
  five omitted guard/unit bits and two blocker bits give 128 complete
  masks. Some has bare drift at least D_0/16384, eligible target/deep
  energy at least D_0/212992, and a coefficient of squared modulus at
  least D_0/429391872. Equivalently, a preselected nonzero finite
  colour in THM-2445 cell i survives in one complete extension with
  loss at most 2^(5-i), and with no loss for i=0. Full-X
  recombination annihilates every THM-2448 transition, after which
  THM-2365 reselects a fresh exact X and a 91-unit m. The old X,m and
  semantic THM-2305 owner/repair alignment are not preserved; no
  scalar row is removed.
source: codex-2026-07-26-aggregate-endpoint-idempotence
depends_on:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2442-delayed-word-septimal-source-completion
  - THM-2445-twenty-four-cell-graft-owner-conditioning
  - THM-2448-right-endpoint-cospan-transition-atlas
related:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
  - THM-2379-anchored-guard-unit-deletion-factor-repair-current
  - THM-2401-common-filter-endpoint-or-first-death-certificate
  - THM-2408-endpoint-prony-resultant-clock-separation-and-shared-node-boundary
  - THM-2419-affine-relation-shell-pushforward-and-observer-homogenization
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2456-two-root-replica-uniform-offset-boundary
  - THM-2457-complete-atom-root-cosupport-graph-and-semantic-word-hostile
script: 04-computation/lrc14_aggregate_endpoint_restoration_thm2452.py
output: 05-knowledge/results/lrc14_aggregate_endpoint_restoration_thm2452.out
script_sha256: 680f68e07d8228c7bc07ebc3ffc6c2a5dba8817c3206dbd4b4f1c7672f8cbf4e
output_sha256: 2f96dcf3165363c8b4b473415238f069757982a07196622c999af653ecd55295
hash_basis: working-tree bytes (LF)
---

# THM-2452 -- aggregate first, then copy the endpoint mask

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2448 preserves a preselected exact `(X,m)` while restoring the
missing right factors. That stringent order produces a finite cospan:
the selected coefficient can move entirely into an off-diagonal mask
transition.

Here the selection order is reversed. Refine the indicator layer to a
complete Boolean atom, copy that same atom to the bare endpoint, and
only then extract an exact Fourier address. Pointwise idempotence makes
the copied endpoint free:

```text
complete mask on the present leg
  -> copy the same mask to the bare leg
  -> the whole nonnegative table is unchanged
  -> reselect a fresh exact X and 91-unit m.                    (1)
```

This removes the omitted-local-mask debt at the aggregate level. It
does not identify the complete literal mask with a semantic THM-2305
owner/repair word.

## 1. Boolean-idempotent cospan slide

Fix one lawful target/deep twist `theta`. Let `E_theta` be the bare
Boolean endpoint envelope and let

```text
F_theta E_theta=F_theta                            (2)
```

be the exact pointwise support identity for the present packet,
including any fixed Boolean terminal word. Let
`Delta_r` be the retained deepest probe. Suppose

```text
{P_omega:omega in Omega}
```

is a finite Boolean partition, co-shifted by exactly the same
THM-2350/2365 target action on both endpoint copies:

```text
sum_omega P_omega=1,

P_omega P_nu=0              for omega!=nu,

P_omega^2=P_omega                                  (3)
```

almost everywhere. Define the present-refined and two-sided tables

```text
H_omega(r,theta)
 =integral_T F_theta P_omega Delta_r,

H_(omega,nu)(r,theta)
 =integral_T
    F_theta P_omega E_theta P_nu Delta_r.           (4)
```

Equations (2)--(3) give the exact boundary identities

```text
H_(omega,nu)=0                     if omega!=nu,

H_(omega,omega)=H_omega,

sum_omega H_omega=H.                              (5)
```

Thus a complete mask already present on the moving leg may be copied
to the bare leg without changing its whole nonnegative table. Every
matched table retains the deepest diagonal zero, target quotient,
graft, arithmetic labels, and any factors already in `F_theta`.

The order is load-bearing. Equations (3)--(5) are indicator-boundary
identities. Each **whole** table in (4) is Poisson--Abel smoothed only
after the Boolean products are formed. No identity
`P_(omega,rho)^2=P_(omega,rho)` is asserted for separately smoothed
factors.

## 2. The complete 128-mask bank

Use THM-2445 notation. For the five guard/unit roles other than `q_*`,
write

```text
g_i=1-d_i,                  1<=i<=5.
```

For the source and other nondeep blocker use the safe/danger bits in
`B_sigma`. A complete local atom is

```text
P_omega
 =product_(i=1)^5 d_i^(epsilon_i)g_i^(1-epsilon_i)
  d_j^(epsilon_j)g_j^(1-epsilon_j)
  d_a^(epsilon_a)g_a^(1-epsilon_a),                 (6)

omega in {0,1}^7.
```

There are

```text
2^5*2^2=128                                             (7)
```

atoms. They satisfy (3) after every lawful target shift. The grafted
`q_*`, deepest `c_3` probe/safe pair, reduced `7`-unit quotient, and
terminal word are not split.

The lexicographic THM-2445 first-failure cells contain respectively

```text
4,64,32,16,8,4                                      (8)
```

complete atoms for labels `i=0,1,2,3,4,5`, counting the four blocker
words. For one fixed blocker word, the numbers of complete guard/unit
extensions are

```text
N_i=1,16,8,4,2,1.                                  (9)
```

## 3. Global energy-first form

Let `mathcal P` be the THM-2365 circulant projection and put

```text
mathcal Q=I-mathcal P,

D_0=||mathcal Q T_0||_2^2>0                         (10)
```

as in THM-2445. Refining the moving endpoint by (6) gives

```text
T_0=sum_omega T_omega.                              (11)
```

Hence

```text
sqrt(D_0)
 <=sum_omega ||mathcal Q T_omega||_2.
```

Some complete atom satisfies

```text
D(T_omega)>=D_0/128^2=D_0/16384.                   (12)
```

By (5), inserting that same complete atom on the partial bare endpoint
does not change `T_omega` or (12).

THM-2365 equation (19) now gives

```text
sum_(a!=0,(b,a+h)!=(0,0))
 |B_omega(a,b,h)|^2
 >=D_0/(13*16384)
 =D_0/212992.                                      (13)
```

There are

```text
12*(13^2-1)=2016                                   (14)
```

eligible colours. Consequently some complete matched table has

```text
|B_omega(a,b,h)|^2>=D_0/429391872.                 (15)
```

These constants concern the bare word `1`. They cannot be copied
unchanged to an arbitrary delayed word.

## 4. Coefficient-first adaptive form

The same mechanism preserves a chosen finite colour more efficiently.
Fix a THM-2445 cell

```text
alpha=(i,sigma)
```

and a nonzero finite target/deep coefficient

```text
B_alpha(a,b,h)!=0.                                  (16)
```

Refine only the tail factors omitted by the first-failure word `S_i`.
The complete extensions `tau` of `alpha` satisfy

```text
T_alpha=sum_(tau extends alpha)T_tau,

B_alpha(a,b,h)
 =sum_(tau extends alpha)B_tau(a,b,h).              (17)
```

There are `N_i` extensions from (9), so some `tau` obeys

```text
|B_tau(a,b,h)|>=|B_alpha(a,b,h)|/N_i.               (18)
```

Copy `P_tau` to the bare endpoint. Equation (5) preserves the same
finite colour exactly. In particular:

- the all-safe ghost `i=0` loses nothing;
- the `D(O)>0` source-owner cell `i=0` loses nothing;
- in the THM-2442 circulant ghost, the source bit is handled separately:
  each lawful phase already inserts its translated source danger on
  both endpoint copies. Apply idempotence only to the residual
  source-neutral five guard/unit bits and other blocker, termwise in
  the source phase, before the finite `C_7` transform. The selected
  nonzero `kappa`, `(a,b,h)`, graft, and word all survive. This
  phase-dependent source-danger mask is not the old stationary
  source-safe atom.

Equation (18) is stronger than the global bound when a useful
THM-2445 colour has already been selected. The energy-first form is
stronger when no prior cell or colour must be retained.

## 5. Delayed-word restoration

Let `Q_term` be a fixed positive rational terminal word, with mean
`mu(Q_term)>0`, and let

```text
V_Q=Var(Q_term),

V_A=max_(r,s,t)Var(A_(r,s,t))
```

for the selected complete bare atom. THM-2365 equation (28) gives

```text
sqrt(D_word)
 >=mu(Q_term)sqrt(D_bare)-V_Q V_A/(12R).            (19)
```

Thus bare positivity implies word-retaining positivity at every
sufficiently large clock. If the explicit stronger threshold

```text
V_Q V_A/(12R)
 <=(mu(Q_term)/2)sqrt(D_bare)                        (20)
```

holds, then

```text
D_word>=mu(Q_term)^2 D_bare/4.                      (21)
```

Combining (12)--(15) with (21) gives the quantitative word-retaining
floors

```text
D_word
 >=mu(Q_term)^2 D_0/65536,

eligible energy
 >=mu(Q_term)^2 D_0/851968,

max eligible |B_word|^2
 >=mu(Q_term)^2 D_0/1717567488.                     (22)
```

Without (20), only eventual nonvanishing is claimed. For the
source-completed ghost, THM-2442 supplies its sharper
same-finite-coefficient qualitative restoration and preserves the
source colour.

## 6. Reselecting an exact address

Let the selected complete matched table have an eligible coefficient

```text
B_tau(a,b,h)!=0,

a!=0,              q=(b,a+h)!=(0,0).                (23)
```

THM-2365 equations (21a)--(23) expand it with absolutely iterated
`m`-then-`X` meaning:

```text
B_tau(a,b,h)
 =sum_(m=a mod 13)sum_X A_(tau;X,m)(q).             (24)
```

Therefore some exact `m,X` satisfy

```text
A_(tau;X,m)(q)!=0.                                  (25)
```

Here `m` is nonzero modulo `13`; the deepest Fourier factor vanishes
when `7|m`. Hence

```text
gcd(m,91)=1.                                        (26)
```

The complete literal endpoint mask, finite colour, word, graft, and
ghost source colour when present all survive. The exact `(X,m)` in
(25) is freshly selected. It need not equal the address chosen before
endpoint restoration in THM-2448.

## 7. Why every transition disappears after full-X recombination

There is an equivalent Fourier statement. Let `L_tau` and `R_nu` be
endpoint functions carrying complete masks `P_tau` and `P_nu`. For
each exact `m` and target twist, endpoint Parseval is absolutely
convergent and gives, up to the fixed convention-dependent modulation
sign,

```text
sum_X
 L_tau_hat(X)conjugate(R_nu_hat(X+m c_3))

 =integral_T
   L_tau(x)conjugate(R_nu(x))e(+-m c_3 x)dx.         (27)
```

If `tau!=nu`, the right side is zero by (3). Thus every THM-2448
off-diagonal transition cancels after the complete `X`-sum, for each
fixed `m`, twist, word, and source phase. Only diagonal complete masks
remain. Equation (24) then reselects one of their exact terms.

The cancellation does not hold term by term in `X`. On the four-point
circle take `L=delta_0`, `R=delta_1`, and right frequency `X+1`.
Every one of the four normalized cross terms is nonzero, while their
sum is

```text
(i-1-i+1)/16=0.                                     (28)
```

The interval version is equally sharp: a danger comb `P=d_L` and its
safe complement `1-P` can have a nonzero fixed-`X` cross coefficient,
although their full-`X` sum is zero. A common translation merely
rotates that coefficient. Pairing orientations or squaring it loses
the target charge.

Thus the theorem's controlled forgetting is necessary:

```text
keep old X,m  -> THM-2448 transition cospan;

release old X,m
  -> complete matched endpoint
  -> fresh exact X and 91-unit m.                  (29)
```

## 8. Repair-current scope

The complete mask records every literal guard/unit safe-or-danger bit
and retains the deepest safe factor. This theorem makes no
THM-2379 repair-current conclusion from that record. Applying
THM-2379 is a separate corollary whose physical carrier hypotheses
must be checked in the chosen branch. Even when that check succeeds,
its repair Fourier selection need not share the target colour `q` or
exact `(X,m)` from (25), and its repair shift is not thereby a lawful
target dipole.

## 9. Exact remaining boundary

The theorem proves a complete **literal local factor mask** on both
endpoint legs. It does not by itself prove:

- that a blocker truth word is the semantic THM-2305
  source/expiration/pure-fork word;
- that a literal guard/unit danger bit has already been routed as the
  particular deletion carrier used in a target-aligned THM-2379
  current;
- that the complete mask has been physicalized as THM-2401's
  canonical `C_13` root section and orientation;
- that the old marked triangle or all-coordinate `91`-unit relation
  address survives.

Conditionally, once the retained `q_*`/`c_3` pair is identified as
THM-2401's disjoint clean pair in one canonical root section, the
complete `P_tau` is exactly a common Boolean filter and the positive
table supplies simultaneous descendant mass. THM-2401 then applies.
That semantic/root-section identification is not proved here.

Accordingly, the omitted local masks are no longer the aggregate
obstruction. The live problem is owner/repair/root semantics, not
endpoint enumeration or fixed-frequency transition cancellation. No
scalar profile is excluded; the ledger remains `165`.

## 10. Exact companion

Run

```text
python 04-computation/lrc14_aggregate_endpoint_restoration_thm2452.py
python -O 04-computation/lrc14_aggregate_endpoint_restoration_thm2452.py
```

The companion:

- enumerates all `128` complete masks and `16,384` ordered mask pairs;
- verifies the lexicographic extension census (8)--(9);
- checks the drift, eligible-energy, coefficient, and delayed-word
  denominators;
- checks every `91`-unit CRT residue;
- verifies an exact four-term fixed-`X` leakage/full-`X` cancellation
  hostile;
- uses explicit `require` checks which remain active under `-O`.

The normal, optimized, and stored transcripts are byte-identical.
LF-normalized hashes are recorded in the frontmatter.

## 11. Audit boundary

Two independent immutable hostile audits accepted the
Boolean-idempotent slide, quantitative floors, aggregate-to-exact
reselection, word threshold, reduced ghost-bank/source-phase typing,
repair boundary, and fixed-`X` hostile. Both independently replayed
the normal, optimized, and stored transcripts and reproduced the
recorded LF hashes.
