# The r=5 contraction distinguishes the two root charts

**Status: VERIFIED-EXACT sidecar to proved THM-2594; not independently
audited and not a bridge/current theorem.**

## 1. Question after the fixed-root hostile

The repaired THM-2594 proves that both the derived-label contraction and the
genuine fixed-absolute-root control are nonzero.  This refutes the claim that
`theta=t-2u` is the unique cause of nonvanishing, but leaves two sharper
questions:

1. does the contraction distinguish the two root charts at all; and
2. after retaining owner root `u`, are there modes that are pure in neither
   chart?

The closest proved mechanism is THM-2594's common-base joint table.  The
canonical hostile is its fixed-root control.  The corrected near miss is
MISTAKE-295: a constant theta column is not an absolute-root test.  The
least-used sidecar is the owner root before the final marginalization.

## 2. The exact shear contrast

Let `J_(u,q,ell,theta)` be the audited nonnegative integer joint table.  Form

```text
R_theta(ell,theta)=sum_(u,q) J_(u,q,ell,theta),
R_t(ell,t)=sum_(u,q) J_(u,q,ell,t-2u).                  (1)
```

Double-center each `7 x 13` response as in THM-2512 and subtract.  The raw
tables differ in `19/91` cells and their centred tables in `28/91` cells.
The fixed primitive coefficient `Psi_(1,1)(1,1)` of the contrast is nonzero
in `Q(zeta_91)`, with `67/72` nonzero reduced coordinates.  More strongly,

```text
all 5,184/5,184 primitive contrast coefficients are nonzero,              (2)
minimum reduced support=60/72.
```

Thus the toothpick contraction genuinely distinguishes “sum in theta, then
forget u” from “return to t, then forget u.”  Equation (2) proves chart
sensitivity, not that either chart is the physical cause of the original
signal.

## 3. Retaining u exposes the shear law

Before marginalizing `u`, put

```text
S_ell(u,theta)=sum_q J_(u,q,ell,theta),
A_ell(u,t)=S_ell(u,t-2u).                                (3)
```

For a primitive 13th root `xi`, use the convention

```text
S_hat(r,s)=sum_(u,theta) S(u,theta) xi^(-ru-stheta).
```

Changing variables in (3) gives the exact Fourier shear

```text
A_hat(r,s)=S_hat(r+2s,s).                                (4)
```

This identity is checked on all `7*13^2` modes at the split prime
`53`, with `xi=16`.

The two exceptional lines now have a precise meaning in the slaved chart:

```text
r=0:    a profile depending only on theta,
r=2s:   a profile depending only on absolute t=theta+2u. (5)
```

Synthetic pure-theta and pure-t controls have zero support off their named
line.  Therefore any nonzero mode with

```text
s!=0,       r!=0,       r!=2s                           (6)
```

is mixed in both charts; no choice between theta and t makes it a
one-coordinate profile.

## 4. Exact verdict

At the certified split embedding, the audited table has

```text
576/1008 nonaxial modes in each chart,
528/924 modes nonzero off both pure-root lines.           (7)
```

The off-both-line census by word slot is

```text
(0,0,132,132,132,132,0).                                 (8)
```

Thus the four interior word slots carry a genuine owner-root/root-label
interaction.  Independently, after integer double centering in `(u,root)`,
every one of the `169` cells is nonzero for slots `1,...,6`; slot zero is
identically zero.  The split-prime statement (7) is a certified lower bound
on exact cyclotomic support: nonzero reduction proves exact nonvanishing,
while modes vanishing at this one embedding are not claimed zero over the
integers.

## 5. Consequence boundary

This repairs the causal picture as follows:

```text
old claim:  nonzero because theta is slaved                 REFUTED,
survivor:   theta and t contractions are both nonzero       PROVED THM-2594,
new fact:   their contrast and off-both-line modes survive  VERIFIED-EXACT.
```

The new modes show that the table is not merely a pure derived-root profile
or a pure absolute-root profile.  They still live in a finite common-base
candidate contraction.  No physical phase/current, THM-2512 generic bridge,
THM-2449 one-point response, U_full ancestry lift, row exclusion, or LRC(14)
conclusion follows.

## Reproduction

```bash
python3 04-computation/lrc14_stage2_root_shear_contrast_probe_20260816.py
python3 -O 04-computation/lrc14_stage2_root_shear_contrast_probe_20260816.py
```
