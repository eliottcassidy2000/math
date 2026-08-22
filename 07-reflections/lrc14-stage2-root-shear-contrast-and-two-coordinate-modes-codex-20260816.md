# The r=5 contraction distinguishes the two root charts

**Status: PROVED + FINITE-EXACT as THM-3505 after an independent full-table
reconstruction; not a bridge/current theorem.**

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

### Restoring the C7 word Fourier axis

The slotwise count in (8) suggests a stronger test than support by individual
`ell`.  At the split prime `547`, the primitive 91st root `64` gives

```text
zeta_7=64^13=81,             zeta_13=64^7=475.          (9)
```

Take the full Fourier transform in word slot, owner root, and root label.
For every nontrivial word mode `beta`, every nontrivial root-label mode `s`,
and every owner mode off the two pure lines, the coefficient is nonzero:

```text
6 * 12 * 11 = 792/792 primitive C7 x C13 x C13 modes fire. (10)
```

Each of the six `C7` modes contributes all `132/132` admissible root pairs.
The complete split-field support ledger has SHA256

```text
51ea8b27b8f14f07ed7099601e80a5b36f18510cf6da9c4f815802e6eb8f05cc. (11)
```

This is a genuine three-axis spectral-closure statement for the realized
candidate table.  It is not the physical `7 tensor 13` bridge sought by
LRC(14): it retains an additional owner-root Fourier coordinate and still
lacks the physical phase/current, common U_full atom key, and Hall-cone
consequence.  Its value is sharper localization: cancellation does not
reappear when the word-slot axis is Fourier transformed.

### The projective slope pencil and the owner quotient

There is a sharper organization of the same transform.  For nonzero
root-label frequency `s`, put

```text
lambda = r/s in F_13,                                      (12)
```

where `r` is owner-root frequency.  The primitive frequencies split into
thirteen affine slope classes, each one a 72-element Galois orbit indexed by
the six nonzero word modes and twelve nonzero root modes.  Exact reduction
modulo `547` gives

```text
lambda=0:  owner marginal                         72/72,
lambda=2:  fixed-absolute-root line                72/72,
lambda other: eleven genuinely two-coordinate lines 11*72/11*72,
all slopes:                                       936/936. (13)
```

Thus summing out the owner does **not** kill the signal: the resulting
`C7 x C13` mixed transform is nonzero in all 72 primitive characters.  This
is the explicit quotient check requested by the spectral-closure question.
For `lambda=0` it recovers, rather than strengthens, THM-2594's already
proved primitive mixed block; `lambda=2` is the corresponding
fixed-absolute-root hostile, and the other eleven slopes are the new
two-coordinate localization.

The full slope-pencil ledger has SHA256

```text
e92e3f1b072db16ada1daa28925803ebd9e11658deb3532680911ed637dee85d. (14)
```

This closes the algebraic primitive-spectrum gate for this one realized
candidate table.  It does not repair the remaining type mismatch: the table
is still a linked-node common-base contraction rather than THM-2449's lawful
one-point response or a positive physical current.  The missing bridge is
therefore realization/descent, not a hidden Fourier zero.

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

## Independent audit and basis boundary

The independent THM-3505 audit imports neither this Fourier helper nor the
THM-2594 constructor.  It reconstructs the full common-base table and matches
its exact digest before transforming.  It then reproduces the complete
`936/936` slope pencil and adds the anchored mixed-Haar basis.  That basis has
only `804/864` nonzero word-Fourier coordinates, despite `864/864` nonaxial
Fourier characters.  A rank-`144` exact change of basis reconciles the two
counts: support is basis-dependent.

The same audit compares the U_full minimum joint-address gate.  A formal
binary checkerboard can be embedded into a selected `(owner,root)` rectangle,
but the actual U_full cells still have no map to the required common base,
roots, sheets, horizons, or address.  Thus the independent audit promotes the
finite spectral statement while preserving the realization boundary.

## Reproduction

```bash
python3 04-computation/lrc14_stage2_root_shear_contrast_probe_20260816.py
python3 -O 04-computation/lrc14_stage2_root_shear_contrast_probe_20260816.py
```
