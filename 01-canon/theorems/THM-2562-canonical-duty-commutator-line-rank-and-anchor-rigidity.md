---
id: THM-2562
title: "Canonical duty-commutator line rank, zero-anchor rigidity, and six-replica covariance debt"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; independent hostile audit REQUESTED.
  The canonical THM-2556 quotient-duty commutator has exact rational ranks
  50, 70, 72, or 74 according to four gain-slope classes.  Its only local
  rank defect is one cyclic puncture syzygy.  More decisively, every
  commutator-null target density vanishes at target zero.  THM-2541's actual
  canonical Abel current has nonzero zero-target aggregate, so its quotient
  duty curvature is nonzero for every gain and equals one forced scalar on
  six consecutive nonzero target fibres.  Exact cancellation would require
  a sixfold covariance replica.  Nonnegative flat-arc hostiles remain sharp;
  no physical common-carrier intertwiner or residual exclusion is proved,
  no row is removed, and LRC(14) remains OPEN.
source: codex-2026-07-27-holotopy
depends_on:
  - THM-2541-canonical-typed-row-full-target-plane-support
  - THM-2544-all-unit-target-projector-kernel-and-lawful-image-intersection-obstruction
  - THM-2556-reynolds-duty-curvature-and-fibre-covariance-mixed-cell
related:
  - THM-2548-seven-step-c91-transfer-and-full-norm-separation
  - THM-2550-canonical-typed-row-double-nondegeneracy
  - THM-2557-doubly-centered-transfer-moment-index-and-hall-cone-boundary
script: 04-computation/lrc14_duty_commutator_line_rank_thm2562.py
output: 05-knowledge/results/lrc14_duty_commutator_line_rank_thm2562.out
---

# THM-2562 -- the canonical duty square forces a six-replica cancellation

THM-2556 decomposes the first non-product transfer/projector square into

```text
quotient duty curvature + within-fibre covariance curvature.             (1)
```

Its exact norm invoice still permits cancellation.  The present theorem
computes the quotient operator itself.  The result is much more rigid than
mere nonvanishing: the actual THM-2541 canonical current forces six identical
uncontaminated outputs for every nonzero gain.

The word *holotopy* is again mnemonic.  The rank defect below is a finite
mapping-cone syzygy, not a claim about a topological homotopy group.

## 1. Canonical duty commutator

Let `G=F_13^2`, write `q=(x,y)`, and use THM-2556's unnormalized canonical
duty count

```text
n(x,y)=2316060
 +210552(1_(x=0)+1_(y=0))
 +12 1_(x+y=0)+19128 1_(x=y=0).                            (2)
```

The normalized duty is

```text
nu=kappa n,       kappa=N_7/F=1439676/27825593350009.       (3)
```

For `g!=0`, define

```text
(tau_g a)(q)=a(q-g),       C_g=sum_(j=0)^6 tau_g^j,
K_g=[M_n,C_g].                                                  (4)
```

Then, exactly,

```text
(K_g a)(q)=sum_(j=0)^6 (n(q)-n(q-jg))a(q-jg).               (5)
```

The normalized THM-2556 quotient curvature is `kappa K_g`.

## 2. Puncture factorization on one gain line

Every affine `g`-line has thirteen points.  Index one such line by `C_13`
and put

```text
(Hf)(i)=sum_(j=0)^6 f(i-j).                                 (6)
```

Subtract the constant value of `n` on the complement of its punctures on
this line.  If `S` is the puncture set, `D` the diagonal matrix of its
nonzero excess weights, `E_S` restriction, and `I_S` inclusion, then

```text
K_(g,L)=I_S D E_S H-H I_S D E_S.                            (7)
```

Thus each puncture contributes at most one incoming and one outgoing rank:

```text
rank K_(g,L)<=2|S|.                                         (8)
```

There is one exact exception to full incidence rank.  For three punctures,
write their cyclic gaps in increasing order.  The two rectangular incidence
matrices `H_(T,S)` and `H_(S,T)` have rank three for all fourteen gap types
except

```text
(1,6,6),                                                    (9)
```

where both have rank two.  This follows by checking one explicit `3 x 3`
minor for each of the other thirteen types; the dependency in (9) is
immediate from the two coincident half-circle boundaries.

For the representative `S={0,6,12}` with weights `(x,y,z)`, the remaining
one-dimensional Schur obstruction is

```text
-z(x-z)/x.                                                  (10)
```

Consequently the block has rank four when the unique adjacent puncture pair
has equal weights, and rank five when those weights differ.  In the equal
canonical case `(210552,12,210552)`, the extra signed circuit restricts to
`(-1,+1)` on the adjacent equal punctures.  Its opposite signs will be
load-bearing for the nonnegative boundary below.

## 3. Complete rational rank spectrum

Let `m=g_y/g_x` when `g_x!=0`.  The three duty walls in (2) are

```text
x=0,                    y=0,                    x+y=0.       (11)
```

If `g` is tangent to one of them, one line has one puncture and twelve lines
have two.  If `g` is transverse, the origin line has one puncture and the
other twelve have three with weights `(210552,210552,12)`.  The exceptional
gap pattern (9) occurs on exactly two affine lines precisely when

```text
(m-1)(m^2-4m+1)=0.                                        (12)
```

Over `F_13`, the quadratic roots are `6,11`.  For `m=1` the adjacent pair
has equal weights and (10) vanishes; for `m=6,11` it has unequal weights.
Therefore:

| gain class | gains | line-block ranks | `rank K_g` | nullity |
|---|---:|---:|---:|---:|
| axes or antidiagonal | 36 | `2,4^12` | 50 | 119 |
| `g_y=g_x` | 12 | `2,4^2,6^10` | 70 | 99 |
| `g_y/g_x in {6,11}` | 24 | `2,5^2,6^10` | 72 | 97 |
| all remaining gains | 96 | `2,6^12` | 74 | 95 |

The table is over `Q`; multiplying by `kappa` does not change it.  The
`(1,6,6)` puncture circuit is the unique local rank loss.

## 4. Zero-target anchor rigidity

The gain line through zero is simpler than the full table.  Its duty value
is constant at every nonzero multiple of `g`, while

```text
d_g=n(0)-n(g)
 =229692   if g_x g_y=0,
 =440232   if g_x g_y!=0 and g_x+g_y=0,
 =440244   otherwise.                                      (13)
```

For every target density `a` and every `j=1,...,6`, the backward half-window
at `jg` contains zero exactly once and otherwise only equal-duty points.
Hence (5) gives the uncontaminated identity

```text
(K_g a)(jg)=-d_g a(0),              j=1,...,6.              (14)
```

In particular

```text
K_g a=0  implies  a(0)=0,
||K_g a||_2>=sqrt(6)d_g|a(0)|.                              (15)
```

This holds for arbitrary real or complex signed densities, not merely for a
positive cone.

THM-2544 equation (27) identifies the actual THM-2334 Abel boundary current
`c` with `Uc=A`.  THM-2541 proves that every one of its 169 target
aggregates is nonzero; in particular `A(0)!=0`.  Therefore the canonical
THM-2556 quotient-duty term is nonzero for **every** `g!=0`, already on the
six nonzero targets `g,...,6g`:

```text
([M_nu,C_g]A)(jg)=-kappa d_g A(0).                          (16)
```

Quotient nonvanishing is no longer a debt on this typed current.

## 5. Exact six-replica covariance invoice

Assume, conditionally as in THM-2556, that a lawful common-carrier
intertwiner identifies the desired physical transfer with `D_v`, where
`Q(v)=g`.  Put

```text
W_g=(R D_v-C_gR)c.                                         (17)
```

If the full mixed face vanishes, equations (1) and (16) force

```text
W_g(jg)=+kappa d_g A(0),              j=1,...,6,             (18)

||pr_{g,...,6g} W_g||_2
 =sqrt(6)kappa d_g|A(0)|.                                  (19)
```

Thus cancellation is not an unspecified phase accident.  The within-fibre
covariance must reproduce one exact complex scalar on six consecutive target
fibres.  Any strict projected residual bound below (19) forces a nonzero
all-unit mixed face.

The common-carrier hypothesis remains open.  THM-2550's non-replica result
uses a different lawful clock from its drift packet and is not yet typed to
`W_g`; it is a promising comparison control, not a proof of (19)'s failure.

## 6. Sharp nonnegative flat-arc boundary

For each gain define

```text
F_g={q:n(q)=n(q+g)=...=n(q+6g)}.                            (20)
```

A coordinate vector `delta_q` is a zero column of `K_g` exactly when
`q in F_g`.  The complete nonnegative kernel is

```text
ker K_g intersection R_(>=0)^G
 =cone{delta_q:q in F_g}.                                  (21)
```

To see this, on every full-rank or rank-five line block the kernel is exactly

```text
f|_S=0,                         (Hf)|_S=0.                  (22)
```

For `f>=0`, (22) kills every backward half-arc ending at a puncture, leaving
exactly the coordinates (20).  The rank-four block has one extra circuit,
but its restrictions to the two adjacent equal punctures have opposite
signs; nonnegativity forces its coefficient to vanish, reducing again to
(22).

The exact flat-arc counts are:

| direction | gains | `|F_g|` |
|---|---:|---:|
| axes or antidiagonal | 36 | 36 |
| slopes `1,3,6,9,11` | 60 | 18 |
| slopes `2,4,5,7,8,10` | 72 | 16 |

Thus every gain has at least sixteen sharp nonnegative fibre-uniform
hostiles: for `q in F_g`, the lift `c=L delta_q` has `Rc=RD_vc=0` and zero
full mixed face.  None contains the origin, consistently with (15).  Every
strictly positive target density is detected.  The condition `a(0)!=0` is
sufficient but not necessary; densities with `a(0)=0` include both live and
flat cases, while signed kernels have the much larger nullities in Section 3.

## 7. Holotopy boundary and next test

The puncture formula (7) is a finite mapping cone: each duty-wall puncture
contributes incoming and outgoing rank, the gap `(1,6,6)` is the sole local
two-syzygy, and unequal adjacent weights contribute a nonzero holonomy which
kills half that syzygy.  This is the first exact mixed-cell spectrum on the
canonical target plane.

It still lives after quotienting to `G`.  The next decisive test is to build
the common chart x semantic x ancestry carrier, evaluate its actual
within-fibre residual on `g,...,6g`, and ask whether one lawful packet can
satisfy the replica law (18).  Failure at one coordinate closes that gain;
success identifies the precise sixfold cancellation mechanism to attack.

No scalar cover, physical later-root map, or Hall arrival follows here.  The
typed row is not a scalar cover, no row is removed, and LRC(14) remains open.

## 8. Exact companion

The dependency-free companion enumerates all `2,184` affine gain-line
blocks and computes their rational ranks independently from the puncture
factorization.  It checks all 286 three-puncture incidence types, every
exceptional slope and Schur boundary, the complete rank/nullity and flat-arc
histograms, the kernel constraints and signed circuit, and all `1,008`
uncontaminated anchor rows.  Run

```bash
python3 04-computation/lrc14_duty_commutator_line_rank_thm2562.py
python3 -O 04-computation/lrc14_duty_commutator_line_rank_thm2562.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_duty_commutator_line_rank_thm2562.out
```

after LF normalization.  Every executable check raises explicitly under
optimized Python.
