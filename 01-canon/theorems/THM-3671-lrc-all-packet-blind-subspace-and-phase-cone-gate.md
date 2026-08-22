---
id: THM-3671
title: "LRC all-packet blind subspace and phase-cone gate"
status: >
  PROVED + INDEPENDENTLY HOSTILE-AUDITED.  For one fixed selected source and
  target pair, the intersection of the zero-target fibres of its 120
  owner-pivot packets is exactly the pure source line when the six unit
  residues have nonzero sum, and gains one all-nonsource diagonal line when
  that sum vanishes.  Four explicit packet charts already cut out the same
  blind subspace.  Off the scalar sum-zero branch, any nonzero residue measure
  whose coefficients lie in one strict origin-centred phase cone and whose
  support has all mod-13 coordinates nonzero has a nonzero target in some
  packet.  Actual LRC currents are complex and neither the nonzero grouped
  measure nor the required phase cone is proved; no LRC(14) conclusion is
  claimed.
source: kps-s193 / all-packet intersection after THM-3666, 2026-08-21
audit: >
  PASS AFTER SCOPE REPAIRS -- Galileo independently verified dim(K)=8, each
  packet-fibre dimension 6, both branches of the all-thirty intersection, the
  four-chart equality graph, the all-unit scalar hostile, and the phase-cone
  pushforward argument.  The audit fixed the selected source and target pair,
  excluded the zero measure, made the half-plane strict and origin-centred,
  and separated termwise all-unit address banks from a nonzero grouped
  residue measure.
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-3666-lrc-owner-pivot-dual-pair-swap-twist-basis
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2325-prescribed-target-gain-full-lattice-91-unit-needle-bank
  - THM-3670-lrc-successor-mass-transfer-and-thirteen-number-gate
---

# THM-3671 -- four pair-swap charts expose every nonsource residue

**PROVED + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem identifies the
exact algebraic object that every owner-pivot packet misses.  It also isolates
the minimal positivity/phase sidecar that would make the all-unit projector
useful.

## 1. Packet zero fibres

Fix one THM-2309 selected source `j` and its target pair `{a,b}`.  Work over
`F=F_13` with their nine owner-pivot labels

```text
U={0,1,2,3,4,5},
B={j,a,b}.                                          (1)
```

The six labels in `U` are units modulo thirteen and the three blocker labels
have zero speed residue:

```text
w_u!=0 for u in U,
w_j=w_a=w_b=0.                                     (2)
```

Let

```text
K=w^perp subset F^9.                               (3)
```

For each ordered pair of distinct grafts `k,l in U`, THM-3666 gives the
target-dual plane

```text
T_(k,l)=span(e_a-e_k,e_b-e_l).                     (4)
```

Its zero-target fibre in the residue space is

```text
L_(k,l)=T_(k,l)^perp intersection K
       ={r in K:r_a=r_k and r_b=r_l}.              (5)
```

The omitted-unit choice does not affect (4), so these thirty planes represent
all 120 owner packets associated with this fixed source and target pair.

## 2. Exact all-packet intersection

Put

```text
sigma_U=sum_(u in U) w_u in F.                     (6)
```

Then

```text
intersection_(k!=l) L_(k,l)
 =span(e_j),                         if sigma_U!=0,

 =span(e_j,v_diag),                  if sigma_U=0,  (7)
```

where `v_diag` is one on `U union {a,b}` and zero at `j`.

Indeed, if `r` lies in every fibre, then for any two unit labels `k,k'`
choose `l` distinct from both and use (5) to get

```text
r_k=r_a=r_(k').                                    (8)
```

Thus all six unit coordinates equal a common value `d`; the same argument on
the second graft gives `r_b=d`.  Since the blocker speed residues vanish,
the remaining equation `r.w=0` is simply

```text
d sigma_U=0,                                       (9)
```

while `r_j` is free.  Equations (7) follow.

Only four charts are needed.  The constraints from

```text
(k,l)=(0,1),(2,3),(4,5),(1,0)                     (10)
```

form a connected equality graph on the eight labels `U union {a,b}`; hence
their common zero fibre already gives (8)--(9).  This four-chart statement is
minimal only as an explicit sufficient bank; no cardinal minimality claim is
made.

## 3. The source and scalar hostiles are real

Every target-dual vector in (4) vanishes at `j`.  Therefore the pure source
character

```text
ell -> zeta_13^(ell_j)                              (11)
```

is constant on every packet plane.  This is the smallest exact witness that
target packets cannot detect source ancestry by themselves.

When `sigma_U=0`, every vector

```text
c e_j+d v_diag                                     (12)
```

lies in all packet zero fibres.  Taking `c,d` nonzero gives a residue whose
all nine coordinates are nonzero, yet which is target-zero for every packet.
Thus support on mod-13 units alone cannot bypass the scalar branch (6).

When `sigma_U!=0`, (7) has no all-coordinate-unit vector: every vector on the
source line has eight zero coordinates.  This sharp dichotomy is purely
mod-thirteen and should not be confused with septimal word charge.

## 4. Positivity or one phase cone would activate the all-unit support

Let `mu` be a nonzero finite complex residue measure on `K`, supported on
vectors with all nine coordinates nonzero.  For a packet `(k,l)`, push `mu`
through

```text
pi_(k,l)(r)=(r_a-r_k,r_b-r_l).                     (13)
```

Assume either

1. every nonzero coefficient of `mu` is a positive real number; or
2. there is one angle `theta` such that every nonzero coefficient `z` obeys
   `Re(exp(-i theta)z)>0`, so all coefficients lie in one strict open
   half-plane through the origin.

In either case, every nonempty fibre aggregate is nonzero: a sum of vectors
in one open half-plane cannot vanish.  If all thirty packet pushforwards were
supported only at target zero, every support point of `mu` would belong to
the intersection (7).  For `sigma_U!=0` this contradicts all-coordinate-unit
support.  Hence

```text
sigma_U!=0 and a common phase cone
  => some packet has a nonzero target aggregate.    (14)
```

By (10), four packet charts suffice in (14).

THM-2325/2331 provide enormous all-`91`-unit term/address banks, whose
reductions have all coordinates nonzero.  They do **not** prove that the
corresponding grouped all-unit residue measure is nonzero, nor that any
aggregate coefficient survives.  Thus neither the nonzero-measure premise
nor the common phase-cone premise of (14) is currently supplied for the
actual marked current.  The coefficients include endpoint, deep and terminal
phases, so (14) is a conditional transfer gate, not a completed current
theorem.  A non-Archimedean alternative would be one uniquely minimal
valuation in a nonzero packet fibre, but no such valuation sidecar is proved
here.

No visible-height or terminal-phase assertion is made, and LRC(14) remains
open.  **QED.**
