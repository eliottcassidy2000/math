---
id: THM-3670
title: "LRC successor-mass transfer and thirteen-number gate"
status: >
  PROVED; PENDING INDEPENDENT HOSTILE AUDIT.  For any owner-pivot packet, a
  nonzero three-site defect in THM-2365's nonnegative successor mass forces
  positive H-drift, hence some exact 91-unit frequency triangle with a
  nonzero target fibre, hence a nonzero physical three-pair-swap defect.
  Across all 120 packet choices the sufficient test depends on only thirteen
  rational masses.  Simultaneous failure on all thirty graft charts is
  equivalent to six equal a-graft masses, six equal b-graft masses, and one
  scalar relation.  The theorem does not retain the preselected THM-2327
  triangle, the all-coordinate 91-unit projector, visible height or terminal
  phase, and does not prove LRC(14).
source: kps-s193 / THM-2365-to-THM-3665 transfer composition, 2026-08-21
depends_on:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-3665-lrc-support-minimal-three-twist-target-detector
  - THM-3666-lrc-owner-pivot-dual-pair-swap-twist-basis
related:
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
  - THM-3669-lrc-typed-control-all-packet-three-twist-defects
---

# THM-3670 -- the covering transfer is a thirteen-number obstruction

**PROVED; PENDING INDEPENDENT HOSTILE AUDIT.**  This theorem composes an older
nonnegative all-frequency observable with the new support-minimal target
detector.  It replaces “control one preselected complex triangle” by a
rational mass test that is allowed to choose its own exact triangle.

## 1. Owner-pivot sign convention and successor mass

Let `a,c` be the two target blockers, with `c` the deepest blocker, and choose
distinct graft units `k,l` in an owner-pivot packet.  THM-3666 gives the dual
dipoles

```text
alpha=e_a-e_k,
beta =e_c-e_l.                                      (1)
```

Use THM-2365's sign convention: its present factors at `(s,t)` are translated
by `-s alpha/13-t beta/13`.  Thus its coordinate plane is the physical twist
plane with ordered basis `(-alpha,-beta)`.

Let

```text
K_(k,l)(r,s,t)
 =integral_T F_(s,t)(u) Delta_r(u) du >=0,          (2)

S_(k,l)(s,t)=sum_(r in F13) K_(k,l)(r,s,t).         (3)
```

Here `F_(s,t)` is the complete present/transported-word product and
`Delta_r` is the shifted deepest-danger factor of THM-2365.  Its exact
successor identity is

```text
S_(k,l)(s,t)
 =integral_T F_(s,t)(u)(2-d(13 c u)) du.            (4)
```

Every quantity in (4) is a nonnegative rational number for the LRC interval
sets.  The word factors are unchanged by `(s,t)` because the clock is
divisible by thirteen.

## 2. Successor variation forces a physical pair-swap defect

Suppose

```text
S_(k,l)(0,0)+S_(k,l)(-1,0)-2S_(k,l)(0,-1) !=0.     (5)
```

Then `S_(k,l)` is nonconstant.  If THM-2365's drift energy `D_K` vanished,
its exact projection theorem would give

```text
K_(k,l)(r,s,t)=G(r-t).                              (6)
```

Summing (6) over `r` makes (3) independent of `(s,t)`, contradicting (5).
Therefore

```text
D_K>0.                                             (7)
```

THM-2365 (19), (21a)--(23) now extracts integers `X,m` and a nonzero target
vector `q` such that

```text
gcd(m,91)=1,
A_(X,m)(q)!=0.                                     (8)
```

This is an exact fixed-frequency THM-2334 target fibre, not merely a formal
colour.  Define its complex twist profile in the same sign convention by

```text
J_(X,m)(s,t)=H_(X,m)(-s alpha-t beta).              (9)
```

Finite Fourier inversion in THM-2334 and THM-3665 say that (8) makes (9)
nonconstant and hence forces some centre `(s,t)` with

```text
J(s,t)+J(s-1,t)-2J(s,t-1) !=0.                     (10)
```

Because the ordered twist basis in (9) is `(-alpha,-beta)`, the three terms
in (10) are precisely the two reversed blocker/graft pair shifts supplied by
(1).  THM-3666 therefore turns (10) into a concrete physical pair-swap
defect.  This proves the transfer

```text
nonzero rational successor defect
  => positive nonnegative H-drift
  => some exact 91-unit nonzero target triangle
  => some physical three-pair-swap defect.          (11)
```

The centre in (10), and the extracted `(X,m)`, need not be the centre or
triangle used in (5).

## 3. All 120 packets reduce to thirteen masses

Let the six guard/unit labels be `U={0,...,5}`.  Put

```text
S_0=S_(k,l)(0,0),
A_k=S_(k,l)(-1,0),
B_l=S_(k,l)(0,-1).                                  (12)
```

The first value is independent of both grafts.  At `(-1,0)` only the
`a/k` dipole moves, so `A_k` is independent of `l`; at `(0,-1)` only the
`c/l` dipole moves, so `B_l` is independent of `k`.  The omitted-unit label
does not occur in either dipole.  Hence the 120 owner packets collapse first
to the thirty legal ordered pairs `k!=l`, and their successor defects are

```text
Delta_(k,l)=S_0+A_k-2B_l.                           (13)
```

The following equivalence is exact:

```text
Delta_(k,l)=0 for every k!=l                       (14)

iff

A_0=...=A_5=A,
B_0=...=B_5=B,
S_0+A=2B.                                          (15)
```

To prove the forward direction, fix two labels `k,k'` and choose `l` distinct
from both.  Subtracting the equations `(k,l)` and `(k',l)` gives
`A_k=A_(k')`; thus all six `A` values agree.  Similarly, for any `l,l'`
choose `k` distinct from both and subtract to obtain `B_l=B_(l')`.  Any one
equation then gives the scalar relation in (15).  The reverse implication is
immediate from (13).

Combining (11) and (14)--(15) yields the operative gate:

> Unless the thirteen rational numbers `S_0,A_0,...,A_5,B_0,...,B_5` have
> exactly the rigid form (15), at least one owner-pivot packet admits some
> exact `91`-unit triangle with a nonzero physical target defect.

Thus a single graft-mass discrepancy, boundary-jump discrepancy, or exact
valuation obstruction to (15) is sufficient.  One isolated drifting graft is
not by itself enough: THM-2367 gives exact positive masks that recirculate one
pair.  The new requirement is simultaneous six-graft recirculation.

## 4. Sharp remaining covering-row lemma

THM-2367 supplies, from a genuine scalar cover in its eligible role branches,
a septimally lawful locally drifting graft.  Its hostile controls show why
that local drift can disappear after multiplication by the other owner
factors.  Equations (12)--(15) identify the exact stronger target:

```text
on a genuine scalar cover with an eligible low target,
exclude the simultaneous thirteen-mass configuration (15).         (16)
```

No current hostile realizes (15) on a genuine cover.  This is a concrete
finite recirculation problem, but it is not solved here.

The transfer (11) retains the delayed word, lawful target pair shifts, a
nonzero deep colour and an exact `91`-unit multiplier.  It does **not** retain
the previously selected THM-2327 triangle, impose all-nine-coordinate
`91`-unit support, land degree-526 visible height, or preserve the terminal
component phase.  Those are separate sidecars.  LRC(14) remains open.

**QED.**
