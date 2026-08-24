---
id: THM-3993
title: "Labelled two-wall polyhedral elimination and local top response"
status: >
  PROVED. A labelled lower envelope with exactly one rising and one falling
  active wall has an exact local elimination formula: the displaced top is at
  (B-A)/(p+q), and its value is the branch-slope-weighted response
  (qA+pB)/(p+q). For fixed slopes the law holds coefficientwise; with variable
  slopes it gives an order-r coefficient after lower offset orders vanish,
  and that coefficient is leading only when nonzero. AP13 has six isolated
  tops and at t=1/14 the weights are 13:1, but THM-2050 proves that invariants
  of the complete unperturbed local germs do not determine the global maximum.
  This is a local polyhedral analogue, not a proof of LRC(14).
source: root + jet_bridge / Hopf-product cross-frontier session, 2026-08-24
audit: >
  PASS (root, 2026-08-24). The proof checks existence and uniqueness of the
  local top, fixed- and variable-slope coefficient laws including a
  cancellation hostile, the exact AP13 equality fibre and one-sided slopes,
  the AP13/loose-lift hostile, a positive nonsmooth control, and all claimed
  consequence boundaries.
depends_on:
  - THM-2050-period14-top-germs-do-not-determine-global-loneliness
related:
  - THM-2047-phase-height-toric-arrangement-for-lrc
  - THM-3910-lrc14-auxiliary-center-erosion-and-t-sheet-variance-response
  - THM-3995-scale-two-parity-hole-support-and-integer-variance-tariff
  - THM-3990-componentwise-harmonic-obstruction-and-repair-quotient
---

# THM-3993 -- labelled two-wall polyhedral elimination and local top response

**PROVED.**  This theorem extracts the exact local max--min response suggested
by the Brendle--Hung normal-minimization mechanism.  It deliberately retains
wall owner, sign, and slope.  Its LRC application stops at a sharp local/global
boundary.

## 1. Exact two-wall elimination

Let `s` range in a neighborhood of zero.  Suppose a finite labelled lower
envelope near `(h,s)=(0,0)` is

```text
F_s(h)=min_j ell_j(h,s).                                  (1)
```

Assume every `ell_j` is continuous near the base point, exactly two walls are
active there, and write those two walls as

```text
ell_+(h,s)=M+A(s)+p(s)h,
ell_-(h,s)=M+B(s)-q(s)h,                                  (2)
```

where

```text
A(0)=B(0)=0,             p(0)>0, q(0)>0.                 (3)
```

Assume `A,B,p,q` are continuous and every other labelled wall has strict
positive slack at `(0,0)`:

```text
ell_j(0,0)>M.                                             (4)
```

Then there are `eta>0` and `s_1>0` such that, for `|s|<s_1`, every other wall
is slack on `|h|<=eta`, and `F_s` has a unique maximizer on that interval.  It
is

```text
h_*(s)=(B(s)-A(s))/(p(s)+q(s)),                           (5)
```

and its exact height is

```text
M_loc(s)=M+(q(s)A(s)+p(s)B(s))/(p(s)+q(s)).               (6)
```

### Proof

Continuity and (3) give `p(s),q(s)>0` for small `s`.  The two affine walls in
(2) meet exactly at (5).  To the left of that point `ell_+` is the smaller of
the two and is strictly increasing; to the right `ell_-` is the smaller and
is strictly decreasing.  Their lower envelope therefore has a unique maximum
at their intersection.  Substitution gives

```text
A+p(B-A)/(p+q)=(qA+pB)/(p+q)
                =B-q(B-A)/(p+q),                         (7)
```

which proves (6) for the two active walls.

The intersection tends to zero as `s` tends to zero.  By continuity, strict
finite slack (4), and finiteness of the wall set, choose a closed rectangle
`|h|<=eta, |s|<=s_1` on which every other wall lies strictly above both active
walls, while `h_*(s)` lies in its interior.  Hence `F_s` equals the two-wall
envelope there and has the asserted unique maximizer.  This proves (1).
The signed slopes and the attachment of each offset to its branch identify the
inputs of (5)--(6).  Owner names are not needed for this bare numerical value,
but they are needed to interpret a later LRC deletion or transport operation.

## 2. Jet response

If the slopes are fixed, `p(s)=p` and `q(s)=q`, and

```text
A(s)=sum_(r>=1) a_r s^r,
B(s)=sum_(r>=1) b_r s^r,                                 (8)
```

as formal or convergent series, then (6) gives coefficientwise

```text
[s^r](M_loc-M)=(q a_r+p b_r)/(p+q).                      (9)
```

For variable slopes, suppose the offsets vanish through order `r-1`:

```text
A(s)=a_r s^r+o(s^r),       B(s)=b_r s^r+o(s^r).          (10)
```

Continuity of the slopes and (6) give

```text
M_loc(s)-M
 =((q(0)a_r+p(0)b_r)/(p(0)+q(0)))s^r+o(s^r).             (11)
```

Thus the order-`r` coefficient is the branch-slope-weighted average, and it is
the leading coefficient if that average is nonzero.  If it cancels, slope
variation can contribute at the next order.  For example,

```text
A=s, B=-s, p=1+s, q=1
```

gives `M_loc-M=-s^2/(2+s)`: the order-one average vanishes, but the first
nonzero coefficient is not determined by nonexistent `a_2,b_2`.  This is
exact local elimination, not a heuristic Taylor fit.

## 3. The AP13 equality fibre and asymmetric wedge

For

```text
C={1,2,...,13},
f_C(t)=min_(v in C) ||v t||,                              (12)
```

one has

```text
M(C)=1/14,
G_(1/14)(C)={a/14:(a,14)=1}.                             (13)
```

Indeed, the fourteen points `0,t,...,13t` have two points at circular
distance at most `1/14`, proving the upper bound.  If every nonzero difference
has distance at least `1/14`, the fourteen cyclic gaps are each at least
`1/14` and sum to one, so all gaps equal `1/14`.  The points are therefore the
full order-fourteen subgroup and `t=a/14` with `(a,14)=1`.  Conversely such a
unit numerator permutes the nonzero residues and attains the bound.  Hence the
residual fibre consists of six isolated points.

At `t_0=1/14`, THM-2050's explicit slack estimate gives, for
`|h|<1/728`,

```text
f_C(t_0+h)=1/14+h       for h<0,
f_C(t_0+h)=1/14-13h     for h>0.                         (14)
```

The active walls are owned by speeds `1` and `13`.  In the notation of (2),

```text
p=1,                    q=13,                            (15)
```

so the local response to offsets `A,B` is

```text
M_loc-M=(13A+B)/14.                                      (16)
```

More generally, at `a/14`, let `r` be the representative in `{1,...,13}` of
`a^(-1) mod 14`.  The active owners are `r` and `14-r`; the weights are
`14-r:r`.  The equality fibre is therefore polyhedral and zero-dimensional,
not a smooth defect manifold with a transverse positive Hessian.

## 4. Exact local/global hostile

Put

```text
C={1,...,13},
L={1,...,11,13,26}.                                      (17)
```

THM-2050 proves

```text
M(C)=1/14,                    M(L)=1/12,                  (18)
```

while `f_C` and `f_L` have identical complete **unperturbed** function germs
on explicit neighborhoods of all six unit points `a/14`.  Therefore every
invariant computed only from those germs--including the active-owner and
slope data used by (5)--(11)--is blind to the global maximum.  THM-2050 does
not compare arbitrary deformation families attached to the two rows.  Any
globalization from the unperturbed carrier needs at least one off-layer
sidecar such as magnitude, first strict-exit denominator, or a proved
gluing/termination certificate.

The topology of the six-point equality fibre supplies no intrinsic edges, so
its ordinary componentwise Laplacian is zero.  A graph corrector capable of
moving response between those points would be additional structure and must
be derived from a lawful LRC operation.  It cannot be imported merely because
THM-3990 has a graph analogue.

## 5. Nonsmooth positive control

Nonsmoothness alone does not prevent a stratified higher-order certificate.
For `s>0`, `|y|<=1`, and `h` real, set

```text
Phi_s(y,h)=min(h+s^2 y^2+s^3,-h+s^2 y^2+s^3)
          =s^2 y^2+s^3-|h|.                              (19)
```

Then exactly

```text
inf_(|y|<=1) max_h Phi_s(y,h)=s^3>0.                     (20)
```

The inner maximum occurs at `h=0`, and the outer minimum at `y=0`.  Thus a
quadratic coefficient can vanish on a residual stratum and a cubic coefficient
can close it even for a polyhedral active coordinate.  What is missing in LRC
is not smoothness; it is a lawful global deformation and a complete labelled
atlas with uniform control.

## 6. Typed comparison with Brendle--Hung

In arXiv:2608.19068v1, the claimed curvature proof eliminates smooth normal
plane coordinates and obtains an effective second and third coefficient on a
compact zero-curvature locus.  Formula (6) is the precise local polyhedral
counterpart:

```text
smooth source        minimize over a nondegenerate normal Hessian
polyhedral target    maximize a labelled lower envelope over phase
preserved            local extremum and first effective response sign
required sidecars    wall owner, sign, slope, slack, and competing tops
lost in transfer     smooth tube, connected residual, legal metric corrector,
                     integer-speed realizability, and global first exit.   (21)
```

There is no map here from a metric tensor correction to a legal deformation
of an integer speed row.  The closest current global LRC responses are
THM-3910's integer `t`-sheet carrier and THM-3995's parity-sensitive support
and variance tariff for its sole scale-two survivor.  Their pairwise
statistics leave a cubic mixed indicator moment.  That use of "cubic" is
response degree in indicator variables, not third Taylor order, and no
identification of the two coefficients is claimed.

## 7. Equality and failure boundary

1. If `p(0)=0` or `q(0)=0`, the active envelope need not have an isolated
   local top; (5)--(6) are outside scope.
2. If a third wall has zero slack, it must be retained.  A two-wall formula
   need not survive a multiwall intersection.
3. Forgetting signed slope, selected side, or which offset belongs to which
   branch destroys (6).  Forgetting owner alone preserves the bare number but
   destroys LRC deletion semantics and owner transport.
4. A positive **leading** local coefficient, after all earlier coefficients
   vanish, certifies improvement only for the nearby top in that lawful
   direction.  It neither constructs the family nor compares remote tops.
5. A vanishing tested coefficient routes the audit to higher order.  A
   negative genuine leading coefficient rules out improvement of that germ in
   the chosen direction for small positive `s`; another stratum or remote top
   can still decide the global family.

Consequently THM-3993 closes no one of THM-3910's remaining 17 conditional
pair types, does not enter its `t<U` branch, and does not prove LRC(14).
