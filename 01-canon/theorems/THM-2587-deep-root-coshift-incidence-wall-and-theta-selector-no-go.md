---
id: THM-2587
title: "Deep-root/translated-gate incidence wall and theta danger-selector no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On the canonical typed sigma={b}, K=2, depth-five packet, an absolute
  deepest-speed digit is not a THM-2365 deep-comb or endpoint co-shift
  label.  Writing its phase as (z+r)/13, the translated danger gate with
  label tau is active exactly at tau=r for z<13/14 and at tau=r+1 for
  z>1/14.  The middle wall has both incidences.  Each theta half of the
  THM-2584 packet has positive low, middle, and high mass, so no
  diagonal-translation-equivariant selector tau=phi(r,theta) can select an
  incident danger translate everywhere.  Retaining the three-state wall gives
  an exact four-edge bipartite toothpick.  The projected danger and safe translate
  tensors have all 169 global, all 1,183 owner-cell, and all 1,183
  owner-centred joint characters nonzero.  This exhausts the deep-coordinate
  projection, but a lawful target co-shift is the paired dipole
  e_(c_3)-e_(k_b); the compensating graft and left endpoint residue are
  absent.  No THM-2365/2334 current, row exclusion, or LRC(14) conclusion
  follows.
source: endpoint-root-comparison-2026-07-28
depends_on:
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2584-b-word-depth-five-absolute-deep-root-tensor
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2563-paired-dipole-deep-target-corner-and-partial-bare-boundary
  - THM-2569-stationary-diagonal-conditioned-paired-corner-and-frozen-future-role-boundary
  - THM-2586-depth-five-arrival-to-future-root-diagonal
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
script: 04-computation/lrc14_deep_root_coshift_incidence_wall_thm2587.py
output: 05-knowledge/results/lrc14_deep_root_coshift_incidence_wall_thm2587.out
script_sha256: 63a83ac103cc50ec8c5bd3cdfe29a188ed918704adcc150f77b2c0e9743f59a4
output_sha256: b65340f6ca063c11ef86b52e721020b250daa201e6f1c89818beb07df4f81779
hash_basis: working-tree bytes (LF)
---

# THM-2587 -- the absolute root sees a two-neighbour translated-gate wall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2584 retains an absolute deepest-speed root on the same positive
`sigma={b}`, depth-five ancestry as the owner-clock host.  THM-2365 uses two
different labels from the same thirteen-element alphabet: a translate of the
deep danger comb and a translate of the safe endpoint.  The labels cannot be
identified with the absolute digit.  The correct comparison is a bipartite
incidence relation:

```text
absolute deep digit -- within-cell wall --> translated danger/safe label.
                                                                    (1)
```

The relation is very small and spectrally complete.  Its failure to be a
function is also sharp: THM-2584's binary rail `theta` does not retain the
within-cell wall needed to choose one incident translate.

## 1. The pointwise incidence law

Put `p=13` and

```text
d_1(q)=1_(||q||<1/14).
```

For a physical deepest-speed phase write uniquely away from cell endpoints

```text
2x=(z+r)/13 modulo 1,

0<=z<1,                         r=floor(26x) mod 13.       (2)
```

For a translated deep-coordinate label `tau in F_13`, define

```text
I_z(r,tau)=d_1(2x-tau/13).                                (3)
```

Then, up to the null boundary points,

```text
I_z(r,tau)=1
 iff [tau=r   and z<13/14]
  or [tau=r+1 and z>1/14].                                (4)
```

Indeed put `q=tau-r` in the representatives `0,...,12`.  For `q=0`,
condition (3) is `z/13<1/14`, hence `z<13/14`.  For `q=1`, it is
`(1-z)/13<1/14`, hence `z>1/14`.  Every other residue stays at circular
distance at least `1/13>1/14`.

Thus the within-cell wall has three states:

| wall state | range of `z` | incident translate labels |
|---|---:|---|
| low | `0<z<1/14` | `{r}` |
| middle | `1/14<z<13/14` | `{r,r+1}` |
| high | `13/14<z<1` | `{r+1}` |

The diagonal translation

```text
(r,tau)->(r+a,tau+a)                                    (5)
```

preserves (4).  This is an intrinsic bipartite relation, not a tournament:
the two vertex classes are absolute root cells and translated deep gates,
the observable is danger incidence, and the middle state has degree two.
Orienting or deleting one middle edge would add a gauge choice.

The same translated family appears twice in THM-2365.  Its `r` indexes the
danger probe `Delta_r`, while its `t` indexes the complementary endpoint
factor `1-d_1(cx-t/13)`.  Equation (4) therefore compares the THM-2584 digit
with either **projected coordinate**, but identifies it with neither label.

## 2. The binary theta rail cannot choose an incident danger translate

In the THM-2584 root chart,

```text
x=(y+v)/13,

theta=floor(2y) in {0,1},       z={2y},

r=2v+theta mod 13.                                      (6)
```

The exact canonical packet has, for each fixed `theta`, positive mass in
all three wall states.  The exact `(s,r,ell)` positive-cell census is

```text
                         low    middle    high

theta=0                   48      154       48
theta=1                   48      154       48.             (7)
```

In both low fibres the only absolute root is `r=0`; in both high fibres it
is `r=-1`; the middle fibres contain both roots.  Every displayed total is
strictly positive before marginalization.

Suppose a selector

```text
tau=phi(r,theta)                                          (8)
```

were compatible with (4) and equivariant under (5).  Equivariance forces

```text
phi(r,theta)-r=c_theta                                   (9)
```

for a constant `c_theta`.  Positive low mass forces `c_theta=0`, whereas
positive high mass forces `c_theta=1`.  Equation (7) gives both constraints
for each `theta`, a contradiction.

This is the sharp selector no-go.  The pure-low packet is a positive control
for `tau=r`; the pure-high packet is a positive control for `tau=r+1`.
Retaining the three-state wall recovers the complete **relation** (4).  If a
single incident edge must be retained losslessly, the middle fibre has two
elements, so one overlap-edge bit there is necessary and sufficient.  A
convention can choose one edge, but the packet supplies no canonical physical
choice between them.

The quantifier "incident" is essential.  Every selector equivariant under
(5) has the form `tau=r+c_theta`.  On a fixed `theta` rail it is safe for
every `z` exactly when

```text
c_theta notin {0,1}.
```

Thus there are eleven uniformly safe offsets per rail, or `11^2` if the two
rails are chosen independently.  These are undistinguished gauge choices,
not a packet-defined target label, and they still translate only the `c_3`
gate.  The no-go above concerns selecting one of the **active danger** labels
in (4), as would be required to identify the absolute digit with the live
translated danger probe.  It does not claim that projected safe translates
do not exist.

## 3. The exact four-edge incidence toothpick

Retain THM-2584 notation

```text
U=P_(13^5)(1_Q P_(13^2)1_E),

V=P_(13^5)1_E,

s=arrival root-source root,                ell in F_7.     (10)
```

Let `H_ell` be the native owner-clock cells pulled to the physical `x`
coordinate as in THM-2584 equation (12).  Define the fine positive tensor

```text
K_ell(s,r,tau)
 =13 integral_T U(x)V(x-s/13)
     1_(floor(26x)=r)
     1_(frac(169x) in win_ell)
     d_1(2x-tau/13) dx.                                  (11)
```

Its absolute-root/translate support is exactly

```text
(r,tau)=(0,0),(0,1),(-1,0),(-1,-1).                       (12)
```

Every one of the four edges is positive in every owner cell.  As a graph it
is the four-edge path

```text
tau=1 -- r=0 -- tau=0 -- r=-1 -- tau=-1.                  (13)
```

This is a second toothpick descendant of THM-2584's arrival/deep-root path.
It comes from the overlap geometry of the `1/7` danger arcs against the
`1/13` digit cells, not from equality of root and target labels.

Define its danger-translate marginal and safe complement by

```text
D_ell(s,tau)=sum_r K_ell(s,r,tau),

S_ell(s,tau)=C^2581_ell(s)-D_ell(s,tau).                   (14)
```

Both are nonnegative rational tensors on one common physical integral.  At
the original translate,

```text
D_ell(s,0)=C^2581_ell(s),          S_ell(s,0)=0.           (15)
```

The exact support counts over `(s,ell)` are

```text
tau                 -1       0       1       all other

D positive           77      84      77          0,

S positive           84       0      84         84.       (16)
```

In particular, the projected safe endpoint is positive in every nonzero
translate, every nonzero collision displacement, and every owner cell.  The
old translate is zero because the `sigma={b}` arrival packet is literally
deep-dangerous there.

## 4. Both projected translate tensors are spectrally saturated

For `k,h in F_13`, put

```text
Dhat_ell(k,h)
 =1/169 sum_(s,tau)D_ell(s,tau)zeta_13^(-ks-ht),

Shat_ell(k,h)
 =1/169 sum_(s,tau)S_ell(s,tau)zeta_13^(-ks-ht).           (17)
```

Let `Dhat,Shat` be their sums over `ell`, and centre each owner cell by
subtracting one seventh of the global value.  Exact reduction modulo
`Phi_13` gives

```text
Dhat(k,h)!=0,              Shat(k,h)!=0
                    for all 169 pairs;

Dhat_ell(k,h)!=0,          Shat_ell(k,h)!=0
                    for all 1,183 triples;

Dhat^c_ell(k,h)!=0,        Shat^c_ell(k,h)!=0
                    for all 1,183 triples.                 (18)
```

The pointwise mechanism is visible in (4): the co-shift character multiplier
is `1`, `1+zeta_13^h`, or `zeta_13^h` on the low, middle, or high wall.
None vanishes because thirteen is odd.  This alone does not preclude
cross-wall cancellation after integration; the exact census (18) proves that
no such cancellation occurs on the canonical packet.  A second finite-field
certificate covers all `5,070` danger/safe global, cell, and centred buckets:
`4,974` are nonzero at the primitive thirteenth root `18 mod 79`, and the
remaining `96` are nonzero at `107 mod 131`.

Thus the missing map is not another Fourier-support theorem for the deepest
coordinate.  That projection is already saturated in both danger and safe
forms.

## 5. Exact loss ledger and remaining target lift

The connection contract is:

| item | exact content |
|---|---|
| source | THM-2584's `(s,r,theta,ell)` depth-five ancestry |
| target | danger/safe translates of the `c_3` coordinate |
| map | the incidence relation (4), refined by the wall state |
| preserved | one common positive Boolean ancestry and literal deep-gate membership |
| root-only loss | `z` relative to `1/14,13/14`; `theta` does not restore it |
| permanent projection loss | the graft coordinate `k_b` and the left endpoint relation residue |
| cheapest next test | insert the oppositely translated `k_b` factor before marginalization on this same fine wall tensor |

A lawful second target direction is

```text
ell=e_(c_3)-e_(k_b),                                      (19)
```

not the raw coordinate `e_(c_3)`.  THM-2563's moving endpoint therefore
translates the `c_3` and `k_b` factors in opposite directions.  Equations
(11)--(18) move only the first factor.  No function of `(r,theta,z)` can
recover the truth value or Fourier index of the separate graft coordinate.

THM-2334 has a further left/right distinction: its target character sees
`ell.(u-v)`, whereas (11) is a positive spatial incidence tensor with no
left endpoint harmonic `u`.  Consequently (18) is neither a target quotient
character nor a fixed-`X` current.  It cannot be inserted into THM-2365's
`H(r,s,t)` or THM-2334's target variance without first adding the paired
graft and endpoint sidecars on the same integral.

The theorem concerns one canonical typed non-cover row.  It proves no
uniform statement on the other `164` rows, excludes no scalar row, and gives
no LRC(14) conclusion.  The scalar ledger remains `165`.

## 6. Exact companion

Run

```bash
python3 04-computation/lrc14_deep_root_coshift_incidence_wall_thm2587.py
python3 -O 04-computation/lrc14_deep_root_coshift_incidence_wall_thm2587.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_deep_root_coshift_incidence_wall_thm2587.out.
```

The first route splits the common root-chart integral into the six
`theta x {low,middle,high}` wall pieces.  The independent route integrates
the physical formula (11) directly against all translated speed-two danger
combs.  They agree on all `15,379` fine cells.  The companion then checks
(7), (12), (15)--(18), the two selector controls, exact rational
normalization, both modular certificates, and the complete loss ledger.

The independent `common-endpoint-seam` hostile audit rederived (4), checked
the sign convention in the danger intervals, reran the exact companion, and
recovered all `15,379` cells, the `48/154/48` census, and the four-edge
support.  It also checked that the safe tensor is literally the same-integral
complement and audited the projected-versus-lawful scope against THM-2365 and
THM-2334.  A third physical-`x` reconstruction, independent of the root-chart
wall split, reproduced the six rational wall masses and checked (4) at all
`123,032` fine-cell midpoints.  No defect was found.

**QED.**
