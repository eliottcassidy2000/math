# The half and `C_221` fibres contain strict deepest-owner cospans, but the present wall deletes every attachment

**Status: FINITE-EXACT stopping/re-root certificate on the canonical typed
row.**  This is a new full-fibre audit of the THM-2698 half carrier and the
transverse `C_221` carrier.  It is not a row exclusion and does not prove
LRC(14).

## 1. The semantic map and the affine maps are different arrows

Use the canonical row

```text
(H,q1,...,q5,c1,c2,c3)
 =(1,14,27,40,53,66,13,13^3,2*13^5).
```

The canonical exclusive sources and THM-2305 terminal words are

```text
E_j=A0 intersect D_(c_j) minus union_(i!=j) D_(c_i),

Q_(j,sigma)
 =A0 intersect intersection_(i in sigma) D_(c_i)
     minus union_(i notin sigma union {j}) D_(c_i)
     minus D_(c_j).
```

Their prescribed clocks are `k_1=2`, `k_2=4`, and `k_3=6`.  In this audit
the semantic chronology is always the ordinary dilation

```text
x in E_j  ->  D^(k_j)x in Q_(j,sigma),
D(x)={13x}.                                               (1)
```

The half arrow `x->{13x+(2k+1)/(2*13^6)}` and the transverse arrow
`x->{13x+m/(17*13^6)}` are separate affine arrows.  No affine phase is
inserted into `(1)`.  This type distinction is load-bearing.

## 2. The scalar `E_3`-to-fork cospan is genuinely positive

The old displayed cycle centres fail as exclusive sources, but that failure
is not uniform on their delayed fibres.  The half fibre contains the strict
point

```text
x_half=54102827/115843416,
D^6 x_half=11/24 in Q_(3,{1,2}).                          (2)
```

The complete source-and-terminal constraint ledger stays strict on the
symmetric physical-`x` interval of radius

```text
1/301082946198216.                                        (3)
```

Likewise the full transverse fibre contains

```text
x_221=13074908/82055753,
D^6 x_221=4/17 in Q_(3,{1,2}),                            (4)
```

with strict symmetric radius

```text
11/853068347561612.                                       (5)
```

Thus both intersections

```text
E_3 intersect D^(-6)Q_(3,{1,2})                          (6)
```

are nonempty open sets in the tested fibres.  The scalar THM-2305 cospan is
positive.  What fails is attachment to the richer local packet cycles.

## 3. Full half-fibre census: the present factor is the first wall

For each half phase `11/24` and `13/24`, the source-one predecessor window
has the same exact unique-centre census:

| stage | `E_1` | `E_2` | `E_3` |
|---|---:|---:|---:|
| raw fibre window | 38,256 | 11,406 | 12,848 |
| nonzero shallow clock | 0 | 11,406 | 12,848 |
| inherited present packet | 0 | positive | 0 |

The `E_1` branch dies at the shallow-clock gate.  Every one of the `12,848`
`E_3` points survives the nonzero shallow clock.  More strongly, all `12,848`
also retain an active rail and a strict private root if the present factor is
temporarily omitted.  Nevertheless:

```text
E_3 intersect union_(r in F_13) Present_(shallow,r) = empty             (7)
```

on both phases.  This is not a wrong inherited present label: none of the
thirteen available present labels contains an `E_3` point.

After the present gate, the full-packet universe contains only source-free
points and `E_2` points.  Its exact candidate multiplicities are

```text
phase 11/24: source-free 139314, E2 17484,
phase 13/24: source-free 154648, E2 21222.                (8)
```

The `E_2` survivors do not repair the semantic leg: their ordinary `D^4`
terminal word is `None`.  The complete packet universe has `417971564`
reciprocal root-glued cycles, but their semantic source pairs are only

```text
(0,0): 312418764,
(0,2):  98768480,
(2,2):   6784320.                                        (9)
```

There is no `E_3` pair even before a root-pair condition is tested, because
the event vertices themselves have already been removed by `(7)`.

## 4. Full `C_221` census: root re-selection reaches the same present wall

On each complete delayed fibre `4/17` and `13/17`, the unique-centre source
census is

```text
E1=313680,       E2=119089,       E3=131752.              (10)
```

The published literal-root subgenerator fixes carry `3` over `4/17` and
carry `9` over `13/17` in order to retain right root `6`.  Those two carry
slices already contain no `E_3` point.  Their final packet candidates again
have only source-free or `E_2` labels, and their `902300` reciprocal pairs
have source pairs

```text
(0,0): 668922,
(0,2): 218306,
(2,2):  15072.                                           (11)
```

The carry/root choice is re-rootable, so the scout also exhausts every
`E_3` carry class before stopping.  They are exactly

```text
phase 4/17:  carry 0 (65653),  carry 6 (66099),
phase 13/17: carry 6 (66099), carry 12 (65653).            (12)
```

The private-root formula then leaves only the extremal choices

```text
carry class one:  left root 1,  right root 0,
carry class two:  left root 0,  right root 12.             (13)
```

Thus roots `1` and `12`, not root `6`, are the complete nonzero re-root
alternative.  After retaining intrinsic nonconstant clocks and delayed
sector zero, the exact survivor counts are

```text
phase 4/17:  53227 and 53590,
phase 13/17: 53590 and 53227.                              (14)
```

Every survivor in `(14)` fails the dynamically typed present factor.  Hence
there are zero present, rail, unit, or reciprocal-pair survivors in the
complete variable-root `E_3` search.  Exact `C_221` microphase transport was
retained in the pair predicate, but it never gets a vertex on which to act.

## 5. What is proved, what is auxiliary, and what remains

The two layers now separate cleanly:

```text
THM-2305 semantic layer:
  E3 + ordinary D^6 + deepest fork Q_(3,{1,2})
  -> strict positive intervals exist;

local-cycle sidecars:
  rail + present + clock + carry + private root + primitive unit
  -> no E3 vertex survives the present factor;

affine half/C221 arrows:
  positive reciprocal cycles exist on other source labels
  -/-> attachment to the positive E3 cospan.               (15)
```

Rail, present, clock, carry, root, and unit data are not hypotheses of
THM-2305.  They are auxiliary coordinates introduced to make a local packet
transition physical and composable.  Therefore the positive scalar cospans
`(2)--(6)` are real, but they do not by themselves inherit the half or
`C_221` cycle.

The strongest survivor is especially precise: the half `E_3` fibres already
have the needed shallow clock, rail, and private root, while the all-carry
`C_221` fibres reach nonconstant clock plus delayed sector zero.  Both stop
at the same present wall.  The next honest construction should therefore be
a relative-present or present-free mapping cone which attaches the semantic
`E_3` current before rebuilding the unit sidecar.  Another phase or private-
root search in the unchanged present category is saturated by `(7)` and
`(14)`.

No endpoint current, uniform row reduction, scalar-row exclusion, or
LRC(14) conclusion follows.

## Reproduction

```bash
python3 04-computation/lrc14_half_c221_semantic_source_fibre_census_20260728.py
python3 -O 04-computation/lrc14_half_c221_semantic_source_fibre_census_20260728.py
```

Both modes byte-match

```text
05-knowledge/results/lrc14_half_c221_semantic_source_fibre_census_20260728.out
```

SHA-256:

```text
script  0fbeb041007fea1b9e14f0ff6e82fc97ebf724b26c2c10ef85732b4c994b94cd
output  f32268a8c255cd5b04e1bd28aef5307e87dbf9a19936ccdde705371bab1e8320
```

The companion uses exact `Fraction` and integer modular arithmetic, explicit
optimized-mode guards, both the inherited candidate generators and an
all-carry extremal-root hostile control.  Normal and optimized transcripts
are identical.
