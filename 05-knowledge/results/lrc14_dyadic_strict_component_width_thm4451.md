# Sharp strict component widths for odd clock-two failure

**Status: PROVED, computer-assisted all-height + VERIFIED-EXACT +
INDEPENDENTLY AUDITED; LRC(14) OPEN.** This report also corrects an
endpoint-topology error in the prior odd-3-unit candidate value `1/49`.

## Theorems

For a positive integer `n`, let

```text
D_n={x in R/Z: ||nx||<1/14},  tau(x)=x+1/2.
```

For three distinct positive odd tails `T={a,b,c}`, define the strict physical
two-sheet failure set

```text
A_T=D_a union D_b union D_c,
F_T=A_T intersect tau(A_T).
```

Write `lambda(T)` for the supremum of the lengths of the actual connected
components of the open set `F_T`.

> **All-odd sharp component theorem.**  For every distinct positive odd
> triple,
>
> ```text
> lambda(T) <= 17/693.
> ```
>
> Equality holds uniquely at `T=(1,9,11)`.

> **Odd-3-unit sharp strict component theorem.**  If additionally no member
> of `T` is divisible by three, then
>
> ```text
> lambda(T) <= 19/1001.
> ```
>
> Equality holds uniquely at `T=(1,11,13)`.

These are statements in physical `x`-space.  Under the quotient coordinate
`y=2x`, antipodal physical components are identified and every component width
is doubled.  The corresponding quotient bounds are therefore `34/693` and
`38/1001`, respectively.  There is no hidden factor-two normalization.

Unlike Haar mass, component length is not invariant under common dilation:
an odd dilation pulls each component back to shorter lifts.  The equality
statements are individual triples, not dilation rays.

## Endpoint correction: why `1/49` is not a connected-component value

For `(1,7,13)`, an a.e.-merged interval routine reports a physical interval
from `3/49` to `4/49`, of length `1/49`.  But its midpoint is
`x=1/14`, where

```text
||x||=1/14, ||7x||=1/2, ||13x||=1/14.
```

No tail is strictly dangerous on the first sheet, so `x` is absent from
`F_(1,7,13)`.  The apparent interval is two open components of length `1/98`,
not one component.  This triple's actual longest components have length
`1/91`; its exact strict width multiset is

```text
4 copies of 1/91 and 8 copies of 1/98.
```

The height-99 exact audit finds 1,540 primitive all-odd rows whose strict
topology differs from the a.e.-merged topology.  Filling isolated endpoint
holes is legitimate for Haar measure and for the a.e.-minimal BV variation in
the exposure proof, but it is invalid for connected-component localization.

The former `1/49 at (1,7,13)` observation is thus **REFUTED as a strict
connected-component equality**.  The all-odd `17/693 at (1,9,11)` observation
survives unchanged because its constituent pair intervals overlap by positive
length.

## 1. Two-lane interval-capacity lemma

For an odd tail `n`, put

```text
C_n=D_n union tau(D_n).
```

The two sets in this union are disjoint.  The set `C_n` is a periodic comb of
period `1/(2n)`, with one open tooth of length `1/(7n)` per period.

Let `J` be an interval of length `L`.  Scaling by `2n`, write

```text
2nL=m+r,  m=floor(2nL), 0<=r<1.
```

The largest possible amount of the unit-period tooth `(-1/7,1/7)` in an
interval of real length `m+r` is `2m/7+min(r,2/7)`.  Therefore

```text
cap_n(L):=sup_J mu(J intersect C_n)
         =[2m/7+min(r,2/7)]/(2n).                    (1)
```

If `J subset F_T`, each point of `J` has a first-sheet owner and a second-sheet
owner.  An odd tail cannot own both sheets.  Hence, pointwise on `J`, at least
two of the six lane indicators are active, and integration gives

```text
2L <= sum_(n in T) mu(J intersect C_n)
    <= sum_(n in T) cap_n(L).                         (2)
```

Define the endpoint surplus

```text
s_n(L)=cap_n(L)-2L/7.
```

From (1),

```text
s_n(L)=[min(r,2/7)-2r/7]/(2n),
0 <= s_n(L) <= 5/(49n).                              (3)
```

The maximum in (3) occurs at `r=2/7`.  Combining (2)--(3), an interval of
length `L` inside a failure component requires

```text
sum_(n in T) s_n(L) >= 8L/7.                         (4)
```

If a component were longer than the claimed bound, it would contain an
interval of exactly that length, so (4) is a valid obstruction to a strict
violation.  Inclusive versions of the bounds below also retain every possible
equality case.

## 2. Finite capacity reduction

Sort the actual tails as `a<b<c` and put `Delta=8L/7`.  Equation (3) first
gives

```text
a <= 15/(49 Delta).                                  (5)
```

After computing the exact `s_a(L)`, it then gives

```text
b <= 10/[49(Delta-s_a)].                             (6)
```

If `s_a+s_b<Delta`, it finally gives

```text
c <= 5/[49(Delta-s_a-s_b)].                          (7)
```

For `L=17/693`, (5) is `a<=1485/136<11`; the exact possible smallest tails
are `1,3,5,7,9`.  Equations (6)--(7) leave only these bounded families with a
possible third tail:

| smallest pair | inclusive `c` bound |
|---|---:|
| `1:9`, `3:9`, `5:9` | `495/8` |
| `1:11`, `3:11`, `5:11` | `45/2` |
| `1:13`, `3:13`, `5:13` | `6435/412` |
| `7:9` | `495/28` |

Every other base pair admits no allowed `c>b` below the inclusive bound in
(7). Usually that bound is at most `b`; the only extra numerical case is
`(7,11)`, whose bound `165/14` is still below the next odd label `13`.

For `L=19/1001` in the odd-3-unit domain, (5) is
`a<=2145/152<15`; the possible values are `1,5,7,11,13`.  The nonempty
bounded families are

| smallest pair | inclusive `c` bound |
|---|---:|
| `1:11`, `5:11`, `7:11` | `715/4` |
| `1:13`, `5:13`, `7:13` | `715/18` |
| `1:17`, `5:17`, `7:17` | `12155/614` |

This leaves a finite exact calculation except when the first two surpluses
already reach `Delta`.

## 3. Geometric closure of the nominally unbounded pairs

Let

```text
P_ab=(D_a union D_b) intersect tau(D_a union D_b).
```

Fill only isolated deleted endpoints of `P_ab` and merge the resulting a.e.
components; call this explicit superset `P*_ab`.  Let `w` be its largest
component width and `g` its least positive circular gap.  Split the full
failure set into pair-owner sets.  Exact odd ownership gives

```text
F_(a,b,c)=P_ab union Sigma_ac union Sigma_bc,
Sigma_ac union Sigma_bc subset C_c.                  (8)
```

Every component of `C_c` has width `ell=1/(7c)`. If `P*_ab` is empty,
(8) gives `F_(a,b,c) subset C_c`, so every component has width at most
`ell`. Otherwise, if `ell<g`, a `C_c` component cannot meet two
components of `P*_ab`. The incidence graph of components in
`P*_ab union C_c` is then a disjoint union of stars, so every component
has width at most

```text
max(ell, w+2ell).                                    (9)
```

For nonempty `P*_ab` this is a superset bound, hence it applies to
`F_(a,b,c)`; the empty-carrier case uses the sharper `ell` bound above.

For the all-odd target, the capacity-unbounded smallest pairs and their exact
geometry are:

| pair | `w` | `g` | safe for every odd `c>=` | finitely checked `c` |
|---|---:|---:|---:|---|
| `1:3` | 0 | 1 | 7 | 5 |
| `1:5` | 0 | 1 | 7 | none |
| `1:7` | `1/98` | `6/49` | 21 | 9,11,13,15,17,19 |
| `3:5` | `1/210` | `5/42` | 15 | 7,9,11,13 |
| `3:7` | `1/98` | `19/98` | 21 | 9,11,13,15,17,19 |
| `5:7` | `1/98` | `1/14` | 21 | 9,11,13,15,17,19 |

For the odd-3-unit target, only `1:5`, `1:7`, and `5:7` are
capacity-unbounded.  Their safe thresholds are respectively 11, 35, and 35;
the audit checks all smaller allowed `c`.  At each displayed safe threshold,
the right side of (9) is **strictly** below the target, so equality cannot hide
in an infinite family.

This geometric step is the missing owner/address sidecar that scalar capacity
discards: `C_c` remembers which narrow third-tail tooth can connect components
of the fixed two-tail carrier.

## 4. Exact residual boxes and equality topology

The capacity and geometric reductions leave exactly:

```text
all odd:          123 distinct triples;
odd 3-units:      209 distinct triples.
```

Every row is evaluated twice with exact rational arithmetic:

1. union the six oriented owner-cross interval families, merging positive
   overlaps but retaining a shared deleted endpoint;
2. partition at every original wall, evaluate the defining strict predicate
   at cell midpoints and shared walls, and join cells only through a live wall.

The implementations agree on every row.  The resulting leaders are:

```text
all odd:
  max 17/693 at (1,9,11), unique;
  runner-up 1/42 at (1,7,9).

odd 3-units:
  max 19/1001 at (1,11,13), unique;
  runner-up 29/1547 at (1,13,17).
```

For `(1,9,11)`, there are eight strict physical components: four of width
`17/693` and four of width `13/1386`.  Its longest component is the positive-
overlap union

```text
(3/77,4/77) from pair 1:11,
(1/21,4/63) from pair 1:9,
```

namely `(3/77,4/63)`, of width `17/693`.

For `(1,11,13)`, there are twelve strict physical components: four each of
widths `19/1001`, `17/2002`, and `3/2002`.  The natural grid is
`2002=14*11*13`; the proof scripts retain the complete exact interval lists.

## Validity and scope checks

- **Physical versus quotient:** all widths above are physical.  Quotient
  widths are exactly twice these values; Haar masses, by contrast, are
  unchanged by doubling.
- **Strict versus essential topology:** `merge_strict` joins only positive
  overlaps.  The separate Haar/exposure code intentionally merges touching
  intervals a.e.; the two conventions are not interchangeable.
- **No height extrapolation:** the height-99 scan discovered the endpoint
  issue but is not used for the theorem.  Equations (1)--(9) give an infinite
  reduction, and only the explicitly bounded 123/209 rows are computed.
- **Parity:** oddness is essential.  For even `n`, `D_n=tau(D_n)`, so the
  distinct-owner inequality (2) fails.
- **No cosmetic tournament:** oriented owner pairs form an intrinsic directed
  relation (sheet-zero owner to sheet-one owner), but the proof uses their
  interval union directly.  A tournament completion would discard ties and
  endpoint truth without preserving component connectivity.

## File and reproduction

Run:

```powershell
python -B 04-computation/lrc14_dyadic_strict_component_width_thm4451.py
python -O -B 04-computation/lrc14_dyadic_strict_component_width_thm4451.py
python -B 04-computation/lrc14_dyadic_strict_component_width_thm4451_independent.py
python -O -B 04-computation/lrc14_dyadic_strict_component_width_thm4451_independent.py
```

The two runs are line-identical to the frozen output and end in `PASS`.
The height-99 discovery scan is provenance, not a theorem dependency.

Raw-LF SHA-256 values:

```text
04-computation/lrc14_dyadic_strict_component_width_thm4451.py
  c8e1a17b87fd0c7b71c4dcb378db42f6ad94fd9a419d5c33bb8832cbb2a4c238
04-computation/lrc14_dyadic_strict_component_topology_thm4451.py
  e6f3dd3933076c05d452211ab6d4d4ef4601746e0735a19e9e5564c48c98a951
05-knowledge/results/lrc14_dyadic_strict_component_width_thm4451.out
  e80aef3abc01c484fd81b5e72d8767c9ef6d063ccc6a7f8e3fcf165579749a19
04-computation/lrc14_dyadic_strict_component_width_thm4451_independent.py
  7c8ebcbc41beee193c422042a7f09b170c7ecb231325628208b7efb02aa2b92d
05-knowledge/results/lrc14_dyadic_strict_component_width_thm4451_independent.out
  772f293b16edefb6750128327bdb7919df0f287dc3c82d5aa3e9e8883771fa9d
```

The clean-room report is
`05-knowledge/results/lrc14_dyadic_strict_component_width_thm4451_independent_audit.md`.

## Next tasks

1. Replace any canonical/component-localization use of a.e.-merged widths by
   strict wall-aware widths; audit especially the quoted `(1,7,13)` control.
2. Feed the physical `17/693` and `19/1001` bounds into the ten-body entry:
   a connected compact body-safe component of either threshold cannot fit in
   the corresponding strict tail component (equality is also excluded by
   compact-in-open topology).
3. Seek a closed-form classification of equality in the capacity inequality.
   The present finite reduction is small, but its residue word may yield a
   paper proof of the 123/209 boxes.
4. Apply the same lane-capacity lemma to the q=4 one-v2 seam after its absorbed
   odd tail is added; retain sheet color rather than collapsing to Haar mass.
