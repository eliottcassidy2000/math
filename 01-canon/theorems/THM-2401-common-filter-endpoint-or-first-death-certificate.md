---
id: THM-2401
title: "Common-filter endpoint or first-death certificate"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING. Two
  disjoint rational Boolean root packets on F_13, subjected to one
  common Boolean filtration, satisfy an exact alternative. Positive
  final two-sided mass gives a nonzero joint coefficient at all twelve
  root colours; if the two packets physicalize in one circle section,
  every colour also has a syndetic common ordinary Fourier gauge. If
  final two-sided mass is zero, the first death has an exact nonnegative
  root-deletion decomposition. Fixed singleton/two-root banks give a
  labelled failed root of mass at least one quarter of the lost
  two-sided mass. Without fixed absolute roots, a relative-displacement
  incidence still gives all twelve charged boundary colours with
  magnitude at least the lost mass divided by 2028. Typed guard/unit
  failures feed the THM-2379 repair current with exact constants; typed
  blocker failures feed a pure-owner/fork word split. The theorem does
  not construct the common terminal filtration, a canonical terminal
  current, or a proof of LRC(14).
source: codex-2026-07-26-common-filter-first-death
depends_on:
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
  - THM-2323-primitive-fixed-colour-cross-correlation-and-same-gauge-word-alignment
related:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
  - THM-2379-anchored-guard-unit-deletion-factor-repair-current
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
  - THM-2397-clean-root-same-parent-charged-role-partition
script: 04-computation/lrc14_common_filter_first_death_thm2401.py
output: 05-knowledge/results/lrc14_common_filter_first_death_thm2401.out
script_sha256: ecfee9ac955fb77cc70759de9809c574fee2500c348d10004ee4f861d3a88e27
output_sha256: 3be5cdc5f0418eb1d0b231c5f6ac9ed5b9318492db13827cae47d1ca01e30e2b
hash_basis: working-tree bytes (LF)
---

# THM-2401 -- common-filter endpoint or first-death certificate

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-2398 reduces terminal survival of a clean charged pair to one
geometric service:

```text
one common endpoint/root gauge
  + disjoint descendants
  + positive simultaneous descendant mass.                    (1)
```

This theorem makes the complementary failure branch equally explicit.
For a finite common Boolean filtration, either (1) survives to the last
stage, or one labelled filter kills co-support.  The co-support loss has
an exact positive decomposition by the failed side and failed root.  It
therefore cannot disappear through phase cancellation.

The result is deliberately conditional on the existence of the common
filtration.  Current LRC canon does not supply that filtration for the
clean-root current and the canonical terminal word.

## 1. Common Boolean descendants

Put

```text
G=F_13,

zeta=exp(2*pi*i/13).                                  (2)
```

Let `Y` be a finite union of intervals with rational endpoints.  Let

```text
A_0,F_0:Y x G->{0,1}                                 (3)
```

be Boolean step functions with rational breakpoints and suppose

```text
A_0(y,r)F_0(y,r)=0                                   (4)
```

for almost every `(y,r)`.  Let

```text
U_1,...,U_n:Y x G->{0,1}                             (5)
```

be Boolean step functions with rational breakpoints.  The word
**common** means that the same `U_j`, in the same root coordinate and
orientation, is applied to both packets.  Define

```text
A_j=A_0 product_(ell<=j) U_ell,

F_j=F_0 product_(ell<=j) U_ell,                      (6)

a_j(y)=sum_r A_j(y,r),

f_j(y)=sum_r F_j(y,r),

M_j=integral_Y a_j(y)f_j(y)dy.                       (7)
```

Assume throughout that `n>=1` and

```text
M_0>0.                                               (7a)
```

Then

```text
M_0>=M_1>=...>=M_n>=0.                               (8)
```

Equation (4) is hereditary under (6).  It is the pinned zero-shift
anchor in the positive branch below.

Use normalized root coefficients

```text
alpha_j(y,k)
 =(1/13)sum_r A_j(y,r)zeta^(-kr),

phi_j(y,k)
 =(1/13)sum_r F_j(y,r)zeta^(-kr).                    (9)
```

## 2. Positive final mass gives all twelve terminal colours

Suppose

```text
M_n>0.                                               (10)
```

Define the unnormalized cyclic cross-correlation

```text
C_n(s)
 =integral_Y sum_r A_n(y,r+s)F_n(y,r)dy.             (11)
```

Every `C_n(s)` is a nonnegative rational number.  Equations (4) and
(7) give

```text
C_n(0)=0,

sum_s C_n(s)=M_n>0.                                  (12)
```

Thus `C_n` is nonzero and nonflat.  THM-2398's prime-cyclic rational
dichotomy gives, for every `k!=0`,

```text
C_hat_n(k)!=(0),                                     (13)
```

where the transform is normalized by `1/13`.  Finite inversion gives

```text
C_hat_n(k)
 =13 integral_Y alpha_n(y,k)
                  conjugate(phi_n(y,k))dy.           (14)
```

Consequently

> **All-colour descendant conclusion.** If `M_n>0`, then
>
> ```text
> integral_Y alpha_n(y,k)conjugate(phi_n(y,k))dy!=0
> ```
>
> for every one of the twelve `k in F_13^*`.         (15)

This is a joint coefficient in one common root gauge.  It is stronger
than separate nonzero spectra of the two descendants.

There is an optional denominator-sensitive floor.  Write

```text
C_n(s)=b_s/D,

b_s in Z_(>=0),

S=sum_s b_s,

g=gcd_s b_s.                                         (16)
```

THM-2398's cyclotomic norm bound and (14) give

```text
|integral_Y alpha_n(k)conjugate(phi_n(k))|
 >=g^6/(169 D S^5)                                  (17)
```

for every `k!=0`.  There is no denominator-free positive floor.

The historical common-filter provenance in (6) is a sufficient way to
retain (4) and one gauge.  Algebraically, Section 2 needs only the final
same-gauge disjointness, rationality, and (10).

## 3. THM-2323 upgrades every colour to a common ordinary gauge

The root packets have a canonical physicalization whenever they truly
live in one circle section.  For the unique representation

```text
x=(y+r)/13,       0<=y<1,       r in F_13,           (18)
```

put

```text
mathcal A_n(x)=1_Y(y)A_n(y,r),

mathcal F_n(x)=1_Y(y)F_n(y,r).                       (19)
```

They are Boolean step functions with rational breakpoints.  Let
`J_A,J_F` be their numbers of nonzero jumps.

For fixed `k!=0`, their root amplitudes are

```text
P_A(y)=sum_r A_n(y,r)zeta^(-kr),

P_F(y)=sum_r F_n(y,r)zeta^(-kr).                     (20)
```

Equations (13)--(14) give

```text
integral_0^1 P_A(y)conjugate(P_F(y))dy
 =13 C_hat_n(k)
 !=0.                                               (21)
```

THM-2323 Sections 3--4 need exactly this nonzero fixed-colour inner
product for their landing half.  Form the two periodic gauges

```text
N_A(y)=exp(-2*pi*i*k*y/13)P_A(y),

N_F(y)=exp(-2*pi*i*k*y/13)P_F(y).                    (22)
```

Their gauge coefficients are

```text
N_A_hat(h)=13 mathcal A_n_hat(k+13h),

N_F_hat(h)=13 mathcal F_n_hat(k+13h).                (23)
```

Parseval and (21) first give some common integer `h`.  The endpoint-
product Vandermonde recurrence from THM-2323 then gives a common index
in every block of `J_AJ_F` consecutive integers:

> **Same-gauge endpoint conclusion.** For every `k!=0` and every
> integer `H`, some
>
> ```text
> H<=h<=H+J_AJ_F-1
> ```
>
> satisfies
>
> ```text
> mathcal A_n_hat(k+13h)
> mathcal F_n_hat(k+13h)!=(0).                       (24)
> ```

In particular, a common nonnegative gauge obeys

```text
0<=h<=J_AJ_F-1,

1<=k+13h<=13J_AJ_F-1.                               (25)
```

This is a direct extension of the landing half of THM-2323: its
short-arc/nesting hypotheses were used there to prove the fixed-colour
inner product nonzero, while (21) now supplies that input by THM-2398.
No nesting or common danger comb is needed after (21).

The two descendants must physicalize through the same map (18).
Separate endpoint sections do not satisfy the hypothesis.  Also, (24)
allows the selected `h` to depend on `k`; it does not give one ordinary
frequency common to all twelve colours.

## 4. Exact loss and the first co-support death

For one layer, put

```text
A_j^-=A_(j-1)(1-U_j),

F_j^-=F_(j-1)(1-U_j),

Delta_j=M_(j-1)-M_j.                                (26)
```

Pointwise in `y`,

```text
a_(j-1)f_(j-1)-a_j f_j

 =|A_j^-| f_(j-1)+a_j |F_j^-|.                      (27)
```

This is just

```text
af-a'f'=(a-a')f+a'(f-f'),
```

but its nonnegative summands retain the failed side and root.

For absolute root labels define

```text
rho^A_(j,r)
 =integral_Y A_(j-1)(y,r)(1-U_j(y,r))
              f_(j-1)(y)dy,

rho^F_(j,s)
 =integral_Y F_(j-1)(y,s)(1-U_j(y,s))
              a_j(y)dy.                             (28)
```

Integrating (27) gives the exact labelled-root partition

```text
Delta_j
 =sum_r rho^A_(j,r)+sum_s rho^F_(j,s).              (29)
```

Suppose the two initial packets use fixed root banks

```text
R_A,R_F subset G,

|R_A|=a,       |R_F|=f.                             (30)
```

Then some failed side and absolute root has mass

```text
rho>=Delta_j/(a+f).                                 (31)
```

For fixed singleton/two-root banks,

```text
a,f<=2

implies

rho>=Delta_j/4.                                     (32)
```

The weight in (28) is allowed to be integer-valued: it collapses the
number of opposite roots.  If a Boolean opposite-root label is also
required, expand instead

```text
Delta_j
 =sum_(r,s) integral_Y [
     A_(j-1)(r)(1-U_j(r))F_(j-1)(s)

    +A_j(r)F_(j-1)(s)(1-U_j(s))
   ]dy.                                             (33)
```

There are at most `2af` nonzero fixed-bank labels, so one Boolean
pair cell has mass at least

```text
Delta_j/(2af)>=Delta_j/8.                            (34)
```

Now suppose `M_n=0` and let `j_*` be the first index with `M_j=0`.
Then

```text
Delta_(j_*)=M_(j_*-1)>0.                            (35)
```

Equations (29)--(34) are the promised first-death certificate: the
filter label `j_*`, failed side, and failed root are fixed on a
positive rational cell.

The constant `1/4` in (32) is sharp for the labelled-root information.
Take fixed banks

```text
R_A={0,1},       R_F={2,3}
```

and four equal rational base cells.  On the four cells, respectively,
use the pre-filter singleton pairs

```text
({0},{2}), ({1},{2}), ({0},{2}), ({0},{3}),          (36)
```

and let the filter kill, respectively,

```text
A-root 0, A-root 1, F-root 2, F-root 3.             (37)
```

Each labelled-root summand in (29) has mass exactly `1/4`.

## 5. A universal relative-displacement charged certificate

Even when the absolute root masks vary over all thirteen positions, the
loss cannot hide its relative root displacement.  For `d!=0`, define

```text
tau_j(d)
 =sum_(s-r=d) integral_Y [
     A_(j-1)(r)(1-U_j(r))F_(j-1)(s)

    +A_j(r)F_(j-1)(s)(1-U_j(s))
   ]dy.                                             (38)
```

The zero displacement vanishes by hereditary disjointness.  Therefore

```text
sum_(d!=0)tau_j(d)=Delta_j.                          (39)
```

Some `d in F_13^*` obeys

```text
tau_j(d)>=Delta_j/12.                               (40)
```

View every incidence in (38) on its own disjoint finite cover of the
base.  A singleton at root `r` has normalized coefficient
`zeta^(-kr)/13`; a singleton at root `s` has conjugate coefficient
`zeta^(ks)/13`.  Thus the derived charged incidence in displacement
class `d=s-r` is exactly

```text
Z_(j,d)(k)
 =tau_j(d) zeta^(kd)/169.                           (41)
```

For every one of the twelve `k!=0`,

```text
|Z_(j,d)(k)|
 =tau_j(d)/169
 >=Delta_j/2028.                                    (42)
```

This all-colour conclusion is independent of rationality.  Its
noncancellation comes from freezing the relative displacement before
summing.

The carrier in (41) is a coefficient-derived pair-incidence current on
a finite cover.  It does not retain the absolute failed root, and it is
not one of the final descendants `A_j,F_j`.  Splitting which side failed
costs a further factor two and gives `Delta_j/4056`.

The denominator `12` in (40) is sharp without fixed absolute root
banks: on twelve equal cells take `A={0}`, `F={d}` for the twelve
`d!=0`, and delete the `A` root.  Every displacement has mass `1/12`.
For fixed banks, `12` may be replaced by the number of realized
differences, at most `af<=4`.

## 6. Quantitative selection before the first death

If `M_n=0`, telescoping gives

```text
M_0=sum_(j=1)^n Delta_j.                             (43)
```

Hence some layer, not necessarily the first-death layer, satisfies

```text
Delta_j>=M_0/n.                                     (44)
```

At that layer:

```text
fixed small root banks:
  one labelled failed root has mass >=M_0/(4n);

arbitrary absolute roots:
  one charged displacement has every colour
  of magnitude >=M_0/(2028n).                       (45)
```

The distinction between (35) and (44) is load-bearing.  The first death
has a canonical stopping time but may retain arbitrarily little of the
initial co-support.  The telescoping layer has a quantitative initial-
mass floor but need not be the terminal death.

## 7. Typed factor-repair and blocker-word routes

The certificate (28) becomes a THM-2379 input only with additional
physical typing.

Fix a labelled-root summand of mass `rho`.  Reparameterize the selected
root by `x=(y+r)/13`, including the Jacobian in the collapsed weight.
Suppose:

1. the failed filter is the safe indicator `u_L(vx-beta)` of an
   ordinary unit (`L=1`) or the guard (`L=2`);
2. failure therefore supplies the danger bit `d_L(vx-beta)=1`; and
3. the deepest safe factor `d_1(cx-alpha)=0` is retained at the same
   physical point.

The collapsed opposite-root multiplicity in (28) is a nonnegative
integrable weight, which is allowed in THM-2379 Sections 1--5.  That
theorem gives a positive deep-by-repair coefficient

```text
ordinary unit:
  Re Khat^+(a,b)>=11rho/24336;

guard:
  Re Khat^+(a,b)>=rho/2704,                          (46)
```

with both displayed probe multipliers coprime to `91`.  Under the fixed
small-bank floor (32),

```text
ordinary unit:
  Re Khat^+(a,b)>=11 Delta_j/97344;

guard:
  Re Khat^+(a,b)>=Delta_j/10816.                    (47)
```

If one complete four-way nondeep-blocker status is required, the
THM-2379 partition costs four in amplitude:

```text
ordinary unit:
  Re Khat^+_sigma(a,b)>=11 Delta_j/389376;

guard:
  Re Khat^+_sigma(a,b)>=Delta_j/43264.              (48)
```

Equations (47)--(48) are coefficient-derived repair currents.  If a
literal Boolean opposite-root cell is required instead of the collapsed
weight, use (34), which can cost an additional factor two.

There is a separate blocker-word route.  Suppose the failed filter is
the safe indicator of a named target blocker `b`, and the certificate
already lies in THM-2305's prescribed-expiration gauge with its source
owner and clock fixed.  To retain physical set mass rather than the
integer-valued opposite-root multiplicity in (28), use a Boolean pair
cell from (34).  The failed root is dangerous for `b`.  Splitting that
cell by the other target-blocker bit `a` gives exactly

```text
{b}        pure owner transfer,

{a,b}      fork hyperedge.                          (49)
```

One word has mass at least

```text
Delta_j/16.                                         (50)
```

Here the factor `8` first selects the Boolean side/root/opposite-root
cell and the factor `2` fixes the other blocker bit.  If the opposite
packet is already a singleton, the sharper root weight may be Boolean
and the cost can improve.  If the source owner was not already fixed,
an additional three-way source partition gives the coarser
`Delta_j/48`.

Equation (49) never guarantees the pure alternative: all mass may lie
in the fork.  Without the prescribed expiration, source owner, root
gauge, and clock, it is only a local blocker-status word and cannot be
cited as a THM-2305 canonical handoff.

The remaining first-death types are honest residuals:

- deletion of the deepest safe factor violates THM-2379 hypothesis 3;
- an arbitrary Boolean root selector need not be a physical factor;
- a composite unlabelled filter has no unique failed role unless it is
  refined into a lawful atomic filtration; and
- the repair current still lacks a lawful target quotient and terminal
  component phase.

## 8. Exact partition costs

The finite costs used above are:

| retained label | number of cells | guaranteed mass |
|---|---:|---:|
| failed side and root, fixed banks `a,f` | `a+f` | `Delta/(a+f)` |
| failed side/root and opposite root | `2af` | `Delta/(2af)` |
| relative displacement only | `12` | `Delta/12` |
| side plus relative displacement | `24` | `Delta/24` |
| arbitrary absolute failed root at `p=13` | `26` | `Delta/26` |
| one additional blocker bit | multiply by `2` | divide by `2` |
| two nondeep blocker bits | multiply by `4` | divide by `4` |
| unfixed source owner | multiply by `3` | divide by `3` |

If a two-root role first ranges over thirteen cyclic translates, fixing
its full translate costs `13`.  A subsequent fixed-bank failed-root
selection costs at most `4`, for total cost `52`; retaining an opposite
root costs at most `104`.  If the full translate is unnecessary, the
direct absolute-root partition can be cheaper.

## 9. Sharp hostiles and exact boundaries

Four small hostiles prevent stronger conclusions.

### 9.1 Separate masses are not simultaneous mass

Split a rational base in halves.  Put an `A` singleton only on the first
half and an `F` singleton only on the second.  Both separate packets are
positive, but

```text
M=0,       C(s)=0 for every s.                      (51)
```

Thus separate terminal survival does not trigger Section 2.

### 9.2 There need not be a common endpoint gauge

On thirteen equal rational strata, put one singleton in each relative
endpoint displacement.  In separately chosen endpoint sections the
apparent correlation is

```text
C(s)=1/13            for every s.                   (52)
```

Every charged mode vanishes.  This does not contradict Section 2:
there is no common pinned anchor `C(0)=0`.  It is the minimal flat
mean-projection hostile from THM-2398.

### 9.3 First-death mass has no floor from initial mass

Take `A={0}`, `F={1}` and split the base into masses `1-epsilon` and
`epsilon`, with rational `epsilon>0`.  Let the first common filter kill
`A` only on the large stratum, and let the second common filter kill
`A` on the remaining stratum.  Then

```text
M_0=1,       M_1=epsilon,       M_2=0.               (53)
```

The first-death mass is `epsilon`, arbitrarily small.  Therefore (35)
cannot be replaced by a positive multiple of `M_0`; use the earlier
telescoping layer (44) when an initial-mass bound is required.

### 9.4 A labelled root deletion is not automatically a factor repair

On one fibre take `A={0}`, `F={1}`, and let the common Boolean filter
delete root `0`.  This satisfies the first-death algebra, but the filter
may be an arbitrary root selector and no deepest-safe factor has been
specified.  Neither THM-2379 nor THM-2305 applies.

These examples also separate three notions which must not be merged:

```text
same-gauge final descendants,

charged first-death incidence,

canonical physical terminal current.                (54)
```

## 10. LRC interface and stopping ledger

Applied after THM-2398, the positive route is:

```text
clean charged parent pair
  + one common rational Boolean terminal filtration
  + positive final co-support
  -> all twelve terminal root colours
  -> a common ordinary Fourier gauge for each colour.           (55)
```

The zero route is:

```text
zero final co-support
  -> exact first-death filter
  -> labelled failed root or charged displacement
  -> typed THM-2379 repair / typed THM-2305 word,
     if the extra physical hypotheses hold.                     (56)
```

This theorem does not prove that the clean-root current and canonical
terminal word are acted on by one common filter.  THM-2370 retains a
common right endpoint only before the fully masked endpoint introduces
cross-layer terms; THM-2380 likewise records common-endpoint service as
open.  Nor does a coefficient-derived displacement incidence in (41)
become a canonical owner current without an endpoint/owner sidecar.

Accordingly, no scalar profile is excluded.  The ledger remains `165`,
and LRC(14) remains open.

## 11. Exact companion

The dependency-free exact companion:

- exhausts all `5,760` common Boolean filter layers on the `180`
  nonempty disjoint two-packet assignments on `F_5`, checking the root,
  pair, and displacement decompositions with explicit raising tests;
- checks all `1,573` disjoint singleton/adjacent-by-two-root packet
  pairs used in the `p=13` all-colour interface;
- constructs the sharp four-cell labelled-root example and sharp
  twelve-displacement example;
- checks every constant in (42), (47), and (48);
- checks rational first-death controls through denominator `101`; and
- records the separate-support and flat-gauge hostiles.

Run

```bash
python3 04-computation/lrc14_common_filter_first_death_thm2401.py
python3 -O 04-computation/lrc14_common_filter_first_death_thm2401.py
```

Both transcripts must byte-match, after LF normalization,

```text
05-knowledge/results/lrc14_common_filter_first_death_thm2401.out
```

All truth-bearing checks raise explicitly, so optimized mode executes
the same audit.
