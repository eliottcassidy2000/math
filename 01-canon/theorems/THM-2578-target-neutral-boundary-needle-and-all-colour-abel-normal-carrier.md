---
id: THM-2578
title: "Target-neutral boundary needle and all-colour Abel-normal carrier"
status: >
  PROVED + VERIFIED-EXACT.  For every positive integer tooth speed k and
  L in {1,2}, the 26k boundary points of the thirteen target translates
  d_L(kx+s/13), labelled by tooth, sign, and target shift, are pairwise
  distinct.  The obstruction is exactly the coprimality of 7 and 13: a
  cross-sign collision would make plus or minus 13L vanish modulo 7.
  Consequently every Boolean trace pattern on this complete boundary orbit
  is realized by one fixed rational target-neutral finite union of intervals
  H.  In particular, H can isolate one chosen boundary x0 belonging to one
  chosen target shift s0.  The lawful disjoint layers
  H d_L(kx+s/13) and H(1-d_L(kx+s/13)) then have THM-2573 handoff measure
  delta_x0 at s=s0 and zero at the other twelve shifts.  At every physical
  offset M their one-sided logarithmic Abel normal has all thirteen target
  characters nonzero, with exact coefficient
  exp(2 pi i Mx0) zeta^(q s0)/(26 pi^2).  This remains nonzero at every
  live M=m c_3 after multiplying by the deepest coefficient because
  gcd(m,91)=1.  Unlike THM-2574's complex component connection, the carrier
  is a positive Boolean aggregate with no relative lift gauge.  It is still
  external: the THM-2569 stationary packet is supported inside the old
  danger gate, so used as a common filter it kills the safe layer rather
  than crossing the boundary.  An inherited target-neutral transverse
  filter retaining word, owner, future root, and deep mode remains open.
  No THM-2334 relation current, row exclusion, or LRC(14) conclusion follows.
source: common-endpoint-seam-2026-07-28-boundary-needle
depends_on:
  - THM-2573-logarithmic-abel-normal-and-common-endpoint-jump-pairing
related:
  - THM-2568-full-x-transition-annihilation-and-refined-pair-drift-boundary
  - THM-2569-stationary-diagonal-conditioned-paired-corner-and-frozen-future-role-boundary
  - THM-2574-oriented-tooth-component-holonomy-and-fixed-frequency-descent
  - THM-2579-socle-flat-target-torsor-and-integral-difference-filling
script: 04-computation/lrc14_boundary_needle_all_colour_thm2578.py
output: 05-knowledge/results/lrc14_boundary_needle_all_colour_thm2578.out
script_sha256: 0f9187a36d8661e46e175cd63c551b84c93b567cb3c45c4347f73d6a45b5f2f3
output_sha256: 67bc42246cf789017a11d6260e0631b8cf57bdb98ec313cb01e8d95597fb2efb
hash_basis: working-tree bytes (LF)
---

# THM-2578 -- one target-neutral boundary needle sees every colour

**PROVED + VERIFIED-EXACT.**

THM-2574 identifies why the bare `k_a` danger/safe normal is target-trivial
in the live bank: equal physical tooth weights cancel off resonance, while
resonance lands at colour zero.  That cancellation is not an intrinsic
limitation of positive Abel normals.  It is exactly the uniform-weight
hostile.

The target translates of a guard/unit boundary form a finite Kakeya-like
needle atlas with a stronger property: no two labelled needles meet.  A
single fixed, target-neutral Boolean filter can therefore prescribe the
trace independently at every boundary in the atlas.  Isolating one boundary
turns the target profile into a delta function and forces **all** thirteen
target characters.

```text
7/13 boundary separation
  -> arbitrary Boolean boundary traces by one fixed rational filter
  -> one-boundary target delta
  -> every C_13 character nonzero.                            (1)
```

This closes the abstract positive-carrier problem for the Abel normal.  The
remaining issue is inheritance: the live common-ancestry packet does not
currently contain such a transverse target-neutral filter.

## 1. The complete target-tooth boundary atlas is separated

Put `p=13`.  For `L in {1,2}`, let

```text
d_L(y)=1_(||y||<L/14).                                      (2)
```

Fix a positive integer `k`.  The two boundaries of the `j`th physical tooth
under target shift `s` are

```text
x_(j,epsilon,s)
 =(j+epsilon L/14-s/13)/k mod 1,                             (3)

j in {0,...,k-1},       epsilon in {+1,-1},
s in F_13.
```

All `26k` points in (3) are pairwise distinct.

### Proof

Suppose two labels give the same circle point.  For some integer `n`,

```text
j-j'-kn +(epsilon-epsilon')L/14 -(s-s')/13=0.                (4)
```

If `epsilon=epsilon'`, multiplication by `13` shows that `13|(s-s')`.
Since both target representatives lie in `{0,...,12}`, one has `s=s'`.
Then `j-j'=kn`; the chosen tooth representatives force `j=j'`.

If the signs differ, divide the equation obtained from multiplying (4) by
`182` by two.  It would imply

```text
7(s-s') plus-or-minus 13L is in 91 Z.                        (5)
```

Modulo seven, (5) says

```text
plus-or-minus 13L=0 mod 7.                                  (6)
```

But `13=6 mod 7` and `L` is `1` or `2`, so (6) is impossible.  This proves
the separation. QED.

Every point in (3) is rational.  Because the atlas is finite and separated,
its minimum circular gap is a positive rational number.

## 2. One fixed Boolean filter realizes every boundary trace pattern

Let

```text
b_(j,epsilon,s) in {0,1}                                    (7)
```

be an arbitrary assignment on the atlas (3).  Choose pairwise disjoint
rational open arcs around exactly those boundary points with `b=1`, each
arc shorter than one third of the minimum atlas gap, and let

```text
H:T->{0,1}                                                  (8)
```

be the indicator of their union.  Then:

- `H` is a rational step function independent of `s`, hence target-neutral;
- `H` is identically `b_(j,epsilon,s)` on a two-sided neighbourhood of each
  corresponding boundary point; and
- no endpoint of `H` is a boundary point in (3).

For

```text
P_s(x)=d_L(kx+s/13),

A_s=H P_s,                         C_s=H(1-P_s),              (9)
```

the two layers are nonnegative and disjoint.  At a selected tooth boundary,
`H` has two-sided trace one, so

```text
-Delta A_s Delta C_s=1.                                    (10)
```

At an unselected tooth boundary, `H` vanishes on a neighbourhood and both
jumps are zero.  At an endpoint of `H`, the gate `P_s` is locally constant;
exactly one of `A_s,C_s` can jump, so the common-jump product is again zero.

Therefore THM-2573's total-layer handoff measure is exactly

```text
nu_s
 =sum_(j,epsilon)b_(j,epsilon,s)delta_(x_(j,epsilon,s)).     (11)
```

Equation (11) proves a full Boolean trace-interpolation theorem:

```text
fixed rational target-neutral filters
  -> every {0,1}-valued weight profile on the complete target/tooth atlas.
                                                                    (12)
```

No complex component character, connection, or relative lift is used.

## 3. A single boundary needle forces all target characters

Choose one label `(j_0,epsilon_0,s_0)` and use (7) with value one there and
zero everywhere else.  Put

```text
x_0=x_(j_0,epsilon_0,s_0).                                  (13)
```

Then (11) becomes

```text
nu_s=delta_(x_0),                    s=s_0,

nu_s=0,                              s!=s_0.                 (14)
```

For every physical offset `M in Z`, THM-2573 gives the one-sided logarithmic
Abel normal

```text
N_s(M)
 =exp(2 pi i Mx_0)/(2 pi^2),          s=s_0,

 =0,                                  s!=s_0.                (15)
```

With `zeta=exp(2 pi i/13)` and the plus-sign target transform,

```text
Nhat(q;M)
 =1/13 sum_s N_s(M)zeta^(qs)

 =exp(2 pi i Mx_0)zeta^(q s_0)/(26 pi^2).                  (16)
```

The right side is nonzero for **every** `q in F_13`, including all twelve
nonzero target colours.  The statement is uniform in `M`; no tooth resonance
is needed because the physical filter has destroyed the equal-component
cancellation before aggregation.

For a live deepest leg,

```text
M=m c_3,                         gcd(m,91)=1.                 (17)
```

The separate deepest danger coefficient is nonzero because `7` does not
divide `m`.  Multiplication by this nonzero scalar preserves all thirteen
coefficients in (16).  This does not by itself attach the remaining LRC word
and owner factors, but it proves that the physical-offset restriction creates
no further analytic cancellation for the needle carrier.

## 4. Lawfulness and the distinction from the oriented component lift

The construction (8)--(9) is a lawful common-target endpoint family:

```text
H                         fixed target-neutral physical filter;

P_s,1-P_s                 the same gate co-shifted on both endpoint copies;

whole Boolean products    formed before Fourier/Abel smoothing.       (18)
```

Thus (16) is an honest target character of a positive whole-layer Abel
normal.  It does not repeat MISTAKE-266's frozen-moving-projector error.

THM-2574 takes a different route.  It Fourier-weights the `k` tooth labels,
retains one complex component character at fixed `X`, and chooses a
connection which trivializes its monodromy.  That reference has a relative
integer-lift gauge.  The present theorem instead changes the **positive
physical boundary weights** with one actual Boolean filter and then takes
the aggregate normal.  The resulting target delta (14) has no component
gauge and needs no fixed-`X` selection.

The two results are complementary:

```text
THM-2574: orient a cancelling component algebraically;

THM-2578: isolate one component physically and positively.            (19)
```

The survivor in (16) is still the singular high-frequency Abel boundary
coefficient of THM-2573, not the ordinary full-`X` sum and not automatically
a THM-2334 factorwise relation-address current.

### Absolute reference, not a relative difference

THM-2579 proves that target differences and complete target Fourier
contractions cannot retain the absolute `13`-primary Cayley class: their
target coefficient sum is zero modulo `13`.  The needle has exactly the
opposite type.  The physical interval `H` is chosen first, independently of
the target variable, and its unique incidence with the boundary atlas picks
one absolute label `s_0` before any target DFT is taken.  Its target profile
is

```text
const * delta_(s=s_0),               target augmentation = const !=0,  (20)
```

not a difference of two profiles.  Changing target character after this
choice only multiplies the already-selected reference by
`zeta^(q s_0)`.  Thus the `q=0` coefficient and every `q!=0` coefficient
survive together.  In particular, the needle is an explicit **ambient
absolute charged reference** for the boundary-normal target torsor, rather
than another relative atlas.

This identifies the correct kind of sidecar requested by THM-2579, but does
not yet transport its integral Cayley class.  No map has been proved from
this ambient interval into THM-2571's carry lattice, and no inherited packet
has been shown to generate the interval while preserving its word, owner,
root, paired-dipole, and deep-mode data.  Consequently the needle can serve
as the missing reference only after that packet-level inheritance and
comparison map are constructed.

## 5. Why the current stationary packet is not the filter

On the THM-2569 common packet, the old target-informed weight satisfies

```text
w_(N,h) P_0=w_(N,h),                                       (21)
```

where `P_0` is the old `k_a` danger gate.  If this one-sided packet is used
as the common filter in (9), then

```text
w_(N,h)(1-P_0)=0.                                          (22)
```

It therefore kills the repaired layer at the base target shift instead of
having positive two-sided trace across the gate boundary.  The total-layer
jump convention of THM-2573 assigns no handoff merely because `P_0` and
`1-P_0` meet underneath a filter that dies on one side.

This is an exact obstruction, not a reason to freeze or discard the packet.
The positive needle theorem changes the remaining obligation to:

> Construct, from inherited live sidecars, a target-neutral rational Boolean
> filter with positive two-sided trace at at least one `k_a` boundary and a
> covariant whole-layer orbit, while retaining the THM-2569 word, owner,
> stationary future root, paired dipoles, and an allowed deepest mode.

The external `H` in Section 2 proves that such a filter exists in the ambient
endpoint algebra.  It does not prove that the scalar-cover row supplies it,
or that arbitrary multiplication by `H` preserves positive packet mass.

## 6. Exact consequence for the live seam

The updated boundary ledger is

```text
complete danger/safe endpoint + ordinary full-X sum
  -> zero;                                                   THM-2568

same endpoint + logarithmic Abel normal
  -> total-layer handoff measure;                            THM-2573

bare equal-weight tooth boundary in the live bank
  -> component cancellation or target-zero resonance;       THM-2574

one fixed target-neutral boundary needle
  -> lawful positive normal with every target colour.        THM-2578 (23)
```

Thus neither Fourier cancellation, physical-frequency resonance, nor target
character support is the remaining analytic problem.  The missing object is
a **transverse inherited filter** on the common semantic packet.  No row is
excluded; the ledger remains `165`, and LRC(14) remains open.

## 7. Exact companion

Run

```bash
python3 04-computation/lrc14_boundary_needle_all_colour_thm2578.py
python3 -O 04-computation/lrc14_boundary_needle_all_colour_thm2578.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_boundary_needle_all_colour_thm2578.out.
```

The dependency-free exact referee checks:

- `1,045,200` distinct labelled boundary points in the complete atlas for
  `k<=200`, together with the symbolic cross-sign congruence obstruction;
- `262,600` exact boundary traces in `200` rational single-needle controls
  and `400` filter-endpoint exclusions;
- `2,080` all-character target coefficients and `12,800` live `91`-unit
  physical phases; and
- the exact three-cell truth table `H<=P => H(1-P)=0` for the one-sided
  stationary-packet hostile.

The all-`k` separation, arbitrary trace interpolation, handoff identity, and
finite-DFT nonvanishing are symbolic proofs above, not finite extrapolations.
**QED.**
