# The fixed-root hostile supplies rank four, but not the `K4` connection

**Status: FINITE-EXACT REPRESENTATION THEOREM CANDIDATE; independent audit
pending.**  The exact companion is
`04-computation/lrc_r5_absolute_root_fourth_channel_probe_20260816.py`.
It reuses THM-2594's one-common-base Boolean table and the proved THM-3514/3515
U_full endpoint bank.  It constructs no U_full ancestry relation, physical
current, absolute `H^1`, grouped coefficient, scalar-row exclusion, or
LRC(14) conclusion.

## 1. Why the fixed-root table was the right hostile to revisit

The folded-`C_7` transporter probe found a sharp dimension obstruction.  The
theta-slaved THM-2594 response has only the three nonempty windows
`theta=0,1,2`, so its seven septimal character rows span a three-dimensional
space.  The U_full chamber Walsh bank spans four dimensions.  No common
drift convolution and fixed channel map can raise rank three to rank four.

THM-2594 already contains a second lawful marginal of the *same* joint table:

```text
A_abs(ell,t)=sum_(u,q)N(u,q,ell,t-2u),                  (1)
```

where `t` is held as an absolute deep-root label rather than the relative
slaved coordinate `theta=t-2u`.  This was introduced as the hostile proving
that affine slaving is not the unique source of nonvanishing.  Because (1)
is formed from the same `N(u,q,ell,theta)` before leaving the common-base
construction, it is the cheapest lawful candidate for the missing fourth
source direction.

## 2. Two three-planes meet in a two-plane

Let `S_b(theta)` be the septimal Fourier rows of the slaved table and
`A_b(t)` those of (1), for `b in F_7`.  Over the same certified split field
used by THM-3514,

```text
rank{S_b:b in F_7}=3,
rank{A_b:b in F_7}=3,
rank({S_b} union {A_b})=4.                              (2)
```

Therefore

```text
dim(span(S) intersect span(A))=3+3-4=2.                 (3)
```

The absolute-root reindex really does supply one new linear coordinate.  It
is neither a copy of the slaved response nor an unrelated rank-seven bank.
The two descriptions are two three-planes sharing a two-plane inside one
four-dimensional response space.

Both sparse physical tables are nevertheless cyclic in every septimal
channel:

```text
slaved raw cells nonzero:        12/91,
absolute-root raw cells nonzero: 18/91,

C7 x C13 Fourier entries nonzero:
  slaved 91/91,                  absolute 91/91.         (4)
```

Thus entrywise spectral support, row rank, and compatibility with a target
connection are three different invariants.

## 3. The exact projective connection equation

Let `X_a in F_P^m` be the source channel vector at drift frequency
`a in F_13`, and let `Y_a in F_P^4` be the four U_full Walsh values at that
frequency.  A fixed channel map `M:F_P^m->F_P^4` followed by one common
circulant drift operator exists only if

```text
M X_a = mu_a Y_a                 for every a,           (5)
```

for scalars `mu_a`.  Eliminating the thirteen scalars gives the homogeneous
wedge equations

```text
Y_(a,i)(MX_a)_j-Y_(a,j)(MX_a)_i=0,
0<=i<j<=3.                                             (6)
```

There are `6*13=78` equations in the `4m` entries of `M`.  Maps annihilating
the whole source span always solve (6), so raw nullity alone is misleading.
The decisive statistic is

```text
excess = nullity(equations)-4(m-rank(source)).          (7)
```

Positive excess is necessary for a nonzero projective transport; excess
zero says every formal solution kills every source frequency.

## 4. Rank four is necessary and still insufficient

The exact systems are

```text
bank          source rank   equation rank   nullity   excess
-------------------------------------------------------------
slaved 7           3             12            16       0
absolute 7         3             12            16       0
union 14           4             16            40       0.   (8)
```

For the union, the forty-dimensional nullspace is exactly the space of
`4 x 14` maps annihilating the four-dimensional source span.  Hence adjoining
the absolute-root table removes the elementary rank obstruction but still
does not produce one shared projective drift connection to the U_full Walsh
curve.

This is stronger than checking one preferred calibration.  It allows an
arbitrary frequency-independent linear combination of all fourteen slaved
and fixed-root character rows.

## 5. Every marked fourth-sidecar test fails

The probe also keeps the folded-`C_7` marking explicit.  Choose one member
of each antipodal class

```text
{1,6}, {2,5}, {3,4}                                    (9)
```

for the three slaved nontrivial channels, and append one named absolute-root
character.  All

```text
2^3*7=56                                                (10)
```

resulting four-row banks have rank four.  Every projective system has excess
zero.

As a second hostile, retain slaved modes `1,2,3` and use every nonempty
binary sum of the seven absolute rows.  Of the `127` choices, `126` have
rank four and all `127` have excess zero.  The sole rank-three mask is `127`,
the sum of all seven absolute characters.  Fourier inversion concentrates
that sum on the `ell=0` row, which vanishes by THM-2594's word-guard
disjointness.  This explains the exceptional mask rather than merely
recording it.

## 6. The full affine torsor gauge does not repair the mismatch

An identification of two regular `F_13` torsors is not canonically tied to
one generator.  Besides translation, it may apply any dilation

```text
t |-> at+c,                a in F_13^*.                 (11)
```

Translation contributes one common Fourier phase and is absorbed by the
circulant multiplier in (5).  The probe separately exhausts all twelve
nonzero dilations, including reflection.  For every one, the union system
still has

```text
(source rank,equation rank,nullity,excess)=(4,16,40,0). (12)
```

Thus the no-go is not an arbitrary choice of drift origin or generator.

## 7. What the fourth coordinate means

The fixed-root hostile answers the previous probe's question with an
important split verdict:

```text
Does a lawful fourth source direction exist?       yes;
Does it make the source and endpoint connections
projectively equivalent through one convolution?  no.                 (13)
```

The remaining debt is therefore not dimension but connection.  This is the
same distinction visible in the Berggren/Fibonacci `T4` atlas: four states
or four independent channels do not determine how those states move around
a loop or drift torsor.  A `K4` carrier is `H^0`-sized data; equation (5)
tests a connection.

The formal channel-dependent circulants from the preceding transporter
remain available, but allowing four unrelated multipliers makes any four
cyclic source/target pairs equivalent and carries no ancestry content.  The
useful target must retain a common operator or a more geometric typed map.

## 8. The next lawful object

The next LRC test should not search for another marginal rank.  It should
keep the coordinate that the two marginalizations destroy:

```text
(u,q,ell,theta) with t=theta+2u,                         (14)
```

and compare it to an ell-independent U_full endpoint address before either
the relative or absolute root is summed out.  Equivalently, the desired map
may be:

- a nonconvolutional address relation on the joint ancestry table;
- a frequency-dependent channel connection with an independently justified
  geometric law;
- a nonlinear Boolean relation before linear Fourier contraction; or
- a typed cocycle whose holonomy explains the failure of (5).

The exact obstruction says that no frequency-independent linear channel map,
even using every slaved and absolute character and every affine root gauge,
can be that relation.

## 9. Connection and loss ledger

| field | exact content |
|---|---|
| source | THM-2594 joint common-base table, viewed through relative and absolute roots |
| target | THM-3514/3515 U_full four-row Walsh drift bank |
| positive invariant | two source three-planes have union rank four and full `91/91` spectra |
| failed map | fixed channel map plus one common `C_13` circulant |
| exhaustive gauges | all fourteen rows, 56 named sidecars, 127 binary sidecars, 12 torsor dilations |
| exact obstruction | every projective nullspace is annihilator-only |
| preserved source typing | common Boolean base for both THM-2594 marginalizations |
| still absent | U_full common ancestry, physical owner/word/current, chronology, absolute flux |

## 10. Reproduction

Run

```text
python -B 04-computation/lrc_r5_absolute_root_fourth_channel_probe_20260816.py
python -B -O 04-computation/lrc_r5_absolute_root_fourth_channel_probe_20260816.py
```

The pinned semantic digest is

```text
6ad55eeab7f8e5479194965e819b030f177f4500be37cb189688262dc8af35da.
```

The primary THM-2594 construction is rerun from its hash-pinned source; the
U_full target is independently reconstructed from the hash-pinned all-role
atom engine.  The result is a finite exact representation obstruction, not a
physical-current or LRC theorem.
