# Desheeting gives a genuine one-base U_full endpoint coupling

**Status: FINITE-EXACT CANDIDATE; independent audit pending.**  The natural
desheeting of the 39 U_full guard atoms produces one common endpoint base,
retains the seven THM-2594 cells, and has complete `7 x 13` mixed spectrum.
It is an alternate diagonal coupling of the actual endpoint integrands, not
the original Cartesian endpoint current.  No THM-2471 collision/horizon
record, source packet, chronology, grouped exact-address coefficient, row
exclusion, or LRC(14) theorem follows.

## Inheritance pass and the missed diagonal

THM-3514 already writes the guard coordinate as

```text
91t=7a+r,                 a in F_13, 0<=r<7,
```

and proves that the endpoint sheet label `a` and THM-2471 root label transform
by the same regular `F_13` action.  The canonical hostiles then tested equality
of endpoint components and equality of the two circle points.  Both miss the
frozen Cartesian bridge.  The corrected near miss is MISTAKE-293: one common
ancestry base is not one circle point.  The least-used sidecar is the residual
coordinate `r`, because it supplies the base rather than another sheet label.

Indeed, set

```text
y=r/7=13t-a,            t=(y+a)/13.                  (1)
```

Equation (1) is exactly THM-2471's owner-node chart

```text
w_a=(y+a)/13.                                             (2)
```

Thus two endpoint points with sheets `a,b` lie over the same proposed base
when

```text
13t_L-a=13t_R-b=y.                                      (3)
```

They are equal as circle points only when `a=b`.  The old point hostile is
therefore the same-sheet sector of (3), not the full common-base relation.

## Exact descent of the word and harmonic

The two endpoint factors have frequencies

```text
X=13,
Y=X+w_c=742599,
w_c=742586=13^2*4394.                                  (4)
```

On the left leg the delayed word is evaluated at `169t`.  Equation (1) gives

```text
169t_a=13y+13a =13y mod 1,                              (5)
```

so the word becomes the same sheet-independent factor `Q(13y)` on every
branch.  The endpoint phase also descends without a sheet residue:

```text
-X t_a+Y t_b
 =(w_c/13)y-Xa/13+Yb/13
 =57122 y + integer.                                   (6)
```

After putting `x=13y`, the frequency is `4394`.  The branch phase is one
because `w_c/13^2=4394` is integral.  These divisibilities are load-bearing:
without them (3) would carry an additional sheet phase sidecar.

For a character `(alpha,beta,tau)`, let

```text
N_(alpha,beta,tau)(y)
 =sum_(a in F_13)
   1_(E_(alpha,beta,tau))((y+a)/13).                   (7)
```

The common-residual endpoint numerator is the one-base integral represented
by

```text
Q(13y) N(y)^2 exp(2 pi i 57122 y).                     (8)
```

Its same-sheet and cross-sheet sectors are respectively

```text
Q(13y) N(y) exp(...),
Q(13y) N(y)(N(y)-1) exp(...).                          (9)
```

The seven half-open THM-2594 cells are inserted on `y` before endpoint
summation.  The computation processes `7,107,008` unguarded intervals,
`7,108,460` atom-split mapped intervals, `1,199,656` active residual
segments, and `186,244` word-active segments.

## Independent controls internal to the candidate

Five character triples are rebuilt through the literal fully guarded E set,
without using the unguarded-plus-safe restoration.  Their seven-cell full and
same-sheet rows agree entrywise.

More decisively, summing the seven cells in the same-sheet sector reproduces
the entire previously audited point-diagonal bank:

```text
sha256=771545a5cb1f0f03459b8d351de668ad950ece5fcb985fa61d599d643de3303f.
                                                               (10)
```

Its inverse values and bridge are exactly

```text
q_H:  633668780131603861,
q_q5: 405160484437854840264,
bridge:167726070588785644466.                           (11)
```

This recovers a disjoint implementation with the correct branch phases,
endpoint convention, guard restoration, and character order.  The full table
then splits exactly as same-sheet plus cross-sheet.

## The cross-sheet sector is nonzero and spectrally complete

The new cross-sheet inverse values are

```text
q_H:  95149748277639510265,
q_q5: 291433430760196195223,
bridge:375969203763952195911.                          (12)
```

Adding (11) and (12) modulo the certified prime gives the full
common-residual coupling

```text
q_H:  95783417057771114126,
q_q5: 124341028951542154618,
bridge:543695274352737840377.                          (13)
```

The value in (13) is not the frozen Cartesian bridge
`389266878372286537904`.  This is expected: a common-base diagonal and a
product of two separately integrated marginals are different currents.  The
result therefore constructs a natural alternate coupling; it does not
retroactively identify the Cartesian endpoint bank with ancestry.

For each of the full, same-sheet, and cross-sheet `7 x 13` cell-by-refined-
residue tables, the split-field Fourier support is

```text
(total,DC,F_7 axis,F_13 axis,mixed)=(91,1,6,12,72).     (14)
```

Doubly centering the full output kills the axes and leaves every one of the
`72/72` mixed modes nonzero.  Since each displayed nonzero reduction has a
unit denominator in a certified cyclotomic embedding, all corresponding
characteristic-zero coefficients are nonzero.

The fixed all-unit relation class `(1,0,6)` reduces to

```text
(289814661037836286866,0,0,0,0,0,0)                   (15)
```

in the seven-cell chart, so all seven of its septimal Fourier values reduce
to the same nonzero number.  The six modular zeros in (15) are reported only
as a feature of this split-field image; one prime does not prove those six
characteristic-zero cell coordinates vanish.  The seven nonzero Fourier
reductions do prove characteristic-zero nonvanishing of all seven modes.

## Connection contract and exact boundary

| field | exact content |
|---|---|
| source geometry | actual U_full E intervals and delayed Q intervals |
| vertices | thirteen endpoint sheets over a residual base `y` |
| map | `t -> (a,y)` from (1), independent of the character |
| pair relation | equal `y`, allowing all ordered sheet pairs `(a,b)` |
| preserved | guard restoration, E/Q interval geometry, seven cells, refined residue transform |
| same-sheet hostile | exactly the old common-point diagonal |
| destroyed/absent | THM-2471 collision root, source and word sheets, horizons, chronological arrival, exact relation address orbit |
| next sidecar | lift (1) through a THM-2471 collision record and test compatibility with the source-aligned `M(omega,nu,ell,c)` tensor |

This closes the cheapest interval-level/Fubini experiment from the preceding
source-aligned reflection: the endpoint factors can be opened into an actual
common-base integrand while retaining the seven cells.  It does not close the
larger temporal theorem, because the resulting base has not been shown to be
the collision base of the canonical LRC row, and the diagonal current differs
from the original Cartesian current.

In tournament language the square `N(y)^2` retains all ordered sheet pairs;
the point hostile kept only loops, and `N(y)(N(y)-1)` is the loop-deleted
directed complete fibre.  That decomposition is intrinsic here, unlike a
cosmetic K4: the observable is literally pair multiplicity over a shared
base.  No orientation, absolute `H^1` flux, or physical current is inferred.

## Reproduction

```text
python -B 04-computation/lrc_ufull_desheeted_common_residual_base_probe_20260816.py
python -B -O 04-computation/lrc_ufull_desheeted_common_residual_base_probe_20260816.py
```

The full table, full spectrum, and semantic SHA-256 values are

```text
64adfc3850b8a8924546a75e051116aa10e3a62661ca409175aad174a9ebb926
398a2c2ccb57714aa9d735d3bca4f24c2212396a8be9e2bd07e8e4ec57b6c05a
97694f50b071dbf875802223dd894d7e1df86ef65a194c22d588505c59eea507.
```
