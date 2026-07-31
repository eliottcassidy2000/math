# A polarized E3/complement frequency seam

Status: `FINITE-EXACT SCRATCH RESULT`; not canon and no LRC(14) row
exclusion.

The 20-cell THM-2847 horn has all ordinary factors at q3, q7, and q11;
q7 differs only by lying in the complementary E3 truth block.  Restricting
the same endpoint pattern to its Boolean complement gives one physical
endpoint interval at q7 with the same two Prony nodes as the q3 and q11 E3
intervals.

Use the three exceptional THM-2868 section labels

```text
(r,q,truth,local offset)
 = (0,3,E3,0), (4,7,not-E3,1), (8,11,E3,2).
```

All six raw multipliers are 91-units.  In both certified endpoint fields,
the actual adjacent samples split and transport to

```text
U_(r+4)=omega^12 U_r,       V_(r+4)=V_r.
```

The labels lie on `q-r=3`.  Therefore the unique THM-2868 invariant channel
`(3,10)` gives

```text
3r+10q=4
```

at all three seam vertices, and

```text
U_r omega^(10q)
```

is literally constant.  The trivial branch `V_r omega^(10q)` is not.
Thus the character-three branch, rather than the unsplit endpoint current,
is the exact coefficient that can be polarized consistently across the
q3/q7/q11 E3/complement seam.

This sharpens THM-2868's numerical seam observation: the two unit-chart
repairs and the physical q-triangle are the three points of one diagonal
character fibre.  It still does not contract the truth blocks.  The q7
chart is physical only in `not-E3`, while q3 and q11 are physical in `E3`;
combining them into one current would still be a type error.  The remaining
object is now precise: a lawful polarized E3/complement morphism preserving
the constant `(3,10)` seam reference.

## The Bockstein hostile: the apparent carry phase is flat

There is an exact tempting coincidence on THM-2851's carry leg.  If

```text
rho(h)=omega^(3h),
```

then `rho(9)=omega`, the same value as the first ancestry-carry character
`chi(T)`.  The whole triangle rejects the identification:

```text
rho(9)rho(8)/rho(4)=1,
chi(T)=omega != 1.
```

Indeed `rho` is an ordinary character of the residue group.  Pulling the
Kummer line back by `r=q-3` makes it depend only on q, so it kills the
normalized carry 2-cocycle.  The target character 10 cancels it precisely
because `10=-3`.  Thus the constant `(3,10)` seam is a flat Cech
coefficient shadow, not the Bockstein transgression.

This is a useful negative result.  The isolated q11-to-q7 phase match must
not be called ancestry carry.  A carry-faithful seam needs dependence on the
ancestry coordinate—equivalently the additional `C13` fibre/nonsplit
`C169` object already invoiced by THM-2851—not another reparametrization of
the frequency index.

Reproduce:

```text
python3 .scratch/lrc_e3_complement_frequency_seam_20260728/probe.py
python3 -O .scratch/lrc_e3_complement_frequency_seam_20260728/probe.py
```
