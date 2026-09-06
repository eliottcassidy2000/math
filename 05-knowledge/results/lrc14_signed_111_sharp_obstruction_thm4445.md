# Signed (1,1,1) sharp obstruction classification

Status: **PROVED ELEMENTARY RELATIVE TO THM-4414 + FINITE-EXACT +
INDEPENDENTLY AUDITED (THM-4445); LRC(14) OPEN.**

Sorted positive signed \((1,1,1)\) is exactly the additive family \(a+b=c\).
The row \((1,4,5)\) has
\[
                     \mu(F)=\min E={1\over28}.
\]
Every other eligible row satisfies the sharp band
\[
 {31\over392}\le\mu(F)\le\min E\le {6\over55},
\]
with lower equality only at \((1,7,8)\) and upper equality only at
\((1,10,11)\). Since \(31/392>6/77\), \((1,4,5)\) is the only row at or
below the old threshold; there is no \(6/77\) equality row.

The complete carrier set is
\[
 \{k(1,1,-1):0<|k|<3c/14,\ 3\nmid k\}.
\]
Its continuum physical bulk is the constant \(9/98\), while its minimum
network bulk ranges up to \(39/392\). The new mechanism is the lower as well
as upper deleted-third discrepancy:
\[
 L_f-{4\over7c}\le {2\over c}\sum_{3\nmid k}f(k/c)
 <L_f+{4\over7c}.
\]
This proves that the additive family is a cofinal body-phase obstruction,
not a finite exceptional list.

The primary and clean-room engines agree on 1,901 raw rows through height
223. The referee imports no candidate or repository geometry code and checks
the modular kernel, generic breakpoint integrals, sharp extrema, all three
parity sheets, and both cutoff endpoints.

```powershell
python -B 04-computation/lrc14_signed_111_sharp_obstruction_thm4445.py
python -B -O 04-computation/lrc14_signed_111_sharp_obstruction_thm4445.py
python -B 04-computation/lrc14_signed_111_sharp_obstruction_thm4445_independent.py
python -B -O 04-computation/lrc14_signed_111_sharp_obstruction_thm4445_independent.py
```

See [THM-4445](../../01-canon/theorems/THM-4445-lrc14-signed-111-sharp-obstruction-classification.md)
for the complete proof and scope.
