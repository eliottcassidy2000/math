# Independent coefficient-box audit for the reserved THM-4434 candidate

**Status: VERIFIED-EXACT SUPPORTING AUDIT; NO THEOREM PROMOTION.**  This note
audits the analytic tail chain and independently reproduces its complete
maximum-coefficient-18 obligation.  It does not consume or alter the status
of reserved THM-4434, and it does not audit a finite speed-height head.

## Verdict

The requested analytic chain passes:

- the error-cube image is a centrally symmetric planar zonotope;
- its slice-width word is even and nonincreasing, including the physical
  endpoint-zero convention for an actual zero coefficient;
- the exact area is `9(a+b+c)/49`;
- the `h=1` and `h=3` rectangle rules give the stated unit and one-zero
  loads and error `2f(0)/3`;
- choosing `M=max_i|v_i|` gives `f(0)<=6c/(7M)`;
- `M>=19` has slope at most `142/931=15/98-1/1862`;
- the count intercept gives the threshold
  `c>=(308/31)S+4312/93`;
- the projected relation lattice has determinant `c`, its `l1` ball has the
  stated area, and planar Minkowski gives `S<4sqrt(c/3)`;
- `56056/93<603`, while the remaining even-norm tail starts with
  `g(58)=3023/372>0` and is increasing.

No analytic hostile was found.  The magnitude pattern `(1,1,2)` remains the
sharp, correctly excluded slope hostile at `2/7`.  The unique coefficient-box
equality is `(1,7,8)` at `15/98`, attained only on the enlarged closed speed
polytope; it is not a physical equality assertion.

## Independent finite calculation

The referee imports no repository implementation.  Instead of the producer's
polygon clipping, it:

1. generates primitive signed coefficient vectors directly from
   `[-18,18]^3`, modulo overall sign;
2. finds error-slice and speed-polytope vertices as exact plane/cube-edge
   intersections;
3. checks every signed coordinate placement before grouping by magnitude;
4. verifies the zonotope determinant sum, even slice word, both lattice-rule
   errors, both residue-load formulas, and the central-section bound at every
   rational speed vertex; and
5. compares all support-two results with the exact rectangle formula.

The exact universe and result are

```text
4,905 signed sectors,
308 magnitude patterns =293 full support +15 support two,
maximum slope 15/98 only at (1,7,8),
norm-four hostile 2/7.
```

All fifteen actual-zero-coordinate patterns are present.  Their largest
slope is

```text
A_(0,11,13)=128/1001<15/98.
```

The complete `(pattern, defect list, slope, intercept)` semantic digest is

```text
3be81c2522a1df6a146e50634754620103f9d2d8840d17f34c9e9a4954e849f7,
```

identical to the existing coefficient-box producer.

## Reproduction

```powershell
python -B 04-computation/lrc14_universal_scale_three_coefficient_box_thm4434_independent_referee.py
python -B -O 04-computation/lrc14_universal_scale_three_coefficient_box_thm4434_independent_referee.py
```

The two runs have raw-LF byte-identical output and execute `437,776` explicit
checks.  The implementation has no optimizable assertions.

```text
source
b326cf2ec6ae35aa4fd4c5434bab9028d17c7b21249047e5e34c01d76174c694

output
7df7679bc41ae62b206b40da774d68ac09595e947cd59b66742a935e747effaa
```

This audit establishes the coefficient-box and analytic-tail obligations
only.  Universal projection closure, theorem promotion, entry,
synchronization, and LRC(14) are separate status decisions.
