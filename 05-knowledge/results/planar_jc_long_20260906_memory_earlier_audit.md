# Independent audit of the complete valuation-eleven response

**Status: ACCEPTED: PROVED FINITE-ROW + INDEPENDENTLY AUDITED.**
This audits [the valuation-eleven theorem](planar_jc_long_20260906_memory_earlier.md).
Both coordinate changes start at row ten; the source deformation has
valuation at least eleven. The coordinate degree caps, source weight cap,
row-fifteen horizon, and actual depth projections are all retained.

The primary program uses 180 raw coordinate coefficients and 301 equations.
This audit first solves every bracket row with a different polynomial
particular solution, leaving 81 complete tangent parameters. It builds
the depth annihilators from literal generators `x^a u^b p^c y^e`, rather
than the primary's closed simplex formulas. The resulting system has 70
source equations and 91 depth equations. The literal-generator helper
is inherited from the independent valuation-twelve audit; no repository
mathematical program or producer output is imported at runtime.

All seven reduced source equations agree exactly. The smaller tangent
matrix has rank 71. A separate fixed `71 x 71` minor has polynomial inverse
`(I-Xi K)B0^(-1)`, with `K^2=0`. Its entire residual source matrix is constant.
A fixed `7 x 7` source minor has the different nonzero determinant

```
1165859945917954076962498289/20006217428243580727910813859840.
```

Thus this independent path certifies every specialization, not just the
generic fibre. Ten columns (one affine point and nine directions) have
polynomial tangent lifts satisfying all 161 equations, and the complete
ten-dimensional homogeneous coordinate kernel is also exhibited over
`Q[Xi]`. Source dimension nine and coordinate fibre dimension ten are
different quantities; the combined relative solution has dimension nineteen.

The independent inheritance checks substitute the sixth background row
into the literal bracket-five equation and predicted source-six equation,
using THM-4308 equations (6), (28), (32), and (33). The Student operator's
injectivity makes this a unique compatible row at each `Xi`. The recovered
boundary graph is checked to factor through the actual `G_m` coordinate,
and this row is verified to vary. The full actual prefix graph remains an
inherited THM-4426 dependency: the separate pinned background reconstruction
checks all 210 raw prefix coefficients and has 276 live gates. It is not
silently replaced by a constant background in this audit.

The complete weight ranks agree; weight 21 fails and weight 22 succeeds.
The audit checks the whole affine family, the inherited valuation-thirteen
transport, and the forced-zero deformation channels `py^5,y^6`.
The producer separately isolates the `py^5` onset between rows thirteen
and fourteen. This audit does not infer that onset from the row-fifteen
result alone.

Analytically, all nonlinear correction terms begin at source order twenty
and bracket order nineteen. The linear system therefore describes the
actual finite changes at the claimed horizon. The converse uses the full
raw equations and depth annihilator, so it is not merely a necessary
response statistic. Neither direction asserts later compatibility,
polynomial termination, arbitrary lower source valuations, or JC(2).

The audit has **5,868** always-active checks. Normal/optimized output
agreement and final source/output hashes are pinned by the manifest.
The relation matrix hash, independent of coordinate choices, is
`85c0211aa3f35ac2548683ad34e041f93e49456dba328e158d1d854a6c62a0c3`.

```
python3 -B 04-computation/planar_jc_long_20260906_memory_earlier_audit.py
python3 -B -O 04-computation/planar_jc_long_20260906_memory_earlier_audit.py
```
