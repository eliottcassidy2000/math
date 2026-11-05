# Independent review: three-parameter ternary-clock family

**PASS**, no corrections requested, 2026-11-05.

Frozen source in `session197-hybrid-constructions/ternary-clock-family/`:

- PROOF.md: bb9249b7a0371ea3c81b24cb5b1de13e8f6b4fac0120d4ac789cdd79f25d9baa
- verify.py: 7733fc0a2eac672eadf55c001244489e64a3a3f383d50bc71874683e23c501d6
- core.py: 3204eeb1e9d2b05cac9bcc21c98f25697d10f925d6ff743393e945ce4e73615c
- Both validation/normal/optimized receipts: 6867645e25d1f028eff48a416f79891ba7074e60b8af3e9a76f22c625dd66bc6

Copied-source normal and optimized runs agree with the frozen receipt. No source file was edited.

## Uniform run bounds and interpolation

The quotient is a valid strong tournament; a->b->U->c->V->R->a is a directed spanning cycle. All positive transitive blowups are therefore strong.

The crossing covers remain bounded independently of all three module sizes:

- U's external crossing edges are a->b, c->b, c->V and R->b. They are covered by the singleton left vertex c and singleton right vertex b.
- V has only the singleton c as an external in-neighbor.
- R's crossing edges are b->U and c->V, covered by singleton left vertices b,c.

Internal edges of another transitive block do not create additional crossing types: that entire block has one fixed incidence type relative to the module under consideration. The covers thus genuinely give run bounds3,2,3 for every positive parameter choice, not merely on the interpolation grid.

For fixed other parameters, the clone-run formula supplies exponential bases1,...,3, then1,2, then1,...,3. Recovering coefficients with the fixed positive-parameter Vandermonde matrices preserves the remaining separate spans. Repeating this recovery proves membership in the18-dimensional tensor span. No simultaneous-independence assumption is hidden in that step. The complete18-point Cartesian grid consequently determines the function at every positive triple. The proposed formula belongs to precisely that span.

The fixed-port return changes only m: R beats a,b,c and loses to U,V, and a newly appended clone loses to the older clones. Reversing the labels inside U and V of length two gives exactly the previously audited seed. Changing p or q is correctly separated from this returning birth.

## Unique cancellation

At p=1 the binary coefficient is6. For p>=2 the denominator2*3^p-3*2^p is positive. A positive q would require the displayed ratio to be at least2. At p=4 the necessary inequality already fails. Writing its positive defect as E_p=2*3^p-9*2^p+6 gives E_(p+1)=3E_p+9*2^p-12>0 for p>=4, proving the exclusion at every larger p. The remaining ratios are4 at p=2 and12/5 at p=3. Therefore p=q=2 is the unique cancellation in this quotient family, giving A=60,C=-1.

This argument does not assert uniqueness among all tournament contexts or supply inverse targeting.

## Independent finite evidence

The new `clock-family-audit/independent.py` constructs the blowup directly from quotient labels and counts all Hamilton paths by literal DFS, without subset DP. It independently verifies all18 grid graphs and all seven off-grid controls.

It then solves the full18-by18 tensor interpolation system by exact rational elimination. All18 recovered coefficients agree with the displayed expression, including its five zero coefficients and13 nonzero coefficients. Quotient-level crossing edges and the fact that all cover vertices remain singleton modules are separately checked. This provides an independent implementation of both the finite certificate and its uniform-scope premises.

The all-order theorem follows from those uniform run bounds plus exact finite-dimensional interpolation. It is not numerical extrapolation from the off-grid checks.
