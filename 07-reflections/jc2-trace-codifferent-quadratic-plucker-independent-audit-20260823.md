# m=2 W2 Plucker probe (scratch; 2026-08-23)

Status: **FINITE-EXACT reduction; full rank-two decision remains open.**

Universe: `W2` is the `QQ`-span of all degree-one and degree-two monomials in
`(R0,R1,R2,U,V,E)` for the `m=2` tower.  All polynomial and linear-algebra
claims below are exact over `QQ`, except the explicitly labelled finite-field
quadratic relaxation.

## Exact reductions

1. The 27 displayed monomials have source-polynomial rank 24.  The constant is
   already in their span:

   `1 = V - R0*E + (3/2)(R2*V - R2)`.

   Hence `W2/QQ` has dimension 23.  Dropping `V` from the printed 24-element
   source basis gives the quotient basis used by the script.

2. On `Lambda^2(W2/QQ)` (dimension 253), the nonconstant coefficient map has
   rank 111 and the full bracket map has rank 112.  Thus the affine fibre
   `{c : J(c)=1}` has dimension 141.

3. A sparse exact point `c0` in that affine fibre has skew rank 6, not 2:

   `R0^E + (9/2)R0^(R1U) + 18 R0^(R2E) - 30 R1^(R1E)`
   `+ (45/8)R1^(U^2) + (9/2)(R0R1)^(UV)`.

   This is the cheapest positive *linear* signal, but it is not decomposable.

4. The fixed-`x` slice is exactly inconsistent: the map `G -> J(x,G)=G_y`
   has rank 21, while adjoining the target constant raises the rank to 22.
   Equivalently, there is no `G in W2` with `J(x,G)=1`.  As a positive control,
   adjoining `y` immediately gives the rank-two solution `x^y` with bracket 1.

5. Every hypothetical solution can be normalized, modulo constants, to
   `F=x+H`, `G=E+K`, where `H,K` lie in the 21-dimensional zero-gradient
   quotient.  This is the exact 42-variable bilinear/Segre chart.  A direct
   total-degree-3 XL test did not produce a unit certificate.

## First structural stopping obstruction

The lowest possible global Pfaffian separation would be a linear functional
on `Lambda^4` that is constant and nonzero on `c^c` throughout the affine
bracket fibre.  Equivalently it would separate `c0^c0` from

`span(c0^ker(J), ker(J)^ker(J))`.

The exact `F_1009` reduction has 10,152 generating rows, span rank 8,178 in the
8,855-dimensional `Lambda^4(F_1009^23)`, and `c0^c0` reduces completely to
zero against that span.  Therefore this prime admits no one-functional
quadratic/Pfaffian separator.  This is a hostile diagnostic, not a
characteristic-zero existence result: a full decision now requires either a
higher-degree elimination certificate or a structural case split in the
42-variable normalized chart.

Replay from the worktree root:

`python 04-computation/jc2_trace_codifferent_quadratic_plucker_independent_audit_20260823.py`
