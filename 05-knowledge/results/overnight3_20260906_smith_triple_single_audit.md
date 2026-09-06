# Independent audit of the ternary four-node residue obstruction

**Status: INDEPENDENT MATHEMATICAL AND EXACT AUDIT PASS.** This audits
[the complete unit-family theorem and diameter claim](overnight3_20260906_smith_triple_single.md).
It does not extend the family to arbitrary nested depths.

The row and degree types agree with the value/first-Hasse-derivative
observer on degree below eight. Clearing the two rows at zero leaves the
six displayed rows on degrees two through seven. A minor has common
factor `3^(2*(sum degrees-number of derivative rows))`; all remaining
terms have exactly the stated polynomial weight. Thus the coefficient
floors give valid lower bounds for arbitrary three-adic units, not only
the finite numerical controls.

The independent [companion](../../04-computation/overnight3_20260906_smith_triple_single_audit.py)
reconstructs **all 923 minors and 5,587 coefficients** using SymPy's
polynomial-domain determinant implementation, without importing the
producer or either of its determinant algorithms. It obtains the same
complete surviving residue table. In particular the rank-four ideal has
the sole possible leading residue `a^2 b(1+a)` at valuation 27. Every
other rank-four minor is divisible by `3^28` when that residue vanishes.
The separately reconstructed minor `3^24 A^4(A-1)^4` has valuation exactly
28 for every unit `a`; consequently higher cancellation in the critical
minor cannot increase the ideal valuation further.

The five other determinantal minima each have an actual unit witness.
Successive differences yield the ordered spectrum
`(0,0,2,6,7,12+kappa,15-kappa,22)`. This checks both lower and upper
directions. No determinant-only or preserved-prefix inference is used.
The intrinsic residue `chi` is well defined at either closest-pair base,
and affine coordinate changes cancel between its numerator and denominator.

The independent companion also uses **integer Smith normal form** on 24
signed unit controls and every one of the 2,600 positive-common-depth
integer configurations with diameter below 81. This differs from the
producer's modular elimination. It finds the same eleven metric trees,
each with one spectrum. Translation to least node zero is integral and
unimodular. At common depth zero, integral CRT splits residue-class
observers of size at most three; the proved three-node theorem determines
each from its metric. Together these arguments cover every smaller integer
diameter, including cases absent from the finite common-depth-positive
enumeration. The pair of diameter 81 therefore attains the stated minimum.

Both hostiles have determinant valuation 64 and largest exponent 22.
Their kernel exponents modulo `3^13` are nevertheless 53 and 54. This
independently verifies an actual observer consequence of the lost residue.
The 24 unit controls also verify the independently incoming
[THM-4435 Hermite formula](../../01-canon/theorems/THM-4435-four-node-metric-blindness-and-universal-hermite-precision.md)
`max_i{2v_3(F'(x_i)),3v_3(F'(x_i))-v_3(F''(x_i))}=22`, with ordinary
second derivative and the zero-derivative case handled separately.
The finite diameter result does not classify all four-node metric trees;
the residue repair is asserted only on the stated full unit family.

Reproduce from the repository root:

```text
python3 -B 04-computation/overnight3_20260906_smith_triple_single_audit.py
```

The [preserved output](overnight3_20260906_smith_triple_single_audit.out)
records 8,231 explicit checks. The symbolic identities and complete finite
diameter reduction carry the proof scope; the check count is only execution
evidence. No external theorem beyond the named current three-node canon is
used in the mathematical audit.
