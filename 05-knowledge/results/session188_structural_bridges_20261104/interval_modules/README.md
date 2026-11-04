# Session 188: arbitrary interval-module reversals

`PROOF.md` proves:

1. Every simultaneous XOR family of original intervals in a transitive block
   quotient weakly increases every ordered path-cover count. Scalar equality
   occurs exactly when the new block quotient is transitive.
2. Interval patterns biject with all named quotient tournaments; exactly k!
   of the 2^(k choose 2) patterns preserve the scalar response.
3. In any fixed ambient quotient, no transitive source module can be replaced
   by any different tournament on the same names invisibly to the full named
   response. This includes arbitrary overlapping and nested original intervals.

`verify.py` is standalone and uses explicit exception gates. Run it normally
and under `python -O`; `normal.json` and `optimized.json` must agree.

The receipts retain complete named difference polynomials, including every
cubic and quartic, for two strong five-vertex ambient controls. The crossing
case changes H by six and has 144 nonzero coefficients. The nested case
preserves H but has a nonzero degree-four contrast with 136 coefficients.
The separate full whole-converse contrast is also retained and directly
verified on every attachment query.

The theorem does not claim a full-response collision or universal named
reconstruction. Nested reversals of nontransitive source modules, and the
full-response effect of transitive quotient permutations with nontransitive
blocks, remain outside its conclusion.
