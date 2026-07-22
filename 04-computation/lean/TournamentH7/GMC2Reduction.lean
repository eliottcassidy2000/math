import TournamentH7.GMC2Formalization

/-!
Compatibility entry point.  The library module now lives at
`TournamentH7.GMC2Reduction` so the NC2/GMC(2) formalization spine can import
it and the default target can typecheck it.  This file retains the historical
single-file invocation `lake env lean GMC2Reduction.lean`.

The concurrent arithmetic-engine theorem names (`wick_expansion`,
`multinomial_dilate_modEq`, `face_sum_frobenius`, `face_sum_ne_zero`, and
`charge_radial`) are retained by the focused library modules imported through
`TournamentH7.GMC2Formalization`.
-/
