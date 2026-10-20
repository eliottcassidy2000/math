# Hollow-frame concentration supplement

Read `HOLLOW_FRAME_CONCENTRATION.md` for the proof and
`CONCENTRATION_REVIEW.md` for the independent proof-review record.

For uniform labelled external frames,

    Var(h)/(k!)^2 = 1/(2k) + O(k^-2),
    h/k! -> 1 in L^2 and in probability.

The proof bounds every summand of the inherited second-moment formula,
including indices growing with k. It also gives the immediate
second-moment maximum lower bound with strict odd rounding.

Reproduce the recorded 16,640 exact term checks:

    python3 -O verify_concentration.py --max-k 256 --output concentration_checks.json

Only the Python standard library is needed. These are analytic-formula
checks, not an external-frame census. No general maximum equivalent or
central limit theorem is asserted. The existing extrema and second-moment
packets are unchanged.
