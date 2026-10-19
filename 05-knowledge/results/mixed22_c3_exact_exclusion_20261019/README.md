# Mixed-(2,2), c=3 exact exclusion

Read `mixed22_c3_exact_exclusion_20261019.md` for the scoped theorem and proof.
It covers all 82,987 degree-valid canonical cores of this morphology, using
51,098 analytic cases and 31,889 residual exact-certificate cases. It does not
claim the c=4 morphology or the general seven-thin conjecture.

Run `python -O verify_mixed22_c3_exact_exclusion_20261019.py` for the independent
standard-library replay. No numerical optimizer is needed. The proof does not
use a short-span/six-thin theorem: all three common neighbors lie between the
roots, so their span is already at least four.
