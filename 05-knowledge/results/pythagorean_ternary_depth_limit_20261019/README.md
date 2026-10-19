# Pythagorean ternary-depth limit and finite-memory brackets

19 October 2026. PROVED, with exact finite controls.

The fraction of depth-n Berggren nodes producing positive AP square-sum triangles converges to the irrational constant

    0.06039052312486240241155650348537631073311780098403...

with error at most 2/3^n. This closes the prior AP square-sum classification note's specific ternary-depth limit gap. It is a different sampling measure from the one-third asymptotic density under a hypotenuse bound.

At suffix depth d>=1, a sound interval classifier marks exactly two of the 3^d suffixes unknown. Every other verdict is correct for every possible earlier prefix. The unknown fraction is 2/3^d under uniform fixed-length sampling. This is an abstention fraction and a suffix-only, arbitrary-earlier-length guarantee; it is not an error rate or a lower bound for every classifier given the total input length.

The proof uses ordered interval cylinders and orientation-aware ternary digits. It does not assume the finite root-node levels are invariant under the wall-swapping symmetry, which sends root slope 2 to 3. Irrationality follows from the quartic boundary and its nonperiodic Möbius address.

Read TERNARY_DEPTH_LIMIT_AND_FINITE_MEMORY.md. Its logical predecessor is the [positivity-language proof](../dynamics_labeling_entry_20261019/BERGGREN_POSITIVITY_NONREGULAR.md), whose shrinking-cylinder argument is cited explicitly. The irrationality direction does not need the external quadratic-sufficiency theorem.

Run verify_ternary_depth_limit.py with Python 3. Its 200,465 checks reproduce the included ternary-depth-receipt.json. Python -O gives byte-identical output. The receipt includes a depth-128 rational interval of width 2/3^128 and a fully integer-certified 50-place decimal truncation. Finite checks audit the formulas and conventions; the text proves convergence and irrationality.

No transcendence claim is made. No Git operation or remote publication is claimed by this staged package.
