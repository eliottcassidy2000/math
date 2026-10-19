# Two complete c=0 morphology exclusions

This folder excludes maximum-span critical-root morphologies (2,4,0) and (3,3,0), and their duals, in a hypothetical seven-thin no-balanced-pair poset. All 9,648 labeled cores and 56,265 residual orders are covered. Peczarski's known six-thin theorem is explicitly used for short spans.

Run `python -O replay_small.py`. The full mathematical replay is standard-library-only; search additionally uses NumPy/SciPy. Read the main proof for the arbitrary-exterior lift, scoped certificate pattern, exact counts and attribution. The four larger preliminary branches are not included in the theorem. No full seven-thin theorem or novelty claim is made.
