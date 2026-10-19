# Mixed(2,3), c≤1 exclusion with known-six-thin dependency

Read `mixed23_c01_exact_exclusion_20261019.md`. All 5,181 cores are covered, with 5,894 candidate orders excluded using Peczarski’s 2008 six-thin theorem and 23,197 remaining orders exactly cone-certified. Independent replay passed.

Run `python -O verify_mixed23.py`; optionally run `python -O verify_obstacle.py` within `unfiltered_obstacle` to replay the unfiltered cone limitation. The external six-thin theorem is an explicit dependency, not reproved here. Cases c=2,3 remain separate. No full seven-thin theorem or novelty claim.
