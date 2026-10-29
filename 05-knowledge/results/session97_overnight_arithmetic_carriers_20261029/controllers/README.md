# Cheaper ternary controller: arithmetic completion and guarded retries

Five proof checkpoints, each with a stdlib exact checker and JSON receipt:

1. `INTEGER_TERMINATION.md` — every positive integer's (3,3) auxiliary tree terminates, via a sharp27/28 potential contraction and an integral branch guard. Checker: `check_integer_termination.py`; receipt:`checks.json`.
2. `FIVE_EXCEPTIONAL_CYLINDERS.md` — complete first-return comparison: five expanding cylinders, equality only at1, descent everywhere else. Checker:`check_exceptional_cylinders.py`; receipt:`exception_checks.json`.
3. `SOURCE_PRESERVING_RETRIES.md` — four repaired half-rays, the3 mod16 repair, exact unbounded-clock compression, and a complete hostile power family. Checker:`check_source_preserving_retry.py`; receipt:`retry_checks.json`.
4. `FAILURE_TRANSDUCER.md` —25 exact failure-to-failure charts versus their rank-one Haar quotient; source carries and clocks retained. Checker:`check_failure_transducer.py`; receipt:`transducer_checks.json`.
5. `PERIODIC_EXCEPTION_CERTIFICATE.md` — all-period arithmetic exclusion: only−1 is an integer periodic failure cycle; only−1,−3 have eventually periodic failure itineraries. Checker:`check_periodic_certificate.py`; receipt:`periodic_checks.json`.

Run each checker with `python3 -B`, then with `python3 -O -B`; outputs agree. No assertions carry verification logic. All claims distinguish native source arithmetic from source-measure calculations.

Still OPEN: whether the remaining aperiodic positive-integer failure itinerary exists; whether retries necessarily descend below their original source. None of the notes claims Collatz convergence, universal ROOT bounds, or an identification with tournament symmetry.
