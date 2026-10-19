# Thue–Morse source exclusion and its arithmetic dichotomy

Read the [elementary positive-integer exclusion](THUE_MORSE_FINITE_SOURCE_EXCLUSION.md), [scalar p-adic extension](THUE_MORSE_SCALAR_MAHLER_EXTENSION.md), and [independent review](review/INDEPENDENT_REVIEW.md). Every fixed non-erasing binary block image, with a finite prefix, has a terminating source-dependent mismatch certificate for positive Mahler integers. The separate cited-input result classifies commuting block images as rational sources and noncommuting images as transcendental sources. These theorems leave arbitrary safe words and the Mahler 3/2 problem open.

From this directory run:

    python3 -B mahler_thue_morse_reset_exclusion.py
    python3 -B -O mahler_thue_morse_reset_exclusion.py
    python3 -B review/verify_independent.py
    python3 -B -O review/verify_independent.py

Match author outputs to mahler_thue_morse_reset_exclusion*.json and independent outputs to review/receipt*.json. The author runs 18,802 explicit checks and the independent verifier 68,642; the scalar theorem follows from the audited cited input and proof, not the finite tests. SOURCE_PROVENANCE.json records reviewed source hashes and publication changes.
