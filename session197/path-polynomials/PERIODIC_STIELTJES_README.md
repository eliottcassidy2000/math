# Repeated physical wrapper words and Stieltjes moments

The reviewed theorem in `PERIODIC_STIELTJES_WORDS.md` proves every Hankel minor nonnegative for any fixed parity-returning word of triangle wrapping and single-vertex transitive padding, repeated from K1 or Tr2. It includes fixed compositions of transitive mirror contexts. The proof factors the exact production matrices into totally nonnegative pieces and proves the output-array/Cauchy-Binet transfer.

Scope is essential: arbitrary initial profiles are excluded, an odd-padding period is covered only after squaring, and a one-time-prefix graph gives the explicit negative determinant -48. No general real-rootedness, strict positivity or arithmetic target coverage is claimed.

Run from this directory:

    python3 verify_periodic_stieltjes_words.py
    python3 -O verify_periodic_stieltjes_words.py
    python3 periodic-stieltjes-review/independent.py

The frozen source receipts and the verbatim independent PASS report/checker are included. The review contains a complete copied source replay under its `source/` subdirectory. `PERIODIC_STIELTJES_DEPENDENCIES.json` binds the previously published dependencies, and `PERIODIC_STIELTJES_MANIFEST.json` binds every new payload file. The broader context-growth investigation is separate and is not included in this release.
