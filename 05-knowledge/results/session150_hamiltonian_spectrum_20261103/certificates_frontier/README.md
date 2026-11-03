# Extending the finite Hamiltonian-path prefix to 90,001

**Conditional finite certificate.** Given the prior certificate that every odd integer H ≤ 84,259 except 7 and 21 is attained, every odd integer through 90,001 except those two is attained. The certificate is a list of labelled strong order-12 tournaments whose path counts cover every prime and every 7 times a prime in the new interval. Any other odd integer there splits as a product of two odd numbers ≤ 84,259, neither 7 nor 21; `check_prefix.py` finds and checks a split for each case. The ordinal sum of tournaments realizes the product of their Hamiltonian-path counts. This is a finite result, not a proof of eventual or global completeness.

The inherited 84,259 prefix follows from the companion `certificates_interval` packet and `01-canon/theorems/THM-4104-selected-order-eleven-strong-ear-solid-interval.md`.

**Encoding.** Each row of `witnesses.txt` is `(H, parent, cut)`. Parent vertices are 0–10; its 55 bits list pairs `(i,j)` with `i<j` lexicographically, low bit first; 1 means `i → j`, 0 the reverse. For added vertex 11, cut bit i set means `i → 11`, otherwise `11 → i`.

**Independent replay.** `verify.py` constructs adjacency from each row and uses start-vertex subset DP to count paths and transitive closure for strongness. `check_prefix.py` verifies each witness and checks every other odd target by direct factor search against the inherited prefix. Run `python3 check_prefix.py` and `python3 -O check_prefix.py` in this directory. Neither uses assertions or the search recurrence.

**Search scope.** The candidate pool used 5,000 fixed-seed pseudorandom order-11 labelled parents (each edge bit set at approximately 25% frequency), tested their nonconstant cuts using path-cover counts and retained candidates only for prime and 7×prime targets in 84,261–90,001. This selection procedure is not a census, and the replay requires none of its formulas.
