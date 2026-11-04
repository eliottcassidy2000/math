# Collatz reset phases: certified local construction, not universal convergence

The packet proves six elementary results:

1. Every finite valuation prefix has an exact arithmetic progression of terminal reset exponents leading to 1. The compiler distinguishes first-hit routes from root padding.
2. Every permissible length-L reset phase has a prefix over {1,2,3,4}. Alphabet {1,2,3} fails already at the phase 14 modulo 18, so four is sharp among initial alphabets.
3. The least direct reset can be exponential in L. Tiny-alphabet prefixes attain the sharp exponent extremes.
4. For a periodic template centered at rho, all repeated-prefix reset phases are truncations of one ternary address j satisfying 4^j=1+3rho. This connects to the S179 binary mismatch parser through an exact identity.
5. This address is rational only for the all-1 and all-2 templates, with addresses 1/2 and 1. Every other template has a non-eventually-periodic ternary reset-digit expansion.
6. Every fixed template has recurring exponential direct-reset cost: the limsup of log_3(b_r)/(rL) is 1. The proof preserves the first-hit convention for the all-2 root template.

`PROOF.md` gives the proofs, inherited work, exact scope, and source links. `reset_phase.py` is standard-library-only.

Run:

```
python reset_phase.py > RESULTS.json
python -O reset_phase.py > RESULTS.optimized.json
cmp RESULTS.json RESULTS.optimized.json
```

The verified outputs match byte for byte: 5,460 exhaustive words, 10,920 literal first-hit root routes, all 59,048 permitted inverse phases through length 10, 150 periodic-center checks, and 213,575 counted explicit checks.

The default `completion(word)` is genuinely compact: it records formulas, the word, and a bounded ternary carry residue, without expanding the full carry or source. `completion(word, expand=True)` deliberately constructs the possibly enormous integers. The binary input lengths of the valuation letters still count toward the compact compiler's work. Tests with a letter of size 10^20 protect both the compact-completion path and the modular periodic-center compiler.

This constructs a source after a word is prescribed. It does not prove convergence for an arbitrary supplied source, improve natural-density coverage, or glue the local ranks into a global rank. Nonperiodicity of the reset digits does not imply randomness, normality, transcendence, or non-automaticity under general input encodings.
