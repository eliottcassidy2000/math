# The sink trick with a larger prefix and a fresh sink before each digit (FINITE-EXACT, negative beyond count 5)

Session klein-2026-11-16-hex, continuing CONTINUATION.md of path_splice_rows_continuation_20261116. Tool
`scripts/allsearch2.c` (all rooted tournaments of order `n`, `k = 1, 2, 3` added vertices with every contact pattern
and every orientation among the new vertices); receipts `as2_5_2.out`, `as2_6_2.out`, `as2_7_1.out`, `as2_4_3.out`,
`as2_5_3.out`.

**Setting.** With a trailing sink the parity split `D = E − O` of the fixed-root continuation vanishes and the
complete continuation is `(X, g) = (E + O, g)`; a fresh sink before each digit therefore means: search the prefix for
three extensions of counts `3H, 3H + 2, 3H + 4` with a common `(X′, g′)`. The previous session found this only for
three prefixes of count 5 (two added vertices, order 6).

## Results

| scan | rooted bases | with a `3H + 2` extension | with a `3H + 4` extension | with a common-continuation triple | counts of those |
|---|---|---|---|---|---|
| order 5, two added vertices | 60 | 50 | 45 | 2 | 5, 5 |
| order 6, two added vertices | 336 | 276 | 246 | 3 | 5, 5, 5 |
| order 7, one added vertex | 3,192 | 1,631 | 1,561 | 2 | 5, 5 |
| order 4, three added vertices | 16 | 16 | 12 | 1 | 5 |
| order 5, three added vertices | 60 | 50 | 45 | 2 | 5, 5 |

- **Larger prefix: none.** Every rooted base of order ≤ 6 (two added vertices), of order 7 (one added vertex) and of
  order ≤ 5 (three added vertices) whose three digits share a common `(X′, g′)` has `H = 5` (pre-registered P-A and
  P-C **confirmed**). The count-5 bases are the three of the previous session plus their order-7 and order-5 relabelled
  forms; with three added vertices the base `(5, 2, 4, 3)` (`D = −2`) also reaches a common complete `(E′, O′, g′)`.
- **The parity split is not the digit's resource.** At order 6, bases with `D ≠ 0` have a `3H + 2` extension in
  217 of 262 cases (83 %) and bases with `D = 0` in 59 of 74 (80 %); for `3H + 4` the fractions are 72 % and 77 %.
  Pre-registered P-B (a factor of two in favour of `D ≠ 0`) is **refuted**: small digits are equally available with
  and without a parity split; what the sink buys is only the coincidence of the continuations, and that coincidence
  happens at count 5 alone in this range.
- By count (order 6, two added vertices): every base of count 3, 5, 9, 11, 13, 19, 23, 25, 29, 31, 37, 43 has both a
  `+2` and a `+4` extension; count 15 has neither (`D ≠ 0`: 24 bases, 0 and 12) or none (`D = 0`); count 17 has `+2`
  only, 27 has neither, 33 has `+2` (30) but `+4` (6), 41 has `+2` only, 45 has 6 of each. So the digit *values* are
  plentiful; the shared continuation is the scarce thing, as CARRY.md and DECODER.md found for the literal sources.

## Reading

The fresh-sink protocol is exactly the weak continuation on the prefix, and the exhaustive scans say that no prefix of
count other than 5 carries a synchronized three-digit step with up to three added vertices. The count-5 prefixes are
the strong 5-tournaments with one 3-cycle removed, i.e. the smallest nontransitive objects; their children 15, 17, 19
are then dead ends (the previous session's second step). Together with wave 45's charge principle (reusable digit ⟺
shared conserved charge) this says the ternary digit system cannot be grown from small seeds by appending vertices with
a sink-erased continuation; the deposit families of fresh_unit_prime_gluing remain the only synchronized digits, and
they are one-step.
