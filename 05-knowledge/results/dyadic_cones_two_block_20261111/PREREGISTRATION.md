# Pre-registration (mac-mini-2026-11-11-rays)

Protocol note. The two runs below were launched at 18:37 and 18:38 (local time), before this file was written;
no output of either had been read when it was written (both output files were empty), and it is pushed before
they are read. The registration is therefore of predictions about runs in progress, not of runs not yet started.

Notation (THM-5126 (A)). A rooted tournament `(L, v)`, `B = L - v`, `O = N+(v)`, has offset
`c0 = H(B) - sum_{pairs of O} W_ef` and response `M = sum_{y in O} U_B(y)`, and `H(L o 0^t) = c0 + M 2^t` for `t >= 1`.
A *cone root* is a root whose response is a power of two; a root with response `c 2^j` (`c` odd) and offset `r`
is a root of type `(c, r)`.

In sample when this is written: all roots of order at most 9 (every `B` among the 6,880 classes of order 8 and
every `O`): cone roots exist for 181 positive offsets, among them every odd `r <= 199` except 27, 79, 99, 155,
195, 199, and for 548 negative offsets, among them every odd `r' <= 199` except 117 and 169; offsets `1 + 2^e`
have cone roots for `e = 1..7` and `e = 9`, not for `e = 8`; table of types `(c, r)` for odd `c <= 63`
(`out/blocks_roots_order9.txt`).

Runs registered: `cone_census 8 1 20000` and `block_census 8 1 2047` on the class list of order 8 (all
one-vertex extensions as bases: every root of order 10).

Predictions:

- (Q1) Among the roots of order at most 10, every odd `r` with `|r| <= 199` is the offset of a cone root.
- (Q2) The offset `257 = 1 + 2^8` has a cone root of order 10.
- (Q3) Every odd `r` with `1 <= r <= 999` is the offset of a cone root of order at most 10.
- (Q4) For every odd `c <= 15` and every odd `r` with `|r| <= 63` there is a root of type `(c, r)` of order at
  most 10.

(Q3) is the boldest: at order 9 the positive offsets below 1000 with a cone root are fewer than 181.

Falsifiers: one missing offset in the stated range. Each prediction is scored separately.

## Outcomes

`out/cones_roots_order10.txt`, `out/blocks_roots_order10.txt`; summaries in `out/summary_cones.txt`,
`out/summary_blocks.txt`. Every cone witness was recounted literally (`verify_cones.c`: 5,177 witnesses of order
10, three tail lengths each, no failure).

- (Q1) **CONFIRMED.** Every odd `r` with `|r| <= 199` is the offset of a cone root of order at most 10; in fact
  every odd `r` with `-791 <= r <= 647`.
- (Q2) **CONFIRMED.** `257` has a cone root of order 10 (response `2^10`); so do `513` and `1025`.
- (Q3) **FAILED.** The positive offsets `649` and `825` have no cone root of order at most 10 (the only two below
  1000; the negative offsets missing below 1000 are `-793, -809, -891`).
- (Q4) **CONFIRMED.** Every odd `c <= 15` and every odd `r` with `|r| <= 295` is a root type of order at most 10
  (the registered range was `|r| <= 63`); every odd `c <= 31` and `|r| <= 275`.

Not predicted: 243055 of the 261888 possible types with odd `c <= 255` and odd `|r| < 2047` occur among the roots of
order at most 10 (92.8 percent).

## Erratum and additions after the audit (2026-11-11)

- **Erratum (audit m7).** The sentence above, "(Q3) is the boldest: at order 9 the positive offsets below 1000 with
  a cone root are fewer than 181", is wrong: they are exactly 181 (every positive cone offset at order 9 is below
  1000; the largest is 597). The sentence is left as registered. The predictions and their scoring do not depend
  on it.
- The audit checked that the registration commit (`18571fc3da`) precedes the output commit (`aaf352fa43`) and that
  the prediction text is unchanged between them, and confirmed the four outcomes with its own census.
- Witnesses for (Q4) and for the wider boxes: the types with `c = 9, ..., 31` that occur only at order 10 had no
  witness on file at the pre-audit checkpoint (audit M4). `out/types_c9_to_c31_roots_order10.txt` now holds them
  (3,598 witnesses, 10,794 literal counts, no failure: `out/verify_types_order10_c9_to_c31.txt`).
