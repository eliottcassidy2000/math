# The path splice (HYP-9862 proved), 3q + 2 / 3q + 4 against Collatz and 3x − 1, rows of doubling chains, and a continuation for a fresh digit

Session klein-2026-11-16-splice (2026-11-16), on the owner's directive: prove HYP-9862 with the splice injection; merge
how `3q + 2` and `3q + 4` relate to Collatz and `3x − 1`; change the continuation so that it can accept another fresh
digit; survey the prior sessions on the `+2` and `×2` recursions; think of the odd numbers in rows of chains growing by
doubling from each odd. Canon: THM-5256 (new); HYP-9862 marked PROVED.

## Status

- **PROVED** ([SPLICE.md](SPLICE.md), THM-5256). For induced `A, B ⊆ U` with `W = A ∩ B`: `H(U) ≥ Σ_σ h_A(σ) h_B(σ)` over
  the linear orders `σ` of `W`, by splicing Hamiltonian paths of `A` and `B` segment by segment between the shared
  vertices (merge lemma) and decoding by restriction. For `|W| ≤ 1` this is `H(U) ≥ H(A)H(B)`: **HYP-9862 is proved at
  every order**, and THM-5312 (A) follows. Two shared vertices buy exactly the order statistics (tight at 5 for two
  triangles on an arc). In the odd-cycle formula the statement is `Mix(U) ≥ (H(A) − H(A − v))(H(B) − H(B − v))`.
  Exhaustive check orders 4–7 (4,143,216 instances at order 7, 0 violations).
- **Merged** ([COLLATZ_ROWS.md](COLLATZ_ROWS.md)). The hosts `3q + 2`, `3q + 4` are `3(q + 1) ∓ 1`; in count coordinates
  the Collatz step is the digit-0 child `3q` and the `3x − 1` step is the digit `−2` child `3q − 4`; in half-counts the
  three children are `x ↦ 3x + 1, 3x + 2, 3(x + 1)`. The odd numbers in rows `L_c = {1 + c 2^s}` (THM-5126's lines):
  Collatz in count coordinates is the Syracuse map on row indices; `{7, 13} ⊂ L_3`, `{11, 21} ⊂ L_5`, 47 and 67 are
  heads; the digit-1 child of `q` lies on the row of the Syracuse image of `q`. The `+2` and `×2` recursions are the
  `B`- and `k`-directions of the one array `h + 2^k B` (THM-4865).
- **Surveyed** ([SURVEY.md](SURVEY.md)): 58 prior items on `+2` moves, doubling lifts, ternary children and rows, with
  their own status words.
- **FINITE-EXACT** ([CONTINUATION.md](CONTINUATION.md)). The fixed-root continuation `(E, O, g)` has `D = E − O ≥ −g`
  for every rooted tournament to order 8, `D = 0` for non-terminal roots, `D = −g` in every literal hit; weakening to
  `(X, g)` changes nothing in the literal families (pre-registered P-3 refuted as stated). Appending a sink erases `D`
  and the three digits of the count-5 sources (15, 17, 19) then share the complete continuation `(40, 40, 9)` — three
  consecutive odd counts with one fixed-root continuation, which DECODER.md forbids for `E ≠ O`. Second step: from the children 15, 17, 19 no further three-digit set with a common continuation exists at two added vertices (P-4 refuted; counts 45–61 reached, continuations differ).
- Pre-registration: [PREREGISTERED.md](PREREGISTERED.md) (P-1, P-2, P-5 confirmed; P-3 refuted as stated, confirmed in
  the modified form; P-4 refuted).

## Files

- `scripts/`: `shuffle_check.py` (the path-splice check), `profile.c` (rooted profiles, extension searches, `check`
  mode against the literal deposit values), `allsearch.c` (all rooted tournaments: `law`, `one`, default two-vertex),
  `common.py` (shared continuations of the three digits), `rows_collatz.py`.
- `receipts/`: `shuffle_check_4..7.out`, `profile_check.out`, `search1_deposit.out`, `search2_*.out`,
  `search2_*_sink.out`, `allsearch_4..6.out`, `allsearch_7_one.out`, `law_678.out`, `anycommon.out`,
  `second_child{0,2,4}.out`, `rows_collatz.out` / `.optimized.out` (byte-identical).
- `prediction_splice.txt`, `data/child{0,2,4}.masks`.
