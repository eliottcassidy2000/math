# 1, 3, 7, 21, 2, 18 in the grammar; neutral choices of the core; the fibres that cannot be changed

Session klein-2026-11-16-core (2026-11-16), on the owner's directive: with Chan–Pak reduced to the lanes `p` and `7p`,
study 1 and 3, their 7-multiples 7 and 21, the even prime 2 between them, and 18; change the core so it gives more
neutral choices without creating another free path run. Concurrent packets on the same words (other machine):
`ready_core_growth_and_mod18_20261117`, `one_port_neutral_continuations_20261117`; read first, cited, not repeated.
Canon: THM-5268 (new).

## Status

- **Arithmetic, typed** ([HOLES.md](HOLES.md)). `7 = 1 + 2·3` and `21 = 3 + 2·9`: the holes are the `k = 1` near-joins of
  the two rigid fibres with the next fibre's response; the pencils `1 + 2^k·3` and `3 + 2^k·9` miss only their first
  terms; `2` is the doubling both need and the centre of the first seam `(1, 3)`; `18 = 2·9 = 21 − 3`, as `6 = 2·3 = 7 − 1`;
  in the lane `7p = p + 6p` the holes are the two members whose shift is twice a wrong-fibre response.
- **PROVED** ([NEUTRAL.md](NEUTRAL.md) §1, THM-5268 (A)). Padding is transparent: the balanced arcs of `M_a ⇒ S ⇒ M_b` are
  the trivial adjacent-rank flips of the paddings and the balanced arcs of `S`. Hence the neutral structure of a fibre
  is that of its strong realisers, and the census below is complete for every `h < 125` (Busch's floor 125).
- **FINITE-EXACT, permanent** (NEUTRAL.md §2). The fibres with no nontrivial count-preserving flip at any order:
  1, 3, 5, 11, 13, 17, 19, 31, 35, 39, 47, 53, 55, 57, 63, 67, 85, 95, 113, 115, 121 — including the near-join
  exceptions, four of the five +2-window failures, the atom 63 and the tower-blocked 53, 67, 113. Order 9: 82,384
  balanced arcs, 27,456 nontrivial; the neutral-richest strong core has count 495 and ten balanced arcs; the richest
  fibre is 547. Spearman 0.613 between neutral edges and `|𝔅(h)|`.
- **Selection rule** ([CORE.md](CORE.md)). Rigid fibres cannot be changed without a free path run, and padding only
  doubles; targets `p`, `7p` should be written over flexible fibres. Not a proof of anything about the lanes.
- Pre-registration: [PREREGISTERED.md](PREREGISTERED.md) — 8 predictions, 2 confirmed in full, 6 with refuted clauses.

## Files

- `scripts/`: `neutral.py` (orders ≤ 8, from `profiles8.npz` of the seams packet), `neutral9.c` + `neutral9.py` (order 9),
  `holes_arith.py`, `nearjoin_check.py` (helper block).
- `receipts/`: `neutral.out` (+ `.optimized.out`), `neutral9.out`, `neutral9_py.out` (+ `.optimized.out`),
  `holes_arith.out` (+ `.optimized.out`).
- `prediction_core.txt`, `SHA256SUMS`.

Reproduce: copy `profiles8.npz` (seams packet, `profiles8.py`) and `counts.txt` (toothpick packet) here;
`python3 neutral.py`; `cc -O2 neutral9.c -o neutral9 && ./neutral9 && python3 neutral9.py` (nauty `labelg`, 10 s + 1 min);
`python3 holes_arith.py`.
