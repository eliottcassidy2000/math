# Toothpick diagrams for A000568, the curried extension grammar, two encodings of {7, 21}, and first creators

Session klein-2026-11-16-tooth (2026-11-16), on the owner's directive: unequal growth and nonuniform contacts; a
toothpick diagram for each term of A000568 showing how the isomorphism classes relate; the first appearances of every odd
count other than 7 and 21 and their relation to each other; a full recursive grammar that naturally encodes the {7, 21}
impossibility; abstract structures as compositions of functions, curried. Canon: THM-5266 (new).

## Status

- **PROVED** ([GRAMMAR.md](GRAMMAR.md)). The count of a one-vertex extension is a quadratic form on the cube,
  `H(T + v_X) = h + ⟨B − e, x⟩ − xᵀ M x` (replayed on 62,422 extensions), so the spectrum obeys the exact recursion
  `S(n+1) = ∪_T H_T({0,1}^n)`; the uniform contacts (sink, source, sole exit, ordinal sum) are closed on `(h, e, B, P)` and
  compose to the near-join `ε_t ∘ σ^k`, the tower and the midpoint; a general contact needs the pair matrix, two of them
  the 3-path tensor ("nonuniform contacts cost state"). **The additive sub-grammar never outputs 7 or 21** because
  `𝔅(1) = {2^j}`, `𝔅(3) = {3·2^j}`, `𝔅(5) = {5, 6, 7, 9, 32}·2^j` and `12 ∉ 𝔅(9)`: 7 is not `2^m + 1`, 21 is not
  `3(2^m + 1)` — the holes are the Mersenne shadows of 9 and 27. The full grammar encodes them through multiplicativity
  and Busch's floors (THM-1370's route).
- **FINITE-EXACT** ([TOOTHPICK.md](TOOTHPICK.md)). The extension graph to order 9 (nauty canonical forms): contacts 2, 7,
  29, 215, 2591, 50357, 1656619; a class of order 9 has 8.65 distinct parents of 9 on average; children per order-8 class
  22–248; the transitive tournament is the most fertile only through order 4; canonical-parent tree branching max 227,
  43 % of order-8 classes childless. Diagrams [figures/toothpick_7.svg](figures/toothpick_7.svg) (532 toothpicks) and
  [figures/toothpick_8.svg](figures/toothpick_8.svg) (7,412): one toothpick per class, generations perpendicular, each
  toothpick crossing its least canonical deletion at the contact, first creators red.
- **FINITE-EXACT** ([CREATORS.md](CREATORS.md)). First creators are 85 % of order-9 classes (1, 1, 6, 33, 328, 5457,
  163676 at orders 3–9; 1,520 values, orphans 3 and the 63 pair as in THM-4922); the genealogy is nearly degree-blind
  (out-degree of the deleted vertex close to binomial, mode at the middle); dead-end appends are 0.37 % of order-9 first
  creators; one of the two first creators of 63 is a dead-end append over an order-7 base.
- Pre-registration: [PREREGISTERED.md](PREREGISTERED.md) — 7 predictions, 4 confirmed, 3 with refuted clauses (the
  genealogy is not dominated by sole exits, which the owner's "grammar" intuition would have suggested).

## Files

- `scripts/`: `quadform.py`, `holes_grammar.py`, `extgraph.py`, `spectrum9.c`, `spectrum9z.c`, `creators9.c`,
  `tree_creators.py`, `toothpick.py`, `sixty_three.py`, `nearjoin_check.py` (helper block).
- `receipts/`: `quadform.out`, `holes_grammar.out`, `extgraph.out`, `creators9.out`, `tree_creators.out`,
  `toothpick_7.out`, `toothpick_8.out`, `sixty_three.out`, with byte-identical `.optimized.out` for the Python scripts.
- `data/`: `counts.txt` (order, count per class in gentourng order; identical from the ascii and digraph6 streams),
  `cparent_3..9.txt` (least canonical deletion of every class, same order).
- `figures/`: `toothpick_7.svg`, `toothpick_8.svg`.
- `prediction_tooth.txt`, `SHA256SUMS`.

Reproduce (needs nauty's `gentourng`, `labelg`): `cc -O2 spectrum9.c -o spectrum9 && ./spectrum9 > counts.txt`,
`python3 extgraph.py` (27 s), `cc -O2 creators9.c -o creators9 && ./creators9`, `python3 quadform.py`,
`python3 holes_grammar.py`, `python3 tree_creators.py` (needs `counts_z.txt` = `counts.txt`), `python3 toothpick.py 7`.
