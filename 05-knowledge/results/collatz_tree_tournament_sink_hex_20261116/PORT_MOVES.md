# The binary port-move search, generalised (FINITE-EXACT, negative to order 6)

Session klein-2026-11-16-hex. Tool `scripts/odometer.c`; receipts `odometer_345.out`, `odometer_6.out`.

Wave 45 of chanpak_waves_20261113 (§7) argued that a fresh digit needs an additive **port move** — an operation
sending a ready `(r, a)` to `(2r + 1, a)` (binary) or `(3r + j, 3a − 2)` (ternary) — and tested one family,
`w45_tworoot.py`: two roots `x_1, x_2`, letter `i` appends a vertex that beats only `x_i` (or its mirror), and the
binary odometer condition `H(w i) = 2H(w) − 1 + a_i` for all words `w`; negative for orders ≤ 5. P107 (a root-moving
move on the ready seeds 13, 23, 41, 59, 75, 95) is queued there and not run here.

**The generalised search.** A letter `i ∈ {0, 1}` appends a vertex `z` that beats exactly an arbitrary subset `P_i`
of the base vertices and either beats (`r_i = 1`) or loses to (`r_i = 0`) every previously appended vertex; wave 45's
letters are `P_i = {x_i}` (or the complement) with `r_i = 0`. The odometer condition is `H(w i) = 2H(w) + d_i` for
every word `w`, with `d_0 ≠ d_1` (the digits need not be `−1 + a_i`). Over every tournament class of the base order
(all labellings are covered because `P_i` ranges over all subsets) and every unordered pair of letters with distinct
single-letter digits:

| base order | classes | letter pairs with distinct digits | surviving all words of length 2 | length 3 | length 4 |
|---|---|---|---|---|---|
| 3 | 2 | 128 | 0 | 0 | 0 |
| 4 | 4 | 1,556 | 0 | 0 | 0 |
| 5 | 12 | 20,952 | 0 | 0 | 0 |
| 6 | 56 | 431,736 | 0 | 0 | 0 |

**No pair survives even the four words of length two.** So a binary odometer made of history-blind appends (contacts
fixed on the base, one rule towards the appended past) does not exist at orders ≤ 6, strengthening wave 45's negative
by a factor of the whole subset lattice. This is consistent with wave 45's charge principle (a reusable digit needs a
shared conserved charge, Hankel rank 2) and with THM-5302's seed monoid: the doubling move `H ↦ 2H − 1 + a` of the
rooted append freezes its charge `a`, and a second letter with a different frozen charge breaks the first one's
rank-2 recurrence at the second step. The remaining candidate is P107's root-moving move, which is outside this
search space (the root changes, so the letters are not history-blind).

**Something better that emerged.** The sink trick (SINK_DIGITS.md) is a ternary port move with the continuation
erased; it works once, at count 5, and does not iterate. The tree tournament of COLLATZ_TOURNAMENT.md is the opposite
regime: a one-node-at-a-time append with history-dependent contacts (the ancestor chain), whose count is 1 or grows
exponentially with no small digits at all. Between the two — history-blind appends (no odometer) and ancestor-chain
appends (no small digits) — the digit programme needs contacts that depend on a bounded part of the history, which
is exactly what the deposit's marked modules are.
