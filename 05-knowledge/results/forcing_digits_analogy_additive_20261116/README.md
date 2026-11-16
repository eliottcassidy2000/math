# Forcing the midpoint but not the +2 child, the tower criterion, the 13/41 analogy, additive gluing, port moves at orders 7 and 8, and the windows around 10007

Session klein-2026-11-16-force (2026-11-16), a long session on the owner's directive: force both the midpoint and the
+2 difference from an arbitrary prime rather than by search; classify the charge-carrying port moves at orders 7 and
8; find an additive construction like the root graft; make the analogy 13, 41 : 7, 21 as powerful as possible (wave 50's
sterile forks 11, 29, 35 and ready seeds 13, 41); explain the windows around 10007, 10009, 30011, 30013, 30021, 30027
in the spirit of approximations that select the integer. Canon: THM-5260, THM-5261 (new).

## Status

- **PROVED** ([FORCING.md](FORCING.md), THM-5260). (A) A history-blind gadget (substitute `T` for a vertex `τ` of `Q`)
  gives a count depending only on `H(T)` iff `τ` is on no cycle of `Q`, and then it multiplies: the midpoint `3p` is
  forced (`C_3 → T`), the `+2` child never is — it must read the realiser. (B) The ordinal-sum laws for the root
  response, end count and charge; a transitive suffix of length `k` sends `B_t ↦ 2^k B_t`, `e_t ↦ 0`. (C) The tower
  criterion: `3h + 2` is tower-forced at some order iff a strong realiser has a vertex with `2^k B_t = h + 1`; a rooted
  doubling child `2h − 1` exists iff `a_t = 0` or `2^k B_t = h − 1`. FINITE-EXACT table for all counts below 125.
- **PROVED** ([ADDITIVE.md](ADDITIVE.md), THM-5261). Gluing without mixed odd cycles: `H(U) = Σ_{X ∩ Y = ∅} H_A(X) H_B(Y)`;
  additivity `p + q − 1` exactly when all cycles of both pieces pass through one shared vertex and the mixed cycles are
  even; the census of pairs on a shared arc.
- **FINITE-EXACT** ([ANALOGY.md](ANALOGY.md)). The all-order charge and digit spectra read on 7, 21, 13, 41, the sterile
  forks and THM-5272's doubling failures; wave 50's P111 tested in the rooted-append sense.
- **FINITE-EXACT** ([PORTS78.md](PORTS78.md)). The charge-carrying port moves at order 7 (two added vertices) and order
  8 (one added vertex).
- **PROVED arithmetic** ([WINDOWS.md](WINDOWS.md)). 30011, 30013 are the twin-prime `+2, +4` children of the hidden
  generator 10003 = 7·1429; 30021, 30027 the forced midpoints of the twin primes 10007, 10009; the Tribonacci row is
  `round(c τ^m)`; the selection principle and why hadamard's placement does not yet select.
- Pre-registration: [PREREGISTERED.md](PREREGISTERED.md).

## Files

- `scripts/`: `module_blind.py`, `padding_check.py`, `digits.c`, `strongB.c`, `forcing_table.py`, `glue.c`,
  `glue_features.py`, `additive.py`, `portmove.c`, `tribonacci.py` (inline in the receipt), `numbers.py` (inline).
- `receipts/`: `module_blind.out`, `padding_check.out`, `digits_4..9.out`, `strongB_3..9.out`, `forcing_table.out`,
  `glue_*.out`, `glue_features.out`, `additive.out`, `portmove_7.out`, `portmove_8one.out`, `numbers.out`, `tribonacci.out`.
- `prediction_force.txt`.
