---
id: HYP-9520
title: "The eternally black values thin out. The odd values that are never the Hamiltonian-path count of a self-converse tournament (THM-4930) have natural density zero, and conjecturally (EXPLORATORY) there are only finitely many. Below 625 there are exactly 61 (with 7 and 21; 377 undecided). Their share falls from 24% below 225 to 16% on [225, 625), and to at most 5.2% of [625, 2001), where only 36 values lack an explicit self-converse witness"
status: >
  OPEN (CONJECTURE (a): density zero). EXPLORATORY (b): finiteness.
  FINITE-EXACT: the eternally black values below 375, and those = 3 (mod 4) below 625 (THM-4930); 377 undecided.
  FINITE: explicit self-converse witnesses for 652 of the 688 odd values in [625, 2001).
  HEURISTIC: Busch windows lengthen.
session: mac-mini-2026-11-08-firstrealizers
source: 05-knowledge/results/eternal_colours_self_converse_values_20261108.md
scripts:
  - 04-computation/experiments/eternal_colours_20261108/sc_holes.py (+ .out)
  - 04-computation/experiments/eternal_colours_20261108/sc_witness.c (+ .out)
related: [THM-4930, THM-4922, HYP-9512, THM-4744 (Busch), THM-791, THM-4875, THM-1370, Chan-Pak Conjecture 6.2]
---

# HYP-9520 — the eternally black values

## Statement

Call an odd `h` **eternally black** if no self-converse tournament has `H = h`, and let `E` be the set of such `h`. By
THM-4930 this happens iff `h ≠ d²e` for every path count `d` and every `e ∈ {1} ∪ Σ_sc^str`.

**Conjecture.**
- **(a)** `E` has natural density zero among the odd numbers.
- **(b)** (EXPLORATORY) `E` is finite.

`E` contains every odd value that is not a path count. So (b) would imply that only finitely many odd values are not
path counts, a weak form of Chan–Pak's Conjecture 6.2, which predicts that the exceptions are exactly 7 and 21.

**Data.** Exact except 377 (THM-4930).

| range of odd `h` | eternally black (attainable) | share |
|---|---|---|
| `[3, 225)` | 27 | 24% of 111 |
| `[225, 625)` | 32 (377 undecided) | 16% of 200 |
| `[625, 2001)` | at most 36 (no witness yet; not decided) | at most 5.2% of 688 |

- Every eternally black value in `[225, 625)` is born at order 8. 49 of the 61 values below 625 are `≡ 3 (mod 4)`.
- The values in `[625, 2001)` without a witness are 671, 683, 707, 719, 723, 755, 785, 795, 815, 851, 869, 887, 907,
  931, 949, 983, 1003, 1009, 1067, 1087, 1115, 1131, 1145, 1231, 1245, 1263, 1391, 1459, 1475, 1535, 1551, 1567, 1583,
  1627, 1699, 1891 (30 of them `≡ 3 mod 4`).

## Heuristic (HEURISTIC)

- **Windows lengthen.** A strong realisation of `e` has order in the window `[n(e), B(e)]`. Here `B(e) ≈ 1.86 ln e`
  (Busch), while `n(e) ≈ ln e / ln ln e` (the least order grows like the inverse of `log(n!/2^(n−1))`). So the window
  contains more and more odd orders.
- **Black values live at the low end.** Below the conjecture's range, every eternally black value sits where every
  available odd order is near the low end of its strong self-converse spectrum, or where only even orders are available
  and THM-791 applies.
- **Large values sit in the middle.** For large `e`, some odd order `m` in the window has `f(m) ≪ e ≪ max H(m)`, i.e.
  `e` lies in the dense middle of order `m`'s spectrum. At order 9, for instance, 1,215 of the 1,679 odd values up to
  3357 are self-converse values.

## Tests

- **Order 12.** Settle the `h ≡ 1 (mod 4)` values in `[375, 625)` that remain open.
- **Witnesses.** Extend the near-transitive search (more reversed orbits, orders 13–15) to clear the values in
  `[625, 2001)` that remain unresolved.
- **The delay.** Is the delay from birth to the first self-converse realisation of an eternally blue value bounded?
  It is 0–3 for every blue value below 625 except 275, where it is 5. 275 is a sandwich around the order-5 atom with
  `H = 11`, the only kind of realisation it has.
