---
id: HYP-9520
title: "The eternally black values thin out. The odd values that are never the Hamiltonian-path count of a self-converse tournament (THM-4930) conjecturally have natural density zero, and (EXPLORATORY) only finitely many exist. Below 1131 there are exactly 83 (with 7 and 21). Their share falls from 24% below 225 to 16.5% on [225, 625) and 8.3% on [625, 1131), and is at most 3.0% on [1131, 2001), where only 13 values lack an explicit self-converse witness. Each block found so far is born at a single order: [225, 625) at order 8, [625, 1131) at order 9"
status: >
  OPEN (CONJECTURE (a): density zero). EXPLORATORY (b): finiteness.
  FINITE-EXACT: the eternally black values below 1131 (THM-4930 statements 3-5; 377 by complete order-12
  enumerations below 625; [625, 1131) by complete order-12 and order-13 enumerations below 1125, statement 5
  PRE-AUDIT, audit K2 pending).
  FINITE: explicit self-converse realisations for 422 of the 435 odd values in [1131, 2001).
  HEURISTIC: Busch windows lengthen.
  Audit K1 (2026-11-08): PASS WITH CORRECTIONS (title typing; cosmetic); MISTAKE-738.
session: mac-mini-2026-11-08-firstrealizers
source: 05-knowledge/results/eternal_colours_self_converse_values_20261108.md
scripts:
  - 04-computation/experiments/eternal_colours_20261108/sc_holes.py (+ .out)
  - 04-computation/experiments/eternal_colours_20261108/sc_odd13.c (+ sc_odd13_1124.out, sc_full1011_le2001.out)
  - 04-computation/experiments/eternal_colours_20261108/sc_witness.c (+ .out), sc_walk.c (+ sc_walk.out)
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

**Data.** Exact below 1131 (THM-4930).

| range of odd `h` | eternally black (attainable) | share |
|---|---|---|
| `[3, 225)` | 27 | 24% of 111 (27 of the 109 attainable values, 25%) |
| `[225, 625)` | 33 | 16.5% of 200 |
| `[625, 1131)` | 21 | 8.3% of 253 |
| `[1131, 2001)` | at most 13 (no realisation known; not decided) | at most 3.0% of 435 |

- **Birth orders.** Every eternally black value in `[225, 625)` is born at order 8, and every one in `[625, 1131)` at
  order 9. Below 225 the birth orders are 6 (3 values), 7 (11) and 8 (13).
- **Mod 4.** 66 of the 83 values below 1131 are `≡ 3 (mod 4)`, where only odd orders can supply a strong centre.
- **The 21 values in `[625, 1131)`** are 671, 683, 707, 719, 723, 755, 785, 795, 815, 851, 869, 887, 907, 931, 949,
  983, 1003, 1009, 1067, 1087, 1115. All were missed by the witness searches (`sc_witness.c`, `sc_walk.c`) and are
  absent from the complete order 11; the order-13 enumeration then showed every one of them black. Above 1131 the
  near-transitive search also missed 1551 and 1627, which the complete order 11 shows blue, so a missed value is not
  automatically black.
- **Without a realisation in `[1131, 2001)`:** 1131, 1145, 1231, 1245, 1263, 1391, 1459, 1475, 1535, 1567, 1583,
  1699, 1891 (11 of them `≡ 3 mod 4`). All are born by order 9. Deciding them needs orders 13-15: a strong centre
  `e ∈ [1125, 1875)` can have order 14, and `e ≥ 1875` order 15.

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

- **Order 12.** Done (audit K1, confirmed by `sc_odd13.c`): 377 is eternally black; every value below 625 is decided.
- **Orders 12 and 13 below 1125.** Done (`sc_odd13.c`): every value below 1131 is decided, and all 21 values left
  unresolved there by the witness searches are black.
- **The values in `[1131, 2001)` without a realisation.** Prediction (EXPLORATORY, from the pattern below 1131): they
  are all eternally black. Deciding them needs complete self-converse inventories above 1124: of order 13 (and 15 for
  1891) for the values `≡ 3 (mod 4)`, and of orders 12, 13 and 14 for 1145 and 1245 (`f(14) = 1125`).
- **The delay.** Is the delay from birth to the first self-converse realisation of an eternally blue value bounded?
  Below 1131 it is 0–3 except for 275 (5) and 711, 855, 1011 (4). 275, 711 and 855 are sandwiches around small
  strong atoms, their only kind of realisation; 1011 is a strong atom that first appears at order 13.
