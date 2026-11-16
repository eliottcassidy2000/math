# Changing the continuation so that it accepts another fresh digit

Session klein-2026-11-16-splice. Tools `scripts/profile.c` (rooted profiles `(H, E, O, g)` and the one-/two-vertex
extension searches), `scripts/allsearch.c` (all rooted tournaments of an order), `scripts/common.py`; receipts
`profile_check.out`, `search1_deposit.out`, `search2_{triple,quad,deposit}.out`, `search2_*_sink.out`,
`allsearch_{4,5,6}.out`, `allsearch_7_one.out`, `law_678.out`, `anycommon.out`, `second_child{0,2,4}.out`.

The continuation of fresh_unit_prime_gluing_20261116 (BRIDGE_RESET.md) is the fixed-root profile `(E, O, g)`: `E`
and `O` are the sums, over root-start paths `P` with an even (odd) number of vertices, of the Hamiltonian-path counts of
the complement (`H(∅) = 1`), and `g = H(T − root)`. Write `X = E + O` and `D = E − O`. The ternary deposit realises
`3H, 3H + 2, 3H + 4` with the common tripled continuation `(3E, 3O, 3g)` on one family, and its note proves that the
nonzero children cannot be sources of the same family (`H = 2E − g` there, the children have defect `d·u`); CARRY.md
and DECODER.md show that no second fresh digit appears with the tripled continuation, and that no three consecutive
odd counts share any `(E, O, 3)` with `E ≠ O`.

## 1. The tool reproduces the literal families

`profile check` recomputes the ternary-deposit source `(163, 84, 89, 5)` and its three children `(489 + d, 252, 267, 15)`
exactly (receipt `profile_check.out`). The one- and two-vertex extension searches (fixed root, every contact pattern,
no old-arc change) on the three literal sources give:

| source (order, root) | profile | `d = 0` hits (strong / tripled `(E,O,g)` / tripled `(X,g)`) | `d = 2` | `d = 4` |
|---|---|---|---|---|
| triple splice (6, 4) | (17, 9, 12, 3) | 28 (12 / 4 / 4) | 32 (24 / 6 / 6) | 0 |
| quadruple splice (6, 5) | (15, 11, 12, 1) | 42 (34 / 6 / 6) | 0 | 14 (2 / 2 / 2) |
| ternary deposit (8, 7), two vertices | (163, 84, 89, 5) | 8 (0 / 4 / 4) | 2 (2 / 0 / 0) | 16 (0 / 2 / 2) |
| ternary deposit, one vertex | | 2 (2 / 2 / 2) | 0 | 4 (4 / 1 / 1) |

In every row the hits with tripled `(X, g)` are exactly the hits with tripled `(E, O, g)`: **weakening the
continuation to `(X, g)` changes nothing for these sources**. The deposit's `+2` child needs its old-arc reversal (the
two-vertex `d = 2` hits have continuation `(329, 338, 9)`), as its note says.

## 2. The parity split is not a free coordinate: `D ≥ −g`

Over all rooted tournaments of orders 4–8 (every root; receipt `law_678.out`, counts by `[strong][root terminal]` of
`D + g < 0, = 0, > 0`):

| order | non-strong, non-terminal root | non-strong, terminal root | strong (root terminal) |
|---|---|---|---|
| 6 | (0, 0, 74) | (0, 28, 24) | (0, 36, 174) |
| 7 | (0, 0, 405) | (0, 116, 200) | (0, 251, 2220) |
| 8 | (0, 0, 3727) | (0, 803, 2446) | (0, 2768, 45296) |

So `D + g ≥ 0` always (FINITE-EXACT to order 8; `D + g = Σ_(|P| ≥ 2) (−1)^(|P|) H(T − P)`, the singleton path giving
`−g`), `D = 0` exactly when the root is not in the terminal strong component (the terminal-component lemma of
DECODER.md: an ordinal suffix after the root cancels in the signed sum), and `D = −g` in a large terminal-root class
that contains every literal hit of §1 (`(84, 89, 5)`, `(252, 267, 15)`, `(9, 12, 3)`, `(27, 36, 9)`, `(11, 12, 1)`,
`(33, 36, 3)`, `(260, 269, 9)`, `(329, 338, 9)`: all have `D = −g`). The pre-registered P-3 is therefore **refuted as
stated**: in the fixed-root literal families `D` is a function of `g`, and the exhaustive scan of all rooted
tournaments of order ≤ 6 with two added vertices finds no source with the three digits under tripled `(X, g)` either
(`allsearch_6.out`: 0 weak triples, 0 strict triples, not even a `{+2, +4}` pair).

## 3. The change that works: make the root non-terminal

Since `D = 0` whenever an ordinal suffix follows the root, **appending a sink** (or any transitive suffix of length
`m`) to a rooted source erases the parity split: `X ↦ 2^m X`, `g` and `H` unchanged, `D ≡ 0`, and the complete
fixed-root continuation becomes the pair `(X, g)`. The DECODER obstruction needs `E ≠ O`, so it no longer applies.

Exhaustively at order 6 (two added vertices) and order 7 (one added vertex) there are sources whose three digits
`3H, 3H + 2, 3H + 4` share a common complete continuation `(E′, O′, g′)` — not the tripled one, a re-encoding — and
all of them have `H = 5` and a non-terminal root (receipt `anycommon.out`). The first:

| source masks (vertex 0 = low bit), root 3 | profile | children (two new vertices `x = 6, y = 7`; `x` loses to every old vertex and to `y`) | profiles |
|---|---|---|---|
| 62, 44, 56, 48, 34, 0 | (5, 6, 6, 3), `D = 0` | `y` beats {4} | (15, 40, 40, 9) |
| | | `y` beats {0, 1, 3, 4, 5} | (17, 40, 40, 9) |
| | | `y` beats {0, 1, 2, 3, 5} | (19, 40, 40, 9) |

Child masks: `254, 236, 248, 240, 98, 192, 0, 80` / `126, 108, 248, 112, 98, 64, 0, 123` / `126, 108, 120, 112, 226,
64, 0, 111`, root 3. The sink `x` is what makes the three continuations coincide: without it the three one-vertex
extensions `T + y` have `(X, g) = (40, 9)` but `D ∈ {−2, 0, 0}`. The same happens for the sources `(5, 12, 12, 3)` of
orders 6 and 7 (continuations `(40, 40, 9)` and `(80, 80, 9)`), and for the source `(5, 2, 4, 3)` with `D = −2` only
after the sink (continuation `(40, 9)` in the weak sense). So three consecutive odd counts can share a complete
fixed-root continuation once the root is non-terminal — the first such triple in the repository — at the price of
strongness (the children are ordinal sums; the digit lives in the prefix component).

## 4. The second step

From the three children (orders 8, root 3, continuation `(40, 40, 9)`) every one- and two-vertex extension was searched
for the next digits `3H′, 3H′ + 2, 3H′ + 4` (receipts `second_child{0,2,4}.out`):

| child | `H′` | `d′ = 0` hits | `d′ = 2` | `d′ = 4` | continuation shared by all three | by `+2` and `+4` |
|---|---|---|---|---|---|---|
| `d = 0` | 15 | 66 | 0 | 26 | none | none |
| `d = 2` | 17 | 60 | 40 | 0 | none | none |
| `d = 4` | 19 | 30 | 28 | 22 | none | none |

So the counts 45, 49 / 51, 53 / 57, 59, 61 are reached, but **no second fresh digit with a common continuation exists at
two added vertices**: the pre-registered P-4 is **refuted** in this range. The sink trick gives one synchronized
three-digit step, not yet an inductive digit system; the next thing to try is a larger prefix (a count-5 source is the
smallest possible) or a second sink after the second step, so that the parity split is erased again before each digit.

## 5. Assessment

- The parity split `D` of the fixed-root continuation is dependent in the literal families (`D = −g`) and vanishes for
  non-terminal roots; the honest continuation is `(H, X, g)` with `D` a module-determined quantity after a bridge step
  (BRIDGE_RESET: `D′ = (E_s − O_s) g` for retention, `(E_s − O_s) H` for reroot).
- A sink after the root converts the deposit's three-coordinate interface into two coordinates; three consecutive odd
  counts then share a complete continuation (count 5 → 15, 17, 19), which the DECODER obstruction forbade for `E ≠ O`.
- What remains for an inductive digit system is §4's question in general: whether the children's `(X, g)` is again a
  three-digit source. The path-splice theorem (SPLICE.md) says that gluing a source into a host at one vertex never
  loses the product, so the lower bound on the next count is always available; the exact digit needs the equality case.
