# Flexible cores composed with the one-port alphabet: invisible flips, delayed choices, permanent pairs

Session klein-2026-11-16-compose (2026-11-16). Owner: *"flexible cores be composed with the one-port continuation
alphabet"*, and *"hold both counts fixed while storing choices that only show up after another wrap"*. The alphabet is the
concurrent `one_port_neutral_continuations_20261117`'s: at a retained root `r` of `T`, with `g = H(T − r)`, `E`/`O` the
weighted even/odd root-start path sums (weight = Hamiltonian count of the complement), **S** appends a sink and **N** a
vertex beating only `r`, with exact transfers `S: (H, g, E, O) ↦ (H, g, E + O, E + O)` and
`N: (H, g, E, O) ↦ (H + 2E, g, E + O − g, E + O)`. In this programme's coordinates `E + O = B_r` (the root response,
THM-4865) and `B_r − e_r = 2E_r` (THM-5260), so the full rooted state is `L3 = (H, g, e_r, B_r)`, with `L2 = (H, g, B)` and
`L1 = (H, g)`. The flexible cores are the fibres with nontrivial count-preserving flips (THM-5268).

## 1. What a count-preserving flip stores ([receipts/compose.out](receipts/compose.out) (1))

Over the 4,378 nontrivial balanced flips of orders ≤ 8 (THM-5268) and every root `r` not on the flipped arc — 18,880
(flip, root) pairs, all with `ΔH = 0` by balance — the change of the root state splits as:

| `Δg` | `ΔB` | `Δe` | pairs | what the alphabet sees |
|---|---|---|---|---|
| 0 | 0 | 0 | 1,040 | nothing, ever: the flip is invisible at this root |
| 0 | 0 | ≠ 0 | 6 | after one **N** (`H + B − e`): a stored choice read by the first wrap |
| ≠ 0 | 0 | 0 | 216 | after two wraps only: `S` and the first `N` give the same counts, the second `N` reads `g` through `E′ = B − g` |
| 0 | ≠ 0 | 0 | 1,092 | after `S·N` or `N` |
| other | | | 16,526 | immediately or after one letter |

At the two endpoints of the flipped arc (6,436 pairs) `B` always changes (the tail loses the paths through the arc, the
head gains them). The 216 **delayed choices** are the owner's "choices that only show up after another wrap" realised
by a change of the core: same `H`, same `e_r`, same `B_r`, different deletion count `g` — invisible to `S`, invisible to
the first `N`, read by the second. The six first-wrap choices all have `Δe = ±2`; after one `N` they are permanent
`+2` pairs, and all three underlying flips (classes of counts 321, 253, 291 at order 8) land on the same pair
`751/753` (centre 752).

## 2. Permanent pairs from the rooted-state census ([receipts/compose.out](receipts/compose.out) (2))

All 52,313 distinct rooted states `(H, g, e, B)` of orders ≤ 8. Two kinds of permanent `+2` pairs (same continuation
under every `S/N` word):
- **K1**: same `(g, e, B)`, counts `H` and `H + 2` — the concurrent packet's rooted pairs (13/15 at centre 14, 41/43 at
  42 are here); 201 centres below 1,000, among them 4, 10, 14, 16, 18, 24, 26, 28, 30, 32, 34, 38, 40, 42, 44, 50, 52,
  54, 58, 60.
- **K2**: same `(H, g, B)`, end-counts `e` and `e + 2` — permanent after one `N`, at centre `H + B − e − 1`; 145 centres,
  95 below 1,000, none below 62.
Together they cover 250 of the 499 even centres in `[4, 1000]` (C-2's 60 % refuted); the misses begin 6, 8, 12, 20, 22,
36, 46, 48, 56, 62, 64, 68, 84, 86, 88 — the hole neighbours again, and 12 (11/13 is not a rooted pair at orders ≤ 8)
and 126 (the concurrent 125/127 switch keeps `(E, O − g)`, which is `N`-permanent but not `S`-permanent, so it is not a
K1 pair). One-contact seams (THM-5267) reach far more centres (979 of 999 below 2,000) but are not permanent; permanence
costs coverage.

## 3. Flexible cores feed the alphabet ([receipts/compose.out](receipts/compose.out) (3), (4))

The number of distinct rooted states `(g, e, B)` a fibre offers correlates with its number of nontrivial neutral edges
(Spearman 0.673 over the 61 odd `h ≤ 125`; C-5 confirmed): the rigid fibres 11, 13, 47 offer 26, 26, 30 states, the
flexible 37 offers 117. Iterating **N** at the newest vertex is the sole-exit tower, tribonacci from the third step
(THM-4865; 1,254 checks, 0 failures): the alphabet reapplied to its own output is the tower, and a flexible core gives
it more starting states than a rigid one of the same count.

## 4. Verdict

Composition works as a selection of starting states, not as new coverage: a nontrivial flip changes the root state at
most roots (16,526 of 18,880 immediately visible), leaves it untouched at 1,040, and stores a delayed choice at 216 —
the deletion count `g`, read only by the second wrap. Permanent pairs exist at half the even centres below 1,000 from
orders ≤ 8. Chan–Pak OPEN.
