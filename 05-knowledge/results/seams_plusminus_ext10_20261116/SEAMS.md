# Seams: placing the centre at a prescribed target ±1 by one port contact or one arc flip

Session klein-2026-11-16-seam (2026-11-16). Owner: *"Work to place the seam's center at the prescribed target ±1 …
Work to change the quotient or the port contacts, not just the suffix word."* Concurrent today:
`center_126_new_interface_20261116` (a two-arc switch 125 ↔ 127 stable under appends; its exact inverse selects the
centres of one pencil and "does not select an arbitrary center, e.g. 30012") and `ready_restarts_rounding_20261116`
(the windows 13/15 and 41/43). THM-5260 showed that suffix words (towers) force the midpoint 3H but never 3H + 2. Here
the centre is moved by the port contacts of one vertex, with the quadratic form of THM-5266 (A), and by single arc flips.

**Definitions.** A *contact seam* at centre `c` is a pair `T + v_X`, `T + v_{X ∪ {x}}` whose counts are `c ∓ 1` (either
order); its *digit* is `Δ_X(x)/2`, where `Δ_X(x) = H(T + v_{X∪x}) − H(T + v_X) = (B_x − e_x) − Σ_{y∈X} (M(y,x) + M(x,y))`.
An *arc-flip seam* is a pair `T`, `T^{a↔b}` (one arc reversed) with counts `c ∓ 1`.

## 1. The arc-flip lemma (PROVED) and why contacts and flips are one thing

    H(T^{a↔b}) = H(T) + M(b, a) − M(a, b),        Σ_{a→b ∈ T} M(a, b) = (n − 1) h.

*Proof.* The Hamiltonian paths of `T` through the arc `a → b` are exactly the pairs `(P_1 ending at a, P_2 starting at b)`
of nonempty paths covering `V`, `M(a, b)` of them; after the flip the paths through `b → a` are the `M(b, a)` pairs
`(P_1 ending at b, P_2 starting at a)` (neither path uses the arc between `a` and `b`, since `a` and `b` lie in different
paths); every other Hamiltonian path survives the flip. Each Hamiltonian path has `n − 1` arcs. ∎ Replay: 10,566 flips on
orders 3–7, 0 mismatches ([receipts/seams.out](receipts/seams.out)).

`T + v_{X∪x}` is `T + v_X` with the arc between `v` and `x` reversed, so **a contact seam over a base of order `n` is an
arc-flip seam at order `n + 1`**, and an arc-flip seam whose arc is incident to a vertex `v` is a contact seam over the
base minus `v`. The two censuses agree where their orders match (the misses below 2,000 of flips at orders ≤ 9 are those
of contacts from bases ≤ 8). The contact digits are the symmetric pair sums `M(y, x) + M(x, y)` and `B_x − e_x`; the flip
digits are the antisymmetric differences `M(b, a) − M(a, b)`; both are read off the pair matrix.

## 2. Census ([receipts/seams.out](receipts/seams.out) bases ≤ 8, Python by the quadratic form; [receipts/seams9.out](receipts/seams9.out) bases ≤ 9, `seams9.c`, 21 s)

| bases | seam pairs | even centres covered | missed |
|---|---|---|---|
| orders ≤ 8 | 53,690 | 979 of 999 in [4, 2000] | 6, 8, 20, 22, 36, 40, 54, 56, 60, 64, 84, 94, 108, 148, 160, 166, 188, 196, 262, 312 |
| orders ≤ 9 | 824,248 | 4,982 of 4,999 in [4, 10000] | 6, 8, 20, 22, 36, 40, 54, 56, 60, 64, 84, 94, 108, 148, 188, 196, 262 |

- **An arbitrary centre is placed by one contact.** From bases of order ≤ 9 every even number in `[264, 10000]` is the
  centre of a one-contact seam; the 17 misses all lie below 263, four of them (6, 8, 20, 22) because a neighbour is a
  hole. Beyond 10,000 the bases of order ≤ 9 thin out (their contact values are counts of order-10 tournaments); a centre
  such as 30012 needs bases of order 11–12 by the same rule, which is the natural continuation of the concurrent
  packet's pencil inverses (those select only their own centres).
- **Residues.** The misses are `≡ 0, 4, 6 (mod 8)` (4, 8, 5 of them) and none is `≡ 2 (mod 8)`; the same held for the
  20 misses from bases ≤ 8. Pre-registered (S-8) and confirmed; no mechanism found — HYP-9863 (b).
- **The quotient barely moves.** Among the 980 first-witness seams at centres ≤ 2,000 (bases ≤ 8) only 2 join
  tournaments with different numbers of strong components: seams live inside the strong class, and changing the
  condensation is neither needed nor typical (S-4 refuted).
- **The +2 window.** For the odd counts `H ≤ 1000` realised as contact values, `(3H + 2, 3H + 4)` is a contact seam at the
  centre `3(H + 1)` for 492 of 498 (bases ≤ 9); the failures are `H = 1, 11, 17, 19, 27, 35`, whose centres 6, 36, 54, 60,
  84, 108 are exactly the missed centres divisible by 3. Where towers never reach `3H + 2` (THM-5260), one port contact
  gives both `3H + 2` and `3H + 4` for 98.8 % of the counts below 1,000; the exceptions are again the small hard values
  (11 sterile-only, 17 and 27 fertile-only, 35 an atom; 19).
- **Digits.** The contact digits `|Δ|/2` take every value in `[0, 3638]` (bases ≤ 9): every multiplier below 3,639 is a
  "±1 with a multiplier" of some contact flip. Digit 1 (a seam) is 0.18 % of all contact flips, digit 0 (a
  count-preserving contact change) 0.51 %. Arc-flip digits take every value in `[−153, 153]` at orders ≤ 8; 3.3 % of
  arcs are balanced (`M(a, b) = M(b, a)`, the count survives the flip) and 2.4 % carry a seam.

## 3. Three ways to move a count

| move | what changes | law | reach |
|---|---|---|---|
| suffix word (THM-5260, HYP-9453) | `k` sinks before the sole exit | `h + 2^k B`: a pencil in `k` | its own pencil; never `3H + 2` |
| quotient (ordinal sum, THM-5264 (E)) | the condensation | counts multiply | multiples only |
| port contacts (here) | the out-set `X` of one vertex | quadratic form; one contact = one arc flip | every even centre in `[264, 10000]` from bases ≤ 9 |

Chan–Pak OPEN; the seam census is finite-exact only.
