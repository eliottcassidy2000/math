# The dead-end descent: a one-outward-edge system on tournaments whose escapes recur to the base case

Owner: *"When it comes to preserving the answer while losing what happens next the takeaway lesson is to arrange our
constructions so that they orient the information they lose in such a way that what does escape from each element in the
set is arranged so that it ultimately recurs back to the base inductive case, like is supposed by the collatz conjecture
and a graph with one outward edge per natural number node."* The repo's earlier reading of this directive is
`collatz_one_out_edge_systems_20260930.md`; here the system is built on tournaments.

## 1. The system

For a tournament `S` with a digit-0 vertex let `δ(S) = S − z`, `z` the digit-0 vertex (unique when `S` is strong of
order ≥ 4, Theorem 5 of [OUTDEG.md](OUTDEG.md); otherwise the one of least canonical label — any fixed rule gives a
functional graph, and every rule terminates). By Theorem 3, `z` is a sink or a sole exit with a dead-end port `t_0`, and

    H(S) = h(δS) + ε(S),     ε(S) = 0 if z is a sink,   ε(S) = B_{t_0}(δS) ≥ h(δS) if z is a sole exit

(THM-5263 (I2)). The **escape** `ε(S)` is what the step loses: the digit of the deleted vertex, a function of the smaller
tournament `δS` and the port alone. Iterating, `S = S_0, S_1 = δS_0, …, S_r` sterile (no digit-0 vertex) or a single
vertex,

    H(S) = h(S_r) + Σ_{i<r} ε(S_i),

every term computed on a tournament of smaller order: the count recurs to the sterile base with the escapes as the carried
digits. This is the owner's arrangement exactly — the lost information (which vertex, which port) is oriented into one
number per step that returns to the base case. Unlike Collatz, termination is trivial (the order drops by one per step),
so the whole content of the system is in (a) which tournaments are sterile and which counts they carry, and (b) which
values the escapes can take — (b) is the near-join / dead-end-append problem of [CONVERSE.md](CONVERSE.md).

## 2. Census

- Strong tournaments of orders 3–8 ([receipts/converse_descent.out](receipts/converse_descent.out) (4)): 6,404, of which
  6,062 are sterile; the longest descents (choosing, at each step, the digit-0 vertex whose subtree is deepest) take
  `{1: 32, 2: 159, 3: 12, 4: 3, 5: 9, 6: 28, 7: 99}` steps (P-4, "≤ 4 steps", refuted: 136 fertile strong tournaments of orders ≤ 8 admit
  descents of 5 to 7 steps).
- One longest chain at order 8 ([receipts/converse_census.out](receipts/converse_census.out) (C)):
  `H = 51` (delete the sole exit, escape 34) → order 7, `H = 17` (sink, 0) → order 6, `H = 17` (sole exit, 16) → order 5,
  `H = 1` (transitive; four sink deletions) → one vertex. So `51 = 1 + 16 + 34`: `17 = 1 + 2⁴·1` is the near-join over a
  single vertex with `k = 4`, and `51 = 17 + 2·17` is the sink trick `(2 + 1)·17` (CONVERSE.md (S)).
- Order 9 ([receipts/sterile9.out](receipts/sterile9.out), `sterile9.c`): 178,133 strong classes, 175,012 sterile; at
  every order 4–9 no strong class has more than one digit-0 vertex (Theorem 5 replayed). Strong floors (least count of a
  strong tournament) 3, 5, 9, 15, 25, 45, 75 at orders 3–9 — Busch's bound `5^k, 9·5^{k−1}, 3·5^{k−1}` for
  `n = 3k+1, 3k+2, 3k` — so a count `v` has strong realisers only at the orders `n` with `f(n) ≤ v`, and the census
  through order 9 is complete for every `v < 125`.

## 3. Three kinds of counts (FINITE-EXACT below 125, permanent where stated)

For an odd count `v < 125` the strong realisers all live at orders ≤ 9, so the sterile/fertile census decides:

- **Sterile-only** (`every` strong realiser is sterile): `11, 13, 47, 67` — THM-5263's permanent near-join exceptions,
  now read as the *sources* of the descent system: they can be bases but are never produced by a dead-end append.
- **Fertile-only** (`every` strong realiser has a digit-0 vertex): exactly `3, 5, 9, 17, 27` below 125 — permanent, since
  their strong realisers stop at orders 3, 4, 5, 6, 7 respectively (floors 5, 9, 15, 25, 45 at the next orders). These
  counts are never bases: `3 = C_3`, `5` and `9 = 1 + 2³` and `17 = 1 + 2⁴` are near-joins over a single vertex,
  `27 = 3 + 2³·3` a near-join over `C_3` (P-13 confirmed: 17 and 27 never acquire a sterile realiser).
- **Mixed** (both kinds of strong realiser): every other odd `v < 125`, i.e. every odd number in `[15, 123]` except
  `17, 21, 27, 47, 67`. The sterile strong counts of orders ≤ 9 are every odd number in `[11, 199]` except `17, 21, 27`.

So the base cases of the descent are the sterile tournaments, which carry every odd count below 125 except
`1, 3, 5, 7, 9, 17, 21, 27`; the sterile-only counts are exactly the near-join exceptions among them; and the counts that
can only be reached *through* the descent (fertile-only) are the five small near-joins `3, 5, 9, 17, 27`.

## 4. What the analogy is and is not (typed)

- **ANALOGY, exact in form.** One outward edge per fertile strong tournament (Theorem 5), terminal nodes = sterile
  tournaments, escape = digit of the deleted vertex, count = base count + escapes. Collatz: one outward edge per natural
  number, conjectured terminal node 1, escape = the parity lost by the forward map.
- **Not an analogy in difficulty.** The hard direction of Collatz (termination) is trivial here; the hard direction here
  (which escape values occur — Chan–Pak) has no Collatz counterpart. The one-out-edge reading of Collatz in
  `collatz_one_out_edge_systems_20260930.md` (functional-graph rigidity, THM-4523) is about the map's graph, not about
  the spectrum of a weight; nothing in this packet moves either conjecture.
- **What it does buy.** A canonical recursion for the count of every fertile strong tournament, and the three-kind
  taxonomy of §3: Chan–Pak for near-join values (HYP-9453) is the statement that the escapes of this system, over all
  bases, realise every odd number except the sterile-only ones and 1, 7, 21.
