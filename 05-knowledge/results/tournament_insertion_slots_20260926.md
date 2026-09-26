# The insertion-slot count is exact, the insertion recursion is not: a tested reading of "T(N) is T(N−1) interspersed along a Hamiltonian path"

**Status: PROVED (one-line slot formula) + REFUTED (the recursion it
suggests) + FINITE-EXACT (exhaustive `n <= 6`, random `n = 7..10`). A probe
of the owner's seed, recorded so the decomposition is not re-proposed as an
identity. Session `collatz-exponent-atlas-20260926` (opus), 2026-09-26.**

Script: `04-computation/experiments/tournament_insertion_slots_20260926.py`.

## The seed

"A tournament of size `N` is the nodes of a complete graph of size `N-1`
converted to edges, strung together in a Hamiltonian path connecting `N`
nodes, and then the smaller graph's `T(N-1)` edge information interspersed
in the `T(N)` edges." Read as mathematics: build the Hamiltonian paths of
`T` on `N` vertices by inserting a new vertex `v` into the Hamiltonian
paths of `T - v`.

## What is true

**Slot formula (PROVED).** Let `P = x_1 ... x_(N-1)` be a Hamiltonian path
of `T - v` and `b_i = [v -> x_i]`. The positions at which `v` can be
inserted to give a Hamiltonian path of `T` are: in front iff `b_1 = 1`, at
the end iff `b_(N-1) = 0`, and between `x_i, x_(i+1)` iff `b_i = 0, b_(i+1) = 1`.
Since `#01 - #10 = b_(N-1) - b_1` (telescoping), the number of positions is

```text
slots_v(P) = 1 + #{ i : v -> x_i  and  x_(i+1) -> v }.
```

Verified with zero failures on all labelled tournaments with `n <= 6`
(every choice of `v`) and on random tournaments up to `n = 10`.

## What is false

**The recursion `H(T) = sum_(P in HP(T - v)) slots_v(P)` fails.** Deleting
`v` from a Hamiltonian path of `T` need not give a Hamiltonian path of
`T - v`: if `v` sits between `x` and `y` with `x -> v -> y` but `y -> x` in
`T`, the shortened sequence is not a path. Minimal witness: the cyclic
triangle `x_1 -> x_2 -> v -> x_1`. `T - v` has the single path `x_1 x_2`,
into which `v` has two slots (`v x_1 x_2` and `x_1 x_2 v`), but `H(T) = 3`:
the path `x_2 v x_1` is the one the recursion misses. Exhaustively, the
recursion is wrong for 3200 of the 5120 `(T, v)` pairs at `n = 5` and for
163,680 of 196,608 at `n = 6`; the missing paths are exactly those in which
`v` "repairs" a backward arc of `T - v`. The correction sum
`sum_P #{10 transitions}` is odd in a large fraction of cases, so it does
not carry Rédei's parity either.

## Reading

This vertex-insertion reading of the seed gives the slot formula: the
actual interior insertion gaps are `01` transitions, whereas one plus
the number of `10` transitions counts all slots, including endpoints.
The separate [marked-spine edge chart](crossroads233_20260926_graph.md)
realizes another exact reading of the seed. The repo's proved Hamiltonian-path machinery (THM-002
OCF, PROP-001 arc-flip identity, LEM-004 odd functions) works with arc
flips and with the cycle-count expansion, not with vertex insertion, and
this probe shows why: insertion is not a bijection between paths of `T`
and decorated paths of `T - v`. The metric-space analogy of the seed
("the triangle inequality is one relation between two points") combines
different types: exactly one arc per distinct pair is the tournament
axiom, whereas triangle inequality constrains distance values. A loop-free
complete weighted graph with edge lengths1,1,3 satisfies edge uniqueness
but violates triangle inequality. The seed's unique-relation portion adds
no condition to the tournament axiom; the metric equivalence is false.
