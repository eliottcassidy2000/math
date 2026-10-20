# The trijection: a tournament with a marked line, cycle-plus-solo, or star is a tournament on one vertex fewer, modulo 1, Z_{n-1}, or S_{n-1}

**kind-pasteur-2026-10-20-S3.** Owner-dispatched topic: tournaments as intransitivity and invariance; a
"trijection" between tournaments and simple graphs, read on the three graphs with 4 nodes and 3 edges
(the line `P4`, the cycle with a solo `K3+K1`, the star `K_{1,3}`), with the Hamiltonian path playing the
line, and its relation to the owner's `F = S + U` triplet `{p, p^3, p^2qr}`.

**Status.**
- **PROVED** (elementary; complete proofs below): Theorem 1 (the three marked bijections and their
  counts for every `n >= 4`; the line and star cases separately remain valid at `n = 3`), Theorem 2 (`n = 4` is the unique self-complementary size of the
  one-less-edge-than-node family), Theorem 3 (the `n = 4` dictionary: `H = 1 + 2c_3`, `c_3 + s = 2`,
  tetrahedral duality is an exact involution with `c_3 <-> s` and `H -> 6 - H`).
- **FINITE-EXACT**: all 38 checks of
  [`trijection_four_nodes_three_edges_kps_20261020.py`](../../04-computation/experiments/trijection_four_nodes_three_edges_kps_20261020.py)
  (`.out` alongside; complete enumeration of all labelled tournaments on 4, 5, 6 vertices, Burnside over
  the full symmetric group; the `n = 7` cycle count by direct Burnside over `Z_6`).
- **ANALOGY** (recorded with its hostile): the link to `F = S + U`. The only exact content shared with the
  owner's triplet principle is the shape "an involution with one fixed point, plus a factor-2 cocycle",
  which this note realises exactly on the graph side (§4). No map from divisors to tournaments is claimed.
- **NOT CLAIMED**: any Collatz statement. The Camion–Busch note already proved that no uniform
  isomorphism between Collatz and a tournament `Pi_1` sentence exists.
- **CONCURRENT WORK, first-pushed today and credited.** The same owner prompt ran in parallel.
  [trijection_reset_20261020](trijection_reset_20261020/README.md) (commit 44ff7e2d37) proves the
  line case of Theorem 1: marked pairs `(T, P)` `<->` simple graphs on the ordered set `[n-1]` `<->`
  tournaments on ordered `[n-1]`, `|M_n| = 2^C(n-1,2)`, with the fibre `H/|Aut|` from
  [LEM-003](../../01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md), the eight `n = 4` rows,
  the free carrier `13,14,24`, and the verdict that the three graph types, the three prime shapes and the
  three converse-merged tournament classes are distinct quotients. [tournament_graph_trijections_20261020](tournament_graph_trijections_20261020/README.md)
  (commit ae23c720da) treats observer charts and cut quotients, a different object, and inherits with
  [defect_boundary_20261020/web.md §4](defect_boundary_20261020/web.md) the `12/4/4` census and the
  distinction between complementation `C` and tetrahedral opposition `J`. This note was written before
  reading them and is kept for what they do not contain: the **cycle and star quotients** of Theorem 1
  with the closed form and the OEIS check, Theorem 2's `C(n-1,2) = n-1` reading and the index-2 chain,
  and Theorem 3's exact-involution law `H(T*) = 6 - H(T)`, `T** = T`, `c_3 + s = 2`. The line case and the
  `12/4/4` facts below are restated as inherited controls, not claimed.

## 0. Inheritance and board

- Closest proved mechanism: tetrahedral duality on 4-tournaments,
  [Camion–Busch note §6.3](camion_busch_gaps_polyhedra_collatz_20261002.md) (`c_3(T) = #sources + #sinks
  of T*`, transitive 24 `<->` strong 24, middle 16 fixed as a class); the fixed-path cube of
  [THM-1430 §1–3](../../01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md)
  with [LEM-003](../../01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md); the Triangle
  foundation ([everything-is-the-triangle](../../07-reflections/everything-is-the-triangle.md): a tournament
  with a linear order is a binary filling of the staircase `delta_{n-2}` above the superdiagonal); and the
  concurrent line-case packet [trijection_reset_20261020](trijection_reset_20261020/graph_audit.md).
- Canonical hostile: `P4` and `K3+K1` have the same independence polynomial `(1+x)(1+3x)`
  (Camion–Busch §1), so equal counts never certify equal structure; and MISTAKE-529 (a repeated numeral
  and a carrier switch were once read as a Lawvere bridge). Both are applied in §4.
- Corrected near miss: every triplet-principle row of Camion–Busch §7 is ANALOGY; this note adds one exact
  instance and keeps the arithmetic link at ANALOGY.
- Least-used sidecar: the opposite-edge map of the tetrahedron (the labelled form of duality on edge sets)
  agrees with complementation on every triangle and star and disagrees with it on every `P4` (§3.3).
- Board: marked structure `X`; its automorphism group `A_X`; the complementary pair set (always of size
  `C(n-1,2)`); the odd-order rule (a permutation fixes a tournament iff it has odd order); the
  orientation reversal (inner for path and cycle, outer for the star).

## 1. Theorem 1 — the three marked bijections

Fix `n >= 4`. A **coherent `(n-1)`-edge structure** in a tournament `T` on `n` vertices is one of

| kind | undirected carrier (`n-1` edges) | coherent orientation | automorphism group `A_X` of the oriented carrier |
|---|---|---|---|
| line | the path `P_n` | directed Hamiltonian path | trivial |
| cycle + solo | `C_{n-1} + K_1` | directed `(n-1)`-cycle, one vertex outside | `Z_{n-1}` (rotations) |
| star | `K_{1,n-1}` | all arcs out of the centre, or all into it | `S_{n-1}` (leaves) |

In each family exactly two of the `2^{n-1}` orientations of the carrier are coherent, and they are
converses of each other. At `n = 4` the carriers are the three graphs with 4 nodes and 3 edges: the line
`P4`, the cycle with a solo `K3+K1`, the star `K_{1,3}`.

**Theorem 1.** Let `lines(n)`, `cycles(n)`, `stars(n)` be the numbers of isomorphism types of pairs
`(T, X)` with `X` a coherent structure of the given kind. Then

```text
lines(n)  = 2^C(n-1,2)                                   = |Tour(n-1)|,
cycles(n) = (1/(n-1)) sum_{d | n-1, d odd} phi(d) 2^(C(n-1,2)/d) = |Tour(n-1) / Z_{n-1}|,
stars(n)  = 2 t(n-1)                                       = 2 |Tour(n-1) / S_{n-1}|,
```

where `Tour(m)` is the set of labelled tournaments on `m` vertices, `t(m)` its number of isomorphism
types (A000568), and the quotients are orbit sets of the rotation group and of the full symmetric group.
Each equality is realised by an explicit bijection.

*Proof.* The complement of each carrier in `K_n` has `C(n,2) - (n-1) = C(n-1,2)` pairs, and `T` is the
carrier's coherent orientation together with an arbitrary orientation of those pairs. `S_n` is transitive
on directed Hamiltonian paths of `K_n`, on directed `(n-1)`-cycles, and on out-stars (and separately on
in-stars), so marked types are orbits of the stabiliser `A_X` on the `2^C(n-1,2)` orientations of the
complement.

*Line (inherited: LEM-003 and trijection_reset_20261020, restated).* The stabiliser of a directed path is
trivial (an automorphism fixing a directed Hamiltonian path fixes every vertex), so `lines(n) = 2^C(n-1,2)`.
Writing the path as `v_1 -> ... -> v_n`, the complementary pairs `{v_i, v_j}` with `j >= i+2` correspond
bijectively to the pairs `{i, j-1}` of `[n-1]`; this is the staircase `delta_{n-2}` read as the upper
triangle of an `(n-1) x (n-1)` matrix, i.e. `Tour(n-1)`.

*Cycle + solo.* Put the directed cycle on `Z_{n-1}` as `i -> i+1` with solo `v`. Define `T'` on `Z_{n-1}`
by keeping every non-cycle pair of `T` and orienting `{i, i+1}` as `i -> i+1` iff `v -> i` in `T`. This is a
bijection from tournaments containing the cycle to `Tour(n-1)`, and it is `Z_{n-1}`-equivariant (rotation
moves `v -> i` to `v -> i+1` and `{i,i+1}` to `{i+1,i+2}`). Hence marked types are `Z_{n-1}`-orbits on
`Tour(n-1)`. Burnside: a rotation of order `d` fixes a tournament iff `d` is odd (an even-order rotation
has a power that swaps the ends of a chord); an odd-order rotation has all `C(n-1,2)/d` pair orbits of
size `d`, each oriented freely. There are `phi(d)` rotations of order `d`.

*Star.* A marked source `v` is the same as the tournament `T - v`, so marked out-stars are `t(n-1)` types,
and in-stars another `t(n-1)`; `S_n` does not identify an out-star with an in-star. `QED`

The three counts are twisted Burnside averages of the same function over the chain
`1 <= Z_{n-1} <= S_{n-1}` acting on the pairs of an `(n-1)`-set, with the odd-order rule
`Fix(g) = 2^{c_2(g)}` for odd-order `g` and `0` otherwise:

| n | lines = 2^C(n-1,2) | cycles (formula, Burnside) | stars = 2 t(n-1) |
|---:|---:|---:|---:|
| 4 | 8 | 4 | 4 |
| 5 | 64 | 16 | 8 |
| 6 | 1024 | 208 | 24 |
| 7 | 32768 | 5472 | 112 |
| 8 | 2^21 | 299600 | 912 |
| 9 | 2^28 | 2^25 | 13760 |

Rows 4, 5, 6 are verified by complete enumeration (checks D2–D4, E2); the cycle column is verified against
direct Burnside over `Z_{n-1}` for `n = 4..7` (E7). The identity `lines = cycles + stars` holds only at
`n = 4` (`8 = 4 + 4`); it fails at `n = 5` (`64 != 16 + 8`) and `n = 6` (E3). Its `n = 4` mechanism is the
accident `2(2^3 + 2 * 2) = 3 * 2^3`, together with the fact that the odd-order part of `S_3` is `Z_3` (E8), so
`stars(4) = cycles(4)`.

**Reading.** "One less edge than node" is exactly the condition that the mark consumes one vertex of
information: the complementary pair set of an `(n-1)`-edge carrier has the size of `Tour(n-1)`. The
trijection is the single statement "tournaments on `n` vertices with a marked coherent `(n-1)`-edge
structure = tournaments on `n-1` vertices modulo the symmetry of the mark", instantiated three times.

## 2. Theorem 2 — why 4 nodes and 3 edges

**Theorem 2.** `C(n-1,2) = n-1` iff `n = 4`. Consequently `n = 4` is the only order at which the
complement of a coherent `(n-1)`-edge carrier is again an `(n-1)`-edge carrier, and there complementation
acts on the three carriers as

```text
P4 -> P4,     K3+K1 <-> K_{1,3}.
```

This is an involution with exactly one fixed type. The tetrahedron is the self-dual polyhedron on these
vertices; the three carriers are the three edge-decompositions of `K4` into two 3-edge graphs:
`K4 = P4 + P4` (a Hamiltonian path and its anti-path, a second `P4`) and `K4 = (K3+K1) + K_{1,3}`
(a cyclic triangle and the star of its solo). Checks A1–A6.

*Proof.* `(n-1)(n-2)/2 = n-1` iff `n-2 = 2`. The three complements at `n = 4` are read off the 20 labelled
3-edge subgraphs of `K4` (12 paths, 4 triangles, 4 stars). `QED`

At `n = 4` the chain `1 < Z_3 < S_3` has index 2 at the top and `Z_3 = A_3`, which is why the three counts
collapse to `8, 4, 4`: the line keeps all `2^3` orientations of its anti-path, the cycle halves them by
rotation, and the star halves again by the index-2 extension to `S_3` and doubles back by its outer
orientation reversal (source versus sink). The reversal is inner for the path and for the cycle (a
vertex permutation reverses a directed path and reverses a directed cycle) and outer for the star.

## 3. Theorem 3 — the `n = 4` dictionary

Write `H` for the number of directed Hamiltonian paths, `c_3` for the number of cyclic triangles, and
`s` for the number of sources plus sinks. These count the coherent copies of the three carriers in `T`.

**Theorem 3.** For every tournament on 4 vertices:

1. `H = 1 + 2 c_3` and `c_3 + s = 2`; the four types carry `(H, c_3, s) = (1,0,2), (3,1,1), (3,1,1), (5,2,0)`.
2. Tetrahedral duality `T -> T*` (orient the dual edge of each arc from the face it runs along to the face
   it runs against, and identify faces with their opposite vertices) is an exact involution on the 64
   labelled tournaments, commutes with the converse, maps every arc to the opposite edge of the
   tetrahedron, and satisfies `c_3(T*) = s(T)`, `s(T*) = c_3(T)`, `H(T*) = 6 - H(T)`. On types it swaps
   `TT4 <-> strong` and preserves the class of the two cyclic-triangle types; the converse swaps those two
   and fixes `TT4` and `strong`. Duality and converse generate a Klein four-group on the four types.
3. On the 20 labelled 3-edge graphs the edge map of duality is the opposite-edge map `J` of
   [defect_boundary_20261020/web.md §4](defect_boundary_20261020/web.md), which already separates `J` from
   complementation `C`; the refinement here is quantitative: `J = C` on each of the 8 triangles and stars,
   and `J != C` on every one of the 12 paths, where `P4^opp` shares exactly two edges with `P4`.

*Proof.* (1) `c_3 = C(4,3) - sum_v C(s_v,2)` (Kendall–Babington Smith), and `C(x,2) - x + 1 = (x-1)(x-2)/2`
equals `[x in {0,3}]` for `0 <= x <= 3`, so `c_3 + s = 4 - sum_v (s_v - 1) = 4 - (6 - 4) = 2`. The odd-cycle
formula [THM-002](../../01-canon/theorems/THM-002-ocf.md) gives
`H = I(Omega(T), 2)`; at `n = 4` two cyclic triangles always share an arc, so `Omega` has no independent
pair and `H = 1 + 2 c_3`. (2) A face is a directed triangle iff all three dual arcs leave it or all three
enter it, i.e. iff the dual vertex is a source or sink (Camion–Busch §6.3); the involution and the
opposite-edge property are read off the construction and verified on all 64 (C1–C8). Then
`H(T*) = 1 + 2 s(T) = 1 + 2(2 - c_3(T)) = 6 - H(T)`. (3) A triangle or a star contains exactly one edge
from each opposite pair of the tetrahedron, so its opposite-edge image is its complement; a path contains
one opposite pair and misses one. `QED`

The line-and-staircase bijection at `n = 4` reads: with the path `0 -> 1 -> 2 -> 3`, the anti-path is
`2-0-3-1` and the eight upset graphs on it give `{} -> TT4`, the two end edges `{02}`, `{13}` give the two
cyclic-triangle types, and the middle edge `{03}` and all four multi-edge sets give the strong tournament,
which indeed has `H = 5` (E5). Complementation on the anti-path pairs the eight marked types as
`(TT4, strong)`, `(C3+source, strong)`, `(C3+sink, strong)`, `(strong, strong)` (E6); it is the involution
"reverse every non-path arc", which sends the transitive tournament to the anti-transitive one with scores
`(1,1,2,...,n-2,n-2)` for every `n`.

## 4. The triplet principle and `F = S + U`: exact content and hostile

The owner's triplet principle (TCPC, Camion–Busch §7) asks for "an involution with one fixed point plus a
factor-2 cocycle", and the owner's arithmetic triplet is the set of shapes with `F = S + U`
([arithmetic_braids_20260917_divisors.md, DB1](arithmetic_braids_20260917_divisors.md)): `F(N)` proper
nontrivial divisors, `S(N)` the squarefree ones, `U(N)` the primes below `N`; equality holds exactly for
`N = p, p^3, p^2qr` with values `0 = 0 + 0`, `2 = 1 + 1`, `10 = 7 + 3` (F1, F2 re-verified below 1000).

Exact content on the graph side (PROVED above):

- the involution with one fixed point is complementation (equivalently the opposite-edge map) on
  `{P4, K3+K1, K_{1,3}}`: `P4` fixed, the pair swapped; on tournaments it is tetrahedral duality,
  `TT4 <-> strong` with the cyclic-triangle class fixed;
- the factor-2 cocycle is the index-2 chain `Z_3 < S_3` with `Z_3 = A_3`, halving `4 -> 2`, compensated by
  the outer orientation reversal of the star, doubling `2 -> 4`; the line, with trivial symmetry, keeps 8.

Connection to `F = S + U` (ANALOGY):

- source: the shapes `{p} u {p^3, p^2qr}`, the squarefree solution and the two non-squarefree solutions of
  `prod(a_i+1) = 2^r + r + 1`; target: `{P4} u {K3+K1, K_{1,3}}`, the fixed point and the swapped pair;
- map: none. The shared datum is the partition shape `1 + 2`. Nothing identifies `p` with the line or
  `p^2qr` with the cycle-plus-solo beyond the owner's word "solo" for the repeated factor;
- preserved predicate: only the partition `1 + 2`;
- lost: every count. Hostile F3: the graph-side triple `(8, 4, 4)` matches no `(F, S, U)`, and `F` never
  takes the value `8`; the owner's "`p^2qr` is `p` and `p^3` combined" has no edge-arithmetic counterpart
  (`3 != 3 + 3`, and the carriers are not unions of one another);
- needed sidecar: an arithmetic statistic on `{p, p^3, p^2qr}` carrying an involution with fixed point `p`
  and a factor-2 cocycle between `p^3` and `p^2qr`. The divisor lattice offers none (Camion–Busch §7:
  no factor-2 cocycle between 3 and 5 at 60);
- cheapest decisive test: exhibit such a statistic, or stop. This note stops; the status stays ANALOGY,
  as in Camion–Busch §7 and as first pushed today by trijection_reset_20261020 ("distinct quotients;
  equality of their counts does not provide an incidence-preserving identification"). What this note adds
  to that verdict is the exact group-theoretic form of the "halving/doubling" on the graph side.

## 5. Open and next

- **Not in the OEIS (checked 2026-10-20, 400519 sequences; the five-term run `4, 16, 208, 5472, 299600`
  matches nothing).** `|Tour(m)/Z_m| = 1, 1, 4, 16, 208, 5472, 299600, 33554432, 7635498336` (`m = 1..9`)
  is the number of tournaments on `m` labelled vertices up to cyclic rotation, equivalently of tournaments
  on `m+1` vertices with a marked directed `m`-cycle, with the closed form
  `(1/m) sum_{d | m, d odd} phi(d) 2^(C(m,2)/d)`. The circulant-tournament sequences A002086 and A049288
  count rotation-invariant tournaments, a different object. Candidate for the OEIS submission pipeline.
- The same equivariant trick gives tournaments with a marked directed Hamiltonian cycle:
  `(1/n) sum_{d | n, d odd} phi(d) 2^((C(n,2)-n)/d)`.
- Tetrahedral duality acts on `H` by `H -> 6 - H` only because `K4` is the tetrahedron. For `n >= 5` the
  natural objects are orientations of polyhedral graphs, not tournaments; a wildcard, not pursued.
- META-PATTERNS card used: "a quotient owes a loss ledger" (the trijection is three quotients of `Tour(n-1)`
  with explicit stabilisers). Candidate counterindication recorded: `lines = cycles + stars` is an `n = 4`
  numeral, not a law.

## 6. Reproduction

```bash
python3 04-computation/experiments/trijection_four_nodes_three_edges_kps_20261020.py
```

38 checks, exact integers only, about two seconds; output in the sibling `.out`.
