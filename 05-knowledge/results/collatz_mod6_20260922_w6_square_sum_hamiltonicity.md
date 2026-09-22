# Square-sum Hamiltonicity: exact horizons 1..40, forced-edge obstructions, the 25-fold blow-up, and what the n >= 25 theorem is not an instance of

**Status: FINITE-EXACT: Hamiltonian path and cycle existence in the square-sum graph `Q_n` for every `1 <= n <= 40`, exact path counts for `n <= 36` and cycle counts for `n <= 40`, all agreeing with OEIS A071983 / A090460 / A071984 (S3); component and leaf tables to `40` (S1); the `n = 15` edge census at vertex `4` (S5). PROVED: the two general obstructions (three leaves; a cut set `S` with `c(Q - S) > |S| + 1`) and the one-round forced-edge certificates for `n = 18, 19, 20, 21, 22` and the two-round certificate for `n = 24` (S4); Gerbicz's 25-fold blow-up step (odd `n`, chain from `1` to `3` in `Q_n` gives a Hamiltonian cycle of `Q_{25n+12}` with the same ends) by a finite junction check, hence the family `(71*25^m-1)/2` (S6, verified explicitly for `887`, `22187`, `554687`). CITED: "every `k >= 25` has a chain, every `k >= 32` a Hamiltonian cycle" (R. Gerbicz, Mersenneforum, 2018-01-21, via OEIS A090461 comment and the archived thread; method summarised in S7, not verified by this lane). REFUTED: the pasted "parity desert 18..22", "N=14 three components join", "N >= 25 conjectured", "4 is the edge to avoid", the "Delta = 4 prime ladder 3,7,11,17", and the pasted Lean `square_sum_graph` (not loopless at vertex value `2`). SCOPE: no map from square-sum Hamiltonicity to Collatz beyond the shape "finite base + inductive extension", and that shape is exactly what the counterexample portrait shows cannot work sheet-blindly (S9).** Session collatz-mod6-20260922, wave 6, lane `square_sum_hamiltonicity`, script `04-computation/experiments/collatz_mod6_20260922_w6_square_sum_hamiltonicity.py`, output `collatz_mod6_20260922_w6_square_sum_hamiltonicity.out` (deterministic; identical under `python3 -O`).

## Inheritance and concept board

The object is `Q_n`: vertices `1..n`, `x ~ y` iff `x != y` and `x + y` is a perfect square (the `x != y` clause matters: `2x` is a square for `x = 2, 8, 18, 32`, S1). The closest proved mechanism in the repository is the swap-fixed operation-fibre closure of canon THM-2422 (the distinct-summand closure of `{2,3}` is `P minus {1,4,6}` with the exact dyadic law `M_t = 27*2^(t-4)+1`), with THM-2433 and THM-362 read for the "one vertex at a time" growth vocabulary; this lane's growth process (adjoin `n`, gain the edges `n ~ k^2 - n`) is a different graph, and the three components `{1,3,6,8,10}`, `{2,7,9}`, `{4,5,11,12}` of `Q_12` are not the summand module `{1,4,6}` (1 and 6 share a component, 4 sits in a third; no map found beyond the cardinality 3, S9 of the `.out`). The canonical hostile is the minus sheet of [counterexample_portrait](collatz_mod6_20260922_counterexample_portrait.md): any "finite residue base plus inductive extension" argument that never uses the sign of the carry is satisfied verbatim by `3n-1`, which has three cycles, so the argument shape borrowed from the square-sum theorem cannot be carried over sheet-blindly (S9 here). The corrected near miss is the pasted "parity desert": the failures at `18..22` and `24` are leaf and forced-edge contradictions, not parity phenomena (S4), and the "unification at `N = 14`" is really the `3 -> 2` merge at `13` (edges `13-3`, `13-12`) followed by the `2 -> 1` merge at `14` (edges `14-2`, `14-11`) (S1). The least-used sidecar is the "nice pair" invariant of Gerbicz's post #22 (a chain of length `n` and one of length `n+1` with every value in same-parity positions), which is the one genuinely inductive ingredient of the `k >= 25` theorem and has no Collatz analogue in the repository. Inherited by path and not re-derived: [summand-graph-fermat-zeckendorf](../../07-reflections/summand-graph-fermat-zeckendorf.md) (three modes `1` unary, `4` binary bridge, `6` ternary; its corrected errors are [arithmetic_braids_20260917_summand](arithmetic_braids_20260917_summand.md) section 8), [arithmetic_braids_20260917_collatz](arithmetic_braids_20260917_collatz.md), the [synthesis](collatz_mod6_20260917_synthesis.md) sections 1-12 and the twenty-two lane notes it indexes. The paste's `196 = 14^2` horizon, `2363 = 17*139` clock, octonion monodromy and "martingale decay per sheet" contain no square-sum content and are not addressed (SCOPE).

Session lead probes (all verified, S1-S2 of the `.out`): 3 components for `4 <= n <= 12`, 2 at `13`, connected for `14 <= n <= 40`; degree-`<= 1` vertices `{16,17,18}` at `18`, `{16,18}` at `19`, `{18}` for `20..30`, none at `31, 32`; the two quoted paths for `15` and `23` are valid; `18`'s only neighbour below `31` is `7`; the "ladder" `3, 7, 11, 17` has differences `4, 4, 6`.

## 1. Components and leaves, `n = 1..40` (S1, FINITE-EXACT)

```text
n   #comp  components (if >1)                     deg<=1 vertices
4   3      {1,3} {2} {4}                            [1, 2, 3, 4]
8   3      {1,3,6,8} {2,7} {4,5}                    [2, 4, 5, 6, 7, 8]
12  3      {1,3,6,8,10} {2,7,9} {4,5,11,12}         [2, 8, 9, 10, 11, 12]
13  2      {1,3,4,5,6,8,10,11,12,13} {2,7,9}        [2, 8, 9, 10, 11]
14  1      -                                        [8, 9, 10]
15  1      -                                        [8, 9]
16  1      -                                        [8, 16]
17  1      -                                        [16, 17]
18  1      -                                        [16, 17, 18]
19  1      -                                        [16, 18]
20..30  1  -                                        [18]
31..40  1  -                                        []
```

The full table is in the `.out`. The three components of `Q_12` are the residue-free classes generated by the squares `4, 9, 16` alone; `13` welds `{4,5,11,12}` to `{1,3,6,8,10}` through `16 = 13+3` and `25 = 13+12`, and `14` welds `{2,7,9}` through `16 = 14+2` and `25 = 14+11`.

## 2. Existence of Hamiltonian paths and cycles, `n = 1..40` (S2, FINITE-EXACT)

```text
n with a Hamiltonian path : [1, 15, 16, 17, 23, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
n without (2<=n<=40)      : [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 22, 24]
n with a Hamiltonian cycle: [32, 33, 34, 35, 36, 37, 38, 39, 40]
```

This is `{1}` union A090461 for the paths (A078107 lists `1` among the impossible `k` by its chain convention; `Q_1` has the trivial one-vertex path) and exactly the support of A071984 for the cycles. The search is a bitmask backtracking with three sound pruning rules (every unvisited vertex keeps an available neighbour; at most one unvisited vertex has exactly one; the unvisited set plus the current vertex is connected), tried from every start vertex in increasing degree order; witnesses are printed in the `.out` (for example `n = 32`: `1,8,28,21,4,32,17,19,30,6,3,13,12,24,25,11,5,31,18,7,29,20,16,9,27,22,14,2,23,26,10,15`, a cycle since `15 + 1 = 16`).

Two general obstructions, both elementary:

**S2a (PROVED). Three leaves.** If a graph on `n >= 2` vertices has three vertices of degree `1`, it has no Hamiltonian path. *Proof.* A Hamiltonian path `P` visits every vertex; a vertex `v` of degree `1` in the graph is joined in `P` only to its unique neighbour, so it has path-degree `1`, i.e. it is an endpoint of `P`. `P` has exactly two endpoints. QED

**S2b (PROVED). Cut criterion.** If a graph `G` has a Hamiltonian path then `c(G - S) <= |S| + 1` for every vertex set `S`, where `c` counts components. *Proof.* Deleting the `|S|` vertices of `S` from the path cuts it into at most `|S| + 1` nonempty sub-paths (one per gap between consecutive deleted vertices, plus the two ends), each lying inside one component of `G - S`; every component of `G - S` contains a vertex of the path and hence meets at least one sub-path. So `c(G - S) <=` number of sub-paths `<= |S| + 1`. QED

## 3. Exact counts versus OEIS (S3, FINITE-EXACT)

Paths are counted up to reversal (this is A071983, which counts each Hamiltonian cycle `n` times, once per cut edge); cycles up to rotation and reversal (A071984); the essentially-different count A090460 is `paths - (n-1)*cycles`.

```text
n   paths(=A071983)  cycles(=A071984)  A090460 = paths-(n-1)*cycles   OEIS A071983 A090460
15  1                0                 1                              1        1        OK
16  1                0                 1                              1        1        OK
17  1                0                 1                              1        1        OK
18  0                0                 0                              0        0        OK
19  0                0                 0                              0        0        OK
20  0                0                 0                              0        0        OK
21  0                0                 0                              0        0        OK
22  0                0                 0                              0        0        OK
23  3                0                 3                              3        3        OK
24  0                0                 0                              0        0        OK
25  10               0                 10                             10       10       OK
26  12               0                 12                             12       12       OK
27  35               0                 35                             35       35       OK
28  52               0                 52                             52       52       OK
29  19               0                 19                             19       19       OK
30  20               0                 20                             20       20       OK
31  349              0                 349                            349      349      OK
32  392              1                 361                            392      361      OK
33  669              1                 637                            669      637      OK
34  4041             11                3678                           4041     3678     OK
35  17175            57                15237                          17175    15237    OK
36  12960            31                11875                          12960    11875    OK
37  (not counted)    20                                               A071984=20   OK
38  (not counted)    25                                               A071984=25   OK
39  (not counted)    50                                               A071984=50   OK
40  (not counted)    64                                               A071984=64   OK
```

At `n = 32` the `392` paths are the `361` essentially different solutions plus the `31` cuts of the unique cycle. The counts confirm, independently of OEIS, that the unique chains at `15, 16, 17` are those of the session lead and that `Q_31` has `349` chains but no cycle (Dobbelaere's "less trivial" case in the A071984 comment; here decided by the exhaustive search).

## 4. Obstruction certificates for `n = 18, 19, 20, 21, 22, 24` (S4, PROVED)

The script runs the following necessary conditions to a fixed point and prints every step; each rule is a theorem about an arbitrary Hamiltonian path `P` of the current graph: (R1) a degree-1 vertex is an endpoint and its edge is in `P`; (R2) both edges of a degree-2 vertex are in `P`; (R3) no vertex has three edges in `P`; (R4) no three leaves (S2a); (R5) a vertex with two forced edges loses its other edges; (R6) forced edges form vertex-disjoint fragments, and a non-forced edge joining the two ends of a fragment on fewer than `n` vertices is deleted (it would close a short cycle). The printed certificates, in prose:

**n = 18 (R4).** Leaves `16` (only `9`), `17` (only `8`), `18` (only `7`): three endpoints, impossible.

**n = 19, 21, 22 (R3 at vertex 7).** `2` has neighbours `{7, 14}` and `9` has `{7, 16}` (degree `2` in all three graphs), so `2-7` and `7-9` are forced; `18` has the single neighbour `7`, so `7-18` is forced. Vertex `7` would carry three path edges. (The adjacency lists are printed in the `.out`; `21` and `22` add only `4-21`, `15-21`, `3-22`, `14-22`, which do not touch `2`, `9`, `18`.)

**n = 20 (R3 at vertex 5).** `4` has neighbours `{5, 12}`, `11` has `{5, 14}`, `20` has `{5, 16}`: the three forced edges `4-5`, `5-11`, `5-20` meet at `5`.

**n = 24 (two rounds, R5/R6 then R4).** Round 1 forces the `21` edges listed in the `.out` (all degree-2 vertices plus the leaf `18`). Saturations: `1` is saturated by `1-8`, `1-24`, so `1-3` and `1-15` are deleted; `5` by `5-11`, `5-20`, so `4-5` is deleted; `6` by `6-10`, `6-19`, so `3-6` is deleted; `7` by `7-9`, `7-18`, so `2-7` is deleted; `14` by `11-14`, `14-22`, so `2-14` is deleted. The forced fragment `1-8-17-19-6-10-15-21-4-24-12` has ends `4` and `12` and covers `11 < 24` vertices, so the edge `4-12` is deleted. Round 2: `2` now has the single neighbour `23`, `4` the single neighbour `21`, `18` the single neighbour `7`: three leaves, contradiction.

**Cut sets (S2b, exact search over `|S| <= 3`).** `n = 18`: `S = {7}`, `c(Q - S) = 3`; `n = 19`: `S = {7}`, `c = 3`; `n = 20`: `S = {5, 7}`, `c = 4`; `n = 21`: `S = {5, 7}`, `c = 4`; `n = 22`: `S = {5, 7, 14}`, `c = 5`; `n = 24`: no cut set with `|S| <= 3`, so for `24` the forced-edge certificate above is the obstruction and the cut criterion (at this size) is not.

**Control `n = 23`.** The same reduction stops at `OPEN` and a path exists; the leaf `18` is an endpoint in every one of the `3` chains (directed counts by start vertex `{2: 1, 9: 1, 18: 3, 22: 1}`), with the other end at `2`, `9` or `22`. The session lead's reading "23 succeeds because `18` can be an endpoint with `7`" is exact.

## 5. `n = 15` and the vertex `4`; the pasted table (S5, S8; FINITE-EXACT)

`4` has neighbours `{5, 12}` in `Q_15` (`4 + 21 = 25` needs `21 > 15`), so both its edges are forced. Exact census: `1` Hamiltonian path up to reversal; through `4-5`: `1`; through `4-12`: `1`; through both: `1`. The degree-2 vertices of `Q_15` are `[2, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15]` and the leaves are `[8, 9]`, which is why the chain is unique: eleven vertices have both edges forced and two are forced endpoints. "4 is the edge to avoid" is REFUTED; `4` is an interior vertex of the unique chain `8,1,15,10,6,3,13,12,4,5,11,14,2,7,9`. The squares occurring are `9, 16, 25` only (`4` cannot occur as a sum of two distinct positive integers `<= 15` other than `1 + 3`, which is used at `n = 16`'s neighbour structure, not here; the `.out` lists the squares used).

The pasted table, entry by entry (S8):

```text
[CORRECT (imprecise)]    N<=13 disconnected            -> 3 components for 4<=N<=12, 2 at N=13
[REFUTED as stated]      N=14 three components join    -> 3->2 at N=13 (13-3, 13-12), 2->1 at N=14 (14-2, 14-11)
[CORRECT]                N=15 first path, squares 9,16,25 -> unique chain (1)
[CORRECT (weak)]         N=16,17 connected             -> connected from 14; one chain each
[REFUTED (mechanism)]    N in [18,22] parity desert    -> leaves (18) and three-forced-edge vertices 7 / 5 (19..22)
[CORRECT, name unearned] N=23 valve                    -> 3 chains, all with 18 at one end
[CORRECT (data)]         N=24 stalled                  -> second-round leaves 2, 4, 18
[REFUTED]                N>=25 conjectured connected   -> connected from 14; PATHS for all N>=25 are a theorem (CITED)
[REFUTED]                4 is the edge to avoid        -> both edges at 4 forced
```

The user's phrase "sequences of length 15 to 25" is also imprecise: the lengths admitting a chain in that range are exactly `15, 16, 17, 23, 25`.

## 6. Gerbicz's 25-fold blow-up (S6): PROVED by a finite junction check

Let `a = (a_1, ..., a_n)` be a Hamiltonian path of `Q_n` with `n` odd, `a_1 = 1`, `a_n = 3` (then `a_1 + a_n = 4` makes it a cycle). For `-12 <= c <= 12` define the block `T(c) = (25 a_k + (-1)^(k+1) c)_k` and `R(T(c))` its reversal. Inside a block, consecutive sums are `25 (a_k + a_{k+1})`, squares. The `25` blocks for `c = -12..12` partition `[13, 25n + 12]` (distinct residues mod `25` because `25` is odd). The glue word of the `.out` (starting `1, T(-1), T(1), R(T(-7)), ...` and ending `..., R(T(8)), 3`) uses each `c` once, each singleton `1..12` once, and its `36` junction sums depend only on `a_1 = 1`, `a_n = 3` and the parity of `n`: a block `T(c)` has first element `25 + c` and last `75 + c`, and `R(T(c))` the reverse. The script lists all `36` junctions with their sums (`25, 100, 144, 49, 100, 144, 36, 81, 25, 9, 16, 25, 100, 169, 49, 144, 49, 81, 16, 36, 100, 81, 16, 36, 144, 16, 25, 81, 36, 144, 100, 49, 100, 121, 169, 36`), all squares, and the closing sum `1 + 3 = 4`.

**S6 (PROVED, construction CITED: R. Gerbicz, Mersenneforum post #19, 2018-01-17).** For every odd `n` and every Hamiltonian path of `Q_n` from `1` to `3`, the glued word is a Hamiltonian path of `Q_{25n+12}` from `1` to `3` (hence a Hamiltonian cycle), and `25n + 12` is odd. Iterating from the base chain of length `35` (`1,8,28,21,4,32,17,19,6,30,34,15,10,26,23,13,12,24,25,11,5,20,29,35,14,2,7,18,31,33,16,9,27,22,3`) gives Hamiltonian cycles of `Q_N` for `N = (71*25^m - 1)/2`, explicitly verified for `m = 1, 2, 3`: `N = 887, 22187, 554687`. *Proof.* The block sums are `25` times a square; the junction sums are the `36` constants above; the multiset of values is `{1..12}` union the `25` residue-shifted copies of `25 a_k`, i.e. `[1, 25n+12]`; the length is `12 + 25n`, odd when `n` is; the first and last letters are `1` and `3`. QED

## 7. Literature status (S7, S7b; CITED)

OEIS A090461 (T. D. Noe, Dec 01 2003) carries, verbatim in the `.out`, the comment chain: the original conjecture "sequence includes all integers k > 24"; Jud McCranie (Jan 11 2018) "25..299"; Robert Gerbicz (dated "Jan 17 2017" in the entry, evidently 2018 from the forum) the `2^20` verification plus the `(71*25^m-1)/2` family; and Robert Gerbicz (Jan 21 2018): the conjecture is proved, every `k >= 25` has a chain and every `k >= 32` a Hamiltonian cycle, with code and a deterministic algorithm on the Mersenneforum topic. A078107 (R. K. Guy) records via Paolo Xausa (May 29 2024) that it has no further terms. A071984's comment (Bert Dobbelaere, Dec 28 2018) gives the degree argument for no cycle at `n <= 30` and states `n = 31` "can be shown by hand". The live forum URL returned `404` on 2026-09-22; the thread is readable on web.archive.org, and its method (post #22, 2018-01-21) is, in this lane's words: blow a "nice pair" (chains of lengths `n` and `n+1` with every value in same-parity positions) up by `49` with shifts `c = -24..24` and glue with `1..24`, obtaining nice pairs for `49n + res`, `res = 24..72`, a complete residue system mod `49`; base table `n = 41..2032`; smaller `n` by lookup; `O(n log n)` time (`n = 10000000` in about `5` seconds). Post #9 (2018-01-12) describes the `2^20 = 1048576` verification by the flip-and-insert extension `F(i, j)` (fresh search last needed at `n = 6109`). The lane verified the post #19 step (S6) and nothing of post #22 beyond reading it: the glue tables live in the linked `squares.c`, which was not fetched. The theorem is therefore CITED, not verified.

## 8. Reproduction block

```text
cd <worktree>
python3 04-computation/experiments/collatz_mod6_20260922_w6_square_sum_hamiltonicity.py \
    > 05-knowledge/results/collatz_mod6_20260922_w6_square_sum_hamiltonicity.out
python3 -O 04-computation/experiments/collatz_mod6_20260922_w6_square_sum_hamiltonicity.py 2>/dev/null | diff - \
    05-knowledge/results/collatz_mod6_20260922_w6_square_sum_hamiltonicity.out   # identical
SQSUM_LIVE=1 python3 ...   # optional: re-fetch A090461 and compare its comments with the embedded copy
```

Pure Python, no dependencies, well under a minute and a few MB.

## 9. Microcosm to macrocosm: what the `n >= 25` theorem is, and the typed non-map to Collatz

**What it is.** Gerbicz's theorem has the shape *finite base* (a verified table of chains, in the final version `n = 41..2032`, with the small cases `25..40` by lookup) *plus an inductive extension* (an explicit embedding `Q_n x {shifts} -> Q_{49n+res}` that multiplies every value by an odd square, uses `25(a + b)` or `49(a + b)` being a square iff `a + b` is, and carries an invariant, the nice-pair parity, that the induction needs). The "fractal repetition in the macrocosm" the user senses is real and completely explicit: the `(71*25^m-1)/2` cycles of S6 are literally `25` scaled copies of the base cycle stitched with the twelve smallest integers, and the general theorem is the same picture with `49` copies and a two-chain invariant.

**Typed analogy (no map found).** Source: Hamiltonian chains in `Q_n`, theorem "finite base + `49`-fold inductive extension". Target: Collatz convergence on the plus sheet. Map: the argument shape "finite set of residue cylinders (base) + inductive extension to all `n`". Preserved: only the shape. Lost: the embedding. The square-sum extension is a constructive self-similarity of the *solution* (scaling a chain by an odd square keeps it a chain); Collatz orbits have no self-similar embedding of orbits into orbits, and the cylinder statements that do extend (Terras words, prefix-descent densities, the gate `q = 1`, the product identity) are proved for both signs, so [counterexample_portrait](collatz_mod6_20260922_counterexample_portrait.md) shows any such "base + cylinder induction" that never uses `n > 0` at `b = +1` is satisfied by the `3n-1` sheet with its three cycles. Sidecar: the nice-pair invariant, the one sign-specific-looking ingredient of the source (it pins positions mod `2`); its Collatz analogue would be an invariant of residue cylinders that is *not* preserved by `n -> -n`, and the repository has no candidate (the portrait lists the additive budget as the only sign-specific object and records that the three minus cycles exhaust it). Test: none available, because the map stops at the shape. The lane therefore records SCOPE: the square-sum theorem is a model of how a "residue base + induction" proof looks when the inductive step is an explicit embedding, and it is evidence *against* expecting such a proof for Collatz without a sign-specific embedding.

## Stopping boundary / next question

Everything decidable at this size is decided and matched to OEIS; the failures are certified; the `25`-fold step is proved and the `49`-fold theorem is cited. The lane stops at the citation boundary: the next cheap, honest step is to fetch `squares.c` (or reconstruct the `49`-fold glue by the same finite search the script uses for `25`) and turn the CITED post #22 induction into a PROVED-by-finite-check statement in the repository, exactly as S6 did for post #19; that would also decide whether the nice-pair invariant can be dropped by using two odd squares (`25` and `49`) instead of one square with two chains. On the Collatz side the lane has no next question: the only thing the square-sum theorem contributes is the shape of an argument the portrait already excludes sheet-blindly.
