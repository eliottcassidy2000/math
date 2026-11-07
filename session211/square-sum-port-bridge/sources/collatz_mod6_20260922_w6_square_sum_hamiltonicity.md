# Square-sum Hamiltonicity: exact horizons 1..40, forced-edge obstructions, the 25-fold blow-up, and what the n >= 25 theorem is not an instance of

**Status: FINITE-EXACT: Hamiltonian path and cycle existence in the square-sum graph `Q_n` for every `1 <= n <= 40`, exact path counts for `n <= 36` and cycle counts for `n <= 40`, all agreeing with OEIS A071983 / A090460 / A071984 (S3); component and leaf tables to `40` (S1); the `n = 15` edge census at vertex `4` (S5). PROVED (classical, not new): the two general obstructions (three leaves; a cut set `S` with `c(Q - S) > |S| + 1`, the path form of the standard toughness-type necessary condition). PROVED (new to the repository, added at the audit, A7): `Q_n` is connected for every `n >= 14`, because a square lies in `[n+1, 2n-1]` for every `n >= 5`. PROVED: leaf/certified-nonendpoint forcing at `n = 18, 19`, and cut certificates for every `n = 18..22` (S4). RETRACTED (2026-09-25): the unguarded degree-2 forcing rule and the resulting claimed forcing certificates for `n = 20, 21, 22, 24`; nonexistence at `24` retains two independent exhaustive searches and now has an independently audited endpoint-aware proof (PROVED; [decoder note](decoder_prime_square_20260925.md), section 4). PROVED: Gerbicz's 25-fold blow-up step (odd `n`, chain from `1` to `3` in `Q_n` gives a Hamiltonian cycle of `Q_{25n+12}` with the same ends) by a finite junction check, hence the family `(71*25^m-1)/2` (S6, verified explicitly for `887`, `22187`, `554687`). CITED: "every `k >= 25` has a chain, every `k >= 32` a Hamiltonian cycle" (R. Gerbicz, Mersenneforum, 2018-01-21, via OEIS A090461 comment and the archived thread; method summarised in S7, not verified by this lane). REFUTED: the pasted "parity desert 18..22", "N=14 three components join", "N >= 25 conjectured", "4 is the edge to avoid", the "Delta = 4 prime ladder 3,7,11,17", and the pasted Lean `square_sum_graph` (not loopless at vertex value `2`). SCOPE: no map from square-sum Hamiltonicity to Collatz beyond the shape "finite base + inductive extension", and that shape is exactly what the counterexample portrait shows cannot work sheet-blindly (S9).** Adversarially audited 2026-09-22 ([audit script](../../04-computation/experiments/collatz_mod6_20260922_w6_square_sum_hamiltonicity_audit.py), [audit output](collatz_mod6_20260922_w6_square_sum_hamiltonicity_audit.out)): every key number recomputed by independent code (union-find, weak-prune DFS, brute cut sets, Gerbicz's own PARI pseudo-code) and confirmed; one presentation error corrected (the `n = 24` forced fragment was printed in DFS order, which is not a path), two labels sharpened, the connectivity theorem added; edits listed in section 10. Session collatz-mod6-20260922, wave 6, lane `square_sum_hamiltonicity`, script `04-computation/experiments/collatz_mod6_20260922_w6_square_sum_hamiltonicity.py`, output `collatz_mod6_20260922_w6_square_sum_hamiltonicity.out` (deterministic; identical under `python3 -O`).

## Inheritance and concept board

The object is `Q_n`: vertices `1..n`, `x ~ y` iff `x != y` and `x + y` is a perfect square (the `x != y` clause matters: `2x` is a square for `x = 2, 8, 18, 32`, S1). The closest proved mechanism in the repository is the swap-fixed operation-fibre closure of canon THM-2422 (the distinct-summand closure of `{2,3}` is `P minus {1,4,6}` with the exact dyadic law `M_t = 27*2^(t-4)+1`), with THM-2433 and THM-362 read for the "one vertex at a time" growth vocabulary; this lane's growth process (adjoin `n`, gain the edges `n ~ k^2 - n`) is a different graph, and the three components `{1,3,6,8,10}`, `{2,7,9}`, `{4,5,11,12}` of `Q_12` are not the summand module `{1,4,6}` (1 and 6 share a component, 4 sits in a third; no map found beyond the cardinality 3, S9 of the `.out`). The canonical hostile is the minus sheet of [counterexample_portrait](collatz_mod6_20260922_counterexample_portrait.md): any "finite residue base plus inductive extension" argument that never uses the sign of the carry is satisfied verbatim by `3n-1`, which has three cycles, so the argument shape borrowed from the square-sum theorem cannot be carried over sheet-blindly (S9 here). The corrected near miss is the pasted "parity desert": failures at `18..22` have leaf/cut certificates; the former degree-2 forcing certificate at `24` was retracted on 2026-09-25 because a degree-2 vertex may be an endpoint (S4), and the "unification at `N = 14`" is really the `3 -> 2` merge at `13` (edges `13-3`, `13-12`) followed by the `2 -> 1` merge at `14` (edges `14-2`, `14-11`) (S1). The least-used sidecar is the "nice pair" invariant of Gerbicz's post #22 (a chain of length `n` and one of length `n+1` with every value in same-parity positions), which is the one genuinely inductive ingredient of the `k >= 25` theorem and has no Collatz analogue in the repository. Inherited by path and not re-derived: [summand-graph-fermat-zeckendorf](../../07-reflections/summand-graph-fermat-zeckendorf.md) (three modes `1` unary, `4` binary bridge, `6` ternary; its corrected errors are [arithmetic_braids_20260917_summand](arithmetic_braids_20260917_summand.md) section 8), [arithmetic_braids_20260917_collatz](arithmetic_braids_20260917_collatz.md), the [synthesis](collatz_mod6_20260917_synthesis.md) sections 1-12 and the twenty-two lane notes it indexes. The paste's `196 = 14^2` horizon, `2363 = 17*139` clock, octonion monodromy and "martingale decay per sheet" contain no square-sum content and are not addressed (SCOPE).

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

The full table is in the `.out`. The three components of `Q_12` are the classes generated by the squares `4, 9, 16` alone (the largest edge sum in `Q_12` is `11 + 12 = 23 < 25`); `13` welds `{4,5,11,12}` to `{1,3,6,8,10}` through `16 = 13+3` and `25 = 13+12`, and `14` welds `{2,7,9}` through `16 = 14+2` and `25 = 14+11`.

**S1a (PROVED, elementary; added at the audit, A7 of the audit `.out`).** `Q_n` is connected for every `n >= 14`. *Proof.* For `n >= 5` there is a square in `[n+1, 2n-1]`: for `n = 5..9` it is `9, 9, 9, 9, 16`; for `n >= 10`, `(sqrt(n+1)+1)^2 = n + 2 + 2 sqrt(n+1) <= 2n - 1` is equivalent to `4(n+1) <= (n-3)^2`, i.e. `n^2 - 10n + 5 >= 0`, which holds for `n >= 10` (the polynomial is `-4` at `9` and `5` at `10`), so the smallest square above `n` is at most `2n - 1` (checked numerically to `200000`). Hence the vertex `n` is adjacent to `k^2 - n in [1, n-1]`, and by induction from the connected `Q_14` every `Q_n`, `n >= 14`, is connected. QED. The pasted "`N >= 25` conjectured connected" is therefore a theorem for all `N >= 14`, and the interesting object from `14` on is Hamiltonicity, not connectivity.

## 2. Existence of Hamiltonian paths and cycles, `n = 1..40` (S2, FINITE-EXACT)

```text
n with a Hamiltonian path : [1, 15, 16, 17, 23, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
n without (2<=n<=40)      : [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 22, 24]
n with a Hamiltonian cycle: [32, 33, 34, 35, 36, 37, 38, 39, 40]
```

This is `{1}` union A090461 for the paths (A078107 lists `1` among the impossible `k` by its chain convention; `Q_1` has the trivial one-vertex path) and exactly the support of A071984 for the cycles. The search is a bitmask backtracking with three sound pruning rules (every unvisited vertex keeps an available neighbour; at most one unvisited vertex has exactly one; the unvisited set plus the current vertex is connected), tried from every start vertex in increasing degree order; witnesses are printed in the `.out` (for example `n = 32`: `1,8,28,21,4,32,17,19,30,6,3,13,12,24,25,11,5,31,18,7,29,20,16,9,27,22,14,2,23,26,10,15`, a cycle since `15 + 1 = 16`).

Two general obstructions, both elementary and classical (textbook facts, recorded here only because the certificates of section 4 use them; the cut criterion is the path form of the standard toughness-type necessary condition for Hamiltonian cycles, whose attribution to Chvátal 1973 is UNCITED-RECOLLECTION):

**S2a (PROVED, classical). Three leaves.** If a graph on `n >= 2` vertices has three vertices of degree `1`, it has no Hamiltonian path. *Proof.* A Hamiltonian path `P` visits every vertex; a vertex `v` of degree `1` in the graph is joined in `P` only to its unique neighbour, so it has path-degree `1`, i.e. it is an endpoint of `P`. `P` has exactly two endpoints. QED

**S2b (PROVED, classical). Cut criterion.** If a graph `G` has a Hamiltonian path then `c(G - S) <= |S| + 1` for every vertex set `S`, where `c` counts components. *Proof.* Deleting the `|S|` vertices of `S` from the path cuts it into at most `|S| + 1` nonempty sub-paths (one per gap between consecutive deleted vertices, plus the two ends), each lying inside one component of `G - S`; every component of `G - S` contains a vertex of the path and hence meets at least one sub-path. So `c(G - S) <=` number of sub-paths `<= |S| + 1`. QED

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

## 4. Obstruction certificates and endpoint correction (2026-09-25)

**Correction lineage.** The original rule R2 said that every degree-2 vertex has both edges in an arbitrary Hamiltonian path. This is false: it may be an endpoint. The 2026-09-22 audit independently reproduced the same false implication, so its earlier certificate confirmation is superseded here. The path/cycle census used separate exhaustive searches and survives unchanged.

**Minimal hostile.** A triangle has Hamiltonian paths; each endpoint has graph-degree 2 but uses only one incident edge. The same error is visible inside this lane's positive control `Q_23`: the valid path

    18,7,9,16,20,5,11,14,2,23,13,12,4,21,15,10,6,19,17,8,1,3,22

ends at `22`, whose graph neighbours are `{3,14}`, but omits `22-14`.

**Repaired forcing rule (PROVED, classical).** Both edges at a degree-2 vertex are forced only when it is certified not to be an endpoint. Two distinct leaves certify the endpoints and therefore make every other vertex an interior vertex. The corrected script uses exactly this certificate. Degree-1 forcing, the maximum path degree of 2, three-leaf contradiction, saturation deletion, and deletion of an edge that would close a short forced cycle remain sound. With fewer than two leaves, unclassified degree-2 vertices are not forced.

**n = 18.** Leaves `16` (only `9`), `17` (only `8`), `18` (only `7`): three endpoints, impossible.

**n = 19.** Leaves `16` and `18` fix the endpoints. Vertices `2` and `9` are therefore interior, with neighbours `{7,14}` and `{7,16}`; their edges to `7` are forced. Together with `18-7`, this gives three forced edges at `7`, impossible. This certificate survives the correction.

**Cut certificates (PROVED; exact search over `|S| <= 3`).** By S2b: `n = 18,19` have `S = {7}`, `c(Q-S) = 3`; `n = 20,21` have `S = {5,7}`, `c = 4`; `n = 22` has `S = {5,7,14}`, `c = 5`. These independently certify all failures at `18..22`. The earlier one-round forcing arguments for `20..22` are retracted because they did not certify the second endpoint.

**n = 24 (PROVED by a new endpoint-aware certificate; former forcing proof RETRACTED).** There is only one original leaf, `18`. The earlier claim of 21 forced edges in round 1, the ensuing seven deletions, and the second-round three-leaf contradiction relied on false R2. They do not constitute a proof. The corrected reduction forces only `7-18` and returns `OPEN`. No cut with `|S| <= 3` was found. Nonexistence is still established in the explicit finite universe by the independent exhaustive searches S2 and A2, both giving zero paths. The later [endpoint-aware proof](decoder_prime_square_20260925.md), section 4, closes the structural slot: the other endpoint must lie in `{2,9,11,22}`; all four cases saturate vertex 5 and force the proper cycle `1-8-17-19-6-10-15-21-4-12-24-1`. Every degree-two use certifies its vertex as interior first. The general reducer still returns `OPEN`; the new endpoint case analysis is implemented in the decoder note's separate script.

**Control n = 23.** Corrected reduction returns `OPEN` and never forces `14-22`. The leaf `18` is an endpoint in all three chains; directed counts by start are `{2:1,9:1,18:3,22:1}`. The triangle control now forces no edges. These controls specifically exercise the missing endpoint coordinate.

## 5. `n = 15` and the vertex `4`; the pasted table (S5, S8; FINITE-EXACT)

`4` has neighbours `{5, 12}` in `Q_15` (`4 + 21 = 25` needs `21 > 15`); the leaves `8` and `9` fix both endpoints, so `4` is interior and both its edges are forced. Exact census: `1` Hamiltonian path up to reversal; through `4-5`: `1`; through `4-12`: `1`; through both: `1`. The degree-2 vertices of `Q_15` are `[2, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15]` and the leaves are `[8, 9]`, which is why the chain is unique: eleven vertices have both edges forced and two are forced endpoints. "4 is the edge to avoid" is REFUTED; `4` is an interior vertex of the unique chain `8,1,15,10,6,3,13,12,4,5,11,14,2,7,9`. The squares occurring along the chain are `[9, 16, 25]`; the edge `1-3` with sum `4` exists in `Q_15` but is unused (the leaf `8` forces `1-8` and the degree-2 vertex `15` forces `1-15`, saturating `1`).

The pasted table, entry by entry (S8):

```text
[CORRECT for 2<=N<=13]   N<=13 disconnected            -> 3 components for 4<=N<=12, 2 at N=13; Q_1 is connected (fails at N=1)
[REFUTED as stated]      N=14 three components join    -> 3->2 at N=13 (13-3, 13-12), 2->1 at N=14 (14-2, 14-11)
[CORRECT]                N=15 first path, squares 9,16,25 -> unique chain (1)
[CORRECT (weak)]         N=16,17 connected             -> connected from 14; one chain each
[REFUTED (mechanism)]    N in [18,22] parity desert    -> leaf/certified forcing (18,19) and cut certificates (18..22)
[CORRECT, name unearned] N=23 valve                    -> 3 chains, all with 18 at one end
[CORRECT (data)]         N=24 stalled                  -> new endpoint-aware proof plus exhaustive searches; old forcing proof RETRACTED
[REFUTED]                N>=25 conjectured connected   -> connected for ALL N>=14 (PROVED, S1a); PATHS for all N>=25 are a theorem (CITED)
[REFUTED]                4 is the edge to avoid        -> both edges at 4 forced
```

The user's phrase "sequences of length 15 to 25" is also imprecise: the lengths admitting a chain in that range are exactly `15, 16, 17, 23, 25`.

## 6. Gerbicz's 25-fold blow-up (S6): PROVED by a finite junction check

Let `a = (a_1, ..., a_n)` be a Hamiltonian path of `Q_n` with `n` odd, `a_1 = 1`, `a_n = 3` (then `a_1 + a_n = 4` makes it a cycle). For `-12 <= c <= 12` define the block `T(c) = (25 a_k + (-1)^(k+1) c)_k` and `R(T(c))` its reversal. Inside a block, consecutive sums are `25 (a_k + a_{k+1})`, squares. The `25` blocks for `c = -12..12` partition `[13, 25n + 12]` (distinct residues mod `25` because `25` is odd). The glue word of the `.out` (starting `1, T(-1), T(1), R(T(-7)), ...` and ending `..., R(T(8)), 3`) uses each `c` once, each singleton `1..12` once, and its `36` junction sums depend only on `a_1 = 1`, `a_n = 3` and the parity of `n`: a block `T(c)` has first element `25 + c` and last `75 + c`, and `R(T(c))` the reverse. The script lists all `36` junctions with their sums (`25, 100, 144, 49, 100, 144, 36, 81, 25, 9, 16, 25, 100, 169, 49, 144, 49, 81, 16, 36, 100, 81, 16, 36, 144, 16, 25, 81, 36, 144, 100, 49, 100, 121, 169, 36`), all squares, and the closing sum `1 + 3 = 4`.

**S6 (PROVED, construction CITED: R. Gerbicz, Mersenneforum post #19, 2018-01-17).** For every odd `n` and every Hamiltonian path of `Q_n` from `1` to `3`, the glued word is a Hamiltonian path of `Q_{25n+12}` from `1` to `3` (hence a Hamiltonian cycle), and `25n + 12` is odd. Iterating from the base chain of length `35` (`1,8,28,21,4,32,17,19,6,30,34,15,10,26,23,13,12,24,25,11,5,20,29,35,14,2,7,18,31,33,16,9,27,22,3`) gives Hamiltonian cycles of `Q_N` for `N = (71*25^m - 1)/2`, explicitly verified for `m = 1, 2, 3`: `N = 887, 22187, 554687`. *Proof.* The block sums are `25` times a square; the junction sums are the `36` constants above; the multiset of values is `{1..12}` union the `25` residue-shifted copies of `25 a_k`, i.e. `[1, 25n+12]`; the length is `12 + 25n`, odd when `n` is; the first and last letters are `1` and `3`. QED

## 7. Literature status (S7, S7b; CITED)

OEIS A090461 (T. D. Noe, Dec 01 2003) carries, verbatim in the `.out`, the comment chain: the original conjecture "sequence includes all integers k > 24"; Jud McCranie (Jan 11 2018) "25..299"; Robert Gerbicz (dated "Jan 17 2017" in the entry, evidently 2018 from the forum) the `2^20` verification plus the `(71*25^m-1)/2` family; and Robert Gerbicz (Jan 21 2018): the conjecture is proved, every `k >= 25` has a chain and every `k >= 32` a Hamiltonian cycle, with code and a deterministic algorithm on the Mersenneforum topic. A078107 (R. K. Guy) records via Paolo Xausa (May 29 2024) that it has no further terms. A071984's comment (Bert Dobbelaere, Dec 28 2018) gives the degree argument for no cycle at `n <= 30` and states `n = 31` "can be shown by hand". The live forum URL returned `404` on 2026-09-22; the thread is readable on web.archive.org, and its method (post #22, 2018-01-21) is, in this lane's words: blow a "nice pair" (chains of lengths `n` and `n+1` with every value in same-parity positions) up by `49` with shifts `c = -24..24` and glue with `1..24`, obtaining nice pairs for `49n + res`, `res = 24..72`, a complete residue system mod `49`; base table `n = 41..2032`; smaller `n` by lookup; `O(n log n)` time (`n = 10000000` in about `5` seconds). Post #9 (2018-01-12) describes the `2^20 = 1048576` verification by the flip-and-insert extension `F(i, j)` (fresh search last needed at `n = 6109`). Post #21 (2018-01-17): the square must be odd because for an even square `c = k^2/2` has the same residue as `-c`; that `9` is too small to glue is henryzz's remark in the same exchange, and Gerbicz's post #22 only says that for `9n + res` "there is nothing on our plate". The audit read posts #9, #19, #21 and #22 directly in the archived page 1 of the thread (posts #1-#25, all dated 2018-01-11 to 2018-01-23) and confirms that the glue word `GLUE` and the base chain `V35` of the script are verbatim transcriptions of post #19 (Gerbicz's `fun_odd`), re-implemented from his PARI pseudo-code in the audit script (A6). The lane verified the post #19 step (S6) and nothing of post #22 beyond reading it: the glue tables live in the linked `squares.c`, which was not fetched. The theorem is therefore CITED, not verified.

## 8. Reproduction block

```text
cd <worktree>
python3 04-computation/experiments/collatz_mod6_20260922_w6_square_sum_hamiltonicity.py \
    > 05-knowledge/results/collatz_mod6_20260922_w6_square_sum_hamiltonicity.out
python3 -O 04-computation/experiments/collatz_mod6_20260922_w6_square_sum_hamiltonicity.py 2>/dev/null | diff - \
    05-knowledge/results/collatz_mod6_20260922_w6_square_sum_hamiltonicity.out   # identical
SQSUM_LIVE=1 python3 ...   # optional: re-fetch A090461 and compare its comments with the embedded copy
python3 04-computation/experiments/collatz_mod6_20260922_w6_square_sum_hamiltonicity_audit.py \
    > 05-knowledge/results/collatz_mod6_20260922_w6_square_sum_hamiltonicity_audit.out   # independent audit
```

Pure Python, no dependencies, well under a minute and a few MB for each script (the numbers of the audit paragraphs, S1a and section 10, are in the audit `.out`).

## 9. Microcosm to macrocosm: what the `n >= 25` theorem is, and the typed non-map to Collatz

**What it is.** Gerbicz's theorem has the shape *finite base* (a verified table of chains, in the final version `n = 41..2032`, with the small cases `25..40` by lookup) *plus an inductive extension* (an explicit embedding `Q_n x {shifts} -> Q_{49n+res}` that multiplies every value by an odd square, uses `25(a + b)` or `49(a + b)` being a square iff `a + b` is, and carries an invariant, the nice-pair parity, that the induction needs). The "fractal repetition in the macrocosm" the user senses is real and completely explicit: the `(71*25^m-1)/2` cycles of S6 are literally `25` scaled copies of the base cycle stitched with the twelve smallest integers, and the general theorem is the same picture with `49` copies and a two-chain invariant.

**Typed analogy (no map found).** Source: Hamiltonian chains in `Q_n`, theorem "finite base + `49`-fold inductive extension". Target: Collatz convergence on the plus sheet. Map: the argument shape "finite set of residue cylinders (base) + inductive extension to all `n`". Preserved: only the shape. Lost: the embedding. The square-sum extension is a constructive self-similarity of the *solution* (scaling a chain by an odd square keeps it a chain); Collatz orbits have no self-similar embedding of orbits into orbits, and the cylinder statements that do extend (Terras words, prefix-descent densities, the gate `q = 1`, the product identity) are proved for both signs, so [counterexample_portrait](collatz_mod6_20260922_counterexample_portrait.md) shows any such "base + cylinder induction" that never uses `n > 0` at `b = +1` is satisfied by the `3n-1` sheet with its three cycles. Sidecar: the nice-pair invariant, the one sign-specific-looking ingredient of the source (it pins positions mod `2`); its Collatz analogue would be an invariant of residue cylinders that is *not* preserved by `n -> -n`, and the repository has no candidate (the portrait's sign-specific objects are the side of `log_2 3` in the sign law, the plus-sheet floor `N >= 2^68`, and the budget-exhaustion discriminator `sum q_L = 3 n_0` verified on the three minus cycles; none of them is an embedding of orbits into orbits). Test: none available, because the map stops at the shape. The lane therefore records SCOPE: the square-sum theorem is a model of how a "residue base + induction" proof looks when the inductive step is an explicit embedding, and it is evidence *against* expecting such a proof for Collatz without a sign-specific embedding.

## 10. Audit record (2026-09-22)

| Claim | Verdict | Note |
|---|---|---|
| S1 components `3 / 2 / 1` at `4..12 / 13 / 14..40`, leaf lists, self-loop values `2, 8, 18, 32` | CONFIRMED | union-find, A1 |
| S2 path set `{1} u A090461`, cycle set `32..40`, session lead paths | CONFIRMED | weak-prune DFS, A2 |
| S3 counts `A071983` (`15..36`), `A071984` (`32..40`), `A090460 = A071983 - (n-1) A071984` | CONFIRMED to `33` (paths) and `36` (cycles) independently; `34..36` paths and `37..40` cycles rest on the lane's search plus the live OEIS match (A2) | |
| S4 certificates `18..22, 24`, cut sets, `n = 23` endpoint counts `{2:1, 9:1, 18:3, 22:1}` | Original confirmation SUPERSEDED on 2026-09-25: cuts/counts survive, unguarded degree-2 forcing is false; see S4 | A3, corrected A4, A5 |
| S5 `n = 15` census, squares `9, 16, 25` | CONFIRMED | A5 |
| S6 25-fold blow-up, `887, 22187, 554687`, `36` junctions | CONFIRMED (independent re-implementation, partition of `[13, 25n+12]` checked) | A6 |
| S7 OEIS comments, dates, authors, `A398909` (Aug 14 2026), the `Jan 17 2017` typo | CONFIRMED against the live entries fetched at the audit | |
| S7b Mersenneforum posts #9, #19, #21, #22 | CONFIRMED against the archived thread; `9 fails` attribution refined | |
| S2a, S2b labelled as new | WEAKENED to classical | textbook facts |
| pasted `N <= 13 disconnected` labelled CORRECT | WEAKENED: fails at `N = 1` | A1 |
| pasted `N >= 25 connected` | strengthened to PROVED for all `N >= 14` | A7 |
| S9 typed non-map | CONFIRMED; portrait wording tightened | |

**2026-09-25 correction.** S4, the status/inheritance/table, the main reducer, and audit A4 now require certified nonendpoints. Fresh script outputs preserve the independent census and cuts. The old forcing certificate is retracted, with the triangle and the explicit `Q_23` endpoint as hostile controls. See [MISTAKES](../../01-canon/MISTAKES.md), entry "Square-sum degree-two forcing forgot the endpoint coordinate".

Edits made at the 2026-09-22 audit (historical): script `forced_reduction` now orders fragments as paths (guarded by forced degree `<= 2`) and prints the forced-edge count `21` at `n = 24` and the junction count `36`; the S8 row for `N <= 13` relabelled; note sections 1 (S1a), 2, 4, 5, 7, 9 and the status line edited as described above.

## Stopping boundary / next question

The finite census remains decided and matched to OEIS; failures at `18..22` have structural certificates, while `24` now has the independently audited [endpoint-aware structural proof](decoder_prime_square_20260925.md), section 4, as well as two independent exhaustive certificates; the `25`-fold step is proved and the `49`-fold theorem is cited. The lane stops at the citation boundary: the next cheap, honest step is to fetch `squares.c` (or reconstruct the `49`-fold glue by the same finite search the script uses for `25`) and turn the CITED post #22 induction into a PROVED-by-finite-check statement in the repository, exactly as S6 did for post #19; that would also decide whether the nice-pair invariant can be dropped by using two odd squares (`25` and `49`) instead of one square with two chains. On the Collatz side the lane has no next question: the only thing the square-sum theorem contributes is the shape of an argument the portrait already excludes sheet-blindly.
