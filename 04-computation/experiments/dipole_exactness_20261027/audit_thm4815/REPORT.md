# Independent audit of THM-4815 (the dipole principle, cut or colour, max-to-min paths, P_7 parity)

**Verdict: PASS WITH CORRECTIONS** (2026-10-27). Every proof (A1, A2, B1-B8 and the parity law of C) is correct, and
every census number in D reproduces exactly with the auditor's own code. The corrections are a missing literature
credit for B7, a false clause in the title, and several factual slips in D/E, the note and the reflection. The auditor
wrote 35 exact edits (`required_edits.json`, built by `make_required_edits.py`; each old string occurs once, and all 35
applied cleanly to in-memory copies); the producer applied all 35. The harness blocked the auditor's own report file, so
this is its report as delivered, saved by the producer. No repository file outside this directory was edited and nothing
was committed by the auditor; all audit code is its own (producer code was read only).

## Major findings

**M1. B7 is classical; Lu was misdescribed; Thomassen was not cited.** The mathematics of B7 is right, the credit is
not.

- Thomassen (Proc. London Math. Soc. (3) 45 (1982) 151-168; Theorem 5.5 of Bang-Jensen and Gutin's survey "Paths, trees
  and cycles in tournaments", whose text the auditor read) characterises exactly when a given maximum-score vertex `x`
  and minimum-score vertex `y` are joined by a Hamiltonian path from `x` to `y`. The only exceptions are:
  - (a) `T - x` ends with the component `{y}` followed by a single vertex `z`, with `z -> x` and `x -> y`;
  - (b) the mirror image of (a) in `T - y`;
  - (c) `R(1,2)` and two small tournaments.

  In (a), `z` also has the minimum score and beats `x`, so the pair `(x, z)` works; (b) is symmetric. B7 therefore
  follows at once ((c) by direct check). The characterisation was cross-checked against an own endpoint computation for
  every (max, min) pair at orders 3-9: the only pairs not explained by (a) or (b) are in `R(1,2)` and in two 6-vertex
  classes.
- Lu (Graphs Combin. 34(4) (2018) 613-618, doi 10.1007/s00373-018-1899-3; abstract only): if a minimum-score vertex beats
  a maximum-score one, a path from the maximum to the minimum always exists, which is exactly B7 case (i); if the maximum
  beats the minimum, such a path exists except in a few cases. THM-4815 said Lu "proves them with a few exceptions"
  (wrong for case (i)) and claimed case (i) and the existential statement as "proved here in full".
- The status block cited a "reversal lemma … (given Lu 2018)" that does not exist in the body, omitted B7 and B8, and
  contradicted "proved here in full".
- B8 (the distance to the regular class is `R` at odd orders) remains a genuine corollary for native moves.

**M2. The title's "strong failures at orders <= 9 all have a cut vertex" is false.**

- Strong 2-connected `L`-failures exist: 1, 5 and 143 classes at orders 7, 8 and 9. Example at order 7:
  `110011111011111100110`, scores `2,2,2,3,4,4,4`, `|Aut| = 1`, self-converse, `lambda = kappa = 2`, `L = 3`, distance 4
  (confirmed by a labelled BFS over all `2^21` tournaments, `kappa` by brute force).
- These are failures but not local minima. Only the strong `L`-local minima have cut vertices (0, 1, 3, 11, 49, 242,
  1,302 classes, all `kappa = 1`), which is what the body proves and counts.
- For the even-order almost-regular target, 30 of the 70 local minima at order 8 are 2-connected, so "cut" cannot mean
  "cut vertex" there.

## Minor findings (all corrected by the 35 edits unless noted)

1. Smallest rigid failure (D, E.2, note section 2, reflection section 5): order 5, not 6. `K1 => S4` (`1111100111`) and
   `S4 => K1` (`1101111111`) are rigid, not self-converse, `L = 1`, distance 2. Strong rigid failures begin at order 6
   (10 of the 20 rigid `L`-failures there). All three strong `L`-local minima of order 6 have scores `1,1,2,3,4,4`: the
   named class, its converse `101011111111110`, and the self-converse `111011111111110`.
2. The forced-arc sentence for the named class was false: "every Hamiltonian path through them contains `5 -> 3 -> 4`"
   holds for only 5 of the 15 paths. The conclusion is right: the only path from 3 ends at 5 and the only path from 5
   ends at 2, so no path joins `{3, 5}` to `{0, 1}`.
3. The A3 table credited the regular class to B6; it is B8.
4. "Order 9 in the note": the note had no order-9 class-graph data ("in progress at filing"). The producer's later
   `exact_n9_summary.out` agrees with the auditor's order-9 BFS.
5. "Two at order 6 … for the regular target": order 6 is even, so that target is the almost-regular class.
6. "960,000 tournaments in all": the output file listed 1,120,000 distinct samples (1,140,000 after two further runs
   were added).
7. "Symmetric" was never defined: it means `|Aut| > 1` or self-converse, and the column counts local minima, not failures.
8. The credit to THM-4818 was overstated: THM-4818 (C2) already uses the endpoint form of the dipole law and the fact that
   every path starts in the top component and ends in the bottom one. New in B1-B3 are the criterion of B2 and the exact
   converse `Pi = C_1 x C_k`. THM-4818's "748 of 904" is consistent: 748 non-strong plus 42 strong rigid local minima make
   790 = 904 - 114.
9. The alternating half of the parity law follows in one line from THM-4810 Theorem P's odd toggle count and is already
   used in THM-4732 Lemma 1; it should be credited, not presented as new.
10. "Session74's source/sink trap" is its "safe pruning lemma"; Session74's "trap" is a different object, a strong
    `L = 1` state at order 11. The auditor checked that state: it is strong, has cut vertices `{2, 3, 7, 9, 10}`, and no
    Hamiltonian path joins a descending pair.
11. Section C was stale: the reflection row and the HYP-9420 (4) body had since been corrected, HYP-9420's title still
    said `d(P_7) = 7`, and MISTAKE-682's main entry said the reflection was "left as written" while its addendum said it
    was corrected (the producer added a note to the addendum).
12. Cosmetic credits: the `c3` law originates in black_descent_20261020c; the note's credit for the identity omitted the
    dyadic-module-spectrum session, which THM-4811 credits; THM-4814 (C) has no per-move tariff, and THM-4814 (D) covers
    regular sources only (the producer amended the A3 rows).
13. The THM-4730 paraphrase should read "palindromic `{K1, C3}` halves with centre in `{none, K1, C3, T5}`"; the two
    readings agree at orders <= 9.
14. B2 should say that the potential "vanishes exactly on the target" (the producer amended B2).
15. Provenance: the stored producer outputs came from older code versions; the "526 arcs" and "67,962 arcs" figures had
    no stored output (both reproduced here); the producer's check C7 tests B7 on strong classes only.

## Proof check

- A1 and A2 are correct: all identities hold on every labelled move at `n <= 7` (165,150,720 moves at `n = 7`). Colour
  is preserved by every move, so the black class graph is undirected.
- The A3 tariffs were re-derived: the `c3` change is exactly `s(a) - s(b) - 1` (all labelled moves at `n = 4-6`), `fas`
  changes by at most `n - 1`, the path norms are as stated (`n = 3-12`), and the unit tariffs hold on every class edge
  at `n <= 8`.
- B1-B6 are correct. B5's claim that optimal matchings are blockwise holds because crossing between separated score
  blocks strictly increases the cost.
- B7 is complete and correct. Cases (i) and (ii) are exhaustive (either some min vertex beats some max vertex, or
  `X => Y`); the path is assembled correctly; `n = 3` (`T'` a single vertex) and `n = 4` (`T'` two vertices, not strong)
  are handled; the counting, the extremal arguments and the `n = 2` contradiction all check. The proof gives more than
  stated: every `y -> x` pair, and every pair when `X => Y`, is joined (exhaustively verified at `n <= 10`).
- B8 is correct: the descending pairs are exactly (start above the mean, end below it), checked on every endpoint pair at
  odd `n <= 9`; an exact BFS gives distance `R` at `n = 5, 7, 9`. With black moves only, the single failure is the
  order-7 class with scores `1,3^5,5` and `|Aut| = 5` (black distance 4 against `R = 2`).
- C is correct: `F' = C(n,2) + n - 1 - F - 2 back(Q)` holds on all labelled moves at `n = 3-7`; the move graph is
  bipartite at `n = 3, 4, 7` and has two parity components at `n = 5, 6`. `P_7` is 7 moves from `TT_id` and 6 from the
  transitive class, the 5,040 transitive labellings splitting 2,520 at distance 6 and 2,520 at distance 7; the six-move
  word is valid and ends at a transitive tournament with `F = 12`. Session74's 15-move `P_11` certificate replays
  (`L` runs 15 -> 0); `F(P_11) = 33` and `F(TT_id) = 55`, so the labelled distance to `TT_id` is even and at least 16.

## Reproduction

| claim | result |
|---|---|
| B3, B4, B5, B7 at `n <= 10` (all 9,733,056 classes at `n = 10`) | 0 failures |
| `delta0` local minima 0, 0, 1, 3, 16, 99, 915 | identical |
| non-strong `L`-local minima 2, 5, 19, 99, 855, 13,303, 376,191 | identical |
| strong `L`-local minima 0, 1, 3, 11, 49, 242, 1,302 | identical, all `kappa = 1` |
| Lu-type failing arcs, `n = 4..10` | 2, 2, 6, 14, 68, 526, 8,212 arcs in 1, 1, 5, 12, 64, 514, 8,156 classes; never untied; `y -> x` arcs never fail (67,962 at `n = 9`) |
| strong 2-connected classes with no min-to-max path (`n = 7, 8, 9`) | 1, 26, 590 (and 19,128 at `n = 10`), all still descend |
| D failures (`L` / `delta0` / `lambda`), `n = 4-8` | 2, 6, 33, 182, 1,431 (symmetric 2, 4, 13, 52, 173) / 0, 0, 1, 3, 16 / 0 |
| D local minima (`L`), `n = 4-8` | 2, 6, 22, 110, 904 (symmetric 2, 4, 10, 32, 114) |
| almost-regular target, `n = 8` | 306 failures, 70 local minima |
| black class graph, `n = 4-8` | one component plus 1, 3, 2, 4, 3 pure-blue classes; colour obstructions 2 (`n = 6`) and 1 (`n = 7`), all self-converse |
| order-9 BFS | `L` failures 17,936 (symmetric 1,062), local minima 13,545; `delta0` failures 99; `lambda` failures 0; distance to the regular class `= R`, any colour and black only (6 pure-blue classes unreachable) |
| random tests of B7 | 0 failures in 1,503,000 tournaments at `n = 11-20`, five samplers |

## Files in this directory

- `required_edits.json` (35 edits) and `make_required_edits.py`;
- the main checkers `audit_endpoints.c`, `audit_classgraph.c`, `audit_bfs_onfly.c`, `audit_labelled.c`,
  `audit_thomassen.c`, `audit_b7_random.c`, each with its `.out` files;
- the brute-force checks `audit_named_classes.py`, `audit_n6_strong_lmin.py`, `audit_session74_trap.py`,
  `audit_a3_tariffs.py`, `audit_p11_replay.py`, with their `.out` files.

Sources consulted by the auditor: Lu 2018 (Springer, abstract), and Bang-Jensen and Gutin, "Paths, trees and cycles in
tournaments" (survey PDF, https://www.cs.rhul.ac.uk/home/gutin/paperstsp/t1.pdf), Theorem 5.5 (Thomassen 1982).
