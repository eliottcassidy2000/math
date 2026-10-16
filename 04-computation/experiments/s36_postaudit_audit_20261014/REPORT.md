<!-- Second independent adversarial audit of session opus-2026-10-14-S36 (the post-audit continuation: HYP-9280 "Post-audit continuation" and section 11 of strong_balance_pisano_golden_eleven_20261014.md). Run by separate agents with their own code in this folder (a first run stalled after items A, B; a second run completed C-H). Saved verbatim from the auditor's final message by the session (subagents cannot write report files). Findings E1-E3 and I1-I11 were applied; see MISTAKE-617. -->

# Referee report: S36 post-audit continuation (HYP-9280 "Post-audit continuation"; results file section 11)

## Scope and methods

**What was checked.**
- The HYP-9280 section "Post-audit continuation (same session)" and its status block.
- Section 11 of `05-knowledge/results/strong_balance_pisano_golden_eleven_20261014.md`.
- The four new session scripts and their outputs.

**Methods.** All counts were recomputed without the session's subset DP:
- explicit DFS over Hamiltonian paths (`pa_core.enumerate_paths_dfs`);
- brute force over all n! permutations (`pa_core.brute_force_perms`, or `itertools` in (G));
- the Landau score test plus the Warshall closure for strongness;
- exact integers and Fractions throughout.

Items (A) and (B) and the core self-test come from the earlier, stalled run. The unfinished (C) script was restructured into short modes.

## ERRORS

**E1. The `0.047` lower bound is false for d = 1, 2.**
- *Where.* HYP-9280, "Where it fails": "`min(E[N] − d/2) >= 0.047` for every `d <= 6` in every sweep". Results file 11.1 repeats it.
- *The session's own output* (`tournament_near_transitive_sweep_20261014.out` prints exactly `E[N] − d/2`) gives smaller values:
  - d = 1: `0.0001` at n = 15 and 16 (feedback 2);
  - d = 2: `0.0057` at n = 16; also `0.0084` (n = 15) and `0.0115` (n = 14, feedback 3).
- *What does hold.* The bound is true only for `3 <= d <= 6`, where the minimum is `0.0471` (d = 6, runs (9,2) and (10,4)).
- *Not uniform either.* The d = 3 minimum falls with n in the feedback-3 runs: 0.1796 at n = 8 down to 0.0731 at n = 14.
- *Fix.* "`E[N] − d/2 > 0` for every `d <= 6`. The minimum over `3 <= d <= 6` is `0.0471`. For d = 1, 2 it tends to 0, but those degrees are proved."
- The qualitative conclusion stands. In fact the least in-degree at a `D_in <= 0` vertex is `>= 9` in every run, not just `>= 7`.

**E2. The status block claims "the sharp T_n minima" as FINITE-EXACT over the feedback-arc families (`<= 2` for n `<= 16`, …, `<= 5` for n = 10). That is false in part and unsupported in part.**
- *Which minima.* The HYP's own list of T_n minima includes "min over vertices of `D_in`, and of `D_out`: `1/H(T_n)`".
- *False part.* In exactly these families `D_in` and `D_out` are negative from n = 11 on: `−47/297` at n = 11 (feedback 2), down to `−6386/8385` at n = 16.
- *Unsupported part.* The session's sweep script never computes the 3-cycle certificate, `Pr[w first]` or `D_out`.
- *What the referee runs supply.* audit_b plus audit_b2 cover all 22 runs:
  - no tournament has a certificate at or below T_n's value;
  - none has `Pr[w first]` at or above T_n's value;
  - the minimum vertex-lemma value in every run with feedback `k >= 2` exceeds T_n's (e.g. n = 16: `3/24579 > 1/16385`).
- *Fix.* Restrict that clause to: F strong, vertex lemma, no majority source, and the T_n minima of the vertex-lemma value, 3-cycle certificate and `Pr[w first]`. State that the `D_in`/`D_out` minima fail there.

**E3. "Not the minimum" is stated for n = 8..12, but it is wrong at n = 8, 9 and unsupported at n = 10.**
- *Where.* HYP, "Conjectured sharp forms": the D_in annealing value `1/(2^(n−2)+1)` for n = 8..12 "was … not the minimum". Results 11.1: the reports at "n = 10, 11, 12" were "a search artefact".
- *n = 8, 9.* It is the minimum. The HYP itself records this as FINITE-EXACT, so this is an internal contradiction.
- *n = 10.* Nothing refutes it. The minimum `D_in` in the feedback `<= 5` runs is `3/387`, `5/517`, `5/517`, `9/681`, all above `1/257 = D_in(T_10)`.
- *Fix.* Only the n = 11, 12 reports were artefacts.

## IMPRECISIONS

**I1. Wording of the PROVED hypothesis.** The text says "o beats every *other* in-neighbour of w", but o is not an in-neighbour of w. It should read "o beats every vertex of `U = V − {w, o}`", equivalently `indeg(o) = 1`. This appears in the status block and in 11.3.

**I2. Scope of the out-degree-1 FINITE-EXACT claims is overstated.**
- *The text.* "Labelled n <= 7 and near-transitive families to n = 16".
- *What the script actually ran:*
  - feedback 2 at even n = 10..16 only;
  - feedback 3 at n = 10..12;
  - feedback 4 at n = 10;
  - n = 8 and 9 were not covered at all.
- *My extension (no failure anywhere):*
  - exhaustive n = 8 (all 1,866,256 strong configurations);
  - odd-n near-transitive families;
  - random n = 6..12.

**I3. "Up to 395 per near-transitive sweep" holds only for the session's own sweeps.**
- The U+ average is already negative at n = 8, in 2,160 of 1,866,256 configurations.
- Other sweeps give larger counts: 644 at (13,3), and 1,125 at (11,4) plus one case where the U+ average is exactly 0.
- The session counted "<= 1/2", but all of its cases are strictly negative (re-checked).

**I4. The G rule ("beats the fewest out-neighbours of w") already fails at n = 8, not first at n = 9.**
- *The counterexample.* Reversed arcs (0,5), (1,6), (5,7), which has feedback arc number exactly 3, with w = 1. The unique G-choice u = 5 has `p = 88/179 < 1/2` (brute force over 8!).
- *Counts.* 9 (T, w) pairs at (8,3) and 16 at (8,4). The session never ran n = 8 with feedback `>= 3`.
- *Tie-breaking.* "Survives" holds only under the most favourable tie-break. A tied G-choice is a non-witness in 1,440 labelled (T, w) at n = 6 (e.g. code 18, w = 1) and in 28,560 at n = 7.

**I5. The W rule is correct but should state the vertex and value.**
- *Vertex and value.* w = 1, value `−656/409² = −656/167281 ≈ −0.003922`, strictly negative. The vertex lemma still holds there (u = 4, `p = 213/409`).
- *Float in the session script.* It computes W in floating point. This is exact at these sizes, but it should be integer.
- *Untested families.*
  - n = 8 was never tested by the session. No failure at n = 8 (feedback `<= 4`) or at n = 9 (feedback 4).
  - There are 4 failures at (10,4). The first is a relabelled copy of the same profile.

**I6. "(SS) at k = 1 cannot be realised by local moves" overgeneralises a seeded sample.**
- *Only two move sets were tested:* single-vertex relocations and 2-block relocations.
- *Not all pairs fail.* Matchings do exist for 23 of 740 pairs (single moves) and 236 of 740 (2-blocks) at n = 8.
- *The obstruction is mostly trivial.* In the n = 6 recount, 447 of the 654 deficient cases contain a path with N <= 1 that has no relocation into `{N >= d−1}` at all.
- *Fix.* Type it NUMERICAL and narrow the wording.

**I7. "A proof … must be global" is a non sequitur.** It appears in the HYP and in 11.4. The failure of five witness rules and one refinement does not exclude every local argument.

**I8. "Fails in 41% of cases at n = 7" mixes up the denominator.** The 41.30% is over triples (T, w, u) with d >= 3 (11,734,800 of 28,413,840), not over (T, w) cases.

**I9. Section 11 refers to a "section 12" that does not exist.** The file ends at 11.4.

**I10. "Dually out-degree <= 4" needs a clarification.**
- The dual statement gives an out-neighbour that follows w with p > 1/2.
- That is not the vertex lemma at vertices of small out-degree, which is exactly where every known (SS) and D_in failure sits (out-degree 1–3).

**I11. Coverage of "feedback arc number <= k" also needs the T_n layer.** T_n is the only strong feedback-1 tournament. It lies outside the 22 runs and is covered by separate checks, which the ranges paragraph should say.

**Remark (no error).** The remark that the swap involution "reduces it to a comparison" is accurate. The exact identity is

    Σ_{u∈U−}(C(u,w) − C(w,u)) = |U−|·N1 + S_N,

- N1 is the number of paths ending at w;
- S_N is the signed U− count over paths that start in U−;
- the swappable class contributes exactly 0.

Verified per path on 1,833,488 exhaustive n = 8 cases and on random and near-transitive families. S_N alone can be negative (4 cases at (12,2), 8 at (14,2), 9 at (12,3)), so the comparison is genuine. The minimum U− average also tends to 0 (0.0117 at n = 16, feedback 2), so there is no uniform slack.

## Findings by item

**(A) Double chains: CONFIRMED** (earlier run, brute force over 11! and 12!, DFS for n <= 18).
- D(11,5): H = 289, w = 0 with d = 9; `p(u before 0)` is 9/17 for u = 1..4, 161/289 for u = 6..9, and 33/289 for u = 10; `E[N] − d/2 = −23/578`; best in-neighbour p = 161/289.
- D(12,6): H = 561, `E[N] − d/2 = −91/561`.
- D(11,6): `−47/594`, the n = 11 minimum (`D_in = −47/297`).
- All D(n, j) are strong. The minimum vertex-lemma value is T_n's, attained at j = 1 and j = n−2 (both copies of T_n). `D_in <= 0` first appears at n = 11 and falls monotonically to −9552/11051 = −0.8644 at n = 18. M is never transitive and F is always strong.

**(B) Near-transitive sweep: CONFIRMED.**
- All 22 session runs reproduced exactly: 12 in audit_b plus the T_n layer for n = 4..20, and 10 in audit_b2. Examples: (14,3): 8,937 strong, 337 with `D_in <= 0`, min `D_in = −1500/3315`, (SS) failures 89/404/587; (12,4): 93,779 strong, 140, `−166/585`; (10,5): 197,076 strong, min vertex-lemma margin `9/681`.
- The total of 401,442 strong arc sets is correct, and the coverage argument is sound.
- New in the referee runs: (SS) k = 2 failures all sit at out-degree-1 vertices with d = n−2 >= 8 (e.g. D(10,5), w = 0, slack −1/153); the certificate is always above T_n's value; `Pr[w first]` is always below T_n's value.

**(C) Out-degree-1 vertices: the PROVED identity is CONFIRMED with no gap.**
- The proof: the successor of w is necessarily o; paths ending at w are `o Q' w`, because o has in-degree 0 in `T[U ∪ {o}]`, giving `H(T[U])` paths with every u before w; all other paths are `A w o B'` with `(A, B')` free, since `U ⇒ w` and `o ⇒ U`; the swap `(A, B') ↦ (B', A)` exchanges "u before w" with "w before u".
- Numerical checks: 403 random shapes, n = 4..12: 0 failures; the decomposition `H = H(T[U]) + Σ H(T[X])H(T[U−X])`: 385 checked, 0 failures; the T_n value holds for n = 4..13.
- Empirical claims (a witness lies in U+; if U− is nonempty its average 2p − 1 is positive; the U+ average can be negative): no failure. Exhaustive n = 8: 2^21 configurations, 1,866,256 strong (equals the condensation count), U− empty in 32,768 and nonempty in 1,833,488, U+ average negative in 2,160. Random samples: 14,923 planted and 10,961 unplanted (T, w) pairs. Near-transitive: 16 sweeps, reproducing every session count, e.g. (12,3): 6,664 / 630 / 6,034 / 395.

**(D) Witness rules: CONFIRMED, counts exact.**
- The W rule at n = 10: H = 409, DFS and brute force over 10! agree. Only w = 1 has `W <= 0`, value −656/167281.
- Labelled sweeps n <= 7: every count matches (n = 6: E 1,440 / 1,440; G all-fail 1,440; n = 7: A 5,040 / 20,160; D 20,160 / 20,160; E 23,520 / 23,520; G 0 / 28,560; P 5,040 / 5,040; W 0). All 10 near-transitive session runs match.
- Explicit n = 7 examples (brute force over 7!; out-sets "v: [x with v → x]"):
  - (i) minimum out-degree, no witness: labelled code 551, out-sets `[[1,2,3,6],[5],[1],[1,2],[0,1,2,3],[0,2,3,4],[1,2,3,4,5]]`, w = 4; in-neighbours 5 (out-degree 4) and 6 (out-degree 5); the unique minimiser u = 5 has p = 20/41;
  - (ii) most frequent predecessor, no witness: labelled code 161, out-sets `[[1,6],[3],[0,1],[0,2],[0,1,2,3],[0,1,2,3,4],[1,2,3,4,5]]`, w = 0; predecessor probabilities 15/41, 13/41, 7/41, 3/41; the unique top predecessor u = 2 has p = 20/41.

**(E) (SS) at k = 1: CONFIRMED, no failure.**
- Algebra (exact LP over all distributions): under k = 0, 1 the minimum of E[N] is exactly 3/2 for d = 3 and 2 for d = 4; parity then gives p > 1/2. For d = 5, k <= 1 only gives 2, so k = 2 is needed.
- Tests: random strong tournaments n = 8..13: 10,175 pairs, minimum relative slack 0.0395; planted out-degree-1 and out-degree-2 vertices n = 8..11: 21,989 pairs, minimum slack 0.0303; a greedy arc-flip search at n = 8..10 converges to T_n with minimum slack exactly `1/H(T_n)`; all 22 near-transitive runs.

**(F) Refinement counts: CONFIRMED exactly, including n = 7.**

| n | `a_u > b_u` | `a_u = b_u` | `a_u < b_u` |
|---|---|---|---|
| 5 | 720 of 2,280 | 1,080 | 480 |
| 6 | 82,320 of 222,720 | 72,480 | 67,920 |
| 7 | 11,734,800 of 28,413,840 | 4,719,120 | 11,959,920 |

**(G) Local moves: CONFIRMED.** Five explicit cases at n = 6 and 7, each with a Hall violator re-checked from the edge lists (e.g. n = 6, out-sets `[[1,3,4,5],[3],[0,1],[2,5],[1,2,3],[1,2,4]]`, w = 1, d = 4: `|X| = 7`, maximum matching 6, the path (5,1,3,2,0,4) has no move into Y; the same T with w = 2: a genuine Hall obstruction with no isolated path, |S| = 5, |N(S)| = 3). The n = 6 part of the session sample (`random.seed(1)`), recounted with independent Hopcroft–Karp code, gives exactly 1,203 / 654 / 99.

**(H) "Nearly transitive" remark: CONFIRMED.** Exact feedback arc numbers: the two n = 9 (SS) counterexamples 2 and 2; D(11,5) 2; the W example 2; the n = 8 G example 3; the n = 7 A and P examples 2 and 2.

## VERDICTS

| Claim | Verdict |
|---|---|
| (A) Double chains, −23/578, −91/561, −0.864 at n = 18 | CONFIRMED |
| (B) 22 runs, 401,442 strong arc sets, all printed numbers | CONFIRMED |
| (B′) "`min(E[N] − d/2) >= 0.047` for d <= 6" | REFUTED as stated (E1) |
| (B″) status line "sharp T_n minima" over the feedback-arc families | CONFIRMED WITH FIXES (E2) |
| Annealing "not the minimum / artefact" at n = 8..10 | REFUTED for n = 8, 9; unsupported for n = 10 (E3) |
| (C) PROVED identity | CONFIRMED (wording I1) |
| (C) empirical out-degree-1 claims | CONFIRMED WITH FIXES (I2, I3); extended to exhaustive n = 8 |
| (C) reduction remark | CONFIRMED |
| (D) W at n = 10 | CONFIRMED: w = 1, −656/167281 (I5) |
| (D) A, D, E, P rules fail by n <= 7; examples | CONFIRMED |
| (D) G rule "n = 9" | CONFIRMED WITH FIXES (I4): fails at n = 8 |
| (E) (SS) k = 1 no failure; algebra | CONFIRMED |
| (F) refinement counts | CONFIRMED |
| (G) matching deficiencies | CONFIRMED; wording CONFIRMED WITH FIXES (I6) |
| "must be global", "41% of cases", section 12 reference | IMPRECISE (I7–I9) |

## Scripts

All in this folder, each with a matching `.out`. Earlier run: `pa_core.py`, `audit_0_core_selftest.py`, `audit_a_double_chain.py`, `audit_b_near_transitive_dfs.py`. This run: `audit_b2_remaining_runs.py`, `audit_c_outdegree_one.py` (modes identity, reduction, random, near, exh8, typeN; the `.out` is their concatenation), `audit_d_witness_rules.py`, `audit_e_ss_k1.py`, `audit_f_refinement.py`, `audit_g_local_moves.py`, `audit_h_fas_examples.py`.
