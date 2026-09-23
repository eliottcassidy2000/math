# Althöfer's 3n±1 game (Conway's Beans-Don't-Talk): no draws below 2^32, exact heights, ascending rays and a log2 3 phase law

**Status.** OPEN: the prize question (no odd start is a draw) is not settled here or, as far as can be checked, anywhere in public. CITED (links in §1): the rules and prize (Althöfer's page, updated 2026-09-21); Hartisch's finite check below 10^6; Pochuev's manuscript "under review" (page, 2026-09-07); a public repository under Pochuev's name that labels its own global proof OPEN; the identity with Conway's Beans-Don't-Talk and Guy's 1996 Problem 42. FINITE-EXACT (bounds stated, §2):
* every odd start n < 2^32 has finite remoteness (no draw): N = 1,122,146,470, P = 1,025,337,178 (n = 1 counted as N);
* the largest height needed below 2^32 is 13,112,245,406,701, attained at n = 3,238,138,315 (a P-position). The height is the least cap C for which n is decided using positions <= C only;
* exact heights are known for every odd n < 10^9 (maximum 1,499,556,851,731 at n = 703,036,003);
* cross-checked by an independent Python retrograde, by OEIS A005694–A005698, and by proof certificates that an independent checker validates.

PROVED (§3.1, elementary, with proofs): the ascending-ray decomposition, the exit-index parity rule and the height bound it gives, the phase identity, and negation equivariance. EMPIRICAL, on stated finite ranges (§3.2–§4): the phase window, the P-density profile, move statistics and the r <-> −r symmetry.

**Reproduction.** Run from the worktree root. Binaries and large files go to `scratch/procgen_game/`.

```bash
E=$PWD/04-computation/experiments; D=$PWD/scratch/procgen_game; mkdir -p $D; cd $D
for p in dense sparse blocks cert struct remote; do cc -O3 -march=native -o game_$p $E/collatz_procgen_20260922_game_$p.c -lm; done
./game_blocks 8589934593 1 4294967297 11000000 400000        # no draws below 2^32   (263 s, 0.85 GB)
./game_sparse 8589934593 1 1000000001 40000000 t1e9.bin      # exact heights < 10^9  (67 s, 0.9 GB)
./game_sparse 8589934593 3233354757 3250131973               # hardest start below 2^32 (61 s)
./game_dense 2097153 0 2097153 dump.bin && python3 $E/collatz_procgen_20260922_game_check.py bfs 2097153 dump.bin
./game_cert 16777217 c.txt 703036003 3238138315 && python3 $E/collatz_procgen_20260922_game_check.py cert c.txt
./game_struct 2147483649 134217728; ./game_remote 268435457; ./game_remote 536870913 rem_268435457.bin
python3 $E/collatz_procgen_20260922_game_phase.py rem_536870913.bin 4 30
```
Collected output: [`collatz_procgen_20260922_game.out`](collatz_procgen_20260922_game.out).

## What this says

1. The prize problem is Conway's *Beans-Don't-Talk*. Guy asked in 1996 whether it has any O-position, that is, any draw, and the question is still open. The only public claimed proof, a repository in Pochuev's name, describes itself as OPEN.
2. The computation stays cheap if it climbs above a cap only where the proof needs it. This settles every odd start below 2^32, 4,000 times Hartisch's bound, in about 4 minutes. The heights needed have a power-law tail, P(h/n > R) ≈ R^{-2.7}, and reach 4,049·n.
3. P-positions are organised along the **ascending ray**, not the descending chain. x is N exactly when the first position on its ascending ray whose descending child is P has even index. The parity of the descending chain predicts the true value only 37% of the time, so it is anti-correlated with it. The winner needs the ascending move at 4.4% of all odd n (first case n = 15: the only winning move is 15 → 23).
4. The strongest visible structure is archimedean. Every move adds log2 3 to frac(log2 n), up to O(1/n), and optimal games end in the universal endgame 19 → 7 → 5 → 1. So the optimal length T* satisfies frac(log2 n + T*·log2 3) ∈ [−0.0251, +0.0457] for every odd 2^9 < n < 7,301,503. The value is the parity of T*, and frac(log2 n) alone predicts it for 84–100% of starts. This explains four observations:
   * the log-periodic P-density, with 1/16-octave bins ranging from 2% to 93% P;
   * the two children of a position agree in value 91.5% of the time;
   * both moves win at 84% of N-positions;
   * the earlier "~28% P" figure was misleading (see 6).
5. The r <-> −r symmetry of the P-density mod 2^k is negation equivariance in the sheet-sorted form (d, u). The map n -> −n swaps the plus and minus sheets, commutes with d and u, and preserves phases and 2-adic exponents.
6. Correction to the session synthesis: "about 28% P-positions" is an artifact of the cap. It counted unresolved positions, and P needs both children resolved. The true P-density is about 0.48 and oscillates log-periodically. Among odd n < 2^32 it is 0.4775.

## 1. Rules, literature and status (CITED, read 2026-09-22)

* **Althöfer, [Collatz prizes](https://althofer.de/collatz-prizes.html)**. The page was last updated 2026-09-21; the Prize 3 section was updated 2026-08-08.
  * Rules: the state is an odd n. The player to move builds 3n+1 or 3n−1 and halves until the result is odd. Reaching 1 wins; otherwise the other player moves.
  * Prize: 500 EUR "for a proof that all odd starting numbers n will lead to 1, if both players act optimally", deadline 2037-12-31. Because the winner of a non-drawn position forces the end in finitely many moves, this is the same as asking that no position is a draw (infinite remoteness).
  * Hartisch's computer analysis: optimal play ends in 1 for every n < 10^6.
  * Pochuev sent a manuscript with a claimed proof; "Sept 07, 2026: Checks are still underway".
  * The game was created in July 2023.
  * These rules match the task statement exactly.
* **Beans-Don't-Talk.** R. K. Guy, *John Isbell's Game of Beanstalk and John Conway's Game of Beans-Don't-Talk*, Math. Mag. 59 (1986) 259–269, [doi:10.1080/0025570X.1986.11977258](https://doi.org/10.1080/0025570X.1986.11977258). Guy, *Unsolved Problems in Combinatorial Games*, Games of No Chance (MSRI Publ. 29, 1996), [Problem 42](https://library.slmath.org/books/Book29/files/unsolved.pdf), asks for the same moves (n -> (3n±1)/2^*) whether there are any O-positions. OEIS [A005698](https://oeis.org/A005698), A005695, A005696, A005697 and A005694 list the positions of remoteness 2–6.
* **Althöfer–Hartisch–Zipproth**, *Analysis of a Collatz Game and Other Variants of the 3n+1 Problem*, ACG 2023, LNCS 14528 (2024) 123–132, [doi:10.1007/978-3-031-54968-7_11](https://doi.org/10.1007/978-3-031-54968-7_11). Paywalled; not read.
* **Pochuev**, [github.com/Grisha-Pochuev/3n-plus-minus-1-game](https://github.com/Grisha-Pochuev/3n-plus-minus-1-game): created 2026-08-06, last push 2026-09-01, preprint [doi:10.5281/zenodo.21844684](https://doi.org/10.5281/zenodo.21844684).
  * Status: the README says **"OPEN … The unconditional theorem is not claimed."** `docs/global-proof.md` says "OPEN — HOSTILE AUDIT FAILED". The remaining gap is a "high-return token-provenance/lifecycle" obligation, Section 136.
  * What is solid there:
    * an exact binary normal form: with n = 2m+1 the children are F(m) and F(R(m)), where F(m) = ⌈3m/2⌉ and R deletes the maximal alternating binary suffix;
    * a Gray-code transducer;
    * Lean checks of the normal form and of finite outcomes;
    * a proof that finite draw traps are detected by cutoff kernels (the finite-trap/unbounded-ray dichotomy);
    * bounded retrogrades up to 10^7 in conjugated coordinates.
  * It records as an "unverified lead" a resolved prefix through n = 24,561,137. The result here supersedes that lead.
  * It is not established that this repository is the manuscript Althöfer is checking.

## 2. Verification (FINITE-EXACT)

**Definitions.**
* Values are the least fixpoint of two rules: n is N (the mover wins) if some child is P; n is P if both children are N. A child equal to 1 counts as P, because the mover has just won. The start n = 1 is N.
* In the full game the positions these rules never resolve are exactly the draws; in a capped computation "unresolved" only means "not decided within the cap". Every label below is produced by these two rules, so every reported N or P has finite remoteness.
* The **height** h(n) is the least cap C such that the retrograde restricted to positions <= C decides n. Equivalently, it is the minimum over finite proof trees of the largest position in the tree.

**Algorithms.** All three share one idea: odd non-multiples of 3 are inserted in increasing order. Multiples of 3 are never children and are evaluated from their children.
* [`_game_dense.c`](../../04-computation/experiments/collatz_procgen_20260922_game_dense.c) is an *incremental bottleneck retrograde*. It uses 2 bits per position, and the time at which a position resolves is exactly its height. At cap 10^8 it reproduces the earlier sweep program exactly (N 22,210,114, P 14,199,339, unresolved 13,590,547, first unresolved 834,437) in 1 s.
  * First unresolved start B for the dense algorithm alone (B ≈ C^0.76):

    | cap | 2^20 | 2^24 | 2^28 | 2^32 | 2^33 |
    |---|---|---|---|---|---|
    | B | 39,111 | 415,803 | 2,924,165 | 18,752,721 | 45,683,717 |
* [`_game_sparse.c`](../../04-computation/experiments/collatz_procgen_20260922_game_sparse.c) *climbs only where needed*.
  * It continues the same fixpoint above the cap. It inserts only positions that are children of unresolved positions reachable from unresolved targets, in increasing order through a min-heap, stored in a hash map.
  * Its resolution times equal the dense heights on all 4,049 overlap targets below 3·10^6.
* [`_game_blocks.c`](../../04-computation/experiments/collatz_procgen_20260922_game_blocks.c) runs the same process in blocks of starts, with flushes.
  * Its labels are sound, but its heights are only upper bounds: Tmax, the largest position ever inserted, bounds every height.
  * It stops loudly, which never happened, in two cases: if the heap empties with a target unresolved (that would be a certified finite draw trap), or if two flushes pass with no progress.

**Results.**
* **All odd n < 10^9, exact heights.** The run from cap 2^33 had 579,249 phase-2 targets, inserted 9.7M positions and took 67 s. The largest height is 1,499,556,851,731 at n = 703,036,003 (N, h/n = 2133). Three starts, 825,978,239, 929,225,519 and 734,202,879, share the height 940,681,949,237: hard starts share hard sub-proofs.
* **All odd n < 2^32** (blocks, run twice with identical output): 28.2M phase-2 targets, 5.4·10^8 insertions, 263 s, 0.85 GB.
  * Counts: N = 1,122,146,470 and P = 1,025,337,178.
  * Height bound: Tmax = 13,112,245,406,701.
  * An exact single-batch rerun of the block containing it gives h(3,238,138,315) = 13,112,245,406,701 (P, h/n = 4049.3). The maximum height below 2^32 is therefore exactly this value. The upper bound is rigorous; the lower bound relies on the exactness of bottleneck-order insertion, which was checked on the 4,049 overlap targets but not proved.
* **Height tail.** For 10^6 <= n < 3·10^6 the dense exact heights give:

  | R | 2 | 4 | 8 | 16 | 32 | 64 |
  |---|---|---|---|---|---|---|
  | P(h/n > R) | 0.146 | 0.028 | 5.7e-3 | 9.2e-4 | 1.0e-4 | 1.5e-5 |

**Independent checks** (see `.out` §4).
1. An independent Python retrograde ([`_game_check.py`](../../04-computation/experiments/collatz_procgen_20260922_game_check.py)) runs a FIFO over remoteness with child counters, on all odd n including multiples of 3. It shares no code with the C programs. At cap 2^21+1 it matches the C table on all 1,048,576 odd positions, with 0 mismatches and the same first unresolved start, 93,165.
2. Its remoteness layers 2–6 equal A005698, A005695, A005696, A005697 and A005694, including even starts, for every listed term. A005698 is also confirmed complete up to its last term, 1,954,687,338,269, through the Jacobsthal characterization.
3. Proof certificates were extracted by [`_game_cert.c`](../../04-computation/experiments/collatz_procgen_20260922_game_cert.c). Each is a DAG in topological order: an N node names a P child, a P node lists both N children. The Python checker verifies the local rules with its own arithmetic.

   | start(s) | value | nodes | largest node |
   |---|---|---|---|
   | 703,036,003 | N | 1,255 | 1,499,556,851,731 |
   | 3,238,138,315 | P | 1,437 | 13,112,245,406,701 |
   | {7, 13, 15, 23, 27, 817,837, 834,437} | — | 952 | — |
   | 64 pseudo-random starts in [10^9, 2^32) | all 64 agree with the block run | 15,752 | — |

   The largest node of each single-start certificate equals the computed height.
4. Block runs with caps 2^27, 2^28 and 2^33 and different parameters give identical P-counts at every common checkpoint.

## 3. Structure

### 3.1 Elementary facts (PROVED)

Notation: for odd n, let s(n) = +1 if n ≡ 1 (mod 4) and −1 if n ≡ 3 (mod 4). Then 3n + s ≡ 0 (mod 4) and 3n − s ≡ 2 (mod 4). The descending move is d(n) = oddpart(3n + s) <= (3n+1)/4. The ascending move is u(n) = (3n − s)/2 > n for n >= 3.

**L1 (ascending rays).** u is a bijection from the odd positive integers onto the odd positive integers not divisible by 3, and u(1) = 1. The odd integers above 1 therefore split into disjoint, strictly increasing rays R_t = {t, u(t), u²(t), …}, one for each odd multiple t of 3.

*Proof.*
* 2u(n) = 3n − s ≢ 0 (mod 3), so 3 does not divide u(n).
* Injective: from m = u(n) we recover s ≡ m (mod 3) and n = (2m + s)/3.
* Surjective: for odd m with 3 ∤ m, take s ≡ m (mod 3) and n = (2m+s)/3. Then n is odd and 3n − s = 2m ≡ 2 (mod 4), which forces s(n) = s and u(n) = m.
* Rays: u^{-1}(m) < m for m > 1, so iterating u^{-1} ends at a multiple of 3. ∎

In Pochuev's coordinates these rays are the orbits of the "A" move.

**L2 (exit-index rule).** For odd x >= 3, set x_i = u^i(x). Suppose d(x_i) is an N-position for i < k and d(x_k) is P or equal to 1. Then x is N if k is even and P if k is odd. If instead d(x_i) is N for every i >= 0, then every x_i is a draw.

*Proof.*
* x_k is N, because the move to d(x_k) wins.
* For i < k the options of x_i are an N-position and x_{i+1}, so x_i is N if and only if x_{i+1} is P. Walk back down from x_k.
* In the infinite case, suppose some x_i has finite remoteness and take one with the least remoteness. A proof of its value must go through x_{i+1} with smaller remoteness, a contradiction. ∎

In other words, the value is the parity of the first index on the ascending ray at which the descending side child is losing. This is the A-ray/side-branch structure of Pochuev's normal form.

**L3 (height bound).** Under L2, every finite proof tree for x contains x_1, …, x_k. Hence h(x) >= x_k >= 1 + (3/2)^k (x − 1).

*Proof.* At x_i with i < k, the only possible winning move is to x_{i+1}; at a P-node both children must appear in the tree. Also x_{i+1} − 1 >= (3/2)(x_i − 1). ∎

**L4 (phase identity).** Write each move as n_{t+1} = (3n_t + ε_t)/2^{j_t} with ε_t = ±1. For any play n_0 → … → n_T = 1:

log2 n_0 + T·log2 3 − Σ j_t = Σ_t log2( 3n_t / (3n_t + ε_t) ).

*Proof.* Multiply 2^{j_t}·n_{t+1} = 3n_t(1 + ε_t/(3n_t)) over t and take log2. ∎

So each move shifts frac(log2 n) by log2 3, up to a correction of about 0.48/n_t. Both children of n have the same phase up to O(1/n).

**L5 (negation, sheets).**
* d(−n) = −d(n), u(−n) = −u(n), U_±(−n) = −U_∓(n), and |3(−n) − ε| = |3n + ε|.
* (d, u) is the sheet-sorted form of the plus/minus pair: on n ≡ 1 (mod 4), d = U_+ and u = U_−; on n ≡ 3 (mod 4), d = U_− and u = U_+.
* Negation swaps the sheets, maps the residue class r to −r mod 2^k, and preserves all move exponents j_t and all correction terms ε_t/(3n_t). ∎

### 3.2 Statistics on exact values (EMPIRICAL; odd n < 2^27 unless stated; `.out` §5–7)

* **P-density.**
  * By octave it rises to 0.491 (2^11–2^14), then slowly falls: 0.4870 (2^20), 0.4813 (2^26), and 0.4773 averaged over [2^27, 2^32).
  * Within an octave it depends strongly on the mantissa n/2^k. In 1/16-octave bins at 2^16 it ranges from 0.018 to 0.899, and at 2^26 from 0.079 to 0.915. The pattern drifts slowly with k.
  * Consequently the cumulative density up to X oscillates with X: 0.4716 at 5·10^4, 0.4892 at 2^16, 0.4723 at 10^9, 0.4775 at 2^32.
* **Residues.**
  * mod 2^k the dependence is strong: the P-fraction ranges over [0.405, 0.535] for k = 5 and [0.232, 0.726] for k = 12. It is symmetric under r <-> −r: the largest z-score of p(r) − p(−r) over all pairs is at most 1.8 for every k <= 12.
  * mod 3^j there is no dependence: every class mod 81 lies in [0.4817, 0.4827].
  * The P-fraction falls with the exponent j of the descending move: 0.498 at j = 2, 0.457 at j = 4, 0.343 at j = 8, 0.151 at j = 16.
* **Which moves win.**
  * Among N-positions, both moves win in 83.6%, only the descent wins in 8.0%, and only the ascent wins in 8.4%.
  * Forced ascents make up 4.36% of all odd n: 15, 27, 47, 49, 53, 59, 83, 107, 117, 131, … 50.00% of them (to four digits) ascend on the plus sheet.
  * The two children of a position have equal values in 91.5% of cases.
* **Down-only game.** Its value is the parity of L(n), the length of the descending chain to 1. It agrees with the true value for only 37.25% of starts. The first disagreements are 15, 27, 47, 49, 53, 59, 63, 65. The P-fraction is 0.25–0.46 when L is even (down-only predicts P) and 0.49–0.66 when L is odd. The simple rule is refuted; L2 is the correct exact rule.
* **Rays.**
  * L2 holds with 0 violations on x < 32,768.
  * The exit index k is 0 with probability 46.7%, 1 with 47.4%, 2 with 3.6%, 3 with 1.6%, and 4 with 0.6%.
  * Side values alternate along rays: P(s_{i+1} = P | s_i = N) = 0.89 and P(s_{i+1} = P | s_i = P) = 0.03. An independent-geometric model would give P-density 0.35, not the observed 0.49.
  * The hardest starts have small exit index (k = 1 for 3,238,138,315). Their height comes from nested side proofs, not from their own ray.
* **Remoteness (game length).**
  * Values are capped remoteness at cap 2^29. On odd n < 2,924,165 they are identical to the cap-2^28 values for all 974,721 non-multiples of 3.
  * The mean is about 2.06·log2 n: 41.2 on [2^20, 2^21).
  * Records (start: remoteness): 47:23, 659:49, 23297:71, 112413:91, 709063:107, 4094795:130, 7279635:132.
  * In optimal play the winner's fastest win ascends at 9.6% of N-positions (8.1% forced), and the loser's slowest loss ascends at 8.9% of P-positions.
  * Hard lines ascend often: 834,437 has remoteness 98 and 37 ascents in its principal variation.
* **Phase law.**
  * For every odd n with 2^9 < n < 7,301,503 (3,650,495 starts), the phase φ(n) = frac(log2 n + T*(n)·log2 3), signed, lies in [−0.02502, +0.04566]. The extremes are at n = 3,719 and n = 4,629, and every octave from 2^10 to 2^21 gives the same range. For [2^4, 2^9) the range is [−0.036, +0.054].
  * The principal variation ends 19 → 7 → 5 → 1, contributing +0.00282 to the phase. With ties broken toward the descending move and 20,000 random starts per octave, the share is 94.1% at 2^10, 98.5% at 2^16 and 99.5% at 2^21. Rarer endings are 77 → 29 → 11 → 1 and 149 → 7 → 5 → 1.
  * Therefore T* is confined to an admissible set whose consecutive gaps are 12, 17 or 29, by the three-distance theorem with 12·log2 3 ≈ 19 + 0.0196, the Pythagorean comma. Its parity is locally a function of θ = frac(log2 n).
  * θ alone, binned into 1024 bins, predicts the value with accuracy 1.000 for k <= 10, 0.911 for k = 16, 0.870 for k = 21, and 0.841 on [2^26, 2^27). Adding n mod 16 improves this by at most 0.4%.
  * The profile drifts because T* grows by about 2.06 per octave while 12·log2 3 misses an integer by only 0.0196.

## 4. Sheets and the r <-> −r symmetry

* **Why the symmetry holds (heuristic, not a proof).** By L5 the move pair (d, u) is negation-equivariant with the sheets exchanged, and by L4 the phase of a position is negation-invariant. The value is statistically a function of the phase together with the 2-adic exponent data along the tree, and both are invariant. The target (1 versus −1) enters only through the endgame, whose phase contributions are themselves invariant. So the limiting P-probability of the class r equals that of −r, as observed. The only thing the sheet labels do is decide which sheet carries the descending move: plus on n ≡ 1, minus on n ≡ 3 (mod 4).
* **Exact mirror at mod 4.** The mirror is exact to four digits:
  * n ≡ 1: only the plus child is P 0.0414, only the minus child is P 0.0436;
  * n ≡ 3: the same numbers with the sheets swapped.
* **Pointwise effect.** Beyond densities there is a pointwise mirror effect. A nearby n' ≡ −n (mod 2^m) agrees with n in 83.8% of cases at m = 10, against 81.7% for a control at the same distance; at m = 16 the figures are 79.9% and 77.1%.
* **Relation to the additive-choice zoo.** In the zoo's language, the one-player sign choice is trivial because d always descends. In the two-player game the ascending sheet matters exactly when the descending child is N, that is k >= 1, at 53% of positions.

## 5. Caveats

* All ranges are finite. A draw above 2^32, or an unbounded draw ray (Pochuev's dichotomy), is not excluded.
* Heights from block mode are bounds. The exact values rest on bottleneck-order insertion, validated against the dense exact heights.
* Remoteness values are capped upper bounds. They are stable between the two caps on the stated prefix.
* The phase window and every percentage in §3.2 are observations, not theorems.

## 6. Files

In `04-computation/experiments/`, each `collatz_procgen_20260922_game_*.c` or `.py`:

| file | purpose |
|---|---|
| `dense` | exact heights and frontier |
| `sparse` | exact heights above the cap |
| `blocks` | no-draw verification to 2^32 |
| `cert` | proof-certificate extractor |
| `check.py` | independent Python retrograde, OEIS comparison, certificate checker |
| `struct` | statistics |
| `remote` | remoteness |
| `phase.py` | phase law |

Output: `05-knowledge/results/collatz_procgen_20260922_game.out`. Large intermediate files are in `scratch/procgen_game/` and are not committed.
