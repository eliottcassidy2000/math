# Machines, automata, FRACTRAN and tournaments: one dictionary; the 25 primes below 100 against the square-sum graphs; the prime square-sum path problem

klein-2026-11-01 (worktree `codex/session-ramsey-20261030`). **Status: DICTIONARY (typed, inherits the exact translations of mac-mini-2026-11-01-lenses, mac-mini-2026-11-01-fractranTM and session135) + PROVED elementary (THM-4860: congruence structure of the prime square-sum graph mod 4, 12 and 36; Chebyshev-bias and 2-matching obstructions) + FINITE-EXACT (scans to 10^5, 7·10^5 and 5·10^6; exact max-flow certificates) + CONSTRUCTION, VERIFIED (THM-4860 (f): the first 50,396 primes arranged with consecutive sums squares) + typed verdicts (NUMEROLOGY for the prime groups against the square-sum graphs, inherited and sharpened) + HYP-9451 (open items).** No open problem of the atlas is closed here; the prime square-sum path problem is new in this packet.

## 0. Directive and inheritance

Directive (2026-11-01, verbatim): "try to draw analogy between turing machine and cellular automata and fractran programs and our tournament structure research as you continue seeing open mathematical questions through these lenses, and try to translate seamlessly between them, see how {2,3,5,7,11,13},{17,19,23,31,47}{29,41,59,71},{37,43,53,61,67},{73,79,83,89,97} relates with the first 25 corresponding graphs for the square-sum sequence, and other instances of 25 we have studied in this repo".

Read newest first; everything below cites and does not redo:

- **mac-mini-2026-11-01-lenses**, THM-4857 (the GF(2) automaton of a tournament) and THM-4858 (exact translations of the Goldbach machines), scripts in [`04-computation/experiments/circulant_fractran_goldbach_20261101/`](../../../04-computation/experiments/circulant_fractran_goldbach_20261101/) (pre-registered, all checks PASSED; pushed 2026-11-01 while this note was written, which is why this packet's theorem is THM-4860 and the Ramsey theorem of 2026-10-30 became THM-4859): the GF(2) tournament automaton `x ↦ Ax` (for a circulant `Circ(Z_q,S)` this is the additive cellular automaton with neighbourhood `S`); converse-blind and relabelling-invariant on all classes `n ≤ 6`; circulant rank `q − deg gcd(a(x), x^q − 1)` on 26,039 connection sets; Paley rank `q` for `q ≡ 3 (mod 8)` and `(q+1)/2` for `q ≡ 7 (mod 8)` with images the Hamming `[7,4,3]` and Golay `[23,12,7]` codes; GF(2)-similarity types equal the merged metagraph nodes for `n = 3, 4, 5` (2, 3, 10) and fall short at `n = 6` (24 types for 34 nodes). Literal translations of Leng's 25-state Goldbach machine: a 53-symbol radius-1 cellular automaton in lockstep for 200,000 steps; the 568-fraction, 349-prime FRACTRAN program LENGGAME whose two stacks agree with the tape for the first 3,000 machine steps (24,405,662 FRACTRAN steps); BB(2) control halts with four 1s.
- **mac-mini-2026-11-01-fractranTM**, [`fractran_turing_circulant_goldbach_20261101.md`](../fractran_turing_circulant_goldbach_20261101.md): THM-4855 (circulant Hamiltonian paths as path words with distinct partial sums; pure-blue `K1, C3, T5` are the circulants of order `≤ 5`), THM-4856 (GOLDBACHGAME, FERMATGAME, LEGENDREGAME halt iff a counterexample exists; "kernel plus controller" as the normal form of `Π⁰₁` prime conjectures), and §2.5: three readings of the owner's five prime groups against the first 25 square-sum graphs, an interleaving obstruction (PROVED) and the verdict NUMEROLOGY.
- **session132**, [`session132_fractran_tournament_recursion_20261031/prime_graphs/`](../session132_fractran_tournament_recursion_20261031/prime_graphs/): the rank-indexed and prime-indexed 25-row tables; the owner's earlier list had 71 twice and no 67, with the conditional repair D: 71 → 67. **The 2026-11-01 directive to klein lists 67 in the fourth group and 71 once (third group): the repair is confirmed by the owner's own wording.** All five groups below use the repaired list.
- **session32**, [`session32_square_monster_20261023/PRIME_BLOCK_CARDINALITY_HOSTILE.md`](../session32_square_monster_20261023/PRIME_BLOCK_CARDINALITY_HOSTILE.md): on the fifteen supersingular primes the square-sum graph is a 14-edge tree with leaves 3, 11, 19, 31, 71. [`monster_levels_square_sum_7_21_20261022.md`](../monster_levels_square_sum_7_21_20261022.md): the split `{2..13} / {17,19,23,31,47} / {29,41,59,71}` is `Fi22 / B∖Fi22 / M∖B` (STRUCTURE, Aricheta); every square-sum link to it NUMEROLOGY.
- **session135**, [`session135_fractran_tm_structure_20261101/`](../session135_fractran_tm_structure_20261101/): exact costs of the FRACTRAN ⇄ register-machine ⇄ single-tape translations (INC/DECJZ `2u+3`, GOTO 1), the 14-fraction square recogniser and the square-sum path/cycle certificate machines (25/36 non-halting states), `lrc_machine` (five fractions accept iff `A ≥ B` and `C ≥ D`).
- Canon used: THM-4849 (circulant tournaments are rational Collatz cycles), THM-4852 (PRIMEGAME is a selector quasi-cellular automaton on base-`29#` digits), THM-4854 (metagraph quotient of the complement-move graph), THM-4505/4510/4850/4851 (square-sum graph: degree formula, onsets, ninja pairs at 42, paths iff `N ∈ {1,15,16,17,23} ∪ [25,∞)`, cycles iff `N ≥ 32`), and klein's 2026-10-27 reflection on halting levels.

## 1. The dictionary: four rewriting systems of one shape

Every system is a state space `X`, a finite *test algebra* (what one step may read), finitely many update maps, and a conflict rule.

| system | state | tests (guards) | updates | conflict rule |
|---|---|---|---|---|
| Turing machine | tape, head, control state | the symbol under the head and the state | write, shift, change state | none (a table) |
| cellular automaton | `Σ^Z` | each cell's radius-`r` neighbourhood, all cells at once | one local table | none |
| FRACTRAN | `n ∈ N`, i.e. an exponent vector | `d_i divides n` — read modulo the remaining denominator lcm `L`; one base-`P#` digit suffices only if those denominators divide `P#` | `n ↦ f_i n`, additive on exponents: a register machine | **first applicable fraction = the transitive tournament `TT_k` on the fraction list** |
| generalised Collatz | `n ∈ Z` | `n mod M` | affine `(a_r n + b_r)/d` | residue classes partition `Z` |
| tournament rewriting | a tournament `T` on `[n]`, matrix `K` | arc tests along a vertex word (a Hamiltonian-path guard is `n−1` arc tests) | THM-4810/4854's complement move `K ↦ −K + 2A_Q`; the linear part is mac-mini's `x ↦ Ax` over GF(2) | the choice of path; by Rédei an **odd** number of guards fire in every state |

**Exact translations (what is preserved) and encodings (what is lost).**

- TM ⇄ CA: radius-1 automaton with head cells `(q, g)` and a HALT marker; one step per step (folklore, Smith 1971; THM-4858 (a), mac-mini T1). Exact; the size changes (25 states ↦ 53 symbols).
- TM ⇄ FRACTRAN: two stacks and one-hot control (THM-4858 (b), mac-mini T2; session135's cost table). Exact as a simulation, polynomially slowed: 3,000 machine steps are 24.4 million FRACTRAN steps.
- FRACTRAN ⇄ generalised Collatz: **exact and forgetful of nothing.** With `M = lcm(d_i)`, the branch of `n` is the least `i` with `d_i | n`, a function of `n mod M`, and on that class the map is `n ↦ f_i n` (all `b_r = 0`). Conversely Collatz-type maps with `b_r ≠ 0` enter FRACTRAN only through an encoding (Conway's), which costs a change of state space.
- FRACTRAN ⇄ CA: THM-4852's selector quasi-automaton reads one primorial digit for the squarefree PRIMEGAME instance; general denominator valuations require the corrected guarded precision rule. Carries still make the update non-local, so this arrow is a dictionary, not a strict CA. See [FRACTRAN_PRECISION P1--P5](../stability_priority_moduli_20261102/FRACTRAN_PRECISION.md). The program `(1/4,3/2)` distinguishes inputs2 and8 although both are2 modulo6. This corrects the incoming generic table's prime-support-only wording; the simulations and their receipts are unchanged.
- tournament ⇄ CA: a circulant's GF(2) automaton **is** an additive CA with neighbourhood `S` (THM-4857, mac-mini A1–A4); the translation forgets the converse (`A ∼ Aᵀ`) and stops separating merged nodes at `n = 6`.
- tournament ⇄ Collatz/FRACTRAN: THM-4849 (a circulant is a rational Collatz cycle of shape `((q−1)/2, q)`), THM-4855 (its path moves are the words over `S` with distinct partial sums), session135's compiled tournament streams. Exact on circulants, by construction.
- tournament ⇄ TM: session135's certificate machines verify a vertex word in `n − 1` arc tests; this arrow is one-directional (checking, not searching).

**Two small facts proved here (elementary).**

- (D1) *FRACTRAN's priority is forced, not chosen.* Replace "first applicable" by an arbitrary tournament `R` on the fractions (apply the applicable fraction that beats every other applicable one). The step is defined exactly when the induced subtournament on `App(n) = {i : d_i | n}` has a source. If the denominators are pairwise coprime, every subset of fractions is realised as some `App(n)` (take `n = ∏ d_i` over the subset), so totality forces every induced subtournament to have a source, i.e. `R` transitive: the only total conflict rules are the orders `TT_k`. For PRIMEGAME's denominators (`91 = 7·13` with `13`, `85 = 5·17` with `17`, …) not every subset is realised and the question of non-transitive total rules is left open (HYP-9451 (3)).
- (D2) *Rédei's theorem is a parity invariant of the tournament rewriting system*: in every state an odd number of Hamiltonian-path guards fire (Rédei 1934). FRACTRAN has no analogue: under PRIMEGAME, `n = 1` fires one fraction (`55/1`) and `n = 2` fires two (`15/2`, `55/1`).

## 2. Open problems through the four lenses: the level is the invariant, the size is a coordinate

An exact simulation is a computable bijection of configurations, so it preserves the logical form of a question (its position in the arithmetical hierarchy) and the identity of a counterexample, while it changes every *size*: states, symbols, fractions, step counts (25 states ↔ 53 symbols ↔ 568 fractions; one machine step ↔ about 8,000 FRACTRAN steps on average in mac-mini's run). The 2026-10-27 reflection made this point for busy-beaver sizes; here is the table the four lenses give for the repository's open questions.

| question | form | level | single "halts iff" program? |
|---|---|---|---|
| Goldbach; Legendre; "no Fermat prime beyond `F_4`" | no counterexample | `Π⁰₁` | yes: GOLDBACHGAME, LEGENDREGAME, FERMATGAME (THM-4856); Leng's machine |
| LRC(14) (lonely runner, 14 runners) | no speed vector with gap `< 1/15`; each instance is a finite periodic check | `Π⁰₁` | yes in principle; session135's `lrc_machine` is the comparison kernel |
| square-sum Hamiltonian cycles for all `N ≥ 32` | each `N` decidable | `Π⁰₁`, **closed** (Gerbicz; THM-4851) | — |
| black-core / metagraph connectivity for every order | each order decidable | `Π⁰₁` | yes in principle |
| Collatz (every `n` reaches 1) | `∀n ∃t` | `Π⁰₂` | **no**: it is the totality of one machine, not the halting of one; "weak Collatz as BB(14)" in the reflection is a `Π⁰₁` fragment |
| every prime has a square-sum prime partner (§3.4) | `∀p ∃s (s² − p prime)` | `Π⁰₂` | **no** (same reason) |
| some `N` admits a prime square-sum Hamiltonian path (§3.5) | `∃N [decidable]` | `Σ⁰₁`, **now PROVED** (THM-4860 (f): `N = 617,059`) | yes, and it halts iff the statement is **true** (the mirror image of the Goldbach programs); it halted |
| ETH is false (the 2026 preprint) | `∃δ ∀A …` | `Σ⁰₄` as stated | no |

The dictionary therefore sorts the atlas into three rows (`Π⁰₁` halting, `Π⁰₂` totality, `Σ⁰₁` search), and "translate seamlessly" means: move the size coordinates freely (a 25-state machine, a 568-fraction program and a 53-symbol automaton are one object) and never move the row.

## 3. The 25 primes below 100 and the square-sum graphs

### 3.1 The first 25 graphs (recomputed; agrees with session132's table)

`Q_n` on `{1..n}`, `i ∼ j` iff `i + j` is a square ([`squaresum_graphs.py`](scripts/squaresum_graphs.py), [output](data/squaresum_graphs.out)): edges 0,0,1,1,2,3,4,5,6,7,8,9,11,13,15,16,17,18,20,22,24,26,28,30,32; one component from `n = 14`; Hamiltonian paths at `n = 15, 16, 17, 23, 25` and nowhere else `≤ 25`; no Hamiltonian cycle `≤ 25` (THM-4851; the path onset 25 is a Pell exit/entry event, `monster_levels` §4).

### 3.2 The groups against the graphs: inherited verdict, two threshold coincidences typed

The three indexings (rank `k ↦ Q_k`, prime `p ↦ Q_p`, path-admissible order) and the interleaving obstruction are session132's and mac-mini-fractranTM's; verdict NUMEROLOGY. Two exact statements in the prime indexing are worth recording because they look like structure and are not:

- `Q_p` is disconnected iff `p ∈ {2,3,5,7,11,13}` = the Fi22 group (session132, threshold 14).
- `Q_p` has **no Hamiltonian cycle iff `p ≤ 31`**, and `{2,…,31}` are exactly the eleven primes that are supersingular *and* form an initial segment: 37 is the first non-supersingular prime, while 41, 47, 59, 71 are supersingular with cycles. Both facts hold because a square-sum threshold (14; 32, equivalently THM-4505's minimum degree `≥ 2` iff `n ≥ 31`) falls into a prime gap, `(13,17)` and `(31,37)`, that happens to bound an initial segment of Ogg's list. Ogg's list is defined by divisibility of `|M|`; the thresholds by Pell/degree arithmetic. **COINCIDENCE of two gaps**, no map.

### 3.3 New: the square-sum graph on the 25 primes themselves (THM-4860 (a)–(c))

Vertices the 25 primes below 100, `p ∼ q` iff `p + q` is a square ([`prime_squaresum.py`](scripts/prime_squaresum.py), [output](data/prime_squaresum.out)).

- 24 edges, one giant component of 22 and three isolated primes 37, 43, 67, all non-supersingular; their first prime partners are 107 (`12²`), 101 (`12²`) and 257 (`18²`). Leaves 19, 31, 73, 79, 89.
- Layer by layer: Fi22 `{2,3,5,7,11,13}` is a **perfect matching** `2–7, 3–13, 5–11` (squares 9, 16, 16); `B∖Fi22` is the path `19–17–47` with 23 and 31 internally isolated (session32); `M∖B = {29,41,59,71}` is a **perfect matching** `29–71, 41–59` with both sums equal to 100; the two non-supersingular pentads `{37,43,53,61,67}` and `{73,79,83,89,97}` are **edgeless**.
- **Mechanism (a), bipartiteness.** Odd `p + q` square forces an even square, so `p + q ≡ 0 (mod 4)` and `{p, q} ≡ {1, 3} (mod 4)`; `2 + q` square forces an odd square `≡ 1 (mod 8)`, so `q ≡ 7 (mod 8)`. Hence `A = {2} ∪ {p ≡ 1 (4)}`, `B = {p ≡ 3 (4)}` is a bipartition of the prime square-sum graph for every bound (brute-force check on all 85,035 edges below `10^5`). **Mechanism (b)**: two primes `≡ 1 (mod 3)` are never adjacent (`p + q ≡ 2` is a non-residue mod 3).
- Consequences. `M∖B` has sides `AABB` (`29, 41 | 59, 71`), so any matching pairs `{29,41}` with `{59,71}`; of the two pairings only `29–71, 41–59` gives squares (the other gives 88 and 112). Fi22 has sides `ABABBA`. On the fifteen supersingular primes the sides are 6 (`2,5,13,17,29,41`) against 9, so **no Hamiltonian path exists by imbalance alone**, re-proving session32's leaf count in one line.
- **Null model for the owner's pentads.** Among the ten non-supersingular primes below 100 the only square-sum edge is `61–83` (`144`); 70 of the 126 splits into two pentads leave both halves edgeless. The size split is one of them. NUMEROLOGY.

### 3.4 Growth, isolated primes, and a `Π⁰₂` question

| bound | primes | edges | isolated | components | max degree | `|A| − |B|` |
|---|---|---|---|---|---|---|
| 100 | 25 | 24 | 37, 43, 67 | 4 | 4 | −1 |
| 1,000 | 168 | 300 | 7 (163, 373, 601, 691, 757, 919, 997) | 8 | 11 | −6 |
| 10,000 | 1,229 | 4,519 | 7 (6301, …, 9769) | 8 | 34 | −9 |
| 100,000 | 9,592 | 85,035 | 1 (57649; first partner 118751) | 2 | 100 | −24 |

Every prime below 100 has a partner below 500 (67 waits for 257). "Every prime `p` has a prime `q` with `p + q` a square" is `∀p ∃s`: `s² − p` prime, a Bunyakovsky-type statement for the family `s² − p`, `Π⁰₂`, and I know no proof; it is the no-isolated-vertex property of the infinite prime square-sum graph (HYP-9451 (2)).

### 3.5 New: the prime square-sum path problem and Chebyshev's bias (THM-4860 (d))

Arrange the primes below `N` in a row so that consecutive sums are squares — the square-sum problem (THM-4851) restricted to primes. Three necessary conditions, from §3.3: no isolated prime; at most two leaves; and, since the graph is bipartite, `| |A| − |B| | ≤ 1`, i.e. **Chebyshev's bias `π(N;3,4) − π(N;1,4)` must be 0, 1 or 2** (a Hamiltonian cycle needs bias exactly 1 and minimum degree 2). [`prime_path_scan.py`](scripts/prime_path_scan.py) ([output](data/prime_path_scan.out)) tracks the three conditions incrementally over all prime bounds `N < 700,000`:

- the bias is in `{0,1,2}` at 400 bounds, but **no bound below 615,731 passes all three**; the first Chebyshev window (around 26,861) is killed by isolated and leaf primes;
- the 323 passing bounds are exactly `N ∈ [615,731, 634,199]`, the second Chebyshev window (Leech's region 616,841–633,798 widened by the allowance of bias 1 and 2); at `N = 615,731` there are 50,283 primes, 930,345 edges, one leaf (582,451, whose only partner below 700,000 is 362,333), four primes of degree 2 and `|A| − |B| = −1`;
- **no bound below 700,000 admits a Hamiltonian cycle** (a leaf prime is present at every bias-1 bound).

So the prime square-sum problem has no solution for any `N < 615,731`, by counting alone (THM-4860 (d), FINITE-EXACT).

**The mod-12 refinement (THM-4860 (d′)).** Mechanism (b) splits each side: class 1 (`p ≡ 1 (mod 12)`, side A) is adjacent only to class 11′ `= {p ≡ 11 (mod 12)} ∪ {3}`, and class 7 (`p ≡ 7 (mod 12)`, side B) only to class 5′ `= {p ≡ 5 (mod 12)} ∪ {2}` (checked on all 85,035 edges below `10^5`). Counting path edges at a restricted class against the capacity of its partner class gives two more necessary conditions, `c_1 ≤ c_11 + 2` and `c_7 ≤ c_5 + 2` (`c_r` = primes `≡ r (mod 12)` below `N`), i.e. two further Chebyshev-type biases must have the right sign. Since 1 is the only square class mod 12, `c_1` is the Chebyshev-depleted class and `c_1 ≤ c_11 + 2` holds generically, while `c_7` against `c_5` is an unbiased fluctuation. [`prime_path_mod12.py`](scripts/prime_path_mod12.py) ([output](data/prime_path_mod12.out)): the three count conditions together hold on the bounds in `(100, 27,031] ∪ [615,731, 634,199]` and nowhere else below `5·10^6`; combined with the isolated/leaf conditions, the 323 admissible bounds below 700,000 are unchanged (at `N = 615,731`: `c_1 − c_11 = −44`, `c_7 − c_5 = −43`, comfortably inside). The refinement does not exclude the window, but it explains the solvers' behaviour: every prime left stranded by the rotation-extension searches (sixty-odd in the diagnostic runs) lay in classes 1 and 7, the classes with half-size neighbourhoods, and a stranded class-1 prime can only be exposed from a class-11 endpoint (class 7 from class 5), while the random rotation walk drifts to low-degree endpoints, which are class 1 and 7 themselves.

**The 2-matching obstruction and the Hall certificate (THM-4860 (e)).** A Hamiltonian path has `n − 1` edges and all degrees `≤ 2`, so `ν₂(G_N) ≥ n − 1`, where `ν₂` is the maximum size of a spanning subgraph with all degrees `≤ 2`; for a bipartite graph `ν₂` is an exact maximum flow (source → A with capacity `min(deg, 2)`, edges with capacity 1, B → sink with capacity `min(deg, 2)`; [`prime_path_2matching.py`](scripts/prime_path_2matching.py)). At `N = 615,731`: `ν₂ = 50,274 = n − 9`, so **no Hamiltonian path exists at the first count-admissible bound**. The minimum cut is a human-readable Hall certificate: `S` = all 4,183 primes `≡ 25 (mod 36)` below `N`, whose neighbourhood is the 4,178 primes `≡ 11 (mod 36)` plus the prime 3, capacity `2·4,179 = 8,358 < 8,366 = 2|S|`.

**The mod-36 partner law (THM-4860 (d″)).** The certificate is explained by one more congruence: if `p ≡ 1 (mod 12)` is adjacent to `q`, then `p + q ≡ 0 (mod 12)`, so the square `p + q` is divisible by 36 and `q ≡ −p (mod 36)`. The three classes `1, 13, 25 (mod 36)` of primes `≡ 1 (mod 12)` are adjacent only to the classes `35, 23, 11 (mod 36)` respectively (and to 3), and the classes `7, 19, 31` only to `29, 17, 5` (and to 2); the remaining edges (classes 5 and 11 mod 12) use squares `≡ 4, 16, 28 (mod 36)` and mix the sub-classes. So six Chebyshev-type pair conditions `c_a ≤ c_{a′} + 2` are necessary ([`prime_path_mod36.py`](scripts/prime_path_mod36.py), [output](data/prime_path_mod36.out); the law checked on all edges below `10^5`). At `N = 615,731` the pair `(25, 11)` has `c_25 − c_11 = 5`; the pairs `(1, 35)`, `(13, 23)` enjoy the bias (1, 13, 25 are the square classes mod 36), the class-7 pairs are fair fluctuations. Below `5·10^6` the mod-4 balance and the six pair conditions hold only on `(100, 27,031] ∪ [617,039, 634,199]`. The exact flow test is finer than the counts (it sees sub-sets and the shared slots of 2 and 3): over the window it excludes 21 of the first 22 count-admissible bounds (`615,731 … 617,051`, gaps `−8 … −1`), the first survivor being `N = 617,053`; see [`prime_path_2matching_window.log`](data/prime_path_2matching_window.log) for the whole window.

**Why the greedy solvers failed.** Five rotation-extension solvers ([`prime_path_posa.py`](scripts/prime_path_posa.py) … [`prime_path_posa5.py`](scripts/prime_path_posa5.py): Warnsdorff extension, forced chains, insertion moves, depth-2 rotation lookahead, class steering) stalled 10–60 primes short at every bound; the stranded primes were always in classes 1 and 7 (mod 12). The two restricted classes leave only `c_{11′} − c_1 ≈ 45` and `c_{5′} − c_7 ≈ 44` free slots, so a Hamiltonian path is essentially a gluing of about 45 alternating segments of each restricted type, and a walk that spends class-11 or class-5 slots on unrestricted neighbours can never recover them. The degree-constrained flow respects every slot globally.

**EXISTENCE (THM-4860 (f)).** [`prime_path_flow.py`](scripts/prime_path_flow.py) takes the flow's maximum 2-matching (one path plus cycles) and merges the cycles into the path by 2-opt exchanges (remove `(c, c⁺)` and `(p, p⁺)`, add `(c, p)` and `(c⁺, p⁺)`, or `(p, c⁺)` and `(c, p⁺)`, successor and predecessor variants) and 3-opt insertions (a cycle between two path cuts with the middle segment reversed), with randomised restarts. At `N = 617,059` it succeeds in 10 seconds: **the first 50,396 primes, 2 through 617,053, can be arranged in a row so that every one of the 50,395 consecutive sums is a perfect square.** The row starts `582,451, 362,333, 293,767, 381,917, 313,639, 391,961, …` (the leaf first) and ends `…, 244,451, 256,813, 351,587, 616,669`; it uses 245 distinct square roots from 6 to 1,110, the root 780 serving 1,021 times. The file [`data/prime_path_N617059.txt`](data/prime_path_N617059.txt) (SHA-256 `f6f520ee0f74371834b9c0de533232764febd35158649e65c1de2243e0f00491`) is checked by the independent [`prime_path_verify.py`](scripts/prime_path_verify.py) ([output](data/prime_path_verify_617059.out)), and the deterministic solver (seed 1) reproduces it byte for byte under `python3` and `python3 -O`. **Minimality (THM-4860 (e′), (f)).** Every bound below `617,053` is excluded ((d), (d′), (e)). The bound `617,053` passes the flow test, yet its 2-matching always leaves the whole `(25, 11)` piece as one unmergeable cycle of about 8,380 vertices, and that is a theorem, not bad luck: inside `X_a` = class `a` ∪ its partner class ∪ the shared prime (3 or 2), the path edges form a forest with at least `2c_a − e_a` edges on `c_a + c_{a′} + 1` vertices, so **`c_a − c_{a′} ≤ e_a`**, the number of global endpoints in class `a`; with `|B| = |A| + 1` both endpoints lie in B and `e_a = 0` for the A-classes. At `N = 617,053`: `|B| = |A| + 1` and `c_25 − c_11 = 1`, contradiction. (`N = 617,039` passes this count but has `|A| = |B| + 1`, which puts both endpoints in A while the leaf 582,451 lies in B: `ν₂ = n − 2`.) Hence **50,396 is exactly the least number of initial primes that can be arranged with consecutive square sums**, and the solution's structure is forced: `|A| = |B|`, the A-endpoint must be a class-25 prime (it is: `616,669 ≡ 25 (mod 36)`), and the `(25, 11)` piece is a single segment threaded through the prime 3. Below `5·10^6` the segment conditions pass only on `(100, 1,423] ∪ [26,821, 27,031] ∪ [617,039, 634,199]` ([`prime_path_mod36.py`](scripts/prime_path_mod36.py) (iv)); the exact flow test over the second Chebyshev window: of the 323 count-admissible bounds in `[615,731, 634,199]`, 118 are excluded by `ν₂ < n − 1` (largest deficiency 8) and 205 pass; the segment condition (e′) then removes those with `|B| = |A| + 1` and a positive class surplus, and the first survivor of everything is `N = 617,059`.

In the level table of §2 the `Σ⁰₁` row has thus been settled by its own method: the search halted.

## 4. The repository's twenty-fives

The typed ledger is [`LEDGER_15_25_MULTIPERFECT.md` §4b](../analogy_reduction_atlas_20261025/LEDGER_15_25_MULTIPERFECT.md) (klein, 2026-10-25); its row on the square-sum onsets now cites the monster-levels Pell law instead of a reflection that never mentions Pell (citation corrected in place). Rows added by the last week, typed:

| twenty-five | where | type |
|---|---|---|
| `π(100) = 25`: the owner's five groups (15 supersingular + 10 others) | §3, session131/132, fractranTM §2.5 | NUMEROLOGY against the square-sum graphs; STRUCTURE as Ogg's list |
| `R(4,5) = 25` as an input of the `R(3,4,5)` recursion `25 + 57 + 77 − 1 = 158` | THM-4859, `ramsey_cascade_barrier…/README.md` §3 | CITED exact value; COINCIDENCE with every other 25 |
| Leng's 25 states ↔ session135's 25 non-halting states of the square-sum path certificate machine | mac-mini T1–T2, session135 §3 | COINCIDENCE of two description sizes (the dictionary says sizes are coordinates: 25 ↔ 53 ↔ 568) |
| strong floor `f(7) = 25` and the path onset 25 | THM-1370, THM-4851, monster-levels §4 | STRUCTURE each; `(15, 25)` against `[15, 25]` NUMEROLOGY (ledger) |
| the two regular groups of order 25 | THM-4764 | STRUCTURE |

Flags from the inventory, not acted on: MISTAKE-694 is used for two entries in `01-canon/MISTAKES.md`; THM-1370 is shared by two canon files (already noted in the repository).

## 5. Verdicts

- Dictionary: the four systems are one shape (tests, updates, conflict rule); exact arrows preserve level and counterexample and move sizes; FRACTRAN's priority is the transitive tournament and is forced by totality on coprime denominators (D1); Rédei's parity is a rewriting invariant without a FRACTRAN analogue (D2). DICTIONARY, not a mechanism between the open problems.
- The owner's groups against the square-sum graphs: NUMEROLOGY (inherited), with the two threshold coincidences typed as prime-gap accidents; on the primes themselves the structure is explained by residues mod 4 and mod 3 (THM-4860 (a)–(c)).
- New object: the prime square-sum path problem, governed by Chebyshev-type biases mod 4 and mod 36, by an exact 2-matching (Hall) test and by the segment condition `c_a − c_{a′} ≤ e_a`; **the least number of initial primes that can be arranged with consecutive square sums is exactly 50,396** (THM-4860 (d)–(f), verified file). Open: Hamiltonian cycles, infinitely many bounds (HYP-9451).

## 6. Reproduction

```bash
cd 05-knowledge/results/machines_automata_fractran_tournaments_dictionary_20261101
python3 scripts/squaresum_graphs.py          # Q_1..Q_25 and the group check
python3 scripts/prime_squaresum.py           # §3.3–3.4 (about 1 min; brute-force mechanisms to 10^5)
python3 scripts/prime_path_scan.py 700000    # §3.5 scan (about 10 s)
python3 scripts/prime_path_mod12.py 5000000  # §3.5 mod-12 refinement (about 1 min)
python3 scripts/prime_path_2matching.py 615731           # §3.5 exact exclusion with Hall certificate (about 10 s; needs scipy)
python3 scripts/prime_path_2matching.py window 615731 634199  # the whole window (about 25 min)
python3 scripts/prime_path_flow.py 617059 1                 # §3.5 construction of the path (about 10 s; deterministic)
python3 scripts/prime_path_verify.py data/prime_path_N617059.txt 617059   # independent check
shasum -a 256 -c SHA256SUMS
```

Receipts in [`receipts/`](receipts/) are byte-identical under `python3` and `python3 -O`.
