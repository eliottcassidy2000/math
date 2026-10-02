# Shaved tournaments: the oriented graphs inside every n-tournament

**opus S15, 2026-10-01.** Worktree `collatz-functional-uniqueness-20261001`; this answers the tournament part
of the owner's fourth prompt of the day.

- Script: `04-computation/experiments/shaved_tournaments_20261001.py` (+ `.out`, ALL CHECKS PASSED).
- C helpers: `shaved_tournaments_20261001_closing.c` (the n = 9 check of Theorem A) and
  `shaved_tournaments_20261001_u8.c` (the n = 8 search of Theorem B).
- Canon: THM-4526.

**Status.**
- PROVED:
  - Theorem A (all n). This is a **classical theorem** (Grünbaum 1971; §1), proved independently here.
  - Theorem C, which is **superseded** by Linial–Saks–Sós (1983; §3).
- PROVED by exhaustive computation: Theorem B (`u(n)` for `n ≤ 8`, hence `κ(7) = 12` and `κ(8) = 17`).
- FINITE-EXACT: the tables, the n ≤ 6 parity, and perfect shavings (checked for n ≤ 60).
- The odd case of Theorem A was found by a proof-search subagent of this session. The orchestrator checked
  it line by line, and the script re-checks it with an independent implementation.
- Independent audit: DONE (2026-10-01, blind re-derivation with the auditor's own code; n ≤ 9 classes, every
  search redone). No claim is unsound. The corrections are applied: the two attributions above, five minimum
  certificates at `n = 7` instead of one, and softened wording (MISTAKE-555). The record is in §8.
- Continued in [`shaved_tournaments_proved_vs_conjectured_20261001.md`](shaved_tournaments_proved_vs_conjectured_20261001.md)
  (thread session thread-2cjcob, 2026-10-01; independently audited): `u(9) = 14`, so `κ(9) = 22` (D70 resolved); the excess
  law fails for every `n ≥ 21` with `n ≡ 0, 2 (mod 3)`; a ledger of what is proved and what is conjectured. One
  sentence of §5 is corrected (MISTAKE-557).

## 0. The owner's object

> for any tournament of size 4, no matter the isomorphism class, it contains edges from A to B, B to C,
> C to D, and A to D. Consider this as a shaved down tournament as fundamental object and study it as the
> size of the tournament/graph increases.

Call this oriented graph `H_4`: the Hamiltonian path `A→B→C→D` plus the arc `A→D` from its first vertex to its
last. Two pairs are missing, `A–C` and `B–D`. Check T2 shows:
- **The 4 completions of `H_4` are the 4 classes of 4-tournaments, each exactly once**: TT4, source + 3-cycle,
  3-cycle + sink, and the strong one. So `H_4` is a perfect shaving: delete two arcs and every class is
  recovered exactly once.
- **The number of copies of `H_4` in a 4-tournament `T` equals `|Aut T|`**: 1, 3, 3, 1. This is the case
  `m = 1` of a general identity. For a fixed labelled shaving `S`, the number of completions isomorphic to
  `T` is `m_S(T) = emb(S,T)/|Aut T|` (orbit counting).

**Definition.** A *shaving* of order `n` (an *`n`-unavoidable digraph* in the literature) is a spanning oriented
graph `S` on `n` vertices that is contained in every `n`-tournament. Equivalently, the `2^k` completions of `S` (`k` = number of missing pairs) meet every
isomorphism class.

This is the repo's flip-rank / transversal-subcube invariant (HYP-3798, mac-mini; seeded by the owner's
"flip only 2 arcs" at `n = 4`), seen from the fixed side:
- `u(n)` := the maximum number of arcs of a shaving;
- `u(n) = C(n,2) − κ(n)`.

There are two ways to let `n` grow:
- **(I)** keep the shape: `H_n` = Hamiltonian path + the arc from the first to the last vertex;
- **(II)** keep the property: the largest shavings at each `n`.

## 1. The shape: path plus first→last arc

A Hamiltonian path (HP) either **closes** (last vertex → first, so path + arc is a Hamiltonian cycle) or is
**non-closing** (first → last). Each Hamiltonian cycle gives exactly `n` closing HPs, and `Aut(H_n) = 1`, so

`#copies of H_n in T = #non-closing HPs = H(T) − n·hc(T).`

**Theorem A (Grünbaum 1971; independent proof below).** Let `n ≥ 2`. An `n`-tournament contains no copy of `H_n` iff it is the 3-cycle `C3` or the
5-vertex tournament `C3[1,C3,1]` (vertex `a` beats a 3-cycle, the 3-cycle beats `b`, and `b → a`). When `n`
is even, every `n`-tournament contains an **odd** number of copies.

Equivalently: **the only tournaments in which every Hamiltonian path closes into a Hamiltonian cycle are
`C3` and `C3[1,C3,1]`.** So `H_n` is a shaving iff `n ∉ {3, 5}`.

*Proof.*

*Even `n`.* By Rédei (1934), `H(T)` is odd, and `n·hc(T)` is even.

*`T` not strong.* Take the strong components `C_1 ⇒ … ⇒ C_k` (`k ≥ 2`). Concatenating Hamiltonian paths of
the components gives an HP from `C_1` to `C_k`, and `C_1 ⇒ C_k`.

*Odd `n`, `T` strong.* Suppose `T` has no non-closing HP.

**Lemma 1 (insertion).** Let `Q = (q_1, …, q_r)` be a non-closing HP of `T − w`. Write `σ_k = o` if
`w → q_k` and `σ_k = i` otherwise. Then `σ = o^j i^(r−j)` with `1 ≤ j ≤ r−1`.

*Proof of Lemma 1.* Each of the following would be a non-closing HP of `T`:
- if `σ_1 = σ_r = o`: the path `(w, Q)`;
- if `σ_1 = σ_r = i`: the path `(Q, w)`;
- if some `σ_k = i` and `σ_{k+1} = o`: insert `w` between `q_k` and `q_{k+1}` (the ends stay `q_1 → q_r`).

A word with none of these three features has the stated form.

**Lemma 2.** A strong tournament on `p ≥ 4` vertices has a vertex whose deletion leaves it strong. This is
classical; it also follows from Moon's vertex-pancyclicity theorem.

Fix a Hamiltonian cycle `c_0 … c_{n−1}` (Camion), with indices mod `n`. Call position `i` **forward** if
`c_{i−1} → c_{i+1}` and **backward** otherwise.

1. *No two consecutive positions are forward.* If `i` and `i+1` were both forward, then
   `(c_i, c_{i+2}, c_{i+3}, …, c_{i−1}, c_{i+1})` would be an HP, non-closing because `c_i → c_{i+1}`.
2. *Two consecutive positions are backward.* Since `n` is odd, the cyclic forward/backward word cannot
   alternate. So some `i` and `i+1` are both backward: `c_{i+1} → c_{i−1}` and `c_{i+2} → c_i`.
3. Let `R = V ∖ {c_i, c_{i+1}}`. Then `Q_1 = (c_{i+1}, …, c_{i−1})` is a non-closing HP of `T − c_i`, and
   `σ_2 = i` there. Lemma 1 forces `σ = o i^(n−2)`, i.e. `R ⇒ c_i`.
4. `Q_2 = (c_{i+2}, …, c_{i−1}, c_i)` is a non-closing HP of `T − c_{i+1}`, and its entry at `c_{i−1}` is
   `o`. Lemma 1 forces `σ = o^(n−2) i`, i.e. `c_{i+1} ⇒ R`.
5. Put `a = c_{i+1}` and `b = c_i`, so `a ⇒ R ⇒ b → a`, with `|R| = n−2` odd. For `s ∈ R` and any HP `Y` of
   `T[R−s]`, the path `(s, b, a, Y)` is an HP, and it is non-closing iff `s → last(Y)`.
   - **`T[R]` not strong** (and `|R| ≥ 2`): take `s` in the initial strong component `R_1`. Every HP of
     `T[R−s]` ends in its terminal component, which lies in `R ∖ R_1`, and `s ⇒ R ∖ R_1`.
   - **`T[R]` strong and `|R| ≥ 4`**: take `s` from Lemma 2. Then `s` has an out-neighbour `t` in `R−s`, and
     cutting a Hamiltonian cycle of `T[R−s]` just after `t` gives a `Y` ending at `t`.

   What is left is `|R| = 1` (`T = C3`) or `T[R] = C3` (`T = C3[1,C3,1]`). ∎

**Remarks.**
- Oddness is used only in Step 2. For even `n`, Rédei's parity does the work.
- The proof is constructive: Check T3b runs Steps 1–5 from every Hamiltonian cycle of every strong class with
  `n = 3, 5, 7` (2529 runs). It fails exactly at `C3` and `C3[1,C3,1]`, and it succeeds on 975 random strong
  tournaments with odd `7 ≤ n ≤ 31`, including the `C3[1,R,1]` family that Step 5 handles.
- Independent exhaustive data (Check T3):
  - no class avoids `H_n` for `n = 4, 6, 7, 8`;
  - every class with `n = 4, 6, 8` has an odd count;
  - `n = 9`: all 1,761,280 one-vertex extensions of the 6880 classes on 8 vertices (this covers every
    9-tournament) contain `H_9`.
- For odd `n` the parity of the count is `1 + hc(T) (mod 2)`, which varies with `T`.
- `C3[1,C3,1]` has scores `(1,2,2,2,3)`, `|Aut| = 3`, `H = 15`, `hc = 3`: all 15 of its HPs close. It is the
  3-cycle with one vertex blown up into a 3-cycle. Blowing up two vertices already breaks this. In
  `C3[X, x, Y]` (`X ⇒ x ⇒ Y ⇒ X`), any split `X = X_1 ⊔ X_2`, `Y = Y_1 ⊔ Y_2` into non-empty parts gives the
  non-closing path `(Y_1, X_1, x, Y_2, X_2)`, so this works once `|X|, |Y| ≥ 2`.
- **Theorem A is classical.** `H_n` is the oriented cycle `D(n,2)`: a directed `n`-cycle with one arc reversed,
  also called a Hamiltonian bypass. Grünbaum (J. Combin. Theory Ser. B 11 (1971) 249–257) proved that every
  tournament of order `n ≥ 3` contains it unless it is `C3` or `C3[1,C3,1]`; see Havet, JCTB 80 (2000) §1.2.
  Thomassen (1980) showed that a strong tournament has at least `n − 5` copies. Benhocine–Wojda, J. Graph
  Theory 7 (1983) 469–473, treat all `D(n,p)`. The proof above is independent; the formula `H − n·hc` and the
  even-`n` parity are at most minor additions. (The first draft said "we found no reference"; see
  MISTAKE-555.)

## 2. The property: the largest shavings

**Theorem B (exhaustive).** `u(n) = 1, 2, 4, 6, 8, 9, 11` for `n = 2, …, 8`, so `κ(n) = 0, 1, 2, 4, 7, 12, 17`.

**(a) `n ≤ 6`.** The maximum shaving is unique up to isomorphism, and it has a unique linear extension:
**the Hamiltonian path plus every span-3 arc `i → i+3`**.
- At `n = 4` this is the owner's `H_4`.
- All spans are odd, so the underlying graphs are bipartite: `P_3`, `C_4 = K_{2,2}`, `K_{3,2}`, `K_{3,3} − e`.
- This is mac-mini's "path + skip-2 diagonal" rule (HYP-3798), now with uniqueness.
- *Parity curiosity (FINITE-EXACT, unexplained).* For `n ≤ 6`, the number of embeddings of
  path + span-3 arcs into **every** tournament is odd. So every class is hit an odd number of times, and
  unavoidability follows from parity, as it does for `H_4`. The multiplicity histograms are `{1}` (`n = 4`),
  `{1, 3}` (`n = 5`) and `{1, 3, 5}` (`n = 6`). At `n = 7` the count is even for 174 of the 456 classes.

**(b) `n = 7`: the Paley heptagon decides.**
- `P_7` (QR_7, `|Aut| = 21`) is the only 7-tournament with no transitive subtournament on 4 vertices
  (classical; re-checked in T5). So every 7-vertex shaving embeds in `P_7`, and the search can be confined to
  subgraphs of `P_7`.
- Under `Aut(P_7)`, the 10-arc subgraphs fall into 16796 orbits, 1792 of them acyclic. **None is a
  shaving**, so `u(7) ≤ 9`.
- Among 9-arc subgraphs, 54 orbits are shavings. So `u(7) = 9` and **`κ(7) = 12`**. This upgrades HYP-3805
  (opus, "very likely") to proved and confirms that the lazy-caterer formula `1 + C(n−2,2)` of HYP-3798 holds
  exactly for `3 ≤ n ≤ 6` (at `n = 2` it gives 1, while `κ(2) = 0`).
- The maximum shavings at `n = 7` form **51 isomorphism classes**. Only 4 of them contain a Hamiltonian path
  (unique linear extension). They are exactly path + span-3 arcs minus one of the four span-3 arcs.
  HYP-3805's span 1 + span 3 − `(0,3)` is one of them.
- The natural continuation, path + all span-3 arcs (10 arcs), is avoided by exactly two classes: `P_7` and the
  `|Aut| = 1` class with scores `(1,2,2,3,4,4,5)`. This matches HYP-3805.
- **Certificate.** To kill every 10-arc candidate you need `P_7` plus 6 more classes, and 6 is the minimum
  (exact branch and bound). The minimum is **not unique**: there are exactly 5 minimum certificates. They
  share five classes (`|Aut| = 9, 9, 7, 3, 3`) and differ in the sixth (four choices with `|Aut| = 3`, one with
  `|Aut| = 5`). So their `|Aut|` multisets are `(9,9,7,3,3,3)` and `(9,9,7,5,3,3)`. Every minimum certificate
  consists of symmetric classes, and `P_7` is needed in each.
- The most frequent blocker is the other vertex-transitive 7-tournament (`|Aut| = 7`, the circulant on
  `{1,2,3}`; it blocks 1417 of 1792 candidates). This is consistent with HYP-3805's heuristic that a class
  with `n!/|Aut|` labellings is hit by few completions. The 13 most frequent blockers are symmetric, but the
  14th has `|Aut| = 1`, so this is a tendency, not an exact law.

**(c) `n = 8`: the excess law's prediction holds.**
- Every acyclic graph on 8 vertices has a topological order, so it suffices to test every set of forward pairs
  `i < j`. The C helper does this for all `C(28,12) = 30,421,755` sets of 12 arcs: **none is a shaving**, so
  `u(8) ≤ 11`.
- Of the `C(28,11) = 21,474,180` sets of 11 arcs, 48,571 are shavings. Up to isomorphism these are **1617
  classes**, each counted exactly once per linear extension (a consistency check). So **`u(8) = 11` and
  `κ(8) = 17`**.
- This is the value klein's HYP-3819 predicted from the excess law
  `excess(8) = κ(8) − ⌈log₂ 6880⌉ = 4` (HYP-3821), which that file called infeasible to verify.
- Only 5 of the 1617 classes contain a Hamiltonian path. Writing the path as `0 → 1 → … → 7`, they add the arcs
  - `{03, 14, 25, 47}`, `{03, 14, 36, 47}` and `{03, 25, 36, 47}`: path + all span-3 arcs minus one interior
    span-3 arc;
  - `{03, 07, 25, 47}`: the owner's `H_8` plus three span-3 arcs;
  - `{03, 07, 16, 47}`.
- Path + all five span-3 arcs (12 arcs) is not a shaving.
- 27 of the 1617 classes have a bipartite underlying graph.

So the unique bipartite pattern "path + span-3 arcs" that governs `n ≤ 6` keeps an echo at `n = 8` (drop one
interior span-3 arc). The number of maximum shavings explodes: `1, 1, 1, 51, 1617` for `n = 4, …, 8`.

## 3. Growth: `u(n) ~ n log₂ n` (Linial–Saks–Sós 1983)

**Theorem C (superseded; see below).** `(1/2 − o(1)) · n log₂ n ≤ u(n) ≤ C(n,2) − log₂ T(n) ≤ log₂ n!`, where `T(n)` is the number of
classes (A000568).

*Proof.*

*Upper bound.* The completions must cover the classes: `2^(C(n,2)−u) ≥ T(n) ≥ 2^C(n,2)/n!`.

*Lower bound.* Two facts:
- **Superadditivity:** `u(a+b) ≥ u(a) + u(b)`. Put two shavings on disjoint vertex sets: any `a` of the
  vertices carry the first, and the rest carry the second.
- **Transitive blocks:** `u(n) ≥ C(k,2) + u(n−k)` whenever every `n`-tournament contains `TT_k`, which holds
  for `k = 1 + ⌊log₂ n⌋` (Erdős–Moser).

Peeling off greedy transitive blocks of size `≈ log₂ n` covers all but `o(n)` vertices with blocks of size
`(1 − o(1)) log₂ n`. That gives `(1/2 − o(1)) n log₂ n` arcs. ∎

**Theorem C is weaker than Linial–Saks–Sós** ("Largest digraphs contained in all n-tournaments",
Combinatorica 3 (1983) 101–104). They prove `n log₂ n − c₁n ≥ u(n) ≥ n log₂ n − c₂ n log log n`. Their `f(n)`
is our `u(n)`, and their upper bound is the same counting argument. So `u(n) ~ n log₂ n`. Their lower bound
uses complete bipartite blocks instead of transitive ones, which is why the transitive peeling here only
reaches `1/2`.

The recursive bound `L(n)` (best of the two rules, seeded with `u(1..7)`), compared with the upper bound
`U(n) = C(n,2) − ⌈log₂ T(n)⌉`:

| n | L(n) | U(n) | 2n−4 | (n/2) log₂ n | log₂ n! |
|---|---|---|---|---|---|
| 7 | 9 (= u) | 12 | 10 | 9.8 | 12.3 |
| 8 | 10 (true value 11) | 15 | 12 | 12.0 | 15.3 |
| 16 | 25 | 44 | 28 | 32.0 | 44.3 |
| 40 | 78 | 159 | 76 | 106.4 | 159.2 |
| 60 | 126 | 272 | 116 | 177.2 | 272.1 |

**Consequence.** The lazy-caterer formula `κ(n) = 1 + C(n−2,2)` (equivalently `u = 2n − 4`) fails **in both
directions**:
- `κ(7) = 12 > 11` and `κ(8) = 17 > 16`;
- `L(n) > 2n − 4` for every `38 ≤ n ≤ 60`, so `κ(n) < 1 + C(n−2,2)` there and, by Theorem C, for all large `n`.
  (An induction with the transitive-block rule extends this to every `n ≥ 38`; see the continuation note §4.)

The counting bound is not tight even at `n = 7` (12 versus 9).

## 4. Perfect shavings end at `n = 4`

A perfect shaving (each class exactly once) needs `T(n)` to be a power of 2. That holds for `n ≤ 4`
(`1, 1, 2, 4`) and fails for `5 ≤ n ≤ 60` (Davis–Pólya formula, Check T7). So the owner's `H_4` is the last
perfect shaving in that range, and conjecturally the last one overall. The maximum shavings at `n = 5, 6`
already hit some classes 3 or 5 times (T8).

## 5. How this meets the Paley and Collatz threads (typed)

- **STRUCTURAL (proved here).** The Paley heptagon is the principal obstruction at `n = 7`: it confines every
  shaving to its subgraphs, and six further classes are needed to kill the 10-arc ones. Lacking `TT_4` confines
  every shaving to its subgraphs. Its maximal symmetry (`|Aut| = 21`) makes it the rarest class (240
  labellings). The rest of the minimum certificate is symmetric too. The same tournament is the extremal object
  of the flip-rank thread (HYP-3805), which is the same problem. *(Corrected 2026-10-01, MISTAKE-557: the draft
  added "and of the LRC thread (HYP-3802)". HYP-3802 only attaches the Paley orientation to the seven roots of
  `−1` of the tight core `{1, …, 13}`, and `{1,2,4}` is not a tight lonely-runner set.)*
- **DICTIONARY, no implication.** `P_7`'s out-neighbourhoods `x + {1,2,4}` are the 7 Fano lines. The set
  `{1,2,4}` is also the parity code of the trivial Collatz cycle (an explained coincidence, nineteenth note), and
  it is the single colour line of a 3-edge-colouring (snarks; twentieth note). None of this transfers
  information between the problems.
- **Selfie tournaments (THM-4524, mac-mini) are a different notion.** There, "shaving" deletes arcs while
  keeping the number of HPs odd. Theorem A is the Rédei parity seen from the endpoints of the paths.

## 6. Directions

- **D69.** Explain the `n ≤ 6` parity of path + span-3 embeddings (a Rédei-type argument?), or show why it
  must stop at 7.
- **D70.** Find `u(9)` (`κ(9)`; the excess law predicts `⌈log₂ 191536⌉ + #{SC classes with |Aut| > 9}`). Which
  9-tournaments obstruct? (The forward-pair search has `C(36,e)` candidates, which needs a smarter search.)
  **Resolved (2026-10-01):** `u(9) = 14`, `κ(9) = 22`, as the excess law predicted. There are 54 classes, all
  rigid, none with a Hamiltonian path. They were found by an orderly search inside a host tournament (continuation
  note §2), and independently by the shave4 lane ([THM-4533](../../01-canon/theorems/THM-4533-redei-graphs-parity-of-shaved-tournaments-and-u9.md)).
- **D71.** Known: `c = 1` (Linial–Saks–Sós 1983). Open: the second-order term, between `−c₁n` and
  `−c₂ n log log n`.
- **D72.** For which spanning oriented graphs is the number of embeddings odd in every tournament? The
  Hamiltonian path (Rédei), `H_n` for even `n`, and path + span-3 arcs for `n ≤ 6` are examples. Compare
  Forcade (Discrete Math. 6 (1973) 115–118) and El Sahili–Abi Aad (Discrete Math. 343 (2020) 111695): every
  antisymmetric Hamiltonian path type occurs an odd number of times in every tournament. Among the maximum
  shavings, the odd-everywhere ones are 2 of 51 at `n = 7`, and none of the 1617 at `n = 8` or the 54 at `n = 9`
  (continuation note §2).

## 7. Reproduction

```bash
python 04-computation/experiments/shaved_tournaments_20261001.py
```

Allow about 15 minutes. The `n = 9` (Theorem A) and `n = 8` (Theorem B) parts compile and run the C helpers with
`gcc -O2 -fopenmp`; without gcc they are skipped and the recorded results are printed.

## 8. Audit record (2026-10-01)

A blind auditor subagent re-derived every claim with its own C and Python code. The record is in the session
scratchpad, `audit20/AUDIT.md`.
- **Classes.** It generated all classes up to `n = 9` itself: 191536 at `n = 9`, with `Σ n!/|Aut| = 2^C(n,2)`.
- **Theorem A.** It counted non-closing paths directly for every class with `n ≤ 9`. It re-implemented Steps
  1–5 from every admissible position (2529 cycle runs, every branch exercised, failing only at the two
  exceptions).
- **Theorem B.** It confirmed `u(7) = 9` by a second route that does not use `P_7`: none of the 352716 forward
  10-sets is a shaving, and 387 forward 9-sets are. It also redid the full `n = 8` search.
- **Verdicts.**
  - A1, A4, A6, A7, A8: SOUND.
  - A2 (Theorem A): SOUND WITH CORRECTION. It is Grünbaum 1971.
  - A3 (Theorem B): SOUND WITH CORRECTION. There are 5 minimum certificates, the lazy-caterer range is
    `3 ≤ n ≤ 6`, and the HYP-3805 heuristic is a tendency, not a law.
  - A5 (Theorem C): SOUND WITH CORRECTION. It is superseded by Linial–Saks–Sós 1983.
- **Agreement.** Every number matches this note's script output line for line.
