# Sumner at `n = 5` on all 6880 eight-tournaments; the three 5-vertex trees against Petersen, `K_5`, `K_{3,3}`; and two critical-path Collatz tasks (the descent set equals the rising-cone union below `2^20`, reduced to a last-dip lemma; HYP-9122 typed as a digit statement)

**Session:** mac-mini (claude) 2026-10-05, worktree `codex/session-trees-tournaments-20261005`, second wave
(first wave: [trees5_tournaments4_erdos_sos_20261005](trees5_tournaments4_erdos_sos_20261005.md)).
**Owner's seed:** run the Sumner `t = 5` test on all 8-tournaments; relate the three unlabeled trees on 5
vertices to the repo's work tangential to the Petersen graph, `K_5`, `K_{3,3}`; identify and pursue the most
critical paths toward a Collatz proof.

**Status.** FINITE-EXACT (nauty `gentourng`, exhaustive): Sumner's conjecture at `n = 5` and the full
unavoidability table of the 27 oriented 5-trees (§1). PROVED (elementary): Petersen/Kneser dictionary (§2);
the Collatz-digraph host census (§1.4); the last-dip reduction and its small-`l` cases (§4). FINITE-EXACT:
`D = ⋃ rising cones` below `2^20` in the last-dip form (§4). CITED: Dross–Havet 2018, Havet–Thomassé 2000,
Grünbaum 1971, Erdős–Ko–Rado. Typed verdicts in §5. Collatz OPEN; no Collatz status changes. Independent
audit OWED on §4's reduction. Scripts in `04-computation/experiments/`:
`sumner_t5_20261005.py`, `trees5_petersen_kneser_20261005.py`, `collatz_last_dip_segments_20261005.py`,
`collatz_descent_set_vs_rising_cones_20261005.py`; outputs `*.out` in this directory (hashes §7).

## 0. Inheritance and concept board

- Sumner's conjecture (1971): every tournament on `2n-2` vertices contains every oriented tree on `n`
  vertices. `PROBLEM-LEDGER.md:379` marks it GREENFIELD, zero work. Nearest canon: THM-4526/4533 (oriented
  graphs unavoidable in every tournament; Rédei parity of embedding counts), THM-4524 (selfie tournaments),
  the Havet–Thomassé exception list reproduced in `tournament_designations_20261001.md:53-54`.
- Petersen in the repo: THM-4529 builds Petersen as `Δ–Y` on the four Pasch lines of the Paley tournament
  `P_7 - 0` (Fano/`Z/7` coordinates, **not** `GP(5,2)` and not Kneser coordinates); THM-261/262 (Petersen =
  `K(5,2)` = `A_4` root orthogonality, partly corrected by MISTAKE-507); the Kuratowski note
  `procgen_kuratowski_20260925_tait_kempe_triples.md` (each Collatz sheet is a planar functional graph; the
  Althöfer union graph `U_N` acquires `K_{3,3}` at `N = 52`, a `K_5` minor at `N = 68`, a Petersen minor at
  `N = 104`; its line 282 types the four 4-tournament classes as the "dual pair + self-dual" shape of Tutte's
  `{F_7, F_7^*} + U_{2,4}`); THM-519 (the graph-to-tournament map sends `P_5` to a tournament with `Ω = K_4`
  and `C_5` to `Ω = K_5`); THM-071 (at `n = 5` the degree-2 Walsh support of `H` is `L(K_5)` = complement of
  Petersen; the degree-4 support is the 60 Hamiltonian paths of `K_5`). The owner's own source table
  (`COLLATZ-KURATOWSKI-2026-09-25-SOURCE.md` §6: `τ(K_5) = 5^3`, `τ(K_{3,3}) = 3^4`, `τ(Petersen) = 2^4 5^3`)
  is already typed "NUMEROLOGY, recorded as hostile".
- Collatz frontier as of today (nine 2026-10-05 notes, read in full by a sweep): every proved coverage family
  has density zero (S10 Theorem CD); S10's next target is the visit measure `π` of the minimum map (needs
  basin densities, HYP-9165); the bounded obligation is the vanishing of the actual no-descent census below
  every `X` (S7, equivalent to Collatz); the "only admissible positive target" is a source-dependent lower
  bound on `λ` across a refuel boundary for an unbounded family of leaves (S5); the E-SCC endgame has Q2
  reduced to the classes `1, 14 (mod 27)` with the `1` class routed through HYP-9122 (Theorem 6.2, `q2_endgame`);
  and the S6/S7 structural question "does the descent set `D` equal the union of rising cones?" is OPEN
  (S6:9-11, S7:244-246; equality was asserted and then withdrawn in the S7 correction).
- Board: {unavoidability number `f(S)`; regular tournaments as the degree wall; `K(5,2)`; cubic hosts;
  excursions above `n`; the positive cycle point `x_w = c_w/(2^A - 3^l)`; convergents of `log_2 3`}.

## 1. Sumner's conjecture at `n = 5` (FINITE-EXACT, exhaustive)

Hosts: all tournament classes on `N = 5, 6, 7, 8` vertices from nauty's `gentourng` (12, 56, 456, 6880;
A000568 confirmed). Guests: all 27 oriented trees on 5 vertices (A000238(5) = 27: 10 oriented paths, 12
oriented forks, 5 oriented stars), generated as the 48 orientations of the three trees modulo isomorphism.
Embedding = subdigraph (not necessarily induced), by bitmask backtracking. `f(S)` = least `N` such that every
`N`-tournament contains `S` (monotone in `N`).

**Theorem 1.1 (Sumner at `n = 5`).** Every tournament on 8 vertices contains every oriented tree on 5 vertices
(185,760 host–guest pairs at `N = 8`, no failure; 203,242 checks in all). CITED: this is also the `n = 5` case of
Dross–Havet's `(21n/8 - 47/16)`-unavoidability (arXiv:1812.05167; `21·5/8 - 47/16 = 10.19`... their bound in
`n` vertices gives order 11; the sharper bound here is the exhaustive one).

**Theorem 1.2 (the unavoidability spectrum at `n = 5`).**

| `f(S)` | count | which oriented trees |
|---|---|---|
| 5 | 14 | 8 oriented paths, 6 oriented forks (all with a `c<-x` leg or balanced legs: see `.out`) |
| 6 | 8 | the two antidirected paths `+-+-`, `-+-+`; 6 oriented forks |
| 7 | 3 | the three mixed stars `out1/in3`, `out2/in2`, `out3/in1` |
| 8 | 2 | the out-star and the in-star |

- At `N = 7` (the `2n-3` layer) **exactly** the three regular 7-tournaments miss anything, and each misses
  exactly the out-star and the in-star: the only obstruction at the Sumner threshold is the degree wall (a
  vertex of out-degree 4 is forced only from `N = 8`: `⌈(N-1)/2⌉ ≥ 4`). This is the directed twin of §5 of the
  first note: at `t = 5` the star is the one tree that sees the regular extremal host, in the undirected
  (Erdős–Sós, `3n/2`) and the directed (Sumner, regular `2n-3`-tournament) settings alike.
- At `N = 6`, 16 of 56 classes miss something; every miss is a star.
- At `N = 5` (spanning), 11 of 12 classes miss something; the transitive `TT_5` contains all 27 (every acyclic
  digraph embeds in the transitive tournament); the regular 5-tournament misses 10 (both antidirected paths,
  four forks, four stars). Among paths this reproduces Havet–Thomassé/Grünbaum: the only 5-vertex exception
  to "every orientation of a Hamiltonian path" is the regular 5-tournament, and it misses exactly the two
  antidirected paths (CITED: Havet–Thomassé 2000; Grünbaum 1971).

**1.4 The Collatz functional digraph as a host (PROVED + FINITE-EXACT to `2·10^5`).** With arcs `n -> C(n)`,
`C(n) = n/2` (even) or `3n+1` (odd): out-degree 1, in-degree 1 or 2; the in-degree-2 ("branch") vertices are
exactly the `n ≡ 4 (mod 6)`; two branch vertices are never adjacent (the preimages `2n ≡ 2 (mod 6)` and
`(n-1)/3` odd are not `≡ 4 (mod 6)`), and branch vertices at distance 2 occur. Hence an oriented tree embeds
iff it is an in-tree (all out-degrees `≤ 1`) with in-degrees `≤ 2` and no two adjacent branch vertices.
**Exactly 5 of the 27 oriented 5-trees embed**: the directed path, the two other in-rooted paths, and the
fork rooted at a short leaf or at the far end of its long leg; the fork rooted at the middle of its long leg
(two adjacent branch vertices) is the one binary in-tree excluded. Compare `TT_5`: 27/27; regular
5-tournament: 17/27; Collatz digraph: 5/27. (Subcubic as an undirected graph, the Collatz graph is
`K_{1,4}`-free like `K_{3,3}` and Petersen, §2.)

## 2. The three trees against Petersen, `K_5`, `K_{3,3}` (PROVED, elementary; FINITE-EXACT census)

**Proposition 2.1 (Kneser dictionary).** Petersen `= K(5,2)`: vertices are the edges of `K_5`, adjacent iff
disjoint, `Aut = S_5`. A 4-subset of `V(Petersen)` is a 4-edge subgraph of `K_5`; the orbits under `S_5` are
six, with induced Petersen subgraphs

| 4 edges of `K_5` | count | induced subgraph of Petersen |
|---|---|---|
| star `K_{1,4}` | 5 | `4K_1` (the 5 maximum independent sets; `α(Petersen) = 4`; Erdős–Ko–Rado for `(5,2)`) |
| fork `S(2,1,1)` | 60 | `P_3 + K_1` |
| path `P_5` | 60 | `P_4` |
| paw | 60 | `K_2 + 2K_1` |
| `C_4` | 15 | `2K_2` |
| `C_3 + K_2` | 10 | `K_{1,3}` |

So the spanning trees of `K_5` (`125 = 5^3`, Cayley; `60 + 60 + 5`) are exactly the 4-subsets of Petersen
inducing `4K_1`, `P_3 + K_1` or `P_4`, and the three unlabeled trees are the three tree-type orbits. The
same table read in `L(K_5)` (= complement of Petersen, THM-071's degree-2 Walsh support) says: stars are the
4-cliques of `L(K_5)`.

**Proposition 2.2 (cubic hosts are star-free).** `K_5` contains all three trees; `K_{3,3}` and Petersen
contain the path and the fork but not the star (maximum degree 3; for `K_{3,3}` also the bipartition
obstruction: the star's colour classes are `(1,4)`, the other two trees' are `(2,3)`). Both cubic graphs have
exactly `3n/2 = ex(n, K_{1,4})` edges, so **`K_{3,3}` (`n = 6`) and Petersen (`n = 10`) are Erdős–Sós
star-extremal witnesses at `t = 5`**; neither is path- or fork-extremal (those need disjoint `K_4`'s, first
note Prop. 4). The Collatz graph is subcubic, hence star-free as well (§1.4), but with average degree `7/3`
it is far from extremal.

**Relation to the repo's Petersen work (typed).**
- THM-4529's Petersen is `Δ–Y` on the Pasch lines of `P_7 - 0` in Fano coordinates; the Kneser coordinates
  here are a different presentation of the same graph (both EXACT). Nothing in THM-4529 concerns 5-vertex
  trees; its Hasse diagram of the Petersen family is a 7-node tree and its Pasch chain `K_6 -> P_10` a
  5-node path (co-location, no map).
- The Kuratowski note's line 282 ("dual pair + self-dual" = `{TT, S} + {±}` ≈ Tutte's `U_{2,4}, F_7, F_7^*`)
  is a **shape** statement; the tree side has no dual pair at all (every tree is reversal-fixed), which is
  exactly the first note's Proposition 1(d). Typed ANALOGY of shapes, consistent with that note's own label.
- The owner's `τ(K_5) = 5^3 = 125` is the labelled count behind Proposition 2.1; `τ(Petersen) = 2000` and
  `τ(K_{3,3}) = 81` are unrelated to 5-vertex trees (their spanning trees have 9 and 5 edges). NUMEROLOGY
  stands as recorded.
- `petersen_metagraph_s287.out` found no Petersen structure on the 10 merged 5-classes; `10 = C(5,2)` is a
  coincidence of counts (same verdict as the first note's census split).

## 3. Collatz task triage (from the frontier sweep)

| candidate task | status | why pursued / not pursued |
|---|---|---|
| A. `D = ⋃ rising cones` (S6/S7 OPEN) | **pursued, §4** | decidable by a finite sweep; reduces to a one-line lemma with a Diophantine core |
| B. HYP-9122 (`minK(s) = ⌊(s+1) log_2 3⌋` for every `s ≥ 2`) | typed, not attacked | the Q2 endgame's Theorem 6.2 assumes it; its own note calls the remaining step a carry-covering statement at the convergent clocks of `log_2 3` with `~2^{0.9L}` heuristic hits, and the synthesis (line ~1313) names the gap "ternary digits of `2^K`". A quick construction is not available; FINITE-EXACT to `s ≤ 190535` already |
| C. S10's visit measure `π` / HYP-9165 basin densities | not pursued | existence of basin densities is itself open; a census extension adds digits, not structure |
| D. S7's census obligation | not pursued | equivalent to Collatz below each `X`; already verified far beyond reach |
| E. S5's "only admissible positive target" | not pursued | needs a new source-dependent inequality; no lever from this session's objects |

## 4. Task A: the descent set and the rising cones

**Objects.** Syracuse map `T(x) = (3x+1)/2^{v(3x+1)}` on odd `x`. `D = {odd n : some odd m < n has T^k(m) = n}`.
A valuation word `w = (a_1, …, a_l)` with `A = Σ a_i` is *rising* if `2^A < 3^l`; its cycle point is
`x_w = c_w/(2^A - 3^l)`, `c_w = Σ_{i=1}^{l} 3^{l-i} 2^{a_1+…+a_{i-1}}`. Shadow theorem (S6 T4, S7 audit): for
rising `w`, `n` has a backward chain along `w` landing on a positive odd `m < n` iff `n ≡ x_w (mod 3^l)` (the
*cone* of `w`). So `⋃ cones ⊆ D`; equality was OPEN.

**Lemma 4.1 (last dip; conjectured in general, FINITE-EXACT below `2^20`).** Let `m* < n` be odd with orbit
`m* = n_0 -> n_1 -> … -> n_l = n` and `n_j > n` for `1 ≤ j ≤ l-1` (the orbit's last visit below `n` before
`n`). Then `2^A < 3^l`.

**Proposition 4.2.** Lemma 4.1 implies `D = ⋃ rising cones`. *Proof.* For `n ∈ D` take any smaller odd
ancestor; along its orbit let `m*` be the last odd value `< n` before `n`. The segment `m* -> n` is an excursion
of the Lemma's type, its word is rising, and the backward chain along it is integral (it is the orbit), so
`n ≡ x_w (mod 3^l)` by the shadow theorem. ∎

**What is proved about Lemma 4.1 (elementary).**
1. `a_1 = 1` and `(2n-1)/3 ≤ m* < n`, equality iff `l = 1`: a first step with `a_1 ≥ 2` lands below `m*`.
2. Every proper suffix word (from `n_j` to `n`, `j ≥ 1`) is non-rising: a rising word has a negative cycle point
   and increases every positive `x`, but the suffix carries `n_j > n` down to `n`. (So in a violating excursion
   the full word is non-rising, with `m*` and `n` both below `x_w`, while `n` lies above every suffix's cycle point.)
3. Identity `(2^A - 3^l) n = c_w - 3^l (n - m*)`. A violation (`2^A > 3^l`) forces `c_w > 3^l (n - m*) ≥ 2·3^l`.
4. Put `t_j = 2^{A_{j-1}}/3^{j-1}` (`t_1 = 1`). From `2^{A_{j-1}} n_{j-1} = 3^{j-1} m* + c_{w_{j-1}}` and
   `n_{j-1} > n > m*`: `t_j < 1 + (3n)^{-1} Σ_{i<j} t_i`, hence `t_j ≤ (1 + 1/(3n))^{j-1}` and
   `c_w/3^l = (1/3) Σ_j t_j ≤ n((1 + 1/(3n))^l - 1) ≤ (l/3) e^{l/(3n)}`.
5. Consequently a violating excursion of length `l` using the word `w` (with `δ_w = 2^A/3^l - 1 > 0`) must have
   **`n < l / (3 ln(1 + δ_w))`**, and from `2 ≤ (l/3) e^{l/(3n)}`: the Lemma holds for `l ≤ 3` at every `n`, and
   for `l ≤ 5` whenever `n ≥ 2^20`.

**FINITE-EXACT (sweep below `X = 2^20`, 1 s).** All 561,601 last-dip segments are rising. The largest ratio
is `2^84/3^53 = 0.997914` (twelve distinct `(m*, n)` pairs, e.g. `93119 -> 93317`, `l = 53`, `A = 84`): the
extremal excursions follow the convergent `84/53` of `log_2 3`. A second, independent sweep
(`collatz_descent_set_vs_rising_cones_20261005.py`, forward orbits of every odd `m < 2^18`, each visited
`n ∈ D` flagged by whether some smaller-ancestor word is rising) finds `|D ∩ [1, 2^18)| = 61,408` (density
`0.468509` among odd, `0.702761` among units; S7's sieve value is `0.468669` at `2^28`), every element with a
rising smaller-ancestor word, and reproduces S6/S7's primitive-cone counts `1,1,0,1,0,2,8,0,28,0,124,602,0,2498,0,12319`
and extends them: `0, 64759, 353912, 0` at depths 17–20. Hence **`D = ⋃ rising cones` below `2^20`**
(FINITE-EXACT), and the two sets' density gap (`0.4687` sieve vs `0.4673` cone series at depth 16) is the tail
of long cones, not a non-cone component.

**Where the question now sits (typed).** The Lemma is a Diophantine statement about near-balanced words: by
item 5 a counterexample at the upper convergent `(A, l) = (485, 306)` would need `n < 9.98·10^4` (excluded by
the sweep); the first convergent not excluded is `(24727, 15601)`, requiring an excursion of 15,601 odd steps
above `n` that returns to `n` exactly, with `n < 2.9·10^8`. Baker-type lower bounds on `δ_w` (`|A log 2 - l log 3|`
`≥ l^{-C}`) bound such `n` polynomially in `l` but do not bound `l`; what is missing is a reason an excursion
cannot hover within a factor `1 + O(l/n)` of `n` for `l` steps while realising a non-rising word — the same
"expense is Diophantine" mechanism as S9's weak-reset records (upper semiconvergents `8/5, 27/17, 46/29, …`).
This is S8's Proposition R ("a non-rising word descends iff `x > x_w`") turned into an obligation: a
violating excursion is one trapped below the positive cycle point `x_w` of its own word. OPEN beyond `2^20`;
not a density question, so S10's "density never forces an atom" does not apply to it.

**Consequence for the frontier.** If Lemma 4.1 is proved, the descent set has an exact arithmetic description
(the union over rising words of the classes `x_w mod 3^l`), its density is the convergent cone series
(`0.4673` at depth 16, `0.468669` sieve), and S6's "replacement argument" is unnecessary. Nothing about
universal entry follows: `D` is the set of sources with a smaller ancestor, not the set of sources with a
smaller descendant.

## 5. Verdicts

| claim / question | verdict | carrier |
|---|---|---|
| Sumner at `n = 5` | FINITE-EXACT (holds); CITED consequence of Dross–Havet | all 6880 hosts × 27 guests |
| the `2n-3` obstruction | FINITE-EXACT: exactly the 3 regular 7-tournaments, exactly the two pure stars | degree wall `⌈(N-1)/2⌉` |
| trees on 5 vertices vs Petersen | PROVED dictionary (Prop. 2.1), stars = EKR independent sets | `K(5,2)` |
| `K_{3,3}`, Petersen vs the three trees | PROVED: star-free cubic Erdős–Sós witnesses | degree 3 |
| `τ(K_5), τ(K_{3,3}), τ(Petersen)` | NUMEROLOGY (as already recorded) | none |
| 4-classes ≈ Tutte's `{U_{2,4}; F_7, F_7^*}` | ANALOGY of shapes (Kuratowski note), no map to trees | — |
| `D = ⋃ rising cones` | FINITE-EXACT below `2^20`; reduced to Lemma 4.1; Lemma PROVED for `l ≤ 3` (`l ≤ 5` with the sweep); OPEN for `l ≥ 6` (Diophantine) | last-dip excursions |
| HYP-9122 | OPEN, typed as a digit statement; not attacked | carry covering at convergent clocks |

## 6. Open items and cheapest next tests

- Lemma 4.1 for `l ≥ 6`: (a) extend the sweep to `2^26` with the same 1-second-per-`2^20` cost scaled (minutes);
  (b) look for a monotonicity argument: show the excursion's running ratio `2^{A_j}/3^j` is pinned by the
  "stay above `n`" constraint so tightly that the final word cannot be non-rising — item 2 (all suffixes
  non-rising) plus item 1 (`a_1 = 1`) are the first two steps of such an induction; (c) the Diophantine side:
  for a non-rising `w` with `n < x_w`, bound the length of an excursion that stays in `(n, (3n+1)/2·…)` using
  the 2-adic rigidity of words (a word of length `l` fixes `m* mod 2^{A}`; staying above `n` for `l` steps fixes
  a residue class whose smallest member grows like `2^{A}`, against `n < l/(3 ln(1+δ_w))`). (c) is the
  promising route: it would make the Lemma a statement that long near-balanced excursions need large sources.
- Sumner at `n = 6` (`N = 10`: 9,733,056 host classes, 91 oriented 6-trees): feasible in C with the same
  bitmask test; the degree wall predicts `f = 10` only for the two pure stars and `f = 9` for the mixed stars.
- Petersen: whether the `n = 5` merged metagraph's 10 nodes admit any natural `K(5,2)` labelling was answered
  negatively by `petersen_metagraph_s287`; Proposition 2.1 suggests testing the **4-subsets** of the 10 merged
  5-classes instead (orbits under the metagraph's automorphisms vs the six Kneser orbits) — a one-script test.

## 7. Reproduction and hashes (SHA-256, raw bytes)

```
python3 04-computation/experiments/sumner_t5_20261005.py > 05-knowledge/results/sumner_t5_20261005.out               # ~2 s, needs gentourng
python3 04-computation/experiments/trees5_petersen_kneser_20261005.py > 05-knowledge/results/trees5_petersen_kneser_20261005.out
python3 04-computation/experiments/collatz_last_dip_segments_20261005.py 1048576 > 05-knowledge/results/collatz_last_dip_segments_20261005.out
python3 04-computation/experiments/collatz_descent_set_vs_rising_cones_20261005.py 65536 12 > 05-knowledge/results/collatz_descent_set_vs_rising_cones_20261005.out   # the 2^18 / depth-20 run quoted above takes 8 min
```

| file | sha256 |
|---|---|
| `sumner_t5_20261005.py` | `dbce989107d7ca09fc86538300ad84dbb776570b231b4fcf318239e7fce0a99d` (first version; the committed version adds §E) |
| `trees5_petersen_kneser_20261005.py` | `25f776350ac7cc2f0c60039d99afce7b41df212ff645c9dd8190df364c034f73` |
| `trees5_petersen_kneser_20261005.out` | `f980a23fc34a78b839304f9b5d61f7002f3d4c2999d1b75181f02fd250209ba2` |
| `collatz_last_dip_segments_20261005.py` | `f14fd4f4980ce96153ba2cf3ba12e09373fe8e1d2df8d162f20898720d885ffe` |
| `collatz_last_dip_segments_20261005.out` | `b21a43759a860bffc9e57c62fdba493d781879a3f0e018bc963ec0d360882a5d` |
| `collatz_descent_set_vs_rising_cones_20261005.py` | `7b326f9466558e8fb89802785475f1f3a9605f2152c99fddc21414e2b65ded1b` |
| `collatz_descent_set_vs_rising_cones_20261005.out` | `b18da999fd9a317ff55bb3323bde3b2aef3368195b74c787fdf5fc873d87ccd3` |

Universe and controls: hosts from `gentourng` with counts checked against A000568 and each line checked to be a
tournament; guests checked to number 27 with the `10/12/5` split; monotonicity of `f` checked; the Collatz host
census checked against the stated in-degree/branch criterion for all 27 guests; the descent-set sweep checked
against "every `n ≡ 2 (mod 3)` is in `D` via the one-step word", "multiples of 3 are never in `D`", and the S6/S7
primitive-cone counts. Hostile controls: the regular 5- and 7-tournaments; the `l = 1` equality case of item 1
(caught by a failed strict check and repaired); the convergent words as the extremal excursions.

References (CITED). F. Dross, F. Havet, *On the unavoidability of oriented trees*, arXiv:1812.05167 (2018).
F. Havet, S. Thomassé, *Oriented Hamiltonian paths in tournaments: a proof of Rosenfeld's conjecture*,
J. Combin. Theory B 78 (2000). B. Grünbaum, *Antidirected Hamiltonian paths in tournaments*, J. Combin. Theory
B 11 (1971). D. Sumner's conjecture (1971) as stated in Kühn–Mycroft–Osthus, *A proof of Sumner's universal
tournament conjecture for large tournaments*, Proc. LMS 102 (2011). Erdős–Ko–Rado (1961). OEIS A000568, A000238.
