# Independent adversarial audit of opus-2026-10-26-S62 (THM-4809 / HYP-9419 / the three-quantities note)

Auditor: independent (this audit modified no repository file; all work in `scratchpad/s62_audit/`).
Audited state: worktree `/tmp/math-wt-s62` at HEAD `64a7985e70` (2026-10-26 18:57:14 -0600, "credit session72 (31002464ae) for the incidence-cycle law and the degree bound"), which is the rebase of the note commit `682ffc3817` / `054847eddb` (18:51) onto the clock-even commits. The worktree moved from `054847eddb` to `64a7985e70` during the audit; everything below refers to HEAD unless stated. One untracked file (`04-computation/experiments/restored_triples_rank_n7_20261026.cpp`, 18:58) is work in progress outside the audited commit.

## Verdict: PASS WITH CORRECTIONS

Every exact number in the note reproduces from independent code (paired cells at orders 4-7, guard ranks over Q / five primes / GF(2), pair-product spans, the complete cover-depth census, the clock census to modulus 120 for all nine sequences). The three theorems are correct. The session's self-correction at `64a7985e70` (Theorems A and B are session72's, published at `31002464ae` before the note) is adequate and is confirmed against session72's own report. The corrections required are: two of the four HYP-9419 candidates are refuted by one-minute computations that the hypothesis text itself prescribes; the "concurrent branch not on main" label is false at HEAD (the packet landed on main at 17:53, before the note, and the cited hash does not exist in the repository); a wrong characterization of the squares' composite oriented moduli; and several attribution and provenance nits.

## Majors

**M1. HYP-9419 (4) is false as stated.** "No sequence given by a polynomial of degree >= 3, nor the Fibonacci numbers, nor the primes, has a tournament pattern at any modulus." Fourth powers `k^4` give the Paley tournament at every prime `p = 3 (mod 4)` (pattern `{1, 2, 4}` at `p = 7`; tournament moduli to 120 identical to the squares' list), `k^6` at `3, 11, 23, 47, 59, 71, 83, 107` (the `3^k` list), `k^8` as `k^4`; even cubics do it: `k^3 + 2` has pattern `{1, 2, 3}` modulo 7 (a tournament), `k^3 + k^2` has `{1, 2}` modulo 5 and `{1}` modulo 3. The Fibonacci / primes / cubes clauses are true but trivial (see Minor 4). Required: retype (4) REFUTED (or restrict it to odd monomials, where it is the trivial symmetry `(-k)^e = -k^e`), and remove "polynomial of degree >= 3" from the HYP title and the ledger.

**M2. HYP-9419 (2), the support-generated rank formula, fails its own first test at every order.** It asserts that the rank of `k`-fold restored guards equals the dimension of the span of the monomials `x_S` over the sets `S` that are unions of at most `k` path supports, "predicting 8, 44, 273 from path supports alone", and prescribes "compute the support-generated monomial span at orders 5, 6 and compare". Distinct monomials are linearly independent, so the prediction is a count. Recomputed: distinct path supports `6, 31, 193, 1390` at `n = 4..7`; distinct subsets of supports (the Walsh characters available, session72's "46 linear-forest subsets" reproduced) `8, 46, 321, 2607`; distinct unions of at most two consistent supports `6, 53, 857, 25766` (all pairs: `6, 53, 889, 26379`). None equals `8, 44, 273` or `8, 64, 1024` beyond the single coincidence `8` at order 4. Session72's report (cited at HEAD) already records `46 != 44` at order 5. Required: retype (2) REFUTED as written; the sharp-depth conjecture (session72's, CITED) is unaffected.

**M3. The clock-even packet is on main, and the cited commit does not exist.** The note (section 0) and THM-4809 (status: "not on main at filing") cite "a branch not yet on main, `clock_even_guarded_composition_20261026` on the worktree `math-wt-normseven-descent-20261026` (commit edd04dad58 ...)". At HEAD the directory `05-knowledge/results/clock_even_guarded_composition_20261026/` is in the tree; its commits `257d004d6c` (17:24), `36b6b63849` (17:53), `d62d2ce890` (18:48), `92e7f12f3b` (18:55) all precede or accompany the note (18:51 / 18:57) and are ancestors of HEAD; `git cat-file -t edd04dad58` fails. The content attributed (C1 arrow trichotomy with Fibonacci modulo 5 as hostile; C2 maps `D_(a,b,q): Z/(qm) -> Z/m`; FOURIER.md F2 "a path legality guard becomes Fourier convolution"; T1 two-stack simulation) is correctly described and present in the main copy. Required: cite the main path and a real commit; delete "not on main".

## Minors

1. `04-computation/experiments/guard_rank_restored_pairs_20261026.out` is 0 bytes in both S62 commits (`31c118ba7e`, HEAD), yet the note ("also computed exactly over Q with fractions") and THM-4809's FINITE-EXACT line cite it for the exact rational ranks at `n <= 5`. The values are right (recomputed exactly: 8; 44, pairs alone 63, singles + pairs 64; and 273 at `n = 6` exactly over Z via fraction-free elimination), but the repository holds no output. Commit the output or drop the "exact over Q" claim to "mod two primes".
2. The "distinct guards" counts `10 / 47 / 280 / 1957` include the identically-zero guard (permutations using a reversed marked arc); the nonzero counts are `9 / 46 / 279 / 1956`. Ranks are unaffected; say so.
3. Section 7, squares at composite moduli "(those whose unit group makes `-1` a non-square)": wrong as a characterization. `-1` being a non-square unit is necessary but not sufficient; 53 moduli below 121 (6, 8, 14, 15, 18, 20, 22, 24, 30, 32, 35, ...) satisfy it and are not oriented, because non-unit squares collide (`3 = -3 (mod 6)`, `4 = -4 (mod 8)`, `9 = -9 (mod 18)`). The exact criterion is Theorem C's own: no nonzero squares with `x^2 + y^2 = 0 (mod m)`.
4. The "never oriented" verdicts are trivial and should be typed as such: cubes contain `1` and `-1 = (m-1)^3`; the Fibonacci pattern contains `F_(-2) = -1` by pure periodicity (checked at every `m <= 120`); the primes attain every unit class (Dirichlet; verified by the first 3000 terms for every `m <= 120`), hence `1` and `-1`; the odd Collatz successors `(3k+1)/2^v` are exactly the odd non-multiples of 3, a negation-symmetric set (verified below 4000). Section 8's reading of the Collatz fact as "division consumes precision" is an over-interpretation of a symmetry triviality. The triangular numbers are also oriented at 15 (pattern `{1, 3, 6, 10, 13}`), omitted from the note's oriented lists.
5. Theorem C's "first sentence is THM-4532's criterion": THM-4532 defines clocks (loops, doubled arcs and missing pairs allowed) but nowhere states "tournament iff `C` and `-C` partition `G \ {0}`" (no greppable form). The criterion is the definition of a circulant tournament (folklore). The `2^k` case of Theorem C is THM-4644 Part 3 verbatim ("if `<2>` has index 2 and odd order `(p-1)/2` it is the subgroup of squares; its two cosets are `QR_p` and `-QR_p`"), which the note does cite.
6. THM-1430 is a duplicated identifier: `THM-1430-explicit-symmetric-keller-counterexample-on-C6.md` and `THM-1430-graph-switching-is-exactly-E-n.md`. The Mallows-Sloane citation is the second file; cite by filename.
7. The twin collision `s(82) = s(43^2) = 44` is the identity `s(2p) = p + 3 = s((p+2)^2)` at the twin pair `p = 41`, `p + 2 = 43`, not "at `p = 43`".
8. Section 0 calls session71's Penrose guards "P1-P5 isomorphic as a named-test system to session70's conjunction family"; the Penrose report's isomorphism is P1-P3 (P4 is a development criterion, P5 the field obstruction), as the note's own section 4 table says.
9. Pre-registration: honest as to outcomes, but `restore_depth_20261026.cpp` was changed after `31c118ba7e` (class 0 "uncoverable" added; the reported statistic changed from the maximum depth to the least depth among covered states). The note should say the "least cover depth = bound" refinement is post hoc and weak: it asserts only that one state exists whose `M` free pairs are covered by `ceil(M/(n-1))` of its paths. The maximum finite cover depth is unknown (`>= 5`) from order 6 on.
10. Index hygiene: the S62 row of `05-knowledge/results/INDEX.md` is the third row, not the top row, and the S62 commit re-inserted a stale duplicate of the mac-mini Berggren row (three copies now, two with older wording). `05-knowledge/hypotheses/INDEX.md` has exactly 120 lines before and after (correct). All nine markdown links in the note resolve.
11. "The minimal complex one degree 8 with Galois group `D_4`": session71's norm audit proves `[K:Q] >= 8`, attained by `Q(i, 3^(1/4))`, and that every degree-8 Galois such field has group `D_4`; the note's wording is acceptable but slightly stronger than the source.
12. The untracked `restored_triples_rank_n7_20261026.cpp` (18:58) tests session72's order-7 case; it is outside the audited commit and must not be reported as done at close-out unless its output is committed.

## Recomputation table

Independent programs: `audit_recompute.cpp` (canonical forms by brute force over `n!` relabelings; guards as subcubes; cover depth by superset-closure bitsets), `audit_exact_rank.py` (+ sympy DomainMatrix over ZZ for order 6), `audit_clocks_fast.py` (modular arithmetic per modulus), `audit_fourth_powers.out`.

| quantity | claimed | recomputed | match |
|---|---|---|---|
| paired cells n=4 (E, src, tgt, V, c, rank, cycle) | 7, 4, 4, 8, 1, 7, 0 | 7, 4, 4, 8, 1, 7, 0 | yes |
| n=5 | 50, 12, 12, 24, 1, 23, 27 | 50, 12, 12, 24, 1, 23, 27 | yes |
| n=6 | 732, 56, 56, 112, 1, 111, 621 | 732, 56, 56, 112, 1, 111, 621 | yes |
| n=7 | 24654, 456, 456, 912, 1, 911, 23743 | 24654, 456, 456, 912, 1, 911, 23743 (E symmetric; 30 diagonal cells) | yes |
| distinct guards / permutations | 10/24, 47/120, 280/720, 1957/5040 | 9+zero/24, 46+zero/120, 279+zero/720, 1956+zero/5040 | yes (zero guard counted) |
| rational rank of single guards | 8, 44, 273 (exact Q at n<=5; two primes at n=6) | 8, 44 exact over Q; 273 exact over Z (DomainMatrix) and mod 1e9+7, 998244353, 3, 5, 7 | yes |
| GF(2) rank of single guards | 8, 43, 272, 1934 | 8, 43, 272, 1934 | yes |
| singles + pair products span | full 8, 64, 1024 at n=4,5,6 | GF(2) rank 8, 64, 1024 (full; certifies Q); mod 1e9+7 8, 64 at n<=5; exact Q 64 at n=5 (pairs alone 63) | yes |
| pairs cannot span at n=7 | <= 32647 = 2^15 - 105 - 15 - 1 | degree-12 dimension 32647 | yes |
| degree-bound dims for singles | 57 (n=5), 638 (n=6) | 57, 638 (also 8 at n=4, 9949 at n=7) | yes |
| dark-pair (uncoverable) states | 1, 13, 68, 492 | 1, 13, 68, 492 | yes |
| cover-depth histograms | n=4: 2/5; n=5: 21/28/2; n=6: 32/593/264/67; n=7: 6806/18690/6780 | identical | yes |
| least cover depth among coverable | 1, 2, 2, 3 | 1, 2, 2, 3 | yes |
| n=8 sample (10,000) | 21 dark; 109 / 3934 / 5936; least 3 | own sample: 27 dark; 82 / 3941 / 5950; least 3 | consistent (different RNG) |
| "only 21 of 64 order-5 states are two-path cubes, pairs span 64" | consistent | depth-2 count 21; pair span 64 | yes |
| squares: tournament moduli | primes 3 mod 4 below 120 (15 of them) | identical list; = primes 3 mod 4 | yes |
| squares: oriented moduli | note's 25 composites + the 15 primes | identical 40-element list | yes (parenthetical wrong, Minor 3) |
| 2^k tournaments | 7, 23, 47, 71, 79, 103 | identical; = primes with ord_p(2) = (p-1)/2 odd | yes |
| 3^k tournaments | 3, 11, 23, 47, 59, 71, 83, 107 | identical (m = 3: pattern {1}) | yes |
| triangular | tournaments 3, 5 | tournaments 3, 5; oriented also 15 | yes (15 omitted) |
| factorials oriented | 6, 20, 24, 40, 45, 60, 90, 112, 120; never tournaments | identical | yes |
| cubes, Fibonacci, primes, Collatz successors | never oriented | never oriented (trivially, Minor 4) | yes |
| dependence on "first 3000 terms" | (not discussed) | no pattern changes between 3000 and 100,000 terms for any sequence, any m <= 120 | stable |
| HYP-9419 (2) support counts | predicts 8, 44, 273 | supports 6, 31, 193; subsets 8, 46, 321; unions<=2 6, 53, 857 | REFUTED |
| HYP-9419 (4) | no degree>=3 polynomial tournament pattern | k^4 Paley at p = 3 mod 4; k^3+2 mod 7; k^3+k^2 mod 3, 5 | REFUTED |

## Attribution table

| claim in the note | source checked | verdict |
|---|---|---|
| (2+i)^3 = 2+11i, (2+11i)^2 = -117+44i, 117^2+44^2 = 125^2 = 5^6; "!5 = 44"; 44 rooted complement targets from the six strong five-tournaments; "no incidence-preserving bijection" | unnamed-carrier CONNECTIONS.md sec. 5 | present verbatim |
| three-layer proposal incidence / transport / realization | CONNECTIONS.md sec. 6 | present |
| Fermat/Dudeney norm obstruction, Q(sqrt 3) conjugate argument, Q(3^(1/4)) | CONNECTIONS.md sec. 4 | present |
| U2 50 free-pair cells at order 5 | census.json (`free_pair_cells` 2, 7, 50 at orders 3-5; `rooted_pair_cells` 64) | present |
| 44 square / 5 triangle cut types; 220 pairs, 64 after cc-diff, 46 after tri-diff | session71-dudeney/REPORT.md (A:14, B:6, C:26) | present |
| angle charge 2a/n - 2b/m in Z; charge group Z/6 (Z/10 for pentagon-square) | session71-dudeney/REPORT.md | present |
| Q(sqrt 3) cannot realize the square; Q(3^(1/4)); degree 8 complex, Galois D_4 | session71-norms/REPORT.md; session71-norm-audit/REPORT.md ([K:Q] >= 8, splitting field of X^4-3, dihedral of order 8) | present (Minor 11) |
| scaling: area 5 and 44 both need four pieces | session71 README ("pairs scaled to either area have the same four-piece optimum") | present |
| (132/5, 176/5, 44) rational right triangle | session71-norm-audit/REPORT.md | present |
| Penrose P1-P5 isomorphism | session71-penrose/REPORT.md: P1 "named combinatorial isomorphism", P2-P3 protocol/realization, P4 development criterion, P5 field obstruction | overstated (Minor 8) |
| session72: Royle-Praeger-Glasby-Freedman-Devillers, automorphism-sign-even, order-4 counts 4/4/3 | session72 README; EVEN_GRAPH_OPERATION_ENVELOPE.md ("at n=4 ... both type counts are four ... degree-even types number three") | present |
| session72 "announced" the incidence-cycle reading | session72-response/REPORT.md (commit 31002464ae): Theorem B rank 2|Q| - c(B), kernel the cycle space, table through n = 6 (732/111/621), degree bound (8), conjecture d_n = ceil((n-2)/2), rank 44, unimodular 64-certificate, order-6 determinant certificate | published, not announced; the session's credit commit 64a7985e70 now states this correctly |
| THM-4532 clocks; "tournament iff C and -C partition" | THM-4532 defines Clk(G, C); criterion absent | partly (Minor 5) |
| THM-474: labelled switching classes and labelled even graphs 2^binom(n-1,2) | THM-474 line 54 | present |
| THM-1430 Mallows-Sloane | duplicated id; graph-switching file | present (Minor 6) |
| THM-4524 E1 (A049313 = sign-even Euler graphs); THM-4531 (tournaments mod U = even graphs in U^perp, no natural bijection) | titles | present |
| THM-1600 derangements 9, 44, 265 | title (!m = 0, 1, 2, 9, 44, 265) | present |
| THM-4644 (index-2 <2>, Paley) and THM-126 | THM-4644 Part 3; THM-126 is Z_7 only | present; THM-126 loosely |
| THM-644 "2/10/44 at order 6" | THM-644 node types 1/1/2, 3/5/4, 2/10/44, 4/84/368 | present |
| HYP-9416 joint targets 1, 5, 44, 547, 9985 | HYP-9416 | present |
| THM-4770 T_d = 24 2^d - 3d - 16, capture family runs forever; THM-4771 Antihydra counter, A4/A11 Hydra | titles / bodies | present |
| THM-3333 Gaussian-square map; THM-3334 sec. 5 torsor 2^(r-1) | THM-3334 sec. 5 (Tripathi 2^(omega(c)-1), factor-choice torsor) | present |
| THM-4808 branch character | title | present |
| THM-4805 B twin collision at p = 43 | identity at p = 41 | index slip (Minor 7) |
| concurrent branch "not yet on main", commit edd04dad58 | on main since 36b6b63849 (17:53); hash absent | wrong (M3) |
| trivial labelled remark (tournaments, graphs on n, even graphs on n+1: 2^binom(n,2)) | elementary; THM-474 states the 2^binom(n-1,2) switching-class count | correct and distinct from THM-474 |
| paper: sigma = 3^(1/4), tau = 2; diameter sqrt(2 sqrt 3) < 2 kills two pieces; contradictions 3 sigma = 3 tau or a sum equal to both; no Fermat/Gaussian/lattice | pinned v3 copy: "set sigma = sqrt(sqrt 3) and tau = 2 so that both areas are sqrt 3"; Lemma 14 (no piece holds two corners of T, so at least 3 pieces); Case B: "3 sigma = 3 tau, which contradicts sigma != tau" and "||z1|| + ||x2|| = sigma and = tau"; zero hits for Fermat, Gaussian, lattice, irrational | correct |
| Waring ladder (Fermat 2 / Legendre outside 4^a(8b+7) / Lagrange 4), typed ANALOGY | classical statements, correctly stated | correct |

## Prior-art findings

- Theorem A is the rank of the incidence matrix of a (bipartite) graph: `|V| - c` with kernel the cycle space over any field (Biggs, Godsil-Royle); equivalently the dimension `|E| - |V| + c` of joint laws with fixed marginals on a bipartite support (transportation polytopes). Classical; PROVED without novelty is the right type. The repository statement is session72's (31002464ae), now credited.
- Theorem B is standard Boolean-function folklore: a subcube indicator of codimension `d` is a multilinear polynomial of degree exactly `d`; products of `k` such indicators are subcube indicators of codimension `<= kd`; the degree-`<= D` space has dimension `sum_(j <= D) binom(M, j)` (O'Donnell, Analysis of Boolean Functions, ch. 1). Session72 states it first for this carrier.
- Theorem C: the criterion is the definition of a circulant (Cayley) tournament (Alspach; Reid-Beineke); "squares give Paley iff `p = 3 (mod 4)`" is the definition of the Paley tournament; the `a^k` clause is the observation that the pattern is the subgroup `<a>`, a skew half-set iff it is the index-2 subgroup with `-1` outside; for `a = 2` it is THM-4644 Part 3. The "pattern-as-connection-set" framing is the only new element, and it is a reading, not a theorem. The census items are elementary (Minor 4).
- HYP-9419 (1) is a research frame (PROPOSED, acceptable); (3) is a question; (2) and (4) are refuted (M1, M2).

## Commands run

```
git log / git show --stat 31c118ba7e 054847eddb 64a7985e70; git ls-tree 31c118ba7e / HEAD; git cat-file -s/-t; git diff 31c118ba7e 054847eddb -- restore_depth_20261026.cpp guard_ranks_paired_cells_20261026.cpp; git merge-base --is-ancestor 92e7f12f3b HEAD
grep over CONNECTIONS.md, census.json, session71 (dudeney/norms/norm-audit/penrose reports, README), session72 README + session72-response/REPORT.md, THM-4532/474/1430(x2)/4524/4531/1600/4644/126/4770/4771/3333/3334/4808/644/4805, HYP-9416, the clock-even packet on main and on the normseven worktree
g++ -O2 -std=c++17 -o audit_recompute audit_recompute.cpp && ./audit_recompute 7 10000      # 16 s for orders 4-7 + the order-8 sample  -> audit_recompute.out
python3 -I audit_exact_rank.py 6                                                            # exact Q ranks at 4-5 (+ sympy DomainMatrix at 6 -> audit_exact_rank_n6.out)
python3 -I  (audit_clocks_fast)                                                             # modular census, 3000 vs 100000 terms -> audit_clocks_fast.out
python3 -I  (fourth powers / cubics)                                                        # -> audit_fourth_powers.out
python3 -I paper_grep.py <pinned 2412.03865v3 html> paper_grep.out; link resolver over the note (9/9 resolve); wc -l hypotheses/INDEX.md (120)
```


## Second check: the order-7 certificate (post-audit; independent verifier with its own code)

Claim checked: products of at most three named guards span all functions on `X_7` (Proposition B' of the note; certificate `04-computation/experiments/restored_products_exact_n7_k3_certificate.txt`). Verdict: **PASS**. The verifier wrote its own C++ program (states built from the 15 bits and the 6 marked arcs; all 5,040 guards evaluated per state; each certificate line's product recomputed by pointwise AND and compared with the listed subcube; GF(2) elimination on the recomputed vectors) and a Python structural cross-check (guards as constraint sets). Exact output (natural order):

```
free pairs (t: i,j): 0:(0,2) 1:(0,3) 2:(0,4) 3:(0,5) 4:(0,6) 5:(1,3) 6:(1,4) 7:(1,5) 8:(1,6) 9:(2,4) 10:(2,5) 11:(2,6) 12:(3,5) 13:(3,6) 14:(4,6)
perm[0] = 0123456, perm[5039] = 6543210
guards identically zero on X_7: 2921 (nonzero: 2119)
distinct nonzero guard functions: 1956
nonzero guard weights 2^k, counts by k: k=9:646 k=10:790 k=11:467 k=12:168 k=13:41 k=14:6 k=15:1
certificate: 2 comment lines, 32768 data lines, 0 malformed lines
lines using 0/1/2/3 guards: 0 / 1 / 6 / 32761
mask popcount histogram (codimension of subcube): 0:1 3:2 5:5 6:9 7:22 8:48 9:191 10:539 11:1245 12:2605 13:7333 14:13962 15:6806
recomputed products that are identically zero: 0
product-vs-subcube mismatches (bit positions, summed over all lines): 0
distinct subcubes described: 32768 of 32768
GF(2) rank of the 32768 recomputed product vectors: 32768 (row XOR operations: 97236)
timings: guards 0.29s, certificate+products 0.94s, elimination 5.43s, total 6.66s
checks: mismatches==0 OK; rank==32768 OK; zeroGuards==2921 OK; distinctNonzero==1956 OK
VERDICT: PASS
```

A second elimination with the rows shuffled (seed 20261026) and the columns processed in reverse order gave rank 32,768 again (49,913 row operations). Negative controls: duplicating one row gives rank 32,767 and FAIL; replacing one guard index of line 0 by guard 0 gives 31 mismatches and FAIL. The verifier agrees with both inferences: full rank over GF(2) of a 0/1 matrix gives an odd, hence nonzero, integer determinant, so rational rank 32,768; and a product of two guards fixes at most 12 of the 15 free pairs, so it lies in the degree-at-most-12 subspace (dimension 32,647) and two guards cannot span, whence `d_7 = 3`.
