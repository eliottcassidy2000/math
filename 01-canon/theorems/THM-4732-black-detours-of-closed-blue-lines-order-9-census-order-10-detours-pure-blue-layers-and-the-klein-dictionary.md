---
id: THM-4732
title: "Black detours of closed blue lines and the Klein dictionary. (1) The class-level line multigraph at n = 9 by a colour-refinement canoniser (2^28 tilings, 314 s): 191,536 classes, every control of THM-643 / THM-4729 / THM-4731 reproduced (SC = 2752, pure-blue 6, black self-lines 1704, merged black self-loops 8998, depth counts = THM-4731 (5)); the line metagraph and the blue spine on the 2752 self-converse classes are connected and the black graph has exactly 1 + pure-blue = 7 components, so all 191,530 non-pure-blue classes are black-connected (HYP-9345 (b) / HYP-9346 (2) at n = 9) and every closed blue line between non-pure-blue classes is black-redundant in every endpoint layer. (2) Lemmas: a detour of L black K moves exists only if E(Q_1) xor ... xor E(Q_L) = E(P) (L odd) or its complement (L even), so one move is impossible and an even number L needs n <= 2L + 2 and C(n-1, 2) even (every black distance is odd for n = 0, 3 mod 4 and at least 3 for n >= 7); black-core connectivity holds iff the pure-blue-deleted line metagraph is connected and every blue line between non-pure-blue classes is black-redundant. (3) Literal black detours (exhaustive over all 2^(e-1) closed blue presentations, n = 5..10): at n = 6, 8, 9, 10 every closed blue line between non-pure-blue classes has three successive black K moves whose literal labelled target is its line partner, in every layer (524,288 lines at n = 10, 321 s); this independently confirms the two concurrent order-10 certificates of the endpoint-three layer (n-4, 3) (every n >= 7 with the lane's n >= 11 theorem; the order-six C3 -> C3 line is pure-blue-touching) and settles every layer, in particular (n-5, 4), through n = 10, and the open sub-case of the (n-5, 4) layer at n = 11 (closed middle with y, c in I: all 196,608 presentations have three-move detours); the only exceptions at n <= 10 are the lines F_5 and F_7 whose open member is the transitive tournament with the arcs d -> 0 and n-1 -> d reversed, d = (n-1)/2 (layer n-3; the open member is path-poor: 6 and 20 black Hamiltonian paths), of black distance exactly 5 by bidirectional search, while F_9 and F_11 have distance 3. (4) For n >= 6 a closed blue line touching a pure-blue class has its pure-blue end open and lies in layer n-2 (bottom block K1) or n-3 (bottom block C3) of the palindromic stack, with sum tc over those stacks as the count (8 at n = 9, 4 at n = 10), so the layers k <= n-4 are pure-blue-free. (5) The Klein dictionary of the trijection reset's 3-4-3 / 4-2-4 chart: kappa (complementation) and kappa sigma (path reversal) are fixed-point-free on tilings at every n >= 3 (the apex tile; THM-849 section 2, cited), so the tiling Klein-set <kappa, sigma> has only free orbits (black twin pairs) and <sigma>-orbits (blue lines); the proper divisors of a non-square N with two equal-exponent primes q, r (Q = their exchange, C = complement) embed <C,Q>-equivariantly into tilings iff some other prime has odd exponent (p^k q r: iff k is odd), the two central CQ-fixed points of the 4-2-4 split being exactly the obstruction (for p^k q r this is, concurrently and earlier by hours, Theorem 1 of divisor_reflection_resolution_20261021); the eight non-central divisors of p^2 q r form a Klein-set isomorphic to the eight tilings of order 4 (C = complementation, Q = converse with path reversal, CQ = path reversal), with exactly 32 equivariant bijections, 4 reading sign f as closedness, 16 reading sign g as chirality, 2 reading both, 2 divisibility-monotone (one of which reads sign g), and none both sign-f-reading and monotone"
status: >
  PROVED: (2) both lemmas; (4) with the pure-blue alphabet and centre law as hypotheses (THM-2444/2454 FINITE-EXACT to
  n = 13, THM-4730 for all n with the audit caveat recorded in THM-4729's status) and THM-4729 L4 (pure-blue nodes are non-strong for n >= 6; conditional on the repaired odd-order part of
  THM-4730) with L3 (non-strong nodes are fully open); (5) the orbit-type classification of both Klein-sets (the tiling side, i.e. the fixed-point-freeness of kappa
  and kappa sigma and the two orbit types, is THM-849 section 2, CITED), the embedding criterion and its specialisation to p^k q r, the isomorphism of the eight non-central
  divisors with the order-4 tilings; the non-existence of an order isomorphism (three minimal divisors, one minimal
  tiling).
  FINITE-EXACT: (1) the n = 9 census (tiling_line_metagraph_n9_20261020.cpp; n = 5..8 agree with S46's independent
  brute-force canoniser class by class in the counts; fibre law tc |Aut| = H with H by subset DP for all 191,536
  classes; the blue-parity, degree and 2^m mass laws); (3) the detour census n = 5..10
  (blue_line_black_detours_20261020.cpp, exhaustive over Q_1, Q_2 for lines without a certificate; literal target
  verified arc by arc) and the black distances of F_5, F_7, F_9, F_11 (blue_line_deep_detours_20261020.cpp,
  bidirectional BFS); the identification of F_5's open member with the tribonacci node T*_5 and of F_7's with a
  tournament of T*_7's score multiset not isomorphic to T*_7 (VF2); (5) the counts 32 / 4 / 16 / 2 / 2 over all 8!
  bijections (trijection_klein_dictionary_20261020.py).
  CITED: THM-643 (frame), THM-644, THM-646, THM-2444, THM-2454 (alphabet and centre law), THM-4566 (idoneal planes),
  THM-4729 (closed/open, L3, Proposition K), THM-4730 (strong alphabet), THM-4731 (twins, depth, tribonacci node),
  HYP-9345, HYP-9346; the endpoint-two and endpoint-three packets (05-knowledge/results/endpoint_two_repairs_20261020g,
  endpoint_three_repairs_20261020h: the (n-4,3) theorem for n >= 11, the "no order-ten census" scope); the trijection
  reset packet (05-knowledge/results/trijection_reset_20261020: the chart, C/Q/CQ, the tetrahedral obstruction); the
  105 guard packets (guarded_children_105_20261018g, native_cancellation_105_20261018g: the 153N identity and the open
  increment bound, untouched here).
  CONCURRENT (same day, read after this session's pre-audit checkpoint): endpoint_three_order10_20261021 and
  endpoint_four_20261020h both closed the order-ten endpoint-three case (affine cells; deterministic bank over the same
  86,016 sources) before this census ran; endpoint_four_repairs_20261021 proved the (n-5,4) layer for n >= 15 except the
  closed middle with I = {y,c,h}, c -> x or (x -> c, h in {x, n-4}) (F_7 lies in that sub-case);
  divisor_reflection_resolution_20261021 (fixed-pair lemma, Theorem 1, twelve-state repair, capacity iff) and
  reflection_context_20261020h state the even-k obstruction; idoneal_information_channels_20261021 gives the exact
  idoneal information law; mac-mini-2026-10-20-debtcones published THM-4733..4736 / HYP-9348/9349 at 18:40-18:41 (pre-audit): THM-4734 has
  the same k = 3, 4 finite results at n <= 10 by the lanes' banks plus search (its "every n >= 6" excludes the order-six
  C3 -> C3 line as the lane's recorded hostile), THM-4733 the same Klein obstruction with a fold of the order-five cube;
  the Codex continuations endpoint_four_endpair_20261021 / endpoint_four_frontier_20261021 (18:43-18:44) repair the
  end-pair branch and the stored n = 15, 18 hostiles by odd-prefix cuts; messaged, cross-cited, not reproduced. This entry cites and does not modify any of them; the endpoint lanes and
  the trijection reset are also same-day packets. No literature-priority claim.
  OPEN: a structural proof of the three-move law for n >= 8 (HYP-9347 (1)) and of the (n-5,4) layer for n >= 11; the
  Collatz increment bound (section 6 of the note); the chirality refinement of the fundamental black graph.
session: opus-2026-10-20-S47
source: 05-knowledge/results/black_detours_closed_blue_lines_klein_dictionary_20261020.md
depends_on:
  - 01-canon/definitions.md
  - 01-canon/theorems/THM-643-blue-black-line-parity-structure.md
  - 01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md
  - 01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md
  - 01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md
related:
  - 01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md
  - 01-canon/theorems/THM-2444-pure-blue-count-refutation-and-rigid-class-law.md
  - 01-canon/theorems/THM-4566-hilbert-class-field-planes-are-the-idoneal-planes.md
  - 05-knowledge/hypotheses/HYP-9345-merged-metagraph-open-laws-after-thm-4729-4730.md
  - 05-knowledge/hypotheses/HYP-9346-self-lines-connectivity-and-growth-laws-of-the-tiling-line-metagraph.md
  - 05-knowledge/hypotheses/HYP-9347-three-move-law-layer-theorems-fundamental-black-graph-and-the-klein-dictionary-beyond-order-4.md
  - 05-knowledge/results/endpoint_three_repairs_20261020h/README.md
  - 05-knowledge/results/trijection_reset_20261020/arithmetic.md
scripts:
  - 04-computation/experiments/tiling_line_metagraph_n9_20261020.cpp (+ _n9.out, _classes_n5..n8.tsv, _classes_n9.tsv.gz)
  - 04-computation/experiments/tiling_line_metagraph_n9_analysis_20261020.py (+ .out)
  - 04-computation/experiments/blue_line_black_detours_20261020.cpp (+ _n5..n10.out, _n11_k4sub.out)
  - 04-computation/experiments/blue_line_deep_detours_20261020.cpp (+ _family.out)
  - 04-computation/experiments/tribonacci_exception_check_20261020.py (+ .out)
  - 04-computation/experiments/trijection_klein_dictionary_20261020.py (+ .out)
---

# THM-4732: black detours of closed blue lines, the order-9 census, the order-10 detours, pure-blue layers, and the Klein dictionary

Conventions: THM-643's tiling model (base path `0 -> ... -> n-1`, tiles = non-path pairs, `kappa` = flip all tiles, `sigma` = staircase reflection, blue = `sigma`-fixed); THM-4729's closed/open (apex tile `(0, n-1)` backward = closed); the endpoint lanes' `K_Q(T) = T^op xor E(Q)`, blue/black words, and the layer `k = d_T(z)` of the closed member (score pattern `(n-1-k, k)`). A **literal black detour** of length `L` is a chain of black words `Q_1, ..., Q_L`, `Q_i` a directed Hamiltonian path of `T_(i-1)`, `T_i = K_(Q_i)(T_(i-1))`, from `T_0 = T` to `T_L = K_P(T)` as labelled tournaments; the **black distance** is the least `L`. The proofs and tables are in the [session note](../../05-knowledge/results/black_detours_closed_blue_lines_klein_dictionary_20261020.md).

## 1. The order-9 line multigraph (FINITE-EXACT)

| n | classes | SC | merged | PB / MX / PK | blue / black lines | self-lines (blue, black) | merged self-loops (blue, black) | reversal-neutral | components (all, blue on SC, black) |
|---|---|---|---|---|---|---|---|---|---|
| 8 | 6880 | 176 | 3528 | 3 / 173 / 6704 | 2048 / 1046528 | (4, 202) | (4, 618) | 416 | 1, 1, 4 |
| 9 | 191536 | 2752 | 97144 | 6 / 2746 / 188784 | 32768 / 134184960 | (0, 1704) | (0, 8998) | 7294 | 1, 1, 7 |

Lines by node type at `n = 9`: blue `(PB,MX) 8`, `(MX,MX) 32760`; black `(MX,MX) 72524`, `(MX,PK) 4,297,276`, `(PK,PK) 129,815,160`. Lines by depth `32768, 2064384, 31457280, 100663296` (THM-4731 (5)). Handshake: 1358 self-converse nodes of odd fundamental black degree, `(T - SC)/2 = 94392`, both even. Single-presentation classes (`tc = 1`, `H = |Aut|`) are exactly the stacks over `{K1, C3}` at every `n = 5..9` (`sum_j C(n-2j, j) = 4, 6, 9, 13, 19`, `H = 3^j`); the palindromic ones are the pure-blue classes (THM-2444's rigidity lemma) and the others (14 at `n = 9`) are exactly the pendant nodes of the fundamental black graph `B` (a merged node of `B` has degree `tc` or `(tc - bm)/2`). Pure-blue classes: the six palindromic stacks `T_9, (1,1,1,C3,1,1,1), (C3,1,1,1,C3), (1,C3,1,C3,1), (1,1,T5,1,1), (C3,C3,C3)`, blue degrees `1, 1, 1, 1, 3, 1`. Closed blue lines by layer `k = 1..7`: `512, 3072, 7680, 10240, 7680, 3072, 512`; every line with two non-pure-blue classes is black-connected; the 8 pure-blue-touching lines sit at `k = 6` (2) and `k = 7` (6).

## 2. The two lemmas (PROVED)

**Lemma 1.** `K_Q` flips the pairs outside `E(Q)`, so `T_L = T xor (A \ E_1) xor ... xor (A \ E_L)` and `K_P(T) = T xor (A \ E(P))`: a length-`L` detour needs `E_1 xor ... xor E_L = E(P)` (`L` odd) or `= A \ E(P)` (`L` even). `L = 1` forces `Q_1 = P`, blue. Even `L` needs `C(n-1, 2) <= L (n - 1)`, i.e. `n <= 2L + 2`, and `C(n-1, 2)` even (the symmetric difference of `L` sets of size `n - 1` has even size for even `L`), i.e. `n = 1, 2 mod 4`; so black distances are odd for `n = 0, 3 mod 4` and at least 3 for `n >= 7`.

**Lemma 2.** All non-pure-blue classes lie in one black component iff (a) the line metagraph with the pure-blue nodes deleted is connected and (b) every blue line between two non-pure-blue classes has black-connected ends. (Necessity: the black graph is a subgraph of (a)'s graph; sufficiency: replace each blue line of a walk by (b)'s black walk.)

## 3. The detour census (FINITE-EXACT)

| n | closed blue lines | three-move detour | pure-blue-touching | non-pure-blue lines without a three-move detour |
|---|---|---|---|---|
| 5 | 8 | 2 | 5 | 1 (`F_5`, layer 2, distance 5) |
| 6 | 32 | 30 | 2 | 0 |
| 7 | 256 | 249 | 6 | 1 (`F_7`, layer 4, distance 5) |
| 8 | 2048 | 2045 | 3 | 0 |
| 9 | 32768 | 32760 | 8 | 0 |
| 10 | 524288 | 524284 | 4 | 0 |

Per layer at `n = 10` (lines / three-move detours): `4096/4096, 28672/28672, 86016/86016, 143360/143360, 143360/143360, 86016/86016, 28672/28671, 4096/4093`. No two-move detour occurs at any `n`. `F_n` is the closed blue presentation with exactly the tiles `(0, d)` and `(d, n-1)` forward, `d = (n-1)/2`; its open partner is the transitive tournament with `d -> 0` and `n-1 -> d` reversed (`T*_5` at `n = 5`); black distances `5, 5, 3, 3` at `n = 5, 7, 9, 11`.

Consequence for the endpoint lanes: the `(n-4, 3)` three-move theorem, proved there for `n >= 11`, holds at `n = 10` (86,016 lines; two concurrent packets certify the same by other methods), hence at every `n >= 7` (at `n = 6` one line touches the pure-blue `C3 -> C3` class); the `(n-5, 4)` layer has the three-move property at `n = 8, 9, 10` and fails at `n = 7` only through `F_7`, which lies in the sub-case left open by the concurrent `n >= 15` theorem; at `n = 11` that sub-case (closed middle, `y, c in I`; `196,608` presentations) has a three-move detour for every presentation (FINITE-EXACT, `blue_line_black_detours_20261020_n11_k4sub.out`).

## 4. Pure-blue-touching layers (PROVED, scoped)

For `n >= 6` the pure-blue end of a closed blue line is open (THM-4729 L4, conditional on the repaired odd-order part of THM-4730, with L3), the last vertex `z` of its presentation has out-degree `delta = 0` (bottom block `K1`) or `1` (bottom block `C3`; `T5` occurs only at the centre by THM-2454), and the closed partner has `d(z) = n - 2 - delta`. So the layer is `n - 2` or `n - 3`, the count in layer `k` is `sum tc(C)` over the pure-blue stacks with `delta = n - 2 - k`, and the layers `k <= n - 4` are pure-blue-free (`k = 3` for `n >= 7`, `k = 4` for `n >= 8`). At `n = 5` the strong pure-blue `T5` (`HC = 2`) gives two closed presentations at `k = 2` and an open one whose partner has `k = 1`.

## 5. The Klein dictionary (PROVED; counts FINITE-EXACT)

(i) (THM-849 section 2, cited) for `n >= 3`, `kappa t = t` is impossible; `kappa sigma t = t` would flip the `sigma`-fixed apex tile; `Fix(sigma)` has `2^e` elements. Orbits of `<kappa, sigma>`: `2^(e-1)` blue lines (stabiliser `<sigma>`) and `(2^m - 2^e)/4` black twin pairs (free). (ii) On `D(N)` (`N` non-square, primes `q, r` of equal exponent, `Q` their exchange, `C(d) = N/d`): `C` is fixed-point-free; `CQ d = d` iff every other exponent is halved (so iff every other exponent of `N` is even) and the `q, r` exponents are complementary. (iii) An equivariant injection `D(N) -> tilings` exists iff some exponent outside `{q, r}` is odd; for `p^k q r` iff `k` is odd; for even `k` the `CQ`-fixed central pair (`f = 0`) is the whole obstruction. (iv) For `p^2 q r`, the eight non-central divisors are the order-4 tilings: `{p, pqr}`, `{p^2, qr}` the two blue lines, `{q, r, p^2 q, p^2 r}` the black twin pair; of the 32 equivariant bijections, 4 read `sign f` as closedness, 16 read `sign g` as chirality (`[tile (0,2) forward] - [tile (1,3) forward]`), 2 read both (`p^2 -> T*_4`, `qr -> T_4`, `p -> (0,2),(1,3)` forward with apex backward, `pqr -> apex only`, `q -> (0,2)+apex`, `r -> (1,3)+apex`, `p^2 q -> (0,2)`, `p^2 r -> (1,3)`, and the same with the blue lines exchanged), 2 are divisibility-monotone (`p -> empty`, `pqr -> full`; one of them reads `sign g`, neither reads `sign f`), none is both `sign f`-reading and monotone, none is an order isomorphism. So the chart's two reflections are complementation and converse-with-path-reversal, with path reversal as their product; the full ten cannot be matched at any order because path reversal never fixes a presentation (the obstruction half is Theorem 1 of the concurrent divisor_reflection_resolution_20261021; the general-N criterion, the counts and the two dictionaries are this entry's).
