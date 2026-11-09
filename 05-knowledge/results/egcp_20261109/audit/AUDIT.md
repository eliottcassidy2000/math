# Audit: THM-5040 (sharp even-graph Chan-Pak), opus-2026-11-09-S89

VERDICT: PASS WITH CORRECTIONS

Independent adversarial audit, 2026-11-09. I read THM-5040 in full (HEAD 7e86bdf8f9 and its predecessors; the file
was being edited during the audit), THM-4996, PROOF.md Part B, HYP-9575, both packet READMEs and every script in
`04-computation/experiments/egcp_20261109/`. All four requested scripts were re-run. Every lemma was re-checked by
hand and by new code in this directory (`indep_lib.py`), which shares nothing with `egcp_lib.py`. It has its own
dominating-set counters (subset DP plus an inclusion-exclusion cross-check) and its own automorphism engine
(individualisation-refinement, with generators along a stabiliser chain), validated against all n! permutations.

**Overall.**
- **The theorem is correct.** E(n) = O(n) for every n. The exclusions (Lemmas 1-3 and {3, 7, 13}), the four
  constructions (Lemmas 5-8), the frontier families (Lemma 9), Theorem T, Lemma 10 and every case boundary of the
  Section 4 induction check out by hand.
- **Independent end-to-end confirmation.** I re-implemented the Section 4 recursion from the text. It built all
  7,934 values that need exactly n vertices for 10 <= n <= 14, and each graph was verified by brute force (order, d
  and evenness): 0 failures.
- **Independent census.** An exhaustive census of all labelled graphs on <= 7 vertices gives E(n) = O(n) for
  n <= 7.
- **Corrections:**
  - two MAJOR items: a missed priority (Brouwer-Csorba-Schrijver 2009, Proposition 1.2, already gives completeness
    *and* conciseness over all graphs), and a false claim in Corollary 4;
  - two MINOR items: Lemma 0 is known; the Corollary 3 arithmetic needs m >= 6 and m >= 4;
  - six NOTEs.
- **None of these touches the main theorem or its proof.**

## Findings

### 1. MAJOR (priority): completeness and conciseness over all graphs are Brouwer-Csorba-Schrijver 2009, Prop. 1.2

**What is known.**
- A. E. Brouwer, P. Csorba, A. Schrijver, *The number of dominating sets of a finite graph is odd*, preprint
  (June 2, 2009), https://aeb.win.tue.nl/preprints/domin4a.pdf. This is the paper already credited for parity.
- Its Proposition 1.2 states that for every odd m with 0 < m < 2^n there is a graph on n vertices with exactly m
  dominating sets.
- The proof is the induction "add an isolated vertex (m -> m) or a vertex joined to all old vertices
  (m -> 2^n + m)". That is THM-4996 Lemma B1 verbatim, including the conciseness ceil(log2(k + 1)).
- pdftotext renders the bound as "2n". I confirmed from the PDF font data (PyMuPDF) that the n is a 7 pt raised
  glyph after a 10 pt "2", so it reads 2^n, both in the statement and in "2^n + m".
- It is restated as Proposition 2.1(iii) of Alikhani-Fatehi-Klavžar, *On the structure of dominating graphs*,
  arXiv:1512.07514, which attributes it to [BCS, Prop. 1.2] ("see also" Alikhani 2013).

**Alikhani 2013 proves something else.** S. Alikhani, *The domination polynomial of a graph at -1*, Graphs Combin.
29 (2013) 1175-1181, is about D(G, -1).
- Its abstract, as restated in Alikhani-Deutsch, arXiv:1305.3734, Introduction, gives for each odd n a connected G
  with D(G, -1) = n.
- The reading "number of dominating sets = n" comes from Alikhani-Fatehi, arXiv:1410.1169. I could not access the
  full text to confirm that the 2013 paper also states the D(G, 1) version.

**Where the canon is wrong.**
- THM-5040, lines 243-246: credits all-graph completeness to Alikhani 2013 and says "the conciseness
  ceil(log2(k + 1)) appears to be new".
- THM-5040, lines 42-43 (Corollary 3): "All graphs attain the lower bound (THM-4996 (1))".
- THM-4996, lines 9-12: the new MISTAKE-747 correction line says "the conciseness ... remains new". Audit P1 of
  THM-4996 missed BCS Prop. 1.2 too.
- HYP-9575, line 35 (status of (B)).
- `05-knowledge/results/egcp_20261109/README.md`, lines 59-60: "only its conciseness part remains new".
- The concurrent THM-5055, line 35, repeats the attribution.

**Fix.**
- Credit BCS 2009, Proposition 1.2 for completeness and conciseness over all graphs (THM-4996 (1)/Lemma B1 is that
  proposition re-proved).
- Credit Alikhani 2013 only as "D(G, -1) for connected G (cited by Alikhani-Fatehi for connected realisers of d)".
- Delete "conciseness appears to be new" everywhere. What is new in THM-5040 is the even-graph statement: E(n), the
  exceptional set {3, 7, 13}, X(n), and the +1 cost.
- Log a MISTAKE: the parity paper's own Proposition 1.2 was not read. This is the pattern of the
  read-cited-papers-fully memory note.

### 2. MAJOR (false claim in a corollary): top-band values ARE realised by disjoint unions

**Location.**
- THM-5040, lines 50-51 (Corollary 4, first bullet): "The top bands T(m) are realisable only on m vertices, and only
  by graphs that are neither unions nor joins (Lemma 0)."
- Repeated in `05-knowledge/results/egcp_20261109/README.md`, lines 57-58.

**Counterexample.**
- P_3 u G_7 with G_7 the census witness for d = 101 (padded to 7 vertices) has 10 vertices.
- d = 5 * 101 = 505, which lies in T(10) = (503, 512).
- It is even (verified: `cor4_unions.out`; also by B2.2, since P_3 and G_7 are even).

**How common.** Nontrivial unions (both factors > 1) realise most top-band values: 1 of 5 at m = 10, all 8 at
m = 16, and 10 of 12 at m = 24 (`cor4_unions.out`).

**What is true.**
- No join realises a top-band value. A join A v B on m vertices has d >= 2^m - 2^a - 2^b + 3 >= 2^(m-1) + 1.
- No graph with an isolated vertex realises one. A graph on m vertices with an isolated vertex realises only values
  of E(m - 1).

**Fix.**
- Replace the bullet with: "only on m vertices, by graphs with no isolated vertex that are not joins (in particular
  not cones); many are disjoint unions, e.g. P_3 u G_7 gives 505".
- Clarify "Joins over them inherit the gap". Presumably B5 cones only over padded graphs, so it cannot cone over a
  top-band graph.

### 3. MINOR (priority): Lemma 0 is known

- **Complement formula.** I. Heinrich, P. Tittmann, *Counting dominating sets of graphs*, arXiv:1701.03453, prove
  D(G, x) + N(comp G, x) = (1 + x)^n. Here N is the neighbourhood polynomial, the face polynomial of Lovász's
  neighbourhood complex plus the empty set, so N(H, 1) = D(H) in THM-5040's notation. BCS 2009 already observe that
  the non-dominating sets form the neighbourhood complex of the complement.
- **Union additivity.** J. I. Brown, R. J. Nowakowski, *The neighbourhood polynomial of a graph*, Australas. J.
  Combin. 42 (2008) 55-68, Lemma 3.
- **C_4-free formula.** Same paper, Theorem 1.
  - Their formula, sum_v (1+x)^deg(v) - x(2m - n) - (n - 1), silently assumes no isolated vertices: it is off by
    one for each isolated vertex, for example K_1.
  - THM-5040's form 1 + #(non-isolated) + sum (2^deg - 1 - deg) is the correct general version. I checked it on
    2,243 random C_4-free graphs.
- **Location:** THM-5040 line 57 ("complement formula; elementary") and line 249 ("Elementary: Lemma 0").
- **Fix:** credit Heinrich-Tittmann (and BCS) and Brown-Nowakowski, and note the isolated-vertex correction.

### 4. MINOR (arithmetic): Corollary 3's count and example need lower bounds on m

- **Count.** Line 45: "These are 1 + floor((m - 2)/2) values per order". The number of k with
  V(k) = ceil(log2(k+1)) + 1 = m is 0, 0, 0, 1 for m = 2, 3, 4, 5, against 1, 1, 2, 2 from the formula. It is
  correct for every m >= 6, checked to m = 20 (`indep_census_7.out`).
- **Example.** Line 46: "2^m - 1 needs m + 1 vertices". This is false for m = 1 (k = 1 needs 1 vertex) and
  vacuous for m = 2, 3 (3 and 7 are unattainable). It is true for m >= 4.
- **Fix:** add "for m >= 6" and "for m >= 4".

### 5. NOTE (Lemma 8): hypothesis of (iii) is automatic; two wording gaps

- **The hypothesis always holds for rho >= 2.**
  - If k >= 2, then deg a = m - 2 >= rho + 1 exceeds every other degree, which is at most rho.
  - If k = 1, then z is universal in H, so no other vertex of H is isolated in H, and every vertex of V(H) - z has
    degree <= rho - 1 < rho = deg a.
  - Exhaustive check: all 5,405 labelled (H, z) with rho <= 5 (`exhaustive_small.out`). The only exception is
    rho = 1, which is irrelevant.
  - "(which holds when rho > h + 1)" (line 117) is correct but understated. The lemma can be stated
    unconditionally for rho >= 2.
- **"U" is undefined** (line 129). Fix: "the set U of the k leaves".
- **"Composing with an even permutation of them"** (line 128) must mean an even automorphism: any permutation of
  non-adjacent twins has edge-sign +1 (B2.1). Read as "alternating-group permutation" it can fail. For example, with
  k = 1 and one universal vertex w != z of H, the only correcting map is the transposition (u w). The argument is
  correct with the edge-sign reading. Fix: "composing with a permutation of them (an even automorphism by B2.1)".
- **(ii) has no written proof.** It is trivial: a misses z, V(H) misses the leaves, and the leaves have degree 1.
- My check of criterion (iii) in both directions passed:
  - exhaustively for rho <= 5;
  - on 2,500 random cases, of which 1,325 had universal vertices of H other than z and 877 had k = 1.

### 6. NOTE (Lemma 5', Lemma 5, Lemma 0, Section 4): implicit side conditions

- **Lemma 5' (line 76).** "a fixed order-preserving isomorphism" presupposes that isomorphic blocks are labelled
  compatibly. Say so; the sign does not depend on the labelling.
- **Lemma 5 needs g >= 1.** For g = 0 the value formula is off by one. Lemma 0's union formula needs both parts
  nonempty. Both always hold where the lemmas are used (g >= 3 in Section 4).
- **Case 2a (line 211).** The graph from Claim_3(g) has at most g vertices and must be padded to exactly g before
  Lemma 5 is applied.

### 7. NOTE (Conclusion wording, lines 225-227)

"The odd numbers in O(n) above 2^n - n - 3 are ... M(n) and the excluded X(n)": X(n) is not in O(n). Suggested
wording: "The odd numbers in (2^n - n - 3, M(n)] are M(n) (n even), and X(n), M(n) (n odd); X(n) is not in O(n)."

### 8. NOTE (cover.py): the "two-hub" family is undefined and unproved in the packet; the top band is injected

- **The family.** `cover.py` (line 63) uses a family with d = 2^(n-1) + 2^f - 2^t + d(G). Its only description is
  "two-hub graphs" (THM-5040 line 239; script docstring).
- **Identification.** I identified it as TH(G, t): two non-adjacent hubs, both joined to all of V(G), with t pendant
  leaves at one hub.
  - Hand derivation: with a: 2^t (2^(f+1) - 1); without a: 2^f + d(G).
  - Brute-force check: `twohub_probe.out`, `twohub_parity.out`.
- **Parity.** TH(G, t) is even iff G is even, except (f, t) = (1, 1). There TH = P_4, which is odd. Proof: for
  f >= 2 the leaves all hang at a and a's unique non-neighbour is b, so both hubs are fixed; for f = 1 and t >= 2,
  a has t >= 2 leaves and g has one. cover.py never uses (1, 1) (n >= 10), so its closure is sound.
- **Top band.** cover.py also adds topband(n) directly, that is, it assumes Theorem T. With `--notop` it still
  closes E(n) = O(n) for every n <= 25, through SR on the census L sets plus TH. So the closure is independent of
  Theorem T, but not of TH.
- **Fix:** document TH(G, t) with this short proof, or describe cover.py's output as corroboration that relies on
  it.

### 9. NOTE (provenance)

- When the audit started, `construct_12.out` held only the n = 10, 11 lines (the run was in progress). It is now
  complete ("TOTAL 706 BAD 0") and committed in 7e86bdf8f9, which matches the status claim "10 <= n <= 12".
- construct.py checks delicate values plus a random sample per order. This audit checked every value that needs
  exactly n vertices for 10 <= n <= 14.

### 10. NOTE (no prior art found for the even-graph statement)

Searches found nothing on dominating-set counts of even graphs (von Brömssen), the exceptional set {3, 7, 13}, the
sharp orders M(n), X(n), or the +1 cost of evenness. Searches covered: BCS 2009; Heinrich-Tittmann; Alikhani 2013 and
the citing papers 1305.3734, 1401.2092, 1410.1169 and 1512.07514; Brown-Nowakowski; RPGFD 2023; web queries on
"number of dominating sets" together with "every odd", minimum order, and even/odd automorphisms. So the core claim
of THM-5040 stands as new. The concurrent independent proof THM-5055 (mac-mini) agrees with the statement.

## Checklist (items 1-12 of the brief)

| # | Item | Verdict | Details |
| --- | --- | --- | --- |
| 1 | Lemma 0 | CORRECT (known; Finding 3) | Complement formula: S misses v iff S is inside N_comp(v). Union: the empty set is the only common set. C_4-free: a set of size >= 2 lies in at most one neighbourhood, and the formula counts singletons once. Brute force: 300 / 200 / 2,243 random cases, 0 failures. |
| 2 | Lemma 5' and Lemma 5 | CORRECT (Note 6) | Aut(comp J) = product of wreath products, generated by block automorphisms and block swaps. A block automorphism keeps every cross edge. A swap reverses s^2 edges between the swapped blocks and 2s\|A_l\| for each block between them. Lemma 5 value: D(tK_2 u comp G) = 2t + D(comp G). If G has a universal vertex, it is a single K_1 co-component, and no E_2 block is isomorphic to it. A join G simply contributes several co-components. Tests: 330 random multi-joins with repeated blocks; 400 Lemma 5 cases (131 with G a cone or a join). |
| 3 | Lemma 6 | CORRECT | D(tK_2) = 2t + 1, so d(CP_n) = 2^n - n - 1. Cones are even by B2.3 (CP has no universal vertex). |
| 4 | Lemma 7 | CORRECT | Value: 2(2^n' - 2^(n'-2)) + d(G). deg a = n' > n' - 1 >= every other degree, since G has no universal vertex. q is a's only non-neighbour. Leaves at a = u plus the isolated vertices of G, all non-adjacent twins. Composing with a twin permutation fixes u, which leaves Aut(G)_q; labelling a first gives sign sgn_G. Exhaustive over all 16,824 (G, q) with n' <= 6; 501 random cases (369 with isolated vertices). A q in a K_2 component never occurs, because K_2 is odd. |
| 5 | Lemma 8 | CORRECT (Note 5) | (i): 2^k (2^rho - 2^h) + 2^rho - D(H) = 2^(m-1) - D(H). (iii): the twin-composition step is right when "even" means edge-sign even. Universal vertices of H other than z are exactly the extra degree-1 neighbours of a, and the composition handles them. B2.5 is the inversion identity. Exhaustive for rho <= 5 (5,405 cases). |
| 6 | Lemma 9 | CORRECT | D values by hand: F-a 1 + (rho - 1) + (rho - 3); F-b 1 + rho + (4 + (rho - 5) + 1); F-c as F-a. Aut(H)_z is trivial: the path flip moves an end, and moves the second vertex once rho >= 5. The spider's possible swaps at rho = 5 (legs 1, 1) and rho = 6 (legs 2, 2) move z. rho > h + 1 holds. Every m from 6 to 17 brute-forced (`indep_check_14.out`, section G). Negative control: F-c at rho = 4 is odd, which is why rho >= 5. |
| 7 | Exclusions | CORRECT | Lemma 1: at most one universal vertex. Lemma 2: H' is 1-regular, so n is odd. Lemma 3: max degree <= 2; pairs per component are r - 2 for P_r, 2 for C_4, k for C_k (k != 4). The classification and co-evenness failures (P_4 flip 3 edges / even permutation; C_4 rotation 2 / odd; P_3 end swap 2 / odd) were verified exhaustively through n = 9. {3, 7, 13}: Ore's bound gives d >= 2^ceil(c/2) for a connected component on c >= 2 vertices, so realisers have a component on <= 6 vertices. My n = 7 census finds exactly 21 / 35 / 210 labelled realisers (K_2, K_3, K_4 - e plus isolated vertices), all odd. |
| 8 | Theorem T and Lemma 10 | CORRECT | Tj(m) \ Tj(m - 2) = {m - 2} (m odd) or {m - 3, m + 1} (m even), checked for m < 200. Note that m - 1 is in Tj(m - 2) but not in Tj(m), which is harmless. SR: 3*2^(m-3) + 2^(m-3) - j = 2^(m-1) - j. Ranges: F-a and F-b need rho = m/2 >= 5, F-c needs rho = (m+1)/2 >= 6 for m >= 10. Bases m = 8, 9 re-verified. All top-band and offset L-graphs for m = 10..15 brute-forced. |
| 9 | Section 4 | CORRECT (Notes 6, 7) | Case 1: j <= n + 1; the missing j is n - 1 (n even) or n (n odd), and that x is M(n - 1). Case 2a: g = n - 2 qualifies; D odd and >= n + 3 forces g >= 3; d_G is odd and lies in [1, 2^g - g - 3]; the fallback value is d_G + 3*2^g - 2 >= 25 (since d_G >= 3 and g >= 3) and <= 2^(g+2) - g - 5; otherwise g = n - 2 and x = 3*2^(n-2) - 2 + d_G. Case 2b: y <= 2^(n-2) - 3; the remaining j lie in [3, n]; the missing j is n - 1 (n even) or n - 2 (n odd), and that y is M(n - 2). The top of O(n) is exactly {M(n)}, with X(n) excluded. |
| 10 | Scripts and independent check | PASS | See below. |
| 11 | Corollaries 2-4 | 2 CORRECT; 3 MINOR (Finding 4); 4 ERROR (Finding 2) | V(k) = least n with k in O(n) and ceil(log2(k+1)) <= V(k) <= ceil(log2(k+1)) + 1, checked for all odd k < 2^18. |
| 12 | Prior art | MAJOR + MINOR (Findings 1, 3) | No prior art found for the even-graph result (Finding 10). |

## Commands run (from `04-computation/experiments/egcp_20261109/` unless noted)

```
python verify_base.py      -> E(n) == O(n) for n = 1..9: True; L(8), L(9) cover T and offsets: True;
                              W witnesses checked: 248 | all base facts: PASS            (run_verify_base.out)
python lemma_check.py      -> complement 300/0, Lemma 5 300/0, Lemma 7 300/0, Lemma 8 400/0 (206 even)
                                                                                         (run_lemma_check.out)
python proof_check.py 20   -> all cases verified for 10 <= n <= 20 (per-n branch counts identical to
                              proof_check_24.out)                                       (run_proof_check_20.out)
python cover.py 20         -> missing E: 0 for every 10 <= n <= 20                       (run_cover_20.out)
python cover.py 25 --notop -> missing E: 0 for every n <= 25 (closure without the injected top band)
python cover.py 16 --noleaf / 14 --notop --noleaf -> misses appear (515, 519, 525, ...; 501-511): SR and the
                              top band are both load-bearing in the closure
```

From `05-knowledge/results/egcp_20261109/audit/` (audit code; nothing imported from egcp_lib):

```
python indep_check.py 14   -> A machinery: 410 graphs, generator evenness == n!-permutation evenness, subset-DP d
                              == inclusion-exclusion d; B-G lemma checks 0 failures (Lemma 9 table m = 6..17);
                              H base: 248 W (set == O(9)), 20 L, re-verified; 72 L-graphs on 10..15 vertices;
                              I end-to-end: all 256 / 511 / 1024 / 2047 / 4096 values needing exactly
                              n = 10 / 11 / 12 / 13 / 14 vertices built by the Section 4 recursion and verified;
                              TOTAL FAILURES: 0                                         (indep_check_14.out)
python indep_census.py 7   -> all 2,097,152 labelled graphs on 7 vertices (and n <= 6): every d odd, E(n) == O(n)
                              for n = 1..7; realised only by odd graphs: {3: 21, 7: 35, 13: 210, 119: 4130,
                              123: 210, 125: 21, 127: 1} at n = 7, {25: 115, 29, 31} at n = 5; Remark "no
                              universal vertex" confirmed for n = 5, 6, 7 (and base_sets.json n = 5..9);
                              corollary arithmetic (Finding 4)                         (indep_census_7.out)
python exhaustive_small.py -> Lemma 8 all 5,405 (H, z), rho <= 5: 0 failures, hub always unique for rho >= 2;
                              Lemma 7 all 16,824 cases, n' <= 6: 0 failures; Lemma 3 shapes through n = 9:
                              P_3 (+ sK_2) for one pair, C_4 / P_4 / 2P_3 (+ sK_2) for two, co-even 0
                                                                                         (exhaustive_small.out)
python cor4_unions.py      -> P_3 u G_7: n = 10, d = 505, even = True; union-realisable top-band counts per m
                                                                                         (cor4_unions.out)
python twohub_probe.py; python twohub_parity.py -> identification of cover.py's two-hub family; TH(G, t) even
                              iff G even except (f, t) = (1, 1)               (twohub_probe.out, twohub_parity.out)
```

## Sources consulted (prior art)

- A. E. Brouwer, P. Csorba, A. Schrijver, preprint 2009: https://aeb.win.tue.nl/preprints/domin4a.pdf (Theorem 1.1,
  Propositions 1.1-1.2, closing remark on the neighbourhood complex of the complement).
- S. Alikhani, D. Fatehi, S. Klavžar, arXiv:1512.07514, Proposition 2.1(iii).
- S. Alikhani, D. Fatehi, arXiv:1410.1169, Introduction (Corollary 2, citing Alikhani 2013).
- S. Alikhani, E. Deutsch, arXiv:1305.3734, Introduction (states Alikhani 2013's result as D(G, -1) = n).
- S. Alikhani, Graphs Combin. 29 (2013) 1175-1181, DOI 10.1007/s00373-012-1211-x (full text not accessible).
- I. Heinrich, P. Tittmann, arXiv:1701.03453 (D(G, x) + N(comp G, x) = (1 + x)^n).
- J. I. Brown, R. J. Nowakowski, Australas. J. Combin. 42 (2008) 55-68 (Lemma 3, Theorem 1).
- G. Royle et al., arXiv:2204.01947 (even graphs; no domination content).
