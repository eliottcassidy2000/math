# Audit E report — THM-4655, THM-4656, HYP-9279, strong-tournament note (mac-mini-2026-10-14-tourn)

Independent adversarial audit (subagent, 2026-10-14). The subagent could not write files outside its scripts, so the coordinator
saved this report verbatim from the agent's final message. All audit scripts and their `.out` files are in this folder. Nothing
was committed by the auditor and no audited file was edited by it. Tools: python3, numpy, sympy, gp, gentourng (installed);
literature by WebSearch/HTML only, no downloads.

Scope: `448c0d367e` plus the updates through `93231731e6`: THM-4655 1(e)–(g), statement 4 at n ≤ 10 and statement 5;
HYP-9279 (B), (B⁺), sharpness, (H), (B_loc) and (Z); note §1.5–1.8. Later note edits were not audited.

## Verdicts

| File | Verdict | Why |
|---|---|---|
| THM-4655 | **CORRECT WITH FIXES** | No mathematical error. One literal gap in 1(d), a wrong threshold in the status line, one unproved "for every m" in 1(e), attribution and wording. n = 10 rests on one program (`tbal2`). |
| THM-4656 | **CORRECT WITH FIXES** | Every formula and list checks out. Prior art uncited for part 1's unit criterion and the part 3 and 4 lists. One parenthetical in the proof of part 3 false in general. Mersenne case misses q = 2. The Fermat analogy needs sharpening. |
| HYP-9279 | **CORRECT WITH FIXES** | All data for n ≤ 9 reproduced. Four statements wrong or misattributed (M1–M4). |
| Results note | **CORRECT WITH FIXES** | Same items, plus one wrong reason, one imprecise claim, untyped rhetoric, and a conjectural finiteness stated as fact. |

## MAJOR findings

- **M1.** HYP-9279(A) is HYP-9277 statement 3 (opus S35) word for word ("every strong tournament on n ≥ 4 vertices has δ ≥ 4/11"). Fix: "Part (A) is HYP-9277(3) (opus S35, evidence n ≤ 8); this file adds the n = 9, 10 evidence and the conjecture that equality holds only at n = 5, 6, 7."
- **M2.** The "no poset analogue" example is wrong: in {a} ∥ {b < c} every element has a balanced partner (P(b before a) = 2/3, P(c before a) = 1/3). Fix: in 3+1 = {x < u₁ < u₂} + {y}, x is incomparable to y but P(x before y) = 3/4, so x (and u₂) has no balanced partner. Four elements are needed: no poset on ≤ 3 elements fails; 48 of the 219 labelled 4-element posets do.
- **M3.** "The tight vertices sit in a cyclic triangle hung at the end of the order by one back arc" is false in general. True for all 32 tight vertices with n ≤ 8: each lies on a cyclic triangle with two of its exact partners. One back arc suffices only for the transitive-plus-one examples (THM-4655's equality word and its reverse), which include every tight vertex at n = 5, 6. At n = 7, 5 of 7 tight vertices need two back arcs; at n = 8, 18 of 20, and the triangle can sit mid-order. Tight counts 1, 2, 5, 17, 82 (n = 5..9) are correct.
- **M4.** HYP-9279(Z): "for k = 2 the best-partner balance of z stayed ≥ 6/17" is contradicted by exhaustive data. Exhaustive minima for k = 2 cut vertices: 5/13, 6/17, 10/29, 18/53, 34/101 (n = 5..9); for n ≥ 6 this is (2^(n−4)+2)/(3·2^(n−4)+5) → 1/3 (witnesses: `z2_k2_witness_check.out`). (Z1), (Z2), (B_loc) hold for all n ≤ 9.

## MINOR findings

THM-4655:
- **T1** 1(d) literally includes transitive tournaments. Use: "Let T be a non-transitive tournament in which some strong component on at least 3 vertices becomes transitive after deleting one vertex. Then δ(T) ≥ 1/3."
- **T2** Status line "analytic bound for m >= 8" must be "m >= 9".
- **T3** 1(e) "two is attained for every m ≥ 3" needs a witness: w = 1 0^(m−1) (H = 2^(m−1)+1) has exactly the balanced cycle arcs (u₁,u₂), P(u₂ before u₁) = 2^(m−2)/(2^(m−1)+1), and (u_m,z), P(u_m before z) = (2^(m−2)+1)/(2^(m−1)+1); the others have P(z before u₁) = 2^(m−1)/(2^(m−1)+1), P(u_(i+1) before u_i) = 2^(m−3)/(2^(m−1)+1) (exact, m = 3..16).
- **T4** "tpo_u2_check.c independently enumerates all valid A" → "enumerates all valid A (independent of the closed-form decomposition, but resting on the bijection 1(b), which audit E checked against DFS enumeration of the actual Hamiltonian paths for all 16,382 words with m ≤ 13)"; same for the note's "brute force for m ≤ 20".
- **T5** Cite: for w = 1^(m−1)0 this is the nearly transitive tournament, H = 2^(n−2)+1 (El Sahili–Ghazo Hanna, arXiv:2207.11510).
- **T6** Statement 4 at n = 10 rests on `tbal2` alone; the audit replicated n ≤ 9 and checked the n = 10 minimiser (δ = 62/157, H = 1413).
- **T7** (optional) all 6 strong tournaments at n = 5 also fail for gaps; at n = 7, 351 of 353 fail for positions.
- **T8** No committed outputs for `trans_plus_one.py`, `rotation_lemma_check.py`, `logconcave_test.py`.

THM-4656:
- **W1** Part 1's unit criterion is KNOWN: Bragman–Rowland 2022 (arXiv:2202.00704, Thm 10), McIntosh–Roettger 2007 (same lesson as MISTAKE-614 E2). "L_p ≡ 1 (mod p²)" is due to Andrejić 2006 and McIntosh–Roettger 2007, not Sun–Sun 1992 (FLT connection). Cite Akin 2018 (arXiv:1806.11241) for the invariant-subset count. Type part 1's first sentence KNOWN (cited); the count transfer is an assembly.
- **W2** Proof of part 3: "−1, the only element of order 2 in ⟨φ⟩ in both cases" fails at p = 11 (π = 10, φ⁵ = (1, −1)). Replace: split case φ^(π/2) = (x, y), x, y ∈ {±1}, xy = (−1)^(π/2) = 1 since 4 | π, so x = y = −1; inert case F_(p²) is a field.
- **W3** Mersenne: "q = 2 or q ≡ 3 (mod 4)" (M₂ = 3 has π = 8). One-line proof: for inert p, φ^(p+1) = N(φ) = −1, so φ^(2^q) = −1 and π(M_q) = 2^(q+1) — the necessity half of Lucas–Lehmer with seed 3.
- **W4** Golden rigidity (p | L_(2^j)) mirrors prime divisors of Fermat numbers (ord_p 2 a power of 2), not Fermat primes; Fermat primes are mirrored by inert p with p + 1 = 2^q. The list 3, 7, 47, …, 524287 is OEIS A074714; type KNOWN. Same fix in note §2.
- **W5** Golden Paley primes = Fibonacci-primitive-root primes ≡ 3 mod 4 (Shanks 1972; OEIS A003147; lists agree below 2000, 52 primes). "Aut acts regularly on the arcs" holds for prime order p.
- **W6** Title: for split primes part 1 counts tournaments on Z_(p²), not Z[μ]/(p²); retitle or adopt S3.
- **W7** `wss_orbits.out` section (3) line "every golden-Fermat prime divides some L_(2^j): False" is a buggy check (factor lists stop at j = 6, wrong fallback); `golden_fermat_check.out` gives 0 mismatches.
- **W8** Checks list omits a = 4: 3; a = 5: 3, 11; a = 7: 5.
- **W9** Cite Freyd–Brown 1992 (AMM E3410) for π(n) ≤ 6n.

Results note:
- **N1** §1.1 add "non-transitive" (as T1).
- **N2** Family minimum not monotone: 8/19 at n = 6 then 2/5 at n = 7.
- **N3** §2, K₆ − P₆: it is P₆ that has independence number 3, so K₆ − P₆ has clique number 3. Conclusion holds (generic rank 9, every edge-deleted subgraph independent: a 2D rigidity circuit; W₄ too).
- **N4** §3: the resistance is a ratio, the projective action of Q² = [[2,1],[1,1]], with a shorter period (5, not π(5) = 20, at n = 5). Type "fixed point of least resistance … fails modulo p" as ANALOGY.
- **N5** §4: "PSL₂ over the residue fields F₄ = Z[φ]/(2) and F₅ = Z[φ]/(√5)"; delete "the projective image of the golden ratio"; type the E₈ / half-spin 128 link NUMEROLOGY.
- **N6** §5: rigidity facts correct, but the "11 coordinates and 55 pairs" reading is ANALOGY, not STRUCTURAL.
- **N7** §6: "finitely many extra cycles of 3n−1" → "the three known cycles (conjecturally all)".

## Suggested strengthenings (proved or exactly checked)

- **S1** δ(T) > 1/3 on the whole strong transitive-plus-one family (m ≥ 3). If w ≠ 1 0 1^(m−3) 0, (z,u₂) is strictly balanced by 1(c). At the equality word with m ≥ 5, P(u₃ before z) = 2^(m−4)/(1+2^(m−3)) ∈ (1/3, 1/2) (valid A containing u₃: 3·2^(m−4); H = 3(1 + 2^(m−3))). At m = 4 (w = 1010), P(u₁ before u₄) = 5/9.
- **S2** HYP-9279(B) is sharp for every n ≥ 5: at w = 1 0 1^(m−3) 0, u₂ has P(u₂ before z) = 2/3, P(u₂ before u₁) = 1/3, P(u_j before u₂) = 2^(m−4)/(3(1+2^(m−3))) < 1/3 (3 ≤ j ≤ m−1), P(u_m before u₂) = 2^(m−3)/(3(1+2^(m−3))) < 1/3; checked m ≤ 16. Within the family exactly this word and its reverse attain it for each m ≥ 5.
- **S3** THM-4656 part 2's formula also holds for split primes on the full ring Z[μ]/(p²): for every p ∤ 2D the number of orbits is (p²−1)(p^(2−ε)+1)/π', no orbit fixed by negation (81 split cases p < 200; split WSS primes 31, 191, 643, 1546463).

## Independent verification (auditor's own code, different methods)

- `tpo_dfs_bruteforce.c`: actual tournaments, DFS over Hamiltonian paths, all 16,382 words m = 1..13 (4,095 strong): 0 mismatches in strongness, the path ↔ A bijection, the H formula and K-form, and 1/3 < P(u₂ before z) ≤ 2/3 with equality only at 1 0 1^(m−3) 0; family minima 2/5 … 400/817 match; argmin of P(u₂ before z) is 1 1 0^(m−2). Every algebraic step of 1(c), including Lemmas A, B, re-derived by hand.
- `tpo_extra_check.c` (m ≤ 14, 8,190 strong words): 1(e)–(g) 0 failures; 2X, 2Y ≤ W_M in all 38,144 cases; the m ≥ 9 tail re-derived.
- `rotation_audit.py`: all cyclic arrangements n ≤ 9; bound attained only by the identity cycle; corollary checked.
- `logconcave_audit.py` (own digraph6 decoder): positions fail 1/1, 4/6, 35/35, 351/353; gaps fail 1/1, 6/6, 35/35, 353/353 (n = 4..7).
- `strong_dfs_balance.c`, `strong_dfs_extra.c`, `strong_dfs_bloc_z.c`: DFS on gentourng `-z -c`; lists complete and irredundant (Σ n!/|Aut T| = labelled strong counts, n ≤ 9). Strong counts 1, 1, 6, 35, 353, 6,008, 178,133 (n = 3..9); min δ 1/3, 2/5, 4/11, 4/11, 4/11, 62/157, 62/157; # minimisers 1, 1, 1, 2, 1, 1, 2; min feedback arc set of minimisers –, 1, 2, 3, 5, 4, 6; transitive-plus-k –, 1, 1, 2, 3, 2, 3; min P(tail first) –, 3/5, 5/9, 7/15, 9/25, 13/45, 17/81; tight (B) counts –, –, 1, 2, 5, 17, 82; balanced graph with a degree-1 vertex –, –, 1, 4, 17, 142, 1,915. (B), (B⁺), (B_loc) hold for n ≤ 9. Hamiltonian cycles: all 133,632 at n = 8 carry a balanced arc (min 1, in 8 tournaments); of 15,013,970 at n = 9 exactly one cycle in one tournament has none. Shorter cycles without a balanced arc: none to n = 6; n = 7: two 4-, four 5-, two 6-cycles; n = 8: two 3-cycles. (Z): 99,695 cut-vertex instances at n = 9 (13,060 with k ≥ 3), no failures; (Z1) minimum exactly 1/3; z always has a partner in D₁∪D₂∪D_(k−1)∪D_k.
- Labelled sweep n ≤ 7 without gentourng: same minima; (B) holds.
- `spot_checks.out`: the n = 9 (H) counterexample matches every stated detail; the n = 10 minimiser has δ = 62/157.
- `component_reduction_check.py`: δ(T) = max over components and H(T) = ∏H(C) on all 33,866 labelled tournaments n = 2..6.
- `cutvertex_chain_audit.py`: THM-4655 (5) re-implemented and compared with DFS on every strong tournament n ≤ 8 with a cut vertex (3,936 tournaments, 5,673 (T, z) pairs): 0 mismatches in H, v-before-z counts and the fair-coin identity.
- THM-4656: Burnside orbit counts (π by direct iteration) in 76 small cases (a = 1..7) plus 1546463 (3,092,928 orbits vs 4), 241 (27,879,360 vs 116,160), 191, 643: 0 mismatches; a-WSS lists via p² | U_(p−(D/p)) identical; golden rigidity via PARI over 78,496 primes < 10⁶: list matches, 0 mismatches against the L_(2^j) criterion; φ^(2^q) = −1 mod M_q for q = 2, 3, 7, 19, 31, 107, 127; Lucas–Lehmer seed 3 agrees with primality for q ≡ 3 (mod 4), q < 300.
- `pisano_audit.c`: max π(n)/n = 6 exactly at 2·5^k (n ≤ 20,000); owner's examples hold.
- `rigidity_audit.py`: K₁₁ minimally rigid, K₁₂ a circuit in R¹⁰; K₆ − P₆, W₄, K₄ are 2D circuits.
- Small facts: Q⁵ ≡ 3I mod 5; ladder resistance F_(2n+1)/F_(2n); triangular Fibonacci numbers only 1, 3, 21, 55; 1² + … + n² square only at n = 1, 24 (to 2·10⁶).

## Prior art

- Tournament 1/3–2/3 for random Hamiltonian paths: nothing found; appears new. The case w = 1^(m−1)0 of 1(b): nearly transitive tournament count, https://arxiv.org/abs/2207.11510.
- WSS ⇔ golden Wieferich: https://arxiv.org/abs/2202.00704 (Thm 10), McIntosh–Roettger 2007; https://en.wikipedia.org/wiki/Wall%E2%80%93Sun%E2%80%93Sun_prime (Andrejić 2006, PrimeGrid 2^64, Pell–Wieferich 13, 31, 1546463). The tournament orbit-explosion reading not found.
- π(p) = 2^k ⇔ p | L_(2^j): classical; list https://oeis.org/A074714; https://en.wikipedia.org/wiki/Pisano_period (π(L_(2k)) = 8k; π(n) ≤ 6n, Freyd–Brown).
- Lucas–Lehmer seed 3: https://en.wikipedia.org/wiki/Lucas%E2%80%93Lehmer_primality_test
- Golden Paley primes: https://oeis.org/A003147
- Luo 1989: https://proofwiki.org/wiki/Triangular_Fibonacci_Numbers
- Gerbicz 2018: https://oeis.org/A090461
- Icosians → E₈: https://en.wikipedia.org/wiki/Icosian
- Wall 1960: only via a search summary of Klaška's survey (PDF not fetched).

## Not checked

The n = 10 exhaustive run (only its minimiser); original texts of Wall, Sun–Sun, Watson, Shanks; the icosian reductions; the motivic-weight remarks; the code of tbal2–tbal7, cutvertex_search and c_counterexample (claims re-derived with own programs instead); the (Z) sampling to n = 46; note edits after `93231731e6`. The (H) counterexample is confirmed and unique at n = 9.

## Claim-by-claim verdicts

- THM-4655 1(a) CONFIRMED; 1(b) CONFIRMED (add T5); 1(c) CONFIRMED (fix T4); 1(d) FIX (T1); 1(e) CONFIRMED, "two attained" FIX (T3); 1(f) CONFIRMED, status FIX (m ≥ 9); 1(g) CONFIRMED; 2 CONFIRMED; 3 CONFIRMED; 4 CONFIRMED n ≤ 9, n = 10 FIX (T6); "HYP-9277 for n ≤ 9" CONFIRMED; 5 CONFIRMED (n ≤ 8); Reading "≈ 3/4" CONFIRMED.
- THM-4656 setting CONFIRMED; part 1 CONFIRMED, cite KNOWN (W1); part 2 CONFIRMED (S3 extends); part 3 CONFIRMED, proof FIX (W2), list KNOWN (A074714); Mersenne CONFIRMED for q ≡ 3 mod 4, FIX (W3); Fermat mirror FIX (W4); part 4 CONFIRMED, FIX (W5); title FIX (W6); checks CONFIRMED, FIX (W7, W8); π(n) ≤ 6n CONFIRMED (cite W9).
- HYP-9279(A) data CONFIRMED, provenance FIX (M1); (B) CONFIRMED n ≤ 9; poset example REFUTED (M2); sharpness counts CONFIRMED; tight-structure description REFUTED (M3); (B⁺), 142 count, (H), (B_loc) CONFIRMED; (Z1), (Z2) CONFIRMED n ≤ 9; "k = 2 ≥ 6/17" REFUTED (M4).
- Results note CONFIRMED except N1–N7 and M3.
