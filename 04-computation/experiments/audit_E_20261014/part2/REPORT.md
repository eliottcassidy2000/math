# Audit E, part 2 — THM-4655 (6), (7) and the new HYP-9279 text (mac-mini-2026-10-14-tourn)

Saved verbatim in substance by the coordinator from the auditor's final message; scripts and outputs are in this folder.
Nothing was committed by the auditor and no audited file was touched; one niced core. The auditor read the on-disk versions
mid-update (before the statement-7 bound was raised from 6 to 7 everywhere).

**Overall:** no mathematical errors. Statement 6's values and statement 7's proof hold up. Two things needed fixing: the sharpness
claim was stated as proved beyond what statement 6 shows, and the block-size bound was inconsistent across files at read time.

## (1) THM-4655 statement 6 (family W_n, renamed V_n after this audit)

- Proof re-derived: vertex 0's only in-neighbour is b, so every Hamiltonian path starts at 0 (3 paths) or contains `b 0`; X, Y are the unique increasing paths on S, S′ in the transitive order `1…t−1, c, a`; the only exclusion is "c ∈ S and a ∉ S"; H = 3 + 3·2^(t−1), each rotation 1 + 2^(t−1) times; p(c before 0) = 2^(t−1)/H, p(c before i) = 2^(t−2)/H, both < 1/3. With z = b, u = 0…t−1, c, a, the family is the transitive-plus-one tournament with word 1 0^(m−3) 1 0, the complement-reverse of 1 0 1^(m−3) 0.
- Independent check (`w_family_check.out`, DFS path enumeration, n = 5..14): strong, H = 3(2^(n−4)+1), uniform rotations, p(c<a) = 2/3, p(c<b) = 1/3, stated p(c<0), p(c<i), c's balance exactly 1/3, W_n − b transitive, both identifications as labelled tournaments.
- **FIX NEEDED (overclaim).** The family proves only that the vertex-wise constant is **at most** 1/3; "equals 1/3" also needs (B), open beyond n = 10. Use: "The vertex-wise constant is at most 1/3 for every n >= 5 (W_n), and equals 1/3 wherever (B) holds, in particular for 5 <= n <= 10 (exhaustive)." Same change in HYP-9279's status and sharpness bullet.
- MINOR: the note uses W_4 for the wheel graph, so the name W_n collides.

## (2) THM-4655 statement 7: CONFIRMED, with consistency fixes

- N and H re-derived from statement 5 (left (1,g,P,Lp), right (Rp,s′,Q,1)): H = Rp + es′ + βQ + γ + gs + P(σQ+α) + Lp; N = Rp + es′ + β(v)Q + γ(v) + P(σQ/2 + α(v)); the expansions of 3N−H and 2H−3N match the script.
- (i), (ii), (iii) correct (new P = β + σP ≥ P; new Lp between old Lp and γ + s + (α+1)P ≤ β + σP; right side symmetric: new Q = σQ + α, new Rp ≤ (β+1)Q + e + γ ≤ σQ + α).
- Corner claim correct: affine in each of g, s′, Lp, Rp with no products, so the box minimum is at a vertex; each corner is aPQ + bP + cQ + d, a = σ/2; substituting P = 1+x, Q = 1+y gives coefficients a, a+b, a+c, a+b+c+d, so the four-condition test on [1,∞)² is exact in both directions.
- Remark: the comparisons β−β(v) ≤ σ/2, α−α(v) ≤ σ/2−1, α(v) ≤ σ/2, β(v) ≤ σ/2−1 hold for every block and give L2, L3, U2, U3 in general; at P = Q = 1 only L4 (worst at g = 1, s′ = 0) and U4 (worst at g = 0, s′ = 1) need the computer, matching `z1_local_slack.py`.
- Independent certificate `z1_cert_audit.c` (every permutation of the block and every cut point, exact integers; all 4 conditions at all 16 corners for both inequalities, plus (i), (ii), the fair coin and the four comparisons): sizes ≤ 6, 2,458 types, 0 failures; sizes ≤ 7, 47,642 types, 0 failures (agrees with `z1_cert_7.out`). Minimum L4/U4 slack by size: 1: 0; 3: 5/3; 4: 2; 5: 26/9; 6: 10/3; 7: 298/69.
- Brute-force validation `z1_bruteforce_random.c` (no chain; exact forward/backward subset DP): 3,000 random strong cut-vertex tournaments up to n = 16, interior blocks of size 1–7, end blocks up to size 7 incl. a large-end-block mode; 11,630 interior vertices checked, 0 violations; minimum exactly 1/3 (block sizes 1 3 1 1).
- FIX NEEDED: block-size bound inconsistent at read time (title/status 7 vs statement/(v) 6) — make consistent; (iii) base case: add "Since T is strong and D_1 has no in-arcs from other blocks, z has an out-neighbour in D_1; dually some vertex of D_k beats z"; Remark wording: "0 exactly for singleton blocks (realised in the equality-word tournament)"; `z1_cone_certificate.py` docstring describes Pmin, cL, Qmin, cR as minima over enumerated small blocks, whereas the proof uses the bound 1 from (ii) for arbitrary end blocks (computed values equal 1, results unaffected).

## (3) HYP-9279

- (Z) paragraph CONFIRMED (bound can become ≤ 7).
- "Beyond n = 10": all values CONFIRMED with own generator and DFS (block δ: S4 2/5, M5/M6a/M6b/M7 4/11, M8 62/157; all 14 pairs for n ≤ 13 enumerated over z-patterns, 31,662 strong tournaments, every stated minimum reproduced incl. 592/1521 (H = 1521) at n = 11, 1052/2721 at n = 12 and the n = 13 values; all seven n = 14, 15 minimisers spot-checked; "deleting suitable vertices keeps 592/1521 at n = 12": 6 of 13 deletions (`glue_delete_check.out`); (B) holds throughout).
- FIX NEEDED (wording): the tower bullet contradicts the n = 13 table entry — use "The extension tower first found an n = 13 tournament with δ = 592/1521 (H = 13689 = 117²); it is M6b + M6a glued at a cut vertex. The glue search then found 1052/2721 at n = 13."; the title says "to n = 13" but the table goes to 14 (plus n = 15).
- Part-1 corrections applied correctly; THM-4656 fixes correctly applied; the Fibonacci-primitive-root converse and the Mersenne q = 2 statement check out.

Not checked: annealing logs beyond quoted values; full z-pattern enumeration at n = 14–15 (minimisers only).

## Claim → verdict

- THM-4655 (6) values, uniformity, probabilities, balance 1/3, V_n − b transitive, reverse of the equality word → CONFIRMED
- THM-4655 (6) "equals 1/3 for every n ≥ 5" → FIX NEEDED (at most 1/3; equality where (B) holds)
- THM-4655 (7) derivation, (i)–(iii), corner reduction, finite check → CONFIRMED (independently, sizes ≤ 7, 47,642 types)
- THM-4655 (7) brute-force validation → CONFIRMED (3,000 tournaments, 0 violations, minimum exactly 1/3)
- THM-4655 (7) Remark → CONFIRMED (size 7: 298/69)
- THM-4655 bound consistency, (iii) base case, Remark wording → FIX NEEDED (minor)
- HYP-9279 (Z) → CONFIRMED; sharpness → FIX NEEDED; glued blocks → CONFIRMED; tower bullet and title → FIX NEEDED (wording); part-1 texts → CONFIRMED
