<!-- Independent adversarial audit of session opus-2026-10-14-S36, run by a separate agent with its own code (this folder). Saved verbatim from the auditor's final message by the session; the auditor could not write files outside its code. Findings E1-E12 and I1-I15 were applied in the session files; see MISTAKE-614. -->

# Referee report: opus-2026-10-14-S36

**Audit code.** Everything was written from scratch in `C:\Users\Eliott\AppData\Local\Temp\math-wt-strong-balance-pisano-36\04-computation\experiments\s36_audit_20261014\`; each `.py` has a matching `.out`. I did not import any session script or edit any session file. The numba caches are deleted.
- **Tournaments:** `tourn_core.py`, `audit_tournaments_labelled.py`, `audit_tournaments_classes.py`, `audit_tournaments_classes_n9.py`, `verify_specific_tournaments.py`, `verify_n9_ss_classes.py`.
- **Number theory:** `audit_pisano_wss_orbits.py`, `audit_golden_psl211.py`.
- **Packing and square-sum:** `audit_erdos106.py`, `audit_square_sum.py`.

**Methods.** I deliberately used methods different from the session's subset DP.
- Hamiltonian paths are enumerated by depth-first search. For n ≤ 6 this was cross-checked against a scan of all n! orders on every labelled tournament, with 0 mismatches.
- Strong connectivity of T uses Landau's score criterion; strong connectivity of F uses a Warshall closure.
- Transitivity of M uses the Kendall–Babington Smith 3-cycle count.
- Isomorphism classes for n = 8, 9 are canonically generated. All counts match OEIS: 6880 and 191,536 classes, labelled totals 2^28 and 2^36, 6008 and 178,133 strong classes, 236,522,496 and 64,026,088,576 labelled strong.
- Every run took under 6 minutes.

**Key independent results.**
- **n ≤ 7, all labelled strong tournaments.**
  - M(T) is never transitive, and F(T) is always strongly connected.
  - Inverted arcs: none for n ≤ 5; for n = 6, 7 exactly 1,440 and 43,680 tournaments have one.
  - Degree types: n = 6 only (1,1); n = 7 has (1,1) 23,520 times, (2,4) 10,080 and (4,2) 10,080. The squeezed-triangle property holds in every (1,1) case.
  - A forward 3-cycle always exists, and the minimum number of 3-cycles of M is n − 2.
  - (SS) has 0 violations. "Best in-neighbour + Pr[w first] ≥ 1" fails at 10,420,928 (T, w) pairs at n = 7, in every strong tournament with 4 ≤ n ≤ 7.
- **n = 8, exhaustive (fills [N8]).**
  - M is never transitive, F is always strong, and a forward 3-cycle always exists.
  - Inverted arcs per labelled tournament: 0 for 231,227,136; 1 for 4,932,480; 2 for 362,880.
  - (SS) holds at every vertex.
  - Extremes are 33/65, 33/65, 1/65, 1/65 and 32/65, each attained only by T_8.
- **n = 9, exhaustive.**
  - M is never transitive, and F is always strong.
  - At most 3 inverted arcs (4 classes).
  - The vertex lemma and D_in > 0 hold everywhere, and the minimum number of 3-cycles of M is 7.
  - (SS) fails in exactly 5 isomorphism classes (1,572,480 labelled tournaments). Every failure is at an in-degree-7 vertex with k = 3.
  - Extremes are 65/129 and 1/129 (and 64/129 for the first-vertex maximum), attained only by T_9.

## ERRORS

**E1. THM-4657 duplicates prior canon THM-4656.** mac-mini reserved THM-4656 ("Wall-Sun-Sun golden orbit explosion") at 2b6d2c3abf, before S36's d00166b48d. Its full text was pushed at 448c0d367e, and mac-mini sent a coordination message to S36 at 0336f594c7.
- THM-4656 Part 1 is THM-4657 Parts 1–2 verbatim, for any metallic unit μ_a.
- THM-4656 Part 2 is the inert case that THM-4657 leaves "not written out".
- THM-4656 Part 4 gives the same 52 golden Paley primes.
- None of the S36 files mention THM-4656.
- **Fix:** cite THM-4656 and either cede the number or keep THM-4657 only as an independent confirmation. Do not claim priority.

**E2. THM-4657 Part 1 ("Wall–Sun–Sun = Wieferich to base φ") is known.**
- Bragman–Rowland, arXiv:2202.00704 (2022), Theorem 10: |F(p−ε)|_p = |φ/ω(φ) − 1|_p, so p is Wall–Sun–Sun iff that quantity is ≤ p⁻². For split p this is equivalent to φ^(p−1) ≡ 1 (mod p²).
- "Fibonacci–Wieferich prime" is the standard alternative name (McIntosh–Roettger 2007; Wikipedia). Wieferich places in number fields are standard (Graves–Weiss 2025; Silverman 1988).
- **Fix:** type Part 1 as CITED. It appears in main note §0.3 and §2.2 and in HYP-9281 item 6.

**E3. Stale search bound.** "PrimeGrid … past 4·10^16 (OEIS A244801)" appears in note §2.2 and the THM-4657 status. The A244801 comment says exactly that, but it is out of date. PrimeGrid's combined Wieferich/WSS project (Nov 2020 – Dec 2022) showed none below 2^64 ≈ 1.8·10^19 (Wikipedia). THM-4656 already uses the correct bound.

**E4. "|G_n| = nL_n" is false for even n** (note §3.1; HYP-9281 item 2).
- The correct order is |R_n| = L_n − 1 − (−1)^n, which is L_n − 2 for even n. For n ≥ 3, |G_n| = n|R_n|.
- So |G_4| = 20, not 28, and |G_6| = 96, not 108. The paper's Table 4 and my run agree.
- Regularity forces |R_n| = 2n + 1. For odd n that means L_n = 2n + 1, so n = 5. For even n it means L_n = 2n + 3, which never happens.
- "Equivalently" in HYP-9281 item 2 is not an equivalence.
- The norm identity, the |R_n| formula and the classification are the paper's own Lemma 5.1, Table 2 and Theorem 6.4. They must be credited, not presented as new.

**E5. "The eleven-square torus is the icosahedron with a marked vertex" overclaims** (§0.5 and §3.2 in bold; HYP-9281 title).
- The colours are described correctly: the paper's Fig. 1 and Fig. 3 show H = {1,3,4,5,9} gold, −H teal and 0 dark gray.
- The paper's "mark" is a boundary point of the flat torus: the canonical lift z ∈ Z_w, which encodes the side T (Prop. 3.4, Thm 3.5, the dotted segment to R_w z). It is not a 12th cell or label.
- The torus cell-adjacency graph is the 6-regular circulant C_11(1,3,4): neighbours ±a, ±b, ±(b−a) carry labels ±1, ±4, ±3. The icosahedron minus a vertex has degrees 4 and 5.
- A_5 does not act on the 11 cells. The only structure matched is the C_5-orbit pattern {0}, H, −H.
- **Fix:** retype as ANALOGY and say exactly what matches.

**E6. §0.1 "Every sharp constant found is attained by T_n" is false.** The session's own `tournament_strong_balance_rules` output gives the least balance δ over strong tournaments: 4/11 at n = 5, 6, 7. My runs give 62/157 at n = 8, 9. None of these extremes is attained by T_n (0 of the extremal tournaments are T_n). **Fix:** restrict the statement to the 3-cycle certificate, vertex lemma, D_in, D_out and first-vertex constants.

**E7. The one-inverted-arc case is overclaimed.** §0.1 says "One inverted arc never suffices to break the conjecture" and §1.4 says "The conjecture's conclusion holds whenever T has at most one inverted arc". The HYP-9280 proof is correct, but it only gives that F ⊇ T − e contains a directed cycle. That is the weak form, which yields a balanced pair. It does not give that F is strongly connected, which is the conjecture itself. **Fix the wording accordingly.**

**E8. §1.3b "That gives E[N] ≥ d/2, and > d/2 by parity" is wrong for even d.**
- H odd only rules out a single p(u,w) = 1/2. It does not rule out E[N] = d/2: when d is even, Σ C(u,w) = dH/2 is possible.
- HYP-9280 states this correctly: max_u p(u,w) ≥ 1/2, hence > 1/2.
- E[N] > d/2 does hold for all n ≤ 9 (minimum D_in = 1/H(T_n)), but that comes from the data, not from parity.

**E9. Several FINITE-EXACT claims are not computed by any shipped script.**
- Not computed: "F(T) strongly connected for n ≤ 7" (§0.1, §1.2, HYP-9280), the inverted-arc degree types, and "millions of violations". The cited `tournament_majority_cycle_20261014.py` returns only M-transitivity, the inverted-arc count and the forward 3-cycle.
- Missing outputs: the annealing script has no `.out`. There is no output for "25,784 n = 8 samples", and no script or definition for "2,427 / 16,426 pointwise tests".
- All of the computed claims are true per this audit, and the n = 8 sample claim is superseded by the exhaustive run.
- **Fix:** ship the checks (or cite the audit `.out` files). The pointwise claim is unverifiable as shipped.

**E10. Packing literature.**
- Chung–Graham's 2020 O(x^(3/5)) argument contains an error, according to McClenagan (arXiv:2602.01484). The 3/5 bound is proved by McClenagan (2026) and, by a different method, by Bui (arXiv:2508.04603).
- The Roth–Vaughan bound "W(x) >> (x − ⌊x⌋)^(1/2)·…" should read W(x) > 10^(−100)·√(x‖x‖) when x(x − ⌊x⌋) > 1/6.
- "W. Wang" should be Shuang Wang (with T. Dong and J. Li).

**E11. Pisano edge cases.**
- §2.1 says inert primes have π(p) | 2(p+1) but not p+1, "checked for every prime to 2000". This fails at p = 2, since π(2) = 3 divides p + 1. The session's script starts at p = 7.
- HYP-9281 item 7 "equality iff n = 2·5^k" needs k ≥ 1, because π(2) = 3 ≠ 12.

**E12. §3.2 "the Paley maximisers of H (A038375: 3, 7, 11) are exactly the Lucas numbers … is exact."**
- A038375 is known only up to n = 11.
- The repo's THM-135 shows P_19 is not even the best circulant on Z_19.
- For q ≥ 23 this is open.
- **Fix:** say "within the known range (q ≤ 19)" and type the coincidence as NUMEROLOGY.

## IMPRECISIONS

- **I1. First-vertex proof wording.** "The end is always valid" is false: it needs x_(n−1) → w. A slot exists anyway: take the end of the maximal run of in-neighbours. Alternatively, reuse the k = 0 injection (move w to just after the last in-neighbour).
- **I2. "Exact factorisation" misused.** "PSL(2,11) = A·B exactly / (exact factorisation)" is not an exact factorisation in the group-theory sense, since A ∩ B = C_5; the exact one is A_5·C_11. This factorisation is classical (it is Galois' 11-point action; factorisations with an A_5 factor are tabulated, e.g. by Liebeck–Praeger–Saxl 1990). §7 item 4 should not claim this reading as new content.
- **I3. Golden rectangles.** The 12 points carry 15 golden rectangles (any two antipodal pairs), forming 5 orthogonal triples. "The three" requires a choice.
- **I4. Galois statement.** "Acts on p points only for p ≤ 11" should be "has a subgroup of index p". The cases p = 2, 3 also have one (and 3 is a Paley prime). "Hence the golden rectangle at 11, and the Klein quartic at 7" is an association, not an inference.
- **I5. Typing in §3.1.** "…all come from the Fibonacci fixed point F_5 = 5" sits under a "(PROVED, elementary)" header but is interpretive. mac-mini's THM-4656 types 55 = F_10 as a COINCIDENCE (Luo). These need reconciling.
- **I6. HYP-9281 predictions.** These confuse L_n with |R_n|. R_4 ≅ F_5 (with G_4 = AGL(1,5)), not a 7-point object, and R_2 = 0.
- **I7. Square-sum credit.** The full result is known (Gerbicz 2018, OEIS A090461): orderings exist for every n ≥ 25 and cycles for every n ≥ 32. The note's "to 40" should cite this.
- **I8. Pisano attribution and wording.** Attribute the 6n bound (Freyd–Brown 1992 / K. S. Brown, OEIS A001175). "Highly composite" is misused.
- **I9. Which unit is invariant.** "Golden-invariant" should mean invariant under the odd part g of φ. φ itself has no invariant tournament when its order is even. At Fibonacci Paley primes, ψ = −1/φ is a primitive root.
- **I10. Erdős 106 credit.** erdosproblems.com credits "Claude Opus 5 (prompted by Silverstein)" and links Sprite143's repository. The rational rotation (143/145, 24/145) is the paper's own refinement. The paper also cites Baek–Koizumi–Ueoro for the axis-parallel case.
- **I11. s(11) status.** Wikipedia reports an AI-assisted optimality proof dated 29 Sep 2026; describe it as "reported, under review".
- **I12. Poset analogue.** Cite the linear-extension-majority digraph literature (Kislitsyn; Fishburn 1974/76; Ewacha–Fishburn–Gehrlein 1990; De Loof–De Baets–De Meyer; Gupta 2026, whose longest LEM cycle at order 14 is 8, consistent with HYP-9280). I found no prior tournament version, so the majority tournament appears new.
- **I13. Cross-reference mac-mini.** HYP-9279(B), "every vertex has a balanced partner", is the session's own R3 check (0 failures for n ≤ 7, unreported). THM-4655 covers the transitive-plus-one family, which contains T_n.
- **I14. Update HYP-9280 ranges.** Change "FINITE-EXACT n ≤ 7" to n ≤ 9 and record that the formula is now verified at n = 8, 9. "D_in > 0 … equality exactly at T_n" is a wording slip: the minimum 1/H(T_n) is attained only at T_n. In §1.3b, "now hinges on" should be "would follow from".
- **I15. (SS) status.** Update to: holds for all strong n ≤ 8, fails at n = 9 in exactly 5 classes. The session's two examples are classes 2 and 4 of the 5.

## VERDICTS

| Claim | Verdict |
|---|---|
| HYP-9280 reductions: 3-cycle of M ⇒ balanced pair in (1/2, 2/3]; F strong ⇒ M strong ⇒ ≥ n−2 cyclic triples (Moon); strong majority ⇒ an arc with p < 1/3 in every 3-cycle; vertex lemma ⇒ F has a cycle; (SS) ⇒ vertex lemma | CONFIRMED |
| First-vertex bound | CONFIRMED WITH FIXES (I1) |
| In-degree-1 lemma; (SS) at k = 0; vertex lemma for in-degree and out-degree ≤ 2 | CONFIRMED |
| At-most-one-inverted-arc case | CONFIRMED in the weak form; note wording REFUTED (E7) |
| n ≤ 7 exhaustive data and (SS) for n ≤ 6 (brute force) | CONFIRMED |
| Sharp T_n family and its uniqueness | CONFIRMED for n ≤ 9 |
| The two n = 9 (SS) counterexamples (H = 85 and 117, E[N] = 307/85 and 419/117, best p = 53/85 and 23/39) | CONFIRMED by brute force over 9! orders |
| HYP-9277 Claim A data (13/45, 17/81; one inverted arc; 10 and 12 three-cycles of M; both strong) | CONFIRMED |
| Annealing claims | NOT CHECKED (no output shipped); consistent with exhaustive n = 8, 9 |
| THM-4657 proof and numerics (40 pairs; the explosion branch tested on 39 artificial Wieferich-type units) | CONFIRMED; novelty and priority REFUTED (E1, E2), so CONFIRMED WITH FIXES |
| Pisano facts and table | CONFIRMED WITH FIXES (E11) |
| No Wall–Sun–Sun prime below 2·10⁶ (two tests, 148,930 primes) | CONFIRMED; literature bound needs fix (E3) |
| 52 Fibonacci Paley primes below 2000 | CONFIRMED |
| HYP-9281 identities: norm, L_n = 2n+1 ⇔ n = 5, C(L_n,2) = F_2n ⇔ n ∈ {2,4,5}, F_n = n ⇔ n ∈ {0,1,5}, Luo | CONFIRMED (credit the paper) |
| \|G_n\| = nL_n | REFUTED (E4) |
| PSL(2,11) facts: 660; Borel = Aut(P_11) regular on 55 arcs; two non-conjugate A_5 ∋ 4x, transitive with stabiliser C_5; A ∩ B = C_5; G = AB; orbitals 12/60/60; icosahedron with ∞ antipodal to 0 and pentagons QR/NQR | CONFIRMED (terminology I2) |
| Torus = icosahedron with a marked vertex | REFUTED as structure; acceptable only as ANALOGY |
| Erdős 106 exact check (independent clipping-area method; side sum 4 + 101/1153304; both contacts exact) | CONFIRMED |
| Square-sum to 40 | CONFIRMED (credit I7) |
| Square-packing literature | CONFIRMED WITH FIXES (E10) |

Sources:
- [erdosproblems.com/106](https://www.erdosproblems.com/106)
- [Wikipedia: Wall–Sun–Sun prime](https://en.wikipedia.org/wiki/Wall%E2%80%93Sun%E2%80%93Sun_prime)
- [OEIS A244801](https://oeis.org/A244801)
- [Bragman–Rowland, arXiv:2202.00704](https://arxiv.org/html/2202.00704v2)
- [Graves–Weiss, arXiv:2503.19144](https://arxiv.org/html/2503.19144v1)
- [McClenagan, arXiv:2602.01484](https://arxiv.org/html/2602.01484v1)
- [Wang–Dong–Li, arXiv:1603.02368](https://arxiv.org/abs/1603.02368)
- [Wikipedia: Square packing](https://en.wikipedia.org/wiki/Square_packing)
- [Gupta, arXiv:2607.23926](https://arxiv.org/abs/2607.23926)
- [OEIS A038375](https://oeis.org/A038375)
- [OEIS A090461](https://oeis.org/A090461)
- [Wikipedia: Pisano period](https://en.wikipedia.org/wiki/Pisano_period)
- [Wikipedia: Wieferich prime](https://en.wikipedia.org/wiki/Wieferich_prime)
