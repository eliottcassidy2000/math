# Independent adversarial audit, wave 2 (audit_wave2)

*Saved by the main session (opus-2026-10-22-ssj). The auditing subagent's own write of this file was refused by its session, so it returned the text and this file records it with light formatting. Its scripts and outputs are in this directory. The downloads stay local (`downloads/` is git-ignored); their sha256 prefixes are at the end. Corrections applied by the main session are listed in the last section.*

**Audited:** `05-knowledge/results/monster_levels_square_sum_7_21_20261022.md` (commit 211ab1a1c1, "pre-audit"), together with `level_ogg_lists.gp`/`.out` and `square_sum_structure.py`/`.out`.

**Auditor:** temporary identity `audit_wave2`, 2026-10-22.

**Rules kept:**
* The note and the scripts were not edited, and nothing was committed.
* PARI stacks were at most 400 MB; `parisizemax` was never set; jobs ran one at a time.
* `mfatkininit` was run only at levels ≤ 150.
* Downloads are data only.

**Method.** Each claim was checked by a route different from the one the note used:
* Riemann–Hurwitz with CM fixed-point counts, in place of newform Atkin–Lehner signs;
* direct computation of the supersingular points of X_0^+(N) mod p, in place of the genus criterion;
* the auditor's own factorisations of the sporadic orders;
* BFS, exhaustive path counts and an exhaustive cycle search;
* a C enumeration of all 7-vertex tournaments;
* the downloaded paper itself.

## Summary verdict: PASS WITH CORRECTIONS

| # | Severity | Location | Verdict |
|---|---|---|---|
| F1 | medium | §0 item 3, §4.1 | ERROR: contradicts its own output |
| F2 | medium | §0 item 3, §4.1, §4.3 | OVERCLAIM: missed inheritance of THM-4510 |
| F3 | medium | §0 item 3(a), §4.2 | ERROR: "k a Pell number" |
| F4 | medium | §0 items 1–2, §1, §3, §6 | WORDING: B is not a subgroup of M |
| F5 | low–medium | §5 row E, §0 item 4 | ERROR: the test sentence is false |
| F6 | low–medium | §5 row C, §0 item 4 | OVERCLAIM: the nesting is automatic |
| F7 | low | §5 row A | ERROR: wrong growth constant |
| F8 | low | §5 row F | WORDING |
| F9 | low | §5 Inheritance | WORDING; THM-1370 number collision in the canon |
| F10 | low | §5 row D | WORDING: the free action holds for every tournament |
| F11 | low | §0 item 1 | WORDING: attribution |
| F12 | low | §1 | WORDING |
| F13 | low | status line, §3 | WORDING: missing or weak hostile tests |
| F14 | trivial | §8, §4.2 | reproducibility nit |

No error reaches a PROVED or FINITE-EXACT result.

## Findings and repairs

### F1. The Pythagorean-size hostile probe contradicts its own output
**Claim.** §0 says "a Pythagorean size triple occurs in no other sum-graph family tested". §4.1 says it "occurs only for squares (once, at n = 12)".

**Finding.** The script's criterion is that some three of the component sizes form a Pythagorean triple. Under that criterion, the script's own `.out` (part 5) lists two more hits, both for twice-squares:
* n = 24, sizes [2,3,4,5,10];
* n = 25, sizes [2,3,4,5,11].

`audit_square_sum.py` §5 reproduces this on 12 families for n < 120. Squares at n = 12 are unique only under the reading "exactly three components".

**Repair.** "Exactly three components of Pythagorean sizes occur, among the families tested, only for the squares at n = 12. For twice-squares, the sizes at n = 24 and 25 contain {3,4,5} among five components." Fix §0 to match.

### F2. "Only structural link" ignores THM-4510
**Claim.**
* §0: "Their only structural content is that the three components are welded by the square 25".
* §4.1: "the only structural link found to the '3, 4, 5 right triangle'".
* §4.3 presents its own extension and loss mechanism.

**Finding.** THM-4510 is cited in §4 as inherited, so its content counts:
* Corollary A3'' proves that the whole first window {15,16,17} is the single tight triple of squares (9,16,25) = (3²,4²,5²).
* A3'' also shows the window closes at 18 because 7 gains a third partner.
* A3''' generalises this to all primitive Pythagorean triples.
* The 15-chain uses exactly the squares 9, 16 and 25 (wave six, S5).

So a stronger 3-4-5 link is already in the repo, and §4.3 restates A3'' without credit. The welding at 13 and 14 uses both 16 and 25 (13+3, 13+12, 14+2, 14+11); 25 only reaches the component {4,5,11,12}.

**Repair.** Cite THM-4510 A3'' and A3''' as the structural content. Write "welded by 16 and 25". Keep NUMEROLOGY for the sizes 3, 4, 5.

### F3. "k a Pell number" is false
**Claim.** §0 3(a): the coincidences are at "k a Pell number". §4.2: "59 and 71 are not Pell numbers".

**Finding.**
* The Pell numbers are 0, 1, 2, 5, 12, 29, 70, 169, 408, 985.
* The values k = 3, 17, 99, 577 are not Pell numbers. They are even-index half-companion Pell numbers, the numerators with p² − 2q² = 1.
* The k-set 1, 3, 5, 17, 29, 99, 169, 577, 985 is OEIS A079496: the "integer square roots of floor(n²/2)+1". That is exactly the coincidence condition, since ⌈(j²+1)/2⌉ = ⌊j²/2⌋ + 1.
* 41 is not a Pell number either, so the §4.2 sentence about 59 and 71 wrongly implies that it is.
* n = 1 is degenerate: no square ever lies in its window. The script's own G ∩ L excludes 1, while its list includes it.

The coincidence set itself is confirmed to 10^7.

**Repair.**
* Write "k ∈ A079496: the odd-index Pell numbers 1, 5, 29, 169, 985 and the even-index half-companion Pell numbers 3, 17, 99, 577".
* Mark n = 1 as vacuous.
* Write "29 is Pell, 41 is half-companion Pell, 59 and 71 are neither".

### F4. B is not a subgroup of M
**Claim.** "[M : B]" is used throughout. §3 says "B is the only subgroup this large". §0 says "C_M(6A) ... contains Fi22".

**Finding.**
* C_M(2A) = 2.B is a non-split double cover.
* Every other maximal subgroup of M is smaller than B. The next is 2^{1+24}.Co1, of order about 1.4·10^26, against about 4.2·10^33 for B.
* So |M|/|B| is a quotient of orders, not a subgroup index. Exactly, |M|/|B| = 2[M : 2.B] = 2·|2A|, with |2A| = 97,239,461,142,009,186,000.
* Fi24' likewise appears as the quotient of the non-split 3.Fi24'. So [M : X] for X = Fi24', Fi22, Co1 should be read as |M|/|X|.
* Th, HN and He are genuine subgroups, as direct factors of C(3C), C(5A) and C(7A).
* From recollection, C_M(6A) ≅ 6.Fi22, so Fi22 is a quotient, not a subgroup.

**Repair.** "|M|/|B| = 2[M : C_M(2A)]"; "2.B is the largest maximal subgroup"; "C_M(6A) involves Fi22".

### F5. Row E's test sentence is false
**Claim.** Row E: "(25, 45) and (45, 75) carry no square-sum event (cycle onset is 32 ...)".

**Finding.** The window (25, 45) contains:
* the Hamiltonian-cycle onset at n = 32: exhaustive search finds cycles exactly for 32 ≤ n ≤ 40, and Gerbicz (cited) gives all n ≥ 32;
* the point n = 31, from which Q_n has no vertex of degree ≤ 1.

**Repair.** "(25, 45) and (45, 75) contain no path onset or loss, but (25, 45) contains the cycle onset at 32. The match is post hoc: NUMEROLOGY."

### F6. The nesting of the lists is automatic
**Claim.** Row C presents the nesting as Fricke-centralizer structure via ATLAS power maps.

**Finding.** By the note's own formula genus(X_0^+(M)) = Σ_{d|M} newplus(d), the genus is monotone under divisibility of M.
* For squarefree N | N' and p ∤ N', genus(X_0^+(N'p)) ≥ genus(X_0^+(Np)).
* Hence L(N') ⊂ L(N) ∖ primes(N').
* This gives L(21) ⊂ L(7) and L(6) ⊂ L(2) ⊂ L(1) with no power maps.
* Checked on 4850 pairs: 0 violations.
* The same monotonicity justifies the script's restriction to Ogg primes.

**Repair.** State the monotonicity lemma as PROVED. Keep STRUCTURAL only for the equalities with sporadic prime sets, which come from Theorem 1.3 plus the ATLAS.

### F7. Row A's growth constant is wrong for N > 1
**Claim.** Row A: the genus of X_0^+(Np) "grows like Np/24".

**Finding.** The genus grows like ψ(N)p/(24·2^{ω(N)}). This equals Np/24 only for N = 1.
* N = 6, p = 106013: genus 13116, against Np/24 = 26503.
* N = 21, p = 121001: genus 39955, against 105876.

**Repair.** Replace Np/24 by ψ(N)p/(24·2^{ω(N)}).

### F8. Row F treats a repeat as a different pattern
The chain 1 | 2 | 10 gives (4, 6, 5), the same multiset {4,5,6} as the owner's chain. So "other chains give other patterns" is wrong for it. The tuple order is also inconsistent: (6,5,4) for the owner's chain, the `.out` order elsewhere. Use one order and note the repeat.

### F9. THM-1370's proof ingredients are misdescribed
* The proof also uses Moon's (n−1)-subtournament theorem.
* The decisive floor values are f(5) = 9 and f(7) = 25.
  * f(8) = 45 is exhaustive in THM-1370.
  * f(9) = 75 is Busch's theorem (cited) and is not used in the proof.
* The auditor's exhaustive check of all 7-vertex tournaments gives f(3..7) = 3, 5, 9, 15, 25. H = 7 and H = 21 are never attained for m ≤ 7, and the maximum H at 7 vertices is 189.
* Two canon files share the number THM-1370 (the h-spectrum one and the elliptic-Jacobian one). This collision predates the note.

### F10. The free action in row D holds for every tournament
H(QR_7) = 189, Aut(QR_7) = F_21, and the action is free with 9 orbits: all confirmed. But the free action is automatic. An automorphism that fixes a directed Hamiltonian path fixes its vertices in order, so |Aut(T)| divides H(T) for every tournament. Checked on all circulant tournaments with 5, 7 and 9 vertices. State this, so that 21 | 189 is not read as a special link.

### F11. Attribution of the group readings
Aricheta's paper names B and Fi24', and HN and He through Nakaya. Nakaya's abstract covers levels 2 and 3 only. The Fi22, HS, M12 and L3(3) readings are the note's own comparisons of the Appendix A lists with ATLAS orders, not statements in Aricheta.

**Repair.** "lists from Aricheta's Appendix A; prime-set equalities checked here; group reading via Theorem 1.3 plus the ATLAS".

### F12. Two loose statements in §1
* **Sylow subgroups.** p² dividing |M| does not by itself make the Sylow subgroup non-cyclic. For p = 11 and 13 the equivalence uses the ATLAS facts that the Sylow 11-subgroup is 11² and the Sylow 13-subgroup is 13^{1+2}.
* **Element orders.** The order-2p proof needs the line "M has exactly two involution classes, 2A and 2B". With it, the proof is complete.

### F13. Missing or weak hostile tests
* The status line promises "each with a hostile test", but rows A and B have none.
* The §3 probe fails every other index only at the exponent of 2. It never tests the post-hoc freedom ("4 read as 13"), which is what makes the cycle numerology.
* The NUMEROLOGY verdict is still right on a priori grounds.

### F14. Reproducibility nits
* Rerunning both scripts gives byte-identical output except for one line: the committed `level_ogg_lists.out` starts with PARI's stack-size warning, which goes to stderr, so the documented command does not reproduce it.
* Runtimes match the claims: 37 s and 127 MB, and 28 s.
* "Brute force to 10^7" (§4.2) is a scan of the defining condition over k, which is adequate.
* Reword "√2-convergent denominators and numerators" to say exactly which convergents are meant.

## Items confirmed

### §2
**The genus formula** for squarefree M is right, with no oldform subtlety.
* For ℓ | M/d, W_ℓ swaps f(τ) and ℓ·f(ℓτ), with eigenvalues ±1.
* For ℓ | d, W_ℓ acts on every f(eτ) by the sign ε_ℓ(f).
* So a newform contributes exactly one W-fixed vector if all its signs are +1, and none otherwise.

Checks of the formula:
* It equals Riemann–Hurwitz at 237 levels up to 1562, with 0 mismatches.
* The CM fixed-point formulas match the Lefschetz identity on 247 pairs, and reproduce the 44 squarefree genus-zero levels N+ of Conway–Norton (from recollection).
* PARI's `mfatkininit` normalising constant was 1 at all 91 levels, so the original trace cross-check is valid.

**Route 1 (Riemann–Hurwitz).**
* All 13 lists match, searching every prime p < 2000.
* Row 2− is {3,5,7,11,23}.
* Aricheta's Table 3 is reproduced.
* The drop-out genera are 1,1,1,1 and 1,1,1,2,3, with 0 at level 2p, as stated.

**Route 2 (direct supersingular points).**
* Explicit Hauptmoduln of X_0^+(N) were built from eta quotients and Eisenstein series. They have integer coefficients and equal T_NA up to a constant.
* X_0(N) was presented as F(s,J) = Π(J − j(Qτ)), with zero residuals to O(q^{100}). This certifies invariance under all the Atkin–Lehner involutions.
* The supersingular j-polynomial came from Deuring's Legendre form, cross-checked by point counting for p ≤ 97.
* For p < 200 every list matches. In every case #non-F_p points = 2·genus(X_0^+(Np)), and every point is F_{p²}-rational.
* X_0(2) gives {3,5,7,11,23}, and Aricheta's Table 2 is reproduced for p ≤ 17.

**Transcription.** Theorems 1.1 (= 3.1) and 1.3 (= 3.3), Remark 1.2 and all thirteen Appendix A rows match the PDF.

### §1
* The auditor's own factorisations of the 11 sporadic orders multiply out to the script's decimals.
* |M|/|B| = 2^5·3^7·5^3·7^4·11·13^2·29·41·59·71.
* 2^5·3^7·5^3·7^4·11·13^2 = 39,046,339,332,000.
* The valuation classes are as stated.
* The prime sets of Fi22, Fi24', HN, He, HS, M12 and L3(3) are correct, and π(Co1) ⊂ π(B).

### §3
The exponent table recomputes exactly.

### §4
* Component counts for n = 1..40 and the n = 12 components are correct.
* s(n), the runs and the boundary law are confirmed for n ≤ 2·10^5.
* The onset lemma's proof is correct.
* Exhaustive counts give the path set {1, 15, 16, 17, 23, 25..30}, with counts equal to OEIS A071983.
* The chains at 15, 16 and 17 are unique, and 18's only partner below 31 is 7.
* The induced edge counts at n = 25 are 5, 3 and 4.
* The hostile-probe numbers are as printed.
* The descriptions of prior work are accurate.

### §5
* THM-200 and THM-079 are the right files.
* The Camion/Busch windows are as stated.
* Row C's lists are correct.
* The 21A power-map recollection is plausible and correctly flagged.

## Independent checks (in this directory)

| Script | Time | What it covers |
|---|---|---|
| `rh_genus_lists.gp` | 17 s | fixed points against Lefschetz; genus-zero levels; all 13 lists for p < 2000; row 2−; Table 3; drop-out genera; monotonicity; growth |
| `direct_rationality.gp` | 2 s | direct supersingular-point test for p < 200; the Theorem 1.1 identity; X_0(2) and Table 2 |
| `formula_vs_rh.gp` | 20 s, 122 MB | the genus formula against Riemann–Hurwitz at 237 levels |
| `audit_group_orders.py` | < 1 s | §1 orders and quotients; \|M\|/\|B\| = 2·\|2A\|; the §3 table; chain sizes |
| `audit_square_sum.py` | 6 s | §4 claims; cycles exactly for 32..40; QR_7 |
| `strong_floor.c` | 5 s | the strong floor f(3..7) by exhaustive enumeration |

All six end with "ALL AUDIT CHECKS PASSED".

## Downloads (kept local; sha256 prefixes)

| file | sha256 prefix |
|---|---|
| Aricheta PDF, arXiv 1809.07421v2 | ce22b774db7c83c0 |
| its pdftotext output | 08522b277b086fe0 |
| Nakaya arXiv metadata (arXiv 1809.02363) | 47fcd67663fbe801 |
| OEIS search for the k-set (found A079496) | e45d9daad651d0f9 |
| OEIS search for 9, 25, 289, 841, 9801, 28561 (no match) | 74234e98afe7498f |

## Corrections applied by the main session

All of F1–F14 are applied; the note's §8 lists them, and the lineage is MISTAKE-668. Script changes:
* `square_sum_structure.py` reports both probe criteria.
* It names the coincidence k-set as A079496 and marks n = 1 as degenerate.
* `level_ogg_lists.gp` checks the monotonicity lemma.
* Both outputs were regenerated with stdout only.
