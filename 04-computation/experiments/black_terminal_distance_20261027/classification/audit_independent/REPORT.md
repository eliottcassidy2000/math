# Independent audit of the black-descent classification (THM-4811 Theorems A, B, C; classification/PROOF.md)

**Verdict: PASS WITH CORRECTIONS** (2026-10-27). Theorems A, B and C are correct. Every step of PROOF.md Parts I-VIII
was checked, and the cited sources (`blue_quotient_20261020b/structure.md` sections 2-3, `black_descent_20261020c/cycles.md`
sections 3-7, the rigidity theorem F6, LEM-003, the prefix-rerouting families) say what the proof uses. One side claim was
overstated (major, status only); seven minor points. The auditor wrote 18 exact edits (`required_edits.json`, each old
string occurring once; built and checked by `make_edits.py`), all applied by the producer. The harness blocked the
auditor's own report file; this is its report, saved by the producer. Nothing in the repository was edited or committed
by the auditor; PROOF.md and classification/ were unchanged since 9245de3e2a when audited.

Checked in detail: Theorem C (the coset step, mu <= 2, the mu = 1 and mu = 2 impossibilities); III.1-III.5 (minimum cuts
never cross, form a chain, and anti-automorphisms reverse it); Theorem D (the stabiliser count); V.1-V.3; the last-block
lemma (the argument leaves exactly m = 6, |B| = 3, and the finite check covers that case and more; (H2) is correctly
derived in both uses); Parts VII and VIII (the descent targets never fall into E or R5 when it matters).

## Major finding (status of a side claim; Theorems A-C unaffected)

THM-4811 (C3) and PROOF.md's Remark said Theorem C "also proves the open-path form of the extremal-predecessor statement"
of `cycles.md` section 7, and the OPEN lists said only "the case of two minimum-score vertices" was open. Both overstated
what is proved. The case analysis never uses lambda = delta; what it proves is C': every nonregular strong self-converse
tournament has a *black* Hamiltonian path starting at a minimum-score vertex. Such a path can still end at a
maximum-score vertex, and Steps 1 and 3 use blueness, not just endpoint scores, so the extremal-predecessor statement
does not follow in general. The open-path form is proved only in case (c1) and when there is a single minimum-score
vertex (mu = 1); the closed form only in case (c1) and when mu = 1 and delta = 1. Both forms stay OPEN; C' and both forms
hold for every nonregular strong self-converse class through n = 10. Fixed by edits 1-5.

## Minor findings

2. Benhocine-Wojda is a CITED external theorem (the repository has only the publisher abstract), not a "PROVED repo
   input", and was missing from the CITED list (edits 6-8).
3. The order-10 and order-11 censuses test the descent property directly only on the nonregular self-converse classes
   (8,599 and 276,971); the rest follow from Part I. "Census of all 884,464,590 classes finds one exception" should say so
   (edits 9-11).
4. "Stays in that regime": from lambda = 2 the target has lambda = 1, so the preserved regime is lambda < delta (edit 12).
5. The title's black-distance clause needs n >= 4 (C3 is strong and has no black move) (edit 13).
6. Section C2 did not say what P_t is (a Hamiltonian path of sigma(S) from sigma(x0)) (edit 14).
7. The status line still said "not yet independently audited" / "under a separate audit" (edits 15-16).
8. Optional: the n = 15 lambda = 2 search is in fact exhaustive, and the new searches can go in the Section 13 table
   (edits 17-18).

## Reproduction (own C + nauty code; nothing imported from tour.hpp)

| Check | Result |
|---|---|
| Theorems A and C, all strong classes n = 3..9 | Only C3 and R5 lack a black descent; Theorem C holds for every lambda = delta class except R5. lambda histograms match THM-4811; no crossing minimum cuts; F10 and F6 (H = \|Aut\| only for C3) hold. |
| n = 10, all 9,355,949 classes read | 8,599 nonregular self-converse (nauty and own backtracking agree). Only E10 lacks a black descent (657 paths, 9 candidates, all blue). Theorem C holds on all 8,595 lambda = delta ones. 16 classes have lambda < delta. |
| Family E built from the definition | All 112 members to n = 18 pass every Part VII/VIII claim, including the two-move certificate. Distinct members 1, 1, 3, 3, 8, 8, 20, 20, 48 (n = 10..18). The n = 10 and n = 11 census exceptions are E_1(S5) and E_2(S5). |
| Last-block lemma, (H1) tested directly | m <= 8 all pairs, m = 9 non-rigid classes: only C3 with A empty. All two-orbit configurations at p = 3, 5, 7 fail (H1). |
| Lemma BL and bottom-layer census, \|S\| <= 9 | Lemma BL holds; Condition A + transitive gives 1, 2, 5, 12, 28, all bottom layers. |
| Exhaustive searches for exceptions | F/G shapes for every strong S with \|S\| <= 8 (n <= 17); every (\|L_0\|, middles) case of lambda = 1 < delta for n <= 17; lambda = 2 < delta for n = 14..17. No exception outside E. |
| Perturbations and random samples | 201,665 sigma-symmetric flips of E members (n <= 16) and 119,525 random self-converse tournaments (n = 11..16): no exception outside E; Theorem C always holds. |
| Exact black-distance BFS on classes, n <= 9 (no use of F1) | d_black = lambda for every strong class except R5. F1 identity checked on 600M+ (move, cut) pairs, no failure. |
| Section 1 facts about E10 | All confirmed. |
| Section 12 strong form at n = 10 | Exactly 2 failures, as claimed. |

Not re-run: the order-11 census (instructed) and the order-11 strong-form counts.
