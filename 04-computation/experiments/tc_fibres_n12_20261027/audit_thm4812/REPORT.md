# Independent audit of THM-4812 (the {7, 21} presentation gaps at every order)

**Verdict: PASS WITH CORRECTIONS** (2026-10-27). THM-4812 is true as stated. The proof chain (a) → (b) → (c) → (d) holds,
and every finite certificate was reproduced exactly with separate code. No tournament with tc in {7, 21} appeared
anywhere. Two proof-text passages were wrong and have been rewritten (M1, M2); neither changes a conclusion.
`required_edits.json` lists the 14 edits, all applied by the producer (with the AUDITED line naming this report and
noting that the (g) extension to v = 23, filed after the audit, is not covered).

Audited: (a)-(f) fully, (g) for v <= 21 (commit 019882e00c). Not audited: the v = 23 extension (commit 44ebc341cc); it is
consistent with (c), since every vertex-transitive tournament other than K1, C3, R5, P7 has tc >= 25. The auditor made
no edits to the repository (the harness blocked the auditor's own report file; this text is its report, saved by the
producer).

## Findings

**M1 (major, text). (c), case m > 401: "with a margin increasing in m" is false.** The margin is
F(m) = k log(0.0458 m^3) - log k! - log 23 - ((m-1)/2) log 3 with k = floor(m/20). For fixed k its slope is
3k/m - (1/2) log 3 < 0, so it decreases; it decreases at 720 of the odd m in [401, 2000), and its minimum over odd m in
(401, 200000] is 25.52 at m = 419 (33 is the value at m = 401). The conclusion survives: the minimum on each block of
fixed k is at m = 20k + 19, where Robbins' bound for k! gives a margin of at least k(2 log k - 4.083) - (1/2) log(2 pi k)
- 13.03, which is >= 22.7 at k = 20 and increasing after (greedy_vt_bound.py). The other two estimates are correct:
c3 - k(D+1) >= 0.0229 m^3 (checked exactly to m = 20,000) and 2^k i_k >= (0.0458 m^3)^k/k!.

**M2 (major, text). (d), index step.** "Odd and at most 6" leaves 1, 3 and 5. The index divides |P|, and every
configuration with r >= 1 has 21|P|/f < 2.6, so 5 cannot occur. Stronger (S1): (b) applied to G itself gives
tc(S) >= [P : G] 3^r 9^p, which forces [P : G] = 1 in all 24 configurations; the 59 index-3 cases are redundant.

**Minor.**
- m1: the FINITE-EXACT greedy certificate for odd m in [19, 401] had no script; recomputed exactly (greedy.py, smallest
  log-margin 0.962, ratio 2.6, at m = 19), now committed as greedy_vt_bound.py.
- m2: the greedy product bound needs c3 > (k-1)(D+1).
- m3: at order 15, K != 1 because 15 divides |G|, and K is normal, so its actions on the blocks are conjugate.
- m4: the relabelling step deserves one sentence (relabel each orbit by an isomorphism onto the standard copy; this
  conjugates Aut(S)|_O into Aut(X_O), so Aut(S) becomes a subgroup of the fixed P, and every subgroup of that P is
  listed, not only conjugacy representatives).
- m5: "fails for every n >= 23" was checked only to 199; f(n) >= 5^(-1/3) 5^(n/3) and |P| <= 1.618 * 3^(n/3) give
  f(n) > 21|P| for n >= 24 and f(n) > 7|P| for n >= 15 (exact check to n = 3000 agrees).
- m6: THM-4744 belongs in depends_on.
- m7: (e) is a routine corollary of the classical step s(qr) = 1 + q + r; worded so.
- m8: the session note's "24.9 million" invariant tournaments is 21,516,352 (18,138,240 strong).
- m9: bookkeeping (the OPEN audit item; HYP-9352 (2) now a theorem).

## Reproduction (independent code, nothing shared with the repository)

**Methods.** Sigma placed on the last points with each cycle reversed (a different conjugate); plain binary enumeration
instead of a Gray code; strong connectivity by Landau's score criterion cross-checked by BFS; Held-Karp in pull form; an
own orbit-stabiliser |Aut| counter (no nauty), run on every strong tournament with no candidate filter.

**(a) Order 12** (scan12/): all seven types agree with the repository in every column.

| Type | Invariant | Strong | Min H | Candidates (repository's filter) | Smallest tc |
|---|---|---|---|---|---|
| (11) | 2^6 | 0 | - | 0 | - |
| (7) | 2^18 | 178,560 | 10,815 | 30,880 | 393 |
| (5)(5) | 2^14 | 9,024 | 14,625 | 2,320 | 585 |
| (5)(3) | 2^18 | 178,560 | 3,645 | 38,400 | 243 |
| (5)(3)(3) | 2^12 | 2,112 | 12,825 | 192 | 285 |
| (3)^4 | 2^22 | 4,126,848 | 20,889 | 742,896 | 2,215 |
| (3)^3 | 2^24 | 13,643,136 | 2,295 | 2,676,384 | 85 |

The candidate |Aut| histograms are identical (e.g. (3)^3: 3: 2,104,704; 9: 515,808; 27: 45,504; 63: 10,368); the
THM-4790 controls (1,088; 3,760; 23,520; 480; 392,000) are identical; |Aut| divides H in every case; no tc below 64
occurs at order 12. The repository's two (3)^3 slices agree because swapping two fixed points commutes with sigma and
reverses exactly one stored orbit, exchanging the parity slices while preserving H and |Aut|.

**(c) Vertex-transitive lemma** (vt/): regular tournaments of orders 5, 7, 9 number 1, 3, 15; the vertex-transitive ones
have tc 3; 9 and 25; 363, 373, 39 (C3[C3], |Aut| = 81). Circulant minima: 1729 (order 11), 94,949 (13, |Aut| = 39),
162,165 (15, R5[C3], |Aut| = 1215), 792,878,221 (17); Z_5 x Z_3 equals Z_15; every Cayley tournament on Z_3^2 is
C3[C3]. G_A (order 405) and G_B (order 75) are transitive without an element of order 15; each fixes exactly 8
tournaments, all with module blocks (tc 162,165 and 515,241). F_3^5 under Z_5 and F_5^3 under Z_3 have exactly four
invariant subspaces each. Citations checked: Dixon (Canad. Math. Bull. 10 (1967) 503-505; via arXiv:1806.11241), Busch
(EJC 13 (2006) N3), Montgomery-Vaughan (Acta Arith. 27 (1975)).

**(d) Window** (window/): own group construction (element lists, index-3 kernels by evaluating every homomorphism to Z_3
on every element), different standard copies (C3 reversed, P7 on the non-residues): 24 configurations, 83 cases
(73 module, 10 non-module), identical to the repository's job by job in strong count, minimum H, candidates and hits;
185,642 strong, no hit, smallest tc 165; |G| divides |Aut| and |Aut| divides H throughout. The subset-run substitution
DP agrees with Held-Karp on 93,646 tournaments of order <= 16 and with 17 spot checks at orders 17-22 (Held-Karp to 20,
inclusion-exclusion at 22).

**(g), v <= 21:** case counts per v (2, 14, 33, 35, 82, 123, 139, 139) and the order bounds (12, 16, 19, 22) identical;
the 569 jobs have 66 distinct contents (37 are (d) jobs; the 29 new ones at orders 11-12 all have smallest tc >= 63);
the strong census of orders 3-10 (9,355,949 strong classes at order 10) identical in every tc <= 61 multiplicity; fibre
sizes 2, 1, 2, 0, 5, 4, 1, 6, 6, 4, 0.

**Break-it** (breakit/): 400 random odd-order groups and 1,600 invariant tournaments at orders 13-18 (1,251 strong): no tc
in {7, 21}, none below 64, no violation of the product inequality, every orbit regular. 51,115 substitutions with
K1/C3/R5/P7 blocks (every strong quotient on 3-7 vertices; P7 and C7(1,2,3) quotients to order 25; C3[C3] quotients to
order 27; C3[C3[C3]] with 128-bit arithmetic: H = 417,382,500,592,116,023,859, |Aut| = 3^13,
tc = 261,792,936,934,433): no hit, smallest tc 125 (order 13), 37,546 Held-Karp cross-checks agree.

**Status typing** is honest except the three points corrected above (the script-less greedy certificate, the false
monotonicity sentence inside a PROVED step, the OPEN audit item).
