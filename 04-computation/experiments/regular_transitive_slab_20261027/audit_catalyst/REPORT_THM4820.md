# Independent audit of THM-4820 (Eisenstein merge debt and the catalyst threshold)

*Independent auditor (subagent of mac-mini-2026-10-27-slab, own code only), 2026-10-27. The auditor could not write
report files from its sandbox, so the producer saved this report verbatim from the auditor's final message. The
scripts, outputs and certificates in this directory are the auditor's own. `./run_audit.sh` regenerates everything in
about 50 s.*

Every number in THM-4820 holds up, including the headline 41 exactly when k ≥ 36. What is wrong is two status labels,
one overclaim and one unproved asymptotic claim. The audit also found that HYP-9426 (3), the claim that
g ≤ ⌊Q/3⌋ + 1 everywhere, is false: g(101) = ⌊101²/3⌋ + 2.

| Claim | Verdict | Basis |
|---|---|---|
| (A1) β = 3u·u′ − \|u\|\|u′\| | PASS | 35,625 pairs, checked against two different ways of computing values |
| (A2) neutral triples free; (A3) f(n,0,0) = n(n−1) | PASS | proof checked; holds on every count vector of three DP boxes |
| (A4) g(n) = min_k f(n+k,k,k) | PASS | proof checked; brute force on small trees; lattice seeded with 3 or 6 units gives identical tables |
| (B) Φ = ⌊Q/3⌋ bound, telescoping, tight merges | PASS | exhaustive on 707,281 pairs |
| (C) norm-form criteria and g(n) ≥ (n²+2)/3 for 11, 23, 25, 29, 41, 47 | PASS | brute force for every odd n ≤ 241 with 3 ∤ n |
| f(46,35,35) = 42 and f(47,36,36) = 41 | PASS | own exact DP, plus a second pass checking the DP recurrence on every pair |
| Certificate `tree_47_36_36.txt` | PASS | own verifier |
| Row f(11+k,k,k) and thresholds k*(n), n ≤ 14 | PASS | own DP up to k = 62; all 14 thresholds agree |
| (D) lattice values and exceptional lists | PASS on numbers | all 9061 points agree |
| (D) calling 15 exceptions "window-dependent" | FAIL (wrong label) | they are exact |
| Reading: "exceptional status is an inertness statement" | FAIL (overclaim) | true only for the (C) class |
| (E) "asymptotically saves two thirds", "order n² triples" | GAP | not proved |
| HYP-9426 (3), downstream | REFUTED | g(101) = 3402 |

## (A)–(C)

- **(A1):** checked on random and exhaustive small pairs, computing values in integer form and in complex form.
- **(A2):** f(u+(1,1,1)) ≤ f(u) on every vector of the boxes (46,35,35), (47,36,36) and (76,62,62).
- **(A3):** f(n,0,0) = n(n−1) for n ≤ 76 in the DP, and on 780 random all-ones trees.
- **(A4) reduction:** replacing negative units by free pairs never raised the optimum, over all 923 six-unit multisets
  with ≤ 6 leaves.
- **DP recurrence:** matches explicit enumeration of all 403,149 trees for multisets with ≤ 7 leaves, and a separate
  labelled-leaf DP for 8–12 leaves.
- **(B):** the floor identity, the bound and the tight-merge characterisation hold exhaustively for |a|, |b| ≤ 14.
  Cost = Φ(root) + Σ slack holds on 3000 random trees.
  - One wording slip: "a lossless paid merge at 60°" is true only for two single units. In the certificate tree, 28 of
    the 39 β = 1 merges are at 60°; the other 11 are at 79°–89°.
- **(C):** the translation identity, both splitting criteria, the odd-exponent form and M ≡ 1 (mod 8) all agree with
  direct enumeration. All six factorisations of M are correct (41; 3·59; 11·19; 281; 3·11·17; 11·67), and none of these M
  is b² + 3t².

## Finite computations

- **Method.** The auditor's DP differs from the producer's on purpose: β comes from Eisenstein values rather than the
  count formula, states are processed in a different order, and all ordered splits are scanned.
- **Certification.** A second pass checked f(u₁+u₂) ≤ β₊ + f(u₁) + f(u₂) for every pair, and that every state attains
  its value. This certifies f exactly:
  - box (46,35,35): 500,209,345 pairs, 0 violations;
  - box (47,36,36): 581,058,361 pairs, 0 violations.
- **Result.** f(46,35,35) = 42 and f(47,36,36) = 41, so the threshold is sharp.
- **Box (76,62,62).** The row matches the theorem for k ≤ 36 and stays at 41 for every 36 ≤ k ≤ 62. The thresholds are
  k*(1..14) = 0, 1, 2, 4, 7, 10, 14, 19, 24, 31, 36, 44, 52, 61.
- **Certificate.** Cost 41 = Φ(11) + 1, with 118 merges, 119 leaves and root value 11. Its only positive-slack merge is
  at the root, (17,10,16) + (30,26,20) with β = 2. The auditor's DP's own optimal tree is the same tree; that is
  expected, since both programs break ties the same way.

## (D) lattice table

- **Agreement.** The auditor used a Dijkstra-style algorithm, different from the producer's repeated relaxation.
  - At W = 2500 all 9061 points agree, including all 1459 of norm ≤ 400.
  - W = 400, 900 and 2500 agree for n ≤ 16, and W = 7500 agrees for n ≤ 86.
- **Exceptional norms.** The twenty exceptional norms are exactly as listed, and every point at each of those norms is
  exceptional. On norms ≤ 2000 the values are only ⌊Q/3⌋ (7069 points) or +1 (198 points).
- **A count.** "328 norms" counts the norms ≡ 1 (mod 3) in (1, 2000]; including norm 1 gives 329.
- **Rational exceptions.** The rational exception lists are confirmed, and 2601 = 51² is the only norm ≡ 0 exception
  below 7500.

**FAIL: "Their +1 is therefore window-dependent."** A one-line consequence of (B) makes the window irrelevant:
- Every subtree of a cost-C tree costs at most C, so every node has Q ≤ 3C + 1.
- Hence the windowed value g_W(x) equals g(x) whenever W ≥ 3g_W(x) + 1. At W = 2500 that covers every point of norm
  < 2500; at W = 7500 it covers every n ≤ 86.
- Sharper: a zero-slack tree never has a node of norm above Q(x) + [3 | Q(x)], and a tree of total slack ≤ 1 never goes
  above Q(x) + 4.
- The set of zero-slack-buildable points was computed directly, without any window, in both C and numpy. It matches the
  windowed exceptional set on 9061/9061, 27,175/27,175 and 38,083/38,083 points.

So all 20 exceptional norms and all 16 rational exceptions n ≤ 86 are exact. The theorem is also inconsistent with
itself here: it calls 22² and 44² window-dependent, then states g(22) and g(44) exactly.

**51** is not exceptional "by propagation": it has no tight split at all. For norm ≡ 0, tight means β ∈ {−2, …, 1},
and none of the four norm equations for 51 has a nontrivial solution (2601 = 3²·17², 2599 = 23·113, 2603 = 19·137,
2605 = 5·521). **59** is exceptional by propagation, as stated.

## Labels and overclaims

1. **(D).** Drop "window-dependent" and "inside the stated window"; the table gives exact values of g.
2. **Reading.** "Exceptional status is an inertness statement" is false in general. 229 is a split prime, 59 has tight
   splits, and 51 = 3·17. Restrict the sentence to the (C) class.
3. **(E).** "Saves two thirds" and "order n² neutral triples" need an upper bound g(n) ≤ n²/3 + o(n²), which is not
   proved. "Order n²" also needs k* = O(n²), which is known only for n ≤ 14. Mark both as conjecture.
4. **"(C) excludes 53, 55, 71".** (C) is only a sufficient condition and excludes nothing; 59 fails it in the same way
   and is exceptional. Better: "(C) does not apply to 53, 55, 59, 71; computation shows 53, 55, 71 are not exceptional
   and 59 is."
5. **Status block.** "f(11+k,k,k) = 41 for all k ≥ 36" could be labelled PROVED, since it follows from (C), (A2) and
   the verified certificate.
6. **Cross-link.** `05-knowledge/results/session75_final_mathematics_20261027/EISENSTEIN_INVOICES.md`, committed 4
   minutes after the theorem, independently reaches the same 41 and k ≥ 36 under the same reading.

## Downstream (HYP-9426)

- **(3) is refuted: g(101) = 3402 = ⌊101²/3⌋ + 2.**
  - **Lower bound, by arithmetic.** A tree costing ≤ Φ + 1 would need a first split with β ∈ {−4, …, 2}. The
    corresponding norm equations for 101 all fail:
    - 10201 has only trivial solutions;
    - 10203 = 3·19·179, 10197 = 3²·11·103, 10207 = 59·173 and 10209 = 3·41·83 each have a prime ≡ 2 (mod 3) to an odd
      power;
    - the remaining two are ≡ 2 (mod 3).

    Brute force agrees: 0 tight and 0 slack-1 splits.
  - **Upper bound.** An explicit tree of cost 3402 (357 distinct values, 8566 leaves), independently verified
    (`audit_value_tree_101.txt`).
  - **Extent.** Exact at W = 10,500: g ∈ {⌊Q/3⌋, ⌊Q/3⌋ + 1} on every norm below 10,201. The same arithmetic gives excess
    ≥ 2 for 78 odd n ≤ 3000, starting 101, 115, 197, 235, 257, 289.
- **(4).** "51 and 59 are the two non-window-free odd exceptions" is wrong on both counts. The predictions 92, 94 and
  102 hold, as do 88 and 100, and 99 is a new odd exception.
- **(5).** "The lower constant 2/9 is PROVED" inherits the gap in (E).

Files in this directory:
- **Scripts:** `audit_algebra.py`, `audit_small_trees.py`, `audit_dp.c`, `audit_tree.py`, `audit_lattice.c`,
  `audit_compare.py`, `audit_excess2.py`, `audit_slack_levels.py`, `audit_value_tree.py`, `run_audit.sh`.
- **Outputs:** the `*.out` files and `audit_lat_W2500.txt`.
- **Certificates:** `audit_own_tree_47_36_36.txt`, `audit_own_tree_47_36_36_lastmin.txt`, `audit_value_tree_101.txt`.
- **Ignored:** `.gitignore` covers the binaries and the large lattice dumps, which `run_audit.sh` rebuilds in seconds.
