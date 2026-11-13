# Independent audit: kind-pasteur-2026-11-13-S30 (THM-5155, HYP-9681, observer lattice)

Independent adversarial audit, 2026-11-13. Every repository file was read-only except this report. Nothing was
stashed, committed or pushed. **Chan–Pak remains OPEN.** Nothing in the packet claims otherwise.

## Verdict: PASS WITH CORRECTIONS

| object | verdict |
|---|---|
| THM-5155 (A) | **PASS.** Correct; one wording nit (m6). |
| THM-5155 (B) | **PASS.** Reproduced exactly with independent code. There is no collision through order 8, and order 9 has exactly the six stated pairs. Every structural claim holds. The same pattern holds at order 10 (new data below). |
| THM-5155 (C) and the README table | **PASS WITH CORRECTIONS.** Every number reproduces. But "scores" in the table and in "(scores, Q)" is the per-vertex fold. From order 7 that is strictly coarser than the score sequence up to converse (m1). |
| HYP-9681 | **FAIL AS FILED.** "Strong module" is false at order 13 (M1). "A converse flip", meaning one flip, is false at order 15 even for arbitrary modules (M2). Its evidence sentences on Session 320's order-24 pair and on the S29 order-12 block reversal are correct, but the latter is misleading (M1). It must be restated before merge. |
| README exposition and scope | **PASS WITH CORRECTIONS.** The "This is forced" paragraph argues about slides incorrectly (m2). The H1 and an index line lack the order bound (m3). |
| Credits | **MINOR corrections** (m4). No prior joint (score, arc-multiset) census was found. |
| Scripts and replay | **PASS WITH CORRECTIONS.** Both scripts replay byte-for-byte to the README table and the six codes. The structural check is not shipped, and the inputs are not pinned (m5). |

**In one paragraph.** Every input of the packet was recomputed from scratch. The class lists were checked with an
own canonical form (individualisation-refinement), H and arc weights by plain DFS over Hamiltonian paths, Q from the
descent distribution, and C, c3, |Aut| and strongness directly. All agree with the author's files on all 198,946
classes of orders 3-9. The observer table and (B) reproduce exactly. The structure of the six order-9 pairs is as
stated, and they are exactly the 6 of the 13 converse pairs of non-self-converse strong 6-tournaments whose score
sequence is self-reversed.

The conjecture fails as filed. The S29 audit's interior block reversal keeps the score sequence once its two cores
have equal order and one score sequence, as in TT2 ⇒ A ⇒ B ⇒ TT1 against TT2 ⇒ B ⇒ A ⇒ TT1 for two strong
5-tournaments with score sequence (1,2,2,2,3). This flips a non-strong module and no strong one. A second witness,
two independent core flips at order 15, is not a single module flip at all.

The full non-strong order-10 census, combined with the S29 audit's order-10 strong-sector result, gives exactly 62
coincidences at order 10. All are self-reversed core flips.

---

## MAJOR findings

### M1. HYP-9681 is false as titled: an order-13 coincidence that is no strong-module flip

**Exact sentences.**
- HYP-9681 title: "Every coincidence of the joint observer (score-sequence pair, arc-usage multiset) between distinct
  merged tournament classes is a converse flip of a strong module whose internal score sequence is self-reversed and
  whose seam weights the flip preserves".
- HYP-9681 status: "The order-12 interior block reversal of the S29 audit (C3 => v => C3 => S4 => w against
  C3 => v => S4 => C3 => w) changes the score sequence, so it is not a coincidence of this observer."

**The problem.** The S29 audit (M2) proved that reversing an interior core-block [A ⇒ B] preserves the arc-usage
multiset. Its order-12 instance changes the scores, as the HYP says (recomputed: (0,2,2,3,3,6,6,6,8,10,10,10)
against (0,2,2,2,5,5,6,6,8,10,10,10); the score pairs differ too). That happens only because |C3| ≠ |S4|. Take two
non-isomorphic strong tournaments of one order with one score sequence, and the block reversal keeps the scores as
well.

**Witness (order 13).** Let A and B be the strong 5-tournaments with codes 209 and 305 in `cls_5.txt`. Both have
score sequence (1,2,2,2,3) and are self-converse and rigid, with H = 11 and 13, and they are not isomorphic.
- T = TT2 ⇒ A ⇒ B ⇒ TT1 has out-masks `8190 8188 8072 8144 8100 8076 8116 4352 5632 6272 4736 5504 0`.
- T′ = TT2 ⇒ B^op ⇒ A^op ⇒ TT1 (≅ TT2 ⇒ B ⇒ A ⇒ TT1) has out-masks
  `8190 8188 8176 8132 8104 8136 8080 7680 5248 6400 6656 4352 0`.
- Both have H = 143, scores (0,2,3,3,3,4,7,8,8,8,9,11,12) and 22 zeros. Their arc-usage multisets are identical
  (positive values 2, 3, …, 143). This was checked by the auditor's DFS (`obs.c`) and independently by a
  forward/backward subset DP with the split sum (`crosscheck_dp.py`).
- They are distinct merged classes. The canonical forms differ. Structurally, the component sequences are (1,1,A,B,1)
  and (1,1,B,A,1), and T^op has (1,B,A,1,1).
- **No strong-module flip gives it.** All 13 nontrivial modules of T, and all 13 of T^op, were enumerated by brute
  force over subsets, flipped and canonicalised.
  - Only two flips reach the merged class of T′: A ∪ B (10 vertices) and {v₁} ∪ A ∪ B ∪ {w} (12 vertices). Neither
    is strong.
  - In general, a strong module lies inside A or inside B. Every proper module of a strong 5-tournament is
    self-converse, so a strong-module flip of T is isomorphic to T.
- The other three equal-score pairs give coincidences of the same kind: (209, 152) and (152, 305) with (1,2,2,2,3),
  and (516, 642) with (1,1,2,3,3).

**Why 13 is the least order for this mechanism.**
- Flipping a module X keeps the scores iff X's internal score multiset is self-reversed, by the argument of THM-5155
  (A). Then X and X^op have the same score multiset. By Landau's equality points (classical), they then have the same
  ordered component sizes.
- For X = A ⇒ B this forces |A| = |B| = m and scores(B) = m − 1 − scores(A).
- A new class needs B ≇ A^op. So there must be two non-isomorphic strong m-tournaments with one score sequence, which
  needs m ≥ 5, since the strong 3- and 4-tournaments are unique.
- The block must be run-flanked (S29 M2), and the runs must differ, because TT_a ⇒ X ⇒ TT_a with X flipped is the
  global converse.
- Hence the order is at least 2·5 + 3 = 13. Orders 11-12 remain open for other mechanisms.

**Fix.**
- Mark the filed statement REFUTED at order 13 (this audit), or restate it together with M2.
- Replace the evidence sentence with: "The S29 order-12 block reversal changes the scores; block reversals with two
  equal-order cores of one score sequence keep them (order 13, audit M1)."
- Say "module (strong or not)".
- "Whose seam weights the flip preserves" is undefined for modules that are not run-flanked cores. Replace it with
  "whose flip preserves the arc-usage multiset".

### M2. One flip is not enough: an order-15 coincidence needs two

**Exact sentences.**
- HYP-9681 title: "… is a converse flip of a strong module …".
- HYP-9681 body: "The conjecture is that this holds at every order, with "module" in place of "run-flanked core",
  once the seam weights are preserved."
- README: "HYP-9681 conjectures that this is the only failure mode at every order: reversing a module whose own score
  sequence is self-reversed and whose seam weights are preserved."
- hypotheses/INDEX.md: "[HYP-9681] (OPEN): this is the only failure mode at every order."

**Witness (order 15).** Let S1 and S2 be the 6-vertex cores of the THM-5155 classes 67436800 and 203522304, with
out-masks `2 40 3 21 7 29` and `2 12 9 49 7 23`. Both have score sequence (1,2,2,3,3,4), are rigid and are not
self-converse, and they come from different converse pairs.
- T = v ⇒ S1 ⇒ v ⇒ S2 ⇒ v and T′ = v ⇒ S1^op ⇒ v ⇒ S2^op ⇒ v. Both have H = 575 and the same score sequence, and
  their arc-usage multisets are identical by both methods. They are distinct merged classes.
- Flipping each of the 16 nontrivial modules of T, and of T^op, never reaches T′. It takes two strong-core flips, and
  the intermediate v ⇒ S1^op ⇒ v ⇒ S2 ⇒ v is itself a coincidence of the observer.
- The same holds with S2 = S1^op. Then T and T′ are both self-converse; each has 19 nontrivial modules, and no
  single flip reaches the other class.

Hence even with "module" in place of "strong module", "is a converse flip" must become "is generated by converse
flips". Through order 10 the two readings coincide.

**Fix.** A restatement consistent with every known example: "Every coincidence of the joint observer between
distinct merged classes is connected by a chain of converse flips of modules (strong or not), each of which keeps
the score sequence (the module's internal score multiset is self-reversed) and the arc-usage multiset." Its evidence:
- orders ≤ 9 (THM-5155) and order 10 (this audit, below);
- the order-13 and order-15 witnesses above;
- Session 319's orders 22 and 42 (U ⇒ TT1 against U^op ⇒ TT1; U is regular, so the scores agree);
- Session 320's order 24.

Update hypotheses/INDEX.md and the README's HYP paragraph to match.

---

## MINOR findings

### m1. The "scores" column is the per-vertex fold, not the score sequence up to converse

**Exact sentences.**
- README table headers: "scores (symmetrised)" and "scores + Q".
- README: "'Symmetrised' scores fold s with n − 1 − s. The exact converse-invariant observer, the pair {scores(T),
  scores(T^op)}, gives the same last column."
- THM-5155 (C): "184 symmetrised score classes … 63686 of (scores, Q)".

**The problem.** The last-column claim is right: 97,138 either way, and 188,584 either way on the order-10
non-strong classes. The other two columns are not the same:

| n | 5 | 6 | 7 | 8 | 9 |
|---|---:|---:|---:|---:|---:|
| scores, per-vertex fold (packet) | 7 | 14 | 33 | 76 | 184 |
| scores up to converse (pair) | 7 | 14 | **37** | **93** | **269** |
| fold + Q (packet) | 9 | 27 | 182 | 2,093 | 63,686 |
| pair + Q | 9 | 27 | **184** | **2,178** | **66,868** |

- The pair counts equal (A000571 + #self-reversed)/2: (59 + 15)/2, (167 + 19)/2 and (490 + 48)/2. The score-sequence
  counts 2, 4, 9, 22, 59, 167, 490 were reproduced.
- Example: at order 7, (0,1,2,4,4,5,5) and (0,2,2,2,5,5,5) fold to the same multiset but are not converse.

**Terminology clash.** Elsewhere in the repository a "folded score profile" is the class of the sorted multiset under
global reflection, which is the pair. See merged_metagraph_endpoint_transport_20261020 §2 and `folded_score` in
incidence_closure_20261113/closure/check.py. The packet's per-vertex fold is a different, coarser observer.

**Fix.**
- Label the column "per-vertex folded scores {min(s, n−1−s)}" and add the two pair rows, or switch to the pair.
- In (C), write "63686 of (per-vertex folded scores, Q); 66868 with the score pair".

### m2. The "This is forced" paragraph

**Exact sentences** (README): "This is forced. Slides shift the core's scores by the number of singletons after it,
so the score sequence sees them. A flip maps the core's internal scores s to |S| − 1 − s, which is invisible exactly
when that multiset is self-reversed."

**The problems.**
- It reads as a general argument. In fact it holds at order 9 only through THM-5154 (B)'s finite classification.
- The slide sentence is false for the converse-invariant observer. Each of the six residual pairs is itself a slide,
  since TT2 ⇒ S^op ⇒ TT1 = (TT1 ⇒ S ⇒ TT2)^op (S29 audit m1), and the score pair does not see it.

**Fix.** "Given THM-5154 (B), this is forced at order 9. The score sequence determines the ordered sizes of the strong
components (Landau's equality points; classical). So the score pair separates two move-related classes unless their
component-size sequences agree up to reversal. At order 9 that leaves only TT2 ⇒ S ⇒ TT1 against TT1 ⇒ S ⇒ TT2 =
(TT2 ⇒ S^op ⇒ TT1)^op. Up to converse this is a core flip, invisible by (A) exactly when S's score multiset is
self-reversed."

This was checked on the 22 order-9 multiset groups:
- the 13 groups of non-self-converse 6-cores are the only ones whose size sequences agree up to reversal;
- the 8 slide groups with self-converse cores and the two-core slide are all separated by component sizes alone.

### m3. Scope

- The THM-5155 H1 reads "score sequence and arc-usage multiset pin the merged classes except symmetric-score core
  flips", and the filename matches. It carries no order bound, and it is false in general. Neither M1 (order 13, a
  block reversal) nor Session 320's order-24 pair (a module inside a strong tournament) is a core flip. Session 319's
  order-22 pair is a core flip, but of an end core, which is not run-flanked. Fix: append "through order nine", as
  THM-5154's H1 does.
- hypotheses/INDEX.md (quoted in M2): fix it as in M2. The results/INDEX.md entry is correctly scoped.

### m4. Credits

- **H column.** The values 7, 19, 77, 320, 1,520 are |S(n)| of THM-4889 (klein-2026-11-06).
- **Other inherited columns.** The merged counts and the Q column are THM-5153 (S28), and the arc-multiset column is
  THM-5154 (S29). Mark these KNOWN in (C). The new columns are the two score columns, (H, C, |Aut|, c3) and the two
  joint columns.
- **Framing.** THM-1965 (kind-pasteur-2026-07-21-S128c142) already maps the tournament-invariant lattice under
  refinement for n ≤ 6. It covers score, cycle vector, H, |Aut| and spectra, and finds that {score, specA} pins the
  class for n ≤ 5. Cite it as the precursor of the observer lattice.
- **Sibling joint observer.** incidence_closure_20261113 closure/PROOF.md §3 has the "augmented Q" (Q, folded scores,
  Aut, SC, d), with 9, 29, 204 cells at n = 5, 6, 7. Cite it as related.
- **Landau / Moon.** Cite them if the m2 argument is adopted. Session 319's order-22 and order-42 pairs are further
  evidence for the HYP.
- **Novelty search.** The canon and results were grepped for score-sequence and arc-usage, arc-weight or arc-multiset
  combinations, "observer lattice", and "symmetrised" or "folded" scores. No prior joint (score, arc-multiset) census
  was found. THM-1370 (the H spectrum omits 7 and 21) is unrelated to observer completeness.
- **What is new.** (B) and its residue analysis are new. The score half of (A) is new but elementary, and the arc half
  is THM-5154 (A)(ii).

### m5. Replay and scripts

- **Replay works.** `observers.py` and `score_arcs.py` were run on the author's inputs. They reproduce every number in
  the README table and the six codes (auditor directory `replay/`).
- **Unshipped structural check.** The THM status claims "runs (2, 1); cores converse but not isomorphic; core score
  sequences …", and (B) says "rigid". No shipped script produces these: `score_arcs.py` prints only the codes, the SC
  flag and |Aut|. Ship a check; the auditor's `structure9.py` is one.
- **Unpinned inputs.**
  - Say where classgen and fprof live: first_creators_genealogy classgen, and kps_s26 scripts/README.
  - Point to kps_s29 README steps 1-2 for arcms.
  - Pin the line counts 2, 4, 12, 56, 456, 6,880, 191,536.
- **Dead code.** `observers.py` computes `info[r]['S']` under the comment "converse-invariant score". Its second
  component, the sorted scores, is not converse-invariant, and it is never used. Delete it.

### m6. Wording

- "Flanked on both sides by nonempty singleton runs" should read "nonempty runs (of singletons)". The residual pairs
  have a run of length 2. The S29 audit (m9) made the same point for HYP-9680.
- "Rigid" in (B) is true at order 9 but incidental: at order 10, one of the 56 seven-core pairs has |Aut| = 3.
  Optionally say so.

---

## Proof check of THM-5155 (A)

- **Scores.** A strong component S is a module, so flipping it changes only the arcs inside S.
  - For v in S, score(v) = s_S(v) + b before the flip, where b is the number of vertices below S. After the flip it
    is |S| − 1 − s_S(v) + b.
  - Every other score is unchanged. Cancelling the unchanged part of the multiset gives the "iff".
  - This half needs no run-flanking and holds for any module.
- **Arc multiset.** This is THM-5154 (A)(ii), which the S29 audit checked.
- **Edge cases.**
  - A self-converse S gives an isomorphic tournament, and its scores are self-reversed automatically.
  - With TT_a ⇒ S ⇒ TT_a, the flip is the global converse. It gives the same merged class, and (A) still holds as
    stated.
  - End cores (one side empty) lie outside (A). In general they change the multiset.
- **Right object.** (A) concerns the sorted score sequence, not the pair. That is the right object: the residual pairs
  have equal, not merely converse, score sequences.

No gap.

---

## Recomputed independently

All auditor code is new. It shares nothing with the packet or with the S29 auditor. It lives in
`C:\Users\Eliott\AppData\Local\Temp\claude\2f1de76a73e81b8a\a7c1f371-6162-4db2-99bd-68121f2beaaf\scratchpad\audit30\`
and was not committed.
- `obs.c` computes, per tournament:
  - the canonical form of T and of T^op by individualisation-refinement (minimum code over all leaves, no pruning),
    with |Aut| equal to the number of minimal leaves;
  - H and every arc weight by plain DFS over Hamiltonian paths, with no split formula;
  - C by DFS from vertex 0, and c3 over all triples;
  - strongness by reachability;
  - Q from the descent distribution A_d (a mask, last-vertex, descent-count DP), as F_k = Σ_d A_d C(n−1−d, k−1−d).
- `analyze.py`, `structure9.py`, `hyptest.py`, `crosscheck_dp.py` (a second arc-weight method) and `crosscheck15.py`.
- `gen10.py`, `an10.py` and `check10.py` for order 10.
- `replay/` holds the packet scripts.

| check | result |
|---|---|
| class lists, n = 3..9 | own canonical forms pairwise distinct on every `cls_n.txt`, with counts 2, 4, 12, 56, 456, 6,880, 191,536 (A000568), so the lists are complete. The converse code of every class is the class of T^op. SC 2, 2, 8, 12, 88, 176, 2,752. H, c3, strongness, the SC flag and \|Aut\| agree on every class |
| `prof_n.txt` | H, C and F_1..F_n agree on every class. F was obtained independently via descents, with A_0 = H and Σ A_d = n! checked |
| `arcms_n.txt` | (H, zeros, nonzero weights) agree on every class. Σ W = (n−1)H, and no non-arc carries weight |
| converse invariance | every observer is equal on T and T^op for every merged class |
| observer table, n = 5..9 | merged 10, 34, 272, 3,528, 97,144. H 7, 19, 77, 320, 1,520. Per-vertex fold 7, 14, 33, 76, 184. Q 8, 24, 156, 1,728, 50,511. Arcs 10, 34, 271, 3,526, 97,121. (H,C,\|Aut\|,c3) 8, 27, 162, 1,633, 29,742. Fold+Q 9, 27, 182, 2,093, 63,686. Fold+arcs 10, 34, 272, 3,528, 97,138. **All as in the README.** Pair-based variants: m1 |
| (B) | (pair, arcs) and (fold, arcs) both give 0 collision groups for n ≤ 8 and 6 pairs at n = 9. The codes are {67436800, 201393152}, {69501184, 134349824}, {75596032, 1140982784}, {134414592, 203522304}, {134482944, 1277266176}, {1077936640, 1140851712}, identical to `score_arcs.py` |
| (B) structure | each pair is TT2 ⇒ S ⇒ TT1 against TT2 ⇒ S′ ⇒ TT1 with S′ ≅ S^op and S′ ≇ S. Core scores are (1,2,2,3,3,4) ×4, (2,2,2,3,3,3) ×1 and (1,1,2,3,4,4) ×1. Every core is rigid and non-self-converse |
| strong 6-tournaments | 35 classes: 9 self-converse and 13 converse pairs. The pairs' score sequences up to reversal are (1,1,2,3,4,4) SR ×1, (1,1,3,3,3,4) ×3, (1,2,2,3,3,4) SR ×4, (1,2,3,3,3,3) ×4 and (2,2,2,3,3,3) SR ×1. So the six residual pairs are **exactly** the self-reversed ones. All 12 such classes are rigid |
| packet replay | `observers.py` and `score_arcs.py` on the author's inputs reproduce the README table and the six codes |
| HYP evidence: S29 P12 | the S29 out-masks and a rebuild from blocks agree. H = 45, the multisets are equal and the score pairs differ, so the HYP sentence is **correct** |
| HYP evidence: Session 320, order 24 | U is vertex-transitive, hence regular (internal score 10), so the flip leaves every labelled score unchanged. The multiset equality is Session 320's theorem and was verified numerically by the S29 audit. The HYP sentence is **correct** |
| M1 and M2 witnesses | as stated above, with two arc-weight methods, canonical forms and brute-force module enumeration |

### New data offered to the author: order 10 (`gen10.py`, `an10.py`, `check10.py`; 45 s)

- **Census.** A non-strong tournament is fixed up to isomorphism by its ordered sequence of strong components. All
  377,107 non-strong classes of order 10 were therefore generated from the verified strong classes of orders 1 and
  3-9. That count is A000568(10) − A051337(10), and the canonical forms are pairwise distinct. They form 188,646
  merged classes, 185 of them self-converse.
- **Arc multiset alone.** There are 188,457 values in 173 groups (157 pairs and 16 triples). This reproduces, for the
  non-strong sector, the S29 audit's full order-10 group count.
- **Joint observer.** (Pair, arcs) and (fold, arcs) both give 188,584 values with exactly **62 collision groups, all
  pairs, all of flip form** TT_a ⇒ S ⇒ TT_b against TT_a ⇒ S^op ⇒ TT_b, with S non-self-converse and its scores
  self-reversed:
  - 6 groups have the same six 6-cores, with runs (3,1);
  - 56 groups have 7-cores, with runs (2,1). These are exactly the 56 of the 137 non-self-converse converse pairs of
    strong 7-tournaments whose score sequence is self-reversed. One of these pairs has |Aut| = 3.
- **Whole order-10 count.** The score sequence decides strongness (Landau), so strong and non-strong classes never
  share a score pair. The S29 audit found that no order-10 multiset group contains a strong class. Together these give
  exactly 62 coincidences of the joint observer at order 10, that is 4,870,920 − 62 = 4,870,858 values. The
  strong-sector input is the S29 audit's; it was not recomputed here.

So THM-5155 (B)'s pattern and HYP-9681, in either reading, hold at order 10. The first known counterexamples to the
filed HYP are at orders 13 (M1) and 15 (M2).
