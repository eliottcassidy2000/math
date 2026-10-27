# Independent audit of THM-4819, parts (A)–(E)

*2026-10-27. Independent adversarial audit by a subagent of mac-mini-2026-10-27-slab, audited at commit `db9fb1c110`,
including (E). The auditor could not write report files from its sandbox, so the producer saved this report verbatim
from its final message.*

**Method.**
- Every proof was checked line by line, including the inherited Session75 §§1–4, the Session74 pruning lemma, THM-4810
  Theorem L, THM-4815 B1/B2/B7/B8/(C), THM-4814 (D) and THM-4812 (c).
- Every number was recomputed with the auditor's own code: nothing from the producer is imported, and no nauty is used
  (the canonical form is the auditor's own). gentourng only lists the classes.
- The scripts and `.out` receipts are in this directory. The binaries were built in a scratchpad and are gitignored.

## Verdict summary

| Part | Verdict |
|---|---|
| (A1) slab inequality and equality case | **PASS** |
| (A2) TT–regular geodesics lie in the slab | **PASS** |
| (A3) (i) ⇔ (ii) ⇔ (iii) | **PASS** |
| (A4) source/sink classes with L ≥ 1 never on a geodesic | **PASS**, one wording slip |
| (B) censuses at orders 5, 7, 9 | **PASS**: every entry reproduced |
| (C) listed instances exactly c moves from TT | **PASS** where checkable; one reproducibility GAP and one descriptive ERROR |
| (D) Lemma B | **PROVED as written** |
| (D) Theorem D | **PROVED as written** |
| (D) remarks | one **ERROR** in a side clause; the data and the cycle-type formula pass |
| (D′) census | **PASS**, reproduced by running the proof's construction |
| (E) Lemma B′ and Theorem E | **PROVED as written for n ≥ 10**; census PASS; minor wording and bookkeeping |

No error affects a proof.

## (A) The slab

- **(A1).** The termwise triangle inequality is right, Σ|i − m| = m(m+1), and the equality case is right. The
  inequality holds on all 12, 456 and 191,536 classes at orders 5, 7 and 9.
- **(A2).** Both scores move by at most one per move, so every step of a word of length c must raise L by one and lower R
  by one.
- **(A3).** All three implications check:
  - (ii) ⇒ (iii) uses R ≤ R(T) + L(T) − j = c − j, which needs T in the slab.
  - (iii) ⇒ (i) holds because ΔR = −1 together with L + R ≥ c forces ΔL = +1.
  - At orders 5 and 7, a full breadth-first search of the class move graph (`slab_audit.c`) confirms that the three
    sets coincide.
- **(A4).** R ≥ m holds, the Session74 pruning lemma holds, and the involution step is valid. The order-5 and order-7
  checks pass.
  - *Wording slip:* "without ever passing through a class with a source or sink" should add "other than TT", since TT
    has both.

## (B) Censuses — PASS

- `slab_closure_audit.c` (own canonical form, forward closure) reproduces the order-9 row exactly:
  - 191,536 classes, 165,342 slab classes, 152,276 on geodesics;
  - all 11 levels and every failure split;
  - 60,056,967 R-descending moves, the same count as the producer's.
- Orders 5 and 7 are reproduced exactly, both by this closure and by the full breadth-first search.
- The single order-7 "other" failure is the class with scores (2,2,2,3,4,4,4), λ = 2 and distance 4 from TT.
- All four observations hold.

## (C) Exact distances

- **P19 — PASS** (`replay_p19.py`). P19 is built from the squares mod 19. All 45 paths are directed Hamilton paths at
  their own stage, L runs 45 → 0, and the end is transitive. (43 moves are black and 2 blue; colour is not claimed.)
- **Logic — PASS.** The producer's replay code is sound.
- **GAP (reproducibility).** Only the P19 word is archived; the search prints a word only in Paley mode. The auditor's
  own search (`geodesic_search_audit.c`, replayed in pure Python by `replay_words.py`; words in
  `geodesic_search_audit.out`) found and validated words for:
  - 40 of 40 regular classes at order 11;
  - all 64 order-13 circulant connection sets;
  - 20 of 20 random regular 13-tournaments;
  - R5[C3], C3[R5] and D(P7) (3 of 3).
  - *Recommendation:* archive the words.
- **ERROR.** The table calls R5[C3] and C3[R5] "THM-4812's vertex-transitive non-Cayley classes". Both are Z₁₅
  circulants (`crt_check.py`):
  - R5[C3] = Z₁₅{1,2,6,7,10,11,12};
  - C3[R5] = Z₁₅{1,4,6,7,10,12,13}.

  So they are already among the 128 circulants. The wording should change to "the G_A- and G_B-invariant classes of
  THM-4812 (c), both also circulants".

## (D) The prescribed black core

### Lemma B — PROVED as written

- **Line by line:**
  - T − v is strong, with every cut ≥ m − 1.
  - σ_t maps v ↔ c_t and c_x ↔ c_{2t−x}, with c_{t+m} fixed.
  - The key identity (★) holds, and is consistent: it gives c_t out-degree exactly m.
  - (a) and (b) follow.
  - The counting closes: for odd m, m − 1 < m; for even m, the mod-4 argument works because 4 divides 2m.
- **Combinatorics** (`startsets_audit.py`, all subsets, m = 2..11). The sets satisfying (a) are empty for odd m and
  exactly the four mod-4 unions for even m. None satisfies (b) for m ≥ 3. At m = 2 the T5 solutions appear, as the
  remark says.
- **Strong form, exhaustive** (`lemmaB_all_cycles.c`). Every vertex and every Hamilton cycle of T − v was tested, and
  (★) was re-checked on every blue rotation with 0 failures. A brute-force cross-check (`lemmaB_crosscheck.py`) matches
  at orders 7 and 9.

| n | tournaments | (T, v, cycle) instances | cycles with all rotations blue | max blue per cycle |
|---|---|---|---|---|
| 5 | T5 | 5 | 5 (expected) | 2 |
| 7 | all 3 regular | 98 | 0 | 2 |
| 9 | all 15 regular | 6,772 | 0 | 2 |
| 11 | all 1,223 regular | 12,366,554 | 0 | 2 |
| 13 | 64 circulant sets | 21,057,972 | 0 | 2 |
| 13 | 100 random regular | 33,135,819 | 0 | 1 |
| 15 | all 128 circulant sets | 1,905,167,640 | 0 | 2 |

### Theorem D — PROVED as written

- **Reuse of Session75 §3.** Everything after the first move uses only three facts about P₁:
  - it is black;
  - it ends at v;
  - its start a is avoided as the start of the second path, which is always possible.

  Nothing depends on how P₁ was chosen, so Lemma B legitimately replaces the pure-blue classification.
- **Gauge and parity.** The conversion between the reversal gauge and native moves is right, colours are preserved, and
  the native word starts at T itself. The terminal is U_m (v a source) when m is even and U_m^op (v a sink) when m is
  odd.
- **Converse.** K_{rev Q}(S^op) = (K_Q(S))^op holds: on Q's pairs both sides are reversed relative to S, and off Q both
  agree with S. So the converse does swap the native source and sink roles.
- **Construction executed** (`construct_DE.c`). The proof was run with random choices, every claim asserted, and the
  resulting word replayed natively from T:
  - orders 7, 9, 11: every (class, vertex, role), 2,100 / 13,500 / 269,060 runs;
  - order 13: 164 tournaments, 85,280 runs;
  - order 15: 228 tournaments, 68,400 runs;
  - 0 failures.

### Remarks on Lemma B

- **ERROR.** "A full (2m+1)-cycle would make T vertex-transitive and hence pure blue" does not follow: P7 is
  vertex-transitive and not pure blue.
  - Without that clause, the alternative route leaves m ≤ 40 ∈ {**4**, 5, 10, 13, 17, 25, 26, 29, 34, 37} unsettled
    (`remark_nocycle.py`).
  - So m = 4 must be added to the printed list, or the clause deleted.
- **PASS:** the cycle-type formula (direct composition, all pairs, m ≤ 40); the data paragraph (14,862 instances adds
  up).
- *Cosmetic:* the `black_endpoint_check.py` docstring says "(C)" where it should say "(D)".
- *Side note for HYP-9426 (7):* the maximum of 2 blue rotations is also attained at order 13 (4,992 circulant cycles);
  its evidence line says 1.

### (D′) — PASS

All 42, 270 and 26,906 instances are confirmed by the construction runs. The auditor also read `prescribed_isolation.c`,
and its logic is sound.

## (E) Even balanced orders

1. **Halving step and the progression claim — PASS.** Exhaustively for m = 3..12, the sets satisfying (a) are exactly the
   2m − 1 translates of 2·{0, 2, …, 2m−4}.
2. **Witnesses — PASS.**
   - Since 2⁻¹ ≡ m, the pair (0, 4) gives −3m ≡ m − 2 and 9m ≡ m + 4.
   - For m ≥ 7 the pair (0, 4) violates (b), checked for m = 7..60.
   - At m = 5 the violating pair is (0, 3); at m = 6 it is (1, 5).
   - (a) and (b) have solutions only at m = 3 and 4, and the listed sets are correct.
3. **Role forcing — PASS.**
   - After j moves, d_j = (−1)^j (j + ½) for a vertex of score m. So it ends as a source iff m is odd, and a score m − 1
     vertex takes the opposite role. The producer's concurrently edited sentence states exactly this.
   - A brute force without pruning never reaches the non-forced role at n = 4 or 6.
4. **Session75 §4 needs only "P₁ black, ending at a score-m vertex" — PASS.** The direction inequality there was only a
   device to get that score. The score m − 1 case via T^op is fine.

**Computations.**
- Lemma B′, strong form, all cycles, v of score m:
  - n = 6: **3 cycles with all rotations blue**, so the lemma fails at 6. This is consistent with the theorem, which
    does not claim n = 6.
  - n = 8: 0 of 4,382; n = 10: 0 of 12,395,604.
  - Self-converse "circulant minus a vertex" tournaments at n = 8, 10, 12, 14: 0 (up to 127,112,392 instances).
- Theorem E construction: n = 8 (all 680 instances × 50 runs), n = 10 (1,333,300 runs), and n = 12–16 samples — 0
  failures. At n = 6 it fails, as expected below the proof's range.
- Pure-Python census (`census_small_even.py`):
  - n = 4: 2 of 4 succeed; the failures are vertex 2 as sink and vertex 3 as source, identical to the producer's;
  - n = 6: 30 of 30; n = 8: 680 of 680.

**Optional strengthening.** Lemma B′ holds at n = 8 by the finite check, and §4 needs only m ≥ 4. So n = 8 is proved by
the proof route, not just the census. Only n = 6 really rests on the census.

**Minor issues.**
- **S4 sentence.** S4 has two score-2 vertices. The sentence is true for vertex 2; vertex 1 has the black path
  (2,3,0,1). Name vertex 2.
- **YAML not updated for (E).** The `status` field omits (E) and still lists HYP-9426 (2) as OPEN, `scripts:` omits the
  even files, and the last sentence of (D)'s "Consequence" is stale.

## Corrections list

1. **(C)** — replace "non-Cayley" with the circulant description above.
2. **(D) remark** — delete or fix the "hence pure blue" clause, and add m = 4 to the list.
3. **(C)** — archive the order-13, order-15 and special-15 words.
4. **(A4)** — add "other than TT".
5. **(E)** — name vertex 2 in the S4 sentence.
6. **(E) bookkeeping** — update the YAML status and scripts, and the stale HYP-9426 (2) references.
7. **Cosmetic** — fix the script docstring and the HYP-9426 (7) evidence line.
