# Independent adversarial audit of opus-2026-11-04-S79 (graded spine grammars G(k, w); THM-4878, HYP-9470)

**VERDICT: PASS WITH CORRECTIONS** (4 majors, 18 minors, 5 upgrades). The mathematics that carries the packet - the
cut-state transfer, the exact ranks, the censuses - is sound and was reproduced with independent code; the corrections
concern credit (same-day and earlier repository prior art contradicting the note's novelty claims), one false
coverage claim inherited by misreading THM-4875 (D), one convention-dependent set of Hankel ranks, and one PROVED-typed
lemma whose written proof contains a false step.

Auditor: independent agent, 2026-11-04 22:30-23:30. Worktree audited: `C:/tmp/wt-s79` (nothing there was modified).
All auditor scripts, binaries and outputs are under
`C:/Users/Eliott/AppData/Local/Temp/claude/800365894db7d24c/f973fbc2-2adf-4ee3-a3b5-b87fe066eb0a/scratchpad/s79_audit/`
(below: `$A`). The packet was still being edited during the audit; a snapshot of the audited files was taken at 22:33
(`$A/snapshot/`, worktree HEAD f0d7bd0c86 plus uncommitted edits to the engine and the note), and the later live
versions (G(0,2) to 8191, G(0,3) to 511, G(1,1) to 2047, exact minimal-arc layers, counting lemma typed PROVED) were
diffed and audited as far as their outputs exist in `C:/Users/Eliott/AppData/Local/Temp/s79/`. Rows marked RUNNING
were not audited.

---

## 1. Majors

**M1. Prior art missed, novelty claims false (credit; sections 0, 1, 3.1, 3.2, 6 of the note; THM-4878 statements 1,
2, 5).** The note states that among the concurrent lanes "none of them grades tournaments by chain cutwidth, and none
computes cut-state automata, Hankel ranks or seam ranks", and calls the G(0,2) completeness "the direction nobody in
the repository had considered". Both statements are contradicted by material that was in the worktree's own history
before the session's first checkpoint (0ac219b0ca, 21:49:37) and long before its SURPRISE commit (f0d7bd0c86, 22:15):

1. `01-canon/theorems/THM-4880-ribbon-tournaments-...md` and
   `05-knowledge/results/ribbon_tournaments_bounded_span_transfer_20261104.md` (mac-mini-2026-11-04-ribbon, commits
   c586c91b31 at 21:40/21:43 and 3fdca9750a at 21:48): a transitive order with back arcs of bounded span; the same
   fragment ("profile") transfer, 7/27/119 states; **Hankel rank** 4 (span 1) and 14 (span 2) by the same
   reachable/co-reachable pairing; an integer representation of rank 62 for span 3 (the same number as S79's three-arc
   seam rank; not explained by either lane); **junction coupling ranks** of the bilinear form (u, v) -> H(u J v),
   with rank 3 identified as the one-arc join (the Lorentzian law) and rank 2 as the near-join - the concept S79 calls
   a seam rank; `{1}^n` (span-1 run) identified with the sole-exit tower and the tribonacci law (note section 3, "The
   one-letter sub-alphabet {1} is the sole-exit tower of THM-4865"); FINITE-EXACT: strong span-2 words attain every
   odd <= 100,000 except 7, 21, 47 and strong span-3 words every odd <= 5,000 except 7, 21; HYP-9473 (span-3
   completeness => Chan-Pak) is the analogue of HYP-9470 (C). Bounded span and bounded cutwidth are incomparable
   gradings (span <= s implies cutwidth <= 2s - 1; cutwidth bounds no span), so S79's G(0,2) result is not subsumed,
   but an apex-free bounded-structure family complete to a large bound was already on the branch, with a larger range.
   The ribbon audit (`ribbon_transfer_audit_20261104.md`, section 10 of the ribbon note) also typed the transfer as
   KNOWN as a method (Knuth TAOCP 4A 7.1.4, frontier-based DP on bounded pathwidth).
2. `05-knowledge/results/overnight_rank_return_20261029/DEFECT_TRANSFER.md` (Codex, 2026-10-29), section D1-D2:
   defines exactly the backward-arc cutwidth of a supplied order, `c_t = #{(u,v) in D: u < t <= v}, c = max_t c_t`
   (D1), and an exact component-pairing transfer processing the vertices in order with fragments kept as ordered
   (head, tail) pairs and endpoint-expiry rules, cost `n * 2^O((k+c) log(k+c+2))` (README.md line 39, MORNING_REPORT
   line 18, RUN_RECURRENCE line 296: "backward-arc cutwidth is at most two, giving O_k(n)"). This is S79's cut-state
   automaton without apexes, with a state bound of the same shape, five days earlier. Also
   `05-knowledge/results/golden_hamilton_transfer_20261103/HANDOFF.md` line 19 ("the supplied backward-cutwidth
   path-forest compiler") and `05-knowledge/results/fibonacci_parity_precision_20261103/COMPLEMENTS.md` lines 4,
   144-156 ("the general bounded-cutwidth query compiler", "backward-arc cutwidth at most 4").
3. `THM-4876 (6)` and `HYP-9469 (3)` (mac-mini, letter MSG-2447 in `agents/opus/inbox/` at 17:35, before the session):
   the one-centre family `TT_m + u` (= the one-apex class) has the closed form `H = 1 + 2 sum_{a in A, b notin A, a<b}
   c(b - a)`, `c(1) = 1, c(d) = 2^(d-2)`, "the family realises about half of the odd numbers (50,006 of the 100,000
   below 200,000; missing 7, 13, 21, 29, 35, 41, 43, 49, ...)". S79's section 5 / THM-4878 (5) presents the half
   density as a new FINITE-EXACT item; it is an extension (to 2^21) of a recorded fact and should cite it.
4. `THM-4731 (6)` (kind-pasteur S11, March 2026): `T*_n` = the transitive tournament with its Hamiltonian path reversed
   (= the converse of the all-spans-1 chain) has `H = A000213(n)`; with THM-4865 (c) and the ribbon note, the "new
   identification of the tower as a cutwidth-1 word" (note 3.2) is not new.

Required: retype THM-4878 (1) as KNOWN as a method (frontier/bounded-cutwidth path-cover DP; repository: DEFECT_TRANSFER
D1-D2, THM-4880 (A)) and PROVED only for the explicit apex-and-slot automaton with its letter alphabet; cite THM-4880 /
HYP-9473 as the concurrent (earlier-committed) bounded-span analogue in sections 0, 1, 6, 7 and HYP-9470 (C); replace
the two quoted novelty sentences; cite THM-4876 (6)/HYP-9469 (3) in section 5 and THM-4731 (6) in 3.2. A ledger entry
of the same kind as MISTAKE-644/699/706 ("same-day lane missed although git log was read") is due: the ribbon commits
precede every S79 checkpoint in the branch history.

**M2. False coverage claim for G(2,1) (THM-4878 statement 2, note sections 0, 1, 3.2, 6; HYP-9470 (B)).** The packet
says "one-apex x one-apex joins attain every odd <= 2^32 except 7, 21 (THM-4875 (D), HYP-9468 CONFIRMED)" and "Hence
G(2,1) attains every odd number <= 2^32 except 7 and 21 ..., which is the strongest finite evidence for any single
grade". THM-4875 (D) says something else: every odd N <= 2^32 except 7, 21 is a count of **F1 or F2**, where F1 = joins
of two one-apex tournaments and F2 = joins of a one-apex tournament with an arbitrary rooted tournament of order <= 7.
The S78 note (`old_mark_bound_lorentzian_gadgets_one_apex_join_coverage_20261104.md`, "The bilinear family alone")
and `04-computation/experiments/chanpak_one_apex_join_cover_20261104.out` record that F1 alone misses 65 of the 512
odd numbers in [2^10, 2^11) (density 0.127), 103 in [2^11, 2^12), ..., 934,107 in [2^31, 2^32) (density 0.00087):

```
F1-only misses per block [2^k, 2^(k+1)), k, count, ratio, density:
  k=10 65 0.0000 0.126953
  k=11 103 1.5846 0.100586
  ...
  k=31 934107 1.5735 0.000870
```

F2's order-<= 7 side is not a G(2,1) word in general, so nothing in THM-4875 (D) gives G(2,1) coverage to 2^32.
Moreover THM-4875's F1 is `T(w1) => T(w2)^op` with the reversed arc `x2 -> x1` **between the apexes**; after
relabelling the apexes it is a two-apex word, i.e. F1 is contained in G(2,0) (the S78 note says so: "F1 lies inside the
two-apex family"), whereas S79's section 3.2 construction joins across a reversed arc between **chain** vertices
`v' -> v` - a different family, whose coverage nobody has computed. Required: delete the 2^32 clause from THM-4878 (2),
the summary table row, section 1, section 7 and HYP-9470 (B); state instead that G(2,0) (hence G(2,1)) contains F1,
which alone has the positive miss density above, and that the finite evidence for G(2,0) is Session 191's 230,519 and
the census to 2047 (reproduced here).

**M3. The apex-free Hankel ranks 5, 18, 80 (and the c = 0 seam rank 2) depend on a convention the packet does not
state, and under the standard convention they are 4, 17, 79 (seam 1).** The engine's initial state (empty cut) is not
accepting, so the series it computes has H(empty word) = 0. The empty tournament has one Hamiltonian path (the empty
path), and H(empty) = 1 is the convention under which the c = 0 cut law is the pure ordered-join product
H(l r) = H(l) H(r) (the note itself describes the extra term as "the degenerate empty-suffix direction"). Recomputing
the pairing from the session's own automaton dumps with my own exact linear algebra, with the empty cut accepting
(`$A/ranks.py --empty1`):

```
dump3_k0_w1.txt: 14 states; reachable dim 8; observable dim 6; Hankel rank 4 (both alphabets); seam c=0: 1, c=1: 3
dump3_k0_w2.txt: 83 states; reachable dim 39; observable dim 25; Hankel rank 17; seams c=0: 1, c=1: 3, c=2: 13
dump3_k0_w3.txt: 521 states; reachable dim 211; observable dim 120; Hankel rank 79; seams 1, 3, 13, 62
```

versus 5 / 18 / 80 and c = 0 seam 2 with the engine's convention (also reproduced by me, section 6 below). The
brute-force Hankel sub-matrices built from H(ur) by subset DP give rank 4 (G(0,1), lengths <= 6) and 17 (G(0,2),
lengths <= 4) under either convention. The grammars with apexes are unaffected (the apex-only tournament has H = 1 in
both conventions): G(1,0) 4, G(1,1) 17, G(2,0) 17, G(1,2) 78 reproduced. All seam ranks with c >= 1 (3, 13, 62, and
13, 61 with one apex) are unaffected. Consequences: the "unexplained coincidence" paragraph (G(0,1), G(0,2), G(0,3)
exceeding G(1,0), G(1,1), G(1,2) by 1, 1, 2) is an artifact; under H(empty) = 1 the three grammars G(0,2), G(1,1),
G(2,0) share Hankel rank 17 and G(0,3) = 79 vs G(1,2) = 78. Required: state the convention in THM-4878 (3) and the
tables; recommended: adopt H(empty) = 1 (rank 4, 17, 79, c = 0 seam 1; G(0,4): 389 expected, see section 6).

**M4. The counting lemma (THM-4878 statement 6, note section 7, typed PROVED) has a false step in its proof and an
imprecise bound.** The proof says "a single arc of span s already gives H >= 1 + 2^(s-1), and the insertion lemma
transfers this to every sub-arc". False: the connected cutwidth-2 word `1,0,1,5` (chain 0..3 with 3 -> 0 and 3 -> 2,
the engine's own first-found witness for 3) has H = 3 although the arc 3 -> 0 has span 3 (1 + 2^2 = 5). The insertion
lemma bounds H below by induced subtournaments, and the induced subtournament on an arc's span contains the other
reversed arcs, so a single-arc bound does not transfer. Also "at most (log_2 X + 1)^(2c) values below X" omits the
constant counting the ways the 2c arc endpoints can be arranged on <= 2c event vertices with cutwidth <= 2. The
conclusion survives with a different proof: in a connected word the vertices between consecutive slot events are
slot-free and inside a pending span, so by the sink-run bound each of the <= 2c + 1 interior gaps is at most
log_2(X - 1), hence at most C(c) (log_2 X + 2)^(2c+1) values below X have c-arc witnesses, with C(c) depending on c
only; this still forces c >= (1 - o(1)) log_2 X / (2 log_2 log_2 X). Required: rewrite the proof, restate the bound
with its constant, or retype the lemma.

## 2. Minors

1. **State bound.** (c+k+2)!(c+k+3)^(c+k) fails at c = k = 0: the dump of G(0, w) has 3 reachable cut states with
   c = 0 (empty cut; the complete fragment (GS, GE); the open fragment (GS, OPEN, mask empty)) against the bound 2.
   The derivation is also wrong in two places: "masks are disjoint subsets of the c + k pending heads" is false for
   apex bits (an apex beats many vertices, so its bit lies in several masks), and the head-group partition `grp`,
   which is part of the state, is not counted. A correct count - for each of the 1 + c + k start labels used or not,
   GE on at most one fragment, the slot bits an injective partial map from slots to fragments, any apex subset per
   fragment, and 2^(c-1) group patterns - is `2^(max(c-1,0)) sum_j C(1+c+k, j) (j+1) (j+1)^c 2^(kj)`, which is below the
   stated bound for every c + k >= 1 and k <= 2 (checked to c = 8); so the statement holds with "c + k >= 1" and the
   proof sketch replaced. Reachable counts per c from the dumps: G(0,2): 3, 11, 69; G(1,1): 15, 99; G(0,3): 3, 11,
   69, 438; G(1,2): 15, 99, 880; G(2,1): 200, 2344; G(0,4): 3, 11, 69, 438, 2863.
2. **Table error.** Note section 6, row "G(0,1) | 2047 | 360 of 1024 | 662 misses": the engine's own output
   `census_k0_w1_X2047.txt` says 358 of 1024 and 666 misses (664 besides 7, 21); my independent census agrees with
   358 (identical value set).
3. **"No bias modulo 16, 3, 5 or 7"** (section 5, THM-4878 (5)) is false for 7: the one-apex miss fraction below 2^21
   by residue mod 7 is 0.455, 0.488, 0.551, 0.581, 0.533, 0.466, 0.443 (residues 0..6; standard error about 0.0013),
   and mod 3 it is 0.501, 0.490, 0.516. Mod 16 and mod 5 are flat (0.501-0.504).
4. **P3 "CONFIRMED" was premature at filing.** At the snapshot only the 13-clause and X = 511 had been run; the
   second clause (at most 10 percent missed below 2000) became testable with `census_k1_w1_X2047.txt` (22:57), which
   confirms it (1022 of 1024). State the time order.
5. **Sink-run lemma, k >= 1 case.** The lemma is stated for a run "followed by a chain vertex that beats some apex";
   the engine also cuts runs that are followed by no such vertex. That is sound for another reason - then every vertex
   of the run is an ordinal-sum cut (prefix and apexes beat it, it beats everything later), so the word's count equals
   that of the shorter word without the run, which is explored - but the lemma should say so. For k = 0 the engine
   applies the cut to every slot-free vertex, which is correct because in a connected word every such vertex lies
   inside a pending span.
6. **Padding lemma wording.** "a first chain vertex that beats every apex is a source" needs "slot-free" (a first
   vertex that opens a slot is beaten by the slot's tail); the engine's `isSource` has it right.
7. **Exact minimal-arc layers (post-snapshot).** The exactness of the minima 1, 8, 75, 133, 36, 1 relies on the
   length cap 64 in addition to "no deduplication"; the needed fact - every word with <= 4 reversed arcs and H <= 511
   has length <= 58 (at most 9 interior gaps r_i with 1 + 2 sum 2^(r_i - 1) <= 511 by the odd-cycle collection
   formula, maximised at seven gaps of 6 and two of 4, plus <= 8 event vertices) - is not stated. (My independent
   exhaustive search to length 11 found no value with fewer arcs than the engine's minimum; the only values whose
   minimal-arc witnesses need length 12, 441 and 445, were verified from the engine's witnesses.)
8. **Validation range inconsistency.** THM-4878 "Validation" says exhaustive agreement for order <= 6; the note and
   commit f0d7bd0c86 say n <= 7 (2,097,152 tournaments).
9. **Alphabet of the series.** The automaton's letters are indices into lettersFor(c), i.e. the same index means a
   different (T, b, A) at different slot counts (in G(0,1): index 1 opens at c = 0 and closes at c = 1). The note
   defines letters as triples. The two alphabets can in principle have different Hankel ranks; I computed both for
   every dumped grammar and they agree, so this is a definitional remark only - but THM-4878 (3) should say which
   alphabet (and which empty-word convention, M3) its "smallest linear representation" refers to.
10. **First-found witnesses of length 64.** 19 of the G(2,0) witnesses in `census_k2_w0_X2047.txt` (values 55, 57,
    75, 81, 95, 135, ...) are words of length 64 consisting mostly of the letter A = 1: those vertices are ordinal-sum
    cuts (beaten by x_1 and the prefix, beating x_2 and the suffix) which the engine's sink rule (A = all apexes) does
    not recognise, so the search pads them to the length cap. The values are correct (verified through the strong
    components), but the "witness" lists are not what a reader expects; the minimal-arc runs fix this for G(0,2).
11. **Closer bound.** The uncommitted engine change at audit time (lines "closer bound: ... the least such count
    bounds H below") is sound: the first future vertex u that terminates a pending slot of the prefix induces, with
    the prefix and the apexes, the prefix extended by the letter (T_u, 0, A_u) with T_u nonempty and admissible, so the
    minimum over such letters is a valid lower bound by the insertion lemma; letters with inadmissible T give an empty
    vector and are skipped correctly. It must be committed together with the 2047/8191/G(0,3) outputs that used it.
12. **Dedupe and sink-run interplay.** The dedupe key is the count vector alone, not (vector, run length); this is
    sound (the vector fixes c and all future values, and a run cut either exceeds X or is padding with the same count
    as an explored shorter word), but the argument is not in the note. The two 64-bit hashes make a false dedupe
    astronomically unlikely but not impossible; a collision could only lose values (never create them), so attainment
    claims are unaffected and the miss claims for G(1,0), G(0,1) are cross-checked by other methods.
13. **Hankel rank of G(2,0) and "rational dimension 28".** Correctly attributed to Session 191
    (`square_sum_spectrum_overnight_20261030/README.md` lines 51-65: Hankel rank 17, reachable 28, observable 20 on
    36 states). Fine; mention that S79's 200-state representation has reachable dimension 125, observable 25.
14. **Section 7's "no general decision procedure"** sentence cites the squares example; a reference for the
    undecidability of range questions for N-rational series would be better than an example.
15. **HYP-9470 (D)** ("universal over all tournament pairs") is a hypothesis about the two-arc composition law; the
    THM-4880 junction coupling ranks (1, 2, 3, 4, 8, 9, 10, 12 for span-2 blocks) are the data to compare with and
    should be cited.
16. **Credit for the all-parts-1 identification** (see M1 item 4) and for the "parts" convention: THM-4744's parts
    p_i are spans (the interval has p_i + 1 vertices), so "all parts 1" means all spans 1; say so.
17. **G(0,4) numbers** (3,384 states, Hankel rank 390, seam 310) are "modular" and were not committed with a dump in
    the packet directory list; the dump exists in the temp directory (`dump3_k0_w4.txt`); my modular recomputation is
    reported in section 6 (pending at the time of writing if not listed there).
18. **Housekeeping.** The note's section 9 "Audit: PENDING" and THM-4878 "AUDIT: PENDING" should be replaced by this
    report's verdict and the corrections applied.

## 3. Upgrades

1. **Seam rank 13 is verified from both sides.** A brute-force Hankel block (prefixes ending at a two-slot cut and
   complete suffixes from it, lengths <= 5 each, 1,592,001 subset-DP counts, 1242 x 1465 matrix) has rank 13 modulo
   two primes (a lower bound), and the exact pairing from the automaton gives 13 (an upper bound). Likewise seam 3 at
   c = 1 (lower bound 3 from lengths <= 5 in G(0,1) and G(0,2)), Hankel 4 for G(1,0) (exact brute force, 127 x 126)
   and 17 for G(1,1) (lower bound from lengths <= 4; automaton 17). "FINITE-EXACT (exact rational arithmetic)" can
   add "independently confirmed by brute-force Hankel blocks".
2. **Coverage reproduced with independent code.** By my own subset-DP searches (no engine code reused): G(0,2)
   attains every odd <= 511 except 7, 21 (exhaustive to chain length 11 plus random words of length 9-12), and every
   odd in [25, 2047] (random words of length 9-13) - together with the exhaustive part, every odd <= 2047 except 7,
   21; G(1,1) every odd <= 511 except 7, 21; G(0,3) every odd <= 511 except 7, 21; G(2,0) every odd <= 2047 except
   7, 21; G(1,0) exactly 560 of 1024 directly and 763 with products below 2048 (complete census: a reduced one-apex
   word of length m has H >= 1 + 2^(m-2)); G(0,1) exactly 358 values below 2048 (identical to the engine's set) and
   181 of 500 below 1000; the one-apex count 521,915 of 1,048,576 below 2^21 with the dyadic miss sequence exactly as
   stated. All 254 minimal-arc witnesses to 511 and all listed witnesses (<= 301) of every finished census file were
   decoded with my own decoder and recounted (0 mismatches), including the 19 length-64 two-apex words (via strong
   components).
3. **Rank coincidence restated.** With H(empty) = 1: rank G(0,2) = rank G(1,1) = rank G(2,0) = 17, and
   rank G(1,2) = rank G(2,1) = 78 (the latter pending my modular run), rank G(0,3) = 79. The pattern "an apex and a
   slot contribute alike once one apex is present" then reads "a second radix of either kind gives 17; a third gives
   78/79".
4. **Exact minimal-arc layers below 512** (1, 8, 75, 133, 36, 1) are consistent with an independent exhaustive search
   to length 11 (which never beats the engine's minimum for any value) and the two length-12 witnesses were verified.
5. **One-apex formula.** The odd-cycle derivation H = 1 + 2 sum_{i<j, b_i = 1, b_j = 0} w(j - i), w(1) = 1,
   w(d) = 2^(d-2) (THM-4876 (6)) reproduces the DP on all 8,190 words of length <= 12 and gives the 2^21 census in
   1.2 s; the note's "every one-apex value <= 2^21 is attained by a word of length <= 24" follows from
   H >= 1 + 2 w(m - 1) for reduced words (m <= 22).

---

## 4. Proof check of Theorem A (note 3.1; THM-4878 statement 1)

Read against the engine `applyLetter` / `applyApex` (snapshot), the transfer is correct; details:

- **Predecessor options** of the new chain vertex v: (a) none (GS, once - `hasGS`); (b) the end of an OPEN fragment
  whose mask is disjoint from the terminated set T (f_r -> v is a forward arc unless v is the tail of a slot with
  head f_r; the mask of an open end is exactly the set of slots with that head plus the apexes beating it, so the test
  `(e & tailMask) == 0` is the arc condition); (c) a later chain vertex, which beats v only through a reversed arc
  with head v, i.e. one of the b new slots - the slot label is the identity of the later tail; (d) an apex in A not
  already used as a start label (`apexUsed`): an apex has one successor. No option is missing: a predecessor is either
  absent, in the prefix (an open end), a later chain vertex (necessarily a tail of a new slot), or an apex.
- **Successor options**: (a') none (GE, once); (b') later (OPEN with mask = new slots + apexes in A: exactly the later
  vertices that beat v, so a later chain vertex u can follow v iff u is not a tail of a slot with head v, and an apex
  can follow v iff it does not beat v - both are checked when u or the apex is placed); (c') the first vertex of a
  fragment whose start label is a slot in T (v is that tail; the arc v -> head is reversed, hence present).
- **No double counting.** Two new slots with the same head v get distinct labels c, c+1; the word fixes which later
  vertex terminates which label (the bits of later T's), so the two predecessor choices are different paths of the
  same tournament. The same tournament is reached by several words (the two labels can be swapped), which costs
  search time but not correctness. Apex entries: at most one fragment per apex label, and the apex step joins the
  apex to that fragment or kills the state.
- **Dead-state rule.** A fragment labelled by a slot terminated at v but not chosen as v's successor can never be
  entered (its only possible predecessor was v): correctly discarded (`ok = false`). The analogous apex rule in
  `applyApex` (a fragment labelled a surviving the step of apex a) is correct.
- **Head-group rule.** Terminating two slots of one group at one vertex means two arcs v -> h: not a tournament; the
  engine returns no transitions (the word has H = 0, consistent with "no such tournament"). Terminating slots of
  different groups at one vertex (v -> h_1, v -> h_2) is allowed. This is exactly the tournament condition.
- **Cycle check.** (b) and (c') naming the same fragment would close a cycle; excluded (`sc == pc`).
- **Apex steps.** x_1 then x_2 after all chain letters: predecessor = GS, or an open end whose mask lacks the apex
  (the end beats the apex), successor = GE, open (mask 0, so that x_2 can follow x_1 and nothing else can), or the
  fragment labelled by the apex. x_2 -> x_1 is impossible because x_1 is placed first with no apex-predecessor option.
- **Converse direction and the final state.** A choice sequence ending in the single fragment (GS, GE) assembles a
  Hamiltonian path using only arcs of T(word); the count of accepted sequences is the matrix product. Correct.
- **Pruning (finishValue > X)**: the prefix chain plus the apexes is an induced subtournament of every completion, so
  the insertion lemma applies; the complete-state count after the apex steps on a prefix with pending slots equals
  H(prefix + apexes) because configurations entered through a pending slot are never complete. Correct.
- **State bound**: see minor 1 (false at c = k = 0; derivation wrong; a corrected count is below the stated bound for
  c + k >= 1).
- **Sink-run bound (3.4)**: for k = 0 the run u_1..u_r lies inside a pending span (head h < u_1, tail u > u_r), and for
  every odd subset S the cycle u -> h -> S (increasing) -> u uses only forward arcs between slot-free vertices and the
  reversed arc u -> h, so it is a cycle of odd length |S| + 2; the 2^(r-1) odd subsets give distinct odd cycles through
  u, and the odd-cycle collection formula (THM-4876, Grinberg-Stanley) gives H >= 1 + 2 * 2^(r-1). For k >= 1 with a
  later vertex u beating an apex x the cycle x -> S -> u -> x works the same way. Correct, with the caveat of minor 5.
- **Saturation**: an entry that contributes to a final count <= X is bounded by that count (its configurations extend
  to Hamiltonian paths under the same suffix), so no contributing entry ever reaches 2^62; a saturated entry can only
  propagate to a saturated (rejected) final sum. `satMul` saturates slightly early (a >= CAP/m) which is conservative.
  Correct.
- **Closer bound** (post-snapshot engine): sound, see minor 11.

## 5. Independent reproductions (commands and outputs)

Tool: `$A/hp.cpp` (compiled to `hp.exe`, later `hp2.exe` with a random-search mode), a brute-force subset-DP
Hamiltonian-path counter with its own letter enumeration/decoder written from section 2 of the note; `$A/decode.py`
(a second, Python decoder plus the witness checker); `$A/oneapex.cpp` (odd-cycle formula). The engine was used only in
its word mode as the thing being tested, and its dumps only as input for the rank checks (which were also done without
them, section 6).

### 5.1 One-apex miss list and 560 of 1024 (task 2a)

```
$ ./hp.exe census 1 0 2047 13 oneapex_2047_L13.txt
census G(1,0) X=2047 L<=13: 3469 DFS nodes, 4889 completed words; odd values attained directly: 560 of 1024
with products: 763 of 1024; 261 misses: 7 13 21 29 35 41 43 49 59 67 71 77 83 89 101 105 109 113 119 131 133 137 ...
misses per dyadic scale: [2^2:1/2] [2^3:1/4] [2^4:2/8] [2^5:5/16] [2^6:10/32] [2^7:21/64] [2^8:38/128] [2^9:62/256] [2^10:121/512]
```

(259 misses besides 7, 21, as the note's table says; first misses 7, 13, 21, 29, 35, 41, 43 as in V2.) The census is
complete: a reduced one-apex word (leading 1, trailing 0) of length m has 2^(m-3) odd cycles through x, v_1, v_m, so
H >= 1 + 2^(m-2) and m <= 12 for H <= 2047; L = 13 was run.

Section 5 of the note (2^21), by the closed formula verified against the DP:

```
$ ./oneapex.exe verify 12
one-apex formula vs subset DP: 8190 words of length <= 12, 0 mismatches
$ ./oneapex.exe census 2097151 24
one-apex census X=2097151 L<=24: 25499247 nodes; odd values attained directly: 521915 of 1048576
with products: 803615 of 1048576
direct misses per dyadic scale [2^j, 2^(j+1)): 1 1 2 7 14 29 57 117 236 495 998 2035 4077 8180 16382 32923 65646 131959 263502
direct miss fraction by residue mod 16: 1:0.5019 3:0.5010 5:0.5020 7:0.5037 9:0.5038 11:0.5033 13:0.5012 15:0.5013
direct miss fraction by residue mod 3: 0:0.5007 1:0.4900 2:0.5161
direct miss fraction by residue mod 5: 0:0.5030 1:0.5010 2:0.5025 3:0.5031 4:0.5017
direct miss fraction by residue mod 7: 0:0.4546 1:0.4875 2:0.5509 3:0.5813 4:0.5326 5:0.4658 6:0.4431
```

### 5.2 Random words decoded from the letter semantics (task 2b)

`decode.py random k w N maxlen seed`: my Python decoder + my DP (`myDP`), my C++ decoder (`hpWord`), and the engine's
word mode `grammar4.exe k w 1000 64 1 0 "indices"` (`engine`); files `$A/random_k*.txt`.

```
G(0,2): 40 random words compared, 0 mismatches     (e.g. word 1,3,2,0,0,2,2,0,5 m=9 arcs=[(1,0),(2,1),(6,5),(8,5),(8,6)] H=15)
G(1,1): 40 random words compared, 0 mismatches
G(2,1): 30 random words compared, 0 mismatches     (e.g. word 4,3,1,0,8,3,5,1 arcs=[(4,0),(6,4)] apex=03100311 H=1273)
G(0,3): 30 random words compared, 0 mismatches
```

### 5.3 The all-parts-1 chain and the sole-exit tower (task 2c)

```
$ for r in 1..16: ./hp.exe chain $r
H = 1, 3, 5, 9, 17, 31, 57, 105, 193, 355, 653, 1201, 2209, 4063, 7473, 13745
```

exactly the note's list (A000213(r+1)). The chain (arcs i+1 -> i, all other arcs forward) is the converse of THM-4731
(6)'s T*_n and is THM-4865's tower over a point (x_1 -> t, x_j -> x_(j-1), everything else -> x_j): the last vertex beats
only its predecessor, so the tower identification is correct - and already recorded (M1 item 4).

### 5.4 H = 13 (task 2d)

```
$ ./hp.exe spine 6 "3>0" "5>1" "5>4"     -> H = 13
$ ./hp.exe spine 5 "3>0" "4>1"           -> H = 13
$ ./hp.exe word 2 0 1,2,0                -> m=3 n=5 apex:120  H = 13
$ ./hp.exe spine 7 "2>0" "6>2" "6>3"     -> H = 35
$ ./hp.exe spine 7 "2>0" "3>1" "4>2" "6>3" "6>4" -> H = 67
```

### 5.5 Census witnesses (task 2e) - all of them, not 10

`decode.py witnesses k w file` decodes every listed witness with my decoder and recounts H as the product of the
subset-DP counts of its strong components:

```
census_k0_w2_X511.txt:        149 witnesses, 0 mismatches; arcs {2:24, 3:41, 4:46, 5:32, 6:4, 7:1}; chain length <= 12
census_k1_w1_X511.txt:        149 witnesses, 0 mismatches
census_k2_w0_X2047.txt:       149 witnesses, 0 mismatches (19 of them words of length 64, minor 10)
census_k1_w0_X2047.txt:        87 witnesses, 0 mismatches
census_k0_w1_X2047.txt:        63 witnesses, 0 mismatches
census_k0_w2_X2047.txt:       149 witnesses, 0 mismatches
census_k0_w2_X8191.txt:       149 witnesses, 0 mismatches
census_k0_w2_X511_minfa5.txt: 254 witnesses, 0 mismatches; arcs {1:8, 2:75, 3:133, 4:36, 5:1}; chain length <= 12
census_k0_w3_X511.txt:        149 witnesses, 0 mismatches
census_k1_w1_X2047.txt:       149 witnesses, 0 mismatches
```

The arc histogram of the first file is exactly the note's "24 values with 2 arcs, 41 with 3, 46 with 4, 32 with 5, 4
with 6, 1 with 7; chain length <= 12".

### 5.6 Independent censuses and witness searches

```
$ ./hp.exe census 0 1 2047 18 g01_2047_L18.txt
census G(0,1) X=2047 L<=18: 40636 DFS nodes; odd values attained directly: 358 of 1024   (engine: 358; sets identical)
  (181 of the 500 odd numbers below 1000: V1 confirmed)
$ ./hp.exe census 0 2 511 9  -> 220 of 256 direct, 238 with products; misses 7 21 199 331 371 397 413 439 449 457 463 479 487 491 497 499 509 511
$ ./hp.exe census 0 2 511 11 -> 6173887 DFS nodes; 254 of 256; misses 7 21; minimal arcs within L<=11: fa0:1 fa1:8 fa2:75 fa3:131 fa4:38 fa5:1
$ ./hp2.exe randsearch 0 2 511 9 12 400000 11  -> found 245 of 256; not found: 3 5 7 9 11 13 15 17 19 21 23   (short values need short words)
   union of the exhaustive L<=9 census and the random search: 254 of 256, missing exactly [7, 21]
$ ./hp2.exe randsearch 1 1 511 5 11 300000 12  -> 254 of 256; not found: 7 21
$ ./hp2.exe randsearch 0 3 511 6 11 200000 41  -> 254 of 256; not found: 7 21
$ ./hp2.exe randsearch 2 0 2047 5 11 300000 51 -> 1022 of 1024; not found: 7 21
$ ./hp2.exe randsearch 0 2 2047 9 13 200000 21 -> 1013 of 1024; not found: 3 5 7 9 11 13 15 17 19 21 23 (all covered by the L<=11 census)
```

Per-value comparison of my L <= 11 minima with the engine's exact minima (`census_k0_w2_X511_minfa5.txt`): my minimum
is never smaller; it is larger only for 441 (mine 4 arcs at length 9, engine 3 at length 12) and 445 (same), whose
engine witnesses were verified.

### 5.7 State bound numerics

From the dumps (section 2 minor 1): reachable states per c versus the note's bound (c+k+2)!(c+k+3)^(c+k): G(0,1):
c=0: 3 (bound 2!), c=1: 11 (24); G(0,2): 3, 11, 69 (2, 24, 600); G(1,1): 15, 99 (24, 600); G(2,0): 200 (600); G(0,3):
3, 11, 69, 438 (..., 25920); G(1,2): 15, 99, 880; G(2,1): 200, 2344; G(0,4): 3, 11, 69, 438, 2863 (1,728,720).

## 6. Ranks (task 3)

Dumps re-generated from the snapshot source (`$A/eng.exe k w 100 -1`) are byte-identical to `dump3_k1_w0.txt`,
`dump3_k0_w1.txt`, `dump3_k0_w2.txt`, `dump3_k1_w1.txt`; `dump_k2_w0.txt` was generated the same way. My own exact
(Fraction) closure `$A/ranks.py` (reachable rows from state 0, observable columns from the complete indicator after
the apex steps, pairing rank; both the index alphabet and the (T,b,A) alphabet; seams by restricting rows to c-cuts):

```
dump3_k1_w0 (G(1,0)): 15 states; reachable 9, observable 6; Hankel 4 / 4; seam c=0: 4
dump3_k0_w1 (G(0,1)): 14 states; reachable 8, observable 6; Hankel 5 / 5; seams c=0: 2, c=1: 3
dump3_k0_w2 (G(0,2)): 83 states; reachable 39, observable 25; Hankel 18 / 18; seams 2, 3, 13
dump3_k1_w1 (G(1,1)): 114 states; reachable 65, observable 25; Hankel 17 / 17; seams c=0: 4, c=1: 13
dump_k2_w0  (G(2,0)): 200 states; reachable 125, observable 25; Hankel 17 / 17; seam c=0: 17
dump3_k0_w3 (G(0,3)): 521 states; reachable 211, observable 120; Hankel 80 / 80; seams 2, 3, 13, 62
dump3_k1_w2 (G(1,2)): 994 states; reachable 505, observable 116; Hankel 78 / 78 (mod 2^31-1); seams 4, 13, 61
dump3_k2_w1 (G(2,1)): 2544 states; reachable 1549, observable 116; Hankel 78 / 78 (mod 2^31-1); seams c=0: 17, c=1: 61
dump3_k0_w4 (G(0,4)): see the addendum (modular runs still in progress when this report was closed)
```

i.e. every number of the note's rank tables that I could recompute is reproduced under the engine's convention; with
H(empty) = 1 (M3) the apex-free ranks become 4, 17, 79 and the c = 0 seams 1.

Brute-force Hankel matrices (`$A/hankel_bf.py`, entries H(u r) by subset DP, independent of the automaton; ranks mod
2^31-1 and 2^31-19, exact where small):

```
G(1,0): 127 prefixes x 126 suffixes (lengths <= 6): rank 4 (exact)
G(0,1): 377 x 376 (lengths <= 6): rank 4 (exact; 4 also with the entry H(empty) := 0)
G(0,1) seam c=1: 88 x 88 (lengths <= 5): rank 3 (exact)
G(0,2) seam c=1: 1019 x 796 (lengths <= 5): rank 3
G(0,2): 477 x 320 (lengths <= 4): rank 17
G(0,2) seam c=2: 1242 x 1465 (lengths <= 5): rank 13
G(1,1): 673 x 672 (lengths <= 4): rank 17
```

Lower bounds from finite blocks meeting the exact upper bounds from the automaton: Hankel 4 (G(1,0)), 17 (G(1,1)),
and - under H(empty) = 1 - 4 (G(0,1)), 17 (G(0,2)); seams 3 (c = 1) and 13 (c = 2) are therefore verified from both
sides. Not verified by me: the seam 62 and 61 (only from the automaton, exact), the G(2,1) Hankel 78 and the G(0,4)
numbers (modular, pending/addendum), the G(1,2) seam 61 (modular only).

## 7. Prior art (task 4) - full listing of the hits that matter

Keyword sweep (`git grep -il`, no truncation; lists in `$A/priorart/*.lst`): cutwidth 25 files, "cut state" 35,
"seam rank" 8 (all S79), "coupling rank" 32, tribonacci 757, Hankel 1029, "sole-exit tower" 34, "upset tournament"
23, "two reversed arcs" 45, ribbon/bounded span 345, "half of the odd" 47. Relevant earlier statements:

| claim in S79 | earlier file | lines |
|---|---|---|
| cutwidth grading / backward-arc cutwidth of an order | `05-knowledge/results/overnight_rank_return_20261029/DEFECT_TRANSFER.md` (D1) | 28-37 |
| cut-state (component-pairing) transfer with a 2^O((k+c)log) state bound | same, D2-D4; `README.md` 39; `MORNING_REPORT.md` 18; `RUN_RECURRENCE.md` 296 | |
| bounded-cutwidth path-forest compiler | `05-knowledge/results/golden_hamilton_transfer_20261103/HANDOFF.md` | 19 |
| bounded-cutwidth query compiler, cutwidth <= 4 targets | `05-knowledge/results/fibonacci_parity_precision_20261103/COMPLEMENTS.md`; `README.md` | 4, 144-156; 30 |
| finite-state transfer for H on a transitive order with bounded back arcs; Hankel ranks 4, 14; rank 62 (span 3) | `01-canon/theorems/THM-4880-...md` (A); `05-knowledge/results/ribbon_tournaments_bounded_span_transfer_20261104.md` section 3 | |
| composition/coupling ranks across blocks, rank 3 = one-arc join, 2 = near-join | ribbon note section 8; THM-4880 last paragraph | |
| tribonacci tower = the span-1 / cutwidth-1 run | ribbon note section 3 ("{1} is the sole-exit tower"); `THM-4731-...md` (6) line 76; THM-4865 (c) | |
| one-apex density one half, closed formula, first misses 7, 13, 21, 29, ... | `THM-4876-...md` (6) lines 107-112; `HYP-9469-...md` (3) lines 40-43; `chanpak_increments_thm013_thm017_merge_20261104.md` 93-95 | |
| apex-free bounded family covers all odd numbers except 7, 21 (to 5,000 / 100,000 except 47) | THM-4880 FINITE-EXACT; HYP-9473 (1)-(2) | |
| seam gadgets with transfer matrices, "seam" vocabulary | `THM-4879-...md`, `05-knowledge/results/seam_gadget_thickness_spectrum_20261104.md` (notes S79 as "concurrent and not read", line 91) | |

No earlier file states the cut-state automaton **with apexes**, the specific numbers 13 / 62 / 310 as seam ranks
across c reversed arcs, or the G(0,2) census; those are new. The one-arc seam rank 3 is the Lorentzian law (session29
ADDENDUM, THM-4866 (B)) as the note says, and also THM-4880's junction rank 3.

## 8. Status typing (task 5)

| item | typed | audit |
|---|---|---|
| THM-4878 (1) transfer | PROVED | correct as mathematics; KNOWN as a method (frontier DP; DEFECT_TRANSFER D2; THM-4880 (A)); state bound false at c = k = 0, derivation wrong (minor 1) |
| (2) dictionary: G(1,0), G(2,0), connected G(0,1), fa <= 2, 13's grades | PROVED | correct (13: verified in G(0,2), G(2,0); absent from G(1,0) by my complete census and from G(0,1) by my exhaustive census agreeing with the engine) |
| (2) "G(2,1) attains every odd <= 2^32 except 7, 21" | PROVED | **false inference** (M2); F1 is in G(2,0) and alone has positive miss density |
| (2) tribonacci identification "new" | PROVED | correct, not new (M1) |
| (3) Hankel ranks | FINITE-EXACT | reproduced; apex-free values convention-dependent (M3); alphabet to be named (minor 9) |
| (3) seam ranks 3, 13, 62; 13, 61; 17, 61 | FINITE-EXACT | 3, 13 confirmed from both sides; 62, 61 reproduced exactly (61 modular); c = 0 values depend on the convention |
| (4) sink-run bound, padding, saturation | PROVED | correct (minors 5, 6) |
| (5) one-apex density | FINITE-EXACT | reproduced; extends THM-4876 (6)/HYP-9469 (3), not cited (M1); "no bias mod 7" false (minor 3) |
| (6) censuses: G(0,2) 511 / 2047 / 8191, G(1,1) 511 / 2047, G(0,3) 511, G(2,0) 2047 | FINITE-EXACT | attainment claims are existential and are proved by the witnesses; to 511 (all grades), to 2047 (G(0,2), G(2,0)) they were reproduced by my own searches; 8191 only through the engine (witness list stops at 301) - type "FINITE-EXACT by the engine, witnesses <= 301 and the minimal-arc list to 511 independently recounted". The "misses" of G(1,0) and G(0,1) are complete by length bounds (one-apex) and by two independent censuses (G(0,1)) |
| (6) exact minimal-arc layers 1, 8, 75, 133, 36, 1 | exact | consistent with my independent search; the length-64 argument is missing (minor 7) |
| (6) counting lemma | PROVED | proof defective (M4) |
| note 3.3 "universal" 13 and 62 | FINITE-EXACT within cutwidth <= 3 | fine as stated (HYP-9470 (D)) |
| validation "480 random words", "n <= 7 exhaustive" | FINITE-EXACT | my own 140 random words agree; THM says order <= 6 (minor 8) |

On the specific question: "G(0,2) attains every odd number <= 511 except 7 and 21" is an existential statement
certified by explicit witnesses (the minimal-arc file lists all 254 with chain length <= 12), so it needs no
completeness device at all; the length-64 / sink-run / saturation / closer devices matter only for the **miss** claims
and for the **exact minima**, where they are sound (section 4) but the length-64 cap needs the gap argument of minor 7.

## 9. Pre-registered predictions (task 6)

| | prediction | outcome recorded | audit |
|---|---|---|---|
| P1 | G(0,2): 10-60 misses <= 300; >= 3 percent missed in [2^11, 2^12) | REFUTED (first clause at the snapshot; both clauses after the 8191 run) | correct: 0 misses <= 300 (my census to length 9 misses only 199 below 300, attained at length 11; engine witnesses verified); the second clause needed X >= 4095 and is refuted by `census_k0_w2_X8191.txt` (0 misses in [2^11, 2^12)), which I could not reproduce independently beyond 2047 |
| P2 | G(0,3): <= 5 misses <= 300; >= 1 miss <= 2000 | first clause CONFIRMED (post-snapshot), second pending | correct: `census_k0_w3_X511.txt` 254 of 256, reproduced by my random search |
| P3 | G(1,1): 13 attained; <= 10 percent missed <= 2000 | CONFIRMED | 13 = word 1,2,0,3 verified; the 2000-clause was confirmed only by `census_k1_w1_X2047.txt` at 22:57 (1022 of 1024), after the note first said CONFIRMED (minor 4) |
| P4 | G(1,2): no miss <= 2000 | PENDING | file empty; not audited |
| P5 | 13 in G(0,2) and G(2,0), not in G(1,0) or G(0,1) | CONFIRMED | correct (section 5.4; one-apex and G(0,1) censuses) |

## 10. Addendum: results that arrived while the report was being written

- **G(2,1)**: my modular closure on `dump3_k2_w1.txt` gives reachable 1549, observable 116, Hankel rank 78, seams 17
  (c = 0) and 61 (c = 1) - the note's numbers (modular both there and here; a rational rank >= 78 is certain, equality
  is as certain as two agreeing primes).
- **G(0,2) to 8191 (post-snapshot claim)**: my random-witness search over words of length 11-16 (150,000 trials,
  119,673 valid) found 4045 of the 4096 odd values; the 41 not found are 3..127 (short values, all covered by the
  exhaustive census and the earlier runs) and the ten values 5297, 6607, 6611, 6629, 7019, 7583, 7867, 7879, 7883,
  8155 near the top of the range, for which the engine's witnesses are not listed (the witness output stops at 301).
  So 4086 of the 4096 are independently reproduced; the claim for the remaining ten rests on the engine alone (a
  second random run with lengths 12-15 was launched; its result, if it finished, is in `$A/rand_g02_8191b.out`).
  Recommendation: print the witnesses for every value (or at least their minimal-arc words) so that the 2047 and
  8191 claims can be certified by recounting, as the 511 claim now is.
- **G(0,4)** (3,384 states): my modular runs under both conventions were still running at close; expected 390
  (engine convention) and 389 (H(empty) = 1); see `$A/ranks_k0_w4_mod.out` and `$A/ranks_k0_w4_mod_empty1.out`.
- **G(2,0) brute-force Hankel** (lengths <= 5): still running at close (`$A/hankel_bf.out`, last line); the exact
  automaton value 17 was reproduced from the dump.

## 11. Files

- `$A/hp.cpp`, `hp.exe`, `hp2.exe` - subset-DP counter; modes batch / word / spine / chain / census / randsearch
- `$A/decode.py` - Python decoder, random-word comparison, witness checker (strong-component products)
- `$A/oneapex.cpp`, `oneapex.exe` - one-apex odd-cycle formula and 2^21 census
- `$A/ranks.py` - exact/modular closure ranks from dumps, both alphabets, `--empty1`
- `$A/hankel_bf.py` - brute-force Hankel blocks
- outputs: `oneapex_2047_L13.txt`, `oneapex_2pow21.out`, `g01_2047_L18.txt`, `g02_511_L9.txt`, `g02_511_L11.txt`,
  `rand_g*.txt`, `random_k*.txt`, `ranks_*.out`, `hankel_bf.out`, `hankel_bf_empty0.out`, `mydumps/`, `priorart/`,
  `snapshot/` (the audited versions of the packet files, 22:33)
