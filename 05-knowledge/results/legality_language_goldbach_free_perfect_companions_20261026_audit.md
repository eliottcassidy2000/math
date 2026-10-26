# Independent adversarial audit of opus-2026-10-26-S61 (commit 52e3aa76e2; strong-form addition c71ab459f0)

Auditor: independent (Claude Fable 5.1), read-only on the worktree `/tmp/math-wt-s61`; all scripts and outputs in this directory. Files audited: the note `05-knowledge/results/legality_language_goldbach_free_perfect_companions_20261026.md`, THM-4806, HYP-9417, the four S61 scripts with their outputs, the pre-registration file, the index rows, the ledger bullet, and (per the coordinator's addition) the strong-form paragraph, THM-4806 R "Strong form", HYP-9417 (1b) and `legality_strong_lemma_20261026.cpp` (its `.out` was still empty when read).

## Verdict: PASS WITH CORRECTIONS

Every proof in the note checks line by line and every number I could recompute reproduces exactly (tables below). The corrections are about inheritance, typing and scope, and two of them are mandatory because they concern the packet's headline labels:

1. The all-order injectivity that the note, THM-4806 and HYP-9417 (1) label **OPEN** was already **PROVED and independently audited in the same tree** 42 minutes before S61 filed (session70, commit a8fcb1c840, 15:18:11). The note quotes a superseded checkpoint-1 README (15:04) as if current.
2. A second concurrent packet on the identical question (mac-mini-2026-10-25-pyth, commit 9f2a9722ef, 15:37:59, `hamilton_path_recovery_20261026`, with THM-4807 reserved) is not cited although its reservation stub sits in the very INDEX line that S61 edited.
3. The strong form (c71ab459f0) is not a hypothesis: it **follows in six lines from the facts already established in the note's own proof of R (i)**, so HYP-9417 (1a)/(1b) are theorems (consistent with session70's closure). Proof in "Majors", item 3.
4. THM-4806 P (iii) ("covers at most four odd prime factors otherwise") overstates what `perfect_companions_shapes_20261026.cpp` enumerates for p = 19 and 31.

## Majors

**M1. Stale OPEN label and stale inheritance (sections 0, 1, 7, 9; THM-4806 status; HYP-9417 (1); ledger bullet).** Timeline from `git log` (all 2026-10-26, -0600):

| time | commit | lane | content |
|---|---|---|---|
| 15:04 | 1988aba492 | session70 | checkpoint 1: "The all-order tournament reconstruction and graph-labeling lanes remain in progress at this checkpoint." |
| 15:18 | a8fcb1c840 | session70 | "[session70-carriers]: close all-order named guard reconstruction"; README line 1 "legal guards reconstruct tournaments"; `session70-guards/REPORT.md` = "Legal re-marking guards reconstruct every finite tournament", PROVED relative to the cited Busch-Jacobson-Reid nontraceable-arc theorem, with an independent guard audit (`session70-guard-audit/REPORT.md`, Lemmas A-C, decoder replay) |
| 15:37 | 9f2a9722ef | mac-mini | `hamilton_path_recovery_20261026/` (reduction to blocks, exhaustive flip-set injectivity test on all tournament classes n <= 8, decoder, stress test to n = 18 with cut vertices) + THM-4807 RESERVED "legal-path reconstruction - the Hamilton paths of a labelled tournament determine it ... settling the all-order legal re-marking guard injectivity left open by vitali-response C5" |
| 15:49 | 0946a65ace | S61 | pre-registration checkpoint |
| 16:00 | 52e3aa76e2 | S61 | note/THM-4806/HYP-9417: all-order injectivity OPEN; section 0 quotes the 15:04 README; section 9: "fetched ... again before each checkpoint ... no packet with Theorem R ... was found by the sweep" |
| 16:12 | 52823f684f | session71 | `INCOMING_RECONCILIATION.md`: S61's "global OPEN wording ... is stale relative to the independently audited Session70 all-order reconstruction theorem" |
| 16:17 | c71ab459f0 | S61 | strong form added, still as a hypothesis |

`git show 52e3aa76e2:.../session70_attachment_guard_alignment_20261026/README.md` already contains the closure text, so the superseded quotation was taken from a file that had changed in S61's own tree. Required corrections: (a) section 0 must quote the current README and cite `session70-guards/REPORT.md` as the all-order theorem; (b) THM-4806 R (iii) and the "Consequence for the carrier" paragraph must be retyped as a *reproduction/corroboration* (labelled, n <= 8) of an already-proved and audited all-order theorem, and as concurrent with mac-mini's class-level n <= 8 check (`injectivity_flips.out`: all nonempty flip sets of unused intra-block arcs, 93 sets at n = 8, zero leave HP unchanged; `verify_reconstruction.out`: decoder exact on 6,880 classes); (c) HYP-9417 (1) cannot stay OPEN (see M3); (d) the ledger bullet "all-order injectivity OPEN" and the hypotheses-INDEX wording must follow; (e) section 9's inheritance claim must be corrected (the sweep used the words "legality language"/"dark arc"; the concurrent packets say "legal guards"/"nontraceable"/"missed arcs"). Note also that the section-0 quotation is a paraphrase in quotation marks (the 15:04 text is given above verbatim); the anchor/niche/wildcard quotations are verbatim.

**M2. Concurrent mac-mini packet uncited (M1, second row).** Its counts cross-validate S61's labelled counts exactly: strong classes with a missed arc 1, 2, 9, 33 at n = 5..8 and 11 missed arcs at n = 7 correspond to S61's 120, 1,440, 38,640 labelled strong tournaments and 48,720 arcs (8 classes with one arc, one class with three: 6 x 5040 + 2 x 1680 + 1 x 5040 = 38,640, arcs 33,600 + 3 x 5040 = 48,720). Its (W) statement "a missed arc u -> v inside a strong T is the chord of a Hamilton path from u to v" is exactly the chord fact used in M3. Cite it (and THM-4807's reservation) in the note, THM-4806 and HYP-9417.

**M3. The strong form is provable; HYP-9417 (1a)/(1b) are theorems.** The implication "strong form for strong tournaments => all-order injectivity" (note, c71ab459f0) is correct; I checked each step: same components by R (ii); differing pairs inside components and dark on both sides; an arc of a component C is dark in T iff dark in C (Hamiltonian paths of T traverse C contiguously, and paths of C extend by concatenation); the avoiding path of C^(e) is a path of C' through the differing pair, extends to T', and contradicts L(T) = L(T'). What the note missed is that the strong form itself follows from its own proof of R (i):

*Chord lemma.* Let T be strong (n >= 3) and u -> v dark. With W = T - {u, v} and components W_1 -> ... -> W_m, the note's Claims give W_1 -> v and u -> W_m. If m = 1 then u -> W -> v and u -> v, so v is a sink, contradicting strongness; hence m >= 2. T strong forces an arc into W_1 from outside; it is not from W_j (j >= 2) nor from v (W_1 -> v), so u -> b for some b in W_1. Symmetrically some a in W_m has a -> v. Concatenating u, a Hamiltonian path of W_1 from b, Hamiltonian paths of W_2..W_(m-1), a Hamiltonian path of W_m ending at a, and v (Camion in each component) gives a Hamiltonian path P = (u = p_1, ..., p_n = v) of T.

*Strong form.* In T^(uv) the word (v, u, p_2, ..., p_(n-1)) is a Hamiltonian path through v -> u whose other arcs lie on P, hence are bright in T; so it avoids every dark arc of T. QED. Combined with the note's implication this is a self-contained proof of all-order injectivity (no Busch-Jacobson-Reid input), i.e. an independent proof of session70's theorem. Computational corroboration (`audit_chord.cpp`): for every dark arc of every strong labelled tournament with n <= 7 and every strong member of X_8, (a) T - {u, v} is non-strong, (b) a Hamiltonian path from u to v exists, (c) the rotated word uses no dark arc: zero violations on all 120 + 1,440 + 48,720 dark arcs of the strong labelled tournaments with n = 5, 6, 7 (none exist for n = 3, 4) and on the 2,195 dark arcs of the 2,032,504 strong members of X_8 (`chord.out`, 1m56s). The strong-form paragraph, THM-4806 R and HYP-9417 (1b) should be retyped PROVED, with the single-pair lemma demoted to a corollary.

**M4. THM-4806 P (iii) scope.** The theorem says the census "covers at most four odd prime factors otherwise"; the note says "within the enumerated bounds" without stating them, and the `.out` does not print them. Mirroring the script's loop bounds (`audit_shape_coverage.py`): stem primes come from `primes_upto(2,000,000)`, the W = 3 shapes run only when (bound/2^a)^(2/3) < 1e9, the W = 4 shape only when (bound/2^a)^(3/4) < 2e8, and odd sources (odd squares) are never enumerated for p >= 17. Coverage of even sources 2^a M', Omega(M') <= 4:

| p | 2^a q (row) | 2^a q^k, 2^a q r, 2^a q^2 r, 2^a q^3 r, 2^a q^2 r^2 | 2^a q r t | 2^a q^2 r t | 2^a q r t u |
|---|---|---|---|---|---|
| 13, 17 | complete | complete | complete | complete | complete |
| 19 | complete | complete | TRUNCATED a <= 4 (r would need to reach 2.7e7/q; primes list ends at 2e6) | TRUNCATED a = 1 | SKIPPED a = 1, TRUNCATED a = 2..4 |
| 31 | complete | TRUNCATED for a <= 20 (q^k, q r), a <= 18 (q^2 r), a <= 16 (q^2 r^2); q^3 r complete | SKIPPED a <= 17, TRUNCATED a = 18..28 | SKIPPED a <= 17, TRUNCATED a = 18..25 | SKIPPED a <= 25, TRUNCATED a = 26..28 |

So "none for p = 19" is exact for row companions (direct formula, all divisors a of 18) but "no companion with at most four odd prime factors" for p = 19 holds only outside the truncated cells and says nothing about odd-square companions; for p = 31 even Omega <= 2 shapes are incomplete for a <= 20. The theorem text, HYP-9417 (3) ("no companion with at most four odd prime factors within the enumerated bounds") and section 4 must state this table (or the script must be rerun with bounds that make the claim true; the W = 2 cells can be closed by solving the quadratic in q instead of looping). The script header also advertises a shape 2^a q^2 r t u that is not implemented, and "four odd prime factors" should say "counted with multiplicity".

## Minors

1. Transitive 4-tournament (strong-form paragraph): with the standard labelling 3 > 2 > 1 > 0 the two flip paths through 0 -> 2 are 3,0,2,1 (uses dark 3 -> 0) and 3,1,0,2 (uses dark 3 -> 1); "every flip path ... must use the dark arc 3 -> 1" is false as written (and for the labelling 2 > 1 > 3 > 0 the strong form even holds for 2 -> 0). The conclusion (strong form fails for non-strong tournaments; 48 violations at n = 4) stands.
2. "all labelled tournaments of order <= 7 (48,720 dark arcs inside strong tournaments ...; 186,480 violations in non-strong ones)": both are the n = 7 figures; the totals over n <= 7 are 50,280 strong dark arcs and 196,008 violations (0, 48, 600, 8,880, 186,480).
3. Theorem P (ii): "a = 1 gives q = N - 3 ... composite for every p <= 31" fails at p = 2 (N - 3 = 3 is prime; it is N itself since a = 1 = p - 1). Say 3 <= p <= 31.
4. "The dependent values below 180" omits 179 (dependent; present in the script's own output).
5. Section 5: "the first with an earlier column address not from its own square root is 16" - 8 = 2^3 already has the column address (2, 3); the statement is true only among squares.
6. Section 6.3 "on {1, ..., N}": the script adds edges {n, 2n}, {n, n^2}, {n, s(n)} with targets up to 3000, so the tested graph is not the induced graph on {1..N}. The thresholds are the same under the induced convention (21, 259, 25, 16; F_N and doubling+squaring planar to 3000), but the N = 21 witness quoted (branch vertices 3, 4, 8, 9, 11, 16) uses the subdivision vertex 22 = s(20); the induced K_(3,3) at N = 21 has branch vertices 1, 3, 4, 8, 9, 16.
7. Section 4: "28 is a hermit because ... the three-prime shapes start above 2N" is loose (2^a q r sources such as 30, 42 are below 54); the reason is THM-4805's exact census.
8. Section 6: THM-4470 states Skolem gracefulness of the halving forest as a citation ("Lee-Shee"), so "(CITED, PROVED there)" should read CITED (Lee-Shee via THM-4470); THM-4470's own statement is the difference-permutation bullet.
9. Section 6: the ear 11 - 25 - 24 used by all ten Hamiltonian paths of Q_25 is true (verified: 25 has degree 2, neighbours {11, 24}, never an endpoint; `audit_q25_ear.py`) but is not stated in the cited files; the "all ten" phrase occurs in session32/session33 SYNTHESIS.
10. "session70 bridge audit sec. 4.2" is `session70-bridge-audit/AUDIT-INITIAL.md` section 4.2 (the GRAPH-SUPPLEMENT has the Petersen material in its section 3); name the file.
11. Hypotheses INDEX line 2: the S60 entry's tail ("ANALOGY: the crevasse dictionary ... [note](aliquot_table...)") is stranded after the mac-mini RESERVED stub (pre-existing from the reservation order, but S61 rewrote the line and left it).
12. The `.out` of `legality_strong_lemma_20261026.cpp` was empty at audit time; the X_9/X_10 samples are NUMERICAL and use an unseeded xorshift (deterministic, fine) - state the seed and the sample counts in the `.out` when regenerated.
13. The hash in `legality_language_reconstruction_20261026.cpp` is a 64-bit commutative sum of per-path mixes; as a certificate this is sound (equal languages give equal hashes, so zero hash collisions proves injectivity). No change needed; noted for the record.

## Recomputation table

| claim (note / THM / HYP) | claimed | recomputed (auditor) | match |
|---|---|---|---|
| n = 3..7 table: tournaments (strong) | 8 (2), 64 (24), 1,024 (544), 32,768 (22,320), 2,097,152 (1,677,488) | same | yes |
| with a dark arc (strong ones) | 6 (0), 24 (0), 360 (120), 6,720 (1,440), 162,624 (38,640) | same | yes |
| dark arcs total, max | 6/1, 72/3, 1,120/6, 22,320/10, 661,248/15 | same | yes |
| pair-lemma violations | 0 at every n <= 7 | 0 | yes |
| language collisions | 0 (64-bit hash) | 0 with EXACT bitsets for n <= 6 and for all of X_7 (2^15), 0 with independent 128-bit hashes for n = 7 (2^21) and X_8 (2^21) | yes |
| X_8 row: 2,097,152 (2,032,504), dark 5,426 (1,925), 17,199 arcs | as stated | same (max dark per member 21) | yes |
| five-vertex example | 9 paths; u -> v dark; w,v,u,c,a a path after the flip | 9 paths (listed in `tournaments.out`), u -> v the only dark arc, flip path exists (13 paths after the flip) | yes |
| strong form: violations among strong, n <= 7 and X_8 | 0 | 0 (n = 3..7: 0/0/0/0/0; X_8: 0) | yes |
| strong dark arcs 48,720; non-strong violations 186,480 | (stated for "order <= 7") | n = 7: 48,720 and 186,480; per n: 0, 48, 600, 8,880, 186,480 violations | yes (n = 7 figures) |
| X_8: 2,195 dark arcs in strong members | 2,195 | 2,195 (1,790 members with one, 135 with three); 1,001 violations among non-strong members | yes |
| n = 5: 120 strong with one dark arc each; n = 6: 1,440 one each | as stated | n = 5: [1 arc: 120]; n = 6: [1 arc: 1,440]; n = 7: [1: 33,600] [3: 5,040] | yes |
| chord fact (W) and the rotated avoiding path | (auditor's lemma, M3) | every dark arc of every strong labelled tournament n <= 7 (120, 1,440, 48,720 arcs) and of every strong member of X_8 (2,195 arcs): T - {u, v} non-strong, a Hamiltonian path u -> ... -> v exists, the rotated flip path avoids all dark arcs; 0 violations in each of the three tests | yes |
| Goldbach-dependent below 10^6 / 10^7 / 10^8 | 56,267 / 462,732 / 4,117,840 | 56,267 / 462,732 / 4,117,840 (own program, loop over (m, p) pairs and even squares; untouchable sanity below 10^6: exactly {5}) | yes |
| decade shares 10^2..10^7 | 0.4378, 0.2584, 0.1525, 0.1067, 0.0903, 0.0812 | 0.43778, 0.25844, 0.15247, 0.10670, 0.09033, 0.08122 (3,655,108 of 45,000,000) | yes |
| dependent values below 180 | list ending 173 | same list plus 179 | yes (list truncated) |
| patterns at 10^7: 2 x pq 21.5% (15.9%), 2 x 3 primes 13.3% (12.7%), 2 x p 11.0% (7.0%), 4 x pq 10.9% (8.2%) | as stated | 0.2153 (0.1592), 0.1333 (0.1272), 0.1104 (0.0697), 0.1086 (0.0815) | yes |
| 2 x prime class at 10^7 / 10^8 | 0.1466 / 0.12781 (383,574 of 3,001,133) | 0.14657 / 0.12781 (383,574 of 3,001,133) | yes |
| patterns at 10^8 | 0.1990, 0.1438, 0.1006, 0.0931, 0.0697, 0.0483 | 819,451 (0.1990), 592,145 (0.1438), 414,386 (0.1006), 383,574 (0.0931), 287,041 (0.0697), 198,787 (0.0483) | yes |
| enrichment of 2 x p at 10^8 | 1.55 | 1.552 | yes |
| pre-registration before the run | 0946a65ace before `_1e8.out` | prereg committed 15:49:38 at 0946a65ace (file unchanged since); `_1e8.out` first appears at 52e3aa76e2 (16:00:23) and does not exist at 0946a65ace | yes |
| section 8 outcomes vs `_1e8.out` | 0.08122, 4,117,840, 0.12781, 0.1990, 1.55 | all present in `_1e8.out` and reproduced | yes |
| P (i): (2^a - 1) \| (N - 1) iff a \| (p - 1) | proof | 1,038 (p, a) pairs, primes p <= 61, a <= 2p + 2: 0 violations; also all integers 2 <= p <= 80: 0 | yes |
| row formula s(2^a q) | stated | a <= 11, q < 200: 0 violations | yes |
| q_a and companions | 652 (p=5,a=2); 10828 (7,2); 2^8 x 33685759 (17,8); 2^6 x 36600682668888127 (31,6); none for 2, 3, 13, 19 | same (sympy isprime on every q_a; factorisations of the composite q_a printed) | yes |
| direct s() checks | s(652)=496, s(10828)=8128, s(10682)=8128, s(44655764)=33550336, s(8623554304)=8589869056, s(2342443690808840128)=2305843008139952128 | all OK (divisor_sigma) | yes |
| a = 1: 5 \| N - 3 iff p = 3 mod 4; composite for p <= 31 | as stated | divisibility yes; "composite" fails at p = 2 (3 is prime) | partial |
| G (i) algebra, G (ii) | proof | 0 violations for p, q < 60, k <= 7; no p^2 q or p^k source of v = 2q' + 1 for primes 5 <= q' < 1000 | yes |
| A: v_2(x) + tau(g(x)) | brute force to 10^4 | own brute force to 2,000: 0 mismatches; 36: 4 two-ray vs 5 combined | yes |
| planarity thresholds | F_N planar to 3000; +doubling 21; +squaring 259; two-ray planar to 3000; Q_N 25; Q_N + doubling 16 | same under both conventions; every witness a K_(3,3) (six branch vertices of degree 3) | yes |
| Q_25: ten Hamiltonian paths through the ear | stated | 10 paths, all with 25 interior between 11 and 24 | yes |
| hypotheses INDEX 120 lines | 120 | 120 at parent, 52e3aa76e2, HEAD and working tree | yes |
| note links | resolve | all 12 link targets exist | yes |

My first Goldbach program had two bugs of my own (odd sources n = p m can exceed 2N, and the even-square loop stopped at sqrt(N)); both were corrected before the figures above, and the discrepancies they caused (e.g. 4,117,859) vanished.

## Attribution table

| item | file checked | finding |
|---|---|---|
| (i) session70 anchor / niche / wildcard | README at 1988aba492, lines 29/31/33 | verbatim |
| (i) "remain in progress" | README at 1988aba492 line 3 | paraphrase presented as a quotation; superseded by a8fcb1c840 (15:18) which is already in S61's tree at 52e3aa76e2 (M1) |
| (ii) strong-atom C5 counts | `strong_atom_continuation_vitali_20261026/README.md` line 246 | "1,1,2,8,64,1024,32768 for n=1,...,7", guard signature injective on X_n: correct; the same README (lines 290-293) now says the all-order claim is PROVED by concurrent continuations |
| (iii) HYP-9416 | title | "fresh (unnamed, multiset-counted) re-marking": correct |
| (iv) Q_n Hamiltonian paths for n in {1, 15-17, 23, 25-40} | `collatz_mod6_20260922_w6_square_sum_hamiltonicity.md` line 36 | correct |
| (iv) Gerbicz | same file, section 6 | correct |
| (iv) planar iff n <= 24 | THM-4505 line 75 -> `procgen_kuratowski_20260925_tait_kempe_triples.md` | correct |
| (iv) ear 11-25-24 used by all ten paths | cited files | true (verified) but not stated there; Minor 9 |
| (v) THM-4470 halving forest Skolem graceful | THM-4470 section 2 | stated as a citation to Lee-Shee; Minor 8 |
| (vi) bridge audit "sec. 4.2" | `session70-bridge-audit/AUDIT-INITIAL.md` 4.2; GRAPH-SUPPLEMENT section 3 | spokes -> K5 with six 2-factors, Petersen - vertex -> K3,3, exactly two of six cycles lift: correct; section reference needs the file name |
| pants_boundary_holonomy_20261011b: 2^g | line 129 | "exactly 2^g gluing classes", g = dim H^1(Gamma; F_2) (cycle rank): correct |
| boundary_support Prop. 0 ranks 0, 2, 9, 33, 271 | `boundary_support_20261020_metagraph.md` lines 103-136 | rank r - c - b, table 0, 2, 9, 33, 271 for n = 3..7: correct |
| THM-2237 sec. 3 | THM-2237 | section 3 "The missing coordinate is top parity": consistent with the row |
| FROBENIUS-ABC F1 | `strong_atom_continuation_vitali_20261026/FROBENIUS-ABC.md` lines 24-46, 158-160 | F1 phase transport; weight law misses the Frobenius degree (m = 9 example): correct |
| THM-4805 companions | THM-4805 line 68 | s^-1(6) = {6, 25}, s^-1(28) = {28}, s^-1(496) = {496, 652}, s^-1(8128) = {8128, 10682, 10828}, s^-1(33550336) = {33550336, 44655764}: correct |
| S57 counts 8,250 / 56,266 / 462,731 / 4,117,839 | `s57_untouchables_root27_20261026/README.md` lines 297-302, AUDIT.md 145 | correct |
| concurrent packets on the same question | session70 a8fcb1c840; mac-mini 9f2a9722ef (THM-4807 reserved) | NOT cited (M1, M2) |

## Prior-art findings

- **Dark arcs are "nontraceable arcs".** Arc-traceable tournaments (every arc on a Hamiltonian path) are the subject of Busch, Jacobson, Reid, *On arc-traceable tournaments*, J. Graph Theory 53 (2006) 157-166, and Busch-Jacobson, *Extremal results on arc-traceable tournaments*, JCMCC 63 (2007); they characterise a nontraceable arc x -> y of a strong tournament by a vertex z such that x lies in the first and y in the last of the (at least four) strong components of T - z. That is the mechanism of the note's five-vertex example and of "dark arcs in strong tournaments from order 5". The note does not cite BJR; session70 does (and builds its proof on it). The repo already cites Busch for the minimum number of Hamiltonian paths (MISTAKE-560 lineage).
- **Single-pair lemma.** I know no published statement in the "flip" form; it is a few lines from Camion and the component chain, and (M3) its strong form is a few more lines. By the repo standard it may be typed PROVED without a novelty claim, but it should be presented as an elementary corollary of the BJR circle of ideas, not as the new content of the packet.
- **"A tournament is determined by its set of Hamiltonian paths."** No published source known to me. Related classical facts: Redei (odd number of Hamiltonian paths; exactly one iff transitive), Camion, Moon's pancyclicity; the digraph reconstruction conjecture (vertex-deleted subtournaments) is a different notion and is false for tournaments (Stockmeyer 1977), so the word "reconstruct" should not be read in that sense. In-repo priority: session70 (15:18, proved relative to BJR, audited), then mac-mini (15:37, exhaustive classes n <= 8, decoder, stress to n = 18), then S61 (16:00, labelled n <= 8, single-pair lemma, self-contained proof now available via M3).
- **Row theorem.** s(2^a q) = (2^a - 1) q + 2^(a+1) - 1 is standard; "(2^a - 1) | (N - 1) iff a | (p - 1)" is a one-line consequence of ord(2 mod 2^a - 1) = a (2^(p-1) = 1 mod (2^a - 1) iff a | p - 1, then N - 1 = 2 - 1 - 1 = 0): folklore-level, PROVED without novelty is appropriate. The companion lists are OEIS territory (A048138 source counts, A057709 hermit numbers); 652, 10682/10828, 44655764 were already in THM-4805. The indexing of row companions by divisors of p - 1 is a neat elementary observation, plausibly unrecorded.
- **Theorem G, Proposition A, planarity thresholds**: elementary / computational; no novelty claimed, none needed.

## Commands run

All in `C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-ephrepos-math\f973fbc2-2adf-4ee3-a3b5-b87fe066eb0a\scratchpad\s61_audit\` (outputs beside the scripts):

```
g++ -O2 -std=c++17 -o audit_tournaments.exe audit_tournaments.cpp && ./audit_tournaments.exe 7 1 > tournaments.out   # 2m40s
g++ -O2 -std=c++17 -o audit_goldbach.exe audit_goldbach.cpp && ./audit_goldbach.exe 1000000 / 10000000 / 100000000  # 0.6 s / 1.6 s / 20 s
g++ -O2 -std=c++17 -o audit_chord.exe audit_chord.cpp && ./audit_chord.exe 7 1 > chord.out
python3 audit_perfect.py > perfect.out                      # sympy
python3 audit_tworay_planarity.py > tworay_planarity.out    # networkx 3.5, 3m05s
python3 audit_shape_coverage.py > shape_coverage.out
python3 audit_q25_ear.py > q25_ear.out
git log / git show / git diff / git cat-file -e on /tmp/math-wt-s61 (worktree read-only): prereg timeline, session70 and mac-mini timelines, INDEX line counts at 65a93d7747 / 52e3aa76e2 / HEAD, link resolution
```
