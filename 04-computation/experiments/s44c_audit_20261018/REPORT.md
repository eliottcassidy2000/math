# S44c independent audit: THM-4710, the Pythagorean/conic/RP merge note, and the post-audit edits

(This file was saved by the session from the referee subagent's final message. The subagent's file tool refused to
write report files. The text is unchanged apart from this note. All audit scripts and outputs in this directory were
written by the referee.)

- **Audited:** HEAD `db42bf9dd3`, which includes the third pre-registration outcome and THM-4710 extended to N ≤ 64.
- **Method:**
  - Everything was checked again with code written from scratch. Nothing imports the session's `vir_ss`.
  - `vir_pbw.py` is an independent Virasoro PBW engine.
  - The Witt recursion was redone exactly over Q and over Q(N).
  - Singular vectors and supports were recomputed for 5 models, including one the session never used.
  - The finite geometry was recomputed independently.
  - Literature was downloaded into the scratchpad folder `s44c_lit_untrusted/`. The decisive passages were re-read by
    the referee.
- **Verdict:**
  - THM-4710's mathematics is correct. It in fact holds for every N, not just N ≤ 64.
  - Every substantive mathematical claim in the merge note is correct.
  - There is one MAJOR finding, a missed citation in the merge note. Everything else is minor.

## Verdict

**THM-4710: correct, and stronger than claimed.**
- **The Lemma** is correct by hand and in 3,900 random-vector tests at five values of c. The top-length part never
  depends on c. It is already stated in Li 2004 Prop. 4.7 and Butson–Nair 2026 §4.3.3.
- **Uniqueness and the closed forms hold for all N.**
  - Over Q(N), the generic solution equals all 11 closed-form coefficients identically.
  - The certifying minors and all denominators have no integer roots ≥ 1.
  - N ≤ 8 was computed exactly.
- **The written proof of part 1 has a gap (T1):** the step from mod P to Q.
- **Part 2 is logically complete,** given the generator positions of HYP-9315(b). That assumption is stated in the
  status line but not in the title.
- **Actual models agree:**
  - (4,7) and (5,6), exactly over Q;
  - (6,7) and (4,13), at two primes;
  - (5,11) (N = 20, never used by the session), at two primes.
- **Novelty.** The rigidity at Li degrees 5, 7, 9 and the closed forms were not found in the literature. G_3 for every
  (3,p') is in AvEH §2.4.
- **The pre-registered test** is legitimate and honestly reported.

**Merge note: one MAJOR literature finding (M1) and several minor slips.** These were rechecked and are correct:
- the Berggren bijection and the central charges;
- Petersen, 5K3, and the K6/syntheme/total dictionary;
- the K5 contraction and the negative K3,3 result;
- Coxeter and L(Heawood);
- the degrees and the triangle rule for p ≤ 31, and Hamiltonicity;
- PSL/PGL by p mod 4 for p ≤ 31;
- the modular images A_5 and PSL(2,7);
- the RP paraphrase, which is faithful (the longest verbatim overlap with the PDF is 6 words);
- the abstract lemma and its sharpness.

**Post-audit edits.**
- The HYP-9315(d) proof is transcribed faithfully, and the weight claim in (iv) is correct.
- The corrected-window statistics reproduce bit for bit.
- "Every same-p pair is sharp" overclaims (W1).
- MISTAKE-640 is accurate.

## MAJOR

**M1. The merge note's §2 conic graphs are in Coxeter 1983, a paper the note already cites.**
- **Where:**
  - §2 "Literature" says the conic-exterior form "was not located";
  - §4 next obligation 2;
  - the ledger types the result "FINITE-EXACT (partly classical)";
  - the OPEN-QUESTIONS S44 item.
- **What is wrong:** Coxeter, *My graph*, Proc. LMS (3) 46 (1983) 117–136, page 1:
  - It defines the conjugacy graphs on the internal and external points of a conic in PG(2,q), with valencies (q±1)/2.
  - At q = 5 it names the Petersen graph and five separate triangles.
  - At q = 7 it names the edge graph of the Heawood graph and the Coxeter graph.
  - This was read in the abstract/introduction text of the OpenAlex record; the full text is paywalled.
  - The note cites Coxeter 1983, via Brouwer, only for the Fano construction. This is the repo's recurring
    cited-but-not-read failure (MISTAKE-580/629/632/640(i)).
- **Repair:**
  - Cite Coxeter 1983 for the four identifications and the valencies, and retype them as CITED, replayed independently.
  - Delete obligation 2 and the OPEN-QUESTIONS item.
  - Keep the note's own additions: the K6/syntheme dictionary, the M_x/K5 view, the negative K3,3 result, the survey,
    and the Proposition.
  - Log a MISTAKES line.

## MINOR

**T1. The mod-P to Q step in the proof of part 1 is incomplete.**
- **The problem.** The conditions involve M_(L−m), which is itself computed mod P. So the bound needs
  M^P ⊇ red(M^Q) at every earlier level. The justification given ("each previous S_(L') is one-dimensional both mod P
  and over Q") fails in two ways:
  - (a) dimension 1 over Q is shown only where a model supplies the generator; the file says exact admissibility "was
    not checked separately";
  - (b) the mod-P rank of d M + S must match the rational rank, which is not checked.
- **The conclusion is true, and the audit proves it.**
  - `witt_closure_exact.py` runs the recursion exactly over Q.
  - `witt_symbolic.py` runs it over Q(N). The generic dims are {1,0,1,1,1,1,1,1,0,1}, and all closed forms match
    identically.
  - `certificate_minors.py` gives the maximal minors, up to constants:

    | L | maximal minor |
    |---|---|
    | 3 | 384 |
    | 5 | (4N+1) |
    | 7 | (N+4)(4N+1)(4N+3) |
    | 9 | (N+5)(N+7)^5(2N+5)(4N+1)(4N+3)(4N+5)(4N²+4N−43) |

  - `denominators_check.py`: every denominator has integer roots only at −7, −5, −4.
  - So dim S_L(N) = 1 for every N ≥ 9, spanned by the closed form. N = 3..8 was computed exactly.
  - `compare_closed_forms.py` checks N = 3..20, 49..64, 70, 100, 150 and 200, exactly with no primes, and all agree.
- **Repair:**
  - State part 1 as PROVED for every N ≥ 3, and part 3 as PROVED for N ≥ L.
  - Rerun and check the three audit scripts first.

**T2.** Add "(given the HYP-9315(b) generator positions)" to the title.

**T3. The novelty line is inaccurate.**
- The first audit never reviewed literature for THM-4710; it predates it.
- AvEH §2.4, second Lemma and Cor. 2.5, give b(p') for every (3,p'). By hand, b(p') = αG_3 + β d³x_0^N, with
  α = −4(N−2)/(3(2N+1)) and β = −5α/(8N(N−1)(N−2)).
- Cite Li 2004 Prop. 4.7 and Butson–Nair 2026 §4.3.3 for the Lemma.
- Add Feigin–Nakanishi–Ooguri 1992 as unread. Its abstract emphasizes (2,q), so it is probably not prior art.

**T4.** "18 generator comparisons" is a miscount. The JSON has 15 at Li 5, 7, 9, and 28 including Li 3. This affects
THM-4710 and Li note §8.

**T5.** "S_L = 0 for N < L" covered cases that were never run. The runs used Lmax = min(9, N+3), so S_7 at N = 3 and S_9
at N = 3, 4, 5 were not computed. The audit computed them exactly, and they vanish.

**P1. The C,A branch has two limits.**
- After an A-move, t → √3 and c → 13−8√3.
- After a C-move, t → 2+√3 and c → −11. Examples: (4,15) and (2911,10864).
- Fix the table and the ledger.

**P2. The p = 7 group sentence contradicts the note's own Proposition.** The group through which the tree acts on the
Coxeter graph is PGL(2,7), of order 336 = |Aut(Coxeter)|, not PSL(2,7). The coincidence with the modular images holds
only at p = 5.

**P3. Dodecahedron and orientation wording.**
- The intended group is (A_5×±1)/±1, i.e. Aut(hemi-dodecahedron), of order 60.
- The odd elements of S_5 = PGL(2,5) ≅ SO(3,F_5) have non-square spinor norm. They swap the two hemi-dodecahedral
  embeddings of Petersen; they are not "orientation reversals".

**P4. RP margin overreach.**
- σ = 1+w and μ = 2m are a faithful reading of the source for its own case (eqs. 30–32, 40–43).
- Not established:
  - "covers stress tensors": the two-form-specific closedness step (Lemma 4.2 / App. B) is not addressed;
  - "m = 6 covers w ≤ 4": the H ≥ |P| cone and the factor e^|η| both come from the 45° reflection.
- Type "m > 1+w" as a READING. The lemma itself is PROVED.

**P5. Credit and wording.**
- clock_duality replayed only Petersen, L(Heawood) and the orders 60 and 336, not "the dictionary".
- On Hamiltonicity, say "the only known".
- The ledger's "Berggren tree = minimal models with pq even" should note that p = 1 is included and that models with
  p, q both odd are absent.

**W1. "The windows are sharp" overclaims** (Li note §3, HYP-9315(f): "every same-p pair first disagrees at 2q_min−1").
- 16 of the 33 same-p pairs have no data at that level; they are untested, not exceptions.
- The 17 pairs that reach the threshold all first disagree exactly there.
- (U) is sharp in all 87 different-p pairs.

## NOTE

- **Pre-registration timing and faithfulness are fine.** The commit was at 15:50:16 and is on origin/main, and every
  output comes later. The stable-range coordinates are exact for N ≥ 9 = Lmax.
- **Outcome nits.**
  - The first invocation crashed at parse time (sympy read `N` as `N()`) before computing anything. This is not
    mentioned.
  - "Validated against the PBW results" in fact compares against the closed forms, which the PBW runs had reproduced.
- **HYP-9315(d)(iv) is correct.** The transcription's "weight 4N−i−1" is equivalent to the first audit's "raises weight
  by 2N−i−1", and every step checks. Suggest defining "creation modes" as L_(≤−2).
- **MISTAKE-640.** Listing THM-4710 under "Survivors" is a category error: it was added after that audit.
- **Stale or loose wording in THM-4710:**
  - "(proved for 6 ≤ N ≤ 48)" on the G_5 line is stale;
  - it is the products A_ν that span, not A_j;
  - "nested" applies to the ideals, not to varieties;
  - −L/(L−1) is extrapolated from three values;
  - the positions are a reading of BGG, not a theorem.
- **The triangle rule can be typed PROVED.** A point is exterior iff its norm is a square, and the three norms of a
  self-polar triangle multiply to −1 modulo squares.
- **H_t, CGS and adjacency.**
  - Codex's H_t exists only for odd prime powers t.
  - The CGS bijection exists and is unique under both A_5 and PGL(2,5): the split attached to P is the pair of
    Stab_A5(P)-orbits on the conic.
  - Petersen adjacency is not intrinsic to K6: all pairs of splits meet as (1,2).

## Reproduction (this directory)

| script | checks |
|---|---|
| `vir_pbw.py` | the engine |
| `lemma_check.py` | the Lemma |
| `witt_closure_exact.py` | exact recursion over Q |
| `witt_symbolic.py` | recursion over Q(N) (about 40 min) |
| `certificate_minors.py` | uniqueness minors |
| `denominators_check.py` | denominators |
| `compare_closed_forms.py` | exact comparison with the closed forms |
| `model_generators_audit.py` | actual minimal models |
| `conic_audit.py` | finite geometry |
| `berggren_audit.py` | Berggren tree and branch limits |
| `modular_images.py` | modular images |
| `rp_check.py` | RP lemma sharpness and copyright overlap |
| `windows_scope_check.py` | windows sharpness scope |

Each script writes matching `.out` and `.json` files.

Literature lookup used:
- Feigin–Frenkel, hep-th/9301039;
- Feigin–Nakanishi–Ooguri, Int. J. Mod. Phys. A (1992) (World Scientific record);
- Coxeter, *My graph* (OpenAlex abstract record).
