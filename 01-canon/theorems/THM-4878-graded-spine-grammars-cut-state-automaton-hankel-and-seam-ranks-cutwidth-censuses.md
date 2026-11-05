---
id: THM-4878
title: "Graded spine grammars G(k, w) (a transitive chain with reversed arcs of cutwidth <= w and k <= 2 apexes): the Hamiltonian-path count is an N-rational series over an explicit finite cut-state automaton (the method is KNOWN: frontier dynamic programming, DEFECT_TRANSFER D1-D2, THM-4880; the explicit apex-and-slot automaton with its state bound is PROVED); the dictionary G(1,0) = one-apex, G(2,0) = two-apex (containing THM-4875's family F1), connected G(0,1) = upset tournaments with the all-spans-1 chain = the tribonacci tower (THM-4865, THM-4731 (6)), fa <= 2 in G(0,2), ribbon words of span s in G(0, 2s-1) (PROVED); exact Hankel ranks 4, 4, 17, 17, 17, 79, 78, 78 for G(1,0), G(0,1), G(2,0), G(1,1), G(0,2), G(0,3), G(1,2), G(2,1) with the empty tournament counting one path (the 17 of G(2,0) reproduces Session 191; the three two-radix grades share rank 17) and seam ranks 3 (one reversed arc: the Lorentzian law), 13 (two arcs), 62 (three), 310 (four) (FINITE-EXACT); the sink-run bound H >= 2^r + 1 (PROVED); the one-apex class attains 521,915 of the 1,048,576 odd numbers below 2^21 (FINITE-EXACT); pre-registered cutwidth censuses (FINITE-EXACT): the apex-free cutwidth-2 grammar G(0,2) attains every odd number <= 8191 except 7 and 21 with at most five reversed arcs per value (against the pre-registered prediction), and no family with a bounded number of reversed arcs can be complete (PROVED, counting)"
status: >
  PROVED: (1) the cut-state transfer theorem and the state bound; (2) the dictionary of grades; (4) the sink-run lower
  bound, the padding lemma and the saturation soundness. KNOWN: the monotonicity of H under induced subtournaments
  (Redei insertion); Moon's product for interval chains (THM-4744 (1)); the tribonacci law (THM-4865); the
  Lorentzian one-arc join law (session29 ADDENDUM, THM-4866 (B), THM-4875 (B)); the odd-cycle collection formula
  (THM-4876 (4), PROP-001); the two-arc formulas (THM-4748); the one-apex four-coordinate transfer (Session 16,
  HYP-9468); the two-apex Hankel rank 17 (Session 191, reproduced independently here).
  FINITE-EXACT: (3) the Hankel and seam ranks (exact rational arithmetic; the modular script agrees on every grammar it
  was run on); (5) the one-apex density below 2^21; (6) the censuses (see the note, section 6).
  NOT CLAIMED: that any grade is onto the odd numbers other than 7 and 21 (HYP-9470; OPEN).
  CONCURRENT: mac-mini-2026-11-04-ribbon (THM-4880, HYP-9473): the bounded-span analogue, committed before this
  session's first checkpoint and missed at first filing (MISTAKE-720).
  AUDIT: PASS WITH CORRECTIONS (independent, adversarial, own code: 4 majors, 18 minors, 5 upgrades, all applied;
  04-computation/experiments/spine_grammar_audit_20261104/AUDIT.md). Majors: prior art (ribbon lane, DEFECT_TRANSFER,
  THM-4876 (6), THM-4731 (6)); the false 2^32 inference for G(2,1) from THM-4875 (D); the empty-word convention of
  the apex-free Hankel ranks; the counting lemma's proof.
session: opus-2026-11-04-S79
source: 05-knowledge/results/spine_grammars_cutwidth_apexes_hankel_seams_20261104.md
related:
  - THM-4744
  - THM-4748
  - THM-4865
  - THM-4866
  - THM-4875
  - THM-4876
  - THM-1370
  - THM-4749
  - HYP-9470
  - HYP-9468
  - 05-knowledge/results/feedback_word_hamilton_spectrum_20261020/README.md
  - 05-knowledge/results/square_sum_spectrum_overnight_20261030/README.md
  - session194/two-apex/PROOF.md
  - session195/combined-growth/PROOF.md
---

# THM-4878 - graded spine grammars: cut-state automaton, Hankel and seam ranks, cutwidth censuses

## Setting

A **spine tournament** is a transitive chain v_1 < ... < v_m with a set of chain arcs reversed; its **cutwidth** is the
largest number of reversed arcs crossing one cut between consecutive chain vertices. **G(k, w)** is the class of
tournaments with k <= 2 apexes x_1 -> x_2 (each beating an arbitrary subset of the chain) over a spine tournament of
cutwidth <= w. Read left to right, each chain vertex is a **letter** (T, b, A): the pending reversed arcs ("slots")
whose tail it is, the number of slots it opens (their head), and the apexes that beat it; slots opened at one vertex
form a head group and a letter terminates at most one slot per group. A **cut state** is the multiset of fragments of
a Hamiltonian path restricted to the prefix, each fragment labelled by how it starts (global start, a pending slot, an
apex) and how it ends (global end, or open with the set of pending slots and apexes whose head is its last vertex).
Full definitions, proofs, programs and tables: the source note.

## Statements

1. **Transfer (PROVED).** For every (k, w) there are nonnegative integer matrices M_ell (one per letter) and
   P_1..P_k on the finite set of cut states, with H(T(ell_1...ell_m)) = e_0 M_(ell_1) ... M_(ell_m) P_1 ... P_k f.
   The entries count the local predecessor/successor choices of the new vertex. The number of cut states with c
   pending slots (fragments with distinct start labels among the 1 + c + k labels, the global end on at most one
   fragment, each slot bit in at most one mask, apex bits free, and the head-group pattern) is at most
   2^(max(c-1,0)) sum_j C(1+c+k, j) (j+1)^(c+1) 2^(kj), which for c + k >= 1 and k <= 2 is below (c + k + 2)!
   (c + k + 3)^(c + k) (the audit's count; the crude bound fails only at c = k = 0, where three states are reachable).
   Reachable counts: G(1,0) 15, G(0,1) 14, G(2,0) 200, G(1,1) 114, G(0,2) 83, G(0,3) 521, G(1,2) 994, G(2,1) 2544,
   G(0,4) 3384. Hence H on G(k, w) is an N-rational series over a finite alphabet: pieces compose by matrix
   multiplication and grow by word length. The method is KNOWN (frontier dynamic programming over a linear layout of
   bounded cutwidth: Knuth TAOCP 4A section 7.1.4; in the repository the component-pairing transfer with the
   backward-arc cutwidth of a supplied order, `overnight_rank_return_20261029/DEFECT_TRANSFER.md` D1-D2, and the
   bounded-span profile transfer of THM-4880); what is proved here is the explicit apex-and-slot automaton with its
   letter alphabet and its validation.

2. **Dictionary (PROVED, elementary).** G(1,0) is the one-apex class and G(2,0) the two-apex language (x_1 -> x_2
   marked). Connected G(0,1) words are the upset tournaments; with all spans >= 2 they are Moon's interval chains
   (THM-4744); the all-spans-1 chain is the sole-exit tower of THM-4865 (H = 1, 3, 5, 9, 17, 31, 57, 105, ...,
   tribonacci; the identification is THM-4731 (6) and the ribbon note's {1}^n, not new). Ribbon words of span <= s
   (THM-4880) are G(0, 2s - 1) words. Every
   tournament with feedback arc number <= 2 in its transitive order (THM-4748) is a G(0,2) word. THM-4875's family
   F1 (one-apex x one-apex joins with the reversed arc between the apexes) lies in G(2,0) after relabelling the
   apexes, and F1 alone has positive miss density (S78 note: 65 misses in [2^10, 2^11), 934,107 in [2^31, 2^32));
   the join of two one-apex tournaments across a reversed arc between chain vertices is a G(2,1) word whose coverage
   has not been computed (the audit corrected an earlier false inference of 2^32 coverage for G(2,1) from THM-4875
   (D), whose second family F2 is not a spine word). The grades are monotone in both arguments and the two axes are
   independent: 13 lies in G(0,2) and G(2,0) but in neither G(1,0) nor G(0,1).

3. **Ranks (FINITE-EXACT; conventions: letter alphabet = the letter indices of the current slot count, the (T, b, A)
   alphabet gives the same ranks; the empty tournament counts one path, so the empty cut is accepting).** Hankel
   ranks (dimension of the smallest linear representation of H on the grammar): G(1,0) 4, G(0,1) 4, G(2,0) 17,
   G(1,1) 17, G(0,2) 17, G(0,3) 79, G(1,2) 78, G(2,1) 78, G(0,4) 389 (the last two modulo two primes).
   With the engine's convention (empty word = 0) the apex-free grades gain one dimension: 5, 18, 80, 390 (audit M3).
   Seam ranks (minimal number of bilinear terms composing a prefix and a suffix across a cut with c pending reversed
   arcs): c = 0: 1 (the ordered-join product; 2 with the engine convention); c = 1: 3 in every grammar without apexes
   (the Lorentzian law; THM-4880's junction coupling rank 3), c = 2: 13 (in G(0,2), G(0,3), G(0,4); confirmed from
   below by the audit's brute-force Hankel block), c = 3: 62 (G(0,3), G(0,4)), c = 4: 310 (G(0,4), modular); with one
   apex 4 / 13 / 61 (G(1,1), G(1,2)); with two apexes 17 at c = 0 and 61 at c = 1 (G(2,1)). The three two-radix
   grades G(0,2), G(1,1), G(2,0) share Hankel rank 17.

4. **Sink-run bound (PROVED).** r consecutive slot-free chain vertices beaten by every apex, followed by a chain
   vertex that beats some apex (or, for k = 0, lying inside a pending slot), force 2^(r-1) odd cycles through that
   vertex, so H >= 2^r + 1 (by the odd-cycle collection formula); if no such vertex follows, every vertex of the run
   is an ordinal-sum cut and the count equals that of the word without the run. Sinks at the end of a word and a
   slot-free source at its start do not change H. With counts saturated at 2^62, any recorded value <= X is exact.

5. **One-apex density (FINITE-EXACT; extends THM-4876 (6) and HYP-9469 (3), which record the closed formula and
   "about half of the odd numbers" with 50,006 of the 100,000 below 200,000).** 521,915 of the 1,048,576 odd numbers
   below 2^21 are one-apex counts; the miss fraction per dyadic scale tends to 1/2; it is flat modulo 16 and 5
   (0.501-0.504) but biased modulo 7 (0.455, 0.488, 0.551, 0.581, 0.533, 0.466, 0.443 for the residues 0..6) and
   mildly modulo 3 (0.501, 0.490, 0.516) (audit).

6. **Censuses (FINITE-EXACT, pre-registered; all chain words of length <= 64, values <= X, ordered joins added).**
   The apex-free cutwidth-2 grammar G(0,2) attains every odd number <= 8191 except 7 and 21 (4094 of 4096), with at
   most five reversed arcs per value; the exact minimal-arc layers below 512 are 1, 8, 75, 133, 36, 1 values with
   0, 1, 2, 3, 4, 5 arcs (the one-arc layer is 1 + 2^k, the two-arc layer is THM-4748; the three-arc layer has 15
   pairwise-configuration types; the exactness of the minima uses the length cap 64, justified because a word with
   at most 4 arcs and H <= 511 has at most 9 interior gaps with 1 + 2 sum 2^(r_i - 1) <= 511, hence length <= 58).
   G(1,1) attains every odd number <= 2047 except 7 and 21 with at most one reversed chain arc beside the apex, and
   G(0,3) every odd number <= 511; G(2,0) every odd number <= 2047 except 7 and 21 (reproducing the two-apex
   witnesses); G(1,0) attains 560 of the 1024 odd numbers below 2048 directly and 763 with ordered joins; G(0,1) 358
   of 1024 below 2048 and 181 of 500 below 1000. Attainment is certified by explicit witnesses; the miss statements
   rest on the length cap, the sink-run cut, the closer bound and saturation (all sound), and were reproduced by the
   audit for G(1,0) and G(0,1). The pre-registered prediction P1
   (that G(0,2) would miss between 10 and 60 odd numbers below 300) is REFUTED; P3 and P5 CONFIRMED; the rest in
   HYP-9470. **Counting lemma (PROVED; proof corrected by the audit):** in a connected cutwidth-2 word with c reversed
   arcs the vertices between consecutive slot events are slot-free and lie inside a pending span, so by the sink-run
   bound each of the at most 2c + 1 interior gaps has length at most log_2(X - 1) when H <= X; with at most 2c event
   vertices arranged in one of C(c) ways, at most C(c) (log_2 X + 2)^(2c + 1) values below X have witnesses with c
   reversed arcs. Hence no family with a bounded number of reversed arcs is complete, and the arc count of any
   complete witness family grows at least like (1 - o(1)) log_2 N / (2 log_2 log_2 N).

## Validation

Exhaustive agreement of the global-identity transfer with the subset dynamic programme on all spine tournaments of
order <= 7 (2,097,152 at order 7) and random orders 9, 11, 13; the automaton engine reproduces the interval-chain,
one-apex and two-apex spectra below 1000/2048/2048 and agrees with an independent subset dynamic programme on 480
random words of six grammars; the Hankel rank 17 of G(2,0) and the one-arc seam rank 3 reproduce Session 191 and the
Lorentzian law. The independent audit reproduced with its own code the censuses to 511 (all grades) and to 2047
(G(0,2), G(2,0)), the one-apex and interval-chain spectra, the one-apex 2^21 count, all listed witnesses, the ranks
from the dumps, and the seam ranks 3 and 13 from brute-force Hankel blocks.

## Programs

`04-computation/experiments/spine_grammar_census_20261104.cpp` (engine: census, word evaluation, automaton dump),
`spine_grammar_ranks_20261104.py` (exact ranks), `spine_grammar_ranks_mod_20261104.py` (modular ranks),
`spine_grammar_validate_20261104.py` (random-word validation), `spine_cutwidth_hamilton_20261104.cpp` (global-identity
transfer and exhaustive check), `spine_cutwidth_census_20261104.cpp` (per-word census, used for the G(0,1) control),
`spine_grammar_prereg_20261104.md` (pre-registration).
