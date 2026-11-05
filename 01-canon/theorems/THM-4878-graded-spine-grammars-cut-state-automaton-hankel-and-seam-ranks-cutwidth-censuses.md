---
id: THM-4878
title: "Graded spine grammars G(k, w) (a transitive chain with reversed arcs of cutwidth <= w and k <= 2 apexes): the Hamiltonian-path count is an N-rational series over an explicit finite cut-state automaton (PROVED, with the state bound (c+k+2)!(c+k+3)^(c+k)); the dictionary G(1,0) = one-apex, G(2,0) = two-apex, connected G(0,1) = upset tournaments with the all-parts-1 chain = the tribonacci tower of THM-4865, fa <= 2 in G(0,2), one-apex one-arc joins in G(2,1) (PROVED); exact Hankel ranks 4, 5, 17, 17, 18, 80, 78, 78 for G(1,0), G(0,1), G(2,0), G(1,1), G(0,2), G(0,3), G(1,2), G(2,1) (the 17 of G(2,0) reproduces Session 191) and seam ranks 3 (one reversed arc: the Lorentzian law), 13 (two arcs), 62 (three arcs) (FINITE-EXACT); the sink-run bound H >= 2^r + 1 (PROVED); the one-apex class attains 521,915 of the 1,048,576 odd numbers below 2^21 (FINITE-EXACT); pre-registered cutwidth censuses (FINITE-EXACT): the apex-free cutwidth-2 grammar G(0,2) attains every odd number <= 8191 except 7 and 21 with at most five reversed arcs per value (against the pre-registered prediction), and no family with a bounded number of reversed arcs can be complete (PROVED, counting)"
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
  AUDIT: PENDING.
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
   pending slots is at most (c + k + 2)! (c + k + 3)^(c + k). Reachable counts: G(1,0) 15, G(0,1) 14, G(2,0) 200,
   G(1,1) 114, G(0,2) 83, G(0,3) 521, G(1,2) 994, G(2,1) 2544. Hence H on G(k, w) is an N-rational series over a
   finite alphabet: pieces compose by matrix multiplication and grow by word length.

2. **Dictionary (PROVED).** G(1,0) is the one-apex class and G(2,0) the two-apex language (x_1 -> x_2 marked).
   Connected G(0,1) words are the upset tournaments; with all parts >= 2 they are Moon's interval chains (THM-4744);
   the all-parts-1 chain is the sole-exit tower of THM-4865 (H = 1, 3, 5, 9, 17, 31, 57, 105, ..., tribonacci). Every
   tournament with feedback arc number <= 2 in its transitive order (THM-4748) is a G(0,2) word. The one-arc join of
   two one-apex tournaments is a G(2,1) word, so G(2,1) attains every odd number <= 2^32 except 7 and 21 (THM-4875
   (D)). The grades are monotone in both arguments and the two axes are independent: 13 lies in G(0,2) and G(2,0)
   but in neither G(1,0) nor G(0,1).

3. **Ranks (FINITE-EXACT).** Hankel ranks (dimension of the smallest linear representation of H on the grammar):
   G(1,0) 4, G(0,1) 5, G(2,0) 17, G(1,1) 17, G(0,2) 18, G(0,3) 80, G(1,2) 78, G(2,1) 78 (the last modulo two
   primes). Seam ranks (minimal number of bilinear terms composing a prefix and a suffix across a cut with c pending
   reversed arcs): c = 1: 3 in every grammar without apexes (the Lorentzian law), c = 2: 13 (in G(0,2) and G(0,3)),
   c = 3: 62 (G(0,3)); with one apex 13 at c = 1 (G(1,1), G(1,2)) and 61 at c = 2 (G(1,2)); with two apexes 17 at
   c = 0 and 61 at c = 1 (G(2,1)). Unexplained: rank G(1,1) = rank G(2,0) and rank G(1,2) = rank G(2,1).

4. **Sink-run bound (PROVED).** r consecutive slot-free chain vertices beaten by every apex, followed by a chain
   vertex that beats some apex (or, for k = 0, lying inside a pending slot), force 2^(r-1) odd cycles through that
   vertex, so H >= 2^r + 1 (by the odd-cycle collection formula). Sinks at the end of a word and a source at its
   start do not change H. With counts saturated at 2^62, any recorded value <= X is exact.

5. **One-apex density (FINITE-EXACT).** 521,915 of the 1,048,576 odd numbers below 2^21 are one-apex counts; the
   miss fraction per dyadic scale tends to 1/2 with no bias modulo 16, 3, 5, 7.

6. **Censuses (FINITE-EXACT, pre-registered; all chain words of length <= 64, values <= X, ordered joins added).**
   The apex-free cutwidth-2 grammar G(0,2) attains every odd number <= 8191 except 7 and 21 (4094 of 4096), with at
   most five reversed arcs per value; the exact minimal-arc layers below 512 are 1, 8, 75, 133, 36, 1 values with
   0, 1, 2, 3, 4, 5 arcs (the one-arc layer is 1 + 2^k, the two-arc layer is THM-4748; the three-arc layer has 15
   pairwise-configuration types). G(1,1) and G(0,3) attain every odd number <= 511 except 7 and 21; G(2,0) every odd
   number <= 2047 except 7 and 21 (reproducing the two-apex witnesses); G(1,0) attains 560 of the 1024 odd numbers
   below 2048 directly and 763 with ordered joins; G(0,1) 181 of 500 below 1000. The pre-registered prediction P1
   (that G(0,2) would miss between 10 and 60 odd numbers below 300) is REFUTED; P3 and P5 CONFIRMED; the rest in
   HYP-9470. **Counting lemma (PROVED):** at most (log_2 X + 1)^(2c) values below X have cutwidth-2 witnesses with c
   reversed arcs, so no family with a bounded number of reversed arcs is complete and the arc count of any complete
   witness family grows at least like log N / (2 log log N).

## Validation

Exhaustive agreement of the global-identity transfer with the subset dynamic programme on all spine tournaments of
order <= 6 and random orders 9, 11, 13; the automaton engine reproduces the interval-chain, one-apex and two-apex
spectra below 1000/2048/2048 and agrees with an independent subset dynamic programme on 480 random words of six
grammars; the Hankel rank 17 of G(2,0) and the one-arc seam rank 3 reproduce Session 191 and the Lorentzian law.

## Programs

`04-computation/experiments/spine_grammar_census_20261104.cpp` (engine: census, word evaluation, automaton dump),
`spine_grammar_ranks_20261104.py` (exact ranks), `spine_grammar_ranks_mod_20261104.py` (modular ranks),
`spine_grammar_validate_20261104.py` (random-word validation), `spine_cutwidth_hamilton_20261104.cpp` (global-identity
transfer and exhaustive check), `spine_cutwidth_census_20261104.cpp` (per-word census, used for the G(0,1) control),
`spine_grammar_prereg_20261104.md` (pre-registration).
