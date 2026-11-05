# Pre-registered predictions: graded spine grammars G(k, w) (opus-2026-11-04-S79)

Written and committed before any census of G(0,2), G(0,3), G(1,1), G(1,2) was run. The engine is
`spine_grammar_census_20261104.cpp` (canonical cut-state automaton over a transitive chain with backward arcs of cutwidth
<= w plus k <= 2 apexes placed after the chain; values <= X over all chain lengths, pruned by the insertion lemma and
deduplicated on count vectors; ordered-join products added afterwards). "Miss" means an odd number <= X that is not the
Hamiltonian-path count of any tournament in the grammar or an ordered join of such tournaments.

Validation targets (not predictions; observed by other code before this file was written):

- V1. G(0,1) attains exactly 181 of the 500 odd numbers <= 1000, first misses 7, 11, 13, 19, 21, 23
  (per-word engine `spine_cutwidth_census_20261104.cpp`; the family is Moon's interval chains, THM-4744).
- V2. G(1,0) attains exactly 560 of the 1024 odd numbers <= 2047, first misses 7, 13, 21, 29, 35, 41, 43
  (one-apex automaton, `oneapex.cpp` run in this session; the one-apex class omits 13 at every order: Session 16).
- V3. G(2,0) attains every odd number <= 999 except 7 and 21 (two-apex witness file
  `square_sum_spectrum_overnight_20261030/two_feedback_100000.json`).

Predictions:

- P1. G(0,2): besides 7 and 21, between 10 and 60 of the 150 odd numbers <= 300 are missed (attained fraction between
  0.60 and 0.93), and the miss fraction does not tend to zero: at least 3 percent of the odd numbers in [2^11, 2^12) are
  missed.
- P2. G(0,3): at most 5 odd numbers <= 300 are missed besides 7 and 21, and at least one odd number <= 2000 is missed
  besides 7 and 21 (no fixed cutwidth without apexes is complete; a structural guess).
- P3. G(1,1): 13 is attained, and at most 10 percent of the odd numbers <= 2000 are missed besides 7 and 21.
- P4. G(1,2): no odd number <= 2000 is missed besides 7 and 21.
- P5. Minimal grades: 13 has minimal grade (0,2) and (2,0) (incomparable), i.e. 13 is in G(0,2) and in G(2,0) but not in
  G(1,0) or G(0,1) (the G(0,2) membership follows from THM-4748's crossing formula 2^(r+2) + 2*3^r - 1 at r = 1).
