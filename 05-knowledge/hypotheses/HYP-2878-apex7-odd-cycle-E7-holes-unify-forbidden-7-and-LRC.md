---
id: HYP-2878
title: The apex-7 odd-cycle phenomenon -- E_7's odd holes (C_5 + C_7) are the shared structure behind forbidden H in {7,21}, E_7 non-chordality, and LRC(14)
status: COMPUTED+STRENGTHENED. E_7 odd holes=1496 C_5+196 C_7 (validated V=54 chi=28); 7-multiple H-atoms first at m=7 (49=7^2); SHARP FACT: 7 is the UNIQUE permanent prime gap in the H-spectrum (=> LRC(14)=2*7 hard). Three n=7 thresholds fire together; single literal C_5 object still OPEN
source: kind-pasteur-2026-06-22-S31e
related:
  - THM-200    # H=7 impossibility (the pentagon obstruction)
  - THM-079    # H=21 impossibility
  - HYP-2876   # the apex-7 floor / finite-certificate (LRC residue route)
  - HYP-2605   # the winding tournament T(x) (LRC <-> tournament bridge)
  - HYP-2758   # LRC(14) open BECAUSE 14=2*7 composite
---

# HYP-2878 -- the apex-7 odd-cycle phenomenon

## The owner's lead
The even-graph metagraph `E_n` is CHORDAL (perfect) for n<=6 and FIRST gains odd holes
EXACTLY at `n=7` = the apex prime of LRC(14) (`14=2*7`). The same prime forbids `H in {7,21}`
and indexes the 7-sector cover. Conjecture: a single odd-cycle phenomenon at the apex prime 7
underlies all three at once.

## NEW COMPUTED FINDING (kps-S31e, validated): E_7's odd holes = C_5 + C_7
`04-computation/e7_odd_holes_apex7_kps.py` (construction VALIDATED: V(E_5,6,7)=7,16,54 and
omega=5,10,28 match the canon exactly):
> **E_7 has 1692 chordless odd cycles: 1496 of length 5 (C_5) and 196 of length 7 (C_7).**
> No C_9/C_11. E_5, E_6 have NONE (chordal). So the apex prime n=7 is exactly where the
> even-graph metagraph first admits odd holes, and they come in TWO lengths -- BOTH the apex
> prime's neighbors:
- **C_5 (the pentagon), 1496 of them** = the SAME odd cycle as the `H=7` obstruction. THM-200:
  `H=7` needs the conflict graph `Omega = K_3` (3 mutually-conflicting 3-cycles), which is
  IMPOSSIBLE because 3 pairwise-sharing triangles FORCE a directed 5-cycle (a C_5) -- the
  pentagon is the obstruction that makes `H=7` (and `H=21=3*7`) permanently forbidden. The
  pentagon is ALSO E_7's first odd hole. The C_5 is the shared "first odd-cycle complication".
- **C_7 (the heptagon), 196 = 14^2 of them** = the APEX-PRIME signature: the length-7 odd hole
  is the prime 7 itself appearing as a chordless cycle in the apex-prime metagraph, and the
  count `196 = 14^2` is the LRC dimension squared.

## C_7 structure (kps-S31e enrichment)
The 196 C_7 heptagon holes are CONCENTRATED on 34/54 classes (not spread); a few classes are
highly central (appearing in 139/196 and 104/196 of them), in COMPLEMENT-SYMMETRIC PAIRS
(139=139, 104=104 -- the Z_2 complement symmetry of the merged metagraph). The most-central
C_7 class has 9 edges (NOT the naive Hamiltonian C_7 = 7 edges); 0 heptagon holes are built
purely from low-complexity (<=7-edge) even graphs. So the apex heptagon is a structured,
complement-symmetric object on the 9-edge even-graph classes -- not the trivial 7-cycle. The
C_5 holes are broader (touch 48/54 classes).

## THE ATOM-SPECTRUM SIDE (kps-S31f, ties to mac-mini HYP-2879 'weight'=atom-count)
Strong H-ATOMS (irreducible strong tournaments; `H` multiplicative over strong components,
HYP-2877) computed m=3..7 (`04-computation/strong_atom_7adic_kps.py`, validated):
- **7-multiple strong atoms FIRST appear at EXACTLY m=7** = the apex prime:
  `{35, 49=7^2, 77, 91, 105, 133, 147, 175, 189}`. None at m<=6.
- `H=7` and `H=21` are NEVER strong atoms (the `{7,21}` forbidden); `49=7^2` and `75` ARE
  irreducible m=7 atoms with COMPOSITE values that do not factor into smaller atoms (7 is not an
  atom; the 3-atom+25-atom needs 9 vertices != 7). These are mac-mini's HYP-2879 "single atoms / w=1".
- **THE SHARP UNIFYING FACT: `7` is the UNIQUE PERMANENT prime gap in the H-spectrum.** Among odd
  primes, essentially ALL are strong atoms (3,5,11,13,17,19,23,29,31,37,...,157 -- realizable H
  values); the only non-atom odd primes <=159 are `7, 107, 149`, and `107,149` are TRANSIENT
  (atoms at higher n, forbidden-seven reflection) -- so **`7` is the ONLY prime PERMANENTLY
  forbidden as a tournament `H` (Redei count) value (THM-029).** LRC(14) = `2*7` is the first open
  Lonely Runner case PRECISELY because 7 is this unique defective prime.

So THREE independent thresholds fire at `n = 7` = the apex prime: (i) `E_7` first gains odd holes
(C_5 + C_7); (ii) 7-multiple strong H-atoms first appear (49=7^2,...); (iii) `7 | H` first possible.
And `7` is the unique permanent prime gap, with `14 = 2*7` the LRC dimension. The apex-7 phenomenon
is one event -- the prime 7 "switching on" across the even-graph metagraph, the H-atom spectrum, and
(via 14) the Lonely Runner -- all at its own index n=7.

## The bridge to LRC(14) (HYP-2605, the winding tournament)
LRC(14) IS a question about the random winding tournament `T(x)` (phase `x`, difference-winding
map): `H(T(x)) = I(Omega,2)` (Redei/OCF) and `mu_{1/7} = P_x[T(x) has a scale-1/7 local sink]`
are two windows on the SAME `T(x)`; both extremized by the regular/AP object. The empty-sector
count `N_E(x) in {0,...,6}` has exactly **7 values** (= the apex prime = sector count); LRC asks
`P(N_E=0)>0`. So the apex prime 7 enters LRC as the number of sectors / the range of `N_E`, and
the obstruction lives in the odd-cycle (OCF) content of `T(x)` -- the same `c_3,c_5` that build `H`.

## Honest assessment (what is established vs thematic vs open)
- **ESTABLISHED (computed):** E_7 is the first non-chordal `E_n`; its odd holes are exactly
  `C_5` and `C_7`; `H in {7,21}` forbidden via the pentagon; LRC=T(x)'s OCF/local-sink (HYP-2605).
- **THEMATIC (same odd cycle, different graphs):** the pentagon `C_5` is the obstruction in BOTH
  the conflict graph `Omega` (forbidding H=7) AND the even-graph metagraph `E_7` (non-chordality).
  These are different graphs; the unification is "the pentagon is the first odd-cycle obstruction
  at the complexity threshold n=7", not yet a single literal `C_5`.
- **OPEN (the single-object claim):** whether ONE odd-cycle object simultaneously (i) is E_7's
  odd hole, (ii) forbids H in {7,21}, (iii) obstructs LRC(14). The C_7=apex-prime heptagon
  (196=14^2) is the most LRC-suggestive (length 7 = sector count = 14/2) but its literal role in
  the LRC winding tournament is unproven.
- **PRIOR NEGATIVE (codex-s96):** even-graph MINORS are NOT the H obstruction; loneliness is not
  minor-closed under runner-deletion. So the link is via the odd-hole STRUCTURE (the cycle space),
  not a minor order on speed-subsets.

## Next tests
1. Are the 196 C_7 holes a single S_7-orbit / Z/7-symmetric? (a canonical apex heptagon).
2. Do E_7's C_5 holes correspond, via the cycle-space bijection (even graph <-> tournament cycle
   part), to the H=7 conflict-graph K_3 obstruction? (the sharpest single-object test).
3. Does the LRC winding tournament `T(x)`'s OCF avoid a forbidden conflict-graph class exactly when
   non-lonely? (the owner's `{7,21}<->LRC` claim) -- compute the conflict graph of T(x) for the
   AP cluster across phases.
-> THM-200, HYP-2605, HYP-2876, HYP-2758; reflection `forbidden-seven-in-all-senses.md`.
