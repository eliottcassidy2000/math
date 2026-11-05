---
id: THM-4880
title: "Ribbon tournaments (back arcs of span <= w on an ordered vertex set). (A) The Hamiltonian-path count is a finite-state series: H(a_1 ... a_n) = final . M_{a_n} ... M_{a_1} . init with nonnegative integer profile matrices on 7, 27, 119 states for w = 1, 2, 3; the minimal dimension is 4 (w = 1) and 14 (w = 2), and the span-2 series has an exact integer representation of rank 14 with entries in [-5, 5] (span 3: rank 62, entries in [-4, 4] after LLL). (B) Strong connectivity of ribbon words is regular: a condensation automaton with 4, 17, 121 states; strong components need not be intervals. (C) The value set of the span-w family is the multiplicative closure of its strong values. (D) End-twin lemma: if b's only out-neighbour is a and c -> a with every other vertex beating c, reversing the arc between b and c keeps H. (E) 47 is not a span-2 count (Busch floor + exhaustive orders <= 9); every realiser of 47 has width 3 and explicit span-3 words realise it. (F) The constant-letter runs satisfy proved linear recurrences: {2}: order 8 (growth 2.288291), {1,2}: order 12 (growth 2.562972), {1}: tribonacci. (G) The span-2 letters generate the full matrix algebra M_13(Q) on the span of nonempty-word states (irreducible). (H) For blocks of length <= 5 and multipliers 2..9 the residual H(uB) - alpha H(u) is constant on the 8,040 prefix states only for the ordinal-sum blocks (then 0) and takes at least 227 values otherwise; on the 14-dimensional quotient the eigenvalue 1 of every block matrix of length <= 4 is at most simple, so no block is unipotent on the reachable space and none inserts a non-constant arithmetic progression. FINITE-EXACT: strong span-2 words attain every odd number <= 100,000 except exactly 7, 21, 47, with the solid prefix [49, d_n] of the length-<= n values growing geometrically (ribbon tiling law)."
status: >
  KNOWN (method; audit correction 4): the profile transfer of (A) is the standard frontier-based / bounded-bandwidth
  transfer-matrix dynamic programme for Hamiltonian paths (Knuth, TAOCP 4A, section 7.1.4, simpath; frontier-based
  search), and the regularity of (B) follows abstractly from the MSO-definability of strong connectivity
  (Buechi-Elgot-Trakhtenbrot); new here are the tournament specialisation (free tails, at most one unfed head,
  7/27/119 states), the explicit merge-lemma automaton, the integer quotient, (G), (H) and the census.
  PROVED (mac-mini-2026-11-04-ribbon, elementary): (A) profile transfer and the exact integer quotient (exact
  arithmetic: Hankel rank, Z-module echelon/LLL bases, verified transfer on 996 + 500 random words), (B) the
  condensation automaton (merge lemma; verified against brute-force reachability on 86,778 words), (C) closure,
  (D) end twins (300 random checks), (E) via Busch's floor (THM-4744, CITED) and exhaustive enumeration, (F) Krylov
  minimal polynomials plus Berlekamp-Massey minimal recurrences verified beyond the proved orders, (G) exact integer
  computation (python-flint), (H) exact finite computation over 8,040 prefix states, 583 blocks x 8 multipliers, and
  the exact characteristic polynomials of the 340 block matrices of length <= 4 on the 14-dimensional quotient (the
  unipotency half was supplied by the audit: the session's rank bound rank(M_B - I) >= 25 did not exclude a single
  nilpotent Jordan block, MISTAKE-719).
  CITED/KNOWN: the insertion lemma H(T) <= H(T + v) (THM-1370), Busch's strong floor (THM-4744), the one-exit law
  THM-4862 (a), multiplicativity over strong components (THM-1862/THM-4094), the tribonacci tower (THM-4865 (c)),
  the one-arc join (SEAM.md (G1), THM-4875 (B)) and the near-join (HYP-9453) which appear here as the junctions of
  coupling rank 3 and 2. FINITE-EXACT: section 7 of the note (ribbon_enum.c to 5,000 with 27 coordinates;
  ribbon_levels.c to 50,000 on the 14 reduced coordinates; probe47.py orders <= 9; width.py orders <= 8). AUDIT: PASS WITH CORRECTIONS (independent, adversarial; 05-knowledge/results/ribbon_transfer_audit_20261104.md):
  one major (the (H2) proof gap, repaired as above) and twelve minors, all applied on 2026-11-04 (file name and index
  wording: the 14-dimensional representation is reducible, irreducibility holds on the 13-dimensional nonempty-word
  space; (A), (B) KNOWN as methods; (H1) stated as the computed fact; the general-w ordinal-sum rule; components
  inside A are singletons for w = 2; width-table entries 0 at orders 3-4; the all-order width argument for 47; DFS
  first-encounter lengths; level-receipt caveats). The audit reproduced every proof and computation with its own
  code, including an independent undeduplicated DFS to 50,000 (misses exactly 7, 21, 47). NOT CLAIMED: Chan-Pak completeness
  (OPEN; HYP-9473 records the ribbon tiling law and the span-3 completeness conjecture that would imply it).
session: mac-mini-2026-11-04-ribbon
source: 05-knowledge/results/ribbon_tournaments_bounded_span_transfer_20261104.md
scripts:
  - 04-computation/experiments/ribbon_spectrum_20261104/ribbon.py
  - 04-computation/experiments/ribbon_spectrum_20261104/condensation.py
  - 04-computation/experiments/ribbon_spectrum_20261104/intrep.py
  - 04-computation/experiments/ribbon_spectrum_20261104/lll_span3.py
  - 04-computation/experiments/ribbon_spectrum_20261104/minrep.py
  - 04-computation/experiments/ribbon_spectrum_20261104/exact_algebra.py
  - 04-computation/experiments/ribbon_spectrum_20261104/recurrences.py
  - 04-computation/experiments/ribbon_spectrum_20261104/minrec.py
  - 04-computation/experiments/ribbon_spectrum_20261104/probe47.py
  - 04-computation/experiments/ribbon_spectrum_20261104/width.py
  - 04-computation/experiments/ribbon_spectrum_20261104/witness47.py
  - 04-computation/experiments/ribbon_spectrum_20261104/digitscan2.py
  - 04-computation/experiments/ribbon_spectrum_20261104/junctions.py
  - 04-computation/experiments/ribbon_spectrum_20261104/ribbon_enum.c (test_s2_misses.txt)
  - 04-computation/experiments/ribbon_spectrum_20261104/ribbon_levels.c (lv2_5e4_levels.txt)
related:
  - THM-1370
  - THM-4094
  - THM-4744
  - THM-4862
  - THM-4865
  - THM-4875
  - THM-4876
  - HYP-9029
  - HYP-9453
  - HYP-9473
  - 05-knowledge/results/seven_return_native_recursions_20261029/SEAM.md
---

# THM-4880 - ribbon tournaments: a finite-state transfer for Hamiltonian paths

**Setting.** A ribbon word of span `w` is a sequence of letters `a_1, ..., a_n`, `a_q` a subset of `{1, ..., w}`; its
tournament has the vertices `1, ..., n`, the back arcs `q -> q - s` (`s in a_q`, `s < q`) and all other arcs forward.
Span-2 letters are written `--`, `1-`, `-2`, `12`. `H` counts directed Hamiltonian paths. Every tournament is a ribbon
word for `w` large enough (its ribbon width is the least such `w`).

**(A) Profile transfer.** Processing the vertices in order, a Hamiltonian path restricted to the processed set is a
set of fragments whose heads (other than the global start) lie in the window of the last `w` vertices and whose tails
outside the window are free. The profile state (fragments as (head age or unfed, tail age or free)) is finite, the
extension rules are local, and `H = final . M_{a_n} ... M_{a_1} . init` with nonnegative integer matrices on 7, 27,
119 reachable states (`w = 1, 2, 3`). The Hankel rank is 4 (`w = 1`) and 14 (`w = 2`). The Z-module generated by the
functionals `u -> H(uv)` has rank 14 for span 2 (generated by suffixes of length `<= 3`), is letter-invariant, and
yields an integer representation `y(u) in Z^14`, `y(empty) = e_1`, `H = y_2`, letters acting by integer matrices with
entries in `[-5, 5]`; for span 3 the rank is 62 and an LLL-reduced basis gives entries in `[-4, 4]`. Fourteen short
suffixes (lengths `<= 3`, listed in the note) form a Q-basis of the observables: the future of a prefix is determined by
its counts under those fourteen extensions.

**(B) Strong connectivity is regular.** Adding vertex `v` with past out-neighbourhood `A` (inside the window) merges with
`v` an initial stretch of the components from the first component meeting `A` onward, leaving after it only components
contained in `A`; a component without window vertices that precedes the first window-containing component is final.
The resulting automaton has 4, 17, 121 states for `w = 1, 2, 3`. Strong components of a ribbon word need not be
intervals (`--.1-.-2.1-`: condensation `{2} < {1, 3, 4}`).

**(C) Closure.** Induced subtournaments and ordinal sums of ribbon words of span `w` are ribbon words of span `<= w`, so
the value set of the family is the multiplicative closure of its strong values (and its primes are its strong primes).

**(D) End twins.** For `X` with a vertex `b` whose only out-neighbour is `a`, adding `c` with `c -> a` and all other
vertices beating `c`, the orientations `b -> c` and `c -> b` give equal `H` (proof in the note, section 5).

**(E) 47.** No span-2 ribbon word has `H = 47` (a strong realiser has order `<= 8` by Busch's floor; the strong span-2
words of orders `<= 9` are enumerated; 47 is prime). Every tournament with `H = 47` has width exactly 3: its strong realisers are the two order-7 classes (the four order-8
classes with `H = 47` are ordinal sums of one of them with a vertex), both of width 3, and an ordinal sum with transitive
pieces keeps the width (at least 3 by the induced-subtournament bound of (C), at most 3 by its closure);
`---.1--.---.-2-.--3.12-.-23` and `---.---.---.-23.1--.-2-.--3` are span-3 words with `H = 47`. Width distribution of
all classes of orders 3-8: note, section 6.

**(F) Runs.** `H({2}^n)`: `1, 1, 3, 5, 11, 29, 65, 149, 343, 777, 1781, 4081, 9329, 21361, 48881, ...`, for `n >= 9`
`H(n) = 2H(n-1) - H(n-2) + 4H(n-3) - 3H(n-4) + 6H(n-5) - 2H(n-6) + 4H(n-7) - H(n-8)`, growth 2.288291;
`H({1,2}^n)`: `1, 1, 1, 5, 15, 41, 105, 261, 659, 1689, 4347, 11165, 28623, 73329, 187875, ...`, order-12 recurrence
with characteristic polynomial `x^12 - 2x^11 - 2x^9 - 2x^8 - 4x^7 - 4x^6 - 4x^5 + x^4 + 2x^3 + 3x^2 + 2x + 1`, growth
2.562972; `H({1}^n)` is tribonacci (THM-4865 (c)). Neither new sequence is in the OEIS.

**(G) Irreducibility.** On the 13-dimensional span of the states of nonempty span-2 words, the four letters - and already
the three nonempty letters - generate the full matrix algebra `M_13(Q)`: there is no common invariant subspace and no
lower-dimensional quotient of the series.

**(H) No-go.** For the 583 blocks `B` of length `<= 5` (three nonempty letters) or `<= 4` (with the empty letter) and the
multipliers `alpha = 2, ..., 9`, the residual `H(uB) - alpha H(u)` over the 8,040 prefix states (7,873 distinct reduced
states) is constant for exactly 22 of the 4,664 pairs, all ordinal-sum blocks with `H(B) = alpha` (then the residual is 0:
`alpha = 3`, a C3 after a sink; `alpha = 5`, the strong 4-tournament after a sink), and takes at least 227 distinct values
for every other pair. On the 14-dimensional quotient the eigenvalue 1 of the block matrix `A_B` has algebraic
multiplicity 0 for 289 and 1 for 51 of the 340 blocks of length `<= 4`, never 2 or more (exact characteristic
polynomials, audit); since the reachable space maps onto `Q^14`, no block is identity-plus-nilpotent on it, and a
non-constant arithmetic progression `H(u B^k) = H(u) + kc` would need a Jordan block of size at least 2 at the
eigenvalue 1. (`rank(M_B - I) >= 25` on the 26-dimensional reachable space is a datum, not a proof: a single nilpotent
Jordan block of size 26 has rank 25.) Hence no bounded-carry digit system and no arithmetic-progression gadget exist
among short blocks; with (G), a proof of completeness for the family must be a collective statement about the
reachable set in `Z^13` (note, section 9).

**FINITE-EXACT.** Strong span-2 words: every odd number `<= 5,000` except 7, 21, 47 (27-coordinate DFS, 6,133,659 strong
words; the audit's own undeduplicated DFS to 50,000 agrees), every odd number `<= 50,000` except 7, 21, 47 by the deduplicated level enumeration on the reduced states
(27,797,197 distinct states at length 17; every value realised at length <= 17), and every odd number `<= 100,000`
except 7, 21, 47 by the same enumeration restricted to the three nonempty letters (strong words of length `<= 16`;
per-length solid intervals `[49,165], [49,257], [129,639], [249,1833], [397,3033], [725,8609], [1281,19097],
[2715,45457], [4681, >= 99999], ...` for `n = 8..16`, consecutive ones overlapping), whose solid prefixes `[49, d_n]` of the values of length `<= n` are `d_8 ... d_16 = 165, 257, 639, 1833,
3033, 8609, 19097, 45457, 49999` (the last capped by the bound). Strong span-3 words attain every odd number `<= 5,000` except 7 and 21, all by length `<= 10` (2,417,506 reduced states at length 10). Junction coupling ranks of blocks of length `<= 4`:
1 (ordinal sums), 2 (near-join then sink), 3 (one-arc joins), 4, 8, 9, 10, 12.
