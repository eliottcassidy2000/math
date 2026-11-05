---
id: THM-4880
title: "Ribbon tournaments (back arcs of span <= w on an ordered vertex set). (A) The Hamiltonian-path count is a finite-state series: H(a_1 ... a_n) = final . M_{a_n} ... M_{a_1} . init with nonnegative integer profile matrices on 7, 27, 119 states for w = 1, 2, 3; the minimal dimension is 4 (w = 1) and 14 (w = 2), and the span-2 series has an exact integer representation of rank 14 with entries in [-5, 5] (span 3: rank 62, entries in [-4, 4] after LLL). (B) Strong connectivity of ribbon words is regular: a condensation automaton with 4, 17, 121 states; strong components need not be intervals. (C) The value set of the span-w family is the multiplicative closure of its strong values. (D) End-twin lemma: if b's only out-neighbour is a and c -> a with every other vertex beating c, reversing the arc between b and c keeps H. (E) 47 is not a span-2 count (Busch floor + exhaustive orders <= 9); every realiser of 47 has width 3 and explicit span-3 words realise it. (F) The constant-letter runs satisfy proved linear recurrences: {2}: order 8 (growth 2.288291), {1,2}: order 12 (growth 2.562972), {1}: tribonacci. (G) The span-2 letters generate the full matrix algebra M_13(Q) on the span of nonempty-word states (irreducible). (H) No block of length <= 5 is a bounded-residual digit block for any multiplier 2..9 except ordinal-sum blocks (residual 0), and no block of length <= 4 is identity-plus-nilpotent. FINITE-EXACT: strong span-2 words attain every odd number <= 100,000 except exactly 7, 21, 47, with the solid prefix [49, d_n] of the length-<= n values growing geometrically (ribbon tiling law)."
status: >
  PROVED (mac-mini-2026-11-04-ribbon, elementary): (A) profile transfer and the exact integer quotient (exact
  arithmetic: Hankel rank, Z-module echelon/LLL bases, verified transfer on 996 + 500 random words), (B) the
  condensation automaton (merge lemma; verified against brute-force reachability on 86,778 words), (C) closure,
  (D) end twins (300 random checks), (E) via Busch's floor (THM-4744, CITED) and exhaustive enumeration, (F) Krylov
  minimal polynomials plus Berlekamp-Massey minimal recurrences verified beyond the proved orders, (G) exact integer
  computation (python-flint), (H) exact finite computation over 8,040 prefix states and 340/342 blocks.
  CITED/KNOWN: the insertion lemma H(T) <= H(T + v) (THM-1370), Busch's strong floor (THM-4744), the one-exit law
  THM-4862 (a), multiplicativity over strong components (THM-1862/THM-4094), the tribonacci tower (THM-4865 (c)),
  the one-arc join (SEAM.md (G1), THM-4875 (B)) and the near-join (HYP-9453) which appear here as the junctions of
  coupling rank 3 and 2. FINITE-EXACT: section 7 of the note (ribbon_enum.c to 5,000 with 27 coordinates;
  ribbon_levels.c to 50,000 on the 14 reduced coordinates; probe47.py orders <= 9; width.py orders <= 8). Not
  independently audited at the pre-audit commit; audit section in the note. NOT CLAIMED: Chan-Pak completeness
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
words of orders `<= 9` are enumerated; 47 is prime). All six classes with `H = 47` of orders 7-8 have width 3;
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

**(H) No-go.** For blocks `B` of length `<= 5` (three nonempty letters) or `<= 4` (with the empty letter) and multipliers
`alpha = 2, ..., 9`, the residual `H(uB) - alpha H(u)` over 8,040 prefix states is bounded only when it vanishes,
which happens exactly for the ordinal-sum blocks (`alpha = 3`: a C3 after a sink; `alpha = 5`: the strong 4-tournament
after a sink). Every block of length `<= 4` has `rank(M_B - I) >= 25` on the 26-dimensional reachable space. Hence no
bounded-carry digit system and no arithmetic-progression gadget exist among short blocks; with (G), a proof of
completeness for the family must be a collective statement about the reachable set in `Z^13` (note, section 9).

**FINITE-EXACT.** Strong span-2 words: every odd number `<= 5,000` except 7, 21, 47 (27-coordinate DFS, 6,133,659 strong
words), every odd number `<= 50,000` except 7, 21, 47 by the deduplicated level enumeration on the reduced states
(27,797,197 distinct states at length 17; every value realised at length <= 17), and every odd number `<= 100,000`
except 7, 21, 47 by the same enumeration restricted to the three nonempty letters (strong words of length `<= 16`;
per-length solid intervals `[49,165], [49,257], [129,639], [249,1833], [397,3033], [725,8609], [1281,19097],
[2715,45457], [4681, >= 99999], ...` for `n = 8..16`, consecutive ones overlapping), whose solid prefixes `[49, d_n]` of the values of length `<= n` are `d_8 ... d_16 = 165, 257, 639, 1833,
3033, 8609, 19097, 45457, 49999` (the last capped by the bound). Strong span-3 words attain every odd number `<= 5,000` except 7 and 21, all by length `<= 10` (2,417,506 reduced states at length 10). Junction coupling ranks of blocks of length `<= 4`:
1 (ordinal sums), 2 (near-join then sink), 3 (one-arc joins), 4, 8, 9, 10, 12.
