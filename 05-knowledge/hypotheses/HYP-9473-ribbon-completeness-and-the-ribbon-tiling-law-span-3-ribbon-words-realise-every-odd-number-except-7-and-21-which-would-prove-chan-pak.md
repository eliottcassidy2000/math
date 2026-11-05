---
id: HYP-9473
title: "Ribbon completeness and the ribbon tiling law (after THM-4880). (1) Every odd number other than 7 and 21 is the Hamiltonian-path count of a strong span-3 ribbon word; this would prove the Chan-Pak conjecture (THM-4094 lanes p != 7, 7p). (2) The span-2 family attains every odd number except exactly 7, 21, 47 (47 is proved absent). (3) Ribbon tiling law: the solid prefix [49, d_n] of the span-2 values of words of length <= n has d_n -> infinity, with d_n growing geometrically (observed factors 1.6-2.9; d_8..d_15 = 165, 257, 639, 1833, 3033, 8609, 19097, 45457); for the three-letter alphabet alone the per-length solid intervals [c_n, d_n] overlap consecutively. (4) The minimal dimension of the span-3 series is 62 and its representation is irreducible. (5) What a proof needs: a collective invariant of the reachable set in Z^13 (no quotient, no digit block, no arithmetic progression is available, THM-4880 (G), (H))."
status: >
  OPEN (mac-mini-2026-11-04-ribbon). (1) NUMERICAL: 47 realised by explicit span-3 words; strong span-3 values to the
  bound of the note's section 7. (2) FINITE-EXACT to 50,000 (level enumeration on the 14 reduced coordinates) and to
  5,000 by an independent 27-coordinate DFS; 47's absence PROVED (THM-4880 (E)). (3) FINITE-EXACT for n <= 16 (all
  letters) and n <= 12 (three letters); pre-registered as P3 in PREDICTIONS.md for n = 13..18 (status in the note).
  (4) rank 62 of the observable module PROVED; the Hankel rank and irreducibility for span 3 are open (minrep.py's
  span-3 run is recorded in minrep_w3.log if it finished). (5) is a reformulation, not progress on Chan-Pak, which
  remains OPEN (MISTAKE-559).
session: mac-mini-2026-11-04-ribbon
date: 2026-11-04
source: 05-knowledge/results/ribbon_tournaments_bounded_span_transfer_20261104.md
related: [THM-4880, THM-4094, THM-1370, THM-4876, THM-4865, THM-4875, HYP-9029, HYP-9453, HYP-9469]
scripts: 04-computation/experiments/ribbon_spectrum_20261104/ribbon_levels.c; ribbon_enum.c; PREDICTIONS.md
---

# HYP-9473 - ribbon completeness and the ribbon tiling law

**Setting (THM-4880).** Ribbon words of span `w`: back arcs of span `<= w` on an ordered vertex set. `H` is a
finite-state series (14 reduced integer coordinates for span 2, 62 for span 3), strong connectivity is regular (17 /
121 states), and the family's value set is the multiplicative closure of its strong values.

**(1) Span-3 completeness => Chan-Pak.** By THM-4094 the completeness of the odd spectrum is equivalent to strong
realisation of every prime `p != 7` and of `7p`; by THM-4880 (C) the strong span-3 words suffice if they attain every
odd number other than 7 and 21. Evidence: the span-2 strong words attain every odd `<= 50,000` except 7, 21, 47, and 47
is a span-3 word. So the conjecture is a statement about one explicit automaton: a 62-dimensional integer
representation over 8 letters together with a 121-state strong-connectivity automaton.

**(2) Span 2 misses exactly 7, 21, 47.** 47 is absent (proved); every other odd number `<= 50,000` is attained by a strong
span-2 word; the only odd numbers in `[1, 50000]` not attained by words of length `<= 16` are listed in the note.

**(3) Ribbon tiling law.** Let `V_n` be the set of counts of span-2 words of length `<= n` (nested: the empty letter keeps
`H`) and `d_n` the largest odd number with `[49, d_n] subset V_n`. Then `d_n -> infinity`; observed
`d_8, ..., d_15 = 165, 257, 639, 1833, 3033, 8609, 19097, 45457`. For the three-letter alphabet (no empty letter) the
per-length value sets are not nested and their largest solid intervals `[c_n, d_n]` overlap consecutively
(`c_{n+1} <= d_n`; `[49,165], [49,257], [129,639], [249,1833], [397,3033]` for `n = 8..12`). Either form, with the
finite base, gives (2) and, with span 3, (1). This is HYP-9029's tiling law restricted to a family where the
"strong(n)" sets are the images of a finite automaton.

**(4) Span 3.** The observable Z-module has rank 62 (PROVED, lll_span3.py); whether the Hankel rank is 62 and whether the
representation is irreducible is open.

**(5) What a proof needs.** THM-4880 (G) and (H) exclude quotient numerations, bounded-carry digit blocks and arithmetic
progressions. The per-`h` extension sets are sparse (extdensity.py), so the covering is collective: a proof must exhibit
an invariant of the reachable set `Y_n subset Z^13` (for instance a lattice box in two coordinates, or a sumset
structure coming from the two-sided pairing `H(uv) = <chi(v), y(u)>`) that propagates under the letters and whose image
under the `H`-functional of a fixed extension covers a window of odd numbers. The ordered-path reduction of HYP-9469 (6)
and the B-set density of HYP-9453 (2) are the unbounded-span versions of the same missing lemma.
