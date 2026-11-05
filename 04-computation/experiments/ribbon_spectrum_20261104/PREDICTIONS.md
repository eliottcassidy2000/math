# Pre-registered predictions (mac-mini-2026-11-04-ribbon), written before `ribbon_enum.c` was run

Setting. A *ribbon word of span w* is a sequence of letters a_1, ..., a_n, a_q a subset of {1, ..., w}:
vertex q has the back arcs q -> q - s for s in a_q (s < q); every other arc of the transitive base order
points forward. H = number of directed Hamiltonian paths. The profile transfer (ribbon.py, verified
against brute force on 996 random words) gives H as a product of nonnegative integer matrices on 27
profile states for w = 2 and 119 for w = 3. A word of length n >= w + 2 is strong iff every cut between
consecutive vertices is crossed by a back arc.

Established before the enumeration (python, exhaustive):
- 3-letter strong alphabet {1},{2},{1,2}, words of length <= 12: every odd number <= 3033 except 7, 21, 47
  is attained; per-length value sets contain solid odd intervals [49,165] (n=8), [49,257] (9), [129,639]
  (10), [249,1833] (11), [397,3033] (12), consecutive ones overlapping.
- 47 is attained by no span-2 word: a strong tournament with H = 47 has order <= 8 (Busch's floor f(9) = 75),
  all strong span-2 words of orders <= 9 were enumerated (probe47.py), and 47 is prime (so the ordinal-sum
  closure cannot produce it). The two strong classes with H = 47 (order 7) have ribbon width exactly 3
  (width.py), as do the four order-8 classes with H = 47.

Predictions (out of sample):
- P1. Strong span-2 ribbon words attain every odd number in [1, 10^6] except exactly 7, 21, 47.
- P2. Strong span-3 ribbon words attain every odd number in [1, 10^5] except exactly 7 and 21.
- P3. For the 3-letter span-2 alphabet and every length n in 13..18, the value set V_n contains a solid
  odd interval [c_n, d_n] with c_{n+1} <= d_n (the ribbon tiling law), with c_n < 2 c_{n-1} and d_n > 2 d_{n-1}.
Any failure is to be recorded as such.

## Correction (same session, 2026-11-04, before any enumeration result was used)

The sentence "A word of length n >= w + 2 is strong iff every cut between consecutive vertices is
crossed by a back arc" above is FALSE, and so is the implicit claim that every word over the
3-letter alphabet {1},{2},{1,2} is strong: strong components of a ribbon word need not be intervals.
Example (span 2): letters {}, {1}, {2}, {1} give the arcs 2->1, 3->1, 4->3 and all other arcs forward;
vertex 2 is a source, so the condensation is {2} < {1,3,4}, although every cut is crossed by a back
arc. The per-length value sets quoted above therefore mix strong words with ordinal sums (that is why
their minima lie below Busch's strong floor), which does not affect the value-set statements but does
affect how the enumeration must separate strong words. The enumerator now carries the condensation
automaton of condensation.py (17 states for span 2, 121 for span 3; verified against brute-force
reachability on 21,844 and 64,680 words). Predictions P1-P3 stand as written: P1/P2 are about strong
words (now correctly recognised); P3 is about all words of length n over the 3-letter alphabet.
Recorded as a MISTAKE entry in the session's close-out.
