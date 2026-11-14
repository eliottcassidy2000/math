# Independent S328 phase-companion audit: PASS

Reviewed `MINIMUM_PHASE_SWITCH.md` at SHA-256
`cf1eee0f24039c2936a29569e7aa6450651aa2b5905e4d458dcd9a9ce2bf5627`;
 `bounded_catalog.cpp` at
`49eadc86d685f3bb0e61ee327007d217a7e98a0ac7e668bf8134a38878f0f6b9`;
and `replay_switch.py` at
`95405595d95750f4e31613bacc3b95f83f6905a444aa16f90aafc98fad1f5bdc`.

- The stabilizer proof correctly reduces odd k in [1,2a−1] by the reflection
  k↦2a−k. A strict odd witness in (2a/k,4a/k) exists as asserted; for k≥3
  it lies below 2a. The two stabilizers and the two flips exhaust the
  unit-action possibilities even for composite odd slopes.
- For n≥2, Q₂ is twice a sum of 2^(n−1) odd products, one per pair of
  complementary subsets; both Rédei oddness and H(empty)=1 are needed.
  Hence Q₂ is divisible by 4. The unique singleton at Q₂=2 cannot supply
  a distinct companion. No fixed-H phase exclusion was made for H=1.
- The deletion injection recovers both the old cover and the selected
  slot, since the slot containing the inserted vertex identifies the
  choice. Successive deletions to the empty tournament or to a directed
  triangle establish the advertised lower bounds. Induction over budget
  and order therefore exhausts the census, with no assumption of
  strong connectivity. Degree-bucket relabeling is a complete canonical
  comparison. Both path recurrence directions and complementary
  subset clock aggregation are sound; n≤12 is comfortably within
  signed 64-bit bounds.
- Independently compiled with `g++ -O2 -std=c++17` and regenerated both
  literal catalogue and log byte-for-byte. Independently ran the Python
  replay, also matching its complete receipt byte-for-byte. This reports
  3,578 classes, 919 distinct profiles, and the pairs (19,217) and
  (25,223) as the only companions at the global first clock Q₂=4096.
- The target's digit is forced modulo 198, and the canonical positive
  digit interval has one representative per odd class. A phase
  contributes only an integer dyadic scale, still a multiple of 4
  whenever H>1. Thus the forbidden shift 2cd excludes the entire
  fixed-other-H 99-digit chart, not all tournament realizations and
  not all other chart sources. The lower target cutoff 405523 concerns
  **complete complementary equal-clock first-phase envelopes**, not
  individual repairs. The two optimal pairs yield distinct targets
  even under independent dyadic retiming, as shown by coprimality and
  the strict leading-term versus offset inequality.

No mathematical corrections required. The result is not an eventual
Chan–Pak threshold or a classification at later phase levels.
