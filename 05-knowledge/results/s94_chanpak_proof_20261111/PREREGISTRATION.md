# S94 pre-registration (opus-2026-11-11-S94), pushed before any S94 lane or census result was read

S94 is a proof-oriented session on the Chan-Pak conjecture (every odd N except 7 and 21 is the number of Hamiltonian
paths of a tournament; OPEN). Nothing below is a claim about the conjecture.

## Known to the coordinator when this file was written (in sample)

- S93 (same day): THM-5141, counting exponent 3/4 (proof-read, PASS WITH CORRECTIONS).
- Two-apex family G(2,0) (spine with apexes u -> v; S93 lane `twoapex`, files in
  `04-computation/experiments/s93_chanpak_wave_20261111/twoapex/`): with N_L = the smallest odd N > 21 that is not
  H(w) for any word w of length at most L, N_12, ..., N_16 = 230521, 776161, 2327675, 7059819, 22156807, so
  N_L / 3^(L+1) = 0.1446, 0.1623, 0.1622, 0.1640, 0.1716; and N_17 >= 2^26 (every odd N < 2^26 except 7, 21 is
  attained with length at most 17).
- H(0 w) = H(w), so the words of length exactly L have the same value set as the words of length at most L.

## Predictions (out of sample; evaluated with `coordinator/census_twoapex_bitmap.py` of this session)

- **Q1.** N_17 lies in [2^26, 72,000,000] (that is, N_17 / 3^18 <= 0.1859).
- **Q2.** N_18 lies in [192,000,000, 221,000,000] (N_18 / 3^19 in [0.165, 0.190]). In particular every odd
  N < 2^27 except 7 and 21 is attained by a two-apex word of length at most 18.
- **Q3 (a prediction about proofs).** No lane of S94 will produce a complete proof of the Chan-Pak conjecture; at
  least one lane will produce a new statement labelled PROVED that survives its adversarial verifier. Recorded so
  that the session note can say plainly whether the aim was met.

Q1 and Q2 are refuted by a value outside the stated interval. No tolerance will be changed after the fact.
