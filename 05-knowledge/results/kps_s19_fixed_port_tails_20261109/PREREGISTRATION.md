# Pre-registration: kind-pasteur-2026-11-09-S19

Pushed before the runs named below. Exploratory results that came before this file are listed under "Already seen".

## Owner prompt (2026-11-09)

"keep going and try the tournament Chan-Pak with offset moves and the 15-matrix orbit proof"

## Concurrent work, read first, and lanes

- **Session236** (906c04d05a, 14:03 MST; session236_chanpak_offsets_orbits_20261109) is on the same prompt. It already
  has, and this session does not redo or re-claim:
  - fixed-port sink padding with a single final bridge, H = C + 2^m A;
  - the positive-offset barrier for that mechanism;
  - the one-apex zero-tail algebra;
  - the periodic-suffix drift obstruction;
  - the all-word tower-ray rigidity for the 57 ternary matrices.
- **This session's lane: multi-bridge fixed-port tails over arbitrary cores.** Take any tournament T0 with a port u.
  Append a word over the letters:
  - S: a sink, beaten by everyone;
  - I_t: t transitive vertices, then a vertex w that beats only u.

  The port stays u throughout, so every bridge returns to u.
- **What the tails add.** They contain Session236's single-bridge rays (one I letter) and the one-apex words (T0 a
  single vertex). With h_u = I(Omega(T0 - u), x) they satisfy the closed form
  h = h0 + (x/2)(F0 N(I) - h_u Q(I)), where:
  - N(I) is the binary number of the bridge positions;
  - Q(I) = sum over pairs of bridge positions p < q of c(q - p), with c(1) = 2 and c(d) = 2^(d-1) for d >= 2;
  - the closed form applies to I_0 letters after an initial sink.
- **Signed offsets.** Prepending sinks transports the offset h0 - (x/2) h_u Q(I), which can be negative.

## Already seen (exploratory, before this file)

- **Transfer check.** 800 random tails of order <= 15 match brute force in bases 2 and 3 (0 mismatches, after fixing
  a checker keying bug).
- **Cores of order <= 6, letters S, I_0, I_1, I_2:**
  - base 2: the misses below 20,000 are exactly 7 and 21;
  - base 2, below 200,000: the misses are 7, 21, 34869, 85885, 103061, 116587, 153451;
  - base 3: the misses below 20,000 are 10, 31, 5401.
- **Cores of order <= 5:** hundreds of misses below 20,000 in both bases.

## Predictions

- **Q1.** Base 2, cores of order <= 7 (all classes, all ports), letters S, I_0, I_1, I_2: the odd misses below 10^6 are
  exactly 7 and 21.
- **Q2.** Base 3, the same family: the misses v = 1 (mod 3) below 10^6 are exactly 10 and 31.
- **Q3.** Base 2, cores of order <= 6: the misses below 10^6 number between 6 and 500 (sparse but persistent), and
  every one of them lies above 30,000.
- **Q4.** Fewer cores suffice. There is a set of at most 40 core statistics (h0, F0, Delta0, h_u) among the order-<= 7
  cores whose base-2 tails already reach every odd value in [1000, 200000].

Outcomes will be appended below, unedited, after the runs.

## Outcomes

Appended 2026-11-09 after the runs. The text above is unchanged.

- **Q1: CONFIRMED.** Base 2, cores of order <= 7 (2,831 distinct core states), letters S, I_0, I_1, I_2: the odd misses
  below 10^6 are exactly 7 and 21. The search used 57,330,333 states (`out/q1_b2_le7_1e6.out`).
- **Q2: CONFIRMED.** Base 3, the same family (2,901 core states): the misses below 10^6 are exactly 10 and 31. The
  search used 34,006,533 states (`out/q2_b3_le7_1e6.out`).
- **Q3: REFUTED as worded; the intended statement is CONFIRMED.** Cores of order <= 6 leave 20 misses below 10^6, which
  is inside the predicted range 6-500. But the clause "every one of them lies above 30,000" fails as written, because
  7 and 21 are among the misses. The intended clause excluded the permanent gaps; the other 18 misses all lie above
  30,000 (smallest 34,869) (`out/q3_b2_le6_1e6.out`). This wording slip is recorded in MISTAKES.
- **Q4: UNDECIDED.**
  - **Upper bound.** A greedy cover of the odd targets in [1001, 200000] uses 81 core states (base 2); base 3 uses 86
    for [1000, 200000].
  - **Lower bound.** Every cover needs at least 17 (base 2) and 20 (base 3): the program finds that many targets whose
    covering-core sets are pairwise disjoint.
  - So "at most 40" is neither confirmed nor refuted (`out/q4_greedy_b2.out`, `out/q4_greedy_b3.out`).

## Addendum (pushed before the runs it names)

These were registered after the Q1-Q4 outcomes above, while the independent audit was running.

- **Q5.** Base 2, cores of order <= 7, letters S, I_0, I_1, I_2: the odd misses below 10^7 are exactly 7 and 21.
- **Q6.** Base 3, the same family: the misses v = 1 (mod 3) below 10^7 are exactly 10 and 31. The previous ternary
  bound was 3*10^6, the S18 audit's extension of the 15 port-chain generators.

### Addendum outcomes

(pending)
