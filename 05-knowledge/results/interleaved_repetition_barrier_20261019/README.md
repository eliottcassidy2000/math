# A repetition barrier for realized Collatz first crossings

19 October 2026. **PROVED, scoped, with a CITED explicit logarithmic-form input.**
The first-crossing bank language uses K,V,F,H from the preceding work and
starts with a nonempty expansion run over {K,V}.
General convergence, arbitrary bank entry, and the full source-sensitive
2-adic boundary remain open.

## Read first

[The interleaved repetition theorem](INTERLEAVED_REPETITION_BARRIER.md) is the
primary result. A fully lawful expansion-first first-crossing word has strict endpoint descent
if it contains a consecutive expansion-bank power U^k whose duplicated bits
(k-1)|U| comprise at least 1% of all its expansion bits. Episode count and word
length are unrestricted. Repetition can occur after earlier resets.

The proof retains the original source through the preceding context, uses
exact repetition fuel to force exponential source height, and compares that
with a polynomial ceiling for any non-descending endpoint. A new exact
short-range theorem supplies the overlap; it applies to all episode counts.

Additional results:

- The normalized carry obeys e/S<N/6 at every expansion-first first crossing.
- Every expansion-first first-crossing bank word of Terras length at most
  2^52 descends,
  regardless of its episode count.
- Actual fixed-displacement endpoint witnesses must avoid the stated
  macroscopic expansion repetitions. Approximations to K-infinity still
  require exponentially long completions in their initial repetition depth.
- The family K^L(FK)^(j-1)F has unbounded expansion episodes and is covered
  for every L. A separate exact control has 1,557 episodes.

[The power-expansion note](POWER_EXPANSION_ZERO_HEAD.md) preserves the simpler
single-episode discovery route. Its all-length endpoint theorem is a
corollary of the stronger interleaved result. Prior one-episode repair already
excluded least counterexamples from that subclass; the new significance is
strict endpoint descent and, in the main note, unrestricted interleaving.

## Verification

Run from this directory:

    python3 verify_interleaved_repetition.py
    python3 -O verify_interleaved_repetition.py
    python3 verify_power_expansion.py
    python3 -O verify_power_expansion.py

Both pairs agree byte for byte with their respective receipts. The primary
checker makes 255,312 explicit checks and replays 228,395 Terras steps. Its
complete depth-eight prefix forest has 523 first-crossing leaves; the 532
uncrossed frontier nodes are explicitly retained. The large separate family
control uses exact independently reconstructed literal coefficients and the
full-word congruence, not a claim of full parity replay at that size.

The finite controls are not the universal quantifiers. Matveev's theorem is
an explicitly cited external input, with its application and constants
checked; no numerical test is offered as a proof of that theorem.

## Independent audit

The [independent proof review](review/INTERLEAVED_REPETITION_REVIEW.md) passes the frozen expansion-first theorem with its cited Matveev input. It also documents the resolved reset-only scope issue. Run `python3 -B review/independent_interleaved_controls.py` and repeat with `python3 -B -O review/independent_interleaved_controls.py`; compare with the included independent receipts. Both runs make 1,225,098 exact checks, including a complete depth-nine prefix forest and the 1,557-episode literal replay. Source and published hashes are recorded in `SOURCE_PROVENANCE.json`.
