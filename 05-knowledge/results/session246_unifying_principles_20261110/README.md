# Session246: unifying physical mechanisms

**SESSION COMPLETE, 10 November 2026. Independently audited. Chan–Pak remains OPEN.**

No larger numerical coverage range is claimed. The result is four rigorous
mechanisms, with exact physical examples and explicit remaining obligations.

1. **Lawful reset transport.** An infinite physical reset family factors every
   future count, even for a strong tournament. Prime-output words must avoid
   its certified patterns. [Proof](bridges/RESET_PRINCIPLE.md), [THM-5120](../../../01-canon/theorems/THM-5120-lawful-rank-one-reset-transport-and-prime-word-avoidance.md).
2. **Physical infinite dynamics.** Bounded deletion means finitely many root
   changes. An explicit binary controller approximates every normalized ratio
   to arbitrary precision,
   yet a full-interval limit set coexists with positive-density integer gaps.
   [Proof](descent/PROOF.md), [THM-5121](../../../01-canon/theorems/THM-5121-bounded-deletion-binary-shadow-and-physical-ratio-control.md).
3. **Shared-core interval selection.** Exact interaction cancellation gives
   additive switches without a repeated-module multiplier. Negative inertia
   bounds capacity; unit-step interactions admit constructive path/cycle
   selection. [Proof](switches/PROOF.md), [THM-5122](../../../01-canon/theorems/THM-5122-shared-core-exchange-intervals-and-negative-inertia-capacity.md).
4. **Universal periodic recurrences.** A common tensor flag bounds every
   periodic fixed-apex block by an order-(k+1)(k+4)/2 recurrence. Mixed periods
   retain polynomial memory, already (4·4^m+6m+5)/9 for one apex.
   [Proof](periodic/PROOF.md), [THM-5123](../../../01-canon/theorems/THM-5123-common-gap-flag-and-periodic-apex-jordan-bounds.md).

[Unifying synthesis and exact missing step](SYNTHESIS.md) ·
[Independent hostile audit](audit/AUDIT.md).

The sampled connections are reproducible and source-aware: modular factorial
resets, Rule30 spatial resets, response tomography and local/global realization.
Inherited theorems remain credited. In particular the source grammar, the
fixed-apex compiler, quadratic ear formula, and deletion-one density gaps
are not presented as new results.

## Reproduction

From this packet directory:

    python3 -I -B verify_all.py > REPLAY.json
    cmp AGGREGATE.json REPLAY.json
    sha256sum -c SHA256SUMS

The aggregate runs four producer and three independent programs in ordinary
and optimized Python, checks stored receipts and source seals, and requires
all seven deliberate-failure controls to fail. Finite checks validate the
implementations; the linked proofs establish the all-order statements.

The all-target obligation is still an exact physical selector or a lawful
well-founded predecessor, retaining absolute count and continuation state.
Session246 has ended.
