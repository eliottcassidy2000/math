# Audit boundary and reproduction ledger

**Same-agent proof review with separate algorithms.** No independent agent,
human referee or Lean proof is claimed. The algebraic implications are given
in full in the topic documents; numerical exploration is distinguished from
universal certificate verification.

## Collatz certificate universe

At E9, the F2 observer retains m mod25, z mod5^9 and the pending first-return
flag. Its39,203,125 valid states have64 checked outgoing letters each:
2,509,000,000 inequalities. Ten possible initial all-one phases each have63
non-one entry letters, giving630 entries. Both tails beyond64 use a uniform
gross-charge bound of9 and the actual minimum/maximum certificate values.
Unused dense slots must equal the declared sentinel. All valid slots must be
finite and nonpositive. All computed successors must remain valid.

`epoch2.cpp` is the max-plus producer. `audit_epoch2.cpp` does no search and
does not use producer-generated transition tables. It enumerates the active
multiplier orbit independently, updates through successive modular halvings,
and sums the density by individual layers. Its explicit error branches work
under optimized builds; there are no C++ assertions.

The accepted binary's exact length and SHA256 are in[binary_receipts.json](binary_receipts.json).
The checker is sufficient: minimality or uniqueness of the potential is not
needed. The unsuccessful F1 extra-observer probe is a failed sufficient
certificate target, **not** a proof that no better F1 potential could exist.
No actual native excursion contradicting the desired fee was found or claimed.

The E6 run supplies a smaller positive control. Corrupting its terminal entry
is rejected by an edge inequality. The independent epoch semantics test uses
literal integer affine coefficients and a separate closed-epoch calculation:
all21,844 positive words of lengths1..7 over letters1..4, the514-letter
resonant hostile at E9, and initial phases of varied lengths. This finite
test supports the implementation; the all-length enclosure proof is(C1) in
COLLATZ.md. The initial phase orbit and final pending flag are explicitly
tested. Both normal and Python `-O` runs pass.

`lift_receipt.py` then checks the650 boundary inequalities, the exact terminal
fee, sharp all-ones clock and positive non-one slack, all nine E11 families,
and the seven-family next-layer reduction. It reads the independent certificate
receipt; it does not turn a producer diagnostic into a proof.

## Lonely runner universe

`phase_cover.py` uses the exact forward phase gap for queries. Its reference
path constructs and intersects every relevant closed comb interval. It checks
all196 primitive triple cores through12, four stated radii and64 new speeds:
50,176 queries. The point-component and closed-equality cases are explicit.
A separate pair-clock maximum calculation checks the two inherited hostile
maxima4/13 and3/13. The100- and1000-digit inputs have direct rational witnesses.

The all-speed bank for each primitive triple at radius1/5 uses a proved
positive-width tail cutoff, not an arbitrary finite stopping bound.
`lrc14_banks.py` applies that same proof to28 explicitly listed twelve-speed
cores at radius1/14. All1,458 below-cutoff cases pass both the new query and
the literal new-comb reference, and every returned witness is checked against
the actual speed vector. The tail proves acceptance of every remaining
positive integer speed for each listed core. Normal and `-O` runs agree.
The infinite arbitrary-core quantifier remains open.

## Reproduction

From this directory:

    c++ -O3 -std=c++17 epoch2.cpp -o /tmp/epoch2
    /tmp/epoch2 9 /tmp/epoch2-9.bin
    c++ -O3 -std=c++17 audit_epoch2.cpp -o /tmp/audit-epoch2
    /tmp/audit-epoch2 9 /tmp/epoch2-9.bin
    python3 epoch_semantics.py
    python3 phase_cover.py
    python3 lrc14_banks.py
    python3 lift_receipt.py

`lift_receipt.py` consumes the committed independent receipt. To audit a new
binary first, compare the auditor's JSON and SHA256 with the committed values;
do not replace a receipt merely because a producer finished. The F1 negative
control is reproducible with `multiplier.cpp 9 2` and `audit_multiplier.cpp`.
The constructor `tilted.cpp` records the prior slope-only probe; it is not a
load-bearing dependency of the E11 proof.

The checked source and receipt hashes are stored in[manifest.json](manifest.json).
External source read scope, unavailable link and non-imported proof claims are
recorded separately in[SOURCES-AND-TRANSFERS.md](SOURCES-AND-TRANSFERS.md).

## Exploratory adaptive simulation, kept outside the E11 dependency graph

`adaptive_probe.py` verifies the five finite-precision graphs listed in
ADAPTIVE.md. Its search and validation passes use different formulas for the
transition and reward, cover all states given by the proved count, and close
the infinite-letter tail. Its receipts are reproducible exact controls for
the proposed simulation. They do not imply feasibility at every precision,
and no part of the audited E11 theorem depends on that extrapolation.
