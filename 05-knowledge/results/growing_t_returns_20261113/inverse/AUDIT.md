# Inverse-lane audit and source pins

**VERIFIED controls; PROVED IN NOTE mathematical results; target return OPEN.**
This lane uses the session base `9b8aa35188`. The parent session owns Git
integration and publication. This child made no Git mutations.

## Reviewed proof scope

The root agent independently reviewed sections 1--5 of
`UNIT_READY_MINIMUM.md`: the strong score-prefix argument, the OCF
two-cover lower bound, separation of all non-singleton cores, exact
scalar minimization, ties, and the resulting order bound. The independent
geometry agent also read those sections completely and accepted all
these arguments, including the canonical equality-type classification.
Neither review substitutes for the stated inherited grammar or OCF
dependencies.

After those reviews, section 5 was sharpened from source orders21..28
to21..25 by evaluating its already-reviewed inequality at the target.
The root separately checked and requested that exact table: nineteen-
vertex cores allow n0..10; twenty allow n0..4; twenty-one allow n0..1;
twenty-two and twenty-three allow n0 only. This is twenty parameter
cells, not an enumeration of the cores in those cells. The other final
change was the explicit source link for the inherited suffix identity.

The final reviewed proof snapshot is:

    UNIT_READY_MINIMUM.md
    SHA256 666062186917059e90b0860ab70138f6afd5a8bb41e2d8f3845497319bca13ab

The root also independently reviewed and accepted the cyclic
correlation formula, sharp bound and equality classification, and the
reindexed-clock factorial/parity argument. Its requested final
precisions were incorporated: the minimum drift is of order2^(-L/2)
with different even/odd limiting constants; the generic repeated-word
formula begins at j1, since the first prefix can have a parity boundary
correction. The positive physical example has compatible prefix parity
and retains its all-j formula. The clock note explicitly keeps both the
marked root and predecessor outside the uniform module in its eventual
readiness argument. Their final proof snapshots are:

    CYCLIC_DRIFT.md
    SHA256 45fa46b1fc46a2e8790403f242d46bb4f74da59098d5c1ab0f4eb8b99f4de444

    REINDEXED_CLOCKS.md
    SHA256 04f885b3543e5f41d501f2cb5601473f179ff0d50042d481c6eddf8dc0a68bf2

`SOURCE_PINS.json` records these files, all four checkers and eight
receipts, and the inherited proof/data dependencies. It records hashes,
not a claim that the inherited census was freshly generated.

## Executed controls

Each command below was executed both normally and with `-O`; all four
normal/optimized receipt pairs are byte-identical. Runtime validity
checks use explicit exceptions, so optimized execution retains them.

```powershell
python -B 05-knowledge/results/growing_t_returns_20261113/inverse/unit_ready_minimum.py
python -O -B 05-knowledge/results/growing_t_returns_20261113/inverse/unit_ready_minimum.py
python -B 05-knowledge/results/growing_t_returns_20261113/inverse/cyclic_drift.py
python -O -B 05-knowledge/results/growing_t_returns_20261113/inverse/cyclic_drift.py
python -B 05-knowledge/results/growing_t_returns_20261113/inverse/reindexed_clocks.py
python -O -B 05-knowledge/results/growing_t_returns_20261113/inverse/reindexed_clocks.py
python -B 05-knowledge/results/growing_t_returns_20261113/inverse/joint_probe.py --order 18
python -O -B 05-knowledge/results/growing_t_returns_20261113/inverse/joint_probe.py --order 18
```

The minimum checker covers every singleton address of bitlength at most19,
all6008 strong order-eight representatives of the inherited6880-class
census, all196417 unique-start parent classes through order15, and
literal extremal graphs through source order15. The physical controls
check full H, strongness, unique predecessor starts and unit readiness;
the order-seven alternatives have distinct degree multisets.

The cyclic checker covers all32766 words of lengths1..14 and four
literal growing examples. The clock checker compares the inherited
complete ready fibres and independently recounts54 actual source/clone
graphs. The targeted joint probe excludes every2,178,309 canonical
inner core of order18 with arbitrary outer-tail length.

These finite universes support the proofs at their stated scope. They
do not prove completeness of the whole metagraph, produce a unit-ready
k17 witness, or establish a uniform target-preserving405 return. Chan--Pak
remains OPEN.
