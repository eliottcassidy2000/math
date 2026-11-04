# Verification receipt

The author ran `python verify.py` and `python -O verify.py` from this directory. Both completed successfully and produced byte-identical standard output and JSON receipts. All correctness gates are explicit exceptions, not removable assertions.

Exact controls: 75 labelled context/port pairs at context orders 0 through 3; 675 literal full tournaments with outer TT1, TT2 or C3; the three seven-coordinate base signatures; the shorter two-vertex entry signatures; strongness by reachability; endpoint and start counts by directed subset recurrence; two-slot coordinates by explicit subset partitions. No random sampling or numerical approximation is used.

A separate structural worker spot-audited the context differences, shorter TT2 signatures, exact uniform locus and two-birth strongness. Two wording clarifications were applied: the empty formal tuple is used only for h/f/j, not as a seven-coordinate identity, and a favorable isolated clock cannot establish uniformity. The theorem's unrestricted quantifiers follow from the proof and inherited exact assembly/attachment identities, not extrapolation from these finite cases. Native Hamilton-complement legality is not asserted.
