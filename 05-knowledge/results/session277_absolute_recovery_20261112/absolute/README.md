# Why more leading-digit accuracy is not integer recovery

The new granularity theorem says that at every n≥64, **each individual
positive Burnside sector already exceeds181**. Leaving even one out prevents
nearest-integer recovery, no matter how accurately the sum predicts relative
size. The exact-engine sibling supplies the missing low-residue coordinate.

For illustration, support corrections through9 leave these errors:

| n | log10 absolute undercount | log10 relative undercount |
|---:|---:|---:|
|14|5.7521|−10.7013|
|32|66.2265|−47.6642|
|64|397.2075|−120.5655|
|100|1127.3712|−204.7573|
|250|8316.0128|−561.0363|

So at250, a relative error about10^−561 still misses an absolute amount about
10^8316. These are rounded displays of exact rational errors, not interval
endpoints. `ERRORS.csv` gives all orders0..250, not just selected successes.

`PROOF.md` proves the all-order obstruction, the exact rational support-tail
bound, and the two-coordinate interval/residue completion principle.
`certified_interval.py` exports `lower_and_tail(n,s)` and
`integer_interval(n,s)` for the low-residue engine. `CONTROLS.json` records56
exact rational interval comparisons against Session270's inherited cache.

This addresses the user's request for direction and magnitude of failure. It
does not turn a finite data table into an all-order extrapolation: the all-order
claims have separate proofs, and exact values above250 are not asserted here.
