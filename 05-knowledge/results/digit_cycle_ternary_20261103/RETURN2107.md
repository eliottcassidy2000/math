# A retained-root return at 2107 with no old deletions

**FINITE-EXACT literal construction.** The eight-arc minimum below is
relative to a specified target bank, not a global native-move minimum.

Continue from the exact source in
[the preceding handoff](../return701_algebra_20261103/HANDOFF.md). Its names
are `[1,2,4,5,6,7,8,13,14,21,22]`; root 21 has the sole in-neighbour 13,
and spare 6 is retained. Its count is 2107=49*43, and its golden chart
is (1,4,3).

The inherited complete one/two-start grammars were scanned through parent
order 10. At that order they supply two 2-ready targets, both of order 11,
with rooted profiles

`(2107,151,153,2109)` and `(2107,395,397,2109)`.

The scan found no smaller target in those two grammars. It is a witness
search, not an all-size classification of 2107-ready sources: larger
parents and the exceptional-family classification are not needed here.

For these two targets, all `2*9!=725760` identifications fixing both the
actual root and its unique old in-port were compared. The minimum arc
discrepancy from the actual source is eight. The saved minimum witness
has a repair order such that every completed arc reversal is strong and
admits the two-complement helper guard proved in the previous packet.

The resulting literal return uses:

* 16 native complements, each with its supplied actual Hamilton path;
* eight temporary helper births/deletions;
* **zero old-vertex deletions**, retaining every old name;
* peak order 12, and returned order 11.

Its actual old root 21 has profile `(2107,151,153,2109)`, so it is
2-ready. The returned golden observer is `(2107,4824,4860)` and its chart
is `(0,2,5)`. This is not a neutral-phase return. The spare's five Q4
presentations and the full eighteen-polynomial two-port observer are
retained in [return2107.json](return2107.json).

The ordinary two-birth reset at that root now produces

`6325=3*2107+4=25*11*23`.

The next root is name 31, the two new names are 31 and 32, and the spare
remains name 6. The next source has order 13, profile
`(6325,2107,2107,6325)` and golden chart `(6,3,1)`.
All rows and source names are saved. This additional step does not prove
uniform closure under either or both ternary offsets.

The target comparisons count labelled identifications of two specified
rooted classes. An independent endpoint DP recounts the returned profile.
The original 81/27 pair is separately reconstructed in the group/port
scripts; neither its order 27 nor the new digit clock is used as a
replacement for an actual orientation guard.
