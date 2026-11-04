# Independent audit of the complete positive joint guard

**PASS — 2026-11-04.** The coordinator accepts
[JOINT_CLASSIFICATION.md](JOINT_CLASSIFICATION.md), including its extension
to all finite tournaments, uniqueness corollary, and complete strong-input
zero-joint table. This classifies a guard, not return or spectrum closure.

The image/complement proof keeps the two free slots of the K term ordered.
The first free slot is the segment cut before a, and only its endpoint
determines whether the join can be restored. The exact surplus is therefore
[E_B W_w A_a]+[E_B P Q_aw], with no missing factor of2. Coincident marks
make the L term zero and Q_aa the singleton a; empty free slots are retained.

The shortest a-to-w path contains only its first exit from B. Every omitted
exit gives a distinct nonnegative K witness, paying d−1 (or d for a=w).
The all-tournament upgrade is essential and valid: when a cannot reach w,
the scalar failed-join path ending in B contains w in its interior. Cutting
after w gives an injective L-surplus cover. Hence tau>=sigma in that case.
Together with sigma=[E_B A_a], this also deals with degree zero. The remaining
degree-one formula gives the claimed iff, without importing strongness.

Two different positive auxiliaries at the same root cannot both have that
root as sole exit: their mutual arc is a second exit for one. For the SCC
sidecar, a spanning two-ending-path cover must put a terminal-SCC vertex in
a path ending in that SCC. Conversely, cut a Hamilton path ending at either
terminal endpoint after the other endpoint. Strongness supplies a path
ending at every specified vertex (equivalently cut a Hamilton cycle). In a
strong T with a sole exit b, b reaches all of T−a and belongs to its initial
SCC. The proof's four-row zero-joint table follows, including the reducible
deletion case with sigma0 and the actual transported-zero hostile. The
positivity qualification is necessary and has not been dropped.

The coordinator read all three independent sources and replayed the exact
controls. The core DFS source ran in normal and optimized Python modes;
both pass **110,927 gates** and match its frozen receipt. The permutation
SCC and weak-input sidecars were independently replayed and match their
frozen receipts at **26,851** and **47,221** gates respectively; the producer
also ran those sidecars in both modes. Their declared universes, reachable
and unreachable cases, coincident marks, nonempty endpoint slots, positive
digits and hostile zero cases were checked against the source loops.

|Standalone source|Receipt|Gates|Receipt SHA256|
|---|---|---:|---|
|[joint guard](../../../../04-computation/nonzero_preparation_20261103/algebra_joint_guard.py)|[JSON](../../../../04-computation/nonzero_preparation_20261103/algebra_joint_guard.json)|110927|148a19cf1fd9bc54a75576782fc700030f3d53b27dcf3febb8c66e289224396a|
|[SCC sidecar](../../../../04-computation/nonzero_preparation_20261103/algebra_joint_scc.py)|[JSON](../../../../04-computation/nonzero_preparation_20261103/algebra_joint_scc.json)|26851|d136b727d13ef78b877e613a9b20a8d02f778cc3cf60dfc6bf6ce44870e40052|
|[weak-input sidecar](../../../../04-computation/nonzero_preparation_20261103/algebra_joint_weak.py)|[JSON](../../../../04-computation/nonzero_preparation_20261103/algebra_joint_weak.json)|47221|3103157b1b6214be3326a83309af6d8d1ac9d1d8aa7c8ff24c60ccf80da19a7f|

Source pins and complete replay commands are in the primary note. The core
contains26705 marked cases through order5; the SCC sidecar checks10650
unordered endpoint pairs and14002 strong marked pairs. The weak sidecar
checks every literal image in the nonreachability injection (3966 images
over3212 unreachable marked pairs) and all312 positive weak guards.
These finite universes test the proof and are not its all-order premise.

The audited prepromotion proof SHA256 is
`49ed5b683d2d7c3900a100b6cc07775966d405d9a2e0a29df66ab4041e7e11c1`.
Status/link-only promotion is permitted. No Lean formalization or external
priority claim is made.
