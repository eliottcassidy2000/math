# Nice pairs and weak-B/ninja pairs are different: exact onsets41 and42

Session125 terminology addendum,2026-10-31. This is a separate packet; the
original-nice-pair publication is unchanged. All new computational claims below
are **FINITE-EXACT**, generated without external code or copied witness data.
The concurrently reserved THM4850 is not a proved dependency of this check.

**INHERITED PROVED existence:** the already published
[Session124 threshold42 packet](../../session124_square_prefix_threshold42_20261031/README.md),
with its [independent audit](../../session124_square_prefix_threshold42_20261031/INDEPENDENT_AUDIT.md),
constructs a weak-B pair at every n>=42. It uses exact bases42..2081 and the
certified49-fold recursion for every n>=2082. Combining that established upper
bound with the exact lower-bound check here gives: **a weak-B pair exists if
and only if n>=42**. The present checker does not rerun those all-order
construction certificates; its new contribution is the exact onset, counts
and terminology separation.

## 1. Do not conflate the two definitions

Positions start at1. Both types use consecutive orders n,n+1 and require every
common label to occupy positions of the same parity.

| Requirement | Original Gerbicz nice pair | Weak-B pair (the ninja terminology in the current discussion) |
|---|---|---|
| Beginning |1|1,3|
| End for even order |8|8|
| End for odd order |3|10|
| Necessary common port colors |1,3 odd;8 even|1,10 odd;3,8 even|
| First pair |41/42|42/43|
| Pairs at first order |6|20|
| Common phase masks at first order |3|6|

The two classes are disjoint at every order: label3 must have odd position in
an original nice pair and position2 in a weak-B pair. The words "ninja" and
"nice" therefore cannot be silently interchanged.

Closing an original pair member gives a square-sum cycle because1+3=4 and
1+8=9. A weak-B **even** member also closes by1--8, but its odd member does not
close by its endpoints:1+10=11 is not square. This distinction matters when
transporting cycle claims. This packet proves only the stated finite onset,
not a new all-order theorem or any reserved generalization.

## 2. Independent weak-B onset check

Let c(x)=1 for an odd position and0 for an even position. A candidate with
prefix1,3 and normalized end8/10 is **eligible** for pairing only if
c(8)=0 and c(10)=1. The conditions c(1)=1,c(3)=0 are already forced by the
prefix. "Raw" below means only prefix and endpoint are imposed; "eligible"
also imposes the two necessary common port colors.

|Order|Raw paths|Eligible paths|Eligible phase masks|
|---:|---:|---:|---:|
|32|0|0|0|
|33|0|0|0|
|34|7|1|1|
|35|0|0|0|
|36|4|0|0|
|37|0|0|0|
|38|2|2|2|
|39|0|0|0|
|40|11|0|0|
|41|1|0|0|
|42|139|49|39|
|43|136|74|30|

Every pair with smaller order32..41 has an empty eligible side. Thus no
intersection of nonempty spectra has been skipped. At42/43 the shared spectra
have exactly six masks, with shorter/longer multiplicities

```
common odd-position mask   shorter paths   longer paths
1585793692427              1               2
1599483269979              1               1
1602705117979              2               1
3372531545867              1               1
3389445068555              4               3
3526867243787              1               2
```

Bit x-1 represents label x. The total number of normalized weak-B pairs is
1*2+1*1+2*1+1*1+4*3+1*2=20. Reversal is not counted separately because the
prefix fixes the orientation.

For n<=30, the even member would give a Hamiltonian cycle of even order<=30,
contradicting the elementary low-degree obstruction. For n31 or32, the even
member would be a G32 cycle containing edge1--3. The independently enumerated
unique G32 cycle does not contain that edge. Thus42/43 is the first weak-B
pair among *all* positive orders, not merely the first in the scanned window.

## 3. The final obstruction is one label, not nonexistence of a raw path

There is exactly one raw order41 path, independently recovered by both
algorithms:

```
1,3,6,30,19,17,32,4,21,15,34,2,7,18,31,33,16,9,40,41,8,
28,36,13,23,26,38,11,5,20,29,35,14,22,27,37,12,24,25,39,10.
```

It puts8 at position21. But the even-order member must end at8 in an even
position. Consequently label8 alone prevents pairing41 with42 (or40).
For these two nonempty raw normalized path sets, this one-label projection
is minimal: projection onto zero coordinates cannot separate two nonempty sets.
It is incorrect to strengthen "no compatible weak-B member at41" into "no
path with prefix1,3 and endpoint10 at41."

Similarly, all raw even-order paths at36 and40 put10 in an even position,
whereas the odd member must end10 in an odd position. These fail the necessary
port coloring before any deeper phase-spectrum test.

Contrast the original nice-pair obstruction at40/41: eligible paths exist on
both sides and four internal coordinates are needed for a separating projection.
That four-label theorem remains in the unchanged original packet. The two
onsets have different failure mechanisms, not just different numerical labels.

## 4. Positive first-pair witness

One weak-B pair independently generated by the degree solver is

```
P42 = 1,3,33,31,18,7,29,35,14,11,38,26,10,6,30,19,17,32,4,5,20,
      16,9,40,24,25,39,42,22,27,37,12,13,36,28,21,15,34,2,23,41,8
P43 = 1,3,22,42,39,25,24,40,9,27,37,12,13,36,28,8,41,23,2,34,15,
      21,43,6,30,19,17,32,4,5,20,16,33,31,18,7,29,35,14,11,38,26,10
```

Both have prefix1,3, correct endpoints, square adjacent sums and matching
position parity for every label1..42. The common odd-position mask is
1585793692427. New label43 occurs in odd position23, as slot counting requires.
The certificate contains a witness for each of the six shared masks.

## 5. Reproduction and exact scope

The self-contained standard-library checker uses two different exhaustive
searches: path growth with necessary residual-degree/connectivity pruning, and
edge-subset degree propagation with connectedness tested only after completion.
They construct square edges independently and agree on every eligible mask
and its multiplicity for each order32..43. The unique raw41 witness is also
independently recomputed with the path port filter disabled. G32 uniqueness
and the absence of its edge1--3 are checked again by the edge-subset solver.
Raw counts elsewhere are exact outputs of the degree solver; the independently
compared quantities are the eligible spectra and multiplicities.

From this directory, or using absolute script paths from any working directory:

```sh
python3 verify_weak_b_onset.py --output /tmp/weak-b-normal.json
python3 -O verify_weak_b_onset.py --output /tmp/weak-b-optimized.json
cmp /tmp/weak-b-normal.json /tmp/weak-b-optimized.json
cmp /tmp/weak-b-normal.json weak_b_certificate.json
sha256sum -c SHA256SUMS
```

The checker rejects duplicate labels, wrong prefix/endpoints and a nonsquare
edge, and explicitly tests both the port distinction and the new-label parity.
It has no external data dependency and does not import the original packet.
Its constants n32..43 specify a finite universe, not an extrapolation. The
claimed first onset combines that finite classification with the stated
small-order cycle obstruction.
