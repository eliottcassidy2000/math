---
id: THM-3266
title: "All common two-face row pairs require one origin bit"
status: >
  PROVED + FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED on the two promoted
  bank-I2 faces. Exactly 24 unordered row pairs cover both the
  support-(1,2) and support-(1,3) faces. At the single physical state (5),
  their fourteen used rows split into seven available only on the small face
  and seven available only on the full face, and every common covering pair
  crosses this bipartition. Hence none of the 24 fixed pairs admits a
  face-blind count-vector selector. With pair identity supplied, one static
  face-origin bit is necessary and sufficient for every pair. The full
  conflict atlas has union size 47 and six states hostile to all 24 pairs.
  No unrestricted 22-row, tree-complexity, dynamic-memory, or FC(3) claim is
  made.
source: root/creative-synthesis-cont/2026-08-03
audit: >
  The exact companion pins the promoted THM-3249 pair atlas, the direct
  THM-3244 full bank, and the repaired THM-3262 joint companion; reconstructs
  all small-face response vectors and full-face cover sets; enumerates all
  231 row pairs and every exclusive-label conflict; verifies one static
  face-aware selection on every pair/face/state case; and contains no
  assertion node or floating literal. Normal and
  optimized runs byte-match the frozen transcript. An independent
  reconstruction regenerated both banks without trusting the atlas
  constants, recovered all 24 pairs, the singleton bipartition, every
  conflict count, union, intersection and digest, and verified the
  face-conditioned selector on all 109,344 pair/face/nonreset cases. The
  audit also repaired the reflection so the singleton is not overread as an
  unrestricted-22-row theorem.
depends_on:
  - THM-3244-unique-reset-exposure-deletion-graph-nonmorse-boundary
  - THM-3249-cross-support-upset-atlas-local-sections-and-no-constant-gauge
related:
  - THM-3254-first-shell-two-row-clutch-and-graded-gauge-no-go
  - THM-3262-joint-two-face-selector-bit-pareto-threshold-policy-and-signed-collision-tax
script: 04-computation/fc3_common_covering_pair_face_blind_conflict_atlas_20260803.py
output: 05-knowledge/results/fc3_common_covering_pair_face_blind_conflict_atlas_20260803.out
script_sha256: c03a837f1ed2fbbadc4c9aaef8609a79b1411a1898a56a57546e5460f1fdca56
output_sha256: 40d074ea26d838bd43edea52b2cdaffeaf27a6a2cdb0dc9abedcd6156eed0e82
hash_basis: LF-normalized bytes
---

# THM-3266 -- one singleton cuts all 24 common pair charts

**PROVED + FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED on two promoted
physical faces.**

THM-3262 proves that the row pair {2,10} cannot be selected lawfully on both
promoted bank-I2 faces after the face tag is forgotten. The obstruction is
not special to that pair. It is universal across every row pair that covers
both faces, and one literal singleton state explains why.

## 1. Common covering pairs

For a face f, row r and nonreset count vector n, say r is available when its
exact response has a strict Q_f-directed one-pole ascent. An unordered pair
e={i,j} covers the face when at least one of its rows is available at every
nonreset state.

The support-(1,2) face has 52 covering pairs; the support-(1,3) face has 31.
Their intersection consists of exactly

~~~text
(2,10), (2,13), (2,19),
(3,9), (3,11), (3,16), (3,22),
(7,18), (7,22),
(10,11), (10,16), (10,22),
(11,13), (11,17),
(13,14), (13,18), (13,22),
(14,19),
(16,17), (16,21),
(17,22), (18,19), (19,22), (21,22).                   (1)
~~~

Thus there are 24 common pair charts.

## 2. The universal singleton bipartition

At the literal state (5), with padded count vector

~~~text
n=(0,0,0,0,1,0,0,0),
~~~

the available rows among all 22 charts are

~~~text
small face: {2,8,9,11,12,14,16,18,22},
full face:  {3,4,5,6,7,10,12,13,17,19,20,21}.          (2)
~~~

Row 12 is the only common available row in (2), but it belongs to none of the
24 pairs in (1). The vertices that do occur split into

~~~text
S_small={2,9,11,14,16,18,22},
S_full ={3,7,10,13,17,19,21}.                          (3)
~~~

Every pair in (1) is an edge crossing the bipartition S_small | S_full.
Consequently, at the same input (5), its S_small endpoint is exclusively
lawful on the small face and its S_full endpoint is exclusively lawful on the
full face.

It follows immediately that no function of the untagged count vector with
values in a fixed pair e from (1) can be lawful on both faces. This conclusion
is independent of the selector's circuit, formula, lookup table or query
language.

## 3. Complete conflict atlas

For a common pair e, define

~~~text
X_e={n: the small and full availability sets inside e are
         nonempty and disjoint}.                       (4)
~~~

The exact conflict counts in the order (1) are

~~~text
(19,13,16,16,14,23,11,11,8,19,8,24,
 13,14,8,16,16,12,19,18,10,17,16,8).                  (5)
~~~

Nine pairs orient small-left/full-right in increasing row order and fifteen
orient the other way. The union of the 24 conflict loci has 47 states. Their
intersection has exactly six:

~~~text
(5), (4,5), (1,3,5), (1,4,5), (3,4,5), (1,3,4,5).    (6)
~~~

For every pair, the singleton state `(5)` is the unique conflict having
pole-cardinality one. The conflict-pair multiplicity histogram is

~~~text
(number of incident pairs, number of states)=
(1,6),(2,5),(3,8),(4,3),(5,5),(6,3),
(7,6),(8,1),(9,1),(11,1),(20,2),(24,6).                (7)
~~~

The canonical complete-atlas digest is

~~~text
5faf9b8fc55cc49c83caa436d78ed98f4253e1c369e7a625ddce91cd69e61b1b.
~~~

## 4. Exact static sidecar cost

Fix one of the 24 pairs and supply its identity. Zero origin bits are
impossible by the identical input (5) with opposite exclusive labels. One
static face bit is sufficient: after reading the face, choose the first row
of the pair that is available in that face's exact cover bank.

The companion verifies this construction on

~~~text
24*(238+4318)=109344                                   (8)
~~~

pair/face/nonreset cases without one unavailable choice. Therefore the
minimum static face-origin sidecar for every fixed common pair is exactly one
bit.

## 5. Scope

Pair identity is supplied and not charged. The theorem does not optimize a
threshold tree, DAG, arithmetic formula, chart switches, routes or
history-dependent memory. The singleton (5) does not itself obstruct a
selector allowed to leave the fixed pair and use any of the 22 rows; no claim
about that larger problem is made here.

Only the two promoted bank-I2 faces are in the universe. Nothing is proved
for other supports or banks, and no Gaussian-moment positivity, FC(3), SFC(3)
or Jacobian consequence follows.

## 6. Exact reproduction

Run

~~~text
python3 04-computation/fc3_common_covering_pair_face_blind_conflict_atlas_20260803.py
python3 -O 04-computation/fc3_common_covering_pair_face_blind_conflict_atlas_20260803.py
~~~

and compare LF-normalized bytes with the declared transcript. The companion
uses exact integer/rational response data, pins every direct dependency, and
checks all covering and conflict claims exhaustively.

QED.
