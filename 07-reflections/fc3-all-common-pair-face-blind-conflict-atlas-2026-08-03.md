# FC3 common-pair conflict atlas: one singleton cuts all 24 charts

**Status:** FINITE-EXACT discovery package; theorem-strength on the two
promoted `I2` faces, but not yet a canonical theorem.

**Scope:** the 24 unordered row pairs that cover both the support-`(1,2)` and
support-`(1,3)` physical faces of THM-3249.  Pair identity is treated as
given.  The result concerns deterministic selectors from padded count vectors
and a static face tag; it does not optimize a decision tree, analyze dynamic
memory, allow a selector to leave its chosen pair, or prove `FC(3)`.

Exact companions:

~~~text
04-computation/fc3_common_covering_pair_face_blind_conflict_atlas_20260803.py
05-knowledge/results/fc3_common_covering_pair_face_blind_conflict_atlas_20260803.out
~~~

## 1. Inheritance pass

The closest proved mechanism is
[THM-3249](../01-canon/theorems/THM-3249-cross-support-upset-atlas-local-sections-and-no-constant-gauge.md):
the two faces have respectively 52 and 31 covering row pairs, with exactly 24
in common.  The hostile example inherited from the repaired row-`(2,10)`
joint controller is the identical padded state `(5)`, where the two faces
require opposite rows.  The corrected near miss is that this was special to
rows 2 and 10.  The least-used sidecar is the static face bit.

The bounded concept board was:

| object | operation | preserved | lost | decisive test |
|---|---|---|---|---|
| tagged state `(f,n)` | forget `f` | pole multiplicities | response availability | intersect the two availability sets |
| common covering pair `{i,j}` | restrict the 22-row atlas | a Q-monotone route exists on each face | other lawful rows | enumerate exclusive labels |
| conflict locus `X_{i,j}` | intersect over pairs | literal hostile states | pair orientation | test the singleton `(5)` |
| face bit | restore `f` | both facewise covers | within-face geometry | choose the first available row |

No tree DP was run: identical-feature conflicts decide the policy-independent
question before any choice of test language.

## 2. Exact conflict definition

For face `f` and a common pair `e={i,j}`, let

~~~text
A_f^e(n) = {r in {i,j}: row r has a strict Q_f-directed one-pole ascent at n}.
~~~

Both facewise covering statements say that this set is nonempty at every
nonreset state.  Pad the small count vector with `n_6=n_7=n_8=0`; all 239
small vectors occur on the large face.  The small reset imposes no transition
condition, leaving 238 common nonreset inputs.  Define

~~~text
X_e = {n: A_12^e(n) intersect A_13^e(n) is empty}.       (1)
~~~

Because both availability sets are nonempty subsets of a two-element pair,
membership in `X_e` is exactly an opposite-exclusive-label conflict.  Hence a
face-blind selector with values in `e` exists if and only if `X_e` is empty.

## 3. One singleton gives a bipartite obstruction

At the literal state `(5)`, with padded vector

~~~text
(0,0,0,0,1,0,0,0),
~~~

the exact available-row sets among all 22 charts are

~~~text
support (1,2): {2,8,9,11,12,14,16,18,22},
support (1,3): {3,4,5,6,7,10,12,13,17,19,20,21}.       (2)
~~~

Row 12 is the sole common row in `(2)`, but it occurs in none of the 24 common
covering pairs.  The vertices used by those pairs split as

~~~text
S_12={2,9,11,14,16,18,22},
S_13={3,7,10,13,17,19,21}.                             (3)
~~~

Every one of the 24 common pairs is an edge crossing the bipartition
`S_12 | S_13`.  Thus at `(5)` its `S_12` endpoint is exclusively lawful on
the small face and its `S_13` endpoint is exclusively lawful on the large
face.  This single calculation proves, simultaneously and without a policy
model, that no common pair admits a face-blind count-vector selector.

In the canonical increasing order of each pair, nine conflicts have
`small=left, full=right` orientation and fifteen have the reverse orientation.
The orientation is only an ordering convention; the bipartition `(3)` is the
intrinsic statement.

## 4. Full hostile atlas

The complete conflict counts, in THM-3249's common-pair order, are:

| pair | `|X_e|` | orientation |
|---|---:|---|
| `(2,10)` | 19 | small left / full right |
| `(2,13)` | 13 | small left / full right |
| `(2,19)` | 16 | small left / full right |
| `(3,9)` | 16 | small right / full left |
| `(3,11)` | 14 | small right / full left |
| `(3,16)` | 23 | small right / full left |
| `(3,22)` | 11 | small right / full left |
| `(7,18)` | 11 | small right / full left |
| `(7,22)` | 8 | small right / full left |
| `(10,11)` | 19 | small right / full left |
| `(10,16)` | 8 | small right / full left |
| `(10,22)` | 24 | small right / full left |
| `(11,13)` | 13 | small left / full right |
| `(11,17)` | 14 | small left / full right |
| `(13,14)` | 8 | small right / full left |
| `(13,18)` | 16 | small right / full left |
| `(13,22)` | 16 | small right / full left |
| `(14,19)` | 12 | small left / full right |
| `(16,17)` | 19 | small left / full right |
| `(16,21)` | 18 | small left / full right |
| `(17,22)` | 10 | small right / full left |
| `(18,19)` | 17 | small left / full right |
| `(19,22)` | 16 | small right / full left |
| `(21,22)` | 8 | small right / full left |

For every pair, `(5)` is the unique conflict of minimum pole-cardinality one.
The union of the 24 loci contains only 47 of the 238 common nonreset vectors.
Their intersection has six states:

~~~text
(5), (4,5), (1,3,5), (1,4,5), (3,4,5), (1,3,4,5).    (4)
~~~

Thus the obstruction is concentrated but robust: six inputs reverse every
common pair, while the individual loci range from 8 to 24 inputs.  The exact
conflict-pair multiplicity histogram is

~~~text
number of pairs : number of states
1:6, 2:5, 3:8, 4:3, 5:5, 6:3,
7:6, 8:1, 9:1, 11:1, 20:2, 24:6.                     (5)
~~~

## 5. Exact selector-bit theorem and its boundary

For each of the 24 pairs separately:

1. zero face bits are impossible, already at `(5)`;
2. one static face bit is sufficient: after reading the face, choose the
   first member of the pair available in that face's exact cover bank; and
3. therefore the minimum static face-origin sidecar is exactly one bit.

The construction was checked on all 238 small nonreset states and all 4,318
large nonreset states for every pair.  The word **static** matters.  This does
not show that a uniformly initialized or history-dependent previous-chart bit
can recover the face, nor does it optimize threshold trees, DAGs, formulas,
queries, switches, or routes.  Pair identity is supplied and not charged as
sidecar information.

The singleton `(5)` does not itself obstruct a selector allowed to use
arbitrary rows among all 22; the present result makes no claim about that
larger selector problem. Nor does it compare other support/bank faces or imply
Gaussian-moment positivity.

## 6. Reproduction and integrity

The companion pins the promoted THM-3249 script/output, the direct THM-3244
large-face script/output, and the repaired joint-sidecar script/output.  It
reconstructs all 239 small response vectors from exact product-Gamma formulas,
replays THM-3244's exact 4,319-state row-cover certificate, checks all 24
facewise covers, enumerates every conflict locus, and contains no `assert`,
float, randomness, or cache.  Normal and optimized Python runs byte-match.

~~~text
python3 04-computation/fc3_common_covering_pair_face_blind_conflict_atlas_20260803.py
python3 -O 04-computation/fc3_common_covering_pair_face_blind_conflict_atlas_20260803.py

script_sha256 = c03a837f1ed2fbbadc4c9aaef8609a79b1411a1898a56a57546e5460f1fdca56
output_sha256 = 40d074ea26d838bd43edea52b2cdaffeaf27a6a2cdb0dc9abedcd6156eed0e82
conflict_atlas_sha256 = 5faf9b8fc55cc49c83caa436d78ed98f4253e1c369e7a625ddce91cd69e61b1b
~~~

The next sharply typed question is whether the six-state universal locus in
`(4)` follows from a symbolic response-sign factorization, rather than from
the present exact census.  That question is open; the finite selector theorem
above does not depend on such a factorization.
