# Independent audit of the four-marker sandwich

**PASS after one corrected degree typo.** This audit reviews the all-order
proof S1–S4 and its stated boundaries, not merely the producer's finite tests.
It supports promotion of these results to PROVED and closure of the specific
`A2*h*A2` all-core question in HYP-9416. Other response-classification and
dynamic-compression questions in that hypothesis remain outside this audit.

Reviewed [sandwich/README.md](../sandwich/README.md), promoted SHA-256
`227b96f05433d43f63f2ade93dd5de6ba9b78712f03b737bb1f15e0e349ea785`.
The inverse status/link-only edit reproduces the fully reviewed corrected
candidate `31c47969a899f887df248d79e359cfe7609a93b7b3d4b66f195e443bff887325`
exactly; the mathematical body and checker are unchanged.
The producer checker SHA-256 is
`9478dce8edda08d1c8acc03f43bcfc108bd13f7e3786a5846ebef11e5cae5dd0`.

## Exact correction and proof review

The initial candidate, hash
`37243dc9580f70097a528786e6b689692328d943d2e40fb42e4b18163a3165e1`,
said the right marker outdegrees were `n+1,n+2`. Direct construction shows
both are `n+2=N−2`, since both have indegree one. The producer replaced
that sentence with `Both right markers have outdegree n+2.` No conclusion
depended on the incorrect smaller value; both values exclude those vertices
from the degree-one set. The corrected mathematical statement was rechecked.

* **S1, all-order reconstruction: PASS.** Formula (1) gives exactly two
  outdegree-one markers, plus the entrance iff it is a sink. The extra case
  is a homogeneous directed triangle dominated by its entire exterior.
  The indegree-one argument is the dual statement. Since `N>=6`, these
  sets cannot meet. All three choices in either exceptional triangle differ
  by an actual exterior-fixing cyclic automorphism. Independent rotations
  give one joint rooted type even when both exceptions occur. The `n=2`
  reverse arc is explicitly covered; strongness is never used. The decoded
  output is a representative up to rooted isomorphism, so the quadratic
  adjacency-operation bound does not claim a general canonical-labelling
  algorithm.
* **S2, automorphisms and composition: PASS.** An unrooted automorphism must
  preserve both intrinsic extremal-degree sets. Its induced action on the
  quotient core fixes the ordered ports. The embedded port-fixing core
  automorphisms fix all exceptional triangle vertices pointwise, commute
  with their rotations, and intersect that rotation subgroup trivially.
  The kernel is exactly `C3^epsilon`, proving the direct product, not merely
  its cardinality. At the reverse two-vertex core the factor is `C3×C3`
  and the core group is trivial. Consequently the encoding is injective
  on isomorphism classes and faithful on maps; fullness requires the stated
  exclusion of both exceptions. Decode–hinge–encode is a conjugated
  operation on those classes. It is not the native unmarked path operation
  on the encoded tournament, and the note preserves this distinction.
* **S3, exact image guard: PASS.** The necessary cardinality tests alone
  are insufficient. A candidate must supply disjoint marker chains,
  distinct ports outside the markers, and the prescribed incidences.
  Exact re-encoding is an explicit sufficient witness. There are at most
  nine candidate pairs, and S1 proves choice independence after acceptance.
  The displayed six-vertex hostile has the claimed degree sets and fails
  because proposed ports are markers. No decoder promise is being silently
  applied to arbitrary tournaments.
* **S4, fixed-context optimality: PASS.** With nonempty contexts and distinct
  ports except in a singleton, the five listed context pairs exhaust total
  size below four. Literal isomorphisms verify every listed three-vertex
  collision, and the two admissible order-two core types remain distinct
  in all five contexts. S1 supplies a working context pair of size four.
  Thus the lower bound is exact in this grammar. It does not extend to
  arbitrary encodings, empty contexts, coinciding core ports, or a scalar
  observer of the encoded object.

The inherited [HYP-9416 Open clause](../../../hypotheses/HYP-9416-fresh-re-marking-two-port-quotient-joint-classes-are-minimal-and-line-responses-classify-strong-atoms.md#open)
was read at base `4950589d056a64f8807774cf7da893dcafcd63a0`, where this
all-core `A2*h*A2` question is explicitly open. The strong-core two-marker
theorem is inherited, not newly proved by the four-marker claim. Source
versions and the status boundary are preserved in [MANIFEST.json](MANIFEST.json).

## Independent finite route

[check.py](check.py) imports no producer or repository code. It constructs
the sandwich by the literal five-factor binary hinge operation with explicit
ports, independently of the proposed degree formulas. All labelled core
orientations at orders two through five are quotiented by internal-vertex
permutations, yielding exactly `2,8,32,176` rooted types.

For each of these 218 types it checks all marker-choice decodings and also a
separate intrinsic quotient: delete the unique two-vertex chains, or contract
the homogeneous extremal triangles. It enumerates every degree-compatible
vertex permutation without using the decoder and compares the **entire
automorphism set** with the proposed direct-product set, including explicit
commutation. The four boundary combinations occur `185,14,14,5` times.

All `32,768` labelled six-vertex tournaments are tested independently against
the literal relabelling orbits of the two encoded order-two cores. Exactly
800 lie in the image; 2,960 pass the degree-cardinality test. The receipt also
stores concrete isomorphism permutations witnessing each smaller-context
collision and checks 100 small decode–operate–encode compositions.

Reproduce from this directory:

```
python3 -B check.py > normal.json
python3 -B -O check.py > optimized.json
cmp normal.json optimized.json
```

Both modes pass **36,688 explicit gates** and produce the same bytes,
SHA-256 `b209187a6b3deee34c68385bcff31c44ee5bdf5a64a51527903b34c58ebcc2d7`.
The producer's 3,878-gate program was also rerun in both modes; both results
match its frozen receipt byte for byte. Those replays are saved separately
as `producer_normal_replay.json` and `producer_optimized_replay.json`.

No uncorrected mathematical or implementation blocker remains in the
reviewed scope. The general theorems rest on the proof above and the producer
note; the finite universe is a hostile control, not their all-order proof.
