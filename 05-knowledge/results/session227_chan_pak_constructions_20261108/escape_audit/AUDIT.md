# Independent audit: exact graft boundary escape

**Verdict: PASS.** The all-order proof is sound, and the independently written,
portable verifier passes normally and with Python optimization enabled. No
mathematical correction was required. This audit does not establish scalar
completeness or a new unrestricted coverage range.

## Audited claims

1. **Inherited defect bound and sharp zero set.** With
   `f=U-D-d+E`, the difference counts exactly the two nonempty path bins on
   `A-a` whose joining endpoints are distinct out-neighbors of `a`. For each
   unordered endpoint pair, the supplied cut-bit inequality telescopes to two
   covers for each Hamilton path on the remaining vertices. Deleting those
   endpoints recovers the cut and underlying word, so the counting injection
   has no collisions. Empty remainders are included. Consequently
   `f >= 2 sum H(A-{a,u,v}) >= outdeg(a)(outdeg(a)-1)`, and `f=0` exactly when
   the root has outdegree at most one.

2. **Attribution checked.** The lower bound is already THM-4875(A1), which
   explicitly credits THM-4873(ii)(D). These are the long-slug theorem files
   cited by the producer. MISTAKE-716 specifically warns against treating
   inherited endpoint bounds as new. The draft accurately marks the bound
   inherited and gives a self-contained proof. Session198's
   `session198/bilinear-descent/PROOF.md`, Sections 1–3, supplies the earlier
   start-one grammar and terminal faces. The present unrestricted last-letter
   observation is immediate from the inherited O/Z/J matrices.

3. **All predecessors, not only a seed orbit.** Every nonempty O/Z/J word
   finishes on at least one of `E=0`, `E=d`, or `f=0`, since its final letter
   forces that identity. This actually holds algebraically for arbitrary
   vectors. Thus the proof genuinely excludes every nonempty word from every
   physical predecessor, even input-dependent words of unbounded length. It
   does not exclude the empty word, and does not classify all points on or
   outside the three faces.

4. **Exact graft criterion.** For every nonempty old rooted tournament and
   every rooted gadget of order at least two:
   - `d_G-E_G >= d_A H(B-r)>0`. The separately counted deleted paths ending
     at the old root cannot be extended to the new root.
   - `E_G>0` iff the old root is a nonsink. The inherited endpoint formula,
     positivity of the marked-bin cover count, and `U-d >= D-E>0` prove the
     forward case; both terms vanish for a sink.
   - `f_G>0` iff the gadget root is a nonsink, because its output outdegree
     is exactly one plus its old gadget outdegree.
   Hence escape is equivalent to both input roots being nonsinks.

5. **Degenerate boundaries.** The rooted singleton old tournament is a sink
   and gives `E_G=0`. Empty induced deletion graphs contribute one Hamilton
   path, as stipulated. A singleton gadget is J and has `f_G=0`; it is
   deliberately outside the two-or-more-vertex theorem. A sink-rooted
   nonsingleton gadget also has `f_G=0`. None of these cases is silently
   dropped by a strict-cone assumption.

6. **Minimal edge and scalar hostile.** The source-rooted directed edge has
   the stated matrix L and sends `(1,0,2,1)` to `(5,1,9,3)`, whose three
   strict gaps are `E=1`, `d-E=2`, and `f=2`. A singleton gadget cannot escape,
   proving minimal gadget order two for this particular state property.
   The other rooted edge is ZJ. L itself is two Z births followed by marking
   the first newcomer. Order four is also minimal among all physical escaping
   states: at orders at most three a root has outdegree at most one, giving
   f=0, or is a global source, giving E=0. The witness is rank three, not full rank.
   Its scalar H=5 already occurs at `t.OJ=(5,1,5,1)` and also at
   `t.JJ=(5,2,6,3)`. Thus state separation cannot be promoted to scalar
   separation. The explicit five-vertex closing-J tournament has state
   `(13,4,14,5)` and is strong; H=13 is not claimed as a new attainable value.

7. **Additional sidecars.** The root-start multiplicative law follows from
   the unique old-root entry seam. It uses `alpha_B S_A` for old order at
   least two and correctly uses `s_B` instead for an old singleton.
   O/Z/J preserve root starts only under the stated nonsingleton hypothesis.
   The L image plane, kernel, nonzero minor, induced determinant -1,
   two-sided integer inverse, cubic matrix identity, and displayed families
   all check. The proof does not infer physical inverse closure from a
   unimodular integer inverse. Strongness follows from the old root reaching
   all old vertices and the cross-edge directions; the edge gadget itself
   need not be strong.

## Independent finite-exact verification

`independent_check.py` imports no producer or repository code and requires
only the Python standard library. It counts Hamilton paths by a memoized
first-vertex deletion recurrence, with a separately written last-vertex
recurrence for endpoints. Literal permutation/split enumeration checks all
small profiles. Explicit exceptions remain active under `python3 -O`.

The deterministic receipts record:

- All 1,099 labelled tournaments of orders 1–5 with **every root**:
  5,405 rooted instances. All profiles and root-start counts are compared
  with literal permutations; exact unmatched-cover identities, degree iff,
  and 7,878 individual out-neighbor-pair lower bounds are checked.
- All 21,845 binary cut-bit sequences of lengths 0–7 verify the elementary
  telescoping identity and the lower bound, including empty remainders.
- All 75 labelled tournaments of orders 1–4 as old inputs against all 75
  gadgets, both with root fixed at label zero: 5,625 pairs. Of these, 75
  have singleton gadgets and 5,550 satisfy the main gadget-order hypothesis.
  Exactly 3,969 escape; all sink combinations and singleton exceptions are
  retained. The root-start product or singleton exception is checked on
  every pair.
- 225 literal O/Z/J controls and 1,875 controls on arbitrary integer vectors
  test the final-letter identities. These are finite controls, while the
  identities themselves prove the all-word statement.
- The source-edge matrix on all 75 old inputs, the sink-rooted edge's ZJ
  state action, the exact H=5/H=13 witnesses, 125 signed integer controls on
  the restricted image plane, the inverse, kernel and cubic recurrence.
- Literal L iterates through `m=3` (old orders 2,4,6,8) and their closing-J
  graphs verify the claimed states, root starts and strongness. Later
  displayed recurrence entries are algebraic, not literal enumerations.

These checks support the reviewed proof; they are not an extrapolation from
finite data. The producer has its own independent append-endpoint DP,
permutation controls, all 5,625 graft pairs, and literal iterates through
`m=4`. Its verifier was also rerun in this audit in both modes.

## Reproduction

From this audit directory:

    python3 independent_check.py > independent.replay.json
    cmp independent.replay.json INDEPENDENT.json
    python3 -O independent_check.py > independent.replay.optimized.json
    cmp independent.replay.optimized.json INDEPENDENT.optimized.json
    cmp INDEPENDENT.json INDEPENDENT.optimized.json

To replay the separately supplied producer checker, substitute its directory
for PRODUCER:

    python3 PRODUCER/verify.py > producer.replay.json
    cmp producer.replay.json producer.normal.json
    python3 -O PRODUCER/verify.py > producer.replay.optimized.json
    cmp producer.replay.optimized.json producer.optimized.json

`AUDITED_INPUTS.json` identifies the exact producer input bytes. The independent
checker has no path dependency on those files. `SHA256SUMS` covers this audit's
publication artifacts, excluding the checksum file itself and startup notes.
