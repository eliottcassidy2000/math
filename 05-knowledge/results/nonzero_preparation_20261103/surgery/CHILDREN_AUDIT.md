# Independent audit of the actual second-round fork

**PASS — full analytic/source audit, independent normal/optimized replay,
and independent reconstruction of every native intermediate.** This audits
[CHILDREN.md](CHILDREN.md) and its exact finite universe of four continuations.
The first-round dependency has separately passed [AUDIT.md](AUDIT.md).
Promotion and Git remain with the parent.

## 1. Four actual marked continuations

The inputs are the actual first children at root8: order13, H4291, F4479;
and order14, H4289, F5769. The source reconstructs them from the pinned
first-round programs rather than substituting equal-count representatives.

I independently constructed all four funnel/converse targets from their
literal cores, entry sets and tail words, then applied the recorded physical
permutations. Every reconstructed adjacency equals the displayed target.
Each map sends the canonical root to old physical8 and its auxiliary to the
stated physical name. All four target orders equal their respective input
orders, and no old name is deleted. The permutations identify target graphs
on the existing name sets; the programs implement every actual changed arc.
They do not relabel the old input for free.

In each target the auxiliary's outgoing row is exactly `1<<8`. The accepted
sole-exit identity therefore supplies `L=F_8(T-a)`,
`K=P_(T-{a,8})(2)` and `G=L+K=F_8(T)`, without requiring strongness.
The source independently counts both deletion graphs. Its target values are

|Input H|Gap|F=G|S|L|K|Target defect|
|---:|---:|---:|---:|---:|---:|---:|
|4291|1|4563|4292|1491|3072|272|
|4291|2|8371|4293|1715|6656|4080|
|4289|1|6983|4290|2887|4096|2694|
|4289|2|12094|4291|3902|8192|7805|

Thus the four C2 children have `(H,F_8)` equal respectively to
`(12875,13689)`, `(12877,25113)`, `(12869,20949)`, and `(12871,36282)`.
The exact scalar and root-response laws are both checked. No inference from
H-neutrality alone is needed or made.

## 2. Every native operation and the weak boundary

I read the complete new source and its pinned dependencies. The four old-pair
schedules have5,8,4,6 flips. Every one uses the inherited actual helper birth,
two guarded Hamilton complements and helper-only deletion. The source checks
both endpoints' nontrivial sparse helper cuts, exact guard names/directions,
named-set replay, retention of all old names, final target equality, and
reversed/repeated-name hostile guards. The imported compiler is pinned by its
full source hash; no unfrozen search or response receipt enters the program.

In addition to both producer replays, a separate verifier with no producer
imports reconstructed all46 Hamilton complements. It used the alternative
implementation “complement all arcs, then restore consecutive guard arcs.”
It checked every helper cut, guard, intermediate adjacency, and helper-only
deletion against the receipt, including that each completed gadget changes
exactly one declared old pair. Independent reachability checks and forward
endpoint dynamic programming verified all69 helper/native states and their
actual root8 H/F/defect observers. The same recurrence checked every completed
stage, each target and the final children. The supplemental verifier passed
401 always-active checks, including all four marked canonical maps.

All69 helper/native states are strong. The only weak completed stage is the
sixth and final flip of the fourth branch: the deletion after `(0,3)`. At
that point physical vertex3 dominates every other retained vertex. Removing
it leaves a strong thirteen-vertex ready core. I independently checked this
exact source-plus-strong-core decomposition at the target and both permanent
C2 stages. The C2 births keep vertex3 as global source and thus remain weak.

The other three branches have every completed adjacency stage, both C2
stages and the final child strong. This verifies the proof's precise scope:
three strong leaves and one explicitly weak leaf. It would be false to
advertise this receipt as an all-strong four-leaf tree, and the primary does
not do so. Root8 is retained at every stage, including helper stages where
an appended array entry is not the root.

## 3. Padding and the retained-order obstruction

The padding argument is exact for any transitive prefix of length l. In a
Hamilton path the prefix must appear first in its unique order, so H is
unchanged. In a two-slot cover counted by S, the slot starting at a cannot
contain prefix vertices, so all of them are forced into the other slot and
S is unchanged. For F, each prefix vertex can independently be assigned to
either slot, in increasing prefix order before that slot's old vertices.
This is a reversible `2^l`-to-one construction, including an otherwise empty
ordinary slot. The same reasoning applies to L and K: their distinguished
A_a or Q_ar slot cannot contain a prefix vertex, leaving two free slots.
Thus `(H,S,F,G)` transforms to `(H,S,2^l F,2^l G)` and zero joint deficit
persists.

This proves readiness of the padded target but not strongness; a nonempty
prefix is a genuine source component. In the actual branch it represents an
existing old name, not an unrecorded birth or a deleted old name. The source
constructs precisely this full fourteen-name graph. The retained-order
obstruction is therefore handled honestly: the thirteen-name scalar witness
can be padded to a weak actual target, but has not been replaced by a strong
fourteen-name ready witness.

The four continuations form a finite depth-two construction only. The note
correctly leaves a strong continuation of the H4289 +4 branch, and an
unbounded returning two-branch class, open. No absence conclusion is drawn
from the bounded enlarged-core scout.

## 4. Independent replay and pins

From the session worktree:

```
python3 -B 04-computation/nonzero_preparation_20261103/surgery_children.py --output /tmp/nonzero_surgery_children_audit.json
python3 -B -O 04-computation/nonzero_preparation_20261103/surgery_children.py --output /tmp/nonzero_surgery_children_audit_optimized.json
```

Both executions reproduce the frozen84,963-byte receipt and191-byte stdout
exactly:634 gates,4 cases,46 complements,69 strong helper/native states.
I rechecked every `CHILDREN_MANIFEST.json` file length and digest. The supplemental401-gate independent verifier is preserved as well:

```
python3 -B 04-computation/nonzero_preparation_20261103/surgery_children_audit.py
```

Its stdout must match `CHILDREN_AUDIT.out`. It reads the frozen child receipt
and imports no producer implementation. The source is5090 bytes, SHA256
`47d7b0dd685a8e5cf076c592f6b281bff4e417ec431e8863937dda7ff6a3f5b3`; stdout is131 bytes, SHA256
`1ccd07e63f05406b33ca6d545ab6c0246c5281c0b8c3f4ea46a488b72fd57c46`.

|Artifact|Bytes|SHA256|
|---|---:|---|
|surgery_children.py|10698|4142db9f3ee60b1b9a6bba2c017d0a38d4994895472f9f294122979c0d12d9a8|
|children.json|84963|001a49cd1517da14752c3ec96559b688f0802393bf7e2fd28667017367456aac|
|children.out|191|bca81dfdbc7dcffd529e323f5448492a0409b1797946eabceafb229f92e814a0|
|CHILDREN.md, audited prepromotion version|7571|b6a4bf810a5e329c12183f559b10eb85ee9a3f132a95fd93116220f66c68e52b|

No correction is required. The accepted result is this complete actual
finite fork with retained old names/root, full changing joint response, and
its explicitly paid strong/weak boundary.
