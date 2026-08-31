# LRC(14): mixed-rank deletion depth and recursive incoming signals

**Session:** opus-lrc14-incoming-20260831
**Canonical outcome:** [THM-4296](../01-canon/theorems/THM-4296-lrc14-mixed-rank-deletion-depth-and-recursive-signature-closure.md)
**Status:** finite exact progress; LRC(14) remains open.

## Inheritance pass

- **Closest proved mechanism:** THM-4254's repair hypergraph. An active
  deletion mask disjoint from a nine-body certifies that body by safe-set
  inclusion.
- **Canonical hostile:** THM-4287's repaired 9,006-mask carrier has exactly
  101 failures at `(100,636)` and `(256,636)` even though its older prefix
  already closes endpoint 637. Endpoint monotonicity is false.
- **Corrected near miss:** comparing cocycle margin numerators from different
  ordered pairs as if they shared a unit. MISTAKE-532 records why signs
  survive but cross-row weakest-margin rankings do not.
- **Least-used relevant sidecar:** the *complete ideal* below an inactive
  signature, not merely the equal-signature fibre at one target. The ideal is
  the natural recursive object because deleting indices exposes bodies
  monotonically while common activity must be intersected across every row.

## Live concept board

| concept | representation | preserved predicate | chief information loss |
|---|---|---|---|
| deletion depth | upset in the Boolean lattice of pool deletions | existence of a lawful disjoint active repair | one scalar depth forgets which mask and which bodies it covers |
| response quotient | bit pattern on the exact failure ledger | labelled body coverage at one pair | forgets activity at other pairs unless carried explicitly |
| inactive exchange | delete strictly inactive carrier masks, append responders | carrier size and all audited row covers | resembles basis exchange but need not satisfy matroid exchange |
| signature ideal | rows with `I_E(p)` contained in a deleted-index set | common-deck surgery after full activity intersection | signature alone forgets appended-mask activity |
| typed proof graph | row sets labelled by certificate consumer | lawful union of proved consequences | mask-level union would falsely merge incompatible decks |

Every useful incoming result changed at least two entries on this board. The
endpoint cover became an exchange only after the inactive-mask sidecar was
added. The rank-eight hostile changed the response quotient into a depth
filtration. The singleton census became a theorem only after its 110 decks
were kept separate in the proof graph.

## Signal 1: the endpoint cover was really an exchange problem

The 101 endpoint-636 failures do not admit a common-active fourteen-mask
solution; that restricted minimum is fifteen. They do admit a fourteen-mask
solution in the union of masks active at the two rows. This distinction is
exactly what a carrier permits: each mask need only be active where it is
used.

The incoming inactive-mask scan then found fourteen old carrier masks that
are strictly negative on every endpoint-637/636 row. Removing them loses no
certificate on this boundary. The response witness can replace them, keeping
the carrier at 9,006. The lesson is reusable:

```text
response cover + globally inert old coordinates = size-preserving exchange.
```

This is only a basis-exchange analogy. The active-mask families vary with the
pair, and there is no proved matroid exchange axiom. A cheap future hostile is
to select two minimal carriers at one row and test symmetric exchange; one
failure would prevent importing matroid language.

## Signal 2: rank eight was a coordinate choice, not a theorem hypothesis

The exact body `1d106401` at `(256,632)` has no active disjoint rank-eight
mask. Optimizing harder inside rank eight cannot help. The useful response to
that negative result was to reread the consequence direction rather than the
search code.

For any deletion set `R`, disjointness gives

```text
B subset P\R
  => G_((P\R) union {q,r}) subset G_(B union {q,r}).
```

The rank never enters. Activity is upward closed under adding deleted labels.
This produces the intrinsic integer

```text
delta_(q,r)(B)=minimum active deletion rank disjoint from B.
```

For the hostile body, the exhaustive rank-eight complement is empty and all
thirteen one-bit extensions of its best rank-eight mask are active, so its
depth is exactly nine. This is more than a successful workaround: it explains
why the previous lane stopped.

Deletion depth behaves like a stopping time or a persistence birth index.
The sublevel sets `{B:delta(B)<=k}` are nested, and exact-rank scans compute
their cumulative filtration. That analogy is lawful because the map and
monotonicity are explicit. It does not imply any topological stability or
probabilistic law.

## Signal 3: dual weights predicted which local repairs would scale

At `(100,629)`, the mixed response quotient has 419 classes. An incoming
independent audit found the exact quarter-unit dual

```text
weights 0:1, 1:1, 3:1, 4:1, 5:2, 7:1, 8:1, 9:2,
        23:2, 24:1, 26:1;
total 14/4=7/2,
every response-class load <=1.
```

Thus three masks are impossible, while four explicit masks cover all 28
obligations. The dual is not only a lower-bound certificate. Its 42 tight
classes identify the exact face on which alternate four-mask witnesses live.
For the 995-body endpoint-626 problem, the first useful computation should be
the maximal response classes and their rational dual before any large branch
search. A dual concentrated on a small body subset may expose a structural
packet; a diffuse dual warns that local mask addition will scale poorly.

## Signal 4: recursive ideals outperformed target-local fibres

There are 3,738 current rows with singleton inactive signatures. Grouping by
the unique inactive index gives 192 disjoint complete ideals. For 110 indices,
one replacement mask is common-active on the entire ideal and covers every
private body of the deleted joint mask. Those 110 decks close 1,219 rows.

The failed index 19 then supplied the next recursive task. It has no single
common responder, but two response classes cover its eight obligations. This
gives a separate 422-mask deck on 36 rows. The successful recursion is

```text
singleton one-replacement census
  -> inspect a failed singleton ideal
  -> exact two-replacement response quotient
  -> promote only after a full common-activity/body replay.
```

The natural next objects are the remaining 81 failed singleton groups after
index 19, followed by small two-index down-ideals. The cheap order is not by
row count alone. Rank each group by `(obligation count, maximal-response
defect, dual lower bound, common-active density)` while keeping all four
coordinates; do not collapse them into a cosmetic score in the proof.

## Typed union rather than false unification

The incoming proof-graph audit was decisive because the three mechanisms have
different semantics:

```text
110 separate 421-mask decks: 1,219 rows
one 422-mask index-19 deck:      36 rows
one 9,019-mask mixed carrier:    70 rows.
```

There are four pairwise overlap rows and no triple overlap, so the typed union
has 1,321 rows. An independent C++ consumer derives the index-19 ideal,
parses the raw 72-row carrier scan, checks that it is the complete inherited
`r>=626` prefix, and obtains the same 21,326-row complement. This is the
correct place for union: consequences, not decks.

A late incoming theorem then supplied a distinct append-only carrier and three
recursive ideals. Treating it as signal required recomputation, not addition
of headline counts. Its ten-row carrier set is contained in the 70-row mixed
carrier set, and its 36-row `H19` equals the index-19 node above. The `H294`
and `H372` decks remain genuinely independent, but 115 of the incoming
theorem's 118 consequence rows were already present. Exact overlap accounting
leaves only

```text
(147,294), (147,590), (372,619).
```

The aggregate typed union is therefore 1,324, not `1,321+118`, and its
complement has 21,323 rows. Independent Python and C++ consumers reproduce
every intersection and the three output ledgers byte-for-byte. This is a
useful concurrent-research rule: preserve distinct Pareto constructions and
certificate identities, but always recompute their consequence join in the
common inherited universe.

Tournament language is not helpful here unless a genuine pairwise observable
is chosen. The intrinsic relation among deletion masks is inclusion, a
partial order with many ties, not a tournament. One could orient two response
classes by strict dominance only after retaining incomparability as a
sidecar; forcing all incomparable pairs to point somewhere would erase the
exact-cover geometry.

## Anchor / Niche / Wildcard next portfolio

### Anchor: the 1,005 endpoint-626 failures

1. Split the 995 and 10 bodies and compute exact depth histograms beginning at
   ranks eight and nine.
2. Freeze maximal mixed response classes before exact cover search.
3. Compute a rational dual and its tight-class face.
4. Test whether a size-preserving exchange exists by scanning carrier masks
   strictly inactive on the complete already-closed 70-row node plus the two
   endpoint-626 rows.
5. Only then append witnesses and rescan the next complete endpoint layer.

The hostile control is coverage collapse from `(6)`: rank ten is easier to
activate but each mask is disjoint from only `C(20,9)` bodies. A small depth
does not guarantee a small carrier.

### Niche: failed signature ideals

For each of the remaining singleton groups, compute the exact common-active
response quotient and minimum replacement number. Then test two-index ideals
selected by shared private-body structure. A positive node must include its
full ideal, common activity, rebuilt deck, strict margins, and proof-graph
overlaps. A local target witness is not enough.

### Wildcard: the deletion-depth profile as a finite sequence

Order bodies by colex only as an address, then study the sequence
`delta_(q,r)(B)` together with response-class and endpoint sidecars. Useful
questions are finite and decisive:

- Does depth depend only on the rank-nine failure atom? Test duplicate atoms
  with different depths.
- Is the depth-`<=k` family compressed or shifted under colex? Test closure
  directly; one witness refutes it.
- Do adjacent endpoint rows have a bounded edit distance between depth
  profiles? Compare exact labelled bodies, not histograms.
- Does the dual-tight body set recur under a signature ideal? Intersect exact
  ledgers before proposing a bridge.

If any answer is positive, it yields a new representation or theorem target.
If negative, the minimal witness identifies the missing coordinate. Either is
more valuable than another untyped carrier search.

## Session boundary

The exact frontier is now

```text
21,323 residual rows,
maximum endpoint 626,
top rows (100,626) and (256,626),
1,005 failures for the current mixed carrier.
```

No physical entry, exclusive owner, semantic arrival, or LRC(14) proof was
obtained. The durable breakthrough is structural: deletion depth is a lawful
rank-free invariant, and recursive signature ideals plus typed proof nodes
give a disciplined procedure for turning incoming computations into the next
finite tasks.
