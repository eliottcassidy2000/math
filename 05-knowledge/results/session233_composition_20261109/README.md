# Session233: nonunit loops compose into consecutive odd windows

**Session233 complete: three independently audited results. Chan–Pak completeness remains OPEN.**
This session develops positive constructions from the earlier symmetry,
endpoint, graft and moving-core ideas. It does not prove that every admissible
odd integer is a complete Hamilton-path count.

## Main advance: restore unit spacing by composing nonunit loops

[THM-5000](../../../01-canon/theorems/THM-5000-coprime-interface-composition-consecutive-odd-windows.md)
proves that, for every L>=337, one full physical interface supports

    377L−126687 consecutive odd Hamilton-path counts
    on 12L+962 strong-tournament vertices.

For L=1000 this gives **250,313 consecutive odd counts on 12,962 vertices**.
The asymptotic interval-width/vertex ratio improves 377-fold relative to the
inherited single twelve-bit moving-core construction. This compares these
explicit constructions, not all possible constructions.

The key is to retain two short loops previously unhelpful for unit steps:
an h=13 core changes the cut coordinate by10; an h=29 core changes it by6.
Together with the inherited h=45 unit core, their actual composed offsets are
2(6525i+1755j+377t). The 377 fixed choices i=0..12, j=0..28 cover every residue
modulo 377. A bounded CRT selector then fills an entire interval, with a real
strong tournament for every selected target in that interval.

[Weighted proof and literal certificates](graft/WEIGHTED_PROOF.md) ·
[Independent supplement audit](graft_audit/SUPPLEMENT_AUDIT.md).
The earlier [two-unit-core theorem](graft/PROOF.md) is retained as a simpler
composition lemma and valid checkpoint, though its fixed overhead is larger.

**Remaining gap:** the interval locations grow exponentially while the widths
grow linearly. The example's lower endpoint has 12,966 bits. There is still no
selector placing one of these windows over an arbitrary admissible target.

## Complete endpoint spectrum

[THM-5001](../../../01-canon/theorems/THM-5001-every-positive-endpoint-quotient-has-a-small-strong-realizer.md):
every positive integer N is an endpoint count and endpoint-stabilizer quotient
of a strong tournament on at most max(5,ceil(log2(N+1))+3) vertices.

A rigid near-transitive base has terminal vector (1,1,2,4,...,1). Binary digits
select the incoming neighbors of one new endpoint. Thus the former isolated
local 7/21 examples become a complete constructive boundary theorem. The
rooted path-cover response obeys a different exact linear law, so selecting
the endpoint count does not independently select the full count or graft
response.

[Proof and examples](endpoint/PROOF.md) · [Independent audit](endpoint_audit/REPORT.md).

## Exact complementary-block transfer and its boundary

[THM-5002](../../../01-canon/theorems/THM-5002-three-state-transfer-for-bridge-linked-transitive-blocks.md)
gives a three-state transfer for the entire odd-cycle packing polynomial of
transitive blocks joined by single forward bridges. Both canonical run-window
endpoints, even closed under ordinal sums, miss 11. Their tight packing-number
coverage therefore cannot by itself prove spectrum completeness.

The enlarged family of arbitrary block compositions attains exactly 1,374
values at most10,000, proved by an exhaustive bounded-value census and the
inherited strong-order floor. This is not a census of all tournaments.

[Proof](transducer/RESULT.md) · [Independent audit](transducer_audit/AUDIT.md).

## Verification and inherited corrections

Every substantive theorem has an independent proof/computation audit. The
final compact replay checks all producer and audit programs normally and
under optimized Python, with deliberate failure controls. Large constructed
graphs are counted by proved transfer identities, not literal enumeration.
See [reproduction](REPRODUCE.md), REPLAY.json and SHA256SUMS.

MISTAKE-744 records the corrected preliminary head-extension transfer, live
optimized checks, and two repaired inherited THM-250 prose errors. None
changes the prior near-transitive total-count formula.

[What remains to prove](SYNTHESIS.md). The first two verified nonforce
checkpoints were 2ecb04e29ec7b0c322a247971a6acdca25febb36 and
64e98b1f6c73599d6174b89453d7bc494aa929e2. Final publication preserves the
concurrent creator, ternary and Paley reservations already on main.

The final integration also preserves the concurrent pre-audit checkpoint67fa9a1f9 on Paley doubles and even-graph domination. Those separate carrier computations are not dependencies of the proofs here.
