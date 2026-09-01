# LRC14 endpoint-592 response/exchange session — 2026-09-01

## Status

The anchor succeeded: THM-4313 closes the 35-row endpoint-592 layer by a
43-addition/43-deletion exchange while preserving the 3,925-mask carrier size
on all 467 rows. The typed residual now starts at endpoint 591 on 13 rows.
LRC(14) remains open.

## Inheritance pass

- Closest proved mechanism: THM-4311's one-response, one-deletion exchange and
  complete singleton quotient.
- Canonical hostile: the THM-4311 carrier fails 2,468 endpoint-592 bodies,
  dominated by 2,101 failures at q=256.
- Corrected near miss: THM-4282's endpoint-650 packing20 has zero literal
  overlap with the new q=256 failure fibre, and its index-256 replacement
  masks are all inactive at endpoint 592.
- Least-used sidecar: exact responder multiplicity and pairwise co-response,
  rather than total response weight alone.

## Anchor / Niche / Wildcard board

1. Anchor — close endpoint 592 at unchanged carrier size.
2. Niche — find a rigorous response-cover packing lower bound.
3. Wildcard — test whether a retained-pool LP dual globalizes to every
   responder.
4. Control — transport old endpoint-650 packing/surgery objects.
5. Next scale — expose the exact endpoint-591 frontier.

The anchor produced THM-4313. The niche produced a global pair-tagged packing8.
The wildcard failed constructively: exact pricing found 20,986 omitted-mask
dual violations. The control failed twice and explained why old surgery masks
live on the wrong activity side.

## Mechanism

The useful separation was

~~~text
complete frozen-failure response census
  -> explicit cover witness
  -> exact one-deletion quotient
  -> deterministic multi-deletion proposal
  -> independent simultaneous raw replay.
~~~

No implication was taken from singleton-safe to simultaneously safe. The
chosen 43 deletions were merely a proposal until the 6,681,439,050-case replay
returned zero failures.

The pair-tagged packing is also genuinely mixed-row: seven mutually
incompatible q=256 bodies extend by one q=105 body. The six bodies with zero
rank-eight response do not form a cover obstruction, because 151 rank-nine
masks respond to all six. Rank scarcity and cover scarcity are different
coordinates.

## Pool/global loss

The 37,497-mask retention rule preserves enough structure to find a 43-cover
and to certify a pool-relative lower bound 36. It destroys global dual
feasibility: 20,986 omitted responders exceed the dual cap, with maximum score
3.849498499. This is a concrete instance of a recurrent repository rule:

~~~text
good primal compression does not imply good dual compression.
~~~

The missing sidecar is exact pricing against omitted columns. Future
retained-pool dual claims should price the complete universe before promotion.

## Next generated tasks

1. Replay the new carrier on the 13 endpoint-591 rows and freeze its exact
   hostile failure universe.
2. Compute complete rank-eight/rank-nine responder counts and per-obligation
   multiplicities for that universe.
3. Search mixed-row response packings before optimizing a cover; q=256 showed
   that a fibre-only packing can often gain one from another row.
4. Reuse the full-467 singleton ledger incrementally: additions cannot create
   inherited singleton obligations.
5. Price every retained-pool dual over the complete response universe before
   treating it as a global lower bound.
6. Keep THM-4306's H_265 deck separate until an explicit transport map
   preserves activity, response, and deletion predicates.
