# THM-4286 signature-response and two-deck-surgery packet

Status: **PROVED RELATIVE TO THM-4282/THM-4283 + FINITE-EXACT + DETACHED
LITERAL-WALL AUDITS PASS.**  This packet supplies two alternate exact
common-deck proof nodes on 72 disjoint post-THM4282 rows.  Current THM-4283
already contains all 72, so this packet contributes zero current ledger
rows.  It proves neither physical entry nor LRC(14).

## Consequence

The old 421-bit inactive signature does not determine any of the four audited
THM-4282 certificate-response predicates: each has mixed signature fibres.
Two refined surgeries nevertheless give explicit body-covering 421-mask
decks.

```text
index-396 deck:
  (E minus E[396]) followed by 042022c9
  common post-THM4282 rows                         36
  overlap with current THM-4283 proof union        36
  current incremental rows                          0

two-mask/shrunk deck:
  (E minus E[107],E[318],E[374]) followed by
  32043014,20807016,128c8012
  common post-THM4282 rows                         36
  overlap with current THM-4283 proof union        36
  current incremental rows                          0
```

The two 36-row branches are disjoint.  Their 72-row union has FNV
`4f8f4c79707540a6` and SHA-256
`253071871ede4041c658ac7705de5283794f1baa230c9df1e22d16e22ac830b3`,
and is a subset of current THM-4283's 691-row proof union.  The current exact
complement remains `22,682` rows, FNV `f7563445f15efebf`, SHA-256
`7d0044bc4c32f08b9d09dca420171a05666d26e03f38fbc48a9baa1fcb027102`,
with top `(100,637),(294,637),(520,637)`.

These are two separate common-deck certificate nodes.  Their row-set union is
not one common deck and is not interchangeable with THM-4283's carrier.

## Exact components

- `signature_response_primary.out` and
  `signature_response_independent.out` are Python and independent Perl
  replays of signature-fibre purity.  There are `17,604` signature fibres;
  the four predicates have respectively `6,15,34,48` mixed fibres.
- `singleton_fibre_cegar.out` exhausts the index-396 private obligations and
  response/activity quotient.  Exactly two of 495 one-mask body responders
  are active on all 36 fibre rows.
- `singleton_fibre_literal.out` independently constructs literal joint walls,
  proves strict activity on all 36 rows, and scans all `14,307,150` bodies for
  the rebuilt 421-mask deck.
- `singleton_fibre_literal_exhaustive.out` independently enumerates all
  `5,852,925` rank-eight masks with literal walls, recovers the `495`
  responders, `83` activity classes and two fibre-wide masks, and finds 35
  row-response profiles with sole collision `(301,366)/(366,547)`.
- `two_mask_fibre_cegar.out` computes the complete 20-obligation response
  quotient and proves exact fibre-wide replacement minimum three.
- `two_mask_fibre_literal.out` independently checks the five-row fibre, the
  422-mask body cover, and the unique extra removable retained index `318`;
  the final 421-mask deck passes the full body scan.
- `full_gain_two_mask.out` and `full_gain_index396.out` audit every old
  signature contained in the corresponding deletion set.  The minimum-ratio
  comparator uses an overflow-free continued-fraction comparison.
- `two_mask_net36_literal.out` independently reconstructs all 421 old margins
  on all 112 candidate rows, checks the frozen signature exactly, proves all
  three appended masks strictly active on 106 rows, intersects the inherited
  residual in 36 rows, and repeats the global body scan.
- `proof_graph_consequence.out` pins all inherited byte identities and proves
  that both 36-row sets are wholly subsumed by current THM-4283.  It records
  the 72-row alternate-node union and the unchanged current residual/top.
- `proof_graph_subsumption_independent.out` is a separately written Perl
  exact-set replay of the current THM-4283 complement, both full
  subsumptions, branch disjointness, the 72-row identity, and the top triple.

## Authoritative ledgers

`index396_fibre36.csv`, `two_mask_net36.csv`, and
`alternate_node_union72.csv` are lexicographically ordered row ledgers.
`index396_rebuilt421.txt` and
`two_mask_rebuilt421.txt` are the two ordered decks.  Their FNVs are
`1b2cd4f2728db49a` and `c9ac86709cda10df` respectively.

`SHA256SUMS` is the authoritative raw-LF byte ledger.  Canonical sources live
in

```text
04-computation/lrc14_signature_response_congruence_thm4286/
```

and `REPRODUCTION.md` gives the complete commands.

## Scope

- The signature-purity census is derived from audited THM-4282 ledgers; it
  does not recompute their carrier/deck certificates.
- All new common-row claims have a primary endpoint-cocycle calculation and a
  detached literal-wall check.  Body coverage is global for each named deck.
- Exact minimum three applies to the named 20-obligation fibre-wide response
  problem.  It is not a global minimum over all possible deck changes.
- The current residual is inherited unchanged from THM-4283 and is still a
  fixed-pool proof graph.  No arrival, owner, phase, physical-entry, or
  LRC(14) conclusion follows.
