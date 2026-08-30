# THM-4283 endpoint-638 carrier-response packet

Status: **PROVED RELATIVE TO THM-4282 + FINITE-EXACT + DETACHED
LITERAL-WALL AUDIT PASS.** This packet proves 63 additional fixed-pool
carrier rows. It proves neither physical entry nor LRC(14).

## Consequence

Append

```text
014c9084 = {15,42,85,120,143,145,176,193}
```

to THM-4282's ordered 8,996-mask carrier. The resulting carrier has
count/FNV `8,997/8e1860a25d0fcf87` and closes every row of
`closed63.csv`. The exact set consists of all post-THM4282 rows at endpoints
`>=639` and eight of the nine endpoint-638 rows. Layer counts are

```text
644:7, 643:14, 642:9, 641:7, 640:14, 639:4, 638:8.
```

Its removal from the 23,373-row THM-4282 residual leaves 23,310 rows, FNV
`89ea6a588da4ba0a`, SHA-256
`75c5881e9de9622c1627de4da07dca1df0be82f3366ef1e7eb78e36ff0f71d14`,
with unique top row `(256,638)`.

## Exact components

- `endpoint638_carrier_response_primary.cpp` audits the two decisive response
  quotients. Before insertion, `(256,644)` has exactly two failed bodies and
  367 active full responders; `014c9084` is the least, so the exact minimum
  is one. After insertion it has zero failures. The first resistant row
  `(256,638)` has 40 failures, 315 complete active-response classes, no full
  responder, and exact minimum nine, with an explicit witness.
- `endpoint638_activity_margins.cpp` computes exact endpoint-cocycle masses
  for all seven former top rows and the resistant row. The repair and all
  nine proposed next-stage masks have strictly positive margins where used.
- `endpoint638_detached_literal_audit.cpp` reconstructs the complete
  direct-wall geometry independently of the primary response quotient and
  scans all 14,307,150 labelled nine-bodies for each of the 63 claimed rows.
  It reports zero failures and ledger `cbe0c99a6d552e23`.
- `proof_graph_consequence.py` pins the inherited carriers, repair support,
  exact closed set, layer counts, complement, and unique new top.

## First resistant control

The row `(256,638)` is deliberately excluded from `closed63.csv`. Under the
8,997 carrier it has

```text
active carrier                  3,304  (FNV b59e6995074df7e1)
active/inactive joint           316/105
exposed bodies                  16,792
failed bodies                   40     (FNV 917d107c4536efc9)
complete response classes       315
one-mask full responders        0
exact replacement minimum       9
```

One exact minimum-nine witness is

```text
02203226 081e1084 08a89440 180a8281 18261042
18a0d040 1a82a200 202a9440 280a0a88
```

This witness is a pair-local continuation proposal. It is not appended to
the global carrier and contributes no edge in THM-4283.

## Authoritative files

`SHA256SUMS` is the raw-LF ledger. The four maintained sources live in

```text
04-computation/lrc14_endpoint638_carrier_response_thm4283/
```

and the consequence-bearing transcripts and exact closed set live beside
this manifest. `REPRODUCTION.md` supplies canonical commands.

## Scope

All conclusions are relative to THM-4282's fixed pool, threshold, carrier,
and exact residual. Carrier closure is not a common-deck theorem, the
minimum-nine claim is only for the named `(256,638)` active-response
quotient, and neither statement proves physical exact-`M=12` entry.
