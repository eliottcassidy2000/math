---
id: THM-4313
title: "LRC(14) endpoint-592 forty-three-response size-preserving exchange"
status: >
  PROVED RELATIVE TO THM-4311 + FINITE-EXACT + COMPLETE TWO-RANK RESPONSE
  CENSUS, EXPLICIT RESPONSE PACKING, FULL-467 SINGLE-DELETION QUOTIENT, AND
  DIRECT RAW REPLAY PASS. The inherited 3,925-mask carrier has 2,468
  endpoint-592 failures. An explicit 43-mask response cover and 43 old-mask
  deletions give another 3,925-mask carrier closing all 467 inherited rows.
  The complete response-cover minimum is between 8 and 43; the stronger lower
  bound 36 is proved only inside a specified 37,497-mask retained pool and is
  refuted as a global consequence by exact complete-universe dual pricing.
  The typed union has 2,087 rows and leaves 20,560, maximum endpoint 591 on
  13 rows. No physical entry or LRC(14) follows.
source: root + endpoint592_capacity + endpoint592_cover_opt + endpoint592_structure / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4309-lrc14-endpoint-595-support-threshold-residual-hypergraph-compression
  - THM-4310-lrc14-endpoint-594-residual-layer-closure-and-single-deletion-boundary
  - THM-4311-lrc14-endpoint-593-one-response-size-preserving-exchange
related:
  - THM-4282-inactive-signature-deck-surgery-endpoint-663
  - THM-4302-lrc14-endpoint-596-response-minimum-four-and-size-preserving-exchange
  - THM-4305-lrc14-endpoint-595-pair-tagged-response-exchange
artifact_root: 05-knowledge/results/lrc14_endpoint592_fortythree_response_exchange_thm4313
artifact_manifest: 05-knowledge/results/lrc14_endpoint592_fortythree_response_exchange_thm4313/SHA256SUMS
artifact_manifest_sha256: 72d6c40ee06439f06ee7fc17a67d724ae47b33d1448189380056765d4e7fccd4
primary_scripts:
  - 04-computation/lrc14_endpoint592_fortythree_response_exchange_thm4313/endpoint592_exchanged_carrier_scout.cpp
  - 04-computation/lrc14_endpoint592_fortythree_response_exchange_thm4313/endpoint592_response_multiplicity.cpp
  - 04-computation/lrc14_endpoint592_fortythree_response_exchange_thm4313/verify_global_packing8.cpp
  - 04-computation/lrc14_endpoint592_fortythree_response_exchange_thm4313/endpoint592_cover_direct_audit.cpp
  - 04-computation/lrc14_endpoint592_fortythree_response_exchange_thm4313/endpoint592_retained_pool_export.cpp
  - 04-computation/lrc14_endpoint592_fortythree_response_exchange_thm4313/certify_pool_dual_integer.py
  - 04-computation/lrc14_endpoint592_fortythree_response_exchange_thm4313/full_universe_integer_dual_pricing.cpp
  - 04-computation/lrc14_endpoint592_fortythree_response_exchange_thm4313/endpoint592_full467_singleton_capacity.cpp
  - 04-computation/lrc14_endpoint592_fortythree_response_exchange_thm4313/endpoint592_full467_exchange_audit.cpp
  - 04-computation/lrc14_endpoint592_fortythree_response_exchange_thm4313/typed_endpoint592_consumer.py
  - 04-computation/lrc14_endpoint592_fortythree_response_exchange_thm4313/verify_endpoint592_packet.py
audit: >
  PASS / ACCEPT, subject to the frozen SHA256SUMS. Baseline, response
  multiplicity, packing, direct-cover, singleton, complete-pricing, and
  exchange outputs agree byte-for-byte at O2/O3. The final carrier is replayed
  directly on all 6,681,439,050 row-body cases. Exact full-universe pricing
  records the 20,986 omitted-mask violations that prevent the retained-pool
  dual lower bound from being misreported globally.
---

# THM-4313 -- LRC(14) endpoint-592 forty-three-response size-preserving exchange

**PROVED RELATIVE TO THM-4311 + FINITE-EXACT + COMPLETE RESPONSE CENSUS,
EXPLICIT PACKING, FULL-467 SINGLETON QUOTIENT, AND DIRECT RAW REPLAY PASS.
LRC(14) REMAINS OPEN.**

## 1. Hostile endpoint-592 layer

Retain THM-4311's exchanged carrier

~~~text
C_old: size=3,925, rank8=3,857, rank9=68,
FNV=c9e5faef52ca5707, all 421 joint masks retained.
~~~

The next typed frontier consists of 35 endpoint-592 rows, ordered FNV
3eb23833c35b9266. Direct O3 and O2 scans of all

~~~text
35*binom(30,9)=500,750,250
~~~

cases agree byte-for-byte and find 2,468 failures on exactly seven rows:

| q | 96 | 100 | 105 | 192 | 210 | 256 | 416 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| failures | 13 | 3 | 48 | 1 | 288 | 2,101 | 14 |

Their ordered obligation FNV is 2209b8d6760280cc.

## 2. Response cover bounds

Complete enumeration of active, disjoint rank-eight and rank-nine masks gives
188,462 rank-eight responders, FNV bb1ff3319125eaae, and 2,205,776 rank-nine
responders, FNV ec506c87a74e64e1. Every obligation has between 423 and 34,654
responses.

The packet contains eight pair-tagged obligations such that no active
rank-eight or rank-nine mask responds to any pair. The verifier scans all
20,160,075 masks and finds zero co-responses on all 28 pairs. Therefore every
complete-universe response cover has at least eight masks.

The explicit addition ledger A has

~~~text
|A|=43, rank8=2, rank9=41, FNV=ca3cb80f471f2e7e.
~~~

A direct arithmetic replay covers all 2,468 failures. Thus the complete
rank-eight/rank-nine response-cover minimum lies between 8 and 43.

A separately specified 37,497-mask retained pool has an exact integer dual
lower bound 36. This is pool-relative only. Complete-universe pricing finds
20,986 omitted masks violating the dual, with maximum score 3849498499 at
mask 20188724 against scale 1000000000. Hence no global lower bound 36 is
claimed.

## 3. Singleton capacity and selected deletion ledger

Use the full preservation target

~~~text
K_467 = 391 inherited prefix rows
        union 25 complete endpoint-594 rows
        union 16 endpoint-593 rows
        union 35 endpoint-592 rows,
|K_467|=467, ordered FNV=2d6aa988098aa5eb.
~~~

After adjoining A, the carrier has size 3,968 and FNV b8221545b2d668d0.
Filtering THM-4311's exact inherited singleton ledger and raw-scanning the 35
new rows gives 1,510 private obligations, FNV 4dfcadbd1d5c0c31. Exactly 412
masks are protected: 369 old carrier masks and all 43 additions. There are
3,135 individually safe old nonjoint masks.

Let D be the first 43 safe masks in the deterministic order by active-row
count, rank, and mask. Then

~~~text
|D|=43, FNV=0dd14ef0fe3eec62.
~~~

The singleton quotient does not imply that D is simultaneously deletable.

## 4. Complete size-preserving exchange

Define

~~~text
C_new=(C_old minus D) union A.
~~~

Then

~~~text
|C_new|=3,925, rank8=3,818, rank9=107,
FNV=a0d08a38c10bdab7, all 421 joint masks retained.
~~~

A separate raw program reconstructs C_new and visits every one of

~~~text
467*binom(30,9)=6,681,439,050
~~~

row-body cases. O3 and O2 transcripts, 467-row pair ledgers, and empty
failure ledgers are byte-identical. They give zero failures and pair-ledger
FNV 333e401fef6c7240. This direct replay, not the singleton quotient, proves
the simultaneous 43-for-43 exchange.

## 5. Typed consequence and scope

Consuming the 35 endpoint-592 rows produces

~~~text
typed union: 2,087 rows, FNV=23e4136827b770a5;
residual:   20,560 rows, FNV=8d797592a729e0b3;
next frontier: endpoint 591 on 13 rows, FNV=fc332c0697c671c7.
~~~

All claims are finite-exact on the fixed labelled 30-speed pool and the fixed
467-row target. The theorem does not prove an exact response-cover minimum, a
globally minimum carrier, a physical entry, an arbitrary-pair statement, a
terminating descent, or LRC(14).
