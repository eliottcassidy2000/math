---
id: THM-4311
title: "LRC(14) endpoint-593 one-response size-preserving exchange"
status: >
  PROVED RELATIVE TO THM-4309/4310 + FINITE-EXACT + COMPLETE TWO-RANK
  RESPONSE CENSUS, FULL-432 SINGLE-DELETION QUOTIENT, AND DIRECT RAW REPLAY
  PASS. The inherited 3,925-mask carrier has one endpoint-593 failure. Adding
  one rank-nine response and deleting one certified-safe old rank-eight mask
  gives another 3,925-mask carrier closing all 432 inherited rows. The typed
  union has 2,052 rows and leaves 20,595, maximum endpoint 592 on 35 rows. No
  physical entry or LRC(14) follows.
source: root + endpoint595_response / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4309-lrc14-endpoint-595-support-threshold-residual-hypergraph-compression
  - THM-4310-lrc14-endpoint-594-residual-layer-closure-and-single-deletion-boundary
related:
  - THM-4302-lrc14-endpoint-596-response-minimum-four-and-size-preserving-exchange
  - THM-4305-lrc14-endpoint-595-pair-tagged-response-exchange
artifact_root: 05-knowledge/results/lrc14_endpoint593_response_exchange_thm4311
artifact_manifest: 05-knowledge/results/lrc14_endpoint593_response_exchange_thm4311/SHA256SUMS
artifact_manifest_sha256: e697c4edc3b6612cc30c7800e005568ed554b6c5ad69762019925c05783e2655
primary_scripts:
  - 04-computation/lrc14_endpoint593_response_exchange_thm4311/endpoint593_carrier_audit.cpp
  - 04-computation/lrc14_endpoint593_response_exchange_thm4311/endpoint593_response_capacity.cpp
  - 04-computation/lrc14_endpoint593_response_exchange_thm4311/independent_complement_response_audit.cpp
  - 04-computation/lrc14_endpoint593_response_exchange_thm4311/endpoint593_singleton_deletion_quotient.cpp
  - 04-computation/lrc14_endpoint593_response_exchange_thm4311/endpoint593_exchange_full432_audit.cpp
  - 04-computation/lrc14_endpoint593_response_exchange_thm4311/typed_endpoint593_consumer.py
  - 04-computation/lrc14_endpoint593_response_exchange_thm4311/verify_endpoint593_packet.py
audit: >
  PASS / ACCEPT, subject to the frozen SHA256SUMS. Baseline, singleton, and
  exchange programs agree byte-for-byte at O2/O3; a complement-subset
  implementation independently reproduces the response census. The final
  carrier is replayed directly on all 6,180,688,800 row-body cases. A
  pre-promotion audit repaired the hostile baseline's exit status and replayed
  both builds without changing any frozen mathematical output.
---

# THM-4311 -- LRC(14) endpoint-593 one-response size-preserving exchange

**PROVED RELATIVE TO THM-4309/4310 + FINITE-EXACT + COMPLETE TWO-RANK
RESPONSE CENSUS, FULL-432 SINGLE-DELETION QUOTIENT, AND DIRECT RAW REPLAY PASS.
LRC(14) REMAINS OPEN.**

## 1. Inherited objects and exact hostile boundary

Retain THM-4309's carrier

```text
C_old: size=3,925, rank8=3,858, rank9=67,
FNV=6fbd0bffcf0ed78b, all 421 joint masks retained.       (1)
```

THM-4310 leaves exactly 16 endpoint-593 rows, with ordered row FNV
`5424c07fa724011f`. Direct O3 and O2 scans of all

```text
16*binom(30,9)=228,914,400                               (2)
```

row-body cases agree byte-for-byte and find one failure:

```text
(q,r,body)=(96,593,34087401),
obligation FNV=643cc006d08e9a83.                         (3)
```

For exchange preservation, form the disjoint row union

```text
K_432 = THM-4309's 391-row target
        union THM-4310's complete 25-row endpoint-594 layer
        union the 16 endpoint-593 rows.
|K_432|=432, ordered FNV=a7ed492c64d1c0d8.               (4)
```

The complete 25-row layer in `(4)` is essential: it includes the three rows
that were already typed before THM-4310. A preliminary 429-row scout omitted
those rows and was quarantined before promotion.

## 2. Exact response minimum one and zero strict capacity

An active mask responds to `(3)` exactly when it is disjoint from body
`34087401`. Complete enumeration gives

| rank | all masks | active | responders | responder FNV | least |
|---:|---:|---:|---:|---:|---:|
| 8 | 5,852,925 | 1,262,283 | 1,636 | `56f82f5dc11db83b` | `0134012c` |
| 9 | 14,307,150 | 7,218,133 | 16,209 | `0f615d806860553f` | `0036092c` |

Because `C_old` fails `(3)`, zero additions cannot suffice; mask `0036092c`
is an active disjoint witness, so the exact rank-eight/rank-nine response
minimum is one. A separate complement-subset enumeration reproduces both
responder sets. Put

```text
A={0036092c}, rank9, singleton FNV=60873ef7a2b4ab90.      (5)
```

No mask of `C_old` is inactive on every row of either the prior 416-row target
or `K_432`. The exact `3,925*432=1,695,600` sign census has no equality and FNV
`c0a847088f7c904c`. Thus the strict-common-inactivity mechanism cannot pay for
the addition in `(5)`.

## 3. Exact one-deletion quotient

Let `C_plus=C_old union A`. It has 3,926 masks: all 421 joint masks and 3,505
nonjoint masks. Hold the joint deck fixed. For every row-body obligation on
`K_432`, enumerate the active disjoint nonjoint witnesses whenever the joint
deck does not already close the body. A nonjoint mask is unsafe to delete iff
it is the unique such witness of at least one obligation. The complete quotient
over all `432*binom(30,9)` cases gives

```text
private nonjoint obligations=1,520, FNV=39fedd8ceb347304,
protected nonjoint masks=412,       FNV=7ee1d68a078d5b65,
safe nonjoint masks=3,093.                                (6)
```

The new response in `(5)` is protected by exactly the hostile obligation
`(3)`, a negative control. The least safe old-carrier deletion is

```text
D={0006e281}, rank8, singleton FNV=4c14214a64ec202c.      (7)
```

O2 and O3 reproduce all 1,520 private obligations, all 412 protected masks,
and `(7)` byte-for-byte. Equation `(6)` is a one-deletion theorem only; it
does not certify simultaneous deletions.

## 4. Same-size exchange and full raw replay

Define

```text
C_new=(C_old\D) union A.
|C_new|=3,925, rank8=3,857, rank9=68,
FNV=c9e5faef52ca5707, all 421 joint masks retained.       (8)
```

A separate program reconstructs `(8)` from canonical component ledgers and
directly scans every one of

```text
432*binom(30,9)=6,180,688,800                            (9)
```

row-body cases. Both builds give zero failures, 1,509,470 joint-exposed
bodies, 50,693,609 nonjoint hit incidences, and pair-ledger FNV
`48ad5090055eeeae`. Hence `(8)` closes the exact row target `(4)`.

## 5. Typed consequence and scope

Consuming the 16 endpoint-593 rows into THM-4310's audited partition gives

```text
typed union: 2,052, FNV=7c3f8bda9c37c5c2,
residual:   20,595, FNV=0a9a532153a5e6dc.                (10)
```

The residual maximum is endpoint 592 on exactly 35 rows, ordered FNV
`3eb23833c35b9266`.

- The response minimum is restricted to rank-eight/rank-nine masks for the
  single frozen failure `(3)`.
- The deletion quotient holds the joint deck fixed and concerns exactly one
  nonjoint deletion from `C_plus`.
- The carrier acts only on the inherited labelled thirty-speed pool. Carrier
  failures are proof obligations, not physical danger witnesses.
- No global carrier minimum, arbitrary-pair theorem, terminating descent,
  physical entry, or proof of LRC(14) follows.

The next exact frontier test is the direct replay of `(8)` on the 35 rows at
endpoint 592. Its response quotient, rather than a forced tournament model,
is the relevant sidecar if failures remain.
