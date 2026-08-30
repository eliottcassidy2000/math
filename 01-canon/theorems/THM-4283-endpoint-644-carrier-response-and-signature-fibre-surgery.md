---
id: THM-4283
title: "One-mask endpoint-644 repair, carrier descent to the endpoint-638 singleton, and the exact nine-response wall"
status: >
  PROVED RELATIVE TO THM-4282 + FINITE-EXACT + DETACHED LITERAL-WALL AUDIT
  PASS. Appending one explicit active rank-eight mask to the inherited
  8,996-mask carrier closes exactly 63 stated post-THM4282 rows: every row at
  endpoints 644 through 639 and eight of nine rows at 638. The exact
  complement has 23,310 rows with unique top (256,638). That row has 40
  carrier failures, no one-mask responder, and exact pair-local replacement
  minimum nine. The nine-mask witness is not appended. No physical entry or
  LRC(14) follows.
source: codex-continuation-frontiers-20260829
depends_on:
  - THM-4282-inactive-signature-deck-surgery-endpoint-663
related:
  - THM-4281-rectangle-common-joint-deck-endpoint-670-bridge
  - THM-4278-endpoint-520-688-minimum-one-atom-carrier-augmentation
artifact_root: 05-knowledge/results/lrc14_endpoint638_carrier_response_thm4283
artifact_manifest: 05-knowledge/results/lrc14_endpoint638_carrier_response_thm4283/SHA256SUMS
artifact_manifest_sha256: 5e4f27ac7854121dcde44f51a1fab076848d8feb15bd2e30a168b295543162d2
primary_scripts:
  - 04-computation/lrc14_endpoint638_carrier_response_thm4283/endpoint638_carrier_response_primary.cpp
  - 04-computation/lrc14_endpoint638_carrier_response_thm4283/endpoint638_activity_margins.cpp
  - 04-computation/lrc14_endpoint638_carrier_response_thm4283/endpoint638_detached_literal_audit.cpp
  - 04-computation/lrc14_endpoint638_carrier_response_thm4283/proof_graph_consequence.py
primary_script_sha256:
  - 1c469ba1ec67f25e982a3c02b0a45fca90d9434169a58f3cff9d1eba7da32e9d
  - d744f71c995490610c16b050f707e5d0673042a84ec281621cffe58db52976cb
  - fe9ca902d74e995292ee0b4b1751cd9fa39c2e8cc2727cdca639d8fb5ca58ae1
  - cb964e4770063f71a9734ee790fd9822a69c38b070a6d31136b1b3c2a6efc3be
primary_output_sha256:
  - ef8036de9cf266d1fa899d857cb90b9031284e4224fe28bbc15623e55fc7bdb4
  - f0f93a5624fd1d3cc1b19fc6d7b12efd9ba5e1d252566d518beb4a59ea325363
  - 55e0d6d44b93a95878ba3fd161e96550e2ed0eab51e3fe2cce3e1048911102aa
  - 418db367b18dc70d246cab8fb83af8b9034284e94b92dc6884e11523d412db2a
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT, subject to the frozen SHA256SUMS. The primary response scan
  computes complete active response quotients before and after the repair and
  proves exact minima one and nine. An independent direct-wall scanner
  reconstructs every active carrier and scans all 14,307,150 labelled bodies
  on each of the 63 claimed rows, with zero failures. Exact cocycle margins
  are strict, and an optimization-safe Python replay pins the carrier, closed
  set, residual hashes, and unique top.
---

# THM-4283 -- one-mask endpoint-644 repair, descent to endpoint 638, and the exact nine-response wall

**PROVED RELATIVE TO THM-4282 + FINITE-EXACT + DETACHED LITERAL-WALL AUDIT
PASS. LRC(14) REMAINS OPEN.**

## 1. Statement

Retain THM-4282's fixed pool, threshold `4/63`, frozen 421-mask joint deck,
ordered 8,996-mask carrier, and exact 23,373-row residual `U_4282`. Append
the rank-eight mask

```text
w=014c9084={15,42,85,120,143,145,176,193}.               (1)
```

The ordered carrier identities are

```text
before w:     count=8,996,       FNV=fd899660f14b311c;
after w:      count=8,997,       FNV=8e1860a25d0fcf87.    (2)
```

The augmented carrier closes every row in

```text
C={(q,r) in U_4282:r>=639}
  union {(q,638):q in {100,282,294,366,416,420,512,520}}. (3)
```

The exact layer counts are

```text
r=644  643  642  641  640  639  638
count   7   14    9    7   14    4    8,                (4)
```

so `|C|=63`. In lexicographic row order,

```text
FNV(C)=3ed7002531effa14,
SHA256(C)=2c48a1b46f3ac2fce4b9fdaf426373288bd2ee15649c228de2a8e4d9fc6f46fd.
                                                               (5)
```

The exact complement has

```text
|U_4282-C|=23,310,
FNV=89ea6a588da4ba0a,
SHA256=75c5881e9de9622c1627de4da07dca1df0be82f3366ef1e7eb78e36ff0f71d14,
maximum endpoint=638,
top={(256,638)}.                                        (6)
```

The unique top row is a genuine local wall for this repair. Under the
8,997-mask carrier, `(256,638)` has 40 failed bodies. Its complete active
response quotient has no one-mask full responder and has exact replacement
minimum nine. One minimum witness is

```text
02203226 081e1084 08a89440 180a8281 18261042
18a0d040 1a82a200 202a9440 280a0a88.                   (7)
```

The masks in `(7)` are not appended to the global carrier and contribute no
new row in this theorem.

## 2. Exact endpoint-644 repair

Before adding `w`, the only resistant former top row is `(256,644)`. Its
complete primary audit gives

```text
active carrier masks       4,124,   FNV c277f65c86b786a6
active/inactive joint      400/21,  inactive FNV 3b44cf197cffa1e4
joint-exposed bodies       235,     FNV 9d3ecac607ebedc9
failed bodies              2,       FNV 936a8b25300381b7
failed masks               14326401, 1c306401.          (8)
```

Complete evaluation of all `5,852,925` rank-eight masks, with activity at
the named pair retained, produces four response classes. Exactly 367 active
masks hit both failed bodies. The least is `w`; therefore the replacement
minimum is exactly one. After appending `w`, the active carrier has
count/FNV `4,125/f07d2c694c9621e1` and zero failed bodies. The other six
endpoint-644 rows already had zero failures under the prior carrier.

The endpoint-cocycle margin of `w` at the decisive row is strictly positive:

```text
D=10525586771134955520,
63 mass(w)-4D=150690435472027584,
63 mu(w)-4=3951953041/276040244480,
mu(w)-4/63=3951953041/17390535402240>0.                 (9)
```

The detached direct-wall convention gives grid
`6712746665264640` and integer margin `96103594051038`; reduction gives the
same rational value in `(9)`. Exact margins of `w` are positive on every row
in `(3)`.

## 3. Detached 63-row carrier proof

The primary response engine uses the frozen joint deck as a fast exposed-body
router. The consequence itself is audited by a separate literal engine. For
each row of `C`, it

1. rebuilds the exact pair wall geometry and all active carrier masks;
2. enumerates all `14,307,150` labelled nine-bodies independently; and
3. tests direct disjoint coverage by the active 8,997-mask carrier.

All 63 rows have zero failures. The literal active counts and FNVs agree with
the primary geometry where both are evaluated, and its complete transcript
ends with

```text
LEDGER cbe0c99a6d552e23 TOTAL_FAILURES 0.                (10)
```

Thus `(3)` is a carrier certificate. It is not promoted to a common-deck
statement.

## 4. The first resistant response quotient

At `(256,638)`, the repaired carrier and joint router give

```text
active carrier masks       3,304,   FNV b59e6995074df7e1
active/inactive joint      316/105, inactive FNV 692415a9cc19d4a3
joint-exposed bodies       16,792,  FNV c04e9eb475997f8e
carrier failures           40,      FNV 917d107c4536efc9
positive hit incidences    624,748
positive hit range         1..240
response ledger FNV        3fafa3082852ffcc.            (11)
```

The complete active response quotient on those 40 labelled obligations has

```text
315 classes,              class FNV a78d55f5a8e5d73,
active-response FNV       9ae5db161dbddd21,
full one-mask classes     0.                              (12)
```

Discarding a response pattern contained in another cannot increase cover
size: its mask can be replaced by the active mask representing the containing
pattern. Breadth-first search over unions of the remaining maximal response
patterns first reaches all 40 obligations at depth nine. Backtracking gives
the explicit witness `(7)`. Hence nine is an exact pair-local minimum, not a
greedy upper bound.

All nine masks in `(7)` have strictly positive exact activity margin at
`(256,638)`; the smallest displayed numerator is
`20822587869504960>0`. This preserves the activity coordinate that a bare
set-cover calculation would lose.

## 5. Exact proof graph and reproduction

The exact closed set is stored at

```text
05-knowledge/results/lrc14_endpoint638_carrier_response_thm4283/closed63.csv.
```

`proof_graph_consequence.py` pins THM-4282's residual and carrier byte
identities, checks the predicate `(3)`, computes `(4)--(6)`, and verifies the
support in `(1)`. The exact identity is

```text
U_4283=U_4282-C.                                         (13)
```

No overlap correction with THM-4282 is needed because `C` is selected inside
the already final post-THM4282 residual.

The complete packet, byte ledger, and commands are in

```text
05-knowledge/results/lrc14_endpoint638_carrier_response_thm4283/
  MANIFEST.md
  REPRODUCTION.md
  SHA256SUMS.
```

## 6. Scope and connection ledger

The typed connection is

```text
source:              THM-4282's ordered carrier and exact residual;
target:              pair-local labelled nine-body cover obligations;
map:                 inactive joint signature -> exposed bodies ->
                     complete active response patterns;
preserved predicate: exact activity and disjoint body coverage;
destroyed data:      physical arrival, common-deck interchangeability,
                     and response validity at another pair;
restoring sidecars:  pair, layer, margins, carrier order, full body ledger;
sharp hostiles:      the two bodies at (256,644) before w and the forty
                     bodies at the new singleton (256,638).               (14)
```

This theorem proves no common deck on `C`, no minimum-nine global carrier
augmentation, no physical exact-`M=12` entry, and no LRC(14) conclusion.
