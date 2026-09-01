---
id: THM-4306
title: "LRC(14) index-265 recursive ideal two-mask replacement"
status: >
  PROVED RELATIVE TO THM-4296 AND THM-4305 + FINITE-EXACT + INDEPENDENT
  DIRECT-ATOM AUDIT PASS. The 367-row singleton-signature ideal H265 has no
  one-mask common-active rank-eight replacement, but has exact replacement
  minimum two. Rebuilding the 421-mask joint deck with two masks closes all
  367 rows. Relative to T4305, the ideal adds 353 typed rows; the combined
  union has 2,014 rows and leaves 20,633, maximum endpoint 594 on 22 rows. No
  physical entry or LRC(14) follows.
source: root / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4296-lrc14-mixed-rank-deletion-depth-and-recursive-signature-closure
  - THM-4305-lrc14-endpoint-595-pair-tagged-response-exchange
related:
  - THM-4300-lrc14-size-preserving-response-staircase-and-index-297-ideal
  - THM-4303-lrc14-endpoint-595-twenty-five-row-carrier-closure
artifact_root: 05-knowledge/results/lrc14_index265_recursive_ideal_thm4306
artifact_manifest: 05-knowledge/results/lrc14_index265_recursive_ideal_thm4306/SHA256SUMS
artifact_manifest_sha256: 89100e9d7ecc74c1f2a5b7fc373d50da8c95759c27d4592190e5e77acd2e67e1
primary_scripts:
  - 04-computation/lrc14_index265_recursive_ideal_thm4306/index265_response_quotient.cpp
  - 04-computation/lrc14_index265_recursive_ideal_thm4306/index265_direct_audit.cpp
  - 04-computation/lrc14_index265_recursive_ideal_thm4306/typed_union_consumer.py
audit: >
  PASS / ACCEPT, subject to the frozen SHA256SUMS. The primary computation
  traverses the complete common-active rank-eight response quotient. A
  structurally independent literal-atom path exhausts the one-mask hostile,
  audits every retained/witness sign, and scans the complete body universe.
  Its O0/O3 transcripts agree. The typed consumer audits the quotient and
  reconstructs every inherited and enlarged row partition.
---

# THM-4306 -- LRC(14) index-265 recursive ideal two-mask replacement

**PROVED RELATIVE TO THM-4296 AND THM-4305 + FINITE-EXACT + INDEPENDENT
DIRECT-ATOM AUDIT PASS. LRC(14) REMAINS OPEN.**

## 1. Inherited object and target ideal

Retain THM-4296's ordered 421-mask rank-eight joint deck `E`, its complete
22,647-row fixed-pool residual universe `U`, and its inactive-signature atlas.
For a row `p`, let `I_E(p)` be the set of inactive joint coordinates.  The
target is the exact singleton-signature ideal

```text
H_265={p in U:I_E(p)={265}}.                              (1)
```

This is the largest failed singleton group in the complete THM-4296 census:

```text
|H_265|=367,
FNV=d422b161d94ebae4,
SHA256=750a31e2f0ebe6573835cf9d2cd43e83403d94fe802616e539ee9d330a3dab65,
maximum endpoint=626.                                     (2)
```

Its inactive deck coordinate is

```text
E_265=20820649.                                           (3)
```

Deleting `(3)` exposes exactly eight labelled nine-bodies:

```text
B0=09392104 B1=0d30e080 B2=0d382104 B3=15386080
B4=186c9080 B5=19786000 B6=1d489080 B7=1f087000,
FNV=4ebad120ec162881.                                     (4)
```

The other 420 joint masks remain active at every row of `(1)`.  A replacement
mask must likewise be common-active on all 367 rows.  For such a rank-eight
mask `R`, define

```text
rho(R)={i in {0,...,7}:R intersect B_i=empty}.             (5)
```

The map `(5)` preserves precisely which newly private bodies are repaired. It
forgets literal margins and first blocker rows; both are retained by the
independent direct audit.

## 2. Complete response quotient

The primary computation intersects exact activity over all

```text
367 * binom(30,8)                                         (6)
```

mask-row cells.  Its complete common-active family and response quotient are

```text
common-active masks=1,494,889, FNV=62bfab5b7eafb1d2,
nonempty responders=57,752,
realized response classes=38,
full responders=0,
maximal response antichain={0x7f,0xa5,0xd0}.              (7)
```

Thus the full-responder-only route genuinely fails.  The complete quotient,
including zero response, counts, least realizing masks, and maximal flags, is
frozen in `results/index265_response_quotient.csv`.

## 3. Exact replacement minimum two

Give unit packing weight to `B1` and `B7`, and zero to the other six bodies.
Every one of the 38 realized responses has packing load at most one. Hence no
single common-active mask can cover all eight obligations.

The independent path makes this lower certificate literal.  Any mask covering
both `B1` and `B7` must be an eight-subset of the seventeen-label complement
of

```text
B1 union B7=1f38f080.                                    (8)
```

It exhausts all `binom(17,8)=24,310` candidates, finds a first blocker in
`H_265` for every one, and finds zero common-active survivors.

The two common-active masks

```text
22020e09 : response 0x7f,
00868489 : response 0xa5                                (9)
```

have response union `0xff`.  The packing lower bound and cover `(9)` prove
that the exact common-active rank-eight replacement minimum is two.

## 4. Rebuilt-deck audit

Replace `(3)` by `(9)`.  The rebuilt common deck has

```text
size=422,
FNV=12980884e0346f9f.                                    (10)
```

The independent program derives `(1)` directly from signature atoms, checks
all `367*420=154,140` retained-mask signs and all `367*2=734` replacement
signs, and finds no equality or inactive cell.  It then scans all

```text
binom(30,9)=14,307,150                                    (11)
```

labelled bodies against `(10)` and finds zero failures.  Its O0 and O3
transcripts agree byte-for-byte.  Thus `(10)` is a separate common-deck
certificate for every row in `H_265`.

## 5. Typed row consequence

Only the proved row set `(1)` is unioned; the deck `(10)` is not identified
with any carrier or another proof node.  The typed consumer independently
reconstructs THM-4303 and THM-4305.  It finds

```text
|H_265 intersect T_4305|=14,
|H_265 \ T_4305|=353.                                    (12)
```

The three endpoint-595 rows added by THM-4305 are disjoint from `H_265`.
Consequently

```text
|T_4305 union H_265|=2,014,
FNV=ee275cd5d460d153,
SHA256=836b7ccd0f93268ae039d749c6d778778dd0fcf706dd58cb6dfbd3c340fcbcd1,

|U \ (T_4305 union H_265)|=20,633,
FNV=3acd694e62cb7841,
SHA256=77e54cc8614f2750a7e3a46fd9ec5a3a9b4f86db38121d62802d7387934e7e7f.
                                                                    (13)
```

The residual maximum is endpoint 594 on exactly 22 rows:

```text
(96,594),  (100,594), (105,594), (147,594), (186,594),
(192,594), (206,594), (210,594), (220,594), (244,594),
(256,594), (260,594), (282,594), (294,594), (308,594),
(313,594), (315,594), (366,594), (416,594), (440,594),
(462,594), (520,594),
FNV=8413f0d2282e4cd6,
SHA256=2a46ac360974ee95b5c468f1f76fb9ddd6b5165fa6e410dd3b6bad02ca93dd54.
                                                                    (14)
```

## 6. Scope firewall

- The minimum two concerns common-active rank-eight replacements on the exact
  singleton ideal `I_E(p)={265}` only.
- The result does not claim a uniform replacement bound for other signatures.
- The rebuilt deck is a separate row certificate.  It is not merged with the
  endpoint carrier; only row consequences are unioned.
- No terminating recursive-ideal descent, physical-entry construction,
  arbitrary-pair theorem, semantic-arrival theorem, or proof of LRC(14)
  follows.
