---
id: THM-4283
title: "Endpoint-644 carrier response and signature-fibre surgery"
status: >
  PROVED RELATIVE TO THM-4282 + FINITE-EXACT + DETACHED LITERAL-WALL
  AUDITS PASS. A one-mask repair closes the exact 55-row endpoint band
  639..644, and an exact minimum-nine boundary repair extends the same carrier
  to all 64 rows at 638..644. All 127 nonempty top-signature deletions,
  supplemented by four target-aware covers, realize the complete 640-row
  signature fibre. Their typed union removes 691 post-THM-4282 rows and leaves
  22,682, maximum endpoint 637. No physical entry or LRC(14) follows.
source: root/cross-frontier-bridge/2026-08-29; strengthens codex-continuation-frontiers-20260829
depends_on:
  - THM-4282-inactive-signature-deck-surgery-endpoint-663
related:
  - THM-4281-rectangle-common-joint-deck-endpoint-670-bridge
  - THM-4278-endpoint-520-688-minimum-one-atom-carrier-augmentation
  - THM-4276-six-atom-endpoint-671-augmentation-and-one-layer-descent
artifact_root: 05-knowledge/results/lrc14_endpoint_carrier_signature_surgery_thm4283
artifact_manifest: 05-knowledge/results/lrc14_endpoint_carrier_signature_surgery_thm4283/SHA256SUMS
artifact_manifest_sha256: 51aab45c9e098740a0f3b9985f8c45cbcdac1a25d8cc9bcc35e4d1281b79d135
primary_scripts:
  - 04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/endpoint_top_band_scan.cpp
  - 04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/endpoint638_exact_response_witness.cpp
  - 04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/all127_response_greedy.cpp
  - 04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/targeted_response_lift.cpp
  - 04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/all127_common_family_audit.cpp
  - 04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/selected12_detached_literal_audit.cpp
  - 04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/carrier_band_detached_literal_audit.cpp
  - 04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/proof_graph_consequence.py
  - 04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/proof_graph_consequence_independent.rb
primary_script_sha256:
  - 597cedcc15c1f4c07536016ec76c06f5fa4ebe8dbb367fe9bc2a44bb9eec9c19
  - 05b45364b6d31f33392def8652dc5444f0e378eb04ccbafc32a025f7663a74a9
  - 47c71bed753220a2151d9877a62b796fc26ea4acad2129f780bc4ecf33698e48
  - 0d37c2536a8d5392c9b41d6557d700f75dc7e01320a79734c985b83a6c27a36e
  - c2f768f02ed9787cc53aa50a7e64614b84ccca7c9c6809b26464c5d8eaca74b9
  - e7f9c3a1e9f27c9fc319c4770d3ad36d571e0aae41240002b7c0df138a80d7e2
  - 555fe42c3e28e8f88d0ec3fd363e08d599bfe9e6dbe0c37c1ee6c2d9f058b6fc
  - 0b1652b47a1b5f5df56f9927ccfac49b714700f9d1c3d3c73213b382bf2165f2
  - d9b019e8a02912e6dacd6cb06bf44f9dfa069b8db4747a965398f1f58d6b9331
primary_output_sha256:
  - 0bc692f631be452fe36532cc8d673491965301a65618eb2e887040c13b7a6bfa
  - 1cac6ed858c194fb0f9fe539de8f5b3c27ea87b288d93b8630091e2dab2a2585
  - 8161345b46f81175643296202d8ae3121e4b65a5dd248b8cb28424e189a845cc
  - c6f8e24e93a741b9ce139489110b7f74eadb13add69fcdf2111931f6d3a5d9a7
  - bef5652b5fe6c0a6822bafde32907a61078d8bc752242c759d69ab6dce23635a
  - 18ce15618511925870ae9047dc4d51f15566b0f3b1f97b7f2cb46a48b7911368
  - ef78488e5d356138289217290fff4c799fa1b5aa652cdfe761495e73dd0febf6
  - e67f8acaa398335e7f8213a2eb294f6b532732198c6e6b0be39abd8514448b29
  - 85e3db1e25f2b62a9b53d3ac3041048e2c65968e13c05df71914d3e9eb0e90d3
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT, subject to the frozen top-level SHA256SUMS. The primary
  cocycle/response route and detached literal-wall route agree on the carrier
  band and its first hostile layer. The 127 Boolean scenarios and four
  target-aware variants have exact body-cover replays; twelve selected decks
  have a detached 12*C(30,9) body scan and 1,537,426 literal activity checks.
  O0, O3, NDEBUG, and UBSan controls agree where stated. Exact Python and
  independently written Ruby proof-graph replays agree byte-for-byte. The
  independently produced sibling packet
  `lrc14_endpoint638_carrier_response_thm4283` corroborates the one-mask
  63-row prefix and the exact nine-response wall by a third implementation.
---

# THM-4283 -- endpoint-644 carrier response and signature-fibre surgery

**PROVED RELATIVE TO THM-4282 + FINITE-EXACT + DETACHED LITERAL-WALL
AUDITS PASS. LRC(14) REMAINS OPEN.**

## 1. Statement

Retain THM-4282's labelled thirty-speed pool, threshold, ordered 421-mask
joint deck `E`, final 8,996-mask carrier `C`, and exact residual `U_4282`:

```text
|U_4282|=23,373,                 FNV=c6ab0ae49ee32273,
SHA256=c3e5bf37887aa57af79cb166fce4a6e933e5daffc26dd8032fdfc52ce31240f3,
max endpoint=644,
top={(220,644),(256,644),(258,644),(294,644),
     (366,644),(416,644),(512,644)}.                    (1)
```

There are two complementary conclusions.

1. The carrier `C` fails on exactly two labelled nine-bodies at
   `(256,644)`. Appending the least exact full responder `0x014c9084` gives a
   nested 8,997-mask carrier `C+` that closes every row of `U_4282` with
   endpoint `639..644`, exactly 55 rows. At endpoint 638, only `(256,638)`
   fails, on exactly 40 bodies. Appending an exact minimum-nine response gives
   a 9,006-mask carrier `C++` that closes all nine rows there, hence the exact
   64-row band `638..644`.
2. Let `J_*` be the union of the seven top rows' inactive signatures for `E`.
   Then `|J_*|=27`, and the exact signature fibre

   ```text
   G={p in U_4282:I_E(p) subset J_*}                    (2)
   ```

   has 640 rows. Twelve explicit rebuilt body-covering decks have common-row
   union exactly `G`.

The carrier band `K` in item 1 has 13-row overlap with `G`. Thus the typed
proof-graph-new union has

```text
|G union K|=640+64-13=691.                              (3)
```

Removing `(3)` from `(1)` leaves

```text
count=22,682,
FNV=f7563445f15efebf,
SHA256=7d0044bc4c32f08b9d09dca420171a05666d26e03f38fbc48a9baa1fcb027102,
maximum endpoint=637,
top={(100,637),(294,637),(520,637)}.                     (4)
```

The twelve common decks, `C++`, and their typed set union are different proof
nodes. Equation `(3)` does not assert that one global deck realizes all 691
rows.

## 2. Exact endpoint-644 carrier response

The inherited final carrier has

```text
|C|=8,996,                   FNV=fd899660f14b311c.       (5)
```

Exact joint-exposure scanning of the seven rows at endpoint 644 gives 336
exposed bodies in total. Six rows have zero carrier failures. At
`(256,644)`, there are 235 exposed bodies and exactly two failures:

```text
14326401,1c306401,           FNV=936a8b25300381b7.       (6)
```

Complete activity evaluation over all `C(30,8)=5,852,925` rank-eight masks at
that pair gives

```text
active masks             1,465,388
active-mask FNV          7ee496dd5e3b67c6
response classes                  4
full responders                 367
least full responder       014c9084
exact replacement minimum          1.                   (7)
```

The lower bound in `(7)` is the nonempty failure family `(6)`; the displayed
mask attains it. It is new relative to `C`. Appending it gives

```text
|C+|=8,997,                  FNV=8e1860a25d0fcf87.       (8)
```

The repair is active at all seven endpoint-644 rows and is disjoint from both
bodies in `(6)`, so `C+` closes the entire top layer.

## 3. Descending carrier boundary and attribution

The exact descending scan completes a whole endpoint layer before deciding
whether to stop. It gives:

| endpoint | rows | pair FNV | failures under `C` | failures under `C+` |
|---:|---:|---:|---:|---:|
| 644 | 7 | `195e5d7d703b7d4c` | 2 | 0 |
| 643 | 14 | `a14eb5b1ee96edb4` | 0 | 0 |
| 642 | 9 | `32318bdcca33a0f4` | 0 | 0 |
| 641 | 7 | `399bd7d3c8d4b81d` | 0 | 0 |
| 640 | 14 | `b548beb345c96f2` | 0 | 0 |
| 639 | 4 | `fc85aa6b3ab0d13f` | 0 | 0 |
| 638 | 9 | `d03393624ca9d75f` | 40 | 40 |

This table fixes the attribution. The old carrier `C` already closes all 48
rows at endpoints `639..643`; `0x014c9084` is redundant there. Its proved
contribution is exactly the two bodies at endpoint 644.

At endpoint 638, the eight rows other than `(256,638)` are closed. The latter
has 16,792 exposed bodies and exactly 40 failures, FNV
`917d107c4536efc9`. The new mask is active but intersects every one of those
40 bodies, hence changes none of them. Complete quotienting of the active
rank-eight universe by its response to the 40 failures gives 315 response
classes, no full responder, and exact set-cover number nine. The breadth-first
union search over all inclusion-maximal response patterns supplies both the
lower bound and a nine-mask witness. No monotonicity in endpoint was assumed.

The exact sorted witness is

```text
02203226,081e1084,08a89440,180a8281,18261042,
18a0d040,1a82a200,202a9440,280a0a88.                  (8a)
```

It has FNV `02b936529030e4bc`. Appending `(8a)` to `C+` gives the
9,006-mask carrier `C++`, FNV `fdc1c57ae4dc1bb6`; exhaustive replay has zero
failures at `(256,638)`. The old masks already close the other eight rows at
endpoint 638, and adding masks cannot destroy coverage. Thus `C++` closes the
whole endpoint-638 layer. This makes no claim about endpoint 637 or below.

Therefore

```text
K={p in U_4282:638<=second(p)<=644},
|K|=64, FNV=c230e22462f7f3ab,
SHA256=36baf650ff0d7f8db8c8fb264693f5668302e61a8bd34dbb29f21767dcd31f00.
                                                               (9)
```

## 4. The seven-signature Boolean lattice

For the ordered joint deck `E`, write

```text
I_E(p)={i in {0,...,420}:E_i is inactive at p}.         (10)
```

The seven top signatures are

```text
S_220={25},
S_256={9,29,32,75,137,139,150,159,174,205,218,
       309,333,347,358,374,394,399,405,416,417},
S_258={412},
S_294={173},
S_366={396},
S_416={236},
S_512={107,374}.                                      (11)
```

Their union `J_*` has 27 indices; index 374 is the only overlap in the
displayed list. Retaining `E\J_*`, in the original order, exposes exactly 401
bodies, FNV `a149cb077a90ef39`.

For each nonempty subset `T` of the seven top rows, put

```text
J_T=union_{p in T} I_E(p).                             (12)
```

By THM-4282's exact surgery lemma, the bodies exposed after deleting `J_T`
are precisely

```text
X_E(J_T)={B:R_E(B) subset J_T}.                        (13)
```

The complete 127-scenario atlas applies `(13)`, intersects mask activity over
every row of `T`, and appends an explicit deterministic greedy cover `A_T`.
Across all 127 cases:

```text
deleted indices         1..27
exposed obligations     4..401
greedy witnesses        1..27
rebuilt deck sizes    418..424.                        (14)
```

Every rebuilt deck covers all `C(30,9)=14,307,150` labelled bodies, and all
members of `A_T` are active at every selected top row. Thus all selected top
rows are common for that deck. The construction makes no minimum claim for
the greedy cardinalities in `(14)`.

For a fixed scenario, retained masks are all active at a residual row `p` if
and only if `I_E(p) subset J_T`. Therefore the exact common family of the
rebuilt deck is

```text
F_T={p in U_4282:I_E(p) subset J_T and
                   every a in A_T is active at p}.      (15)
```

Exact evaluation of `(15)` over all 127 scenarios yields

```text
|union_{nonempty T} F_T|=636,
FNV=1c8b9758d5b215da,
SHA256=189309567e3259c4bec755a3e5b210b87680d864f3a4de4136e2245f164375b3.
                                                               (16)
```

The old-signature candidate fibre `(2)` has 640 rows. The four rows not in
`(16)` are exactly

```text
(206,263), (250,256), (256,394), (256,400),
FNV=a28328c12584aa2.                                   (17)
```

This is a failure of the first greedy witness choices, not of the signature
deletion mechanism.

## 5. Target-aware lifts and full-fibre closure

For a target `t` in `(17)`, intersect the response universe with activity at
`t` as well as at the selected top rows. Use the inclusion-minimal scenario
whose deletion contains `I_E(t)`. The exact lifts are:

| target | top scenario | deleted | obligations | allowed active responders | witnesses | deck / FNV |
|---|---:|---:|---:|---:|---:|---:|
| `(206,263)` | 2 (`S_256`) | 21 | 235 | 1,093,321 | 20 | `420 / 7af5e8fba00c61e9` |
| `(250,256)` | 10 (`S_256 union S_294`) | 22 | 275 | 1,095,771 | 21 | `420 / 736207470a8cefd7` |
| `(256,394)` | 2 (`S_256`) | 21 | 235 | 1,063,617 | 20 | `420 / 410dfffaa2e31dbf` |
| `(256,400)` | 2 (`S_256`) | 21 | 235 | 1,063,997 | 21 | `421 / 78280ad1f7d48696` |

Each appended family is active at its target and selected top row(s), and is
an explicit cover of every obligation. These are upper bounds only. Adding
the four target-aware families to the atlas changes `(16)` to the exact
identity

```text
union common families = G,
|G|=640,
FNV=45e9ecdf240e6417,
SHA256=3246d76e82e9e19d07e5851810da3107ad8fe98a1dfbd087edb2b9c5d8b27fa1.
                                                               (18)
```

A deterministic maximum-gain scenario cover followed by reverse deletion
extracts twelve of the explicit decks whose common families still have union
`G`. This is an irredundant selected cover, not a minimum-twelve theorem.

## 6. Independent validity gates

The primary exact common-family audit rederives the 401 global obligations,
checks the body response of all 127 greedy decks and all variants, and tests
236,156 needed witness-activity cells. There are zero equalities; the least
positive evaluated margin is `6619628848896`, while the closest negative
margin is `-8071850407824`.

The detached twelve-deck program imports no project source. It rebuilds
literal joint walls from the thirty speeds and each claimed pair, then:

```text
body universes                        12*C(30,9)
body short-circuit checks             4,892,555,059
claimed literal activity cells        1,537,426
activity equalities                             0
least positive literal margin          165,556,944.     (19)
```

The O0, O3/NDEBUG, and UBSan builds produce byte-identical stdout and literal
ledgers.

Independently, the detached carrier program imports no primary cocycle code.
It rebuilds literal walls for all 64 residual rows with endpoint at least 638,
classifies the 8,997 nested-carrier masks, and exhaustively scans every body
using the active carrier. It performs 29,962,550,588 short-circuit checks,
encounters zero equality, and reproduces exactly the primary 42 base and 40
nested failure ledgers. It also independently verifies all 61 incidences of
the nine-mask witness, least positive margin `291404330910`, and zero repaired
failures at `(256,638)`. O0, O3/NDEBUG, and UBSan outputs agree byte-for-byte.

The independently produced sibling packet
`05-knowledge/results/lrc14_endpoint638_carrier_response_thm4283/` supplies a
third implementation of the carrier prefix. Its primary quotient, exact-margin,
detached all-body, and Python consequence programs agree on the one-mask repair,
the 63 rows closed before appending `(8a)`, and the 40-obligation response
quotient with exact minimum nine. That packet remains a frozen precursor audit;
the present proof graph is the first one to append `(8a)` and use the resulting
64-row carrier band.

## 7. Typed proof graph and final residual

The carrier family `(9)` and common family `(18)` overlap in

```text
|G intersect K|=13,
FNV=c9171a79d21e375d,
SHA256=4660a98fc1c7ca28b2f895cb06cd7269728237a443249662284fd70cdaa67f67.
                                                               (20)
```

Equations `(9)`, `(18)`, and `(20)` give `(3)`. Its union ledger is

```text
count=691,
FNV=4b299b49d107a139,
SHA256=c5646e81b3815bdef5168e36bcd76174065ed21339a5d8853d9efddc8fa3efae.
                                                               (21)
```

Exact Python set replay and independently written Ruby replay agree on
`(20)`, `(21)`, and the residual `(4)`.

## 8. Mechanism, boundary, and scope

The structural mechanism is a response sheaf over the Boolean lattice of
inactive-signature unions:

- the signature records exactly which old masks must be deleted before a row
  can become common;
- the exposed-body response hypergraph records the exact body-cover debt;
- activity at the selected top rows is the fibre constraint;
- activity at a desired downstream row is an additional coordinate imposed
  during witness selection.

The four rows in `(17)` demonstrate why an arbitrary greedy cover is not a
complete downstream invariant. Their repair shows the strongest survivor:
every row in the global top-signature fibre has some explicit rebuilt common
deck. The carrier branch changes scale at `(256,638)`: one response no longer
suffices, but the exact nine-response witness closes the entire layer. The
next unscanned carrier boundary is endpoint 637.

This theorem proves neither that twelve common decks are minimum, nor that the
nine-mask boundary repair descends further, nor that carrier failure is LRC
danger. It proves no physical entry, no continuous reduction, and no
LRC(14). The latter remains open.
