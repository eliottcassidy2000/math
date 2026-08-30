# THM-4283 endpoint carrier and signature-fibre packet

Status: **PROVED RELATIVE TO THM-4282 + FINITE-EXACT + DETACHED
LITERAL-WALL AUDITS PASS.** This packet removes an exact 691-row typed union
from the post-THM-4282 fixed-pool residual. It proves neither physical entry
nor LRC(14).

## Consequence

The proof graph has two distinct certificate types:

```text
full top-signature fibre G                  640 rows
repaired-carrier band K                     64 rows
G intersect K                               13 rows
typed union                                691 rows
```

The union has FNV `4b299b49d107a139` and SHA-256
`c5646e81b3815bdef5168e36bcd76174065ed21339a5d8853d9efddc8fa3efae`.
Its exact complement in the 23,373-row post-THM-4282 residual has 22,682
rows, FNV `f7563445f15efebf`, SHA-256
`7d0044bc4c32f08b9d09dca420171a05666d26e03f38fbc48a9baa1fcb027102`,
maximum endpoint 637, and top rows

```text
(100,637),(294,637),(520,637).
```

## Carrier branch

`results/carrier/top_band_scan.out` reconstructs THM-4282's ordered
8,996-mask carrier (FNV `fd899660f14b311c`). At `(256,644)` it leaves exactly
the two bodies `14326401,1c306401`, FNV `936a8b25300381b7`. Complete exact
enumeration of all `C(30,8)=5,852,925` masks gives four active response
classes, 367 full responders, and exact minimum one; the least full responder
is `014c9084`. Appending it gives an 8,997-mask carrier, FNV
`8e1860a25d0fcf87`.

The nested carrier closes every post-THM-4282 row with endpoint `639..644`,
exactly 55 rows. The inherited 8,996-mask carrier already closes endpoints
`639..643`; the appended mask is responsible only for the two endpoint-644
failures. At endpoint 638, only `(256,638)` fails, with 40 bodies. The appended
mask is active there but covers none of those bodies. The complete active
response quotient has 315 classes, no full responder, and exact replacement
minimum nine, so the one-mask augmentation stops there.

`results/carrier/endpoint638_response_witness9.txt` freezes the exact witness:

```text
02203226,081e1084,08a89440,180a8281,18261042,
18a0d040,1a82a200,202a9440,280a0a88
```

Its FNV is `02b936529030e4bc`. Appending it to the nested carrier gives a
9,006-mask carrier, FNV `fdc1c57ae4dc1bb6`, with zero failures at
`(256,638)`. The older carrier already closes the other eight endpoint-638
rows, so this proves the exact 64-row band `638..644`. No endpoint-637 or
lower-layer claim is made.

`results/carrier_detached/detached_literal_audit.out` independently rebuilds
literal joint walls for all 64 rows of the `638+` band and exhaustively scans
every labelled nine-body using the active carrier. Its 29,962,550,588
short-circuit checks reproduce the primary 42 base and 40 nested failures,
with zero activity equalities. O0, O3, and UBSan outputs agree byte-for-byte.

## Signature-fibre branch

The seven endpoint-644 inactive signatures are

```text
(220,644): {25}
(256,644): {9,29,32,75,137,139,150,159,174,205,218,
            309,333,347,358,374,394,399,405,416,417}
(258,644): {412}
(294,644): {173}
(366,644): {396}
(416,644): {236}
(512,644): {107,374}.
```

Their union has 27 deck indices. Retaining the other 394 masks exposes
exactly 401 labelled bodies, FNV `a149cb077a90ef39`. For every one of the 127
nonempty subsets of the seven top rows, `results/all127` freezes a deck made
by deleting the corresponding signature union and appending a deterministic
active greedy cover. The exact obligation counts range from 4 to 401,
witness counts from 1 to 27, and resulting deck sizes from 418 to 424. These
are explicit upper bounds; no greedy minimum claim is made.

The 127 decks are common on a 636-row union. Exactly four rows in the full
global-signature candidate fibre remain:

```text
(206,263),(250,256),(256,394),(256,400).
```

`results/target_lifts` intersects the response universe with activity at each
target and its inclusion-minimal top scenario. Four alternative explicit
decks close those rows. Consequently the exact fibre

```text
G={p in U_4282 : I_E(p) is a subset of the 27-index top union}
```

has 640 rows, FNV `45e9ecdf240e6417`, SHA-256
`3246d76e82e9e19d07e5851810da3107ad8fe98a1dfbd087edb2b9c5d8b27fa1`,
and is the common-row union of twelve explicit rebuilt decks. The exact
primary audit evaluates 236,156 needed witness cells with zero equalities.
The detached audit then scans `12*C(30,9)` bodies (4,892,555,059
short-circuit checks) and 1,537,426 claimed literal activity cells, again with
zero equalities. O0, O3, NDEBUG, and UBSan controls agree byte-for-byte.

## Proof-graph replay

- `results/proof_graph/run.out` is the exact Python typed-set replay.
- `results/proof_graph/ruby.out` is an independently written Ruby replay.
- `results/proof_graph/carrier_band64.csv`,
  `results/common_families/full_signature_fibre_common_union.csv`, and
  `results/proof_graph/carrier_common_overlap.csv` keep the certificate types
  separate.
- `results/proof_graph/proof_union.csv` and `final_residual.csv` are canonical
  convenience ledgers, not a claim that one deck realizes the union.

## Canonical sources

```text
04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/
  carrier_scan_support.cpp
  endpoint_top_band_scan.cpp
  endpoint638_exact_response_witness.cpp
  all127_response_greedy.cpp
  targeted_response_lift.cpp
  all127_common_family_audit.cpp
  selected12_detached_literal_audit.cpp
  carrier_band_detached_literal_audit.cpp
  proof_graph_consequence.py
  proof_graph_consequence_independent.rb
```

`carrier_scan_support.cpp` is a mechanically inherited, source-pinned copy of
THM-4282's primary carrier scanner with its old main renamed. The two detached
auditors import no project source.

## Scope and exclusions

- The carrier, each rebuilt common deck, and the proof-graph union are
  different typed objects. They are never substituted for one another.
- The only new minimum claims are the exact one-mask response at `(256,644)`
  and the exact nine-mask replacement boundary for the 40 failures at
  `(256,638)`. Greedy and target-aware decks are explicit upper bounds.
- The 9,006-mask carrier is not claimed to close any endpoint-637 row. Failure
  of an earlier carrier is not an LRC danger witness.
- Binaries, timing streams, and exploratory discarded covers are excluded.
- The root `SHA256SUMS` is the authoritative raw-LF byte ledger.
