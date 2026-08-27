# THM-4282 final8996 detached literal carrier audit

Status: **FINITE-EXACT / HOSTILE-AUDIT PASS** for the complementary carrier
claim on the canonical post-surgery residual set.  This is not a claim that
the carrier is the THM-4281 common deck.

## Scope and identities

- Reconstructed base carrier: 8,951 ordered masks, FNV-1a
  `188f82ab9dd1695a`, byte SHA-256 `a5dac3c7...`.
- Frozen joint deck inside that carrier: 421 masks, FNV-1a
  `20d63dd42fe8150e`.
- Canonical post-(520+367)-surgery residual: 23,637 ordered pairs, FNV-1a
  `e8b363d2b3d9ba6a`.
- The audit derives, rather than accepts as a separate pair input, the exact
  90-row slice with endpoint `645 <= r <= 663`; its FNV-1a is
  `942995bee7469430`.
- Additions: 45 ordered rank-eight masks, FNV-1a `ec083b65cc8c34e3`.
  Appending them to the reconstructed base gives 8,996 masks with ordered
  FNV-1a `fd899660f14b311c`.

## Independent method

`source/detached_literal_continuation_audit.cpp` uses the direct joint-wall
geometry in `source/lrc14_endpoint_cascade_direct_wall_body_audit.cpp`.  It
does not include or call `response_pattern_atlas.cpp`, `PrimitivePair`, or
`build_active_universe`.  For each of the 90 pairs it reconstructs exact
literal activity, enumerates all `C(30,9)=14,307,150` bodies, partitions the
base carrier into active joint/non-joint masks, and checks exact body cover at
cuts 0, 1, 2, 11, 42, and 45 additions.  Final failures are zero.

The complete response controls enumerate all `C(30,8)=5,852,925` rank-eight
masks using an independent literal zeta construction.  Exact OR-cover search
certifies:

- `(256,660)`: five base failures; least full responder `084a0a81`;
- `(256,657)`: one failure after the first seed; least full responder
  `0016580c`;
- endpoint 650 small fibres: minima `1,4,3,1` at q values
  `294,366,416,512`, respectively, with the claimed nine masks;
- `(256,648)`: six failures after cut 42 and exact minimum three, closed by
  `00524a81 0128d00c 04410e81`.

At endpoint 650, cut 2 has exactly the five resistant rows/counts
`256:368, 294:2, 366:18, 416:14, 512:1`; cut 11 leaves only q=256 with 358
failures; cut 42 closes endpoint 650.  Immediately before the final three
masks, the only resistant row anywhere in the audited 90-row band is
`(256,648)`, with six failures.

## Consequence

`source/derive_final8996_consequence.py` validates the audited TSV against the
slice independently re-derived from the 23,637-row residual, removes exactly
those 90 pairs, and freezes the 23,547-row result.  Its FNV-1a is
`e07243b2f73f978f`, maximum endpoint is 644, and the top fibre has seven rows:
`(220,644) (256,644) (258,644) (294,644) (366,644) (416,644) (512,644)`.
Accounting is `24,223 - 586 - 90 = 23,547`.

## Controls and limitations

- `source/reconstruct_final8951.py` independently replays all carrier stages;
  its frozen output reports byte identity with the dormant base carrier.
- The final literal scanner was built with Apple clang 17, C++20, `-O3
  -DNDEBUG -pthread`; all checks are explicit `require` gates, not assertions.
- Two separately compiled versions produced byte-identical 90-row TSV output
  SHA-256 `d8dfdb69...`.
- The progress file is non-semantic and thread-order dependent; reproducibility
  compares semantic stdout and TSV only.
- A stale scratch 42-addition identity `35f344...` came from the earlier greedy
  ordering.  The current ordered cut-42 carrier identity is
  `05ff3d7ecbaa740c`; final cut-45 remains `fd899660f14b311c`.
- The parent primary scanner originally lacked internal pair/addition pins.
  `controls/primary_scanner_hostile_summary.txt` records that finding.  The
  detached program closes that gap by deriving and pinning the 90-pair slice
  and pinning all 45 additions.

Run `REPRO.md` from this packet directory.  `SHA256SUMS` is the byte ledger.
