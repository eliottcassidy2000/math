# Index-367 singleton-fibre surgery promotion packet

Status: **PROVED RELATIVE TO THM-4281 + FINITE-EXACT + DETACHED LITERAL
AUDIT PASS.** This packet does not prove a physical entry or LRC(14).

## The result

On the exact 24,223-row post-THM-4281 residual, evaluate every repair in the
frozen 421-mask joint deck rather than stopping at the first inactive repair.
The resulting 421-bit inactive-signature census has:

- 24,223 rows and 864,784 inactive incidences, with zero equalities;
- 17,604 distinct signatures and Hamming-weight range 1 through 421;
- 5,038 weight-one rows distributed over 215 deck masks;
- 421 distinct mask columns, no nontrivial equal columns, and 11 strict
  pairwise inactivity implications;
- top signatures
  - `(256,663)`: weight 7;
  - `(366,663)`: exactly `{index 367 = 0x02188125}`;
  - `(520,663)`: weight 5.

No inactive mask is shared by all three top rows.  Only indices 107 and the
pair `(256,663),(520,663)` meet among the three pairwise intersections.

The weight-one rows force 215 masks into any subdeck that certifies every
post-THM-4281 row as noncommon.  They cover 24,222 rows.  The sole uncovered
row `(70,302)` has signature `{78,368}`, so there are exactly two minimum
hostile classifier bases, each of size 216: the forced 215 plus index 78 or
index 368.  This is an order-free exact lower bound and cover.

The complete signature fibre of `(366,663)` has 26 rows, all with inactive
signature exactly `{367}`.  Deleting mask `0x02188125` exposes exactly three
private labelled nine-bodies:

```text
05646408,25065480,35a24080
```

Their ordered FNV is `22bd212c1ffec6e2`.  At `(366,663)`, complete enumeration
of all `C(30,8)=5,852,925` repairs finds 3,095,104 active repairs and 61 active
repairs disjoint from all three private bodies.  Therefore one replacement is
necessary and sufficient.  The least pair-local replacement is `0x02108325`.
It is the one-element exchange

```text
remove pool label 145 (bit 19), add pool label 63 (bit 9).
```

On the detached literal grid at `(366,663)`, with common denominator
`1112710724405280`, the deleted and replacement signed margins are exactly

```text
0x02188125: -8854821928638, normalized -19164539/2408245840;
0x02108325:  14468348998404, normalized 57414083327/4415518747640.
```

The pair-local least replacement is not active across the other 25 fibre rows.
Testing all 61 complete responders on all 26 rows leaves exactly two
fibre-wide replacements, `0x0a188803` and `0x0a188a01`.  The least,
`0x0a188803`, removes pool labels `{15,30,60}` and adds `{10,84,264}` relative
to the deleted mask.  Appending it after the retained 420 masks gives a
421-mask deck with FNV `87b42cf8a2069177`.  A complete labelled body scan has
zero failures, and detached literal evaluation proves all 421 masks strictly
active at every one of the 26 fibre rows.  The weakest replacement gap on the
fibre is at `(227,370)`:

```text
2213184690321570 / 153207497939015520 > 0.
```

Thus all 26 fibre rows are proof-graph-new.  Deleting them from the exact
post-THM-4281 residual leaves:

```text
count   24,197
FNV     c9746d48cbf3fc37
SHA256  8d59477a996503c338434dda77ae6ad897637d9d090101be118d179631c73e09
top     (256,663),(520,663)
```

## Proof architecture

1. `full421_inactive_signatures_primary.cpp` uses the inherited exact
   endpoint-cocycle prefix engine and tests all 421 masks at all 24,223 rows.
2. `analyze_full421_signatures.py` performs the order-free fibre, weight,
   implication, nearest-neighbour, and exact hostile-basis analyses.
3. `index367_private_bodies.cpp` independently exhausts all 14,307,150
   labelled nine-bodies and isolates the exact private fibre of mask 367.
4. `endpoint663_index367_response_atlas.cpp` exhausts all rank-eight repairs,
   forms the exact three-bit active response quotient at `(366,663)`, proves
   the one-replacement minimum, and rescans the resulting pair-local deck on
   all labelled bodies.
5. `index367_singleton_fibre_literal_audit.cpp` uses a detached literal-grid
   engine, rechecks the complete original signature at all 26 rows, tests all
   61 full-response replacements, and rescans the fibre-wide rebuilt deck on
   all labelled bodies.  O3 and UBSan outputs agree byte-for-byte; sanitizer
   stderr is empty.
6. `proof_graph_consequence.py` derives the 26-row fibre directly from the
   signature CSV, checks the exact one-mask deck surgery, and emits the new
   residual ledger.

## Scope and exclusions

- The fibre-wide deck is used only on the 26 stated rows.  No claim is made
  that it preserves the original THM-4281 148,063-row global common set.
- The pair-local one-swap deck and fibre-wide deck are distinct; the one-swap
  mechanism must not be attributed to `0x0a188803`.
- The 216-mask hostile bases classify failure of the frozen deck on the
  post-THM-4281 residual.  They are not body covers and are not claimed to be
  common-activity proof decks.
- No endpoint scan below 663, physical-entry claim, or LRC(14) consequence is
  included.
- Binaries and deterministic progress stderr are excluded.

`SHA256SUMS` is the authoritative raw-LF artifact ledger.
