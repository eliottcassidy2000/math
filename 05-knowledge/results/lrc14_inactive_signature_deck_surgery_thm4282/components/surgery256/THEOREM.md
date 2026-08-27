# Solver-free delete-seven/add-eight common-deck consequence

Status: **FINITE-EXACT / promotion candidate relative to the frozen
THM-4281 and sibling proof-graph ledgers.** This packet makes no physical-entry
or LRC(14) claim and reserves no theorem namespace.

## Explicit deck witness

Start with the ordered 421-mask THM-4281 joint deck, FNV
`20d63dd42fe8150e`. Delete the zero-based indices

```text
75,107,139,374,394,405,417
```

whose masks are

```text
21948006,128c8900,10259240,20027016,10a41016,2051d200,01188550.
```

Retain every other mask in original order, then append, in this order,

```text
0110a550
04871108
10241207
1042d088
12848902
21141284
30249140
31202206
```

The resulting 422 distinct rank-eight masks have FNV
`1c97b54ece61b351` and text SHA-256
`9484ab76c89fd9f18cfd5f9d5bf996a412f7b94e6d1770286380ef3b6238f4ed`.

An exhaustive scan of all `C(30,9)=14,307,150` labelled nine-bodies finds
zero failures after 405,775,694 short-circuit checks. The eight appended masks
cover the exact 71 obligations exposed by the seven deletions, with 82
incidences and between one and two hits per obligation.

This is an explicit cover witness only. **No minimum-cardinality claim is
made.** The exploratory HiGHS result used during discovery is noncanonical and
is deliberately excluded from this packet and its proof graph.

## Exact current-residual common family

On the frozen 24,223-row post-THM-4281 residual, the primary exact activity
route finds precisely 188 rows on which all 422 masks are active:

```text
count   188
FNV     6588121dbec57bcb
SHA256  5946ff653c51a74eba09a14430c494074e53b5aba87c3159bd17bafbe9e605d5
```

There are 395 rows whose old 421-bit inactive signature is contained in the
deleted seven-set. Testing the eight replacements on all 395 gives 3,160
exact cells, no equality, and the 188 accepted rows. A detached literal-wall
implementation independently verifies the explicit deck, body coverage, and
the exact 188-row activity family without importing the primitive-cocycle or
atom/zeta activity code.

## Proof-graph union

The sibling `(520,663)` plus index-367 certificates supply a 586-row deletion
ledger, FNV `9769742754535d84`. Its complementary residual has 23,637 rows,
FNV `e8b363d2b3d9ba6a`. Inside that residual, the inherited carrier descent
removes exactly the 90 rows with second endpoint at least 645, FNV
`942995bee7469430`.

The three ledgers have the following exact overlaps:

```text
surgery188 ∩ prior-gain586     9   FNV 1f049d3f1c65c5c9
surgery188 ∩ carrier90         5   FNV f3cb9446a58b5543
prior-gain586 ∩ carrier90      0
triple intersection            0
```

Thus the surgery contributes 174 new rows and the exact union contains

```text
188 + 586 + 90 - 9 - 5 = 850
```

rows. Its ledger has FNV `8f595510210a5785` and SHA-256
`7ad581bccd253e1778b972e8a303207da44534e6b995fa3ba15bd34b2801505b`.
Removing it from the post-THM-4281 residual leaves:

```text
count   23,373
FNV     c6ab0ae49ee32273
SHA256  c3e5bf37887aa57af79cb166fce4a6e933e5daffc26dd8032fdfc52ce31240f3
max r   644
top     (220,644),(256,644),(258,644),(294,644),
        (366,644),(416,644),(512,644)
```

The Python proof-graph replay and an independently written Ruby set audit
agree byte-for-byte on every output ledger.

## Serialization and scope

Mask FNV-1a-64 serializes each ordered mask as one little-endian `u64`. Edge
FNV serializes `q` then `r`, each as little-endian `u64`, in lexicographic
ledger order. Edge SHA-256 is over raw ASCII `q,r\n` bytes.

The conclusion is finite and conditional on the frozen input ledgers. It does
not claim that the eight-mask append is minimum, that the 422-mask deck
replaces THM-4281 globally, or that any endpoint below 644 is closed.
