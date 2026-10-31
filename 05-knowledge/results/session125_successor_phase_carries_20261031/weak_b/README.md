# Weak-B ninja pairs: radix-25 path-bank carries

**Separate addendum; not a replacement for the original-port Q49 theorem.**
Weak-B pairs start 1,3, end at 8 for even order and 10 for odd order, and share
common-label position parity. Read [the proof](WEAK_B_CARRY.md) and
[independent audit](audit/AUDIT.md).

- No fixed-phase weak-B family spans all orders **42–47**: all six onset
  phases have a later internal vertex with only one compatible neighbor.
- Permitted head cuts preserve the exact carry law
  `H(25r+12)=25H(r)+K`, with K≥0.
- Every selector in the stated constructor's two-representative **path bank**
  has limsup deletion ratio at least **4/93**, for any valid base choices.
- The owned A/B certificates have K=0 in both parities. Their seed-45 ray has
  local two-step maximum edit-rate limit **28/91**. No global maximum is claimed.

These are path-bank results, not restrictions on arbitrary graph families.
Odd weak-B paths cannot close directly to 1 because 1+10=11 is not square.
The separate cycle-consumer maps have not been given an edit-cost theorem.

## Owned dependencies and portable replay

Python 3.10+ standard library only; **no external C source**. Both checkers locate
all inputs relative to their own published source location, independently of
the working directory:

- `../../session124_mixed_radix_square_rays_20261031/` — owned Q25 A/B words;
- `../../square_sum_B_base_interval_20261019/` — manifest and compressed base chunks;
- `../../session124_square_prefix_threshold42_20261031/SMALL_BASES.json`;
- `../../session125_origin32_numeration_20261031/pair_taxonomy/weak_b_certificate.json`.

[DEPENDENCIES.json](DEPENDENCIES.json) lists pins and source-relative defaults.
Use `--operator-dir`, `--base-dir`, `--small-bases`, and `--weak-b-catalogue`
when those inputs are elsewhere. No data bundle is duplicated in this addendum.

From this directory, with the published siblings present:

```sh
python3 check_weakB_carry.py --output producer.normal.json
python3 -O check_weakB_carry.py --output producer.optimized.json
cmp producer.normal.json producer.optimized.json
python3 audit/independent_audit.py --output auditor.normal.json
python3 -O audit/independent_audit.py --output auditor.optimized.json
cmp auditor.normal.json auditor.optimized.json
sha256sum -c SHA256SUMS
```

The included receipts pass and agree in normal/optimized modes. They cover
1,020 bases, 100 paired states, 12 bounded ray samples, all six phase
obstructions, and bank-selector controls. The independent audit constructs the
phase census separately, decodes bases sequentially, and extracts boundary
edges from expanded paths rather than importing the producer's symbolic builder.
[RELEASE_CHECK.json](RELEASE_CHECK.json) records default-layout relocation,
negative controls, the strict whitelist, and preservation of the frozen parent
packet. No new mathematical search was performed during curation.
