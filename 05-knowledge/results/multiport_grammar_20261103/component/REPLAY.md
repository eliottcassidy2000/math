# Frozen component replays

2026-11-03. All three producer sources and receipts were independently replayed;
normal and optimized modes agree byte for byte. Exact source and receipt pins,
universe definitions and independent permutation controls are in [AUDIT.md](AUDIT.md).

| Receipt | Scope | SHA256 |
|---|---|---|
| `component_inverse.json` | Four-bank levels 0–4, 108 inverse comparisons, 63 literal signatures, 14 attachments | `a4415b55139346e6facc7e611f5a0070966357510687fc34c3c433acceab4145` |
| `component_inverse_five.json` | Five-bank levels 0–4, 108 inverse comparisons, 108 literal signatures, 20 attachments | `225e228394ac04faea87ec663fb1391ecb1256a27eae61b44cbb4d5db99978f6` |
| `spectrum_acceptance.json` | Every odd target through 2000, 1,999 inverse comparisons, 72,190 candidates, 997 literal attachments | `cb2c54346f46f719ee2d3fa946873a1c59c18c0d9f3c3fcdcc97f5c84c77027e` |
| `audit_component.json` | Independent permutation/cut, marked-isomorphism and first-step closure checks; 18,770 gates | `b1719d7a7e787ef2947dfe8cddc4efb34fb2eb1eb5add80b1a78e2f9d3581d93` |

A preliminary left-inverse-only scalar scan through 10,000 also missed only
7 and 21 (878,679 candidates). It was a direction-finding experiment and is
not the independently audited receipt. The retained executable can reproduce
and fully audit that larger range with `spectrum_acceptance.py --bound 10000`
and a separate output path. No infinite conclusion is inferred from either
finite range.
