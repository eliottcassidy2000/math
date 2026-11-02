# Release correction: the full shortcut hostile fibre

2026-11-02. Status: CORRECTED; no theorem or sharp-clock conclusion changed.

The original manuscript accidentally displayed `T^(-3)(7)={56,9}`.
The correct fibre is `T^(-3)(7)={56,18}`:

```
T^(-1)(7) = {14}
T^(-2)(7) = {28,9}
T^(-3)(7) = {56,18}.
```

The wrong display retained the preceding level's multiple-of-three node `9`
instead of its doubling predecessor `18`. Both are excluded by the unit filter,
so the original independent inverse enumeration and test that the surviving
3-unit fibre equals `{56}` were already correct. Nevertheless that filtered
check did not validate the full fibre printed in the manuscript. The verifier
now additionally checks and emits both full hostile fibres explicitly, making
this distinction regression-visible. No other proof statement was modified.

- Prior `PROOF.md` SHA256: `8384004d7e04695b4534dc4eac99988187e6403a8946dedb29c80ded43ee8b65`.
- Corrected `PROOF.md` SHA256: `e742e475d53cd2bd8532dc3aaef7ace817ede5fed9803b46477260f9f8bfc07c`.
- Prior mistaken line: section 3, Shortcut companion, displayed fibre.
- Correction identified by the root/release review.

Normal and optimized verifier outputs are regenerated and compared byte-for-byte;
`SHA256SUMS` records the release payload. `PUBLIC_FILES.json` selects portable
relative paths only; paths are resolved against this packet's directory.
