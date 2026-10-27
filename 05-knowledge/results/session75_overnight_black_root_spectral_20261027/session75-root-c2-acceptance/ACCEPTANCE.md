# C2 final-byte acceptance — native ROOT packet

**ACCEPTED WITH PRESERVED QUALIFICATIONS.** Final producer report reread; final checker and receipts inspected. Fresh normal and `-O` replays reproduce the producer receipts byte for byte. All 12 payload checksum entries and all seven source-ledger hashes match. No mismatch was found; producer bytes and the earlier independent audit remained unchanged.

## Exact reviewed SHA-256 inventory

- `REPORT.md`: `553963c69c948b315c055473c826e14a77f31cb51ae779cfb24458d7b55ab5aa`
- `check.py`: `e45112beb8ed4d035d7c21bcbdb42301f39519be0f648f0de571cb7cf360eb34`
- `normal.json`: `51c0e70941c3163492d6d2b238a7e87dd6a9d05e94f1f87a9d0a1367d883e8bf`
- `optimized.json`: `51c0e70941c3163492d6d2b238a7e87dd6a9d05e94f1f87a9d0a1367d883e8bf`
- `SOURCE_PINS.json`: `d8d2c7712d3862d9789a71494faf01ef6b29513c34110dd70ef191ba7b0c778d`
- `SHA256SUMS`: `2eeb0ffd5c0346ce65efad647f8a5079b0de107376731f0521caff054f5603a3`

The JSON acceptance binds all 13 producer files, including the seven source files and SHA256SUMS itself.

## Qualifications retained

1. Positive-b f=0 in e_next=b*e+f is not signed-c F=0 in e_next=epsilon*b*e+F. The former yields terminal B_s*A, the latter epsilon_s*B_s*A. The producer invoice is accepted in its explicitly defined signed-c convention; this acceptance does not identify the forcing variables.
2. A single global homogeneous profile with fixed ROOT value z is distinct from separate source-labelled paths with terminal epsilon_s*A_N(s). Initial inflation A_N(s)/B_s and its energy invoice are accepted sourcewise, not as a construction of one common ROOT boundary value for all prescribed sources.

This is provenance closure, not new research or a broader theorem promotion. Existing frozen audit: `/tmp/session75-root-audit`; neither that packet nor the producer packet was edited.
