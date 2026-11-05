# Singleton reassociation: exact defect and physical equality fiber

`PROOF.md` gives the exact N/M scalar defects when transferring one singleton
birth across a source-pair seam, the complete physical N equality fiber at the
smallest source, and minimal physical counterexamples. It does not exclude
all target-dependent source replacement or solve exact odd-target coverage.

Reproduce with Python standard library only:

```
python3 verify.py > VALIDATION-normal.json
python3 -O verify.py > VALIDATION-optimized.json
cmp VALIDATION-normal.json VALIDATION-optimized.json
```

The source universe is the complete forward grammar through U=501; algebra is
checked on all basis pairs and on every physical ordered source pair through
U=41. Minimal source orders are independently enumerated, and hostile grafts
are counted from actual adjacency rows. The all-order fiber proof additionally
uses the frozen physical source bounds and exact last-step theorem.
