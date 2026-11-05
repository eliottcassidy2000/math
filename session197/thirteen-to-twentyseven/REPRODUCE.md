# Reproduce the reviewed induced13-to27 obstruction

The frozen source and independent review are copied unchanged. There are no external runtime or theorem dependencies: the all-order reduction uses SCC products, induced path monotonicity, the proved strong-deletion lemma and an independently regenerated empty order-eight frontier.

Run from this directory with standard-library Python:

```sh
sha256sum -c RELEASE_FILES.sha256
python verify.py > /tmp/s197-thirteen.json
python -O verify.py > /tmp/s197-thirteen-optimized.json
cmp /tmp/s197-thirteen.json VALIDATION.json
cmp /tmp/s197-thirteen-optimized.json VALIDATION_OPTIMIZED.json
python audit/closed_frontier.py > /tmp/s197-thirteen-independent.json
cmp /tmp/s197-thirteen-independent.json audit/closed-frontier.json
```

The independent checker regenerates all1,572 extension candidates using literal Hamilton words and cycle canonicalization. Its only file read is audit/final-source/VALIDATION.json, an exact copy of the main receipt used for the final bank comparison. All graph generation and counts are recomputed independently before that comparison. The main checker uses no data file.

Original research paths in the proof and review are provenance references. PUBLICATION.sha256 covers every delivered file except itself. The known low-count gaps are a consequence of this closed bank, not a newly claimed spectrum theorem beyond27. The H15-to-H29 positive extension is retained as a hostile to overgeneralizing this result.
