# Reproduce the reviewed retained-unit-clock-reduction packet

The source, original six-file manifest, independent review and its checkers are copied without changing their frozen bytes. Original research paths in those documents are provenance references, not runtime requirements.

From this directory:

```sh
sha256sum -c RELEASE_FILES.sha256
python verify.py > /tmp/s197-source.json
python -O verify.py > /tmp/s197-source-optimized.json
cmp /tmp/s197-source.json VALIDATION.json
cmp /tmp/s197-source-optimized.json VALIDATION_OPTIMIZED.json
python audit/independent.py > /tmp/s197-independent.json
cmp /tmp/s197-independent.json audit/independent.json
```

The source uses only the included core.py, controls.json and the Python standard library. The independent checker is standalone. Exact named mathematical proof snapshots are in dependencies/ with their original hashes, as bound in the proof and review. These copies document the inherited theorems without making the source verifier depend on their code or data.

The control bank supports the all-order structural proof; census completeness is not a premise of either reduction. This is not a scalar coverage theorem. PUBLICATION.sha256 covers all delivered files except itself.
