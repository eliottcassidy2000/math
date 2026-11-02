# Verification and limitations

```sh
python3 -B /absolute/path/to/this/packet/replay_checkpoint1.py
cd /absolute/path/to/this/packet
sha256sum -c SHA256SUMS
```

The replay relocates only the selected frozen scientific files into a fresh
temporary tree and executes them from a separate unrelated working directory.
Ordinary and optimized(`-O`) Python must reproduce the frozen output bytes.
All gates use explicit exceptions, not assertions disabled by optimization.
The exact configuration and outcomes are in REPLAY_CONFIG_CHECKPOINT1.json
and REPLAY_CHECKPOINT1.json.

Scientific hostiles change a small strong spectrum, a polynomial determinant
claim and an independent exact identity. A same-size proof-byte corruption
must fail the source-integrity check in both modes. These prove the named
gates fail closed, not that every conceivable implementation bug is caught.
Both Git blob and SHA256 identities of selected source snapshots are verified
offline; the packet seal covers every published file except its own checksum
file. Independent review is a separate human-readable proof/status audit.

Infinite statements depend on their written proofs and cited hypotheses;
finite calculations are corroboration. The exhaustive n<=6 Boolean census
has33867 labelled inputs. The structural exclusion cover has272 overlap
presentations, not272 representatives of the entire labelled universe.
Arithmetic target-monoid tests do not establish tournament realization.

Documentation checks are deliberately scoped to startup budgets, unchanged
startup truth surfaces, owned-file conflict markers and curated publication
links. The full repository checker was **not run**, because it reads the
central mistakes ledger. This is not a full-repository documentation pass.
No Lean/Isabelle build, global spectrum proof, unrestricted construction
classification or independent rerun of the large inherited prefix is claimed.
