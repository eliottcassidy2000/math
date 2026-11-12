# Release notes

Research publication only; Chan–Pak remains open. Independent review, exact proof controls and explicit scope boundaries accompany the results.

Verified earlier checkpoints:
- `34b3747a33c7099a1ccd22066eb733378394befe`: session plan.
- `620dc6bde65cf48bbc27a9e08c76cd15db458cbd`: initial complete central-fiber and semigroup result, before the additional reviewed sectors and reset decoder.

The final publication extends those checkpoints. It is based on the then-current remote main tree with a nonforced expected-head lease, retaining concurrent changes. Only this session's result directory is included; pre-existing unrelated working-tree changes are not staged or published.

`git diff --check` for the session packet passed. The repository-wide `agents/check_docs.py` check still reports two pre-existing unrelated startup-budget failures: META-PATTERNS has401 maintained lines against400, and startup surface has180862 bytes against180000. No global documentation clean-pass is claimed, and those unrelated files were not changed for this research release.

Reproduction reports and source hashes describe the mathematical artifact, not a global completeness certificate.
