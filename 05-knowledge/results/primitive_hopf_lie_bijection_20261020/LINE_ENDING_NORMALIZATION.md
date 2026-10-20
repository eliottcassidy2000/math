# CSV line-ending normalization

Packaging repair, 20 October 2026. The repository requires LF for text files. The initial approved CSV contained 62 CRLF line endings, so Git reported a clean-filter mismatch even though the fetched worktree and committed bytes were identical. This checkpoint removes only those 62 carriage returns and updates SHA256SUMS. Every CSV field and record is unchanged.

- Original approved research commit: `3072387c2551e6d235aeec693fc47966163f40af`
- Original CRLF bytes: 58009; SHA-256 `beaaa259a3f14597b9cbc0d62032d1ea08976c50ae5c735a5fc1273c66f4d10c`
- Current LF bytes: 57947; SHA-256 `100c00966207b06e44d7f8384ea08b7e8d143c47a280576f031c2eee6fb8cec4`
- Reversibility check: replacing each LF in the normalized file with CRLF reproduces the original bytes exactly.

The original publication manifest and original committed file remain historical evidence. Current checksum validation uses the normalized CSV and updated checksum list. No mathematical claim, computed value, or verifier was changed.
