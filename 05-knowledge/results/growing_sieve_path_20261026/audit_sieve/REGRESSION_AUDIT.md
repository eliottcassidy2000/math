# Bounded audit of the incoming S60 regression repairs

**PASS, 2026-10-26; bounded review, not a full THM-4805 audit.**
The record `../REGRESSION_CORRECTIONS.json`, SHA-256
`5fc1c9092928605c900766aed2d7af497aeb20b5642025d7e19bcc7f67e85759`,
contains exactly 18 replacements in six maintained files, based on incoming
revision `2c85cab14a`. The auditor independently reconstructed every file from
`git show source_revision:path`, required every old anchor to occur exactly
once at its step, and checked the recorded before/after hashes and final bytes.
All six reconstructions pass.

The reviewed mathematical changes are restricted to the established repairs:

- The Mersenne singleton conclusion keeps S57's finite range through k=26,
  distinguishes the later audit's **reported finite** partitions through 80,
  and retains distinct-prime Goldbach for the all-exponent conclusion.
  Calling the result CITED does not discharge that hypothesis.
- The duplicated older A'(ii) statement is removed while the repaired scoped
  statement remains. The obsolete duplicate results-index entry is removed
  while the newer entry remains and receives the same conditional scope.
- The first-handful source-shape prose separates the 0 and 1 fibres and
  restricts its square-only claim to targets 2,...,5. The broader three-shape
  regime is restricted to 2,...,14.
- The positive abc domains n>=2 and k>=1 are restored consistently.
- The mistake entry records merge regression and provenance without deleting
  the earlier correction or MISTAKE-675.

One intermediate replacement described a Mersenne value as having a “unique
row-of-1 source”. Literally all Mersenne values have a unique source *within
that row*. Before freeze, root replaced this by the precise statement that its
**entire aliquot fibre** consists of that row-of-1 source. The final replayed
record contains the repaired wording. No unrelated numerical result or S60
claim was re-audited; the new reported evidence through 80 remains cited at
its inherited finite scope.
