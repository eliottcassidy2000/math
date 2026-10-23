# Final Singer publication-scope verification

Verdict: PASS after the characteristic-seven wording correction below. Verified 2026-10-23 05:03 UTC.

Reviewed the complete final root note, companion program, and exact JSON receipt. The note distinguishes a Lie-automorphism action from a subgroup of the unipotent Heisenberg group; states the standard SL8 branching input; restricts the E7 conclusion to the distinguished augmentation character; and does not confuse an unavailable university PDF attempt with the independently recovered author-uploaded full text.

One correction was requested and has been applied: characteristic7 must not be grouped with3 and19 as a failure of this permutation module's semisimplicity or splitting. Neither57 nor171 is divisible by7, and the form's polynomial gcd remains x-1 over F7. The final note now confines its theorem to characteristic zero and flags3 and19 without making a positive-characteristic E7 claim. No remaining correction is required.

## Final artifact SHA256 pins

- session35-root/SINGER_HEISENBERG_E7_OBSTRUCTION.md
  377096e575dd62395aaa8842231417eb9576d4fb1526fa084c2650781deb334a
- session35-root/check_singer_heisenberg.py
  14e2559615ab9a936e9dce71e5453b9fecf0501e7b9a21eb8c688d0fc0b0b1b5
- session35-root/singer_heisenberg.json
  1d9d32aa3b7742439fb4721c8c585406f327b5fcd1dbb5092ec48649358ad093

Normal and python3 -O replays both exited0 and matched the receipt byte-for-byte. Replay copies: root_singer_replay.json and root_singer_optimized.json in this audit directory, both with the same receipt SHA256 above. The exact Bezout identity, allowed SL8 triples, cyclic difference-set property, and Frobenius coefficient invariance were therefore rerun rather than inferred from saved output.

## Verified source pin

Authors: Robert L. Griess Jr. and A. J. E. Ryba.
Title: Classification of finite quasisimple groups which embed in exceptional algebraic groups.
Document date visible in the author upload: October5,2004.
Stable source URL:
https://www.researchgate.net/publication/243106241_Classification_of_finite_quasisimple_groups_which_embed_in_exceptional_algebraic_groups
Locator: printed page39, “Rows for G2(3) in Table QE,” order-three trace discussion. Full text was opened and the values {-7,-25,2,20} on the56-dimensional module were inspected. The page identifies its full text as uploaded by coauthor Alex Ryba. This is a URL/edition/passage pin, not a purported binary-PDF hash: the source was read through the web text extractor, not downloaded as a local PDF.

No Git operations or repository edits were made by this auditor.
