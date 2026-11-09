# Session231 strongest-bank arithmetic lane

No Chan–Pak proof or F2-cofiniteness proof was found. `PROOF.md` recovers the exact applicable Session220 obstruction and gives a verified, magnitude-preserving arithmetic inverse for the exact S83 family. Its low-bit-cylinder strategy is inherited from Session217, specialized here without a source-length loop. The mathematical contribution is a clean one-candidate certificate; exact branch occupancy remains open.

Files: `PROOF.md`, `inverse.cpp`, `verify.py`, `bank.txt`, `CASES.txt`, `CHECK.json`, `CHECK.optimized.json`. The binary `inverse` is generated and need not be published. `bank.txt` is copied from the parent's frozen repository HEAD's rooted-triple bank, not newly generated.

Independent audit: ../source_audit/REVIEW.md and ../source_audit/ROUND2_REVIEW.md. The published scripts and receipts include only the verified source scope.

Round2 is in `ROUND2/`: an arbitrary-context two-clock bank-rigidity proof, independently audited, exact half-word coordinates/decoder for the inherited centered-fan class, literal graph checks and scoped discovery experiments. Read its README for provenance and reproduction. No completeness claim.
