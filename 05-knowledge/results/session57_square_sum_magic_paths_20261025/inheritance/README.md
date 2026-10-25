# Session57 inherited connection board

Read CONNECTION_BOARD.md for the eight ranked source/target maps and the exact limits of each transfer. This is a research board, not a solution of the integer magic-square-of-squares problem.

The added elementary theorem is the short-gap/AP-petal bijection; the added universal negative control is that magic opposite pairs are square-sum nonedges. The proof is in the board. The Python certificate independently compares direct opposite-square enumeration with pair-root AP enumeration for all center roots1..500, checks68 positive petals, and checks explicit gap stocks at5,41,65,1105. The finite test supports but does not replace the proof.

Run python3 verify.py and python3 -O verify.py from any working directory. Output should equal certificate.json. All substantive checks use explicit exceptions. The verifier also validates the14 source snapshots against SOURCE_PINS.json; no network or original checkout is required.

Publish only the files named in FILES.json. Other source-recovery scratch files in the working directory are not part of this packet. SOURCE_PINS.json provides immutable commit URLs and SHA-256 hashes. External citation: https://arxiv.org/abs/1908.03236v2 (Onno Cain, 2019, version2).
