# Quaternary unit-clone obstruction

The reviewed all-order conclusion is A>=60 for any eventual fixed-root transitive-clone identity H_m=A*4^m+1 or A*4^m-1. There is no symmetry assumption and no claim of existence at the bound. The proof also establishes the general factorial and neutral-SCC reduction for unit clocks of any base.

Read `PROOF.md` and its independently hash-bound `REVIEW.md`. Run `python3 verify.py` and `python3 -O verify.py`; they compile the two included C++17 counters and reproduce the receipts. The complete strong H<=217 bank is included and hash-pinned to the previously reviewed equality classification.

A third independent method is in `audit/root_revisit.cpp`. Compile it with a C++17 compiler and supply CSV lines `n,H,row0,...` on standard input, obtained from `strong-bank.jsonl`. Its deterministic output is `audit/receipt.json`. It uses revisited root placeholders rather than either source counter.

The source proofs/review retain their original named scratch-era paths for hash stability. In this release, all four source/code files and the bank/receipts are adjacent, and the independent counter is under `audit/`. The published equality dependency is linked directly in the proof. `RELEASE_FILES.sha256` pins every released file except itself.
