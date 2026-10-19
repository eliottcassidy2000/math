# Sharp induced-book height: proof and reviews

The [full proof](BOOK_HEIGHT_ASYMPTOTIC.md), [arithmetic audit](review/ARITHMETIC_REVIEW.md), and [independent coordinating review](review/ROOT_PROOF_REVIEW.md) establish the sharp leading asymptotic

log H(m) ~ log(m) log log(m) / log 2.

The equivalent book capacity is log b(N) ~ log 2 · log N / log log N. The elementary construction works for every number of pages. It rules out a uniform polynomial label bound for all finite two-degenerate graphs, while leaving their proved induced realizability intact. The theorem is not a sharp worst-case height theorem for the whole class.

All twelve frozen source files are published byte-for-byte. The original [README](README.md) gives scope and reproduction commands; FROZEN_SHA256.txt can be checked from this directory. Both proof audits and both normal/optimized checker pairs are included.

Publication replays use python3 -B and python3 -B -O on verify_book_height.py and review/independent_controls.py. All four outputs exactly match the frozen receipts: 4,081,492 primary and 105,750 independent checks in each mode. The proofs establish the asymptotics; finite tests do not replace them.
