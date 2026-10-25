# Session56: overlap-safe collision-signature queries

Read THEOREM.md for proofs and scope, certificate.json for exact finite checks.

Reproduce from any directory:

    python3 /path/to/this/verify.py
    python3 /path/to/this/verify.py --verify

Python standard library only. Run without -O: assertions implement the verification gates, and an explicit startup guard rejects optimized mode. The certificate is deterministic. The frozen Session55 source under sources supplies both pivot enumeration and an independent line-equation support enumerator. SOURCE_PINS.txt gives inherited source hashes; SHA256SUMS covers this packet.

Core result: a set family is idempotent under meet convolution precisely when it is meet-closed; a finite reduced antichain is idempotent precisely when empty or singleton. Ordinary overlapping sparse tables therefore cannot automatically aggregate CRT alternative signatures. The actual ring Z/30Z, components [2,5,3], supplies a false-positive distinctness witness if overlap is treated as independence. Disjoint sparse tables and segment trees avoid the error.

Tests: 256 families on three bits; 81 integer histograms on two bits; 65,536 compression pairs; 8,000 antichain associativity triples; 45,760 order-sensitive queries; 31 exact CRT intervals; 28 fixed-mask overlap controls; duplicate-identity rejection. Weighted histograms count residue arrays, not choices of square roots. Sources and bounds are separated from original proofs.

The packet is owned research output only; no remote writes were performed. Initial workspace write failed because the overlay was full; complete artifacts were regenerated in /tmp, a separate filesystem. Targeted processor startup reported no .machine-id; no registration was invented. Parent publication handles current branch integration and remote checkpoints.
