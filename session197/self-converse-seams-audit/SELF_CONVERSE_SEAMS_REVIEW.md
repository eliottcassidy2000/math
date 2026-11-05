# Independent audit: self-converse palindromic seam extremizers

Status: **PASS, no correction required.** The following frozen sources were read without modification:

- `session197-lowcycle/SELF_CONVERSE_SEAMS.md`: SHA256 `eb47d9a4f5eb86c7f34480711c5306eecda081f262baf9fea1537952252a6995`.
- `session197-lowcycle/verify_self_converse_seams.py`: `b5cb9d98efea19629c2d3bfe8045d1f329147d5047a08c0009ef966e6c2a7fc8`.
- `session197-lowcycle/self_converse_seams_certificate.json`: `a0ef61aaf5ab043162eac8545b464070ae79d0dfb5481f66e4ffb4a6de6d6c43`.

The named `ROOTED.md` dependency was read to check the one-seam response conventions. For the transitive sink gadget of order d>=2, its constants are indeed g=1, t=0, E=Q=1 and L=B=2^(d-1). Thus the invariant R=H and the factor 1+2^(d-1) follow exactly. The root-start count remains one. The unique reverse cross arc gives a bijection between Hamilton cycles and the unique old root-start path followed by the increasing gadget path ending at the new root. This is valid beginning with a singleton old graph because total order is at least three.

The reflection proof checks every exceptional seam and every ordinary arc. The three displayed palindromic families cover exactly all n>=3 not divisible by six and attain the corresponding exact Busch lower bounds. Hence the stated self-converse strong minima follow at all those orders, without any classification of minimizers.

For n divisible by six, the even self-converse congruence rounds the ordinary lower bound up by two. The n=6 singleton-gap construction attains 17. For n>=12 the asserted upper value is explicitly a family optimum, not the global minimum. The side-gap and central-gap reductions are valid and strictly improve the stated large gaps; splitting a side gap four preserves the product. The residual two integer optimizations give respectively 405*25^(t-2) and 425*25^(t-2), with the center-five case alone feasible at t=1. The 27/25 comparison is exact. No gap in the reduction argument was found.

## Replays and independent physical checks

The source verifier was copied to an audit-only directory because it writes its certificate beside itself. Normal and optimized runs there reproduced the frozen certificate byte for byte. This replays all 143 physical graphs through order twelve, 4,178 palindromic lists through order 32 and the independent family-minimum dynamic program through order 300.

A separately written C++ checker, `check_self_converse_seams.cpp`, constructs the literal adjacency arrays and performs direct depth-first enumeration of Hamilton vertex words, without importing source code or using a subset recurrence. For all 143 gap compositions it independently verifies H, the root-start count, the root-end count and the Hamilton cycle count. It also verifies the physical anti-isomorphism in all 31 palindromic graphs.

Independent checker SHA256: `3d5a2251110ea2672e3501b582764060aedfcd525396584d98537e8911019d67`.
Independent receipt `seams-independent.txt` SHA256: `06ce3f638e926173ab748b4a14b2b74dea77373da5c2753612dd3f4ffa3cdc16`.

The finite controls agree with the all-order proofs. The unrestricted self-converse minimum at orders 6t>=12 remains open in this packet.
