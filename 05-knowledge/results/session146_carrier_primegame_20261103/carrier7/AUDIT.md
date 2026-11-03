# Independent order-seven fixed-attachment rank audit

The frozen bank assembled from `contexts7.part00` and `contexts7.part01` supplies 2,668 exterior coefficient-extraction contexts. The independently generated representative catalogue and recomputed matrices have the following ranks and determinants:

| Field | Rank | Determinant |
| --- | ---: | ---: |
| GF(101) | 2,668 | 48 |
| GF(103) | 2,668 | 66 |

Replacing the last row with the first gives rank 2,667 over GF(101). Deleting the last row gives rank 2,667 over each field. The bank SHA-256 is `ee0c332d4dd580eb3073cf95eba3ca098d005cd1f3f26e14c0aaaa9c6671baef`. This is a nonzero rational minor because the coefficients have only powers of two in their denominators. Together with Session145's all-order pooled-parity upper bound, it proves that the unnormalized linear response dimension is 2,668 at order seven. On normalized laws the corresponding number of nonconstant coordinates is 2,667. The target is the full language of fixed source-oblivious attachments. This result does not establish the same rank for S/B-only operations or arbitrary nonlinear protocols.

Each exterior coefficient row is a finite rational linear combination of evaluations at legitimate fixed-cut attachments: for exterior y_e=±1/2, Walsh extraction multiplies the average of H_N times the selected signs by 2^|G|. It is not one single fixed-cut response. The context validator rejects old-old edges, loops, reversed duplicate edges, out-of-range endpoints, missing endpoints, and trailing tokens. The independent verifier enumerates path-block splittings of all permutations to reconstruct and order all 2,668 signatures and least-mask representatives, checking both masks and canonical signs against the supplied catalogue. It evaluates each matrix entry via the even-linear-forest Hamilton coefficient formula, detects cycles with disjoint sets, and performs modular elimination. Direct permutation expansion tests nine edge sets on eight vertices, including nonforest and odd-path cases, and an exact Walsh extraction checks a nonzero coefficient -2 on four vertices. Controls run identically with Python normal and -O.

The verifier also independently reconstructs the previous n6 certificate and its determinants 40 mod 101 and 84 mod 103; its duplicate-row control yields rank 570. The n6 contexts derive from Session145's independent audit. The current audit inherits, rather than reproves, Session145's general closure/compression theorem and coefficient derivation.

Reproduce from this directory:

```sh
cat contexts7.part00 contexts7.part01 > contexts7.txt
g++ -O3 -std=c++17 -o verify verify.cpp
./verify 6 rank_contexts6.txt reps6.txt
./verify 7 contexts7.txt reps7.txt
./verify 7 contexts7.txt reps7.txt --delete-last
python3 -B controls.py > controls.normal.txt
python3 -B -O controls.py > controls.optimized.txt
cmp controls.normal.txt controls.optimized.txt
sha256sum -c SHA256SUMS
```

If distributing the order-seven bank in numbered chunks, concatenate them in the fixed prescribed order into a temporary file, check its SHA-256 against the value above, and pass that file as the second argument of `./verify 7`. `reps7.txt` is checked against the verifier's independent construction; it is not trusted as a source of representative correctness. Bank provenance: `generate_rank7.cpp`, deterministic seed 146146, trial 36,425. The frozen replay runs in a temporary directory.
