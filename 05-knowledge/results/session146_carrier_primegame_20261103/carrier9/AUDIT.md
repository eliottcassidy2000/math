# Independent order-nine endpoint-block audit

The 9-vertex no-isolate signatures are partitioned by endpoint set E. For a representative endpoint set E={0,...,2r−1}, all 2r vertices are touched by each banked exterior context. Independent finite-field calculations give:

| r | Block size | det mod 101 | det mod 103 | deletion and duplicate ranks |
|---|---:|---:|---:|---:|
| 1 | 1 | 41 | 43 | 0 on deletion; duplicate not applicable |
| 2 | 6 | 66 | 74 | 5 |
| 3 | 60 | 90 | 64 | 59 |
| 4 | 420 | 6 | 65 | 419 |

Each block has full rank in both fields. The verifier generates endpoint-pair paths using ordered interior subsets, combines paths by an exact-cover recursion, and independently matches all masks and signs in the producer's catalogue. It checks every context's endpoint-set saturation, edge bounds, lack of loops, duplicate edges, old-old edges, truncations, and trailing tokens. Deletion and duplication checks run over both fields. Hostile controls include deliberately removing every external edge incident to an endpoint; normal and Python -O runs match. See `verify.log` and `hostiles.normal.log`.

A full-support forest has degree one at its endpoints and degree two at all other old vertices. If an exterior context touches every vertex in E, a coefficient vanishes on any full-support signature whose endpoint set does not contain E: an E vertex outside that set would have selected degree at least three. Relabelling the independently verified representative blocks to every E changes gauges only by nonzero signs. Ordering endpoint sets by inclusion makes the full matrix block triangular with these nonsingular diagonal blocks. Its rank is therefore 36·1 + 126·6 + 84·60 + 9·420 = 9,612. The certified order-eight result and lifting argument supply all 56,296 signatures with at least one isolate, so all 65,908 pooled-parity coordinates are observable. Together with the established upper bound, the order-nine unnormalized linear dimension is 65,908, and the normalized-law dimension is 65,907. This concerns all fixed source-oblivious attachments; it does not prove S/B-only minimality or nonlinear observability.

The proof inherits Session145's coefficient formula and upper bound, and the separately certified order-eight minimum and lifting. Exterior coefficient rows are finite rational combinations of fixed-cut evaluations rather than individual attachments.

Reproduce from this directory:

```sh
g++ -O3 -std=c++17 -o verify9 verify9.cpp
for r in 1 2 3 4; do ./verify9 "$r" "block$r.contexts" "block$r.reps"; done
python3 -B hostiles.py > hostiles.normal.log
python3 -B -O hostiles.py > hostiles.optimized.log
cmp hostiles.normal.log hostiles.optimized.log
sha256sum -c SHA256SUMS
```
