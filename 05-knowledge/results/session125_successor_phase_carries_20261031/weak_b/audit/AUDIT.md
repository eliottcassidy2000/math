# Independent audit: weak-B / ninja Q25 path carries

Session 125, 2026-10-31. **PASS.** Separate from the frozen original-port Q49
packet. Only this auditor's scratch directory was edited. No producer module
was imported; no external C was read or executed for this extension. No
permission denial or unresolved blocker occurred.

## 1. Correct mathematical object

A weak-B path starts `1,3`, ends at 8 for even order and 10 for odd order.
A pair at orders n,n+1 additionally agrees on all shared-label position
parities. The necessary port phases are 1 and 10 odd, 3 and 8 even.

These are **path banks**, not direct cycle banks: `1+10=11` is not square.
The original-port Q49 cycle theorem and its six-defect bound are not silently
transferred to this different convention. In particular, any permanent weak-B
phase already differs from numerical parity at labels 3 and 10.

## 2. Independently rebuilt finite phase obstruction

The standalone auditor constructs square-sum edges by square enumeration and
enumerates all paths with prefix 1,3 by endpoint extension. Residual degree
and connectivity pruning are necessary conditions, so completeness is
preserved. The eligible spectra independently give:

| Order | Paths | Masks |
|---:|---:|---:|
|42|49|39|
|43|74|30|

The six common masks and shorter/longer multiplicities are:

| Mask | Shorter | Longer | Obstructed order | Only allowed neighbor |
|---:|---:|---:|---:|---:|
|1585793692427|1|2|46|18|
|1599483269979|1|1|46|35|
|1602705117979|2|1|47|34|
|3372531545867|1|1|46|18|
|3389445068555|4|3|47|34|
|3526867243787|1|2|47|34|

Thus there are exactly 20 compatible path pairs. The independent results
agree with the separately owned onset catalogue. Its positive witnesses
are revalidated, but completeness here does not rely on importing its solver
or accepting six examples as exhaustive.

If phases persist, slot counting forces every new label v>42 to have phase
`v mod 2`. Vertex 46 has arithmetic neighbors 3,18,35; vertex 47 has neighbors
2,17,34. Filtering by each fixed mask leaves the sole displayed neighbor.
The arriving vertex is internal, so a path through it would need degree two.
All six masks therefore fail by order 47. **No claim is made that each listed
failure is the first failure of its mask.** No larger path census was needed.

Hostile controls reject a duplicate, wrong prefix, nonsquare edge, and two
individually valid paths with incompatible phase masks.

## 3. Protected head and the 61-edge boundary

Let `R(P)=E(P)\{1--3}` and `T={1,3,8,10}`. The only square-sum pairs with
both endpoints in T are 1--3 and 1--8. The first is removed from R(P); the
second is absent because the start vertex 1 already uses its unique path edge
1--3. Hence no edge of R(P) has both endpoints in T.

At signed offset c, the Q25 transform sends a label x in position j to
`25x+(-1)^j c`, with positions starting at 1. Splitting only the protected
head edge changes no lifted R(P) edge. Nearest-center decoding
`floor((z+12)/25)` recovers the input label; singleton labels 1..12 decode
to 0. Every connector or retained head carrier has both decoded endpoints
in `{0,1,3,8,10}`, whereas no lifted R(P) edge does.

Thus moving-body edges and boundary edges remain disjoint even when input
phases, head-cut sets, or A/B word choices differ. With h head cuts there
are `25+12+h` atoms, hence `36+h` connectors, and `25-h` uncut head carriers:
the boundary has exactly **61 edges**, independently of h.

The auditor verifies this separation differently from the producer. Rather
than rebuilding symbolic endpoints, it expands actual paths and extracts
the boundary as precisely those edges whose decoded endpoints are both in
`{0,1,3,8,10}`. All 200 outputs from the 100 paired states give 61 boundary
edges. This also independently determines every boundary deletion constant.

## 4. Exact edit laws and finite constants

The member threshold changes exactly one offset between the shorter and
longer output. Shared input phase preserves the map on common labels; the
common 1--3 edge contributes no difference. Therefore

`d(n)=d(parent(n))+K_pair`, with `0<=K_pair<=61`.

Noncarry representatives have identical moving bodies, so `H(t)<=61`.
At a carry `t=25r+12`, r>=43, the full signed-offset union of a parent body
edge `{x,y}` is

`{{25x+c,25y-c}: -12<=c<=12}`.

It is independent of the input phase. Unique decoding separates different
parent edges, and boundary disjointness prevents cancellation. Since removing
the common head edge leaves H unchanged,

`H(25r+12)=25H(r)+K_carry`, with `0<=K_carry<=61`.

The audit sequentially decodes the owned compressed root bytes instead of
using the producer's random-access index formula. All 1,020 base pairs at
42..1061 pass coverage, square-edge, prefix, endpoint and shared-phase tests.
All 100 A/B paired states are independently expanded and checked. Results:

- Maximum base deletion count: 813 at order 1054.
- A within-pair K range: 4..22; B range: 3..22.
- A noncarry maximum: 37; B maximum: 36.
- Both A/B carry boundary sets are **equal** in both parent parities, so
  their K_carry values are exactly zero.
- Base gaps H(43),H(44),H(45),H(46) are 0,9,14,14.

Both homogeneous operators give the independently expanded rays:

| Seed | Orders | Gaps |
|---:|---|---|
|44|44,1112,27812|9,225,5625|
|45|45,1137,28437|14,350,8750|

At 28436→28437 the canonical path deletion counts are 8787 for A and 8769
for B; the corresponding within-pair counts are only 42 and 24. All eight
possible bank-selector triples around each sampled carry were also checked.

The general operator-validity theorem is the inherited owned construction;
this audit independently checks its finite data and the new edit identities.

## 5. Universal and pinned limits: quantifiers checked

If H(43),H(44),H(45),H(46) were all zero, equality of fixed-start paths would
glue Pair(42)..Pair(46) into the impossible coherent family 42..47. Thus one
of these seed gaps is positive. Two Hamiltonian paths with the same endpoints
cannot differ by a single deleted edge: equal degrees would force its
replacement to be the same edge. One seed gap is therefore at least 2.

Along `t_k=25^k(t0+1/2)-1/2`, the carry law implies `H(t_k)>=2*25^k`.
Within-pair deletions are O(log n). Directed set-difference triangles give

`H(n+1)-d(n)-1 <= delta(S_n,S_(n+1)) <= H(n+1)+d(n)`.

For a selector choosing either S_n or L_(n-1), at a carry r the middle
choice S_r forces an incoming cost at least `H(r)-d(r-1)-62`; choosing
L_(r-1) forces an outgoing cost at least `H(r)-d(r)-62`. The neighboring
representation gaps are at most 61; reversing an adjacent-order comparison
accounts for the additional 1. Hence every such path-bank selector has

`limsup delta(P_n,P_(n+1))/n >= 2/(46+1/2)=4/93`.

For either pinned homogeneous operator, seed 45 gives
`H(t_k)=14*25^k`, `t_k=(91*25^k-1)/2`. The canonical predecessor-step
ratio tends to **28/91**. For every two-bank selector, the maximum of its
two adjacent step costs around t_k, divided by t_k, also tends to 28/91:
the preceding case split supplies the lower bound, while the uniform upper
bound `H(n)+d(n)+H(n+1)` supplies the matching upper bound.

**Scope:** this is a ray limit, not a global maximizing-seed computation;
the universal lower bound concerns the stated Q25 head-cut architecture,
not all weak-B path families or all cycle families. No cycle-consumer
transport theorem, extra constructor search, or unbounded census is claimed.

## 6. Replay and receipts

```sh
python3 independent_audit.py \
  --operator-dir /path/to/session124_mixed_radix_square_rays_20261031 \
  --base-dir /path/to/square_sum_B_base_interval_20261019 \
  --small-bases /path/to/SMALL_BASES.json \
  --weak-b-catalogue /path/to/weak_b_certificate.json \
  --output normal.json
```

Run from this audit directory, or use the script by absolute path. All four
input options have source-relative published defaults and may be omitted in
the repository layout. Repeat with `python3 -O` and `--output optimized.json`, then compare the
receipts. Checks use explicit exceptions, never optimization-sensitive
assertions. Receipt `source_sha256` records portable result-relative names for every used owned input, including
all 100 paired certificates and the compressed chunks. Compressed and decoded
chunk pins, headers, lengths and manifests are checked. The final producer
snapshot is recorded separately in `AUDIT_INPUTS.json`; `SHA256SUMS` covers
the frozen audit artifacts. The earlier Q49 audit is unchanged.
