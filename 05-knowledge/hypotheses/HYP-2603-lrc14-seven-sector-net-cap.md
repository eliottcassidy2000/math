# HYP-2603: LRC(14) Seven-Sector Net-Cap Route

**Status:** open proof program with exact bounded evidence; not a proof of LRC(14).

Script: `04-computation/lrc14_mu17_sector_cover_codex.py`  
Output: `05-knowledge/results/lrc14_mu17_sector_cover_codex.out`

## Claim

For `8 <= k <= 13`, let

`N(E) = {x in [0,1) : all cyclic gaps of {e*x mod 1 : e in E} are <= 1/7}`.

Recent work reduced the S3 residual to proving

`meas(N(E)) <= cap_k = min_{|P|=13-k} meas(G_P)`.

A sufficient, weaker target is the seven-sector cap:

`meas(S7(E)) <= cap_k`,

where `S7(E)` is the set of `x` for which the points `{e*x mod 1 : e in E}`
hit all seven fixed sectors `[j/7,(j+1)/7)`, `j=0,...,6`.

The inclusion is elementary up to boundary measure zero:

`N(E) subset S7(E)`.

Indeed, if one fixed sector of length `1/7` contains no point, then the
configuration has an empty arc of length at least `1/7`; strict boundary
exceptions have measure zero.  Therefore a global sector-cover cap would imply
`mu_{1/7}(E) >= 1 - cap_k`, and hence would close the HYP-2602 endpoint.

## Evidence

The exact bank scanned primitive normalized shapes `0 in E` using exact
rational sector breakpoints.

| k | bank | max `meas(S7(E))` | maximizer | `cap_k` | margin |
|---|------|-------------------|-----------|---------|--------|
| 8 | `max(E)<=14`, `3431` rows | `481/1470` | `(0,1,2,3,4,5,6,7)` | `2243/5880` | `319/5880` |
| 9 | `max(E)<=14`, `3003` rows | `2447/5880` | `(0,1,2,3,4,5,6,7,8)` | `1979/4004` | `65669/840840` |
| 10 | `max(E)<=15`, `5005` rows | `8899/17640` | `(0,1,2,3,4,5,6,7,8,9)` | `55/91` | `22913/229320` |
| 11 | `max(E)<=14`, `1001` rows | `3419/5880` | `(0,1,2,3,4,5,6,7,8,9,10)` | `66/91` | `10993/76440` |
| 12 | `max(E)<=14`, `364` rows | `11381/17640` | `(0,1,2,3,4,5,6,7,8,9,10,12)` | `6/7` | `3739/17640` |
| 13 | `max(E)<=14`, `91` rows | `12023/17640` | `(0,1,2,3,4,5,6,7,8,9,10,12,14)` | `1` | `5617/17640` |

No sector-cap violations appeared.  The margins are much larger than the
corresponding net-measure margins, which suggests that fixed-sector cover may
be a forgiving certificate even when `mu_{1/7}` itself is delicate.

## Important Negative Result

Naive one-step compression is not a proof lemma.

For `mu_{1/7}`, compression can increase `mu`, for example:

`(0,1,2,3,4,6,7,13) -> (0,1,2,3,4,5,6,12)`

raises `mu` from `62113/63063` to `1`.

For sector cover, compression can decrease the proposed upper-bound object:

`(0,1,2,3,4,5,6,10) -> (0,1,2,3,4,5,6,9)`

lowers `meas(S7)` from `4/21` to `79/420`.

Thus the proof should not try to prove monotonicity under a local gap move.
The plausible route is global: sector occupancy, relation-height splitting, or
a finite low-height reduction plus high-height discrepancy tail.

## Relation to Other Agents' Work

HYP-2602 asks for a direct bound on `N(E)`, or equivalently a lower bound on
`mu_{1/7}`.  HYP-2603 weakens that endpoint by replacing the exact net set
with a larger fixed-sector cover set.

HYP-2600/HYP-2601 and HYP-2599 suggest how to attack the global part:
high relation height should make fixed-sector occupancy close to an iid
sector-hit model, while low relation height leaves finitely many affine
patterns to certify exactly.

HYP-2593/HYP-2595 remain the integer-placement layer after the continuous
reservoir is known.  The sector route is continuous and LRC-free; it still
needs to be paired with the colored CRT layer before it becomes a runner proof.

## Tournament Analysis

Vertices were proof routes, not runners:

`sector_cover_cap`, `consecutive_mu_min`, `gap_max_stratification`,
`relation_height_tail`, `one_step_compression`.

The observable was current evidence/risk score.  The switch oriented stronger
routes toward weaker routes, with the listed order as the tie Hamiltonian path.
The tournament was transitive:

- score histogram: one vertex at each score `0,1,2,3,4`;
- directed `3`-cycles: `0`;
- SCC sizes: all `1`;
- Hamiltonian path count: `1`;
- leader: `sector_cover_cap`;
- sink: `one_step_compression`.

Assumption challenged: the quotient preserves the predicate "a 1/7-net must
hit every fixed sector" and destroys cyclic adjacency, which runner creates
which boundary point, and equality-wall timing.  Those destroyed data may still
be needed for a final proof, but they are not needed for this sufficient cap.

## Next Lemmas

1. Prove or refute the global sector cap
   `meas(S7(E)) <= cap_k` for all normalized `E`, `8 <= k <= 13`.
2. Split by affine relation height.  High-height rows should obey a sector
   discrepancy tail; low-height rows should reduce to finitely many patterns.
3. Try a dual formulation: if fixed sectors are all hit too often, extract a
   small affine relation or a repeated sector-word obstruction.
4. Reuse HYP-2601's high-relation-height certificate machinery, but with
   seven sector indicators replacing the smooth empty-arc minorant.
