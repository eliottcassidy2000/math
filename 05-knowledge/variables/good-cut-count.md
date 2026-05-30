# Variable: good-cut count

**Symbol:** `g(τ)` or `goodCutCount(τ)`
**Type:** integer statistic of a base-path tiling
**Defined in:** THM-330, THM-336, `TournamentH7.GoodCuts`

## Definition

For a tiling `τ` on vertices `{0,...,n-1}`, a tile `(hi,lo)` is upward when
`lo -> hi`. The tile crosses every cut

`k ∈ {lo+1, ..., hi}`.

The good-cut set is

`G(τ) = { k ∈ {1,...,n-1} : some upward tile crosses k }`.

The good-cut count is

`g(τ) = |G(τ)|`.

## Values

Exact tiling counts by `g`:

| n | counts |
|---|--------|
| 3 | `{0: 1, 2: 1}` |
| 4 | `{0: 1, 2: 2, 3: 5}` |
| 5 | `{0: 1, 2: 3, 3: 10, 4: 50}` |
| 6 | `{0: 1, 2: 4, 3: 15, 4: 101, 5: 903}` |
| 7 | `{0: 1, 2: 5, 3: 20, 4: 153, 5: 1816, 6: 30773}` |
| 8 | `{0: 1, 2: 6, 3: 25, 4: 206, 5: 2739, 6: 61671, 7: 2032504}` |

Sources: `05-knowledge/results/goodcut_bucket_merged_s13.out` and
`05-knowledge/results/goodcut_bucket_n7_fast_s1.out`; recurrence extension from
`05-knowledge/results/goodcut_interval_union_s15.out`.

## Equations

**Structural interpretation (THM-354):** for any Hamiltonian base path `P` of a
tournament `T`,

`g_P(T) = n - #SCC(T)`.

Equivalently, bad cuts are exactly the boundaries between consecutive strong
components in the condensation order.

- `g(τ) = 0` iff every tile is downward.
- `g(τ) ≠ 1` for every tiling. Each upward tile spans at least two cuts.
- `g(τ) ≤ n-1`.
- `g(τ) ∈ {0} ∪ {2,...,n-1}`.
- `g(τ)=n-1` iff every legal cut is good.
- `g(reflect τ) = g(τ)`, where grid reflection maps cut `k` to cut `n-k`.
- THM-330: `τ` is strongly connected iff `g(τ) = n-1`.
- Monotonicity: if every upward tile of `τ` is also upward in `τ'`, then
  `g(τ) ≤ g(τ')`.
- Interval form: `k` is good iff `k` lies in the interval `{lo+1,...,hi}`
  of some upward tile `(hi,lo)`.
- THM-349: the full bucket distribution is the interval-union recurrence
  `B_N(x)=B_{N-1}(x)+Σ_{L=2}^N c_L x^L B_{N-L-1}(x)`, where `c_L`
  counts covers of one connected run by intervals of length at least 2.

## Relationships

- Related to `delta_proj`: `g` is a cheap tiling-coordinate projection whose
  residue is the strong-component defect `n - #SCC(T)`.
- Related to projection-defect polarity: S14 found that single-tile lines with
  `|Delta g|>0` are never even-only through n=6; even-only defects live inside
  `g`-neutral lines.
- Related to `H(T)`: bucket extremes are transitive (`g=0`, `H=1`) and strongly connected (`g=n-1`, broad H range).
- Related to `G_n/Z_2`: THM-354 proves that `g` is constant on each tournament
  isomorphism class and on complement-merged classes, because SCC count is
  preserved by isomorphism and by taking the opposite tournament.
- Related to `B_N(x)`: the polynomial coefficient `[x^g]B_{n-1}` counts all tilings with `g(τ)=g`.

## Lean Formalisation

`04-computation/lean/TournamentH7/TournamentH7/GoodCuts.lean` proves:

- `StTiling.goodCuts_empty_iff_all_down`
- `StTiling.goodCutCount_ne_one`
- `StTiling.goodCutCount_reflect`
- `StTiling.isGoodCut_iff_exists_upward_tile_interval`
- `StTiling.mem_goodCuts_iff_exists_upward_tile_interval`
- `StTiling.cutInterval_subset_goodCuts_of_upward_tile`
- `StTiling.goodCutCount_mono`
- `StTiling.goodCutCount_bucket_bounds`
- `StTiling.goodCutCount_eq_top_iff_all_cuts_good`

These build without project-specific axioms; the audit shows only Lean foundations.

## Verification

`04-computation/goodcut_scc_defect_s354.py` exhaustively verifies THM-354 for
all labeled tournaments through `n=6`, all fixed-base tilings through `n=7`,
and sampled labeled/fixed-base cases at `n=7,8`.

Result: `05-knowledge/results/goodcut_scc_defect_s354.out`.
