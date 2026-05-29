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

Source: `05-knowledge/results/goodcut_bucket_merged_s13.out`.

## Equations

- `g(τ) = 0` iff every tile is downward.
- `g(τ) ≠ 1` for every tiling. Each upward tile spans at least two cuts.
- `g(reflect τ) = g(τ)`, where grid reflection maps cut `k` to cut `n-k`.
- THM-330: `τ` is strongly connected iff `g(τ) = n-1`.

## Relationships

- Related to `delta_proj`: `g` is a cheap tiling-coordinate projection that appears to descend to merged tournament classes through n=6.
- Related to `H(T)`: bucket extremes are transitive (`g=0`, `H=1`) and strongly connected (`g=n-1`, broad H range).
- Related to `G_n/Z_2`: HYP-1764 conjectures that `g` is constant on each merged tournament class.

## Lean Formalisation

`04-computation/lean/TournamentH7/TournamentH7/GoodCuts.lean` proves:

- `StTiling.goodCuts_empty_iff_all_down`
- `StTiling.goodCutCount_ne_one`
- `StTiling.goodCutCount_reflect`

All three build without project-specific axioms; the audit shows only Lean foundations.
