# HYP-2594: LRC(14) S3 Colored Discrepancy Boundary

**Status:** OPEN proof program with a proved interval-component bound and
strong exact-grid evidence.

Building on HYP-2593, write the exact `q=14V` witness count as colored grid
hits in

`G_P cap C_b(E)`, with `x=(b+14t)/(14V)`.

Let

`Sigma(P,E)=sum_b meas(G_P cap C_b(E))`

and let

`K(P,E)=sum_b #components(G_P cap C_b(E))`.

For each color, the grid is a shifted mesh-`1/V` progression.  A union of `m`
open intervals of total length `L` contains at least `V*L-m` points of such a
grid.  Therefore

`actual_count(P,E,V) >= V*Sigma(P,E) - K(P,E)`.

Consequently every fixed shape `(P,E)` is certified at `q=14V` for all

`V > K(P,E)/Sigma(P,E)`.

## Evidence

Script: `04-computation/lrc14_colored_discrepancy_bound_codex.py`.
Output: `05-knowledge/results/lrc14_colored_discrepancy_bound_codex.out`.

Named hard shapes have small rigorous cutoffs:

- `via_zero_k7`: `Sigma=1063/780`, `K=124`, cutoff `V>=91`.
- `via_zero_k9`: `Sigma=25369/28665`, `K=104`, cutoff `V>=118`.
- `broad_1_90`: `Sigma=17239/17640`, `K=122`, cutoff `V>=125`.
- `near_via_min`: `Sigma=57089/64680`, `K=114`, cutoff `V>=130`.
- `quarter_min`: `Sigma=13067/17640`, `K=100`, cutoff `V>=135`.

In the structured bank of `123405` cases, the worst component cutoff was

`K/Sigma = 9729720/39911 ~= 243.785`,

attained at

`P=(1,2,3,7,8,9,10,11,12,13)`, `E=(0,1,6)`,

with `Sigma=39911/60060`, `K=162`.  Thus this entire structured bank is
formally certified for `V>=244` by the crude interval-component lemma.

The global minimum `Sigma` case from HYP-2593,

`P=(1,2,3,5,7,8,9,11,12,13)`, `E=(0,2,3)`,

has `K=112`, so the positive-mass floor is not accompanied by large boundary
complexity there.

## What Broke

The same crude component bound is too pessimistic for very wide random tails.
In a random large-spread sample, the worst row had `K=8226` and cutoff
`4128`, about `9.34*(maxE+14)`.  This does not close the full unbounded
tail by itself.

However, the exact `q=14V` grid discrepancy appears much smaller than `K`.
On `160` random primitive covering tail rows:

- zero actual witnesses: `0`;
- minimum `actual/(V*Sigma) = 61776/98035 ~= 0.630142`;
- maximum positive deficit `V*Sigma-actual ~= 71.867`;
- `55` rows had negative deficit, meaning actual count exceeded `V*Sigma`.

So the next sharpened conjecture is that the true colored-grid discrepancy is
bounded by a constant, or at least by a small function of `k` and the small
part `P`, rather than by the raw number of continuous interval components.

## Proof Program

1. Use the proved `V*Sigma-K` lemma to close bounded-spread banks above their
   explicit cutoffs, with finite exact checks below.
2. Replace `K` in the large-spread tail with a residue/Fourier discrepancy
   bound that sees the arithmetic alignment of endpoints with the grid.
3. Combine with the HYP-2593 uniform `Sigma` floor.

This is not yet a proof of LRC(14).  It compresses the finite CRT placement
problem to a boundary-complexity theorem and exposes the raw component count
as an overcounting proxy.

## Tournament Analysis

Vertices are phase colors.  The pairwise observable is boundary efficiency:
aggregate colored mass divided by aggregate component count.  The tournament
is again transitive (`score_hist=0..13`, no directed `3`-cycles, one
Hamiltonian path).  The leaders are colors `3` and `11`, unlike the pure
mass tournament of HYP-2593 where `1` and `13` led.  This says the most
massive colors are not the most boundary-efficient colors; boundary geometry
is a genuinely new quotient.

See also `HYP-2593`, `HYP-2592`, `HYP-2591`, THM-527/528, and `OPEN-Q-108`.
