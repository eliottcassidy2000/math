# HYP-2593: LRC(14) S3 Phase-Color Reservoir

**Status:** OPEN proof program with a proved exact reformulation and exact
case-bank evidence.

For an S3 set written as

`S = P union {V-e : e in E}`,

the finite CRT witness at modulus `q=14V` is not merely an uncolored max-gap
problem.  A residue `a` has a phase color `b=a mod 14`.  Define

`C_b(E) = {x : ||e*x - b/14|| >= 1/14 for every e in E}`.

Then, for `a=b+14t` and `x=a/(14V)`,

`a/q` is a level-`1/14` witness for `S`

if and only if

`x in G_P cap C_b(E)`.

Thus the exact finite object is the 14-colored reservoir

`Sigma(P,E) = sum_{b=0}^{13} meas(G_P cap C_b(E))`,

sampled by the colored grids `x=(b+14t)/(14V)`.

## Evidence

Script: `04-computation/lrc14_phase_color_reservoir_codex.py`.
Output: `05-knowledge/results/lrc14_phase_color_reservoir_codex.out`.

The script verifies the exact identity

`actual CRT witness count at q=14V = colored grid hit count`

with `0` mismatches on the named hard lifts tested.

In the structured bank of `123405` cases, the colored continuous reservoir was
large:

- Global minimum `Sigma = 14249/28028 ~= 0.508384`, at
  `P=(1,2,3,5,7,8,9,11,12,13)`, `E=(0,2,3)`.
- Global minimum of the largest single color measure was
  `2941/63063 ~= 0.046636`, at
  `P=(1,2,3,5,7,8,9,11,12,13)`, `E=(0,1,6)`.
- Named via-max-zero shapes are not thin in this colored sense:
  `via_zero_k7` has `Sigma=1063/780 ~= 1.362821`;
  `via_zero_k9` has `Sigma=25369/28665 ~= 0.885017`.

Random named lift stress over `250` covering lifts found `0` zero-witness
rows at `q=14V`; the worst observed ratio

`actual_count / (V*Sigma)`

was `38808/57089 ~= 0.679781`.

## Why This May Be The Last Lemma

The previous uncolored slack-reservoir target still had a placement issue:
once `a mod 14` is fixed, each color grid has spacing `1/V`, not `1/(14V)`.
The phase-color formulation removes that mismatch.  It says the expected
number of actual CRT witnesses is `V*Sigma(P,E)`, up to discrepancy.

So the remaining proof can be stated as:

1. Prove a uniform positive lower bound for `Sigma(P,E)` over admissible S3
   shapes, or reduce large-spread shapes to the relation-lattice Fourier
   tail.
2. Prove a colored Erdos-Turan/Koksma bound showing the `14` progressions
   `x=(b+14t)/(14V)` cannot all miss `G_P cap C_b(E)` when `V*Sigma` is
   sufficiently large.
3. Finish the remaining small-`V` cases by exact finite enumeration.

This reframes the binding denominator obstruction: covering would have to make
all colored grid hits disappear despite a large continuous colored reservoir.

## Tournament Analysis

Vertices are phase colors `b=0,...,13`.  The pairwise observable is aggregate
colored reservoir measure over the structured bank.  The tournament is
transitive (`score_hist=0..13`, no directed 3-cycles, one Hamiltonian path),
with colors `1` and `13` leading symmetrically and color `0` the sink.  This
matches the geometry: color `0` puts the max runner at the observer, while
colors near `1` and `13` sit on the threshold edge.

## Caveat

This does not prove LRC(14).  It turns the residual into a sharply stated
colored discrepancy theorem plus finite base check.  The structured bank is
evidence for the `Sigma` floor, not a proof over all relation-lattice shapes.

See also `HYP-2592`, `HYP-2591`, `HYP-2589`, THM-527/528, and `OPEN-Q-108`.
