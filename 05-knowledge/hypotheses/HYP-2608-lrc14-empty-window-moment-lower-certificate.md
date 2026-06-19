# HYP-2608: LRC(14) Empty-Window Moment Lower Certificate

**Status:** open proof program with a proved per-set lower-bound certificate and
exact bounded evidence.  This does not prove LRC(14).

Script: `04-computation/lrc14_empty_window_moment_lower_codex_0618.py`  
Output: `05-knowledge/results/lrc14_empty_window_moment_lower_codex_0618.out`

## Claim

Let `A={j/14 : j=0,...,6}` and let `W_j(E)` be the set of times `x` for which
the open arc `(j/14, j/14+1/7)` is empty of all points `{e*x mod 1 : e in E}`.
The previous EWLB route asks for a lower bound on

`EWLB_A(E) = meas(union_j W_j(E))`.

This hypothesis replaces the raw union by a low-degree moment lower certificate.
Set

`R(x) = # {j : x in W_j(E)}`,

and define empty-window factorial moments

`T_s(E) = E_x[C(R(x),s)]`.

If a polynomial in the binomial basis

`h(r)=sum_s z_s C(r,s)`

satisfies `h(0)<=0` and `h(r)<=1` for `r=1,...,7`, then

`EWLB_A(E) = P(R>0) >= Phi_z(E) := sum_s z_s T_s(E)`.

This is the exact lower-dual analogue of THM-534's moment majorant for the
seven-sector cover set.  The useful certificates are:

| row | certificate | `h(0..7)` |
|---|---|---|
| `R=1` | `z=(0, 1/7)` | `0,1/7,2/7,3/7,4/7,5/7,6/7,1` |
| `R=2` | `z=(0, 2/5, -1/10)` | `0,2/5,7/10,9/10,1,1,9/10,7/10` |
| `R=3` | `z=(0, 5/7, -3/7, 1/7)` | `0,5/7,1,1,6/7,5/7,5/7,1` |
| `R=4` | `z=(0, 16/21, -11/21, 2/7, -2/21)` | `0,16/21,1,1,20/21,20/21,1,1` |

The first degree that clears the consecutive/AP threshold is:

| `k` | degree | `Phi(AP_k)` | threshold | margin |
|---:|---:|---:|---:|---:|
| 8 | 4 | `19633/30870` | `3637/5880` | `431/24696` |
| 9 | 3 | `127/245` | `2025/4004` | `1769/140140` |
| 10 | 3 | `10765/24696` | `36/91` | `12937/321048` |
| 11 | 2 | `21673/70560` | `25/91` | `29749/917280` |
| 12 | 1 | `89533/543312` | `1/7` | `11917/543312` |

Therefore a sufficient endpoint for HYP-2602 is the scalar rearrangement

`Phi_k(E) >= threshold_k` for the row-specific certificate above.

## Evidence

Exact bounded primitive sweeps give zero threshold failures:

| `k` | certificate | bank | minimum `Phi_k` | minimizer | below threshold |
|---:|---:|---:|---:|---|---:|
| 8 | degree 4 | `maxE<=14`, `3431` rows | `19633/30870` | `(0,1,2,3,4,5,6,7)` | `0` |
| 9 | degree 3 | `maxE<=13`, `1287` rows | `127/245` | `(0,1,2,3,4,5,6,7,8)` | `0` |
| 10 | degree 3 | `maxE<=13`, `715` rows | `10765/24696` | `(0,1,2,3,4,5,6,7,8,9)` | `0` |
| 11 | degree 2 | `maxE<=13`, `286` rows | `21673/70560` | `(0,1,2,3,4,5,6,7,8,9,10)` | `0` |
| 12 | degree 1 | `maxE<=13`, `78` rows | `263533/1681680` | `(0,1,3,4,5,6,7,8,9,10,11,13)` | `0` |

The dangerous rows `k=8..11` are minimized by the consecutive cluster in these
exact banks.  The `k=12` minimizer is non-AP, but the row has enough slack that
the degree-1 average still clears the threshold.

Small exact random stress through speeds up to `35` found only large positive
margins.  The worst reported random margins are already much larger than the AP
margin for `k=8..11`, consistent with the thesis that AP-rich coupling is the
hard case.

## Proof Route

The proof route is now:

1. The polynomial minorants above prove `EWLB_A(E) >= Phi_k(E)` for every `E`.
2. Prove the scalar rearrangement `Phi_k(E) >= threshold_k`.
3. Combine with the already-proved reduction `mu_{1/7}(E) >= EWLB_A(E)`.

For `k=8..11`, the exact banks suggest the sharper statement

`Phi_k(E) >= Phi_k(AP_k)`.

For `k=12`, AP extremality is not needed; a weaker first-moment lower bound
`T_1(E)/7 >= 1/7` plus slack, or a coarse large-spread/finite-residual split,
should suffice.

## Relation To THM-534

THM-534 upper-bounds the bad cover event `P(N=0)` where `N` is the number of
missed fixed sectors among `{1,...,6}`.  HYP-2608 lower-bounds the good
empty-window event `P(R>0)` where `R` counts anchored empty arcs of length
`1/7`.  The two objects are dual in spirit:

- THM-534 uses a polynomial majorant of `1[t=0]`;
- HYP-2608 uses a polynomial minorant of `1[t>0]`.

Both reduce LRC(14)'s continuous sector problem to a low-degree scalar moment
functional.  HYP-2608 is especially close to the user's region-first picture:
the vertices are loop regions/empty windows, not runners.

## Tournament Analysis

Vertices are the seven anchored empty-window regions `W_0,...,W_6`.

The observable is the certificate load assigned to each region by equal
Shapley splitting of `h(R(x))` among active empty windows at time `x`.
The switch orients `i -> j` when region `i` carries at least as much load as
region `j`; exact ties follow the Hamiltonian path `0 -> 1 -> ... -> 6`.

For every AP row `k=8..12`, the region-load tournament is transitive:
`score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}`, directed `3`-cycles `0`, singleton
SCCs, and Hamiltonian path count `1`.

This quotient preserves which fixed loop regions carry the lower certificate.
It destroys runner labels, boundary timing, and which pair of runner phases
creates a given empty arc.  The challenged assumption is that LRC tournament
vertices must be runners; in this route the natural vertices are fixed regions
or proof obligations.

## Caveat

The script proves the lower-bound implication `EWLB_A(E) >= Phi_k(E)` and
verifies the scalar threshold in bounded banks.  It does not prove
`Phi_k(E) >= threshold_k` for all integer sets, nor the full LRC(14).
The next target is a rearrangement/three-distance proof that AP-rich rows
minimize the low-degree empty-window certificate.

See also HYP-2603, HYP-2604, THM-531, THM-532, THM-534, and OPEN-Q-108.
