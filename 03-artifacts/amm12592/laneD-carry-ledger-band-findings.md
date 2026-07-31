# Lane D findings: carry-chain construction for gamma < 1 (AMM 12592 / HYP-9061)

Referee: `laneD_referee.py`; output: `laneD_referee.out` (this directory).
All asserted items ran in exact integer/Fraction arithmetic. Frame = THM-2966
spine normal form (cited): depth law `d_m = floor(gamma*m) + D0`, envelope
`T(m) = m + 1 + d_m = (1+gamma)m + O(1)`, i.e. `C = 1 + gamma`. Cell (m,k)
on side 0 owns monomial `p^{m+d_m-k} q^{k+1}`; side 1 mirrored; all row-m
cells sit on the anti-diagonal `z + o = A_m = m + d_m + 1 = T(m)`. Doubled
deficit `delta` at each cell is an integer with `|delta| <= binom(d_m,k)` and
`delta == binom(d_m,k) (mod 2)`; feasibility of a fair extractor at this
deadline is equivalent to `D_M(p) = (1/2) sum delta * p^z q^o -> 0` on (0,1).

## 1. Ledger geometry (PROVED; verified exactly for m <= 400, 9 depth laws)

On anti-diagonal `A_m`, the cells cover positions (q-exponents)
`o in [1, d_m+1] union [m, A_m - 1]`; the **band** `[d_m+2, m-1]`
(directions theta = o/A in (gamma/(1+gamma), 1/(1+gamma))) contains **no
cell of any row, ever**. Consequences, all verified:

- **Band birth**: the band is nonempty iff `m >= d_m + 3`, i.e. from
  `m* ~ (D0+2)/(1-gamma)` on. For gamma = 1 the band never opens
  (this is the structural home of the ratio-2 / dyadic construction).
- **Pairing horizon**: two cells (from the two sides) share a monomial iff
  `m <= d_m + 1`, i.e. `m <~ (D0+1)/(1-gamma)`. Beyond it no two cells of
  the whole scheme ever share a monomial on the deep side.
- **Corner isolation**: the corner cell `k = d_m` (monomial
  `(z,o) = (m, d_m+1)`, capacity `binom(d_m,d_m) = 1`, odd) is attained by
  no other cell of any row once `m >= d_m + 2`. So beyond the pairing
  horizon **every row emits a forced, uncancellable-in-place half-token**
  at its corner; its value must be cancelled collectively by other rows
  through Pascal descent (`p^z q^o = sum_i binom(t,i) p^{z+i} q^{o+t-i}`),
  whose kernel is positive and drifts toward the band center theta = 1/2.
- Per-row exact cancellation is impossible for every depth law (odd cells
  exist in every row by Lucas); backlog between rows is unavoidable. The
  gamma = 1 vs gamma < 1 dichotomy is exactly: at gamma = 1 all backlog
  stays forever editable (no band); at gamma < 1 value that descends into
  the band can never again be touched by any cell.

## 2. Forced-parity closed form (PROVED; asserted at every step of every run)

The parity of the doubled homogenized ledger is choice-independent:

    beta_M(x) = (1+x)^{A_M} + (1 + x^{M+1}) (1+x)^{d_M}   (mod 2).

Derivation: `2 D_M = 2 S_M - 1 + p^{M+1} + q^{M+1}` and `2 S_M` has even
integer coefficients; homogenize to degree `A_M`. Corollary: whenever
`binom(A_M, o)` is odd at a band position `o in [d_M+2, M-1]` (e.g. all o
when `A_M = 2^r - 1`), **every** scheme has a half-odd homogenized
coefficient in the untouchable band at time M. Amortization, if possible
at all, must come from sign oscillation, never from magnitude decay.

## 3. Exact carry-chain simulation (Task 2) — EVIDENCE

Greedy policy (each cell zeroes its own position within box+parity),
referee-checked at every step; independent value recomputation from the raw
choice list agrees exactly with the descended ledger. Runs (>= 30 shells;
M up to 150 rows):

| gamma | D0 | M   | verdict | frozen |D_M(1/2)| |
|-------|----|-----|---------|---------------------|
| 1     | 0  | 60  | amortizes (band never opens)   | 6.0e-36 and falling |
| 1/2   | 0  | 80  | FROZEN  | 3.33e-3 |
| 1/2   | 6  | 80  | FROZEN  | 2.05e-10 |
| 3/5   | 0  | 80  | FROZEN  | 2.97e-3 |
| 3/5   | 6  | 80  | FROZEN  | 2.17e-18 |
| 3/4   | 0  | 80  | FROZEN  | 1.07e-5 |
| 3/4   | 6  | 80  | FROZEN  | 3.48e-20 |
| 9/10  | 0  | 120 | FROZEN  | 8.54e-11 |
| 19/20 | 0  | 150 | FROZEN  | 1.44e-20 |
| 99/100| 0  | 60  | still decaying — band unborn (m* = 200): mimics gamma=1 |

Hostile lazy control (never absorbs): residual 0.155 — greedy does real work.

**Frozen-residual law (EVIDENCE):** the residual value freezes at the band
birth scale and satisfies `log2 |D_inf(1/2)| = -Theta(m*)`,
`m* = (D0+2)/(1-gamma)` (measured -8.2 vs predicted -(1+gamma)m* = -6 at
(1/2,0); -33.4 vs -38 at (9/10,0); -66 vs -78 at (19/20,0); D0=6 rows
overshoot the naive slope but keep the linear-in-m* form). Exact zero
requires m* = infinity, i.e. gamma -> 1 or D0 -> infinity. The **empirical
critical gamma is 1**: no gamma < 1 amortizes once M exceeds the band
birth; every gamma < 1 mimics gamma = 1 exactly until m*, then strands a
nonzero residual. This extends the single-jump slack lemma
(`amm12592_single_jump_routing_slack_deathstar.py`, gamma = 1/2) to the
full greedy multi-hop ledger and all tested gamma.

**Necessary envelope (PROVED, used as the failure criterion):** any feasible
scheme has `|D_M(p)| <= (p^{M+1}+q^{M+1})/2` for all p, M (boundary shares
`|G - 1/2| <= 1/2`). The frozen greedy trajectories violate it at explicit
small M (e.g. gamma = 3/5, D0 = 0: residual 3e-3 > 2^{-M-1} from M = 11 on),
so each greedy prefix is certified dead, not merely slow.

**Local sign-search control (EVIDENCE against sign-oscillation rescue):**
on (3/5, 0), flipping 29 small-capacity cell deficits (m <= 36) tuned at 5
biases reduced the tuned residuals by 1-2 orders but left 6 held-out biases
at (or above) the frozen scale (one worsened 8x). Consistent with the
mechanism: finitely many bumps can interpolate zeros, not kill the function.

**Mechanism verdict:** at the level of the carry-chain design space (shells
of ratio rho < 2 = depth slope gamma < 1, any transport policy realized as
per-cell absorptions; greedy + local-search probes), C = 1 + gamma < 2 with
constant D is REFUTED empirically with a structural explanation: the
untouchable middle band opens at m* and permanently strands value; the
per-row forced corner half-token stream (Sec. 1) keeps feeding it. This is
not a proof: a policy that keeps the band value identically zero forever
(anticipatory, cross-row-coordinated) is not excluded — but it must beat
the PROVED constraints: corner isolation, forced band parity (Sec. 2), and
the necessary envelope, simultaneously.

## 4. Realization / brute-force control (VERIFIED-EXACT)

For the gamma = 1 control (M = 4): decomposed the merged deltas back into
per-cell `w_{m,k}, v_{m,k}`, built the literal stopping rule, enumerated
every stopped word (all words of length T(m) = 2m+1 per class), and verified
the exact polynomial identity
`P_H(p) = 1/2 - (p^{M+1}+q^{M+1})/2 + D_M(p)` plus the pathwise deadline.
This closes the loop scheme <-> ledger exactly. (Task 3's length-24
brute-force for a gamma < 1 winner is moot: no winner emerged.)

## 5. Certificate (27) shape comparison (Task 1) — OPEN, with exact anchors

Re-verified exactly: `(2457/6592) L(t_B) - U(t_A) = 1/25 + F, F > 0`, hence
`(2457/6592) log(p_B/q_B) - log(p_A/q_A) > 1/25`.

- Shape match: our analysis shows single-bias tests can never obstruct
  (at p = 1/2 the row sums are unconstrained-even; the LP relaxation is
  feasible with all deltas 0), so **any true lower-bound certificate must
  couple >= 2 biases with the integrality/parity structure** — precisely
  the form of (27): two log-likelihood ratios (rapidities) with rational
  weights and a rational margin. On shape grounds (27) reads as the numeric
  gate of a two-bias dual/transport certificate (lower-bound step), not as
  an upper-construction gate.
- If gamma = 2457/6592 (C = 9049/6592 ~ 1.3727): band edges give binding
  biases {6592/9049 ~ 0.72848, 2457/9049 ~ 0.27152}; p_B = 0.74840 is 2.7%
  from the upper edge; p_A = 0.58918 matches neither. Alternative reading
  gamma' = q_B/p_B = 0.33619 (C' = 1.336) puts p_B exactly at the band
  edge by construction but leaves p_A unexplained.
- Near-miss flagged as probable numerology (repo has prior coincidence
  retractions): `-log p_A / log q_A = 0.594678` vs
  `2457/(6592-2457) = 2457/4135 = 0.594196` (0.08% apart) — the equal-mass
  ray slope at p_A vs a lattice slope from the certificate. Not asserted.
- 2457/6592 as gamma directly: no ledger quantity computed here selects it.
  Status of (27)'s role: OPEN.

## 5b. Odd-ladder addendum (coinC2 session, greedy control at 1/3 and 1/5)

THM-2976 (binary-clock parity) singles out the odd unit fractions
`gamma = 1/J` as the corner-clocked rates, `1/3` the top rung. Greedy runs
at `gamma = 1/3, 1/5` (D0 = 0, M = 120, same referee asserts): **frozen**,
residuals `|D_inf(1/2)| = 4.179e-2` and `7.295e-2` — *larger* than
`gamma = 1/2` (3.33e-3), consistent with the frozen-residual law
`-Theta(m*)` since band birth is earlier (`m* = 3` resp. `2.5`). Verdict:
corner clocking does **not** rescue greedy transport; the clock's live role
is extraction-side (lower bounds), not construction-side, unless an
anticipatory policy (lane F wave 2) exploits it in a way greedy cannot.
Driver: `04-computation/amm12592_oddladder_greedy_freeze_deathstar.py` over
this file's Ledger (verbatim import); envelope violations from `m = 10` on
at all three rates.

## 6. Next obligations

1. Attempt a theorem: "band value is asymptotically conserved" — quantify
   the positive-kernel drift (band-to-cone leakage is a binomial tail,
   exp(-c(theta) m) per row) and combine with corner isolation + forced
   band parity to prove no scheme with gamma < 1, D0 finite exists
   (would give C* = 2). The three PROVED ingredients are in Sec. 1-2.
2. Exact beam/branch-and-bound over deficits against the PROVED envelope
   `|D_M(p)| <= (p^{M+1}+q^{M+1})/2` at several biases simultaneously for
   gamma = 1/2, D0 = 0, small M: could yield the first unconditional
   finite-M infeasibility (lower bound C > 3/2 for that depth law).
3. Reverse-engineer (27) against a two-bias dual of the band-transport
   problem at gamma = 2457/6592 (weights (-1, 2457/6592) on rapidities);
   needs the dual formulated with parity, not the (trivially feasible) LP.
