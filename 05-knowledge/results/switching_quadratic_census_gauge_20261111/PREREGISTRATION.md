# Pre-registration (mac-mini-2026-11-11-gauge)

Written and pushed before the runs of block R1 and before any attempt to score block R2.
Statements typed PROVED in the note (the endpoint-graph theorem, the quadratic law of `H mod 4` on a switching
class) are not registered: their checks (`qf_exhaustive.c`, exhaustive on every labelled tournament of orders 3-7;
random checks to order 10) are controls of a proof, and were run before this file was written.

## R1. The 2-adic derivative law on a switching class (out of sample)

Notation. `T` a tournament on `V`, `T^S` its switching at `S`. For `U` a set of vertices,
`c_U(F) = sum_{S subset U} (-1)^(|U|-|S|) F(T^S)`.

Predictions (conjecture HYP-9669 of this session):

- (P1) `v_2(c_U(H)) >= ceil(|U|/2)` for every non-empty `U` (`H` = number of Hamiltonian paths);
- (P2) `v_2(c_U(hc)) >= ceil(|U|/2) - 1` for every non-empty `U` (`hc` = number of Hamiltonian cycles);
- (P3) the binary digit `k` of `H` has `F_2`-degree at most `2^k` on the switching cube, for `k = 1, 2, 3`
  (`k = 1` is the proved theorem; `k = 2, 3` follow from (P1) by the digit lemma of the note).

What is already known when this is written (in sample): (P1), (P2) hold for every `U` in the random tournaments of
orders 6-10 examined (`val_family.py`, `digit_degree.py` in the scratch runs; 2 to 8 tournaments per order), and
for one random tournament of order 8, one random, one circulant and the transitive tournament of order 9
(`two_adic_cube.c`). `|U| <= 4` is proved. `c_U(H) = 0` for `|U| > 2 floor(n/3)` is KNOWN (THM-4524) and is a
control, not a prediction.

Runs registered (program `two_adic_cube.c`, arguments `n seed mode`; mode 0 = random tournament from the seed,
1 = transitive, 2 = circulant with a connection set drawn from the seed):

    11 1 0     11 2 0     11 3 0     11 7 2
    12 1 0     12 2 0     12 0 1
    13 1 0     13 4 2
    14 1 0

How the test objects were chosen: random tournaments (the least symmetric objects available) at every order, and
the two most symmetric kinds (transitive, circulant) once each at the largest affordable orders; no selection
after looking. Orders 11-14 are outside every run made so far.

Falsifiers: one `U` with `v_2(c_U(H)) < ceil(|U|/2)` or `v_2(c_U(hc)) < ceil(|U|/2) - 1`; a digit degree above
`2^k`. Any falsifier refutes HYP-9669 as stated (the proved case `|U| <= 4` cannot fail; if it does, the
program is wrong).

## R2. Forecasts for order 12 (cannot be scored today; recorded so that they can be)

Data used (THM-4889 / HYP-9481, klein, 2026-11-06; A038375): least unattained odd value other than 7 and 21,
`L(n) = 611, 2883, 14685, 80453` for `n = 8..11`; number of distinct values `|S(n)| = 320, 1520, 7663, 43911`;
maximum `661, 3357, 15745, 95095`; mean `n!/2^(n-1)`.

- (F1) `L(12)`, the least odd number other than 7 and 21 that is not `H` of a 12-vertex tournament:
  point forecast 488,000; interval [476,000, 500,000].
  Method: `L(n) 2^(n-1)/n! = 1.940, 2.034, 2.072, 2.064` for `n = 8..11`; forecast ratio `2.085 +- 0.05`.
- (F2) `max H(12)`: point forecast 545,000; interval [531,205, 552,000] (the lower end is the proved lower bound
  THM-329). Method: for even `n`, `e - max H(n) 2^(n-1)/n! = 0.718, 0.620, 0.497` at `n = 6, 8, 10`, successive
  ratios `0.864, 0.801`; forecast next ratio in `[0.72, 0.82]`. (The unsourced value 545,835 flagged in
  MISTAKE-743 lies in the interval; the forecast is not evidence for it and does not use it.)
- (F3) `|S(12)|`, the number of distinct values of `H` on 12-vertex tournaments: point forecast 262,000; interval
  [250,000, 272,000]. Method: `|S(n)| / (max H(n)/2) = 0.968, 0.906, 0.973, 0.924` for `n = 8..11` (even orders
  0.97), applied to (F2).
- Typing: NUMERICAL extrapolations from four or five data points; no model of the upper tail is claimed.

## R3. Exploratory runs (no prediction registered)

Corank statistics of the Hessian `M` of `H mod 4` on random switching classes at orders 9, 10, 11
(`rank_stats.c`). In sample (orders 6, 7, 8, seeds shown in `out/`): the uniform alternating-matrix law that
describes the vertex-cube form of THM-5127 does NOT describe `M` (order 6: corank 1 always; order 7: corank 0 in
62 percent; order 8: `M = 0` in 1.7 percent), so no random-matrix prediction is registered.

## Outcomes

(filled after the runs)
