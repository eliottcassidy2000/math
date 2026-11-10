# Pre-registration - reach of the sign-flip towers (kind-pasteur-2026-11-10-S22)

Written and pushed BEFORE the runs E1 and E2 below. Status of everything in this file: PRE-REGISTERED PREDICTIONS,
not results. Owner directive of this session: construct a seed for every odd `c` (HYP-9658 (A)).

## What was already seen before this file (not predictions)

- Notation of THM-5126: a pair `(B, O)` has offset `r = H(B) - W_B(O)` and slope `M = U_B(O)`; the closed tournament
  `B => TT_t => v`, `v -> O`, has `H = r + 2^t M` (KNOWN cut form). A *seed* of the line `1 + c 2^s` is a pair with
  `r = 1`, `M = c 2^j`; an *anti-seed* (offset `-1`) gives the ray `c 2^s - 1`.
- Working claim (derived by hand, checked by brute force on 400 random base / gadget pairs, 0 failures, script
  `towers.py`): for a pair `(B, O)` with state `s = (r, M, h, E)` (`h = H(B)`, `E` = Hamiltonian paths of `B` ending
  in `O`) and a *gadget* `(X, P)` put at the tail (`B+ = B => X => v`, `v -> O`, new out-set `P` inside `X`), the new
  state is `G(X, P) s` for an explicit 4x4 integer matrix. For `X = TT_a` and `P` two vertices at distance `delta`:
  `r+ = (1 - 2^(delta-1)) r + 2^(a-delta) (2^delta - 1 - 3^(delta-1)) M`, so the offset is exactly reversed
  (`r+ = -r` for every base) iff `delta` is 1 or 2, the two solutions of `2^delta - 3^(delta-1) = 1`.
- Exploratory scan already run (`neutral.py`): among all gadgets `(X, P)` with `X` of order <= 7 the coefficient of
  `M` vanishes only for `|P| = 2`, and then `r+ = -H(X) r`; with `H(X) = 1` these are exactly the transitive pairs at
  distance 1 or 2.
- Two-module families of kernels of order <= 8 (`mfan.py`): the offset `-1` families seen there have odd parts
  3, 5, 9, 13, 15, 17, 27, 35, 39, 41 and several affine forms; no family with odd part 1.
- Hand computation of levels 1 and 2 from transitive bases (families `3 2^j + 5 3^m`, `15 3^m`, `30 3^m - 1`, ...).
  No coverage count of any kind has been computed.

## The state recursion that the runs iterate

With `g = h + r`, `alpha = 2^delta + 1`, `gamma = 3^delta + 2^delta` (`delta` in {1, 2}), `x >= 0` (number of sinks
added to the base plus the number of tail vertices below the pair), `m >= 0` (tail vertices above the pair):

    r+ = -r,   g+ = 2^(x+1+m+delta) M,   M+ = alpha 2^m g + 2 gamma 2^x 3^m M - [x = 0] E,
    E+ = alpha 2^x M + [m = 0] (2 r - g).

## E1 - pure towers

Start states: `(r, g, M, E) = (1, 2, 2^k, [k = 0])` for `0 <= k <= 10` (transitive base, one out-neighbour) and
`(1, 2, 0, 0)` (transitive base, empty out-set). Box: `delta` in {1, 2}, `0 <= x <= 10`, `0 <= m <= 4`, depth <= 4.
A state of even depth has `r = +1`; `c` = odd part of `M`. `N_d` = number of odd `c <= 4097` occurring at even
depth <= d (the value `c = 1` of depth 0 included).

- **R1.** `300 <= N_4 <= 1400` (of the 2,049 odd numbers up to 4097).
- **R2.** `N_4 >= 1.5 N_2`.
- **R3.** No state of odd depth <= 3 has `M` a power of two (no ray `2^s - 1` from towers).
- **R5.** The least odd `c >= 3` that is not counted in `N_4` is at most 13.

## E2 - towers fed by the census

Start states: all pairs `(B, O)` with `B` a tournament of order <= 7 (one representative per isomorphism class, every
vertex set `O`) and `r = +1` or `-1`. Same box, at most two levels. Seeds = states with `r = +1`.

- **R4.** At least 3 of the 33 odd `c <= 4097` listed in HYP-9658 (B) (no seed among unit-offset roots of order <= 10
  and doubling vertices of order <= 9) acquire a seed.

## Decision rule (G3)

- `N_4 >= 1900`: the towers are close to complete on this range; the next step is a completeness proof for the
  recursion (a covering argument for the odd parts).
- `300 <= N_4 < 1900`: the towers are one mechanism among several; report the exact families, the uncovered `c`, and
  what the uncovered seeds of the census look like.
- `N_4 < 300`: the towers are thin; the write-up is restricted to the structure theorem and the explicit families.

## E3 - depth six (added and pushed before the run; E1 and E2 were already run when this was written)

Known when this was written: E1 gave `N_2 = 805`, `N_4 = 1612`; the odd `c <= 129` without a pure-tower seed of depth
<= 4 in the E1 box are 5, 11, 23, 29, 73, 85, 91, 95, 127 (29 is reached from the transitive pair start
`(r, g, M, E) = (-1, 0, 3, 1)`, which E1 did not include). Nothing of depth six has been computed.

Run: pure towers from the E1 start states, box `delta` in {1, 2}, `0 <= x <= 5`, `0 <= m <= 2`, depth <= 6.
`N_6'` = number of odd `c <= 4097` at even depth <= 6 in this smaller box.

- **R6a.** At least 4 of the eight numbers 5, 11, 23, 73, 85, 91, 95, 127 occur at even depth <= 6.
- **R6b.** 5 does not occur.
- **R6c.** `N_6' >= 1500`.
