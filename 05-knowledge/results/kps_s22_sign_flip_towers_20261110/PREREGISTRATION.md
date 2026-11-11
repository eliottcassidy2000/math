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

## Outcomes of E1 and E2 (written after the runs; nothing above this line was changed)

Runs: `runs/e1_d2.err`, `runs/e1_d4.err`, `runs/e1_d4.txt` (E1); `runs/e2_states_order_le7.txt`, `runs/e2_d2.err`,
`runs/e2_d2.txt` (E2).

- **E1.** `N_2 = 805`, `N_4 = 1612`; 1,330 odd `c <= 4097` occur with offset `-1` at depth <= 3.
  - R1 (`300 <= N_4 <= 1400`): **FAILED**, the towers reach more than predicted.
  - R2 (`N_4 >= 1.5 N_2`): **CONFIRMED** (1,612 against 1,207.5).
  - R3 (no power of two as slope at odd depth <= 3): **FAILED**. 15 states of depth 3 have offset `-1` and slope
    `2^23`, `2^24` or `2^25`; all 15 were rebuilt as tournaments and recounted by the transfer-matrix count
    (`scripts/mersenne_towers.out`): `2^N - 1` paths on `N >= 24` vertices. (That `2^k - 1` is a path count for
    `k >= 4` is KNOWN, THM-4748 (4); what the prediction missed is that a tower can have a power of two as slope.)
  - R5 (least uncovered `c >= 3` is at most 13): **CONFIRMED**, it is 5. Uncovered up to 129: 5, 11, 23, 29, 73,
    85, 91, 95, 127.
- **E2.** Seeds for 1,872 and anti-seeds for 1,782 of the 2,049 odd `c <= 4097`.
  - R4 (at least 3 of the 33 exceptional `c`): **CONFIRMED**, 32 of 33; the missing one is 3773.
- **G3.** `300 <= N_4 < 1900`: the towers are reported as one mechanism among several.

Not pre-registered (done after E1 and E2, reported as exploratory in the note): the transitive pair starts
`(-1, 0, 3, 1)` and `(-1, 0, 5, 1)` (`runs/e1x_d3.txt`); starts of order 8 (`runs/e2o8_d2.txt`) and three levels
(`runs/e2_d3.txt`); pairs of order 9 for the certificates. With them every odd `c <= 4097` has a seed and an
anti-seed (`certs_plus.txt`, `certs_minus.txt`, each recounted).

## Outcome of E3 (written after the run)

Run: `runs/e3_d6.err`, `runs/e3_d6.txt` (about two and a half hours of one core). `N_2' = 614`, `N_4' = 1357`,
`N_6' = 1592` odd `c <= 4097` with offset `+1`; 1,554 with offset `-1` at depth <= 5.

- R6a (at least 4 of the eight numbers): **FAILED**. Three occur at depth 6: 5 (slope `5 * 2^32`), 73 (`73 * 2^31`),
  85 (`85 * 2^30`). 11, 23, 91, 95, 127 do not occur.
- R6b (5 does not occur): **FAILED**. So the line `1 + 5 * 2^s` has a pure-tower seed with seven closers.
- R6c (`N_6' >= 1500`): **CONFIRMED** (1,592).

Reading: both failures go in the same direction as R1 and R3 (the towers reach more than expected), and depth six
adds 235 values of `c` to the 1,357 of depth four in this box; the reach is still growing with depth.

## E4 - sources for the 9-grid from the gadget law (added 2026-11-11, pushed before the run)

Context (owner's continuation: reusable source interfaces and the domination-defect inequality of Sessions 260-261).
Session261 realises `N(9, r, s) = 1 + 2^s + 9 * 2^(r+s)` for all `r >= 1` and `s <= 6` from singleton-marked sources
`(D, U)` = (offset, slope) with `D = 1 + 2^s` and odd part of `U` equal to 9, and leaves `s >= 7` open.

Derived by hand before any run (THM-5130 (A) with a one-vertex out-set in a transitive tail): from a seed
`(r, M, h, E) = (1, 2^k, h, E)` the gadget `(TT_a, {p})` gives the source

    D = 1 + 2^(a+k),      U = 2^(a-p) (1 + h) + 2^(p+k) 3^(a-p) - [p = 1] E.

With `p = a`: `U = 2^s + h + 1` (`s = a + k`), so the odd part of `U` is 9 exactly when `h + 1 = 2^(s-3)`, or
`5 * 2^(s-2)`, or `7 * 2^(s-1)`; since `M >= 2 H(B)` only the first is possible, with `k = s - 2` or `s - 1`.
Nothing below has been computed: no state `E` of the pairs with `(h, M) = (15, 2^k)`, no census of singleton sources
of order 8 or 9, no pair with `H(B) = 31` and a power of two as slope.

- **G1.** Some pair of order at most 8 with `(r, h, M) = (1, 15, 64)` and one gadget `(TT_a, {p})` give a source with
  `D = 129` and odd part of `U` equal to 9 (level `s = 7`), confirmed by a literal recount.
- **G2.** Among pairs `(B, {o})` with `B` of order at most 9 there are sources of the 9-grid for `s = 7` and for
  `s = 8` (`H(B) = 129`, `257` and odd part of `U_B(o)` equal to 9).
- **G3.** No pair of order at most 9 has `r = 1`, `H(B) = 31` and slope 64 or 128 (so level 8 does not come from
  this mechanism at these orders).
- **G4 (domination ray).** For every pair `(B, O)` with `B` of order at most 6, `M >= sum over nonempty S in O of
  2^(number of out-neighbours of S outside O)`, with equality for transitive `B` and one-vertex `O`. (If my
  computation of the dominating-set count of the closed tournament is right this is a consequence of the
  inequality of Session261; the run checks the computation.)

## Outcome of E4 (written after the runs)

Scripts and data: `scripts/sources_check.py` (+ `.out`), `scripts/singles.c`, `scripts/src10.c`,
`runs/grid_levels_singletons.txt`, `runs/src9grid.txt`, `runs/src43grid.txt`.

- G1: **FAILED**. The pairs of order at most 9 with `(r, h, M) = (1, 15, 64)` have `E = 6`, so the gadget
  `(TT_1, {1})` gives slope `144 - 6 = 138`; no gadget `(TT_a, {p})`, `a <= 6`, on a seed of slope `2^k` gives a
  9-grid source at level 7. The gadget law does produce the sources of levels 2, 4 and 5.
- G2: **CONFIRMED**. One-vertex sources of order at most 9 exist at levels 7 and 8 (and at 9 and 11; not at 10,
  where a two-vertex selected set of order 9 works).
- G3: **CONFIRMED**. No pair of order at most 9 has `r = 1`, `H(B) = 31` and a power of two as slope; the base
  counts of seeds with a power of two as slope are 1, 5, 9, 15, 19, 27, 51, 99, 351, 709, 1297, 1467, 1779, 1903.
- G4: **CONFIRMED**. `M >= delta_dom` for all 3,978 pairs of order at most 6 (and 32,767 labelled pairs of order at
  most 5), with equality exactly for transitive `B`; the ray form of the dominating-set count and the inequality of
  Session261 hold on every tested case.

Not pre-registered (found while running E4, reported as exploratory in SOURCES_AND_DOMINATION.md): a port with slope
family `16 + 172 * 2^a` (the 43-grid), and the sources that complete that grid.

## Audit notes and errata (2026-11-11, after audit/AUDIT.md; nothing above this line was changed)

- (m9) In "What was already seen": the coefficient of `r` is `-1` for `delta = 1` and `1 - 2^(delta-1)` for
  `delta >= 2`; the displayed formula is right only for `delta >= 2`. The conclusion (reversal iff `delta` is 1 or 2)
  is unaffected.
- (m8) Outcome of R3: "`2^N - 1` paths on `N >= 24` vertices" holds for the states with slope `2^23` and `2^24`; the
  three states with slope `2^25` give `2^(N+1) - 1` paths on `N >= 25` vertices.
- (M3) E4: the sentence "Nothing below has been computed" was wrong for G1 and G3. The file
  `runs/e2_states_order_le7.txt`, committed with that pre-registration, contains the state `1 16 64 6`
  (`(r, g, M, E)` of the pairs with `H(B) = 15` and slope 64), and the unit states of orders 8 and 9 already on disk
  show no `H(B) = 31` among the seeds whose slope is a power of two. G1 and G3 therefore do not count as blind
  predictions. G2 and G4 were blind. (MISTAKE-759.)
- (m25) The label G3 is used twice: for the decision rule of E1/E2 and for a prediction of E4. The message of commit
  b5912dbd67 says "outcomes of E3 and E4"; the outcome of E3 was already in the previous commit.
- (m10, m11) `scripts/domination_ray.py` (added after the audit) reproduces the 3,978 pairs of G4; `scripts/src10.c`
  now prints offsets from 3 on, and `runs/src9grid.txt`, `runs/src43grid.txt` were regenerated with it.
