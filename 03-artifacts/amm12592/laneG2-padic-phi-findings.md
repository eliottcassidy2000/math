# Lane G2 findings: the 2-adic coupled functional Phi — invariance verified, one-bit collapse, and the product-formula wall settled

**Session:** death-star-2026-07-30-coinC2, lane G2.
**Script (all asserted computations exact int/Fraction):**
`04-computation/amm12592_laneG2_padic_phi_deathstar.py` (exit 0, "ALL LANE-G2
CHECKS PASS"; runtime ~1 min).
**Frame:** THM-2966 spine normal form; laneD ledger conventions
(`laneD-carry-ledger-band-findings.md`). Certificate biases (context
CONFIRMED): `p_A = 1285/2181` (`b_A - a_A = 896 = 2^7*7`, `s_A=7, u_A=7`),
`p_B = 8847357/11821757` (`b_B - a_B = 2974400 = 2^6*46475`, `s_B=6`,
`u_B = 46475` odd).

**Conventions.** `N_X := 2 D_M(p_X) b_X^{A_M} = sum_cells delta a^z (b-a)^o
b^{A_M-z-o}` (integer; per-term `v_2 = v_2(delta) + s_X o`). Work with the
integer `Psi := 2 Phi = x N_A/2^{s_A} + y N_B/2^{s_B}`, `x = 1`, `y` = the
odd solution mod `2^7` of the unit congruence
`u_A a_A^{A-1} + y u_B a_B^{A-1} == 0 (mod 2^7)` (`A = A_M`). "Phi mod 64"
below = "Psi mod 128". Depth laws `d_m = floor(gamma m) + D0`.

## 1. (D1) VERIFIED-EXACT and PROVED: Phi mod 64 is choice-independent; K = 6 exactly

- **Invariance (VERIFIED-EXACT + PROVED).** For `(gamma, D0, M)` in
  `{(1/2,0,8), (1/2,0,11), (1/2,0,14), (3/5,0,10), (3/5,0,14), (1/2,4,7)}`,
  33 schemes each (all-w=0, all-w=box, canonical-parity, 30 random in the
  boxes with correct forced parity): `Psi mod 128` **identical across all
  schemes** in every case. Mechanism (PROVED, as in (D1)): `a == b mod 2^s`
  makes `a^z b^{A-1-z} == a^{A-1} mod 2^s` (z-independence), so every `o=1`
  shift changes `Psi` by `2[x u_A a_A^z b_A^{A-1-z} + y u_B a_B^z
  b_B^{A-1-z}] == 2[x u_A a_A^{A-1} + y u_B a_B^{A-1}] == 0 mod 2^7` by the
  unit congruence; `o>=2` shifts have `v_2 >= 1 + min(s_A, s_B) = 7`.
- **The maximal modulus is K = 6 exactly (PROVED + exhaustive check).**
  Per-cell shift spectrum at `(1/2, 0, 12)`: min `v_2(Psi-shift) = 7`
  attained BOTH at `o=1` cells in the odd `z`-parity class (because
  `v_2(b_B - a_B) = 6` exactly, so `a_B^z b_B^{A-1-z} mod 2^7` depends on
  the parity of `A-1-z`) AND at side-B `o=2` cells
  (`v_2 = 1 + 2 s_B - s_B = 7`, **independent of y** — so `K <= 6` for
  every odd `y`, a proof, not a search). Exhaustive scan over all odd
  `y mod 256`: max over `y` of min shift-`v_2` = 7. Odd unit rescaling of
  `(x,y)` preserves all valuations, so `x = 1` is WLOG. The prediction
  "maybe `K = 7` because side-A `o=1` shifts are z-independent mod `2^7`"
  is **REFUTED**: side B binds through two independent channels, both at
  the `2^6` (Phi-units) level. Random-scheme differences realize `v_2 = 7`
  (sharp) in every tested case.
- **Single-bias baseline (PROVED + VERIFIED-EXACT).** At one bias the
  choice-independent information is exactly `N_X mod 2^{s_X+1}` — **one
  parity bit** — and it equals the parity of `M+1`: homogenized, `c_0 = 0`
  (every cell has `o >= 1`) and `c_1 = sum of o=1 deltas == M+1 mod 2`
  (the `o=1` cells are exactly the `k=0` cells of rows `1..M` plus the row-1
  1-side corner, all with box 1). Hence `v_2(N_X) = s_X` iff `M` even
  (asserted for all schemes in T1, both biases). Corollary (cute, PROVED):
  `Psi` is always even, i.e. **Phi is always an integer** — the two odd
  leading bits cancel in the coupled sum.
- **What coupling buys.** On the integer `x N_A + 2 y N_B` the forced
  modulus is `2^{14}` versus `2^8` (bias A alone): the congruence-tuned
  coupling is a real 2-adic boost — but bounded (Sec. 4), and its content
  collapses to one bit (Sec. 2).

## 2. The forced value: closed form and ONE-BIT collapse (PROVED + VERIFIED-EXACT)

Since the all-`w=0` scheme is a valid box point and makes `S_M = 0`, the
invariant equals its boundary-term value, exactly:

```text
Psi_forced = sum_X wgt_X (a_X^{M+1} + (b_X-a_X)^{M+1} - b_X^{M+1})
                          * b_X^{d_M} / 2^{s_X}
```

(asserted `== Psi(all-w=0)` exactly, and `== Psi(random) mod 128`, for
M <= 14 over three depth laws). Expanding mod 128 with
`a_X^{M+1} - b_X^{M+1} = -2^{s_X} u_X sum_i a^i b^{M-i}` and the unit
congruence gives the **one-bit law** (M >= 2; verified for
`M in [2, 3000)` across six depth laws including `gamma = 2457/6592`):

```text
Psi_forced == 64 * bit(M),   bit(M) = (M+1) d_M + floor((M+1)/2)  (mod 2).
```

(`M = 1` special: the extra term `(b_B-a_B)^2 b_B^{d}/2^6 == 64` flips the
bit.) Tables (M = 1..40): `gamma=1/2`: nonzero exactly at `M == 1 mod 4`,
`M >= 5`; `gamma=3/5` and `gamma=2457/6592`: denser patterns, 19-20/40
levels nonzero; `D0 = 4` does not change the `1/2` pattern (`d_M mod 2`
unchanged).

**Interpretation (the decisive structural point).** The invariant depends
on the depth law **only through `d_M mod 2`** and on `M mod 4`. It is pure
boundary-telescoping bookkeeping — the 2-adic shadow of
`2D_M = 2S_M - 1 + p^{M+1} + q^{M+1}` — and is **completely blind to the
band**: two depth laws agreeing in `d_M mod 2` get identical invariants
regardless of interior/band parity structure. Answering the assigned
question: yes, the forced value is nonzero at levels where (D3)'s timing
puts deep forced-odd corners — e.g. at the certificate slope
`gamma = 2457/6592, D0 = 0`, both deep-corner clock levels `M = 2`
(`A=3, d=0`) and `M = 5` (`A=7, d=1`) carry `Psi_forced = 64` — but by the
one-bit law this is boundary arithmetic (`bit(2)=1` since `d_2=0`,
`bit(5) = 6*1+3 = 1`), **not** band information, and by Sec. 4 it cannot
be turned into an obstruction. The nonzero forced value is consistent with
(indeed forced upon) every scheme, feasible or not.

## 3. (D3) timing scan: the deep-corner clock is NOT tuned to the certificate slope

Scan `M <= 8192`, `gamma in {1/2, 3/5, 2457/6592}`, `D0 = 0..8` (T4):

- All-odd-band levels (`A_M = 2^r - 1`) occur for every law with geometric
  spacing (0-11 hits), EXCEPT `gamma = 3/5, D0 in {4,7}` where `A_M`
  provably skips all of them in range (parity of the `floor` steps) — a
  reminder that band-parity pressure is depth-law-sensitive.
- Deep-corner simultaneity (`A_M = 2^r - 1` AND `d_M = 2^rho - 1`,
  equivalently `M = 2^r - 2^rho - 1` with the floor-timing condition):
  fires only at **small sporadic M** everywhere: `(M,r,rho) = (1,2,1),
  (3,3,2), (7,4,3), (15,5,4)` for the appropriate `D0`, plus
  `(2,2,0), (5,3,1)` at the certificate slope `D0 = 0`.
- **PROVED (mechanism):** infinitely many deep-corner firings force
  `|2^rho - 1 - D0 - gamma(2^r - 2^rho - 1)| <= 1` for infinitely many
  `(r, rho)`, i.e. `2^r |(1+gamma) 2^{rho-r} - gamma| = O(1)`; since
  `2^{rho-r}` is a power of two this requires
  `gamma/(1+gamma) = 2^{-j}` exactly, i.e. **`gamma = 1/(2^j - 1)`**
  (j=1: the classical `gamma=1, C=2`; then `1/3, 1/7, ...`).
  `2457/6592` has `gamma/(1+gamma) = 2457/9049 ~ 0.2715 != 1/4`, so the
  clock fires only finitely often there. The (D3) speculation "the
  certificate rates are tuned to the deep-corner clock" is **REFUTED in
  its literal form**. (Side observation, SPECULATION: the natural clock
  slopes `1/(2^j - 1)` bracket the certificate slope `0.3727 in (1/3, 1)`;
  if the deep-corner clock drives the obstruction, the certified slope
  sitting strictly above the accumulation slope `1/3` reads naturally as
  rate slack, not as a clock resonance.)

## 4. (D2) The product-formula wall is REAL — every evaluation variant fails

The contradiction template needs a derived integer forced nonzero mod
`2^{K+1}` with archimedean envelope `< 2^K`. Using the PROVED necessary
envelope `|D_M(p)| <= (p^{M+1}+q^{M+1})/2`:

- **Single bias, single level (PROVED, one line).** The envelope of `N_M`
  is `E = (a^{M+1} + (b-a)^{M+1}) b^{d_M}`, and

  ```text
  E >= (b-a)^{M+1} b^{d_M} >= (b-a)^{A_M} >= 2^{s A_M} >= 2^{s(M+1)},
  ```

  while the available forced modulus is `2^{s+1}` (Sec. 1; shifts at `o=1`
  have `v_2 = s+1` exactly). So `E / 2^{forced} >= 2^{sM - 1}`: the wall
  holds for **every** bias `p = a/b` and every `M >= 1`, with margin
  exponential in M. The middle inequality also kills the "deep forced-odd
  cell" hope: `2^{s o} <= (b-a)^o <= E` for every `o <= A_M`, so **no
  per-cell forced valuation, however deep (including the 1-side cells with
  o up to A_M - 1), ever exceeds the envelope**. Verified exactly for
  `M <= 40` at both certificate biases: the log2 gap is >= 13 (bias A,
  M=1, the global minimum) and grows to ~1409 (bias B, M=40).
- **Two-bias coupling (PROVED for scalar weights).** Envelopes ADD under
  additive coupling; the forced modulus is capped at relative `2^6` =
  `min(s_A, s_B)` by the side-B `o=2` shifts, whose valuation no odd `y`
  can raise (Sec. 1). Coupling buys a bounded 2-adic boost and pays an
  added envelope: hopeless against a `b^{A_M}`-sized bound.
- **Row polynomials.** Feasibility forces nothing on a single row beyond
  the box (`0 <= W_m <= 1` is automatic); there is no decaying row
  envelope to race against. Dead on arrival.
- **Differences `D_{M2} - D_{M1}`.** `|D_{M2}-D_{M1}| <= |D_{M1}|+|D_{M2}|
  ~ max(p,q)^{M1+1}` and clearing denominators still costs `b^{A_{M2}}`:
  same chain, same wall.
- **Weighted sums over M (telescoping) — dollar-for-dollar (PROVED +
  VERIFIED-EXACT).** For `Psi_c = sum_M c_M N_M` (single bias), a shift in
  row m changes it by `2 a^z (b-a)^o b^{A_m-z-o} T_m` with
  `T_m = sum_{M>=m} c_M b^{A_M - A_m}`; invariance mod `2^{1+s+j}` forces
  `T_m == 0 mod 2^j` for all m, a triangular system whose last equation is
  `T_{M0} = c_{M0}`, so **`2^j | c_{M0}` exactly**: the envelope inflates
  by the same `2^j` that the forced modulus gains. T6 demo
  (bias B, `gamma=1/2`, `M0=6`, `j = 0,4,8,12`): forced `K = 7+j` on the
  nose, and `log2(E/2^K) = 225` **independent of j**. The `b^{A}` factor
  is untouched by any telescoping.
- **General multi-bias, multi-level functionals (PROOF SKETCH).** Any
  ledger functional `L = sum lambda_cell delta_cell` invariant mod `2^K`
  has `2 lambda_cell == 0 mod 2^K` at every cell (cells shift
  independently by 2 inside their boxes). If `L`'s boundedness on feasible
  schemes is derived from the proved envelope at finitely many
  `(M_i, p_j)`, then `lambda` is a combination of evaluation kernels
  `p_j^z q_j^o`, and the `o=1` cells of the deepest row pin
  `v_2(weights)` up to cross-bias cancellation; cancelling the
  z-DEPENDENT part mod `2^{s+j}` (a combination of J geometric sequences
  `(a_j/b_j)^z` with ratios `== 1 mod 2^{s_j}`) supports at most
  ~`J` congruence conditions, so the forced-modulus gain is bounded by
  `O(J max_j s_j)` while every participating bias contributes envelope
  `>= 2^{s_j(M_j+1)}`. Bounded gain, exponential cost: the wall stands for
  every fixed finite bias set as `M -> infinity`, and already at `M = 1`
  for the certificate pair (gap `>= 2^{13}`).
- **Scope note.** The wall applies to archimedean envelopes cleared at
  rational biases. A genuinely 2-adic envelope (`v_2(N_M) -> infinity`
  forced) would be an equality-type statement, not a size bound — nothing
  in the proved toolkit supplies one.

**Which (M, A, K) would be needed:** `K > log2 E ~ A_M log2 b - (M+1)
log2(1/max(p,q))`, and since `max(p,q) >= 1/2` and `b >= 3`, the needed K
grows at least like `s(M+1)` (indeed like `A_M log2 b`), while the
available K is `<= 1 + s + O(1)`. Needed grows linearly, available is
constant. **No variant beats the wall.**

## 5. Verdict for the certificate (27) — stated boldly

**The evaluation reading of (27) is DEAD.** No single-bias evaluation, no
congruence-tuned two-bias coupling, no row/difference/telescoping variant
can force a nonzero integer below its envelope: the forced 2-adic
information available from the ledger is a bounded number of bits (one bit
per bias, one extra coupled bit for this pair), and that bit is **pure
boundary bookkeeping** (Sec. 2), while every denominator-clearing costs
`b^{A_M}`. Therefore, given the PROVED wall, **certificate (27) must be
the numeric gate of a rate/entropy (tropical) dual** — a comparison of
exponential rates (transport capacity along the two critical rays vs
forced deficit mass) in which the two rapidities `log(p_X/q_X)` enter as
Legendre-dual slopes (HYP-9061 sec 2d) and `alpha = 2457/6592` as a ray
weight; the `1/25` is a rate margin per unit M. This is consistent with
laneD Sec. 5 (single-bias tests cannot obstruct; the LP without parity is
feasible) and with the no-near-cancellation anatomy (the linear form sits
above a floor — a rate gap, not a small linear form). What the 2-adic
engineering of the biases (`s_A = 7` vs `s_B = 6`) is actually for, under
this reading: it is the minimal misalignment that makes the coupled
one-bit invariant exist at all — a **checksum** aligning the two
evaluations' forced parity layers, in the spirit of THM-2225's dyadic
checksum shares, not a smallness mechanism (this last role assignment:
SPECULATION).

## 6. Status ledger

| claim | status | where |
|---|---|---|
| (D1) Phi mod 64 choice-independent | PROVED (mechanism) + VERIFIED-EXACT (33 schemes x 6 laws) | T1 |
| K = 6 exactly; side-B binds twice; y-scan cannot help | PROVED + exhaustive y mod 256 | T2 |
| single-bias forced info = one bit = parity of M+1; v2(N_X)=s_X iff M even; Phi always integer | PROVED + VERIFIED-EXACT | T1 |
| forced value = boundary closed form; ONE-BIT law `64[(M+1)d_M + floor((M+1)/2)]` | PROVED (M>=2) + VERIFIED-EXACT (M<3000, 6 laws) | T3 |
| invariant blind to band (depends only on `d_M mod 2`, `M mod 4`) | PROVED (corollary of one-bit law) | T3 |
| deep-corner clock infinite iff `gamma = 1/(2^j - 1)`; cert slope excluded | PROVED (mechanism) + scan M<=8192 | T4 |
| wall, single bias: `E >= 2^{s A_M}`, gap `>= 2^{sM-1}` | PROVED + VERIFIED-EXACT (M<=40) | T5 |
| wall, weighted sums: dollar-for-dollar | PROVED + VERIFIED-EXACT (j=0..12) | T6 |
| wall, general multi-bias functionals | PROOF SKETCH (finite-difference bookkeeping open) | Sec. 4 |
| (27) is a tropical/entropy dual, not an evaluation bound | EVIDENCE (forced by the wall within the evaluation class) | Sec. 5 |
| 2-adic engineering = checksum alignment | SPECULATION | Sec. 5 |

**Next obligations (for other lanes):** (i) lane G's reconstruction should
now target ONLY the entropy/max-flow dual: weights `(-1, alpha)` on the
two rapidities as a two-ray capacity-vs-mass rate inequality, margin
`1/25` absorbing `O(D0)`; (ii) the one-bit invariant is available as a
cheap consistency checksum for any candidate construction at
`gamma = 2457/6592` (its value is forced; a construction violating it at
any M is wrong); (iii) the finite-difference bookkeeping in the general
wall sketch is a clean small theorem if anyone wants the wall airtight
beyond the certificate pair.
