# AMM 12592 / THM-3029 toolchain repair: broken imports from a deleted session worktree

**Session:** boxeph, 2026-08-03.
**Scope:** reproduction-chain repair only; no mathematical claim of THM-3029 is
changed. All re-verification below is exact (int / Fraction arithmetic; the
only floats are in the gamma_c *table*, unchanged from the committed
`amm12592_finite_R_capacity_threshold.py`).

## The hazard

`04-computation/amm12592_floor_rate_attained_thm3029.py` — the referee of
canon theorem THM-3029 — began with

```python
sys.path.insert(0, '/tmp/math-wt-coinC2/04-computation')
from liftrate import prof, lift_block, admissible, epoch_identity, eff
from gammac import gamma_c, GSTAR
```

`/tmp/math-wt-coinC2` was the death-star session worktree of 2026-08-01. It
has been deleted, and `liftrate.py` / `gammac.py` were **never committed**.
Consequence: the committed referee of a PROVED+VERIFIED-EXACT canon theorem
could not be run at all (`ModuleNotFoundError`), i.e. the reproduction chain
of THM-3029 was broken at the toolchain level. (The committed
`amm12592_profile_monotonicity_thm3029.py` had the same dead `sys.path` line,
though it only needs the beam module, which IS committed.)

This is exactly the failure mode the repo contract warns about: referee
scripts must be self-contained in the repository, never dependent on session
scratch directories.

## The repair

1. **`04-computation/liftrate.py` (reconstructed).** Definitions taken
   verbatim from the committed `amm12592_profile_monotonicity_thm3029.py`:
   - `prof(R, g1, g2, D0)[i] = (g1*(R+i))//g2 + D0`, `i = 0..R-1` (floor
     profile with slack `D0`);
   - `admissible(delta, d)`: Lucas-box conditions in the Bernstein basis
     `B_{d,k}(x) = x^{d-k}(1-x)^k`: `|delta_k| <= binom(d,k)` and
     `delta_k == binom(d,k) (mod 2)`;
   - `lift_block(delta, d, d')`: convolution with `[binom(d'-d,k)]_k`, the
     admissible block representing the constant 1 (THM-3026 (L));
   - `epoch_identity(R, sol, d)`: the THM-3002 normal-form identity (*)
     `sum_{i<R} x^i Delta_i(x) == (1-x)^{R-1}` in `Z[x]`, with
     `Delta_i = sum_k delta_{i,k} B_{d_i,k}`;
   - `eff(R, d) = max_i d_i/(R+i)` as an exact `Fraction` (guards the
     THM-3029 sec. 2 slack trap).
   The module keeps the original's module-level demonstration (the
   gamma = 1/2, D0 = 3 solve at R = 32 and its lifts), which prints lines
   1–7 of the committed expected output at import time.
2. **`04-computation/gammac.py` (reconstructed).** `GSTAR = 2 log(phi)/log 5`
   and `arch_margin` / `gamma_c` verbatim from the committed
   `amm12592_finite_R_capacity_threshold.py`, with one historical difference
   recovered from the expected output: the default bisection bracket is
   `lo=0.5`, so the module-level demo table clamps `gamma_c` at 1/2 for
   R = 8, 16, 32 (lines 12–20 of the .out show 0.5000000000 there), while the
   THM-3029 referee passes `lo=0.02, hi=0.75` explicitly and reproduces the
   unclamped 0.375 / 0.4412 / 0.5 values of
   `amm12592_finite_R_capacity_threshold.out`. Both demo prints were kept at
   module level because the committed expected .out contains them as
   import-time side effects (lines 1–25); that is also documented in each
   module's docstring.
3. **`amm12592_floor_rate_attained_thm3029.py`**: the single `sys.path` line
   now points at the script's own directory
   (`os.path.dirname(os.path.abspath(__file__))`), so the referee is
   worktree-relocatable. Nothing else was changed.
4. **Same one-line repair** applied to the three sibling THM-3029 scripts that
   carried the identical dead path: `amm12592_profile_monotonicity_thm3029.py`
   (needs only the committed beam module),
   `amm12592_effective_rate_thresholds_thm3029.py` (needs `gammac`, hence was
   equally unrunnable until now), and `amm12592_sub35_direct_beam_thm3029.py`.
   Not rerun here (long beam sweeps; their committed .out files are search
   logs, not exactness certificates).

   Out of scope but recorded: several `sk_*.py` scripts (S_k PSLQ lane,
   THM-3012 family) also reference `/tmp/math-wt-coinC2` **data/output paths**
   (`Sk_520.json`, result files) — a separate instance of the same hazard in a
   different lane, left untouched.

## Reproduction status

Rerun of the repaired referee vs the committed
`05-knowledge/results/amm12592_floor_rate_attained_thm3029.out`:

- **R1** (gamma_c table, R = 8..1024): reproduced digit-for-digit.
- **R2** (profile coincidence + eff rates 0.583333 / 0.592593 / 0.596774 /
  0.597938): reproduced.
- **R3** (R = 32 lift: source SOLVED+verified, pointwise domination,
  admissibility, exact epoch identity): reproduced, all True.
- **R4** (gamma* floor profile closes at R = 8, 16 direct and R = 32 by
  lifting; C = 1 + gamma* = 1.5979874356654402 attained for n <= 63):
  reproduced.
- Header lines 1–25 (import-time demo output of the two modules): reproduced.

DIFF RESULT: **BYTE-IDENTICAL** (`diff` empty, exit 0, stderr empty) between
the rerun of the repaired referee and the committed 55-line
`amm12592_floor_rate_attained_thm3029.out`. In particular the lo=0.5
default-bracket reading of the lost `gammac` module is confirmed exactly by
the clamped header table, and the reconstruction is validated end to end.

## Witness extraction

`04-computation/amm12592_floor_witnesses_extract.py` regenerates the three
concrete witnesses and verifies each EXACTLY (every block admissible at its
profile degree; epoch identity checked as an integer polynomial identity)
before writing

`04-computation/amm12592_floor_witnesses_R8_R16_R32.json`

a JSON array of `{"R", "profile", "blocks", "verified"}` with:

- `R = 8`, profile `[4,5,5,6,7,7,8,8]` (direct beam solve, eff rate
  7/12 = 0.583333);
- `R = 16`, profile `[9,10,10,11,11,12,13,13,14,14,15,16,16,17,17,18]`
  (direct, eff rate 16/27 = 0.592593);
- `R = 32`, profile `[19,19,20,...,37]` = `floor(gamma*(32+i))` (lift of the
  gamma = 1/2, D0 = 3 solution, eff rate 37/62 = 0.596774).

Log: `05-knowledge/results/amm12592_floor_witnesses_R8_R16_R32.out`.

## Standing hazards re-affirmed

- Beam-search negatives remain SEARCH ARTEFACTS (THM-3029 sec. 1); nothing
  here bears on infeasibility questions.
- Session-worktree paths in committed scripts are a repo-level MISTAKE
  pattern: any `sys.path.insert` pointing outside the repository should be
  treated as a broken reproduction chain and repaired as above.
