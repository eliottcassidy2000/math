# AMM 12592 — ANGLE 5: hybrid endgame algebra (S_L interval decode) — boxeph

Session 2026-08-03. Code: `04-computation/amm12592_r128_endgame_algebra_boxeph.py`
(+ `.out` in this directory). Companion engine: `amm12592_r64_floor_solver_boxeph.py`.

## What this delivers

An EXACT characterization of the residuals the last L rows of an epoch can
absorb, wired into the winning l1deg beam as (i) a complete 2-row endgame that
replaces the coarse target-grid completion, and (ii) an optional
distance-to-absorbable beam score.

Notation: closure at profile (d_i) is sigma_{-1} = q^{R-1},
p sigma_i = sigma_{i-1} - Delta_i, each Delta_i admissible (box+parity) at d_i.
S_L := set of sigma absorbable by exactly the last L rows.

## The algebra (all verified exactly in --selftest)

1. **Dual formula.** For deg s <= d the Bernstein-d coordinates are
   `b_m = sum_j s_j C(d-j, m)` (from x^j = sum_m C(d-j,m) B_{d,m}); the decode
   is forced, so S_1 is decidable in O(d^2). (base.final_decode = same thing
   computed triangularly.)

2. **Band reduction (L=2).** sigma = A + xB at degrees (da, db), e := db-da.
   Fold the forced top cell a_da = sigma(0) = +-1; let beta = dual transform of
   (sigma - sigma(0)(1-x)^da)/x at db. Then for the unknowns a_0..a_{da-1}:

       b_m = beta_m - sum_k C(e+1, m-k) a_k ,   m = 0..db,

   i.e. the decode column of B_{da,k}/x collapses by Vandermonde
   (sum_u (-1)^u C(k,u) C(N-u,m) = C(N-k,m-k), N = e+1+k) to an (e+2)-band.
   2-row absorbability = staircase integer program of bandwidth e+1 with box
   `|a_k| <= C(da,k)`, `|b_m| <= C(db,m)` and ALL parities forced.

3. **e = 0: complete O(d) decision + construction** after the transform.
   Chain DP over parity-normalized intervals
   `U_m = A_m ∩ (beta_m - B_m - U_{m-1})`. Completeness: parities are all
   forced, so progression-endpoint normalization makes interval arithmetic
   exact; the Minkowski sum of two step-2 progressions is the full step-2
   progression of endpoint sums; constraint m touches only (a_{m-1}, a_m), so
   the forward intervals are exactly the reachable projections and backward
   choice a_{m-1} in U_{m-1} ∩ (beta_m - B_m - a_m) never fails. Verified
   against brute-force enumeration at (da,db)=(5,5): 60/60 agreement
   (30 feasible / 30 infeasible), every feasible instance constructed.

4. **e = 1: sandwich.**
   - *Necessity*: the same DP with window intervals (U_{m-1}, U_{m-2}) is a
     relaxation (tracked intervals contain the true reachable sets), so its
     first gap or a parity-precheck failure exactly refutes membership.
   - *Sufficiency*: (a) downgrade both rows to (da,da), solve by the complete
     chain DP, lift the second block da->db (multiplicative lift, preserves
     polynomial + admissibility); (b) budgeted DFS over the relaxed intervals
     with exact look-back constraints; (c) steered-grid + exact decode.
     Every success is verified independently (poly identity + box + parity).
   - Selftest at (5,6): all 30 brute-feasible instances constructed
     (29 dp-dfs, 1 downgrade-lift), 0 misses, necessity never violated.
   - Upgrade path to full e=1 completeness: the exact reachable sets of
     (a_{m-1}, a_m) are lattice points of a convex polygon on the fixed parity
     sublattice (the discarded coordinate enters only through a 1-D
     fiber-interval nonemptiness condition), so a piecewise-linear polygon DP
     closes the gap; not needed for the hunt.

5. **L = 3**: semi-decision both ways — necessity via the L=2 relaxation
   after one forced step; sufficiency via exhaustive steered grid on the first
   row followed by the L=2 solver (completion stage C2 below).

## Validation on real data (R = 64)

- sigma_61 = Delta_62 + x Delta_63 of the committed R=64 floor witness
  (degrees (75,75), e=0): absorb2_solve reconstructs an admissible pair
  (dp-e0); splicing it back verifies the full epoch identity.
- Full pipeline (l1deg beam 400 + S_2 completion) closes R=64 from scratch:
  `SOLVED-C1(dp-e0)`, exact verify passes.

## R = 128 (profile d_i = floor(gamma*(128+i)), 76..152, 76 unit gaps)

STATUS: no witness from this angle in 9 hunts (A..I; beams 400-1200, ranks
l1deg / a2-deficit, rand_frac, E-prune, E-branch, corridors); every endgame
refutation was an EXACT certificate (parity precheck or first interval gap of
the S_2 relaxation) and the 'unknown' (relaxed-pass, unconstructed) count was
0 across ~200k C2 leaves -- the sandwich never even entered its ambiguous
zone.  These are search negatives, not infeasibility evidence (THM-3029);
R=128 is independently CLOSED by three concurrent-angle witnesses
(R128_direct / R128_lp / R128_rule), all re-verified exactly by this
session's referee.  Final hunt (I: a2-from-64 + E-branch, beam 1200) was
captured mid-search by a fat low-deficit family (L1 ~ 1.9e72): the a2
deficit's degree-overflow term misranks states while deg >> db+1; the fix
(future work) is to score a2 only once deg <= db+3, i.e. from ~row 118, with
E-branching enabled from ~row 100.

### Full-scale positive controls

Synthesized true pairs at the actual endgame degrees: 12/12 reconstructed at
(151,152) (all via dp-dfs), 8/8 at (151,151) (complete regime). The machinery
is sound at scale; refutations below are genuine certificates.

### Hunt ledger (beam rows 0..125, then C1 = S_2 on all terminal states,
### C2 = wide grid on row 125 then S_2; all endgame refutations exact)

- HUNT-A l1deg beam 400 (exact R=64 winning recipe): exhausted; C1 400/400
  refuted, C2 15,680 leaves refuted, 0 unknowns.
- HUNT-B a2-deficit rank from row 64, beam 400: exhausted; C2 34,400 leaves
  refuted, 0 unknowns.
- HUNT-E l1deg + hard E-prune, beam 400: DIED at row 116 -- the corridor
  |E_i| <= R-1-i empties because delta_{i,0} is greedily clamped, i.e. the
  E-walk is NOT steerable in the stock expansion.  (Diagnosis, not
  infeasibility.)
- Fix: E-BRANCHING.  Both signs of the forced k=0 cell are always admissible;
  flipping delta_{i,0}: v -> -v only shifts the residual coefficient of
  x^{d_i} by 2v (ns[d_i - 1] += 2v) -- verified exactly.  This is the unique
  control surface of the E-walk.  Real witnesses calibrate the corridor:
  R=64 slim witness max|E| = 10, fat max|E| = 4 (E-walks meander near 0).

### Terminal-state anatomy across ranks (all at (d_126,d_127) = (151,152))

- l1deg beam 400 (HUNT-A): E in [-36,-22], 0 E-viable, min def2 974,
  gaps at m=0 AND m=1.
- l1deg + rand_frac 0.25 beam 800 (HUNT-D): E in [-8,14], 147 states with
  |E| <= 2, but fat -- min def2 482, median 3.5e6 (gaps beyond m=0).
- a2-deficit rank beam 400 (HUNT-B): produces SLIM states (L1 down to 236)
  whose ONLY gap is m=0 with gap = |E| - 2 (E = -10, deficit 8).  For
  a2-steered slim states the E-wall is the ENTIRE remaining obstruction.
- Corollary attempts: E-necessity prune alone (HUNT-E l1deg beam 400) DIES at
  row 116 (E not steerable in the stock expansion); adding E-branching with a
  hard corridor (HUNT-F ecap 12; HUNT-G necessity-only; HUNT-H a2 + both):
  the beams then keep 100% E-viable terminal states (|E| <= 2) but they are
  catastrophically fat (L1 ~ 2e72) -- the E-drift of SLIM lineages is
  systematic, and pruning it away leaves only hopeless fat states.  The
  corridor and slimness are in structural tension under low-coefficient
  steering.
- HUNT-C (a2 rank from row 90, beam 800, no E-machinery): 261/800 terminal
  states E-viable spontaneously (the deficit's m=0 term steers E), min
  deficit 42030 with gap profile {m=1: 26, m=2: 42004} -- the wall is a
  LADDER of dual moment functionals b_m = sum_j s'_j C(db-j, m): fixing m=0
  (E-wall) exposes m=1 (first moment, nearly closed at 26), then m=2.  The
  a2-deficit rank squeezes the ladder level by level; the steering grid
  controls low-x coefficients while the ladder lives at the TOP coefficients,
  reachable only through delta_{i,0} flips (exact, +-2 on one coefficient)
  and even delta_{i,1}, delta_{i,2} shifts -- the identified control surface
  for future work.

### Refutation anatomy -> the E-wall (the structural find)

HUNT-A (exact R=64 winning recipe, l1deg beam 400) exhausted: C1 refuted all
400 terminal states, C2 all 15,680 grid leaves, with ZERO 'unknown' outcomes
-- every refutation was an exact necessity certificate, never a search-budget
failure. Anatomy of the 400 C1 refutations: 400/400 pure interval gap,
0 parity failures, 0 degree overflows; deficits 974..~5e5, and for the best
state the gap concentrates entirely at m = 0 (gap 22) and m = 1 (gap 952) --
the tiny-box low-k cells of the final rows.

Decoding m = 0: b_0 = sum_j s'_j = sigma'(1) is the EVALUATION functional.
Since the top cell delta_{i,0}... (k=0 cell, B_{d,0} = x^d) has box C(d,0)=1
and odd parity, delta_{i,0} = +-1 is forced EVERY row, and
E_i := sigma_i(1) obeys E_{-1} = 0, E_i = E_{i-1} - delta_{i,0}: a forced
+-1 walk that must return to 0 at closure. Exact necessity:

    |E_i| <= R - 1 - i        (parity is automatic),

an O(1)-per-state hard prune the l1deg rank is blind to. Measured: HUNT-A's
entire terminal beam had E in [-36, -22] (need |E| <= 2) -- the beam was
E-dead long before the endgame; this is WHY the R=64 recipe fails at R=128
at beam 400, and quantifies the m=0/m=1 low-cell wall (cf. the epoch-level
evaluation wall, THM-2977) at the residual level. The m=1 row of the dual
transform gives the analogous first-moment envelope (step size <= d_j per
row) as the next-order obstruction.

Beam telemetry (l1deg run): the S_2 deficit against the next two profile
degrees sheds like the L1 (e.g. row 123: L1 ~ 3.0e5, def2 = 27280; row 125:
def2 = 1010), confirming it as a usable steering functional; mid-search
deficits are astronomically large exactly per the effective-capacity picture.

Independent re-verification of the concurrent-angle witnesses
(`--verify-existing`): `R128_direct`, `R128_lp`, `R128_rule_boxeph` all pass
admissibility + epoch identity under this session's referee.

## Hazards respected

Search negatives are never treated as infeasibility (THM-3029); the only
exact refutations claimed are the L=2 necessity certificates (parity precheck
/ first interval gap), which are theorems, not search outcomes. All arithmetic
integer-exact; every witness re-verified before `verified:true`.
