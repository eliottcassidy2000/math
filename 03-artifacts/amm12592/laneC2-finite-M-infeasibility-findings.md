# Lane C2 findings: exhaustive finite-M feasibility of gamma < 1 depth laws (AMM 12592 / HYP-9061)

Session: death-star-2026-07-30-coinC2, lane G3 (= repo lane C; resumed).
Script: `04-computation/amm12592_laneC2_finiteM_bnb_deathstar.py` (exact
integer DP + targeted DFS modes dfs/corridor/pin/witness/extend; floats
appear in logs only). Raw run logs archived in
`05-knowledge/results/amm12592_laneC2_finiteM_bnb_laneG3.out`.

Frame = THM-2966 spine normal form with laneD ledger conventions: depth law
`d_m = floor(gamma*m) + D0`, envelope `T(m) = m + 1 + d_m`; row-m cells on
anti-diagonal `A_m = m + d_m + 1`; doubled deficits `delta` integer with
`|delta| <= binom(d_m,k)`, `delta == binom(d_m,k) (mod 2)`;
`N_M(p) := 2 D_M(p) = sum_{m <= M} sum_cells delta p^z q^o`.

## 1. Method: exact DP whose emptiness is a theorem

**Necessary envelope (PROVED; THM-2966 eq (4) + laneD Sec 3).** Every fair
extractor with deadline `T(m) = m+1+d_m` satisfies, for every `M >= 1` and
every `p in (0,1)`:

    |N_M(p)| <= p^{M+1} + q^{M+1}.                              (ENV_M)

(Mechanism: `N_M(p) = 2 p^{M+1}(1/2 - G_{0^{M+1}}) + 2 q^{M+1}(1/2 -
G_{1^{M+1}})` with the boundary shares `G in [0,1]`.)

**Exact mass bookkeeping (PROVED, one-line).** Row m's total cell mass is
`sum_k binom(d_m,k)(p^{m+d_m-k}q^{k+1} + p^{k+1}q^{m+d_m-k}) = p^m q + q^m p`,
so the total mass of all rows `> M` is exactly `p^{M+1} + q^{M+1}` — the
same quantity as the envelope. Hence (ENV_M) says precisely: "the remaining
rows still have enough mass to cancel the current partial sum". The trivial
(mass-only) relaxation is therefore never infeasible; any finite-M
infeasibility is forced by the integrality + parity structure alone.

**DP.** State = vector of exact integers `(N(p_j) * b_j^{A_m})_j` at rational
biases `p_j = a_j/b_j` (all row-m monomials have level `z+o = A_m`, so this
scaling is exact). Cells are processed row by row (within a row:
capacity-descending); a state is kept iff at every bias

    |S| <= R := (mass of not-yet-assigned cells of the current row)
                + (p^{m+1} + q^{m+1}),                          (PRUNE)

all scaled to integers by `b^{A_m}`. States are deduplicated; at row
boundaries, if the bias set is closed under `p -> 1-p`, mirror orbits are
canonicalized (exact: see Lemma below).

**Prune-safety Lemma (PROVED).** Fix an assignment of rows `1..M` satisfying
(ENV_1)...(ENV_M) at bias p. Then every mid-row partial sum passes (PRUNE).
*Proof.* Let `S_part` be the partial sum inside row m and `S_m` the full sum
through row m. `|S_part| <= |S_m| + (mass of the row-m cells not yet added)
<= (p^{m+1}+q^{m+1}) + (remaining row mass) = R`. QED.

**Exactness Corollary (PROVED).** The DP state set at the end of row M equals
exactly the set of evaluation vectors of assignments of rows `1..M`
satisfying (ENV_1)...(ENV_M) at the chosen biases (prunes remove nothing
that satisfies the boundary envelopes, by the Lemma; and the final cell's
(PRUNE) of each row *is* (ENV_m), so survivors satisfy all boundary
envelopes inductively). Consequently:

- **If the state set empties at row M**: no assignment of rows `1..M`
  satisfies all boundary envelopes at the listed biases; since (ENV) is
  necessary, **no fair extractor with this depth law exists — a theorem**,
  with exact minimal infeasible M (states at M-1 were nonempty).
- The boundary state set is order-independent and mirror-symmetric for
  conjugation-closed bias sets, so orbit canonicalization is exact.

**Selftest (VERIFIED-EXACT).** `--selftest` compares the DP boundary set
against unpruned brute-force enumeration filtered by the boundary envelopes:
gamma=1/2 M=3 biases {1/2,2/3} (61 = 61 states, sets equal); gamma=1 M=2
biases {1/2,3/4} (246 = 246, equal); mirror-canonicalization consistency at
{1/2,1/3,2/3} (36 = 36 orbits, raw sets equal). All PASS.

**Parity-relaxed control (PROVED, no computation needed).** If the parity
constraint `delta == binom (mod 2)` is dropped (keeping the box), the
all-zero assignment satisfies every (ENV_M) with `N_M == 0`. So the
relaxed problem is feasible for every depth law at every M: **any finite-M
infeasibility found below is caused by parity (Lucas-forced odd cells),
not by magnitudes.** Likewise the real-LP relaxation is trivially feasible
(laneD Sec 5). Parity is the entire obstruction content.

## 2. Full-frontier DP runs (exhaustive corridor, small M)

### Run B — gamma=1/2, D0=0, biases {1/2, 2/3, 3/4}

    --gamma 1/2 --D0 0 --biases 1/2,2/3,3/4 --Mmax 24 --state-cap 1000000

RESULT: FEASIBLE_THROUGH_5 (state cap 1.65e6 > 1e6 hit inside row 6).
Exact boundary sets: 3, 12, 69, 3154, 47146 states at M=1..5. Exact
corridor at p=1/2 (scaled by 2^{A_M}): minN=-envN, maxN=+envN, min|N|=0 at
every M<=5 — the corridor FILLS the envelope box at 1/2.

### Run C — gamma=1/2, D0=0, biases {1/2, 1/3, 2/3} (mirror-closed)

RESULT: FEASIBLE_THROUGH_5 likewise (orbit counts 3, 9, 36, 1378, 16527);
same full corridor at 1/2. The DP's full-frontier enumeration is the wrong
tool beyond row 5 (state explosion ~ 1.6e6); targeted DFS below supersedes
it. (Run A/D superseded before launch by the G3 modes.)

## 3. Lane G3 targeted experiments (same script, --mode dfs/corridor/pin/witness)

All runs use the 9-bias set
`P9 = {1/2, 1/3, 2/3, 2/5, 3/5, 1/4, 3/4, p_A=1285/2181, p_B=8847357/11821757}`
(mirror-closed core + both certificate biases), exact integer arithmetic at
scale `b^{Lmax}`, `Lmax = A(Mtry)`. First-witness DFS with the sound prune
of Sec. 1; INFEASIBLE = exhausted search = theorem for that (M, biases);
FEASIBLE witnesses are certified by independent re-evaluation (assert).
Referee cross-check: corridor mode at M=5, biases {1/2,2/3,3/4} reproduces
the DP's exact corridor (achievable = even values 0..8 = full envelope box,
odd = infeasible) — PASS.

### 3a. Depth-law sweeps (task a): gamma = 3/5 and 3/4, D0 = 0, parity ON

    gamma=3/4: FEASIBLE at every M = 1..12 (nodes: 3 -> 16653; < 0.1 s each).
    gamma=3/5: FEASIBLE at M = 1..9 (nodes: 3 -> 11628), then M=10
               UNRESOLVED: abs-order CAP 9.0e7 nodes / 440 s, desc-order
               (saturation-first) CAP 4.8e7 / 190 s.

No finite-M infeasibility at these rates through the exhaustible range:
the 9-bias boundary-envelope relaxation stays satisfiable at gamma = 1/2,
3/5, 3/4 well past band birth (m* = 4, resp. 5, 8 at D0=0). EVIDENCE that
the necessary-envelope-at-finitely-many-biases relaxation does NOT detect
the gamma < 1 obstruction at small M — the greedy freeze (laneD) is a
policy failure, not a value-space failure at these depths. Search-hardness
onset (first M needing > 4e7 nodes without verdict): M=10 at 3/5, M=11 at
1/2 (main session + desc both CAP), M=11 at 1/3 feasible but 7.4e7 nodes
(main session), > 12 at 3/4.

### 3b. Parity OFF controls (task b) — PROVED trivial, verified

With `delta == cap (mod 2)` dropped, all-zero satisfies every ENV_M
(Sec. 1); the DFS confirms instantly (node count = cell count) at
gamma = 3/5, 3/4, M <= 12. Parity is the entire finite-M content.

### 3c. Corridor geometry at (1/2, 0) with all 9 biases (task c)

Exact achievable SET of final `N_M(1/2) * 2^{A_M}` (eq-target DFS per value,
global-negation symmetry gives the negative half):

    M=4: envN=8   achievable = {0,±2,±4,±6,±8}      (all even values; FULL box)
    M=5: envN=8   achievable = {0,±2,...,±8}         FULL
    M=6: envN=16  achievable = {0,±2,...,±16}        FULL
    M=7: envN=16  achievable = {0,±2,...,±16}        FULL
    M=8: envN=32  achievable = {0,±2,...,±32}        FULL

Odd values are infeasible by FORCED PARITY (proved: scaled sum == sum of odd-
weight caps mod 2, no search needed); every even value in [-envN, envN] is
achieved. **The 9-bias envelope system does not narrow the corridor at 1/2
at all through M=8**: the corridor width equals the envelope width exactly,
including sitting ON the boundary (|N| = envelope), and min |N| = 0.

At M=9 (envN=32): achievable v>=0 = {0,2,...,28} at cap 6e6, plus v=30
FEASIBLE on retry (2.9e7 nodes); only the exact-envelope point v=32 stays
undecided (CAP at 4.5e7 nodes/200s). So the M=9 corridor CONTAINS all even
values up to (15/16) envelope — the narrowing reported in an earlier
checkpoint of this file was search hardness, not certified absence. At
M=10 (envN=64) the descending scan v=64,62,...,40 is ALL undecided at 8e6
nodes each (job timeout at v=40), while the certified M=10 witness
realizes v=30: corridor max in [30, 64]. **Certified statement: through M=9 the corridor at 1/2 loses
AT MOST its extreme boundary point; no interior hole exists anywhere in
the tested range.** The exhaustible corridor questions end at M~9-10:
past it only witnesses (corridor lower bounds), never certified absences,
are affordable at ~1e8 nodes.

### 3d. Witness anatomy at (1/2,0), M=10 vs in-lane greedy (task d)

Witness = first ENV-feasible assignment found by abs-ordered DFS (1562
nodes), certified by independent exact re-evaluation. Log:
`g3_witness_g12_M10.log` (scratchpad; verbatim cell list).

- **Envelope-riding is what survival looks like.** The witness's row-end
  max-over-biases ratio rho_m = max_j |N_m(p_j)|/env_m(p_j) climbs
  0.00, 0.50, 0.67, 0.89, 0.91, 0.96, 0.99, 0.9990, 0.9984, 0.9977 for
  m=1..10: from row 4 (= band birth m* = 4) on, the surviving trajectory
  sits ESSENTIALLY ON the proved envelope. Since the DFS tries small |delta|
  first, every smaller-|delta| completion of the row-1..3 prefix was
  exhausted or soundly pruned: **survival to M=10 requires committing
  near-maximal deficits (pre-loading)** — the opposite of greedy's
  zero-chasing.
- **Saturation profile.** Rows 8-10 of the witness take delta = +cap on
  nearly every cell (e.g. row 9: all ten cells saturated +cap; row 10: 11
  of 12 saturated, one cell 8 of cap 10); rows 1-3 carry small mixed signs;
  the trajectory N_m(1/2)*2^16 runs -8192, -4096, -2048, -1024, -640,
  -384, -160, -32, +30 — a controlled geometric decay toward the boundary
  crossing, not a hover at 0.
- **In-lane greedy dies at M=4.** The policy "each cell takes the parity-
  allowed value minimizing |N(1/2)|" holds N(1/2) = 0 forever but violates
  the proved envelope at other biases from M=4 on (first violation exactly
  at band birth). 48 of 70 cells differ from the witness. This is the lane-D
  freeze seen value-side and certified: greedy is dead by row 4, while
  feasible assignments exist through row 10+ — the gamma<1 obstruction at
  small M is a POLICY failure (myopia), not a value-space failure.
- **Forced-parity referee PASS (PROVED cross-check).** Homogenizing both
  assignments to level A_10 = 16 reproduces beta_M(x) = (1+x)^16 +
  (1+x^11)(1+x)^5 mod 2 coefficient-for-coefficient. Band positions
  o in [7,9]: witness coeffs (-120, 84, 228) vs greedy (-1098, 0, 1098) —
  same forced parity, but survival keeps band magnitudes ~5-10x smaller.
- **The binding biases are a two-bias pinch that includes p_A.** Per-bias
  rho_m of the witness (exact, floats for display): argmax bias runs
  1/2 (m=2), 3/5 (m=3-5, up to 0.914), p_A (m=6-7, 0.957/0.988), then
  1/3 (m=8-10, 0.998-0.999) with p_A second at 0.979-0.997. The certificate
  bias p_A = 1285/2181 = 0.589 is (near-)binding at EVERY level m >= 3,
  and the late-time pinch pairs it with a q-heavy bias (1/3) on the other
  side of 1/2 — a two-ray shape consistent with HYP-9061's Legendre/dual
  reading of (27). EVIDENCE (one witness at one depth law), mechanism:
  exact rho recomputation from the certified witness.

### 3e. Local rigidity at the M=11 frontier, (1/2, 0)

Extend mode (clamp solved rows, free the last `back` rows; every
INFEASIBLE is an exhausted search = certified for that clamp):

    M=11 back=1 (row 11 free):        INFEASIBLE (91 nodes)
    M=11 back=2 (rows 10-11 free):    INFEASIBLE (227 nodes)
    M=11 back=3 (rows 9-11 free):     INFEASIBLE (259 nodes)
    M=11 back=4 (rows 8-11 free):     INFEASIBLE (5.2e6 nodes, exhausted)
    M=11 back=5,6:                    node/time-capped (undecided)

So the M=10 witness is DEAD as a prefix: any ENV-feasible assignment at
M=11 must differ from it somewhere in rows <= 7 (PROVED for this witness).
Meanwhile the main-session first-witness DFS at (1/2,0) M=11 node-capped
at 1.5e8 nodes (their out file), and a saturation-first (desc value-order)
full DFS at M=11 also capped (4.6e7 nodes / 210 s): M=11 at (1/2,0) is the
live frontier — either infeasible or requiring a globally different
prefix, unresolved by three independent search orders.

### 3f. Deep-row pinning at (1/2,0), M=10 (task c, full table for rows 9-10)

Per-cell feasible |v| sets (testing v >= 0; sets are sign-symmetric by
global negation; "undecided" = clamp search node-capped at 4e6):

    row 9:  (k=0) |v|=1 PINNED   (k=1,side0) {4} conf, {0,2} undecided
            (k=1,side1) {0,2,4}  (k=2,side0) {6} conf, {0,2,4} undecided
            (k=2,side1) {0,2,4,6}  (k=3,side0) {2,4} conf, {0} undecided
            (k=3,side1) {0,2,4}  (k=4) |v|=1 PINNED
    row 10: (k=0) |v|=1 PINNED   (k=1,side0) {3,5} conf, {1} undecided
            (k=1,side1) {1,3,5}  (k=2,side0) {6,8,10} conf, {0,2,4} und.
            (k=2,side1) all      (k=3,side0) {4..10} conf, {0,2} und.
            (k=3,side1) all      (k=4) {1,3,5} both sides
            (k=5) |v|=1 PINNED

Structural signal: **the 0-side deep cells are forced toward saturation
(only near-|cap| values certified feasible), while every 1-side cell keeps
its full value range.** Mechanism: the bias set is mirror-closed except for
the two certificate biases, which both have p > 1/2 and hence weight the
0-side monomials (z large) — the certificate-side pressure is what pins.
EVIDENCE grade (undecided small-value entries are capped searches, not
certified infeasible).

### 3g. The zero-slice at 1/2: exact balance is possible but hardens at M=9

Mode `--zero-bias 1/2` replaces the envelope at 1/2 by the far stronger
demand `N_m(1/2) == 0 at EVERY row end` while keeping ENV at the other 8
biases:

    M=1..5: FEASIBLE (<40 nodes)     M=6: FEASIBLE (407)
    M=7: FEASIBLE (8.6e3)            M=8: FEASIBLE (5.5e5, 2.3 s)
    M=9: CAP (1.5e7 nodes / 60 s)

So perfect balance at the symmetric bias is NOT what kills greedy (a
non-myopic scheme holds N(1/2) == 0 through M=8 with the other envelopes
intact); but the zero-slice hardness onset (M=9) precedes the general
hardness onset (M=11) — the first structural pinch appears exactly where
the corridor's extreme point becomes unresolvable. EVIDENCE via certified
witnesses M<=8; M>=9 undecided.

### 3h. Dyadic clock cross-checks (D3 of the session hand-derivations)

Exact computation of beta_M for (1/2, 0), M <= 16: beta_M == 0 precisely
at M = 1, 3, 7, 15 (i.e. M+1 = 2^r; forced-odd weight 2-8 otherwise) —
the claimed checkpoint-vanishing clock is VERIFIED-EXACT here. And the
deep-corner condition (d_M = 2^rho - 1 all-ones AND A_M = 2^r - 1) demands
2^(r-1) = 3*2^(rho-1) at gamma = 1/2, which has no solutions — PROVED
one-line: the deepest-corner clock NEVER fires at gamma = 1/2, consistent
with lane G2's restriction of the infinite clock to gamma = 1/(2^j - 1).
Suggestive timing: the free levels M = 7 (beta_7 == 0) and M = 8 (row 8's
caps binom(4,k) with 2^2-block structure) are exactly where full-box
corridor persists; the first never-free stretch M = 8..14 is where
hardness sets in (M=9 zero-slice, M=11 general). SPECULATION as stated.

## 4. Interpretation (checkpoint; final numbers in Sec. 5)

1. **No finite-M infeasibility theorem yet** at gamma in {1/3, 1/2, 3/5,
   3/4}, D0=0, with the 9-bias envelope relaxation, through the exhaustible
   range. But the search-hardness onset is sharply depth-law-dependent:
   first M needing > 5e7 nodes is M=11 at gamma=1/3 (feasible, 7.4e7),
   M=11 at 1/2 (UNRESOLVED at 1.5e8), M=10 at 3/5 (UNRESOLVED at 9e7),
   > 12 at 3/4 (1.7e4 nodes at M=12). Hardness arrives when the feasible
   set thins to envelope-riding trajectories, a few rows after band birth.
2. **Parity is the entire content** (PROVED trivially feasible without it),
   and within parity the corridor at 1/2 is EXACTLY the full even box
   through M=8, still >= (15/16) envelope at M=9 (only the on-envelope
   extreme point unresolved), no interior hole anywhere tested. Through
   M=9 the finite-bias envelope necessary condition has essentially zero
   discriminating power at 1/2 beyond the mod-2 class — consistent with
   THM-2977's evaluation-wall verdict (single evaluations cannot obstruct).
3. **What survival requires (task d):** pre-loading to the proved envelope
   (rho -> 0.99+ from band birth on), near-saturation of deep rows
   (0-side pinned by the p > 1/2 certificate biases), and a two-bias pinch
   ({1/3, p_A} at late rows) — while zero-chasing greedy is certified dead
   at band birth (M=4). The gamma < 1 obstruction, if it materializes at
   finite M, will be the closing of the envelope-riding corridor, not the
   greedy freeze.

## 5. Bottom line and next obligations

- **No unconditional finite-M infeasibility theorem was reached** in the
  exhaustible range (M <= ~10-12 depending on gamma); the assigned
  falsifier target ("certified infeasible M") remains open, now with the
  precise frontier localized: (1/2,0) M=11 and (3/5,0) M=10, both
  unresolved at ~1e8 nodes under three different search orders. These are
  the natural targets for a stronger prover (dedicated C search, better
  bounding, or an LP/parity-hybrid dual) — NOT for more Python DFS.
- **Structural yield (certified):** parity-only content; full corridor
  through M=8 and no interior hole through M=9; witness envelope-riding
  with the certificate bias p_A near-binding at every level and a late
  {1/3, p_A} two-bias pinch; M=10-witness prefix rigidity at M=11
  (rows <= 7 must change); 0-side saturation pinning; zero-slice at 1/2
  feasible through M=8, hardening at M=9; beta-clock M+1 = 2^r verified,
  deep-corner clock provably never fires at gamma = 1/2.
- **Reading for HYP-9061:** everything seen is consistent with (and
  sharpens) the rate/entropy-dual reading of certificate (27): the finite
  system's binding structure IS two-bias (one ray near p_A, one q-heavy),
  and the obstruction to gamma < 1 — if real — lives in the closing of an
  envelope-riding corridor at rates too slow for any finite-bias
  evaluation test to certify cheaply, exactly the regime a tropical/dual
  certificate with a per-M rate margin (the 1/25) would govern.

Raw logs archived verbatim in
`05-knowledge/results/amm12592_laneC2_finiteM_bnb_laneG3.out`
(dfs/corridor/pin/witness/extend/desc + the Sec. 2 DP logs); all
regenerable from the mode/flag combinations recorded above with the
script at head of file (final `--selftest`: PASS).
