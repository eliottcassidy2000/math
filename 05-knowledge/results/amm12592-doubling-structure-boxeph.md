# AMM 12592: doubling structure of the gamma* floor-profile epochs — backbone/Lucas reduction, a deterministic doubling map, and the FIRST R = 64 floor closure

boxeph, 2026-08-03. Session deliverable for the "uniform scheme — doubling structure" lane.
Referee: `04-computation/amm12592_doubling_structure_boxeph.py`
Output: `05-knowledge/results/amm12592_doubling_structure_boxeph.out` (27/27 checks PASS)
New witness: `04-computation/amm12592_floor_witness_R64_boxeph.json` (verified exact)

**Outcome class: (iii) + a hard new closure.** No all-R construction yet; instead
(a) a structural reduction of every witness to backbone + Lucas polynomial + small
corrections, (b) a deterministic doubling map that **closes the gamma\* floor profile
at R = 64 with D0 = 0** — extending THM-3029's (A) from n <= 63 to **n <= 127**, on
the epoch the beam could not crack — and (c) a precisely quantified obstruction to
iterating the map, with the decisive next test isolated.

## 1. Headline result: R = 64 closes at the floor profile

The profile `d_i = floor(gamma*(64+i))`, D0 = 0, closes: all 64 blocks admissible
(Lucas box + parity), epoch identity `sum_i x^i Delta_i = (1-x)^63` exact over Z.
Effective rate 0.597938 (the 3/5 profile has 0.600000). Hence

```text
C = 1 + gamma* = log_5(5 phi^2) = 1.5979874356654402   is ATTAINED for n <= 127,
```

matching the proved archimedean floor (THM-3009/3017/3024) on that whole range and
superseding `C <= 8/5` there. The direct beam (`amm12592_gamma35_beam_deathstar.solve`)
FAILS at this profile at every width tried (250–900, ctrl 2–3, span 2–3, and an
attractor-scored variant) — always "exhausted at final row" — yet the profile closes:
one more demonstrated instance of THM-3029 sec. 1 (beam negatives are artifacts).

Provenance of the witness: **not search** — a deterministic doubling map applied to the
committed R = 32 floor witness (below), then verified exactly from scratch.

## 2. Structure found in the witnesses (task 1)

Every committed floor witness (R = 8, 16, 32) has the same three-phase anatomy:

* **Backbone.** `Delta_i = (p - q) + c_i` for `1 <= i <= R-2`, where `(p-q)` is the
  ballot block `[binom(d-1,k) - binom(d-1,k-1)]_k` (admissible at EVERY degree >= 1),
  and `Delta_{R-1} = -1` (minus the full box). The corrections `c_i` are **even in
  every cell** (parity class of the backbone) and supported in a short high-degree
  window that marches linearly with the row.
* **Endgame attractor (exact).** `E_m := -1 + x + ... + x^m` satisfies
  `E_m - (2x-1) = x E_{m-1}`; once the residual hits some `E_m` the epoch closes
  deterministically by emitting `(p-q)` blocks down to `E_0 = -1` = last row. All
  three witnesses enter the attractor (after rows 3 / 8 / 28 resp.). This is the
  closed-form endgame the beam keeps failing to find at R = 64.
* **Slimness.** max|c| = 8 / 3432 / 184756 at R = 8 / 16 / 32 — corrections are
  tiny against the boxes (~1e-4 of capacity at R = 32).

**Reduced identity (new, exact).** Summing the backbone geometrically, the epoch
identity is EQUIVALENT to

```text
q * C = p^R + q^R - p(p-q),      C := Delta_0 + sum_{i=1}^{R-2} p^i c_i .
```

`p^R + q^R = L_R(pq)` is the Lucas polynomial in u = pq (p + q = 1); its coefficients
`(R/(R-m)) binom(R-m,m)` have Fibonacci mass `sum_m binom(R-m,m) = F_{R+1} ~ phi^R`.
**This is where phi enters the construction side**: the corrections must carry a
Lucas/Fibonacci-weighted deficit, against Lucas-box capacities — the same
phi-vs-capacity balance whose threshold is gamma* (THM-3002 sec. 4). And the Lucas
doubling identity

```text
L_2R(u) = L_R(u)^2 - 2 u^R        (verified exact, R = 8, 16, 32)
```

is precisely THM-2160's middle pair `-2(pq)^R` at epoch scale: the epoch-doubling
of the reduced problem is "square the corrections, then absorb one middle-pair term".
This gives HYP-9061's two-ray/middle-pair reading a construction-side mechanism.

## 3. The doubling map (tasks 2–3)

THM-3026 leaves the naive square `q^{2R-1} = q(sum_i p^i Delta_i)^2` with row pile-up
~ pair count (overshoots 49/276/1541). Two exact moves collapse it:

* **q-split differencing.** With `D_j := sum_{i+i'=j} Delta_i Delta_{i'}` (cell
  convolution = product, by (M)), take `E_j := D_j - D_{j-1}`, lifted to the 2R floor
  profile. Then `sum_j p^j E_j = q^{2R-1}` EXACTLY; degrees fit at ANY rate by floor
  superadditivity (no +1 loss: the q is realized as the difference); backbone products
  telescope (`Delta_{i+1} - Delta_i = c_{i+1} - c_i`). Measured worst cell overshoot
  of the carve: **7 / 9 / 25** at 8->16 / 16->32 / 32->64 (naive lifted `q*D_j` rows:
  6 / 11 / 31; THM-3026's 49/276/1541 measured raw row sums without the lift-smoothing).
* **Cell-space carry sweep.** Division by x is EXACT in cell space: a block at degree
  d with top cell 0 re-reads at degree d-1 with the same cells. So sweep rows
  j = 0..2R-1: X := carve row + carry; emit the parity-correct cellwise clamp of X
  into the box; the leftover becomes the carry one degree down. The only forced cell
  is the top (must be +-1), steered by a tau-aware boundary recursion
  `rho(j, m) = f(rho(j+1, m-1))` prescribing the top W carry cells (W = 10–14).

**Results.** 8->16 closes (18 clamp fallbacks), 16->32 closes (22), and 32->64 closes
with **ZERO fallbacks** — at 32->64 the map is a completely deterministic algorithm:
carve, prescribe, clamp. Output verified exactly; that is the R = 64 witness of sec. 1.

**Rate cost of naive doubling (task 2 answer).** Naive squaring needs
`D0(2R) >= 2 D0(R) + O(log R)` to absorb pile-up by slack alone — D0 grows linearly
in R, i.e. a CONSTANT rate penalty (useless). The carve+sweep needs no slack at all
where it closes (D0 = 0 preserved through 32->64). Also note the softening: since
`C* = inf{C : exists D}`, per-epoch slack `D0(R) = o(R)` suffices for
`C* = 1 + gamma*`; bounded D0 is NOT required. Every closure here nevertheless has
D0 = 0.

## 4. The obstruction to iterating, quantified (the honest negative)

64->128 under this sweep does NOT close, and the failure is understood:

* **Fat output.** The map's R = 64 blocks have max|c| = 3.4e21 (box-scale corrections
  at the clamped cells) versus 1.8e5 for the R = 32 input — the sweep saturates boxes
  where the carve overshoots, so output slimness is destroyed even though the epoch
  closes.
* **Transport runaway.** With fat input, the 64->128 carve seeds carries that grow
  ~x2.6/row (lift spread 2^{1+gamma*}) against box growth ~x2^{gamma*}/row; the carry
  hits ~1e26 by row ~122 and the forced top cell leaves {+-1}. Slack does not help
  (D0 up to 12 tried — caps x4096 — still runaway), nor do wide prescriptions
  (W = L up to 100; deep chains amplify through the binomial cross-terms).
* **Multi-folding does not bypass slimming.** Carving 2^s-fold products of a slim
  seed with (2^s - 1)-th differences is exact at every level, but worst overshoot
  grows 25 -> 341 -> 5739 -> 1.4e6 for (32,s=1), (16,s=2), (32,s=2), (16,s=3):
  differencing gains factorial-in-s smallness only on smooth mass, and the c-windows
  are not smooth.
* Per THM-3029 sec. 1, ALL these negatives are artifacts of THIS transport scheme —
  none is evidence that the R = 128 floor epoch fails to close.

**Precise obstruction statement.** The carve+sweep doubling map closes R -> 2R
whenever the input is slim (measured: corrections ~1e-4 of box at R = 32 suffice;
box-scale corrections at R = 64 do not). The map does not preserve slimness. So the
missing piece of the all-R induction is exactly a **slimming pass**: any R = 64 floor
solution with witness-grade corrections closes R = 128 by the same mechanism that
closed 32->64.

## 5. Decisive next tests

1. **Slim the R = 64 witness** (same-scale transport: re-distribute the saturated
   cells across neighbouring rows using the even-correction freedom, objective
   max|c| -> witness scale ~1e7), then run the existing map for 128. This is now a
   well-posed finite optimization, not a search in the dark.
2. **Emission policy that preserves slimness**: emit `ballot + bounded-c` instead of
   the full clamp, accept a larger (but structured) carry, and analyze the resulting
   carry recursion — the tension slim-blocks vs small-carry is the mathematical crux
   of the induction, now isolated in one dial.
3. **Reduced-problem doubling**: work at the level of `q C_R = L_R(pq) - p(p-q)`
   directly; `C_{2R} = q C_R^2 + 2p(p-q) C_R - p(p-q)(3p+q) - 2 p^R q^{R-1}` (exact),
   so the doubled corrections are squares/shifts of the slim `c_i` plus ONE explicit
   middle-pair term `-2 p^R q^{R-1}` — the pile-up never touches the backbone.
4. For the write-up of THM-3029's successor: the endgame attractor lemma explains the
   universal "exhausted at final row" beam mode; any future solver should target
   `sigma = E_m` explicitly.

## 6. Files

* `04-computation/amm12592_doubling_structure_boxeph.py` — referee (27 exact checks).
* `05-knowledge/results/amm12592_doubling_structure_boxeph.out` — expected output.
* `04-computation/amm12592_floor_witness_R64_boxeph.json` — the R = 64 gamma* floor
  witness (profile, 64 blocks, verified flag); re-verified inside the referee run.
* Exploration scripts (scratchpad, not committed): carve/sweep prototypes.

Coordination note: another agent attacks R = 64 directly by search; this witness is
independent of that lane and settles the epoch either way. The all-R statement
(hence `C* = log_5(5 phi^2)` with THM-3024's floor) remains OPEN, blocked exactly on
the slimming pass of sec. 5.
