# Lane G5 findings: the construction gate, taken seriously — exact checkpoint-closure program

**Session:** death-star-2026-07-30-coinC2, lane G5.
**Script:** `04-computation/amm12592_laneG5_construction_gate_deathstar.py`
(modes: `parity`, `dfs`, `depth0`, `threshold`, `solve`, `abeldini`).
All asserted computations exact integer arithmetic; floats only orientation.
Frame: THM-2966 spine normal form, laneD ledger conventions, THM-2976
binary-clock parity (T1), laneF2 policy findings absorbed.

## 0. Summary of what this lane establishes

1. **(PROVED, sec. 1) The checkpoint-closure reduction.** Closing the books
   *exactly* (`D ≡ 0` as a polynomial identity) at every dyadic checkpoint
   `M = 2^r - 1` is a *sufficient* program for `C* <= 1 + gamma`: if the
   base system (rows `1..2^{r0}-1`) and every epoch system (rows
   `2^r..2^{r+1}-1`) admit box+parity deficits summing to the zero
   polynomial, the concatenation is an exactly fair extractor with pathwise
   deadline `T(m) = m + 1 + d_m`. The proved necessary envelope holds
   automatically. No approximation, no tail analysis left.
2. **(PROVED, sec. 2) Parity is free exactly at dyadic cuts** — for ANY
   per-row depths, by a two-line Frobenius argument; and non-dyadic cuts
   are parity-obstructed in all 18 hostile controls. This is (D3)/THM-2976
   T1 sharpened into the exact reason the checkpoint program is the right
   normal form for a construction.
3. **(VERIFIED-EXACT, sec. 3) The gate is real, and it is NOT at the small
   epochs:** every tested rate `gamma` (including `2457/6592, 1/3, 3/8,
   4/11`) closes the base `1..7` and epochs `[2,3], [4,7]` exactly; the
   `gamma -> 0` limit (all depths 0, `C -> 1`) is exhaustively INFEASIBLE
   from epoch `[4,7]` on — consistent with (and a micro-proof mechanism
   for) THM-2967.
4. **(sec. 4) The threshold `gamma*(r)`** = least rate whose epoch
   `[2^r, 2^{r+1}-1]` closes: measured per `r`; comparison with
   `alpha = 2457/6592` is the construction-gate test of certificate (27).
5. **(sec. 5) Abel-Dini hunt** for `(896, 1285)`, `(2974400, 8847357)` as
   consecutive partial sums / structured products of the ledger's
   combinatorics.

## 1. The checkpoint-closure reduction (PROVED)

Let `d_m >= 0` be any depth law. For a block of rows `B = [m_lo, m_hi]`
define the **closure system** `(E_B)`: integers `delta_{m,k}` (both sides,
`|delta| <= binom(d_m,k)`, `delta == binom(d_m,k) mod 2`) with

```text
sum_{cells of B} delta * p^z (1-p)^o  ==  0   identically in p.     (E_B)
```

**Theorem G5-1.** If `(E_B)` is solvable for `B_0 = [1, 2^{r0}-1]` and for
every epoch `B_r = [2^r, 2^{r+1}-1]`, `r >= r0`, then an exactly fair
extractor with pathwise deadline `T(m) = m + 1 + d_m` exists. In
particular, with `d_m = floor(gamma*m) + D0`, `C* <= 1 + gamma`.

*Proof.* Set `w_{m,k} = (binom(d_m,k) + delta_{m,k})/2`; box+parity make
these integers in `[0, binom(d_m,k)]` (THM-2966 (2)). Then

```text
sum_m [p^m q W_m + q^m p V_m]
  = (1/2) sum_m (p^m q + q^m p)           [box midpoint: each row's box
                                            sums to p^m q resp. q^m p]
  + (1/2) sum_r [closure sum of B_r]      [= 0 term by term]
  = (1/2) * 1  =  1/2   for all p in (0,1),
```

with absolute convergence (each row's total mass is `p^m q + q^m p`, sum
1). THM-2966's (S)⇒extractor direction finishes. The proved necessary
envelope `|D_M| <= (p^{M+1}+q^{M+1})/2` is automatic: for the assembled
scheme `D_M = -(1/2) * (tail)`, and the tail is bounded by the tail box
mass. QED.

Two remarks that kill tempting shortcuts:

- **Exact closure cannot be done monomial-by-monomial** beyond the pairing
  horizon: distinct cells share a bivariate monomial `p^z q^o` only when
  `m <= d_m + 1` (laneD sec. 1), so a *bivariate* zero forces all deltas
  zero (parity-impossible). Closure lives in the *univariate* ring: the
  Pascal identity `v(z,o) = v(z+1,o) + v(z,o+1)` (`v(z,o) := p^z(1-p)^o`)
  is the whole cancellation mechanism. The classical ratio-2 construction
  is the degenerate case where bivariate pairing never dies (`gamma = 1`
  keeps `m <= d_m + 1` forever).
- **Depth freedom costs nothing:** a row stopped at effective depth
  `e_m < d_m` refines to the full-depth box (THM-2966 budget-box
  exactness), so `(E_B)` at full depth is the most general closure system
  for the deadline; separately, T1's checkpoint vanishing is depth-blind,
  so the parity consistency below holds for any per-row depths.

## 2. Parity is free exactly at dyadic cuts (PROVED + VERIFIED-EXACT)

Solvability of `(E_B)` mod 2 requires the parity vector
`sum_{odd cells} v(z,o) mod 2` to vanish. Summing a row's odd cells gives
`(x + x^m)(1+x)^{d_m}` at level `A_m` (marking `o` by `x`); summing the
block and homogenizing to level `L`, the box-sum telescopes to

```text
hom_2 = (1 + x^{m_lo})(1+x)^{L - m_lo} + (1 + x^{m_hi+1})(1+x)^{L - m_hi - 1}
```

so the block parity vanishes identically whenever both `m_lo` and
`m_hi + 1` are powers of two — Frobenius `(1+x)^{2^r} = 1 + x^{2^r}` —
i.e. for all dyadic epochs `[2^r, 2^{r+1}-1]`, bases `[1, 2^r - 1]`, and
unions of adjacent epochs (PROVED). Conversely, comparing `(x+1)`-adic
valuations forces `(m_hi + 1) - m_lo = 2^{v_2(m_hi+1)} - 2^{v_2(m_lo)}`
as a necessary condition (PROVED), and the full converse holds in every
tested case. **Script `parity` mode: PASS for 9 rates x 2 D0 x r < 6, and
parity-INCONSISTENT in 18/18 non-dyadic hostile controls.**
This is the exact sense in which (D3)/THM-2976-T1 is a construction-side
gift: at dyadic cuts the *entire* forced-parity debt of an epoch is
internally cancellable; between them it never is.

## 3. Small systems: every rate closes; the gamma→0 limit dies (VERIFIED-EXACT)

Exhaustive DFS (class-sequential in `z`, capacity-pruned; every witness
re-verified by an independent exact polynomial-sum assert):

- **FEASIBLE with exact witnesses** for all of
  `gamma in {1, 1/2, 1/3, 2/5, 3/8, 4/11, 3/5, 3/4, 2457/6592}` (D0=0):
  base `[1,1]`, base `[1,3]`, epoch `[2,3]`, epoch `[4,7]`, base `[1,7]`.
  Node counts 2–171617. In particular **the certificate slope
  `alpha = 2457/6592` closes its books exactly through row 7.**
- **Depth-0 anchor (the `C -> 1` limit):** with all `d_m = 0`, epochs
  `[1,1]` and `[2,3]` close (the `[2,3]` witness is the antisymmetric quad
  `p^2q - p^3q - pq^2 + pq^3 = pq[(2p-1)+(1-2p)] = 0`, THM-2160's middle
  pair in miniature), while epochs `[4,7]`, `[8,15]`, `[16,31]` are
  **exhaustively INFEASIBLE** (7 / 71 / 12871 nodes). So the
  checkpoint-closure gate genuinely closes somewhere between `gamma = 0`
  and the small rates — the program sees THM-2967's wall from below, as an
  integer infeasibility, with a 7-node proof at `[4,7]`.

## 3b. Depth-0 rigidity is Prouhet–Tarry–Escott-type, in closed form (PROVED)

**Theorem G5-2.** With all depths 0, the epoch closure `(E_{[R,2R-1]})`
(signs `a_m, b_m in {+-1}`, `(1-p) sum a_m p^m + p sum b_m (1-p)^m = 0`)
is solvable iff `R <= 2`, and at `R = 2` the only solutions are the two
middle-pair sign patterns.

*Proof.* Valuation at `p = 0` forces `g(u) = sum b_m u^m` to vanish to
order `R - 1` at `u = 1`; since `g = u^R h` with `deg h = R - 1`, this
forces `h = ±(u-1)^{R-1}`, i.e. `b_m = ±binom(R-1, m-R)` up to sign —
compatible with `b_m in {±1}` iff every `binom(R-1, j) = 1`, i.e.
`R <= 2`. Mirror for `f`. Verified exhaustively: `R = 2` has exactly the
2 middle-pair solutions; `R = 4` has 0 among all `2^8` sign patterns.
QED.

So THM-2967's wall appears inside the closure frame as binomial-coefficient
rigidity: the `C -> 1` limit fails because the checkpoint books can only be
balanced by `±(u-1)^{R-1}`, whose coefficients explode. The rigidity is
robust to coarser checkpoints: double-epoch blocks `[2,7]`, `[4,15]` are
parity-consistent (dyadic endpoints; verified PASS for 5 rates x r <= 3)
yet remain exhaustively INFEASIBLE at depth 0 (21 / 925 nodes) — more
rows do not help when the boxes cannot hold binomial heights. The general-`gamma`
gate is the same tension with budget boxes `binom(d_m, k) ~ 2^{gamma m}`
against required binomial corrections — a log-race, i.e. a rate
inequality. This is the mass balance of the lane brief in its exact form.

## 4. The threshold hunt gamma*(r): the gate curve (VERIFIED-EXACT verdicts)

Definition: `gamma*(r)` = least depth rate whose epoch
`B_r = [2^r, 2^{r+1}-1]` (D0 = 0) closure is solvable.

**r = 2 (epoch [4,7]) — all verdicts EXHAUSTIVE (proofs):**

```text
INFEASIBLE: gamma = 0, 1/20, 1/10, 1/8, 1/6, 1/5
FEASIBLE:   gamma = 1/4, 2/7, 3/10, 1/3, 4/11, 3/8, 2457/6592, 2/5,
            5/12, 3/7, 1/2, 3/5, 3/4, 1
```

`gamma*(2) in (1/5, 1/4]`; the transition is exactly where the epoch's
first row gains depth (`floor(4 gamma) >= 1`).

**r = 3 (epoch [8,15]):**

```text
INFEASIBLE (exhaustive proofs): gamma = 0, 1/20, 1/10, 1/8, 1/6, 1/5
                                (nodes 71..131381)
FEASIBLE (exact witnesses):     gamma = 1/2 (support 90/104; sampled DFS)
                                [more as runs land]
```

**Headline so far:** `gamma = 1/2` closes rows `1..15` COMPLETELY (base
`[1,7]` exhaustive witness + epoch `[8,15]` witness, both re-verified by
the independent exact polynomial-sum assert): an exactly fair extractor
prefix with `D_15 == 0` *identically* under deadline `T(m) = m + 1 +
floor(m/2)`, i.e. `C = 3/2` behavior through two full dyadic epochs —
strictly past laneC2's former M=10/11 frontier and laneF2's beam-death
rows, with EXACT zero residual rather than envelope-surfing.

(r = 3 boundary refinement and r = 4 results checkpointed below.)

## 5. Abel-Dini integer hunt: NEGATIVE on telescoping, two structural facts

Full scan (`abeldini` mode, exact, exit 0). Factorizations (all exact):
`x_A = 389` prime; `b_A = 2181 = 3*727`; `x_B = 5872957 = 19*103*3001`;
`b_B = 11821757 = 31*381347` (381347 prime);
`S_B = 8847357 = 3*2949119`; `alpha = 2457/6592 = 3^3*7*13/(2^6*103)`.

1. **No consecutive-partial-sum telescoping found.** Scanned: cumulative
   checkpoint budgets `sum_{m<=2^r-1} 2^{d_m}` over ~700 rates x 7 offsets;
   binomial partial sums `sum_{k<=K} binom(n,k)`, `n < 80`; cumulative
   forced-odd counts / cell counts / corner counts at checkpoints. **No
   sequence produced `(896, 1285)` or `(2974400, 8847357)` as consecutive
   values.** The Abel-Dini construction-side reading (the biases as mass
   ratios of a ledger cascade) finds no support in the ledger's natural
   counting sequences. (VERIFIED-EXACT, negative.)
2. **`1285` is a forced-parity product numerator (exact).** At bias
   `p in {1/3, 2/3}`, `prod_{i in bits(d)}(p^{2^i} + q^{2^i})` for
   `d = 10 = 1010_2` has numerator `(1^2+2^2)(1^8+2^8) = 5*257 = 1285` —
   the certificate's `S_A` is literally the Lucas forced-parity row mass
   numerator of a depth-10 row at the simplest nontrivial bias. `896`, by
   contrast, only appears in trivial parametrizations (it *is* `q_A`'s
   numerator `2^7*7`). Suggestive, not decisive. (VERIFIED-EXACT; the
   *interpretation* is SPECULATION.)
3. **Prime-103 note + a near-alpha cluster (non-decisive).** `103` divides
   both `x_B = 5872957` and `alpha`'s denominator `6592 = 2^6*103`.
   Cumulative forced-odd cell counts hit exactly `6592` at
   `(gamma, D0) = (16/43, 1), (8/21, 4), (17/43, 4), (19/48, 2)` —
   note `16/43 = 0.37209` vs `alpha = 0.37272` (0.17% apart). With ~10^5
   scanned values this is within multiple-comparison noise; recorded, NOT
   asserted (repo has prior numerology retractions).

## 6. The mass-balance mechanism: closure as a binomial-height race (analysis; mechanism PROVED at gamma=0, shape EVIDENCE otherwise)

Generalizing G5-2's proof: in epoch `[R, 2R-1]` (`R = 2^r`), orders
`p^1..p^{R-1}` of the closure polynomial receive *only* 1-side
contributions, so, substituting `u = 1 - p`, the 1-side polynomial
`V(u) = sum_m u^m S_m(u)` (with `S_m` in the row-`m` Bernstein box of
degree `d_m`, parity-forced) must vanish to order `~R` at `u = 1` while
having `u`-valuation `>= R`. Hence

```text
V(u) = u^R (u - 1)^{R'} Z(u),    R' ~ R,  deg Z <= ~(1 + 2gamma)R - R'.
```

At `gamma = 0`: `deg Z = 0` forces `V = ±u^R(u-1)^{R-1}`, whose
coefficients `binom(R-1,j)` overflow the ±1 boxes — G5-2's rigidity. At
`gamma > 0`: box heights grow like `2^{d_m} <= 2^{2 gamma R}` while the
*minimal* coefficient height of `(u-1)^{R'} Z` over integer `Z != 0` of
degree `D = c gamma R` is a Chebyshev/lattice extremal quantity
`H_min(R', D)`; closure requires roughly `H_min <= 2^{2 gamma R}` — a
log-race between the forced binomial mass of the vanishing condition and
the budget entropy, i.e. exactly the "emitted corner mass vs absorbable
budget" rate inequality of the lane brief, now in exact polynomial form.
Its Legendre dual is a two-ray comparison (the binding `u`-scales are the
maximizing `j/R` of `binom(R', j)`-vs-budget, one on each side of the
band): **the (27) shape** — two rapidities with rational weight and a
per-`R` margin. Direction of the certificate along this reading: OPEN;
what this lane adds is that the gate function is concretely computable
(sec. 4's curve) and that its `r -> infinity` limit is *the* threshold
`gamma*` of the checkpoint-closure class.

## 6b. Canon-candidate statement (for the main session to referee/canonize)

**"Dyadic checkpoint closure" (proposed THM).** Fix any depth law
`d_m >= 0`. (i) For a block `B = [m_lo, m_hi]`, the parity vector of the
closure system `(E_B)` vanishes whenever `m_lo` and `m_hi + 1` are both
powers of two, and dyadic-difference timing
`(m_hi+1) - m_lo = 2^{v_2(m_hi+1)} - 2^{v_2(m_lo)}` is necessary
(Frobenius + `(x+1)`-valuation; sec. 2). (ii) If the base and all dyadic epochs are
closable, an exactly fair extractor with deadline `m + 1 + d_m` exists
(sec. 1; with `d_m = floor(gamma m) + D0` this gives `C* <= 1 + gamma`).
(iii) At `d == 0` the epochs `[R, 2R-1]` are closable iff `R <= 2`
(binomial rigidity `V = ±u^R(u-1)^{R-1}`; sec. 3b). (iv) The gate
function `gamma*(r)` (sec. 4) is monotone-bounded data whose limit
`gamma*(infinity)` is exactly the threshold of the checkpoint-closure
class; `gamma*(2) in (1/5, 1/4]` exhaustively.

Referee hooks: `parity` mode (i); `dfs`/`depth0` modes (iii) + witness
verifies (ii's hypotheses at small r); the G5-1 proof is 6 lines from
THM-2966. Note (ii) needs no envelope side condition: closure implies the
envelope pointwise (tail bound).

## 7. Relation to prior lanes (consistency audit)

- **laneC2** ("no finite-M infeasibility theorem yet" at these rates) and
  **laneF2** (beam witnesses to M = 20-26; M = 10/11 walls were
  search-order artifacts) are CONSISTENT with, and superseded by, the
  exact closure through M = 15 at `gamma = 1/2`: closure witnesses
  automatically satisfy every intermediate proved envelope (`D_M` equals
  minus the epoch tail, bounded by tail box mass — no envelope check
  needed, sec. 1). LaneF2's observation that beam corridors reset `D = 0`
  exactly at rows 3 and 7 is explained: those runs were *discovering* the
  checkpoint-closure corridor.
- **laneD's greedy freeze** is a policy fact, not a class fact — already
  known from laneF2; here sharpened: the closure class does not even see
  the band freeze (no intermediate constraint exists in the program).
- **THM-2976/T1**: used as the parity-consistency engine (sec. 2); the
  A-clock (all-band-odd at `A_M = 2^r - 1`) imposes NO constraint on
  closure systems (intermediate homogenized states are unconstrained) —
  the A-clock death law of laneF2 is therefore a *policy* phenomenon tied
  to envelope-tracking searches, not to the closure class. The
  anti-checkpoint `m_ac = 20` sits inside epoch `[16,31]`; whether that
  epoch closes at `gamma = 1/2` is exactly the live r = 4 question.
- **Witness anatomy**: the `[8,15]` closure witness saturates the deepest
  rows' budgets (row 14-15 deltas at full caps 21/35/35/21) — the budget
  HEIGHT of the deepest rows is the resource paying for the high-order
  vanishing, as the sec. 6 height-race mechanism predicts. Witness stored:
  `laneG5-witness-epoch8-15-gamma-half.json` (this directory), reproduced
  deterministically by `solve_system(seed=12345)`.

## 7b. Reproduction

```text
python3 amm12592_laneG5_construction_gate_deathstar.py parity
python3 amm12592_laneG5_construction_gate_deathstar.py dfs        # sec. 3
python3 amm12592_laneG5_construction_gate_deathstar.py depth0     # sec. 3/3b
python3 amm12592_laneG5_construction_gate_deathstar.py threshold --spec 2
python3 amm12592_laneG5_construction_gate_deathstar.py threshold --spec 3
python3 amm12592_laneG5_construction_gate_deathstar.py abeldini   # sec. 5
python3 amm12592_laneG5_construction_gate_deathstar.py solve \
    --spec "1,2,0,8,15"                                           # sec. 4 witness
python3 amm12592_laneG5_construction_gate_deathstar.py inspect    # anatomy
```

Every FEASIBLE verdict is accompanied by an exact witness re-verified by
`ClosureSystem.verify` (independent polynomial-sum assert over Z); every
INFEASIBLE verdict is an exhausted class-sequential DFS with sound
capacity pruning (prunes only provably-dead branches). The
`[8,15] gamma=1/2` witness is archived as JSON in this directory.

## 8. Verdict state and next obligations

- Checkpoint-closure program: sound and sharp (G5-1); parity-free at
  dyadic cuts only (sec. 2). This is the construction-side normal form:
  any gamma < 1 fair extractor OF CHECKPOINT-CLOSURE TYPE exists iff the
  gate function's limit allows gamma — and the class is closed under
  everything the certificate could gate.
- Gate curve: `gamma*(2) in (1/5, 1/4]` exhaustively; `gamma*(3) <= 1/2`
  with `INFEASIBLE <= 1/5` exhaustively (boundary refinement + r = 4 runs
  recorded below as they complete).
- Telescoping reading of the certificate integers: NOT SUPPORTED by the
  scans (sec. 5.1); the one structural hit is 1285 as a Lucas
  forced-parity product numerator (sec. 5.2).
- **Next obligations** (for the parent session):
  1. Referee + canonize G5-1/G5-2 (script modes `parity`, `dfs`,
     `depth0`; witnesses re-verifiable in seconds).
  2. Decide the r = 4 epoch `[16,31]` at `gamma = 1/2` — the single most
     informative open bit: feasible extends the exact frontier to M = 31
     and strengthens the construction gate; exhausted-infeasible would be
     the first closure-class obstruction above the trivial rates and the
     right dual target to compare with (27).
  3. If r = 4 closes, attack an INDUCTION: closure of `[2R, 4R-1]` given
     closure machinery at `[R, 2R-1]` (e.g. via the doubling map
     `u -> u^2` on the sec. 6 factorization, which is exactly the
     Frobenius that makes checkpoint parity free). A proof for all r at
     any gamma < 1 gives `C* <= 1 + gamma < 2`, refuting `C* = 2`.
  4. The dual side: formulate the exact lattice extremal problem
     `H_min(R', D)` (sec. 6) and evaluate its two-ray Legendre form
     against `(r_A, r_B, alpha, 1/25)` — the surviving candidate home for
     certificate (27) inside this lane's frame.
