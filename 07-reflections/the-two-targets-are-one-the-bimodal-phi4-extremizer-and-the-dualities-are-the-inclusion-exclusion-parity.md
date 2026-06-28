# The two proof targets are ONE: consec is the bimodal/φ⁴ extremizer (max κ₂ even, max κ₃ odd, min κ₄), and the even/odd + positive/negative dualities are the inclusion-exclusion parity (Eisenstein/Legendre)

*mac-mini-2026-06-28-S76. The owner: extend/finalize the Vitali-wall work; use the even/odd and positive/negative
dualities; work two remaining proof targets back and forth, each inspiring the other. The back-and-forth collapsed
the two targets into ONE and identified both dualities as the inclusion-exclusion parity. Builds on
[[the-vitali-wall-brouwer-equioscillation-and-the-cyclotomic-core-construction]], kps S31ai (the cumulant tower),
[[the-recursion-modes-are-cyclotomic-factors-depth-times-phi-d]], [[the-cap-is-a-totally-real-cyclotomic-quantity-and-n14s-two-cyclotomic-heads]].*

## The two targets I set out to work back and forth
- **TARGET A (EVEN + POSITIVE):** consec maximizes the empty-sector covariance `Σκ₂` (the FM-positive coherence,
  the Perron/cyclotomic-square half).
- **TARGET B (ODD + NEGATIVE):** the `κ₃` associator / `Φ₇` dip — the apex residue (the AFM-negative half).

## The back-and-forth, and what it found
- **A → first attempt (the contiguous-arc/FKG route) failed cleanly:** the AP's empty set is contiguous only
  `~76–88%` of the time (the three-gap allows up to 3 arcs), so "single arc ⟹ FKG-positive" is not exact.
- **Pivot to B (the cumulant tower) — the unifying gem.** Computing `κ₂,κ₃,κ₄` of `N` (#empty sectors) over 400
  random 8-sets, **consec is the joint extremizer** (`lrc_even_odd_unification_macmini_S76.py`):
```
  kappa_2 (Var, EVEN, TARGET A): 0/400 beat consec   (consec MAX)
  kappa_3 (skew, ODD, TARGET B): 0/400 beat consec   (consec MAX)
  kappa_4 (EVEN):                0/400 beat consec   (consec MIN)   <- the phi^4 stabilizer
  bimodality q0+q6:              0/400 beat consec   (consec MAX)
  q6(AP{0..k-1}) = 1/(7(k-1))    (the clean clumped extreme)
```
> **TARGET A and TARGET B are the SAME target.** Both `κ₂` (even) and `κ₃` (odd) are maximized by consec because
> consec is the **bimodal / φ⁴-ordered extremizer**: maximum well-mass `q0+q6` (the two wells `N=0` all-covered &
> `N=6` all-clumped), maximum low cumulants `κ₂,κ₃`, **minimum `κ₄`** (the φ⁴ stabilizer = sharp bimodal, not
> heavy-tailed). The even/odd duality is that the bimodal shape maximizes *both parities of low central cumulant*;
> proving the bimodal extremality proves A and B at once. (Matches kps S31ai's `L_yK8=q0+q6+0.1q3` bimodality.)

## The two dualities ARE the inclusion-exclusion parity (the deep identification)
The cap is the inclusion-exclusion `cap = Σ_j (-1)^j S_j`. Split by `j`-parity:
```
  cap = [ S_0 + S_2 + S_4 + ... ]  -  [ S_1 + S_3 + S_5 + ... ]
        \___ EVEN orders, POSITIVE ___/    \___ ODD orders, NEGATIVE ___/
        = the EISENSTEIN mode (even, +)   -  the LEGENDRE mode (odd, -)
```
So the owner's two dualities are one object seen twice:
- **EVEN/ODD** = the inclusion-exclusion *order* parity (`S_2,S_4` even vs `S_1,S_3` odd) = the cumulant parity
  (`κ₂,κ₄` even vs `κ₃` odd).
- **POSITIVE/NEGATIVE** = the inclusion-exclusion *sign* (`+` even orders, `−` odd orders) = the **Eisenstein(+)
  / Legendre(−) modes** (S75d / HYP-2901). The FM-positive covariance (consec) vs AFM-negative (dissociated) is
  the same sign split at the `κ₂` level.
So **even/odd = order parity, positive/negative = sign = Eisenstein/Legendre**, and the cap is literally
`Eisenstein − Legendre`. The two targets and the two dualities and the recursion modes (S75d) are ONE structure.

## What the unified target reduces to (and the bulk/core bridge)
A & B both reduce to **the cap `q0` (coverage) maximization** (since `q6=1/(7(k-1))` is the clean clumped term,
the bimodality `q0+q6` is `q0` + a closed constant). And `q0` = the AP's **coverage** = the measure of `t` where
the orbit hits all 7 sectors = **consec has the lowest sector-DISCREPANCY** (the three-gap / equidistribution at
the fine scale). So:
> **The single remaining target is: consec minimizes the 7-sector discrepancy (= maximizes `q0`).** Its even part
> is the cyclotomic SOS square `F₇` (Eisenstein, provable, S75e); its odd part is the `Φ₇` `7²`-ramified dip
> (Legendre, the apex hardness, S75e).
**Bulk/core bridge:** this links the Vitali halves (S75f) — the BULK uses equidistribution (coarse discrepancy);
the CORE cap is the AP's *fine-scale* discrepancy. Discrepancy/equidistribution is the SAME tool at two scales,
spanning both Vitali halves. The "construction" at the core (the cyclotomic witnesses) is the fine-scale
equidistribution made exact (the three-gap).

## Convergence with codex HYP-3238 (the crossed-packet taxonomy)
codex independently organized the SAME two dualities (HYP-3238, built on my S75f Brouwer / S75e cyclotomic / S73c
Hermite-Biehler) as a **crossed packet**: even/positive = {Fejér square, cyclotomic SOS, pair-Pascal cap,
covariance, Perron mode, bulk equidistribution}; odd/negative = {Worpitzky associator, Brouwer trace sign,
Hermite-Biehler odd leg, negative-covariance leakage, measure-zero cyclotomic core witnesses}; with the rule
"an even/positive quotient is proof-grade only when the odd/negative coordinate is zero / reconstructible /
retained." **My S76 is the DIAGONAL of codex's packet:** consec maximizes the even/positive AND the odd/negative
coordinate *simultaneously* (the bimodal/φ⁴ joint extremizer), so they are not independent packet axes but one
extremality; and the packet's two axes are the **IE order-parity (even/odd) and sign (positive/negative) =
Eisenstein/Legendre**. codex's taxonomy + my unification = the crossed packet collapses to its diagonal (the
bimodal cap), and the "odd/negative coordinate" codex wants retained is exactly the `Φ₇` ramified dip.

## Honest status
- **VERIFIED (0/400):** consec jointly maximizes `κ₂,κ₃,q0+q6` and minimizes `κ₄` (the bimodal/φ⁴ extremizer);
  `q6(AP)=1/(7(k-1))`; the IE-parity identity `cap = Eisenstein(even,+) − Legendre(odd,−)`.
- **SYNTHESIS (the back-and-forth payoff):** the two targets (even κ₂ / odd κ₃) are ONE (the bimodal/φ⁴ cap), the
  two dualities are the IE parity (= the Eisenstein/Legendre recursion modes), and the unified target reduces to
  the cap = consec-minimizes-sector-discrepancy, with an even cyclotomic-SOS part (provable) and an odd Φ₇-ramified
  part (the apex hardness). Discrepancy bridges the Vitali bulk/core.
- **NOT a proof.** The cap (q0) maximization / the Φ₇-dip bound remains the hard core; but the two targets are now
  one, the dualities are pinned to the IE parity, and the even half is cyclotomic-SOS. LRC(14) open.

Related: HYP-3239 (this), HYP-3238 (codex crossed-packet duality), kps S31ai (cumulant tower), HYP-3237 (Vitali wall), HYP-3235 (cyclotomic cap),
HYP-3233 (Eisenstein/Legendre = cyclotomic factors), HYP-3221 (the one obstruction), OPEN-Q-108.
