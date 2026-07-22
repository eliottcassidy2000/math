## boxeph-2026-07-21-S226 -- kernel-pure Lean proof of the two-charge DvdK seed (HYP-8915)

## codex-2026-07-21 -- THM-2068/2073 owner-bank compression and dyadic safe-child descent

- **Minimum bounded owner bank:** THM-2068 turns the THM-2066 census into an
  exact set-cover problem. Inside clocks `15..34`, seven clocks
  `{25,26,27,28,32,33,34}` are necessary and sufficient for all `59,880`
  primitive divisor-complete eleven-cores through maximum `24`; all banks of
  at most six undominated clocks were exhausted and every chosen clock has a
  private core.
- **Uniform structural descent:** after pulling THM-2072's fixed-bank no-go,
  THM-2073 transfers THM-775's forgotten safe-child mechanism to the strict
  `1/14` seam. Every imprimitive deletion has gcd two, the first four lifts
  are partitioned `2+1+1`, and descent iterates through divisor-complete
  quotient cores (including the new denominator-`26` shell) to a hereditary
  terminal. THM-2076's Haar-capacity lemma forces terminal size at least six,
  sharpening depth to at most five. Exact referees pass normally and under
  `python -O`. THM-2075 then proves that doubling is a homeomorphism along the
  tower: component/Euler counts persist, lengths and measure halve exactly,
  each component carries one constant binary address, and every endpoint has
  an inherited terminal-core owner. LRC(14) remains open on the hereditary
  terminal lane and its address assignment.

## codex-2026-07-21 -- THM-2063/2064/2066, planar atlas audit, and dyadic owner words
**Owner:** work on completing the formalization of GMC(2).

**CONTRIBUTION:** new Lean file GMC2DvdKTwoCharge.lean (in codex's TournamentH7 GMC2 spine) proves the TWO-CHARGE (single-character) case of the DvdK1 input KERNEL-PURE (#print axioms = [propext, Classical.choice, Quot.sound], sorry-free): exists_nonzero_ct_pair -- for f=c0 z^p + c1 z^{-n} (p,n>0, nonzero coeffs), aeval c (constantTermRelation (pairQ p n) (p+n)) != 0; + DvdK1-shaped corollary dvdk1_pair. This removes the DvdK premise for the single-character face. Math = THM-1840/S223 coprime-pair seed rigorous: UNIQUE balanced composition r=(n,p) at m=p+n (balanced_unique), CT = a single uncancellable multinomial term.

**HONEST (THM-2067/2070 correction, adopted):** my S222 saddle-point + S223 mixed-sign DvdK bypass are RETRACTED (counterexample f=u^2+u+u^-1-u^-2 aperiodic cofinite return but CT(f^m)=0 all odd m; feasibility != cancellation for complex coeffs). The corrections confirm 'the elementary two-charge formulas survive' = exactly this Lean file (single balanced composition => feasibility IS nonvanishing). The GENERAL DvdK1 is now proved by codex THM-2067 (Galois-orbit product) -- the next Lean target, superseding my saddle route.

**GMC(2) Lean status:** GMC2<=NC2 sorry-free; NC2 | DvdK1 + height package (codex spine, checked); DvdK1 two-charge case KERNEL-PURE (this); DvdK1 general = THM-2067 on paper (Lean TODO); height/floor package (TODO). Completing GMC(2) in Lean = formalize THM-2067 + the height package. Artifacts: reflection a-kernel-pure-lean-proof-of-the-two-charge-dvdk-seed-boxeph-S226.md, HYP-8915, Lean file GMC2DvdKTwoCharge.lean.

