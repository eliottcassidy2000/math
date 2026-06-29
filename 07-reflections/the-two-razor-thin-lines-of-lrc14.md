# The two razor-thin lines: what a disproof of LRC(14) would require

*mac-mini-2026-06-29-S15. Understanding the disproof to sharpen the proof, with an Explore-agent map of the floor machinery (THM-523/580, HYP-3129/3415, OPEN-Q-108) and fresh boundary computation. This reframes — and honestly re-places my own recent work.*

## What a disproof is, exactly

A disproof of LRC(14) is a **primitive covering 13-set** `S` — 13 distinct positive integers, gcd 1, containing a multiple of every `q ∈ {2,…,14}` — with `M(S) = max_t min_{v∈S} ‖vt‖ < 1/14` strictly. Equivalently, the lonely set `L(S) = {t : min_v ‖vt‖ ≥ 1/14}` is **completely empty** (the danger combs, each of measure exactly `1/7`, cover all of `[0,1)`). The covering reduction (THM-523) is why "covering": any *non*-covering set omits all multiples of some `q ≤ 14`, and then `t = 1/q` is lonely — so non-covering sets are closed by an elementary witness.

## The locus I had been staring at is off the path

Last sessions I dwelt on the extremal `{1,…,13}`: `M = 1/14` exactly, lonely set the six units mod 14 `{1,3,5,9,11,13}` = three antipodal pairs `{1,13},{3,11},{5,9}` (= the saddle index 3), each touched by a multiplicative-inverse runner pair, and `(ℤ/14)^* ≅ (ℤ/7)^*` the apex-7 cyclotomic group. All true and beautiful — and **`{1,…,13}` is non-covering** (no multiple of 14), so it is closed by the `t=1/14` witness and **off the critical path**. The apex-7/units structure is the *tight locus of the whole problem*, not the locus where a disproof could live. Worth knowing precisely because it is *not* where the danger is.

## Two razor-thin lines, at different scales

The disproof must be a covering set, and there the boundary splits into two lines that are thin in different coordinates:

**Line 1 — the gap `M`.** Covering sets cluster *far* above `1/14`. I verified the tightest known, `{1,…,11,13,84}` (`84 = lcm(12,14)`), at `M = 7/89 = 0.07865` exactly — a **+10.1%** margin. Scanning the single-large family `{1..11,13,L}`, `L=84` is the minimum; larger `L` only loosens (`168→14/173`, `252→13/159`, …). The gap is *quantized*: at the optimizing `t^*` two runners bind and `M = (v_a q - v_b p)/(v_a ± v_b)`, a fraction whose denominator is a binding-pair sum/difference, so the achievable gaps form a discrete lattice that the covering constraint pushes well clear of `1/14`. **In gap-value the disproof is nowhere near.**

**Line 2 — the measure floor `R′`.** `meas(L(S)) = R′ · (m_R m_Q)`, `R′ = 1 + SPEC/baseline`, `SPEC = Σ_{n≠0} \hat c(n)\overline{\hat g(n)}`. `R′ > 0 ⟹ meas(L) > 0 ⟹ M ≥ 1/14`. The certificate `R′ ≥ 0.642` is proved for **bounded** covering families (speeds `≤ 84`, `r = |Q| ≤ 6`) — a 60% margin — but the **unbounded** case is open. This is the thin line, and the single gatekeeper (OPEN-Q-108): prove `ρ_j ≥ c > 0` uniformly, where the 2-adic descent (THM-580) peels even speeds level by level into 2-sheet per-level floors.

So the disproof cannot hide in the gap (covering sets are 10% safe) — it can only hide in the **unbounded measure floor**, in a covering set whose even-speed structure drives some descended level "lonely-poor" enough that the per-level Cauchy–Schwarz bound goes non-positive. That is a short, sharp list of necessary conditions for a counterexample: primitive, covering, unbounded, even-heavy, with a resonant descended level.

## Where my own work actually sits (the honest placement)

The binding obstruction is **2-adic** (even speeds pull the witness off `t=1/2`, where they vanish), and the floor is **σ-even** (existence of a lonely time is invariant under `t ↦ 1-t`; the measure can be bounded by SOS without Borsuk–Ulam). My recent thread — the metagraph antipodal spectrum, the Ky-Fan oriented-path parity, the signed cycle index, the perfect-number/octonion apex — lives on the **σ-odd witness side** (Rédei, the odd index, klein's R-odd Betti numbers). That side is **orthogonal to the floor gatekeeper**: it certifies that a lonely set, *once known nonempty, has odd-indexed structure*, but it does not by itself produce the floor `R′ > 0`. Saying this plainly is the point of the exercise — the beautiful tournament structure is not on the critical path; the critical path is the even-category measure inequality `ρ_j ≥ c`, and it belongs to the descent/SOS machinery.

The one place my apex-7 work *does* touch the path: the 2-adic descent bottoms out, after peeling all the even speeds, on an **all-odd face**, and that face is the apex-7 cyclotomic object where the gentlest-cyclotomic Heegner SOS (`F_7`, `ℚ(√−7)` class number 1) is the natural certificate for the base-level `ρ`. The three pillars (S14) are not three independent attacks — they are the three stages of *one* descent: **2-adic descent (Mersenne) peels → per-level SOS (Heegner) certifies → Borsuk–Ulam (3 mod 4) is the witness backstop at the all-odd base.**

## Reorganized proof targets

1. **Primary (the gatekeeper).** Prove `ρ_j ≥ c > 0` uniformly by applying HYP-3129's exact-low + Parseval-tail estimate *per 2-sheet descent level*. The agent's read, which I share: this is mechanically *simpler* than the monolithic 14-sheet version (2 sheets, smaller sets), so it is the lowest-hanging rigorous close.
2. **Exploit the 10% gap as slack.** Because covering sets sit at `M ≥ 7/89`, the floor need not be sharp — any `R′ > 0` suffices, and the gap margin can absorb loss in the Cauchy–Schwarz step. Target `R′ > 0`, not `R′ ≥ 0.642`.
3. **The fragile point has a silver lining.** The dangerous descended sets are the *lonely-poor* ones — i.e. the highly **structured/resonant** ones. But structure is exactly where a cyclotomic SOS is *strongest* (resonances live on a sublattice where `F_7`-type positivity is sharpest). So the worst case for Cauchy–Schwarz is the best case for SOS; the two tools should be deployed on complementary strata, not uniformly.
4. **The gap-quantization route (alternative).** Prove `M ≥ 1/14` directly: bound the binding-pair denominator `D = v_a ± v_b` for covering sets so that `M = j/D` with `14j ≥ D` forced. This converts the analytic floor into a number-theoretic statement about speed differences (= the tournament arcs), where the tools of the recent thread (Rédei/OCF on the binding graph) might finally bear on the critical path rather than beside it.

## The one-line summary

A disproof must be a primitive, covering, **unbounded**, even-heavy 13-set with a resonant descended level — it cannot hide in the gap (10% safe) but only in the un-certified tail of the measure floor. So the whole proof reduces to a single even-category inequality `ρ_j ≥ c` on 2-sheet descents, the gap margin gives slack, and the apex-7/Heegner SOS is the base-case certificate at the bottom of the 2-adic descent. Everything σ-odd — Rédei, Ky Fan, the metagraph — is the witness scaffolding around that inequality, not the inequality itself.

See [[the-dihedral-recursion-existence-is-even-witness-is-odd]] (THM-581, the even/odd split), [[two-order-two-structures-parity-and-descent]] (THM-580), [[twentyeight-the-octonion-apex-and-the-three-pillars]] (HYP-3547, the descent as one chain). New: HYP-3548. The Explore-agent map underlies the floor-machinery citations.
