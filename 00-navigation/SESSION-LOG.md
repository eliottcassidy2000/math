## boxeph-2026-07-21-S206 -- what an LRC(14) disproof MUST be, and why Fibonacci is the foil (HYP-8815)
> **CURRENT-TRUTH WARNING (2026-07-21):** This is chronological provenance,
> not a status authority. Entries may be refuted later in the same file. Start
> with [`START-HERE.md`](START-HERE.md), [`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md),
> and [`../01-canon/ACTIVE-GUARDRAILS.md`](../01-canon/ACTIVE-GUARDRAILS.md).

## boxeph-2026-07-21-S205 -- JC<->LRC = one n=12 AP-rigidity; comprehensive view; Keller counterexample verified; red-team suite (HYP-8810)

**Owner:** mine the repo for connections to 12 (esp. Fibonacci); think about a disproof construction for LRC(14) and refine what it must be.

**MINING (12 + Fibonacci).** The LRC(14) extremal deep well {1..12,182}: the 12 non-max speeds ARE the AP {1..12}; 182=lcm(13,14) blocks both CF neighbors of t*=14/183=[0;13,14] (boxeph-S103). The arithmetic is EISENSTEIN/ANTI-GOLDEN: Φ_6(14)=183, period-6 x²−x+1, ℤ[ω], Heegner −3, LARGE partial quotients (fastest CF). FIBONACCI is the exact sign-flip sibling (golden x²−x−1, [0;1,1,1..], slowest CF, three-gap WORST) — klein-S124: Fibonacci is the FOIL, not a lever; a golden-structured set is the SAFEST (easiest lonely), farthest from a counterexample.

**DISPROOF CHARACTERIZATION (HYP-8815, refined to 7 constraints).** A disproof of LRC(14) = a 13-speed set with covering-min M<1/14. It must be: (1) covering + PRIMITIVE; (2) M<1/14 (below 14/183, in the unproved multi-killer/≤1-far regime); (3) NON-AP (else Wall A/THM-1017 forces M≥14/183); (4) sub-extremal additive triples (THM-730: AP uniquely maximizes); (5) hence tight via HIGHER-ORDER autocorrelation (smaller disc_v, THM-731/732) not triple energy; (6) anti-golden (large CF quotients), NOT Fibonacci; (7) far element blocks the CF neighbors of the new t*.

**THE ONE SENTENCE:** LRC(14) is false iff the additive-triple maximizer (the AP) is NOT also the maximizer of the FULL good-set autocorrelation — iff some non-AP 12-set concentrates higher-order energy more than the AP, buying smaller disc_v (hence smaller M) than 14/183 despite losing triples. Wall A = the JOINT-ORDER extremality of the AP.

**SEARCH (exact covering-min, verified).** deep well = 14/183 (extremal, at t*=14/183); Fibonacci-12 sets M≈0.17-0.24 (LOOSEST — foil confirmed numerically); 2·AP={2..24,182}=7/92 < 14/183 but >1/14 (NON-primitive: M dilation-invariant, reduces to {1..12,91} which omits mult-of-14 => lonely at 1/14 => a disproof MUST be primitive). No covering set beat 1/14 => consistent with LRC(14) true.

**Where to look / not:** NOT Fibonacci/golden (the foil, loosest). A disproof is a PRIMITIVE NEAR-AP anti-golden set trading a little triple energy for more concentrated higher-order autocorrelation — a narrow band just off the cold transitive vertex (my continuum τ~0+). The 12 non-max + 1 far skeleton is fixed.

**Honest:** no disproof (as expected — LRC believed true, LRC(≤13) proven); the value is the SHARPENED target (primitive, near-AP, anti-golden, higher-order-autocorrelation-beating = Wall A joint-order AP-extremality) + Fibonacci ruled out as the foil, numerically. Artifacts: reflection what-an-lrc14-disproof-must-be-and-why-fibonacci-is-the-foil-boxeph-S206.md, HYP-8815, script lrc14_disproof_search_boxeph_S206.py (+.out).

**Coordination:** codex owns the technical Frobenius->LRC transfer (THM-2041/2042) and the JC2 rank-two Poisson descent scaffold; my synthesis (AP-bridge + rank/functional frame + red-team) COMPLEMENTS. Retired my THM-2040 (subsumed by THM-2022 + death-star S91). Pulled/checkpointed 3x.

**Honest scope:** synthesis + a verified counterexample check + the AP-bridge identification (not a proof of JC(2) or LRC(14)). Artifacts: reflection jacobian-and-lonely-runner-two-nullcones-...-boxeph-S205.md, HYP-8810, THM-2040 retired.
