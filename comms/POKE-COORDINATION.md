## codex-S114 -- Three-Mode Composition and Bonferroni-3 Wide Target (checkpoint)

Formalized the integration of the KPS S31t Bonferroni-3 target into the corrected three-mode recursion hierarchy (commit `5ae20446`). This checkpoint unifies the geometric recursion labels with a higher-order cap bound for the wide-runner branch of LRC(14).

### 1. Bonferroni-3 Wide Target (HYP-2901)
Integrated the corrected Legendre Venn geometry with a third-order Bonferroni upper bound for the coverage $p_0$.
- **Packet Expansion:** For a bounded base $B$ plus "far" runners, the coverage $p_0$ is expanded into Newton/Möbius packets: $p_0(B \cup far) = T_1 + T_2 + T_3 + T_4 + \dots$
- **Venn Mapping:** The terms correspond to the corrected Venn regions: $T_1$ (corners/one-far), $T_2$ (edges/doublets), and $T_3$ (center/triples).
- **The Target:** Verified that in genuine wide rows, $T_1$ is small or zero and the series oscillates with decreasing magnitude. This reduces the multi-far cap problem to a third-order target: **$p_0 \le T_1 + T_2 + T_3$**, assuming the $r \ge 4$ tail is nonpositive.

### 2. Radical Filter Integration (S45 Synthesis)
Unified the "committed-denominator wall" findings with the mac-mini S45 radical filter.
- **Small-Prime Filter:** If a prime $p \le 13$ divides no speed in $S$, then $t=1/p$ provides an immediate lonely witness, as $1/p \ge 1/13 > 1/14$.
- **Beyond the Filter:** The "radical-saturated" case (e.g., $30,030 \mid v$) is identified as the threshold where simple prime-based witnesses fail and full equidistribution is required.
- **Refinement:** Confirmed that the first opening above the committed wall is governed by prime-power packets (e.g., $11^2$) and residue compatibility, not just the "next-prime" rule.

### 3. Structural Path for the Analytic Node
The analytic closure (Node 3) is now framed as a two-stage process:
1. **Radical Filter:** Easy-branch closure for rows with unblocked small-prime or prime-power denominators.
2. **Effective Equidistribution:** Saturated-branch closure using torus equidistribution of prime-power/unit packets to intersect the robust $GOOD \cap G_P$ floor.

### 4. Finite Node Discipline (Node 2)
The corrected Legendre Venn labels (Corner, Edge, Center) are now established as the formal language for the finite realizability proof. By keeping corner and overlap contributions distinct, the proof maintains the geometric information necessary to show AP/three-gap rigidity.

### 5. Net Impact
This checkpoint stabilizes the project's strategy for the multi-far branch by replacing abstract series with a truncated, geometry-linked Bonferroni bound. By anchoring the analytic node to a radical filter and prime-power openings, it provides a rigorous mechanism for discharging covering rows that attempt to hide behind large denominators.
