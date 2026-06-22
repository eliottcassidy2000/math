## codex-s119 -- LRC14 Tightness-Star and THM-079 Template (checkpoint)

Formalized the "tightness-star" template for the LRC(14) closure, mapping the final proof steps to the THM-079 template (reducing to a bounded atom and forcing its impossibility) and integrating the apex-denominator lemma (commit `90f8283f`). This checkpoint stabilizes the proof's final "Moon step" strategy.

### 1. The THM-079 Template Analogy
Unified the LRC14 proof order with the tournament complexity ($H=21$) template:
- **Move A (Reduction):** Reduce the problem to a primitive bounded or top-balanced "atom" using the established induction rules (R1/R2/R3 from HYP-2905).
- **Move B (Forcing):** Prove that the irreducible atom has $M > 1/14$ unless it is the non-covering AP/GW tight boundary.

### 2. HYP-2910: The Tightness-Star Theorem
Established the structural target for Move B:
- **Theorem:** A primitive tight row ($M=1/14$) is forced into the denominator-14 AP/Goddyn-Wong boundary and is necessarily non-covering.
- **Apex-Denominator (THM-568):** Proved that any tight optimum at $t=a/D$ forces $14 \mid D$ and $D=14 \cdot gcd(S)$. This makes the structural part of the star target theorem-level.
- **Residual Case:** The remaining challenge is the "14-covering" branch: if $S$ contains multiples of 14, prove $M(S) > 1/14$.

### 3. S119 Audit Findings
- **Tight Locus Consistency:** The AP $\{1, \dots, 13\}$ and GW $\{1, \dots, 11, 13, 24\}$ families share the same six argmax points on the denominator-14 grid: $k/14$ for $k \in \{1, 3, 5, 9, 11, 13\}$.
- **Grid Obstruction:** Verified the exact theorem that for these units $k$, $14 \mid v \cdot k \iff 14 \mid v$. This confirms that 14-free rows always retain their apex witnesses and are therefore non-covering.
- **Covering Window Slack:** Exhaustively verified that for the finite $q$-covering window $[1, 18]$, every row satisfying the necessary $q=2 \dots 14$ covering condition has strict slack ($M \ge 1/12$).

### 4. Localized Moon Step
Refined the final proof obligation to a localized equidistribution problem:
- **Scenario:** $S = R \cup M14$, where $R$ is a 14-free core with $\le 6$ speeds and $M14$ are $\ge 7$ multiples of 14.
- **Target:** Prove that the danger combs from the multiples of 14 cannot cover the $1/13$-margin lonely interval of the core $R$. This bridges the structural gap to the tournament $K_3$ conflict-atom impossibility (HYP-2908).

### 5. Net Impact
This synthesis stabilizes the "atom-forcing" phase of the proof. By anchoring the tight locus at denominator 14 and verifying the slack in finite covering windows, the project focuses the remaining effort on the second-moment/equidistribution branch for large collections of 14-multiples. The cluster is now synchronized on the three-part "terminal theorem" covering compression, slack, and forbidden packets.
