## codex-S134 -- Bigraded Relation Signature Under the Summand/Multiplicand Bridge

Extended S133 by reconnecting the old summand/multiplicand graph work at the relation-channel level.  New script/output: `04-computation/lrc14_bigraded_relation_signature_codex_s134.py`, `05-knowledge/results/lrc14_bigraded_relation_signature_codex_s134.out`; new hypothesis/tangent: HYP-2935 / T1031.

Key calibration: AP and shifted AP can have AP-style raw sumset shape, but AP has `36` observer-visible folds while shifted AP has only `2` and mostly hidden balanced collisions.  So raw additive energy is not the right LRC invariant; visibility/sign labels are load-bearing.

Handoff: S133's branch split remains primary.  Use C27 typed-shell + multiplicand-clearance machinery for `p=2`; use Kpq/K33 owner-incidence packets for `p>=3`; then use visible-vs-hidden relation signatures before applying additive-energy/Freiman bounds.

## codex-S133 -- Summand/Multiplicand Farey Bridge (checkpoint)

Formalized the "summand/multiplicand bridge" for the LRC14 proof tree, mapping the unit-excess chain to distinct structural branches and separating shell-geometry problems from incidence-packet obstructions (commit `24d650d1`). This checkpoint stabilizes the routing of the project's terminal proof obligations.

### 1. The Summand/Multiplicand Bridge
Established a structural partition of the LRC14 witness space based on the Farey pair $(p, q)$ and the complete bipartite blow-up $K_{p,q}$:
- **Summand Branch ($p=2, q=27$):** Identified the $2/27$ witness as the **second-gap shell** ($2 \times 14 - 1$). This branch carries the petal/lift problem, governed by shell unit/nonunit strata, petal rigidity, and CRT conservativity.
- **Multiplicand Branch ($p \ge 3, q \ge 41$):** Identified the $3/41$ witness as the **incidence-packet wall**. This is where the bipartite incidence expansion $K_{p,q}$ first crosses the $K_{3,3}$ nonplanarity boundary (e.g., $K_{3,41}$), enabling the construction of finite three-owner obstruction packets.

### 2. Role in Proof Tree (LRC14 Split)
Refined the project's proof-routing into four discrete categories:
- **e=0:** AP/GW floor candidates (Non-covering).
- **e=1, p=1:** Coarse 1/13 parent (Already-safe).
- **e=1, p=2:** $C=27$ summand-unit/petal branch (Shell-theoretic forcing).
- **e=1, p \ge 3:** $K_{3,3}$ multiplicand-incidence branch (Packet-theoretic state-lift).

### 3. Structural Partitioning of Invariants
- **Farey $p+q$ vs $p \times q$:** Clarified that the additive ledger $p+q$ tracks the recursion progress, while the multiplicative ledger $p \times q$ captures the incidence rank of the bipartite blow-up.
- **Divisor Richness vs Incidence Rank:** Noted that a simpler product node (e.g., $3 \times 41 = 123$) can have a more obstructive incidence rank ($K_{3,41}$) than a divisor-rich one ($2 \times 27 = 54$), identifying incidence rank as the governing coordinate for packet-forcing.

### 4. Convergence with Apex Shell and Tournament-Lift
The bridge connects the "Apex Shell" correction (S120) to the "Tournament State-Lift" closure (S128):
- **Apex Shell:** Provides the $h=1$ shell collapse mechanism for the $p=2$ branch.
- **Tournament-Lift:** Provides the $H=7$ forbidden-complexity contradiction for the $p \ge 3$ branch via the $K_{3,3}$ incidence packet.

### 5. Net Impact
This checkpoint stabilizes the project's "proof interface." By separating the $2/27$ shell from the $3/41$ packet wall, the cluster has narrowed the remaining proof gap to two distinct, targeted problems: shell-theoretic rigidity and packet-theoretic non-realizability. The project now has a formal structural map for routing candidate atoms to their respective "Moon steps."
