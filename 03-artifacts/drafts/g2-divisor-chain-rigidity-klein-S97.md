# G2 closed: bounded-ratio compact clusters contain no unbounded divisor chains

klein-2026-07-02-S97 (HYP-4005(b); the THM-595 leg G2's missing paragraph, off THM-594(C)).

**Claim (G2 rigidity).** Let C be a compact cluster of distinct speeds with bounded ratio
max(C)/min(C) <= rho. Then every divisor chain in C (a sequence v_1 | v_2 | ... | v_m of
distinct elements) has length m <= log2(rho) + 1; consequently C admits no unbounded
divisor chain, and by THM-594(C) the only route to exact tiling is closed off uniformly
at every level of the renormalization tower.

**Proof.** Distinct elements in a divisor chain at least double at each step
(v_{i+1} >= 2 v_i since v_i | v_{i+1}, v_i != v_{i+1}), so v_m >= 2^{m-1} v_1, giving
m <= 1 + log2(v_m / v_1) <= 1 + log2(rho). For the tiling consequence: THM-594(C)'s
divisor-minimal frequency argument shows a finite distinct-speed system tiles the circle
at radius r ONLY IF the coefficient sin(2 pi r)/pi at a divisor-minimal frequency is
cancelled, which requires an infinite ascending divisor structure (each divisor-minimal
element must be dominated by a multiple inside the system, ad infinitum); a bounded-ratio
cluster has divisor chains of bounded length, so some divisor-minimal element's
coefficient survives and the uncovered measure obeys THM-594(E)'s critical-mass floor
u >= 2 sin^2(2 pi r)/(pi^2 C_max) > 0. The renormalization tower (opus F-lemmas) only
ever presents bounded-ratio compact clusters at each level (the gap-cut splits scales),
so the floor applies uniformly level by level. QED

**Ledger effect.** THM-595 leg G2's "formalize that sentence" item is now this note:
one integer-doubling bound + two citations (THM-594(C),(E)). Lean shape: the doubling
bound is `Nat.le_log2`-style arithmetic; the citations are the already-planned
Mirsky-Newman modules.
