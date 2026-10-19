# Correction proposal for THM-1870

Apply an unmistakable correction banner immediately after frontmatter and replace its misleading title/status, preserving its historical body below a historical-boundary marker.

Suggested title: `Cycles of lengths three through five are single traces; six-cycles are nonspectral at every ambient order at least six`.

Suggested status: `CORRECTED / PROVED: c_k=tr(A^k)/k for k=3,4,5 in every loopless digraph without digons. FINITE-EXACT: six-vertex tournament masks73/83 have equal characteristic polynomial and c6=0/1. PROVED: adjoining universal sources propagates this discrepancy to every ambient order n>=6. REFUTED: the historical conjecture that all sub-Hamiltonian cycle counts c_k, k<n, are spectral. No complexity claim is inferred.`

Suggested banner:

> **Correction, 2026-10-19.** The all-order Hamiltonian-length extrapolation below is false. Co-spectral six-vertex tournament masks73 and83 have `c6=0` and `1`; adjoining a universal source preserves every cycle count and multiplies the characteristic polynomial by `x`, giving co-spectral seven-vertex tournaments with different sub-Hamiltonian `c6`. This extends to every `n>=6`. Conversely, `c3,c4,c5` are exactly `tr(A^k)/k` at every order: a repeated vertex in a closed walk of length at most five would split it into two closed walks of lengths at least three. The old “NOT single traces” claim and the complexity inference are withdrawn. See the squarefree-cycle bridge result for proofs and exact witnesses. The finite six-vertex census observation survives.

Historical-body warning:

`## Historical text below: retain for provenance; the correction above controls all general claims`

Suggested mistake entry, with number assigned centrally:

- **Claim:** ambient Hamiltonian length is the first nonspectral cycle length; sub-Hamiltonian `c_k` is spectral for all orders.
- **Witness:** lexicographic tournament masks73/83, order6, common `charpoly=x6-4x3-3x2-2x`, `(c3,c4,c5)=(4,3,2)`, `c6=0/1`; adjoin a universal source to reach order7 while retaining the `c6` split.
- **First invalid implication:** finite coincidence of `k=n=6` was promoted to a universal dependence on ambient size. A target-preserving padding operation changes ambient size without restoring a forgotten invariant.
- **Repair:** `c3,c4,c5` are single traces at all orders; `c6` is nonspectral for all orders `>=6`. No claim about every `k>=7` follows.
- **Additional correction:** a co-spectral split is an information-loss witness, not a proof of a complexity-class distinction.
- **Mechanism:** characteristic polynomial's cycle-cover term `-c6+disj_33` cancels changes in individual six-cycles against disjoint triangle pairs. Marked square-zero trace-log restores the cycle inventory.
- **Reproducer:** `verify_squarefree_bridges.py` and independent `review/check_spectral_hostile.py`.

Also consider a short warning at THM-1926's `THM-1870` references: its closed-walk zeta facts survive, but the cited ambient-length spectral extrapolation is corrected. No change to the valid zeta identity is needed.
