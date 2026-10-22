# Modular-doubling corollary audit

PASS, subject only to interpreting base 0 as a formal fixed point when discussing cycle entry and period. The actual constraint skeleton omits its loop, so its root component terminates at the base.

Reviewed DOUBLING.md SHA256: eabbb0db8a1ffd34d667f56ac626ad6acb18d08de4dbbdcbffcbc9e7b5ceff8e.

The original repair index j=k-a yields pair {k-a,r-a} and third index k-2a. Ownership gives 2a when 2a<M, zero when 2a=M, and 2a-M otherwise. This proves the claimed doubling conjugacy, including the even-modulus base boundary.

CRT gives periodic residues with zero 2-primary coordinate, and depth max(0,s-v2(a)), with a=0 assigned infinite valuation. Each odd-factor additive order d contributes phi(d)/ord_d(2) cycles. Omitting the zero fixed loop leaves C nonzero cycles and 1+C weak components. Expanding each contracted pair adds one vertex and one edge, preserving cycle rank. Hence the expanded graph has cycle rank C. Its F2 difference constraints, when consistent, have one free bit per connected component and exactly 2^(1+C) solutions. Conversion from auxiliary bits to actual path reversals is a fixed coordinate-wise XOR, so the count is unchanged.

The every-order iff is valid for the actual fixed incidence labels, not merely for some labels: if C=0 the graph is a forest and every order is consistent. If C>0, start with all nonactive vertices before all active vertices and move one cycle spare between its two active endpoints. That operation toggles the selected cycle obstruction and no other cycle obstruction. Thus at least one of those two orders is inconsistent. This establishes the necessity, including the parallel-edge 2-cycle case.

Independent control.py checks parent indexing, CRT counts, divisor formula and depths for all 2<=M<=200, plus the stated cycle counts C(21)=5, C(63)=12 and C(341)=37. It exhausts all path orders and all reversal-bit assignments for M=2,3,4 with fixed arbitrary incidence labels: M2 has 6 compatible orders and 2 solutions each; M3 has 80 compatible and 40 incompatible orders, with 4 solutions on each compatible order; M4 has all 5040 orders compatible with 2 solutions each. No producer code or repository imports are used.

No identification with Collatz, Hamiltonian-path counts, distinct tournament counts, or quotient metagraphs is supported or needed.
