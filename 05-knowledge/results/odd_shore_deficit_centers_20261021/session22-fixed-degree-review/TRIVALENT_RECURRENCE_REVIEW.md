# Exact trivalent recurrence: independent recovery audit

Verdict: PASS. This is a fresh algebra verification, not reliance on recovered Session21 computations. The sole new classification input is the separately reviewed bound deg R_6<=47. This audit does not reprove that center bound.

Put p=1-z-z^3, q=p(z²)=1-z²-z^6. The inherited trivalent formula has denominator D=p²q² and numerator

    N = z^5 q²/2 - z^5 p²q/2 + R_6 p²q
        + z^10 p²(1+z+z³).

At a root alpha of p, q(alpha)=2 alpha^4. Since alpha is nonzero, N(alpha)=alpha^5 q(alpha)²/2 is nonzero. Thus N shares no root with p, independently of R_6.

Modulo q one has p(1+z+z³)=-2z^4. Hence at a root beta of q, neither p(beta) nor 1+beta+beta³ vanishes, and beta is nonzero. Consequently N(beta)=beta^10 p(beta)²(1+beta+beta³) is nonzero, again independently of R_6.

Both p and q are squarefree and they are coprime, verified by exact polynomial gcds. Therefore N/D is reduced for every polynomial R_6. All nine distinct denominator roots occur as double poles. Its exact minimal eventual homogeneous constant-coefficient recurrence order is deg D=18. This minimality remains true if any finite prefix is changed, since subtracting a polynomial cannot remove these poles. It is not merely an upper bound or a generic-R_6 assertion.

If deg R_6<=47, then deg(R_6 p²q)<=59. The other three displayed summands have degrees 17, 17, and 19. Thus deg N<=59. For b_n=[z^n]B_3 and D=sum_(j=0)^18 d_j z^j, coefficient comparison gives

    sum_(j=0)^18 d_j b_(n-j)=0 for every n>=60.

The coefficients d_0,...,d_18 are

    1, -2, -1, 2, 1, 2, -4, 2, 0, 0, -1, -4, 3, -2, 3, -2, 2, 0, 1.

The recurrence start n=60 is safe, not asserted earliest. The actual numerator degree and earliest onset can be smaller. For a generic degree-47 R_6, the z^59 coefficient of N is minus the z^47 coefficient of R_6, confirming the bookkeeping but not asserting that this coefficient is nonzero for actual tournaments.

Exact symbolic controls: trivalent_recurrence_check.py uses a generic 48-coefficient R_6, verifies the two quotient-ring identities, both generic numerator residues, three squarefree/coprime gcd conditions, numerator and denominator degrees, top coefficient and recurrence normalization. All 13 checks passed normally and under python -O. Results are in trivalent_recurrence_receipt.json. This separate addition does not modify the frozen main proof or either earlier checker.
