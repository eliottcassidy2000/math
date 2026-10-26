# Independent audit of the unique-six theorem and small aliquot fibres

2026-10-26. **PASS.** No mathematical correction required.

Reviewed `../aliquot/README.md`, SHA256
`b86f1eb6e7f671370bd12fac0253a7ca3e3af7d9fb166d9c2283c1947b676063`,
and `../aliquot/check_aliquot.py`, SHA256
`61324c5a5c1be9b3094006a00372710c41fb9a6aee205d696deb2e54ac6afd06`.
A later status-and-audit-link promotion is covered by this review; substantive
mathematical edits require renewed review.

## Quantifiers, direction and prime guards

Theorem A first uses s(p²)=p+1 only with p prime. The p=2 case gives3,
which is not perfect; every other possible perfect target is even. This
justifies applying Euclid–Euler without an assumption about odd perfect
numbers. Its Mersenne exponent r is prime; after the r=2 exception every
r is odd. The modulo3 calculation gives N−1 divisible by3 and greater
than3, excluding the proposed prime p. The converse25→6→6 is direct.
The classical characterization was independently read in
[King, Theorem3, pages1–2](https://www.jonathanking.org/PDF/mersenne-perfect.pdf).
The proof neither classifies arbitrary square incoming edges nor arbitrary
prime-power incoming edges.

For the inherited fibre bound, n composite and p its least prime divisor
give n=pq with p<=q and the distinct proper divisors1,q. Equality in both
n<=q² and q<=m−1 forces p=q and m=p+1. If n is not a prime square,1,p,q
are distinct; their sum and AM–GM give4n<=(m−1)². Thus the finite source
bound covers every source for targets2..60, including prime-square equality
sources above60. The separate infinite fibre at1 and terminal fibre at0
are handled correctly.

The distinct-prime Goldbach implication has the correct direction and the
correct threshold. I independently read
[Pomerance–Yang, §1](https://math.dartmouth.edu/~carlp/uupaper3.pdf),
which explicitly makes the distinct-prime hypothesis. Equal primes do not
give the semiprime formula. No failed Goldbach search is promoted to an
empty aliquot fibre, and prime-untouchable uniqueness stays open.

## Independent computation

The producer was rerun normally and with `-O`: both match its frozen
19,801-gate receipt byte for byte, SHA256
`9ef780d23709b11196065167dac972ed32d60434fecd086d79c474547a7d3103`.

The separate [audit_aliquot.py](audit_aliquot.py) imports no producer code.
It adds each proper divisor to all its strict multiples, independently
marks composites, and reconstructs every target fibre2..60 from sources
through3481. It also checks both inverse bounds, their equality case, and
the missing-late-source and repeated/composite-prime hostiles. Its16,005
gates pass normally and under `-O`; both independent receipts have SHA256
`df07d0ec96ecb1ff43368075c456bcdeaed835298f77027b4f89ee1cba25128e`.
Reproduce from this directory with
`python3 -B audit_aliquot.py` and `python3 -B -O audit_aliquot.py`.

The finite tests corroborate the implementation. The unique-perfect-target
classification and completeness of the finite fibres follow from the
quantified proofs above, not extrapolation from those controls.
