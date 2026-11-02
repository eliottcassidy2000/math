# Weighted stability, arithmetic clocks and Hamiltonian witnesses

Research date: 2026-11-02. Starting revision:
`291b4d9757842c02aa451beb65acfade2848609b`.

**PROVED:** the sharp weighted tournament determinant-gap coefficient is
now established for every n>=466, including its fair-pair equality cases.
This closes 1,034 of the formerly open orders 10..1499. Orders 10..465
remain open; the previously proved orders 3..9 are inherited.

The [proof](STABILITY_466.md) retains the expected binary gap and allocates
disjoint interaction contributions once. The standard-library
[checker](verify_stability.py) passes 3,475 exact gates, including native
Hamiltonian/determinant calculations and a polynomial identity check.
The infinite range follows from the proof, not finite sampling.

The ongoing portfolio is: weighted stability (anchor); priority-preserving
FRACTRAN clocks (niche); marked Hamiltonian witnesses and old odd-cycle
carriers (wildcard). All are compared on a five-concept board: native
observable, retained baseline, overlap, priority, and lost quotient data.

See [HANDOFF.md](HANDOFF.md) for the current continuation boundary. This
packet makes no claim to solve Collatz, RH, full Hamiltonian-value spectrum
classification, or the remaining weighted orders.
