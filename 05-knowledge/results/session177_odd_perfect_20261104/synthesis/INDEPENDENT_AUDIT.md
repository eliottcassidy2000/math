# Independent synthesis audit

Status: PASS. Reviewed 2026-11-04.

Reviewed NOTE.md against all four source lanes: mod12, oddperfect (including the LCM closure refinement), networks, and adjacent amicable construction. No mathematical defect or overstated open-problem conclusion found.

Checks:

- Tournament realization remains explicitly induced, not closed; the common-exponent and valuation-lower-bound claims match the proof.
- Euler-aware graph roles, strong-connectivity necessity, weaker odd-order sink filter, and exact simultaneous realization match the support packet. Order-one denominator exceptions are retained.
- LCM closure is correctly stated as a fixed-point criterion, with the empty positive-even-exponent case supplied by the companion proof. Orders 5,3,15 give the stated conjunction forcing 4561.
- Exact cut ratios, coprimality of outgoing/incoming products, residual divisibility, and fourth-power clauses match their source assumptions. Exact components are not silently replaced by exponent lower bounds.
- The amicable factorization preserves distinct primes and coprimality, correctly attributes the classical Euler construction, and distinguishes existence of odd–odd pairs from the unresolved infinitude problem.
- All five relative document links resolve to existing sibling files. No synthesis-specific executable exists.

Replayed source checkers successfully:

- check_odd_order_support.py: 648 realized support/Euler cases and 2,952 direct sigma checks.
- check_lcm_closure.py: 4,032 closure subsets and the conjunction witness.
- verify_networks.py: 500 direct sigma checks, eight unitary-cut controls, eleven conditional clauses.
- verify_amicable.py: ten fixed-kernel certificates, four odd kernels, independent bounded comparisons.
- mod12/check.py: 1,099 labelled tournaments and 5,580 residue checks.

Primary-source URLs were opened and checked: Euler translation arXiv:2409.08783v2, Bibby–Vyncke–Zelinsky arXiv:1908.09420, Zhu–Li Zenodo22005604, Pomerance's dated manuscript, and Pollack arXiv:2610.05661v1. Their cited titles, author identities, dates, and restricted attribution match the synthesis. The Zenodo result is correctly labelled unreviewed and not imported as a bound or nonexistence proof. MIT classical inputs and the Fletcher–Nielsen–Ochem paper were also checked in the source-lane research. This audit does not establish literature novelty or independently audit the external papers' full proofs.

No edits to the synthesis were required. The remaining questions are accurately labelled research targets rather than completed results.
