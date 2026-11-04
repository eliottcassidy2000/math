# Independent prime-flow audit

Verdict: PASS after checker robustness correction, 2026-11-04.

Read the complete PROOF.md and checker independently. The exact LTE law retains the order-one denominator exception. For d>1, r is prime to both d and ell=lcm(b,d), so the logarithmic capacities are valid. The exponent bound follows from the maximum column balance; the stated B condition implies B>s-1 and strict monotonicity of x-C-(s-1)log_3(x+1) for x>=B. The finite exponent-box conclusion and support-fixed exact decision procedure are sound. They do not bound the varying prime support.

The cut identity correctly assigns the single factor of 2 to the Euler source. Proper unitary blocks have abundancy strictly between 1 and 2, giving both claimed strict cut intervals and all four rejection inequalities. Strong connectivity is used only as a necessary condition. LCM conjunction pruning retains exact source exponent constraints and does not assume arbitrary independently selectable target edges.

Independent controls in independent.py import no producer code and compute sigma by direct summation. Both normal and -O runs give byte-identical JSON:
- 82,248 LTE comparisons for distinct odd primes below 100 and exponents 1..149.
- Both maximal graphs on S={3,5,7,11,13} are strongly connected.
- Recomputed C,B=(3,16) for q=5 and (4,32) for q=13.
- Direct maxima over every role-compatible exponent in each box reproduce the column capacities. For q=5, the ordinary exponent at 13 has an empty surviving domain. For q=13, surviving upper bounds are (4,2,4,4,1), and every possible component at source 5 has no outgoing prime in S.
- Checked the exponent-bound inequality by exact integer powers for B<=A<=512; the all-height justification is the separate analytic proof.
- Twelve strict cut controls from six even perfect numbers, using reduced exact fractions. These corroborate the algebra only; no odd perfect positive control is claimed.

Producer replay independently performed after its correction: normal and -O stdout agree byte-for-byte; regenerated RESULTS.json agrees byte-for-byte. Producer totals: 34,200 LTE checks, 58,482 capacity checks, 600 closure checks, eight cut controls.

Correction requested and applied: original checker used assert, which optimization strips. Author replaced all assertions with explicit require(...) guards; final inspected code has no assert statements. There was no mathematical correction.

Audited hashes:
- PROOF.md: da41c74502773a9d9d1af8762b64380d842c8cd6cc5e506c1bed9fd6dab5aa6a
- check_primeflow.py: 1a068600f8307892de8f4a92b16253d18bed7cdeefaebb30c806780090c45b22
- RESULTS.json: 95281d7e26515002470fc855aaf0dcb6a5b800c2d0f0ee781000249481d6a91a

Final added corollary independently checked: a proper sink has O_T=1; the two cut cases force respectively I_T<1 or 1<I_T<2, both impossible for a positive integer. This correctly recovers strong connectivity. Final proof hash above includes this paragraph.
