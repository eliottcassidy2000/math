# Independent degree-local enumeration audit

Status: mathematical kernel APPROVED; finite-exact controls PASS. Scope: `local_scale.py` SHA-256 b76167ef82dda0b134681abce33e159d5dd80d15e92edb9e31ab9b490f6939ad. No source or repository edits were made by this reviewer.

## Exact mathematical justification

Write O(d)=odd(d!), and store O(d) times the rational aggregate of histories at degree d. A transition with actual multiplicity k at odd cycle length ell, and substitution excess e, reaches q=d+ell(k+e). For k=0 only e=0 is allowed. The quantity

q! / (d! ell^k k!)

is integral: choose ell*k of q labels, put them into an unordered collection of k directed ell-cycles, and order the remaining q-ell*k labels modulo permutations of d distinguished positions. Equivalently it factors as binomial(q,ell*k) times (ell*k)!/(ell^k*k!) times (q-ell*k)!/d!. Its odd part is precisely O(q)/(O(d) ell^k O(k)), because ell is odd. Hence the transfer precomputation's checked divisions are always exact, independently of response or signs. Multiplication by the integer substitution coefficient (-1)^e binomial(k+e-1,e) preserves integrality.

The shift is nonnegative: v2(k!) <= k(k-1)/2, with k=0,1 handled directly; ell>=1 and response>=0. The rational state invariant follows by induction, and normalization commutes with merging because merged histories share degree. Final division must use O(d), not O(bound), as implemented.

## Future-response and signed typing

Actual cycle counts, not substitution excess, contribute gcd response. At stage ell a future response coordinate j is sum of k_L*gcd(L,j) over already processed actual cycles. This is at most total actual used degree, which is at most used substitution degree d <= bound. Width max(1,bit_length(bound)) therefore prevents carries in response updates. A coordinate j above remaining degree bound-d is unaffordable forever because degree never decreases; after processing ell, lengths above ell-2 are likewise unavailable. The mask preserves exactly the future-relevant coordinates. Stage, degree, bound and mode determine the retained length, so omitted leading zero digits cause no state alias. Scaling by a nonzero integer depending only on degree preserves cancellations, zero pruning, state graph and transitions.

Core coefficients and intermediate weights are signed cycle-index aggregates, not counts of objects at intermediate stages. The no-fixed observer remains a rational mass: coefficient 3 in core/no-fixed mode is 2/3. The all-fixed-point history proves O(d) is minimal for unrestricted raw histories, but not for restricted observers or canceled aggregates. The kernel does not discover a smaller future-response quotient or recover class representatives, automorphisms or adjacency.

## Independent tests

Run `python session22-enumeration-review/check.py` and `python -O session22-enumeration-review/check.py` from the workspace root. Both produce identical 35,704-check PASS receipts. No checks depend on assert. The checker does not import the producer's controls.

- All bounds 0..20, all four core/no-fixed modes: independent unmerged odd partitions, exact rational weights and direct polynomial substitution.
- Same range: frozen incoming baseline coefficients and stage-by-stage state/transition counts.
- Literal permutation action on directed edge orbits, with direct homogeneous-pair rejection, for raw and pair-free Burnside coefficients through order 6.
- 17,712 (d,ell,k,e) integral transfers and nonnegative shifts, including ell=1 and k=0,1.
- Literal integer bounds, Boolean flags, rational observer hostile.

These finite tests supplement the proof, not replace it. The small Burnside checker enumerates the declared odd-cycle types; even cycle types have no invariant tournament because an even cycle power reverses an edge. No external published sequence or library oracle is needed for these controls.

## Benchmark assessment

The producer's receipt is an honest alternating three-repeat end-to-end comparison of two actual functions with full-vector equality checked. It is not a pure normalization ablation: the baseline collects extra per-transition coordinate telemetry and per-state tuple/bit telemetry that the candidate omits. Consequently an assertion that every observed speedup comes solely from the new arithmetic would overstate evidence. The reviewer supplies a stripped-telemetry sensitivity control in `benchmark.py`, preserving the baseline mathematical kernel, with local source and hash receipts. Remaining algorithmic differences include precomputation, loop organization and integer operand sizes; these are part of the new implementation, not separately isolated causal effects. Shared-host timing is not universal performance or an optimum claim; memory is unmeasured.

Matched-telemetry samples (median ratios, three alternating pairs): raw n80 1.833x, n100 1.971x, n150 2.158x, n200 1.855x; core n100 1.326x. Exact full coefficient vectors and per-stage mathematical state/transition counts agree at all five settings. n80/150 state checks were integrated into the timed-pair driver; n100/200 state checks were subsequently rerun separately (their original timed-pair driver checked coefficients only, and is now augmented with state checks for reproduction). The normalized baseline has SHA-256 1573f6915057357ad8b57efcd6e8369547384ad6a745dd8aab3bdc94ba1c042a; `telemetry.patch` records the complete difference. The auxiliary zero-valued telemetry fields in that control are intentionally placeholders and must not be interpreted as genuine telemetry.

The producer's division-count upper bounds also pass review: at a degree/stage with m=floor((N-d)/ell), raw has m+1 entries and core has 1+m(m+1)/2. Summing over all degrees and odd stages gives O(N² log(N+2)) and O(N³), respectively. These do not bound response-state transition work or bit complexity.

Conclusion: approve mathematical identity and implementation at the frozen hash; accept the finite tests and properly qualified performance measurements. No source change is required. Keep the telemetry caveat and signed/rational typing explicit when publishing.
