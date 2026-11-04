# Independent audit of the Session 178 synthesis

## Verdict

**PASS after correction.** Mathematical synthesis agrees with the four sibling packets. Both precision corrections below were applied by the author and rechecked:

1. In section 4, explicitly require the observation pi:V→Y to be surjective, or set Y=pi(V). This is an actual hypothesis of the source theorem and is needed for the displayed |V|−|Y| rank formula. Unused observation labels cannot be counted as fibres.
2. In next target 5, replace “length-15 path” by “15-vertex path.” The inherited path has 15 vertices and 14 edges.

The residue counterexample is in the context of doubling; explicitly saying “false doubling pair” would remove a possible ambiguity. This is optional clarification, not a mathematical blocker.

## Scope and consistency checks

- Detector: exact spiral convention, eight neighboring arrows, shell-six completeness argument, ordinary rather than shortcut Collatz map, and nonlocal-site-conjugacy scope agree with `../detector/detector.md`.
- Local verifier: fixed carry alphabets, carry and branch consistency, O(log max(n,m)) strip size/reading time, endpoint-data binding, and communication costs are retained. No bounded-radius geometric oracle or one-cell unbounded-integer storage is claimed.
- Honeycomb: dyadic chart occupies only the specified even-row half-strip; odd steps are jumps; finite wires preserve marked arrows after subdivision, not native one-edge time. Wire counts/lengths and lack of nested layouts agree with the independently audited construction.
- Observer theorem: edge identities, weak connectedness, signed integral chains, additive-linear rank restriction, disconnected codomain correction, and trajectory-order caveats are preserved. The corrected text now explicitly includes surjectivity.
- Prime flow: LTE denominator exception, exponent box, both numerical certificates, strict cut directions, common-branch bounds, exterior-prime checking, and Euler 2-adic balance agree with the prime-flow proof. The result remains a fixed-support finite decision procedure; no global support bound or odd-perfect-number conclusion is asserted.
- Arithmetic zero holonomy is correctly separated from universal vertex-boundary closure, including N=6 and the 2→4→8 example.
- All five next targets are explicitly proposals rather than completed research. They respect current missing interfaces and do not reuse inherited results as new theorems.

## Inheritance and links

Independently fetched and checked the exact linked repository files:

- Session 140 `collatz_hex/PROOF.md`: unconditional global honeycomb subdivision with cited infinite-grid-minor input; exactly five ordinary steps for binary predecessors surviving prime-to-three; exponential inverse growth and unbounded stretch in every component. The synthesis preserves those hypotheses and limitations.
- Session 140 `principles/PRINCIPLES.md`: triangular-site conjugacy and first-positive-return suspension were already present and are properly credited.
- Session 167 `honeycomb/NOTE.md`: the 704/694 count, arithmetic deck subgroup, cube kernel 2Z², primitive two-lap synchronization, and extra square-target chord 1–3 all match the synthesis.
- Session 177 `synthesis/NOTE.md`: order graphs, LCM closure, support connectivity, exact cuts, and induced tournament realizability are properly inherited rather than treated as closed perfect-number constructions.

All four relative sibling links resolve locally. All four external GitHub source links returned their referenced contents successfully. No forbidden path was read or scanned, no installs were performed, and no historical-priority claim was inferred.

This is a synthesis consistency/source audit, not a new experiment or an independent proof of the cited infinite-grid theorem.
