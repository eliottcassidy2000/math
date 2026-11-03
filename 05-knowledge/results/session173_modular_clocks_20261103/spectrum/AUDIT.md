# Independent audit of Session 173 spectrum packet

**Verdict: PASS. The requested scope clarification has been applied and verified.**

Reviewed NOTE.md, check.py, common.py and the saved certificates. Ran normal and optimized main checker; both outputs are byte-identical to certificate.json. Main file hashes match those printed in NOTE.md.

## Mathematical review

- Uniform-module substitution: sound. A cover in labelled slots has a unique maximal-run skeleton, and each module contributes exactly k!p(k) fillings. Every module's size is recovered from its full profile. The argument includes empty slots and prohibits equal adjacent module names.
- Arbitrary-r connectivity: sound under the strict inequalities on both old-side sizes. Any deletion of at most r vertices except precisely the feedback module leaves all three cyclic sides nonempty. Removing the feedback module leaves a nonstrong ordinal sum, proving exact connectivity and uniqueness of its minimum cut.
- Canonical recovery: sound for the unmodified stored graph. The unique minimum cut identifies the feedback set, and the ordered strong components of its complement identify all seed types. This argument does not require treating the originally chosen labels as canonical.
- Hamilton skeleton coefficients: correct. The two labelled feedback vertices give the stated seven coefficient pairs. Each old side has at most three runs.
- Basis transformation and storage/read formulas: correct. Independent exact matrix multiplication verifies the q-to-F transformation; literal subset recursion checks all addressed responses at r=1,2,3,4, and arithmetic checks both closed read formulas for M=2..30.
- Positive arbitrary-r gap: sound. Falling-factorial coefficients remain nonnegative under ordinal multiplication, and the specific all-A/all-L/transitive-feedback skeleton supplies a gap at least 2·9^(M−1). The proof does not extrapolate from the two-feedback formula.
- Arithmetic clocks: exact v3=M and the stated mod8 values follow from the displayed three-exponential expression. No spectrum completeness conclusion is claimed.
- Tested scope: accurate. The main test uses full profiles only at M=2, Hamilton counts at three M=3 instances, and separate resilience tests for r=1,3,4. The arbitrary-r and arbitrary-M results have combinatorial proofs.

## Independent controls

`independent_audit.py` uses permutation descent enumeration rather than the supplied least-vertex partition recurrence to obtain seed/read profiles. It verifies substitution in every one of the eight labelled three-position tournament contexts with a five-vertex seed and two singleton modules. Its separately implemented incoming-subset recursion checks both read addresses at every r=1..4. It also independently verifies the coefficient-basis matrix and closed formulas.

Two hostile controls sharpen scope:

1. Deleting local vertex0 or2 makes B0 andB1 literally identical. Therefore connectivity resilience is not data-erasure robustness. The note now explicitly states in its opening Resilience scope paragraph that the M bits belong to intact storage graphs; a deletion may erase a bit even though the remaining tournament stays strong. This clarification was verified in the final text.
2. At r=5 with two five-vertex old sides, deleting an old side gives another five-element separating cut, alongside the feedback cut. Thus the strict side-size hypotheses cannot simply be replaced by nonstrict ones in the unique-cut claim.

Run `python3 independent_audit.py` and `python3 -O independent_audit.py`; their receipts agree. The independent audit imports no supplied recurrence or wrapper, and reads no repository ledger.
