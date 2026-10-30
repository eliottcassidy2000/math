# Session 117: Ramsey structure and unrestricted Boolean circuits

**PROVED scoped results; no new Ramsey-number bound or explicit-NP circuit record.**
Research session of 30 October 2026. The challenge uses all 16 binary Boolean gates,
unrestricted DAG sharing, and actual input-bit length. Recent preprint theorems
and provisional leaderboard claims are distinguished in [SOURCE_MODELS](SOURCE_MODELS.md).

## Results

1. **Exact B2 classification.** For n>=3, an affine Boolean phase plus a single-point
   defect has cost n-1 when its linear support has weight k<=1, and n+k-2 otherwise.
   See [proof and explicitness](circuits/POINT_DEFECT_CIRCUITS.md) and
   [independent raw-DAG audit](audit_circuits/AUDIT.md). Thus the previous extreme-Walsh
   family has at most 2n-2 gates and cannot improve the challenge benchmark.
2. **Complete affine Ramsey geometry.** We classify and count every monochromatic
   affine flat of those functions, with colors and low-dimensional exceptions.
   One GL(n,2) orbit has identical colored-flat geometry and Walsh multiset while
   its exact circuit costs range from n-1 to 2n-2. Consequently an orbit-only lower
   bound on this orbit cannot exceed n-1. [Proof](affine_ramsey/AFFINE_RAMSEY_GATE_GAP.md),
   [independent audit](audit_affine/AUDIT.md).
3. **Uniform tournament family closure and enumeration.** Every twin-containing
   transitive subset of a skew double projects to a transitive base subset, with
   exactly k(k+1)/2 lifts per transitive k-set. This yields the derivative formula
   for the full enumerator and an O(poly(n)3^n) algorithm versus naive 4^n subset
   enumeration. The alternating family has t(D(A_l))=2l+1 for every l>=1, with a
   unique maximum for l>=2 and twelve maxima at l=1. [Proof](ramsey/ramsey_twin_projection.md),
   [independent audit](audit_ramsey/AUDIT.md). Not a best-known-algorithm claim.
4. **Complete local-potential frontier.** For the specified provisional elimination
   template, the feasible coefficient/penalty pairs are exactly
   0<=c<=4 and b>=max(c,2c-2,4c-7). All admissible weights are classified.
   Coefficient four and penalty nine force the normalized weights (1,5,3).
   [Proof](challenge_audit/POTENTIAL_FRONTIER.md), [bounded source audit](challenge_audit/AUDIT.md).
   This is a theorem about the row inequalities, not certification of a four-coefficient
   circuit theorem or an upper bound on circuit complexity.
5. **Classical bilinear counterbenchmark.** Inner product has perfectly flat Walsh
   spectrum, full exponential sign-matrix rank, sharp small monochromatic rectangles,
   and classified maximal affine flats, yet exactly n-1 B2 gates.
   [Proof and checks](bilinear/BILINEAR_BARRIER.md). This prevents another false transfer.

## Inheritance and research direction

[INHERITANCE](INHERITANCE.md) recovers the precise old Ramsey/tournament work and
separates it from historical gate-complexity analogies. [SYNTHESIS](SYNTHESIS.md)
records preserved information, losses, concurrent contributions and the next decisive
structural test. No challenge account, submission, or external proof code was run.
No literature-priority claim is made for these elementary deductions.

## Verification

Run `python3 verify.py` from this folder, or use its repository-relative path.
It runs eight standard-library Python checkers both normally and under -O,
requires identical output, and validates the sealed file hashes. The optional
four-input C++ search is reproduced with `python3 verify.py --cpp` (requires g++).
Both Python modes and the C++ exclusion were run during the session. Finite tests
corroborate the displayed all-parameter proofs; they are not Lean certification,
SAT proof logs, or independent human review.
