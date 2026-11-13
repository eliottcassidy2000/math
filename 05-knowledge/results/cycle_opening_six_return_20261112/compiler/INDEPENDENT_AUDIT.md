# Independent audit of the six-coordinate compiler result

**Status: SOUND within the stated scope.** Reviewer: six/dynamics lane,
2026-11-12. Audited `SIX_STATE_MINIMALITY.md`, `six_state_audit.py`, the inherited
`funnel_six_state.py`, and the definitions of ordered path-cover bins and
literal funnel/birth graphs. The new infinite statements have been read as
proofs, not inferred from their finite controls.

The linear lower bound is valid because the genuine past/future observation
table factors through any proposed linear hidden state. The affine lower bound
is also valid: subtracting the seed's future response removes all observation
constants, and the difference of two states evolves under the linear part of
an affine transition. The determinant811008 therefore forces six affine
coordinates, rather than merely five coordinates plus a homogenizing one.
The scope remains rational linear/affine realizations on all continuations of
the stated reachable family. It gives no lower bound on nonlinear encodings.

The mixed E/B periodic obstruction is valid for arbitrary nonempty W containing
E. The displayed conjugacy gives operator column norms at most2 and3, so the
upper eigenvalue modulus is below3^m. The lower determinant has3-adic valuation
b<m. An integer eigenvalue3^m would divide that integer determinant, which is
impossible. Cayley–Hamilton then rules out a nonzero geometric observed
sequence of that ratio; seven matching consecutive macro-time values already
suffice for the contradiction.

The pure B case is also valid. Nilpotence modulo3 survives every positive
power, so all nonleading characteristic coefficients of N_B^m vanish modulo3.
At argument3^m its determinant term has valuationm while all other terms have
larger valuation. Therefore the lower block has no such eigenvalue. The upper
count becomes exactly geometric after one step, leaving a nonzero lower
geometric sequence if the claimed full ready count were geometric. Positivity
is not a heuristic: E(2)>=H(T−r)>=1 by using singleton{r} as the first bin and
a Hamiltonian path of the remaining tournament as the second bin. Thus the
cancellation coefficient cannot vanish. This is a simpler positivity witness
than requiring an entire Hamiltonian path ending at r.

The empty-tail legality guard is necessary and retained. From a nontrivial
strong parent, the new root reaches its unique port, that port reaches every
old vertex through the parent, and an old nonport reaches the new root. Thus
an empty-tail step remains strong. Nonempty legal tails end in1 and retain the
inherited strongness guarantee. The unique-start condition is inherited from
the audited funnel grammar. Nothing in the new algebra asserts that arbitrary
six-tuples represent tournaments.

Independent computation used two methods different from the original receipt:

* Fraction-based Gaussian elimination gave determinants552960 and811008.
* Direct permutations with weakly increasing cut positions enumerated the
  six cover coordinates for every one of the seven literal control parents.
  A separately constructed ready child was counted by direct permutations.
  All six coordinates and all seven full ready counts matched. No compiler,
  subset-DP counter, or convolution helper was imported for these checks.

Reproduction:

    python -B -O 05-knowledge/results/cycle_opening_six_return_20261112/dynamics/audit_compiler.py

The receipt is `../dynamics/compiler_independent_receipts.json` and reports
PASS. The original checker was also rerun and reports PASS for510 periodic
programs and seven literal controls. These are finite checks; the proof's
quantifier over all E/B periods comes from the norm/divisibility argument.

The corrected sibling reference is `../dynamics/PROOF.md`. Its order-six
carrier is a different assertion from the rank-six information bound here.
Neither result closes the repeatable ready-source construction or Chan–Pak.
