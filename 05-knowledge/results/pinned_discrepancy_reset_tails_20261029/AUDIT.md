# Scope and independent-path verification

2026-10-29. Same-agent proof review; independent mathematical/computational
paths where specified. No independent-agent audit is claimed.

## Proof checks

- The carrier requires every nonidentity multiplier to have h-1 invertible,
  and odd subgroup order. Composite Fourier modes are retained. The column
  identity `C^t C=Nr I` and row identity `sum C_(k,i)^2=N` have different
  roles and are not interchanged.
- The pin LP has a box constraint, and its dual has an l1 constraint.
  The dual objective subtracts the free-column support function. No pin
  can move during rounding. The fractional criterion is an iff only for
  families with uniformly bounded r and with an O(sqrt N) allowance.
- Lovett–Meka Theorem4 was checked at the primary source. Its arbitrary
  starting-point lemma is essential; an incidence-set theorem alone would
  not justify real Fourier entries. Compactness gives the exact existence
  statement; the packet does not claim an implemented edge-walk algorithm.
- The aligned hostile subtracts **two** free-column costs: removal from
  the aligned total and potential opposition. The correct expression is
  N/(2r)-4rm. An arbitrary pin assignment need not admit small norm.
- First reset means all k valuation-one steps followed by the first
  valuation a>=2. Both the actual odd-step clock k+1 and total division
  clock k+a are retained. ROOT j=1 is killed; it is not an active loop.
- The inverse-fibre calculation uses the original w at each original
  preimage. The elementary dominated-convergence bounds are supplied
  separately for a and k. Endpoint residue masses are not replaced by
  equal thirds. The two forbidden/unequal phase controls are exact.
- The word theorem at depth b uses only Fatou's lower-bound direction.
  It does not infer an upper bound by interchanging an infinite word sum
  and a limit. The original cylinder density and endpoint slope both
  contribute to the coefficient 3^L/4^A.
- The independent finite-depth power envelope degrades with depth. Uniform
  tightness would already imply boundedness of every source, since every
  source has positive prior mass; neither boundedness nor universal ROOT
  entry has been proved.
- The coherent ramified product has a regular kernel; both invariant
  subspaces are necessary for its norm formula. This obstruction does not
  apply to arbitrary noncoherent lifts, which retain additional variables.

## Computational controls and boundaries

1. `check_pins.py`: numerical LP/search only proposes witnesses. Rational
   primal and dual vectors and integer interval Fourier arithmetic certify
   the reported bounds. All legal completions are checked for the N=43
   cases. Hostiles include almost fully pinned aligned banks with free
   signs still available. Trigonometric error is accumulated over **all**
   columns in both primal and dual evaluations.
2. `check_resets.py`: literal odd-map trajectories verify 20,000 first
   resets, 1,002 inverse fibres and the terminal-family iff. One million
   original weights have explicit rounding and omitted-tail errors. The
   infinite terminal-family mass has its own geometric remainder proof.
3. `check_words.py`: 819 bounded-bank words have independently tested
   affine source cylinders, native endpoints and one-bit hostile changes.
   Four-block prefixes on 10,000 actual sources preserve the original
   label. The actual ring-difference tournament at625 vertices is compared
   with the tensor formula; both integer invariant-space identities pass.
4. Every check uses explicit exceptions, not removable `assert` statements.
   Normal and optimized executions must give identical receipts.

The first100-bit interval calculation rejected its own excessive Horner
width; increasing internal precision to192 bits repaired the calculation
without weakening the requested error bound. This was a detected numerical
implementation issue, not evidence for any mathematical claim.

## Resulting live targets

- Bounded-r logarithm-free existence: **PROVED from CITED input**, with
  the explicit fractional pin criterion. Better constants, a practical
  deterministic constructor and growing-r families remain **OPEN**.
- Independent first-reset residual bound and exact terminal family:
  **PROVED**. Every finite reset depth also has a proved weaker envelope
  and a critical lower-tail obstruction. Uniform-depth positive or signed
  residual control remains **OPEN**.
- Fully coherent ramified tower quasirandomness: **REFUTED by exact quotient
  amplification**, a scoped deduction from the incoming attachment theorem.
  Noncoherent low-norm towers are a separate **OPEN** construction target.
