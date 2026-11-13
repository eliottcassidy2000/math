# A fixed arithmetic orbit, actual uniform hardness, and physical descent

2026-11-13. **Chan–Pak and ordinary Collatz remain open.** There are three
different results: a logical classification of fixed Chan–Pak, an exact
arithmetic encoding, and a tournament-native inverse on an infinite
restricted grammar.

## Logic and arithmetic dynamics

- [Logic and bounds](LOGIC_AND_BOUNDS.md): the SCC product and Camion
  already bound the order needed for a target by the target itself.
  Busch sharpens this to an integer inequality 3^n<=h^3 after stripping
  singleton SCCs. Exact ordinary target membership is decidable.
  Chan–Pak is a universal Pi_1 statement; a false instance can be
  refuted by a bounded computation, so independence from PA would
  imply standard truth. Neither truth nor independence is established.
- [FRACTRAN orbit](FRACTRAN_ORBIT.md): compile the bounded counterexample
  tester into a prime-state generalized Collatz map. A designated orbit
  reaches 1 iff Chan–Pak is false, and otherwise diverges. A second
  **positive-slope** piecewise-affine map has global convergence to 1
  iff Chan–Pak is true. It uses bounded, restartable epochs and an
  odd reset lane, so arbitrary starting states are handled.
  [Epoch escape](EPOCH_ESCAPE.md) proves that every nonterminating
  orbit of this second map actually diverges, not cycles.
  Neither is the standard 3x+1 map.
- [Uniform hardness](UNIFORM_HARDNESS.md): varying the admission
  validator makes coverage genuinely undecidable even along one
  fixed Fermat target ray. The Pi_1-complete construction uses
  hereditary, ordinal-sum-closed, two-way padding-invariant
  tournament classes. One-way singleton padding can store time
  and raises the uniform problem to Pi_2-complete. A consistency
  search gives actual Gödel-II instances for the constructed
  classes, not an independence result for the ambient spectrum.
- [Guard closure](GUARD_CLOSURE_AND_HALTING.md): exact native
  applicability sets for FRACTRAN are a closure system under
  divisibility. Every finite closure system has a squarefree
  realization. Cyclic tournament priority has a sharp shielding
  criterion; it neither implies computational universality nor
  turns nontrivial fractional rules into autonomous residues.

## A native constructive step and its boundary

- [Paired-end descent](PAIRED_END_DESCENT.md): extend Session 322's
  finite factor bank to **all** strong one-apex dyadic reservoir
  blocks in ordinal combination. The target fixes its dyadic
  shell. A two-ended base-four carry filters the possible cut
  words; a divisor recursion then chooses actual components.
  Every successful deletion of a component or singleton more
  than halves the absolute Hamiltonian target, leaving the same
  physical cut. This is a native descent, but entering this grammar
  is not proved for arbitrary odd targets.
- [Scalar descent barrier](SCALAR_DESCENT_BARRIER.md) and
  [strong version](STRONG_SCALAR_BARRIER.md): a uniform homogeneous
  substitution that reads only H must be a monomial C h^k,
  even when the input is promised strong. A nontrivial marked
  SCC inevitably reads a higher cover coordinate. Any actual
  scalar analogue of 2h+1 or 3h±1 via these operations therefore
  needs more than the count as its state.

Replays check the literal bounded checker, compiler instructions,
an arbitrary-state toy epoch wrapper, all one-apex cuts through
chain length 14, the finite guard closure controls, and the Fermat
ray separation. The all-order proofs do not assume the conjectures
or infer infinite behavior from these finite checks.
