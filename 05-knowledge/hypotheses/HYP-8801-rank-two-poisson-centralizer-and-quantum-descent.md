---
id: HYP-8801
title: "Complete and quantize the rank-two Poisson centralizer without doubling"
status: >
  OPEN / active. THM-2042 constructs the canonical pair (R,D) and an explicit
  joint-invariant scaffold. The target is a polynomial canonical pair (T,S) in
  the joint centralizer, followed by an ordering/correction scheme whose Weyl
  commutators are exactly canonical in A_2. A Poisson completion alone does not
  imply DC(2), and its base shadow does not imply planar JC.
source: codex-2026-07-21-DC2-JC2
related:
  - THM-2042
  - THM-1345
  - THM-1300
---

# HYP-8801 -- the two gates below the rank-two Poisson witness

The owner-supplied abstract announces a polynomial completion `(R,T,D,S)` of
the THM-2042 pair satisfying two canonical bracket relations and all cross
brackets zero. Once the missing formulas are available, the first task is an
independent exact audit of all six brackets and of the claimed three-point
fibre.

The constructive route from the proved scaffold has two separate gates.

1. **Poisson centralizer completion.** Replace the elementary joint invariants
   `I=q^3L^2` and `J=(3xq-1)^3L`, whose bracket has the residual density in
   THM-2042(7), by polynomial `T,S` with `{S,T}=1`.
2. **Quantum descent.** Choose orderings and lower-order corrections
   `Rhat,That,Dhat,Shat in A_2` so their six commutators are exactly canonical.
   Principal symbols only prove the top-order Poisson identities; Moyal/PBW
   correction terms can survive. The cotangent/Hamiltonian-dual construction
   avoids this issue by doubling to four Weyl pairs, but that proves a statement
   about `A_4`, not `A_2`.

For the planar Jacobian conjecture the useful positive target is different:
classify how a smooth non-coordinate polynomial such as `R` can acquire a
Hamiltonian mate after symplectic stabilization, and prove that no such mate
exists without the extra pair. This converts the proposed counterexample into
a candidate obstruction invariant for de-stabilization.
