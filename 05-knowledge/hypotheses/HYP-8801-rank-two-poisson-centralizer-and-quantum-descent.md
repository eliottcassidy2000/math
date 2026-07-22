---
id: HYP-8801
title: "Complete and quantize the rank-two Poisson centralizer without doubling"
status: >
  RESOLVED at the Poisson gate by THM-2044; quantum gate moved to HYP-8802.
  The required second polynomial canonical pair exists after changing the
  preliminary Hamiltonian mate by an explicit polynomial connection. Direct
  Weyl-symmetric quantization has nonzero ordering anomalies.
source: codex-2026-07-21-DC2-JC2
related:
  - THM-2042
  - THM-1345
  - THM-2044
  - THM-2045
  - HYP-8802
  - THM-1300
---

# HYP-8801 -- the two gates below the rank-two Poisson witness

## Resolution update

THM-2044 reconstructs and proves the complete Poisson endomorphism. The
elementary invariants in THM-2042 were the wrong centralizer coordinates for
the final connection; importing the sheared THM-1300 coordinates and solving
one cubic coefficient identity produces the required `(T,S,D)`.

The second gate remains open and is now HYP-8802. The computation there proves
that naive Weyl ordering fails but repairs the entire `R`-column exactly. The
coupled `T`-column and flatness relation remain.

The owner-supplied abstract announced a polynomial completion `(R,T,D,S)` of
the THM-2042 pair satisfying two canonical bracket relations and all cross
brackets zero. THM-2044 reconstructs explicit formulas independently, audits
all six brackets exactly, and transports the exactly-three fibre of THM-1300.

The constructive route from the original scaffold had two separate gates.

1. **Poisson centralizer completion (closed by THM-2044).** Replace the elementary joint invariants
   `I=q^3L^2` and `J=(3xq-1)^3L`, whose bracket has the residual density in
   THM-2042(7), by polynomial `T,S` with `{S,T}=1`.
2. **Quantum descent (open as HYP-8802).** Choose orderings and lower-order corrections
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
