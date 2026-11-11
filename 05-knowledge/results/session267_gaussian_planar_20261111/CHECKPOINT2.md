# Checkpoint 2: from a scalar bridge to all-row and recursive controls

The first checkpoint established the exact Student/Gaussian rotation
identity. Three stronger results now make the research program concrete.

1. **A complete rational-carrier test.** For reduced P(t,x)/Q(t,x) with
   Q(0,x)=1, every coefficient has weight 2j−m≤B if and only if
   w(Q)≤0 and w(P)≤B. At B=5, this is exactly all-row rigid-debt
   vanishing in the inherited chart. It requires the full actual carrier,
   not a fitted finite prefix. [Proof](planar/EXTENSION.md).
2. **A precise substitution law.** For full imprimitive wreath actions,
   N_fix(H wr K)=H^m semidirect N_fix(K). The fixed-sheet gate can hide
   arbitrary inner groups while seeing the outer quotient. A nontrivial
   nilpotent outer block quotient is excluded, but cyclic inner blocks
   need a separate affine-source recognition theorem.
   [Proof and exact controls](topology/boundary_recursion_bridge.md).
3. **Two kinds of complexity.** Arbitrarily high polynomial degree can
   occur in one unchanged covering-topology orbit. Meanwhile an exact
   universal algebraic row test requires at least floor(m/2)−2 scalar
   equations even after paying the source scalar. A growing marked
   Gaussian polynomial stores those debts without losing them.
   [Degree](degree/PROOF.md) · [State lower bound](degree/ALGEBRAIC-STATE.md).

This does not close JC(2). It replaces vague high-degree/topological
concerns with explicit obstructions: a forbidden rational pole, a lost
source-depth coefficient, an unrecognized affine intermediate chart, or
a missing descent/termination rule. The Gaussian null-packet lane and
final synthesis are being integrated next.
