# Mid-session integration: common response quotient, different exact carriers

While the local proofs and checks were running, the session fetched
`fe848abe68` and `602286dc6f`. Their mathematical content was inspected
before rebasing, not inferred from commit titles.

## What overlaps, and what extends it

The incoming [live-divisor transfer](../defect_boundary_20261020/boundary.md)
uses \(M_d=\#\{\text{old cycles divisible by }d\}\), discarding a coordinate
once no smaller cycle can inspect it. It also explicitly states the more
general future-response criterion in terms of
\(J_f=\sum\gcd(f,\ell)\). That is the principle independently developed
in this packet's Proposition 1. It must not be advertised here as a new
principle beyond the incoming work.

The implementations differ:

- Incoming: descending divisor bank through order100, then unit-tail
  reconstruction; 3,663 peak stored states and 66,404 transitions in
  the nontrivial-cycle transfer.
- Here: direct future-response bank, additionally deleting coordinates
  whose lengths exceed the remaining *state-specific degree budget*;
  2,953 peak states and 90,219 transitions through100, including the
  fixed-point stage. The transition totals count different stages,
  so they are not a direct runtime comparison.
- Here: signed substitution excess is carried separately from actual
  cycles, computing ordinary pair-free core counts and fixed-point-free
  core mass through100. The incoming note already proves the underlying
  signed ordered-list transform, including its phase factors.
- Here: A000568 is extended through200. The output through100 agrees
  with the independent recorded hash cited by the incoming packet.

The difference is thus a refinement and an independent check, not two
unrelated solutions of an allegedly unsolved numerical recurrence.

## Structural results imported after reading the proofs

The incoming [sparse theory](../defect_boundary_20261020/sparse_theory.md)
proves, at feedback budget \(k\), moved support at most \(3k\), pair-free
core order at most \(4k+1\), exact weighted-block straightening, and
weight capping at \(k+1\). It compiles finitely many capped weighted
templates into an eventual polynomial count. Its top degree is \(2k\)
and leading coefficient \(1/(2^k k!)\).

The [next-coefficient proof](../defect_boundary_20261020/next_coefficient.md)
classifies the next layer by a shortened gap, a missing gap, or a
shared defect endpoint, obtaining
\[
a_{n,\le k}=
\frac{n^{2k}-3k n^{2k-1}}{2^k k!}+O_k(n^{2k-2}).
\]
The exact budget-two polynomial is already supplied there. Re-proving
only a fixed-budget finiteness statement here would be weaker.

Our run result preserves the *chosen path*, rather than only the
minimum-feedback class. Proposition 6, obtained after reading the
incoming viewpoint, is the new connecting calculation:
\[
f(L_P(T^*_n))=r-1\qquad\text{when all }a_i\ge2.
\]
The proof has both an explicit optimal order and edge-disjoint
directed-triangle witnesses. The checker includes an independent
subset-DP optimum through order10. This illustrates why raw defect
count and minimum feedback must be distinguished.

The incoming [one-chord grammar](../defect_boundary_20261020/grammar.md)
is a complementary extension of the original source family. It changes
one distance-two chord, retains two tail lengths, proves a rank-three
bilinear transfer, and identifies the actual partner
\(TT_m\to C_3\to TT_r\). Our construction instead varies all marked
Hamiltonian paths of the unchanged tribonacci source. These parameters
and coverage claims must remain distinct.

## Small graphs and weighted boundary responses

The incoming [signed-interaction web](../defect_boundary_20261020/web.md)
adds information beyond our signless-incidence discussion: signed
adjacent Hessian entries reconstruct tournament edges up to converse;
at four vertices their signs distinguish a vertex star from a face.
It carefully separates graph complementation, tetrahedral opposition,
and characteristic-two cut/cycle intersection.

Our four-vertex kernel calculation is complementary: all three graph
types have signless rank three, but the free coordinate is spread over
all vertices on a tree and is localized at the isolate on a triangle
plus isolate. Neither calculation is replaced by just noting that
both line graphs of a claw and a triangle are \(K_3\).

The incoming [boundary-support asymptotics](../boundary_support_20261020_asymptotics.md)
separates a node count, a presentation-fiber size, and the actual
matching of line endpoints. Its strip example has equal support at
two powers but different nonproportional weights. This is the same
information discipline needed by the signed recurrence: preserving
possible continuations is weaker than preserving their coefficients.
No Hamiltonian-path extremal estimate from that note is needed here.

The incoming [105 dictionary](../defect_boundary_20261020/incoming.md)
also links budget-four leading matching templates with the prior
four-chart Hamming atlas while retaining endpoint order and binary
addresses. That is a genuine marked pairing map. It strengthens the
earlier number-105 archive synthesis; no bare equality of counts is
used as a substitute for that map.

## A second incoming wave during integration

The final fetch brought `b57c87b03d`, `bccc56795e`, and `fcf72e76dd`.
Their proofs were read before resolving the shared navigation and correction
ledger conflicts; both independent correction entries were retained.

The new [uniform prefix obstruction](../pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md)
addresses the earlier affine-only palindrome gap in a different
way. A rerouting supported strictly in one half of a path has a reflected
copy with disjoint support. If both path reversals were anti-automorphisms,
their product would have the rerouting's order. An even-order prefix
therefore supplies an obstruction. The proof also gives explicit violated
arcs, so it tests actual anti-automorphisms rather than constant partner sums.
Three forbidden sign patterns followed by unit reindexing by `2` or `-2`
exclude all normalized odd circulants of order at least seven.
We checked the displayed increment formulas and the small `m=4` exception.

The companion
[strong-core reduction](../merged_metagraph_purity_and_strong_core_20261020/PURE_BLUE_STRONG_CYCLE_REDUCTION.md)
uses Camion Hamiltonicity and the explicitly cited Lu crossing theorem to
reduce strong pure-blue atoms to this odd-circulant case. The new prefix
proof is a repair route for the classification, while the original
THM-4730 affine-only implication remains false. No part of our recurrence
or run theorem depends on either proof.

There is a provenance boundary on further corollaries. The incoming
all-order node-count/branch note invokes companions named
`FINITE_ATOMS_AND_RIGID_CLOSURE.md` and
`endpoint-transport-and-complement-leakage.md`. Neither filename was present
in the fetched `origin/main` tree at `fcf72e76dd`. We therefore do not
promote those dependent all-order stack/branch corollaries into this
packet. Their missing dependencies are a retrieval limitation, not a
counterexample to their statements. The explicit prefix obstruction is
independent of them.

The newly available
[Petersen support reduction](../boundary_support_20261020_synthesis.md)
is directly useful for our continuation distinction. Four claw modules
can be replaced by two while retaining all six ordered ports: the active
three-state transfer blocks change from `I+5J` to `I+J`. Their nonzero
supports agree, so every exterior has the same extension-existence
predicate, but the multiplicities differ. This reduction covers the
specified claw-ring family, not all snarks.

That gives a precise hierarchy of our interfaces:
retaining support suffices for a Boolean existence test; the Burnside
recurrence must retain numerical continuation weights; marked line
transport must additionally retain the actual path and ports.
The same quotient need not serve all three purposes.
