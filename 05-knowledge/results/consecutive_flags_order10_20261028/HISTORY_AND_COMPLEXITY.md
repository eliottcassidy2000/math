# Historical recovery, external sources, and ETH scope

**CITED** primary literature / **PROVED** elementary output bounds /
**CONDITIONAL** complexity transfers / **OPEN** proposed applications.

## Months-old enumeration work and what changed

1. March 7 and March 16 prototypes (`04-computation/a000568_speedup.py`,
   `a000568_v3_s90dm.py`) evaluate the odd-cycle Burnside sum. Recent
   [Session87 repairs](../kernel125_closed_carriers_20261028/HISTORY.md) correct
   an even recursive restart and an omitted odd-cycle restriction. Those scalar
   computations do not materialize tournament representatives or attachment data.
2. April approximation proposals confuse small relative error with exact
   integer rounding. The [later audit](../strong_attachment_tail_compiler_20261028/HISTORY.md)
   preserves the valid relative asymptotics and rejects the rounding claim.
   Today the exact scalar channel already runs through 300; extending this
   session's *bank* to ten is a different kind of progress.
3. The May 5 reflection
   [Two Models, One Staircase](../../../07-reflections/two-models-staircase-recursion.md)
   supplies the useful fixed-Hamilton-path strip representation. Its dimensional
   comparison is valid, but its description of the spaces as asymptotically
   identical needs qualification: their dimension ratio tends to one, while the
   fixed-spine space is only a 2^(-(n-1)) fraction of all labelled tournaments.
   At order ten that fraction is 1/512. A warning and mistake-lineage entry now
   repair the inference without deleting the historical discussion.
4. The October 24
   [parent-nine packet](../session51_structural_errors_relational_packets_20261024/parent-nine/README.md)
   already used parent-only canonical growth and a safe singleton shortcut. It
   predicted the order-ten scalar count using a separately justified correction
   to a refinement trace. It did not supply this complete order-ten bank.
5. The October 27 [nauty packet](../order10_prime_types_n11_ledger_spectral_blackcore_20261027/README.md)
   did enumerate all order-ten classes, including the prime census. Our claim
   is a different refined constructive implementation, compressed consecutive
   dictionary, and full alternate audit, not historical first discovery.
6. Current [root-port work](../root_ports_closed_carriers_20261028/CARRIERS.md)
   supplied the stable-cell refinement integrated here. Concurrent
   [Session89](../session89_consecutive_tournament_pairs_20261028/generation/ALIGNED_GENERATION.md)
   gives exact independent-pair and overlapping-stream gluing formulas. Its
   n=5 pair bank has 279,968 classes and only five concrete residual mask-action
   recipes. Reusing those recipes and reusing deletion parents are compatible
   optimizations, but their timings cannot be multiplied into a claimed speedup.

## Reproducible curiosity pulls

Using `Random(20261028).sample(sorted(candidates),3)`, where candidates are
`07-reflections/*.md` with matrix, tournament, graceful, square, staircase or
a000568 in the filename, the bounded header pulls were:

- [June dominance/condensation](../../../07-reflections/dominance-dodge-is-tournament-condensation-and-it-is-orthogonal-to-L-kps.md).
  Its peeling analogy suggests separating a cheap transitive shell from a
  difficult strong core. **Stopping boundary:** no explicit equivalence between
  LRC runner dominance and tournament composition was imported from the prose.
- [July prime-grid/sector carries](../../../07-reflections/lrc14-prime-grid-cover-sector-carry-and-threshold-tournament-codex-20260714.md).
  It asks for scale, carry and microphase data beyond an orientation. This is
  an orthogonal example of why a representation's query language matters.
  **Stopping boundary:** no new LRC bound follows from a tournament bank.
- [August D=0 square face](../../../07-reflections/d-zero-square-face-critical-value-collapse-and-source-row-invoice-codex-20260831.md).
  It distinguishes a critical-value coordinate from actual source-row existence.
  **Stopping boundary:** the coincidence of the word square supplies no map to
  square-sum labelling. A transfer would need an explicit operator and preserved
  feasibility predicate, unlike the signed-incidence map proved in THEORY P5.

These are provenance pulls, not truth-source promotions. They motivate the
concept board but are not dependencies of any new proof or numerical receipt.

## A concrete shared linear operator

The strongest additional connection is the signed root-incidence operator in
THEORY P5. Its forward map sums competitive e_i-e_j and collaborative e_i+e_j
columns into tournament-pair scores. Its transpose evaluates differences and
sums of vertex labels, respectively the fixed graceful and square-sum edge
equations. This gives a precise connection with Coxeter tournaments and the
existing LCA affine-transport solver.

The decisive hostile is the kernel: scores do not recover orientations, and
mod-two columns do not recover integer signs or magnitude. Future operations
testing those coordinates require a sidecar. Feasibility of the fixed linear
system remains distinct from choosing the edge lengths, signs, squares and
permutation bounds. This is the point where the general graceful and square-sum
existence problems remain open.

The primary literature checked for this session includes Eplett (1979),
[structural representation](https://doi.org/10.4153/CMB-1979-004-6), and
Kolesnik--Mitchell--Przybylowski (2025),
[Coxeter interchange graphs](https://doi.org/10.1007/s00026-025-00768-9).
Classical representation and root-system mechanisms are credited, not claimed
as new. The exact local contributions are the executed finite audits, the
rooted reset/deletion carrier experiment, and the normalization-specific
signed-score/deletion formulas and boundaries.

## ETH: distinguish the reported result from its conditional consequences

As checked on 2026-10-28, OpenAI's primary
[formalization scope document](https://github.com/openai/math/blob/main/lean/docs/378.md)
reports deterministic exact fixed-width #k-SAT in subexponential variable time
and polynomial space. Its stated Lean bound is

    (L+1)^O_k(1) * 2^O_k(n (log log n)^2 / sqrt(log n)),

where L is canonical input length. It separately reports a sharper October 25
manuscript bound with one power of log log n, and a sparse-clause improvement.
This session did **not** build or audit that Lean development or independently
verify its manuscript. The repo's
[trusted-advances intake](../../reference/TRUSTED-ADVANCES-2026-10-27.md)
classifies it as a preprint claim; this packet retains that status.

The following deductions assume a concrete uniform 2^o(n)-time SAT algorithm
with a polynomial input-length factor. The reported counting theorem is stronger
than a decision-only negation of ETH; it should not be silently inferred from
decision alone. Likewise one must specify uniform algorithms when discussing
arbitrarily small exponential rates.

- **Parameterized reductions:** an O(n)-variable, polynomial-size fixed-width
  SAT encoding transfers the assumed subexponential bound to n. An n^2-variable
  encoding gives 2^o(n^2), not automatically 2^o(n). Tseitin variables, clause
  width and arithmetic bit length must all be included. The foundational source
  is Impagliazzo--Paturi--Zane,
  [Which Problems Have Strongly Exponential Complexity?](https://cseweb.ucsd.edu/~paturi/myPapers/pubs/ImpagliazzoPaturiZane_2001_jcss.pdf)
  (2001), which develops reductions preserving subexponential complexity.
- **Search:** decision self-reduction can recover one satisfying witness using
  polynomially many calls with fixed assignments. Counting requires a counting
  theorem; listing every witness requires time for every emitted witness.
- **P versus NP:** subexponential time need not be polynomial time. No such
  collapse is deduced here. An ETH-based lower bound loses its premise; that
  alone does not construct a practical solver for every problem to which it was
  applied.
- **Research opportunities:** compile the remaining sign/length choices of
  the existing affine tree solver into a counted constraint system, or compile
  port compatibility and symmetry-breaking into SAT. The missing work is an
  explicit encoding with controlled variables and a useful runtime, not a
  claim that a general existence conjecture has been solved.

## Unconditional barriers survive

Let a_n be the number of ordinary tournament classes. Every S_n orbit has size
at most n!, so

    a_n >= 2^(n(n-1)/2)/n! = 2^(n^2/2 - O(n log n)).

Any algorithm explicitly listing one record per class requires at least a_n
output steps, independently of ETH. It therefore cannot run in 2^o(n^2) total
time under that output contract. The aligned square-bit pair census similarly
has at least 2^(n^2)/n! outputs. Output compression, as in our parent/address
dictionary, is a different contract; expansion still incurs the output cost.

The scalar a_n is much easier to specify. The classical odd-partition Burnside
sum has exp(O(sqrt(n))) terms and polynomial bit work per term using integral
factorial-scaled accumulation; its output has O(n^2) bits. This already
separates exact scalar counting from full representative enumeration without
any SAT hypothesis. A000568 recreation should therefore keep separate goals:
scalar terms, representative banks, full automorphism actions, and metagraph
transition outputs.

Finally, no algorithmic breakthrough reconstructs data that a representation
has erased when two inputs give the same retained value and different requested
answers. The 81/27 group difference, the interior-deletion forest and the
anti-involution centralizer gap are such unconditional witnesses. A smaller
carrier is justified by an exact quotient for the stated operations, not by
assuming that future computation will recover an omitted bit.
