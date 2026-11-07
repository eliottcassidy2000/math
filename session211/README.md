# Session 211: local path data and global tournament structure

**Status: four independently audited packets of elementary structural results and exact finite controls. The Chan–Pak Hamiltonian-path spectrum conjecture remains OPEN.** No literature-priority claim or extension of the previous finite coverage range is made.

The common thread is precise: local path-cover and comparison data control what survives at larger scales. Matching small numbers, matching vertex counts, and even matching an entire cover polynomial do not identify the underlying class of objects. The packets recover earlier permutohedron and recursion work, the September square-sum constructions, the October arithmetic reconciliation, and Session210's recognizable gate. Their source manifests retain the individual provenance and corrections.

## 1. Square-sum onsets and the 3/10 obstruction

[Square-sum port bridge](square-sum-port-bridge/PROOF.md).

Keep the axes separate: `7,21` are forbidden tournament **path counts**; `3,10` are their normalized indices `(H-1)/2`; `15,25` are square-sum **vertex-order milestones**, with a third onset at 23.

The inherited square-sum counts, independently replayed, are `p(Q23)=3` and `p(Q25)=10`, where reversal identifies undirected chains. Therefore no universal assignment to tournaments can preserve `(H-1)/2=p(Q)`: it would realize the forbidden counts. This rules out that exact scalar transport, not every useful encoding.

The positive correspondence retains named ports. At the three onsets 15,23,25, deleting the new degree-two vertex produces respectively 1,3,10 spanning two-path forests, with designated ports `{1,10}`, `{2,13}`, `{11,24}` in different components. Nine additions share old count zero and new degree two but have different outputs, so those two scalars cannot determine the next count. The complete local census runs through order 25; separate finite controls replay 22 inherited Pythagorean constructions through order 3600.

## 2. What `n+2` and `2n` actually do

[Tournament recursion bridge](recursion-bridge/PROOF.md).

Source/sink padding `Sigma(T)` changes order by two and preserves `H`. Reversing just its apex arc gives the strong gate `Gamma(T)`. For a nonempty core, every old Hamilton word survives and exactly `P_T(2)` new words appear, where `P_T(2)` counts ordered two-path covers, allowing empty paths:

    H(Gamma(T)) = H(T) + P_T(2),
    HC(Gamma(T)) = H(T).

In the coherent permutohedron this is exactly the creation of `P_T(2)` sinks by reversing one parallel edge class. The gate's path-count minimum at core order n is `1+2^n`; among nontransitive cores it is `3+3*2^(n-1)`, with the proved equality cases. The gated path values at most 31 are exactly `{3,5,9,15,17,27}`, a statement about this restricted class.

Skew-Sylvester doubling has a different geometry. The collection of reflected doubled words spans an affine diagonal permutohedron section, with feedback energy `binom(n,2)+2E`. The section's edges lift through square faces; it is generally not itself a face. Endpoint distributions count two exact one-switch sectors, while paths changing sheets more often remain a separate contribution. Equal final vertex counts do not make padding, gating and doubling commute, as sink and gate degrees prove.

## 3. A common algebra does not repair missing comparisons

[Arithmetic and cover bridge](arithmetic-cover-bridge/PROOF.md).

The ordered-cover polynomial and gate operator extend to every loopless digraph:

    P_Gamma(D)(q) = q[P_D(q+1)+qP_D(q)-P_D(q-1)].

For nonempty cores the gate still transports directed paths to directed cycles. But the number of missing unordered comparisons and the number carrying both directions are unchanged. No gate iterate turns a nontournament into a tournament.

The literal directed arithmetic graph `x->x+2`, `x->2x` on `1,...,N` has no Hamiltonian path and exactly four ordered two-path covers for every `N>=3`. Its gate is therefore strong with four Hamiltonian paths and no Hamiltonian cycle. Symmetrizing the same arithmetic graph gives positive even path counts; completing it to the numeric transitive tournament gives one. These operations preserve different information.

The distinction cannot be recovered from aggregate covers alone: a minimal order-three digraph with one missing and one reciprocal pair has the same polynomial `P(q)=q^3` as a transitive tournament. All nested gate cover and path-count observations agree forever, although only one tower consists of tournaments.

## 4. Scores preserve the macro decomposition, but not the gate

[Score permutohedron and gate boundary](score-permutohedron/PROOF.md).

Tight Landau inequalities recover the ordered strong-component partition and its product face from the labeled score vector. Adding the gate moves the entire padded core score face into the larger score polytope's relative interior, with exact unnormalized facet slack one at the completed gate. This score geometry is distinct from the full coherent permutohedron orientation that remembers every tournament arc.

A score-preserving triangle reversal can erase the gate and changes the cycle count from `H(T)` to `H(T-v)`. A four-cycle reversal gives a two-component cover formula with specified endpoints. The explicit equal-score pairs have cycle counts `11->7` and `37->21`; the resulting 7- and 21-cycle tournaments have no gate. These are **Hamiltonian-cycle counts**, never the forbidden Hamiltonian-path counts. The triangle-switch response remains a tournament path count and therefore cannot reach 7 or 21; the four-cycle response counts endpoint-constrained two-path covers and can. Thus the score projection loses the microscopic condition defining the faithful gated carrier, even at the exact same labeled score point.

## The combined micro/macro picture

The square-sum births and the padding-to-gate flip use the same organizing idea: two-path states invisible to a one-path count can assemble into new global paths. The square-sum map requires specific endpoint ports; the gate sums over all ordered two-bin covers. The four-cycle score switch reaches the same necessary refinement: specified two-path data can give cycle counts excluded from the one-path spectrum. This is a shared mechanism, not a count-preserving map between the square graphs and those tournament witnesses. None of these correspondences justifies replacing its full input by the old count.

The geometric connections are equally typed. Padding supplies a score face; the gate enters the interior; doubling supplies a diagonal section and square-face lifts. These are different operations on different faithful or lossy representations. Score data, aggregate cover data, and named-port data are complementary observations, not a claimed hierarchy in which one always determines the others.

The sharp boundary is **realizability**. Counts can be transported beyond tournaments, but the comparison alphabet and the recognizable gate must remain attached if tournament-specific theorems are to survive. The legitimate cross-scale arithmetic mechanism here is the inherited Pythagorean/reflection construction of square-sum chains, not a numerical identification of `{3,10}`, `{15,25}`, and `{7,21}`.

## What remains open

- Realizing every positive odd tournament path count other than 7 and 21 is still open. The restricted gate bounds and finite controls do not establish completeness or extend the inherited coverage range.
- A count-controlling growth construction must return to a valid marked tournament class. The two order recursions, the common cover polynomial, and the score projection do not by themselves provide it.
- The full contribution of multiply interleaved sheet paths is not determined by the two computed one-switch sectors.
- No independently classified number-theory family with precisely the missing values 3 and 10 was found, and no Collatz or Lonely Runner consequence is claimed. The eventual square-sum threshold is inherited as a citation, not reproved in this session.

## Verification

Each packet includes its verifier, exact outputs, source manifest, hashes and independent review:

- [Square-sum review](square-sum-port-bridge/audit/REVIEW.md): complete count census through 25, all fourteen onset paths, and the separately scoped macro controls.
- [Recursion review](recursion-bridge/audit/REVIEW.md): all 1,099 labeled tournament cores through order five, plus independent folded-square and sheet-sector checks.
- [Arithmetic-cover review](arithmetic-cover-bridge/audit/REVIEW.md): all 4,166 labeled loopless digraphs through order four and an independent smaller full-permutation implementation.
- [Score-gate review](score-permutohedron/audit/REVIEW.md): exact score geometry, triangle and four-cycle changes, and independent counts of the literal cycle witnesses.

Normal and optimized Python replays agree. The finite checks support the reviewed all-order proofs; they are not substitutes for them.
