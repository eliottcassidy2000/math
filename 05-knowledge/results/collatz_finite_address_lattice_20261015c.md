# Finite witness lattices, diagonal address folds, and the cost of changing bases

2026-10-15. **PROVED:** finite attainment for fixed-source signed graph
quotients; the address maps and basis-cost statements below.
**FINITE-EXACT:** 750 toy graph/order checks and the stated algebra controls.
**CITED:** four primary-source results with their scope retained.
**OPEN:** any all-scale Collatz payment or ROOT theorem.

This is the structural companion to
[annular weighted controllers](collatz_annular_overlap_20261015c.md).
[Control program](../../04-computation/experiments/collatz_address_controls_20261015c.py),
[output](collatz_address_controls_20261015c.json).

## 1. Recover the existing logic work before using the papers again

The four requested papers were already studied in
[the October12 logic-lenses packet](logic_lenses_jc2_collatz_hn_20261012.md),
with its corrected canon THM-4620 through THM-4626. In particular,
[THM-4623, functional-graph-order-topology-endomorphisms-and-koopman-torsion](../../01-canon/theorems/THM-4623-collatz-functional-graph-order-topology-endomorphisms-and-koopman-torsion.md)
already gives an orbit-order reformulation of divergence and a Koopman
kernel dictionary. Re-proving those weaker reformulations would not advance
the current target. Here the object is instead the **ordered family of
finite weighted certificates for one fixed finite source set**.

The relevant primary sections were checked again:

| source | checked result | use here / limit |
|---|---|---|
| Pinsky, [2305.11382v1, §2, Lemma2](https://arxiv.org/html/2305.11382v1#S2) | rigid connected pointed graphs encode colours and directions while preserving an automorphism group | retain typed source and terminal ports; automorphism-group preservation does not itself preserve native step counts or Gram weights |
| Xu–Ji, [2609.19612v1, Theorem1.1 and §3](https://arxiv.org/html/2609.19612v1) | countable meet-continuous **complete lattices** are Scott sober; directed tests detect Scott openness | ask whether a limiting certificate is finitely detectable; neither an orbit poset nor the countable poset of finite edge sets automatically meets those hypotheses |
| Fernandes–Mezabarba–Rodrigues, [2609.20140v2, Theorem1.2 and §2](https://arxiv.org/html/2609.20140v2) | existence of bases for all characteristic-zero vector spaces implies choice even without Foundation | use the explicitly indexed finite source coordinates; an abstract basis neither supplies native receipts nor preserves the terminal cube |
| Chini, [2609.20184v1, §2.2, Definitions2.3–2.4](https://arxiv.org/html/2609.20184v1#S2.SS2) | polynomial-kernel constraints can be organized by kernel configurations | track the finite relation subspace of source labels; the o-minimal/model-companion conclusions are not asserted for Collatz |

These are lenses, not imported hypotheses establishing our inequalities.
The proofs below are finite linear algebra and literal graph reasoning.
The Pinsky HTML conversion ends during §2; the claim used here is the fully
displayed Lemma2, not an unchecked later construction. No novelty is claimed
for finite-dimensional stabilization or the powerset Scott topology.

## 2. A fixed finite source has a finitely attained certificate optimum

Fix finitely many original source ports, including ROOT, and finitely many
real weighted rows on them. Each observed edge is an authenticated native
step, carrying the minus sign, with observation stopped at ROOT. In the
union of observed edges allow only these deductions:

1. a finite path to ROOT fixes the label with its actual parity;
2. two finite paths to a common vertex identify their source labels up to
   their clock difference;
3. a finite path to a witnessed non-ROOT cycle fixes that source label to
   zero, because the label represents ROOT support.

Missing edges and an unwitnessed infinite tail impose no zero label.
Every deduction keeps both integer clocks; a residue or common symbolic
tail is not a receipt. This is the same rule set as
[G1](collatz_gram_slack_20261015b.md), extended here from finite graphs to
an arbitrary union of observations.

The deductions form a signed partition of the finitely many source ports,
with a distinguished ROOT block and a zero block. Within any other block,
choose its least indexed source as representative; this is only a sign
gauge, not a declaration that the source is terminal. Push the original
rows to these blocks and form G. Its entrywise norm does not depend on
that sign gauge or on renaming the blocks.

**Theorem F1 (PROVED).** For any observed-edge union E, finite or infinite,
there is a finite subset E0 of E which induces the same signed partition
on the fixed source ports. Consequently

    ||G(E)||_1 = min_(finite F subset E) ||G(F)||_1.            (F1)

The minimum is attained even at equality with a proposed budget; a strict
margin is unnecessary in this fixed-source representation.

**Proof.** There are finitely many source ports. For every nontrivial block,
choose one finite common-endpoint receipt linking each source to its chosen
representative. For the ROOT block choose finite ROOT paths; for each zero
block choose paths to its witnessed non-ROOT cycle. These receipts are finite
by the definition of the deductions. The union of their edge sets is finite
and induces every relation on the source ports. It cannot induce any new
relation not valid in E. Thus the signed partitions agree. Any finite subset
of E induces a coarser partition, and signed identification/zero deletion
cannot increase the Gram norm. This proves (F1). ∎

Equivalently, in a finite-dimensional space with basis e_s, the relation
subspace is generated by e_s-epsilon e_t and e_s for witnessed zero labels.
A finite basis of that subspace has finitely many native witnesses. This
version exposes the finite-dimensional stabilization mechanism, but does
not license arbitrary replacement of the remaining terminal coordinates.

For an increasing observation sequence starting with r free terminal
classes besides ROOT, there can be at most r strict increases of the source
relation subspace, and hence at most r strict decreases of the Gram norm.
Moving a lone terminal along its path only changes a common sign gauge and
pays nothing. An arbitrarily long plateau may precede the next decrease.
The proof gives **no computable stopping test** announcing that the final
partition has arrived, no bound on the lengths of its receipts, and no
uniform bound as the source window grows. These are essential limitations.

This strengthens finite-graph closure only in its limiting scope. The
[October7 terminal-basis packet](collatz_terminal_basis_20261007.md) already
proved finite closure and a minimum for its frozen ROOT-word grammar;
neither statement is replaced by F1.

### Scott-open acceptance, on the correctly typed space

Let the native edge universe be countable, and order its **full powerset**
by inclusion. This is an algebraic complete lattice; its finite subsets
are compact. It is generally uncountable. The poset of finite subsets alone
is countable but is not directed-complete.

For a fixed source Gram and budget B, the predicate ||G(E)||_1<=B is an
upward set. If a directed union satisfies it, F1 supplies a finite witness
E0, and directedness places all of E0 in one member of the family. That
member already satisfies the predicate. Thus acceptance is Scott open,
including its non-strict numerical boundary.

This is a direct proof for our certificate lattice, **not an application**
of the countable-lattice theorem to the Collatz orbit space. Sobriety here
does not force any source to ROOT. Our new point is the finite witness for
the weighted acceptance predicate, not a second reformulation of convergence.

## 3. Put 10, 14, and 16 in one explicit address diagram

The recovered divisor work is
[the October1 triplet analysis, §4.3](procgen_tcpc_20261001_tournament_clock_prime_collatz.md)
and [the divisor-balance lattice reading](collatz_mod6_20260917_divisor_balance_family.md).
The owner's original prime/squarefree-composite/square-containing partition
of p^2qr is 3–4–3. The proper-divisor rank partition is also 3–4–3 but
exchanges p^2 and pqr; only the latter respects complementation as rank
reversal. This distinction is inherited, not newly discovered here.

Here is a precise common construction for the three requested counts.
It is an offered interpretation of 4–6–4 and 4–8–4; the targeted repo search
recovered no earlier maintained definition of those two phrases.

Start with four binary addresses (a,b,c,d) in {0,1}^4, and fold

    pi(a,b,c,d)=(a+b,c,d).

The target is {0,1,2} x {0,1}^2, the exponent box of p^2qr.

| observable / operation | exact count |
|---|---|
| keep all four address bits | 16 |
| group the 16 addresses by a+b=0,1,2 | 4–8–4 |
| remove empty/full addresses; group by total bit rank 1,2,3 | 4–6–4, total14 |
| fold a,b to their sum | 12 target addresses |
| remove the target bottom/top; group by total exponent rank | 3–4–3, total10 |

Every target with first coordinate1 has **two** preimages; the others have
one. Retaining this multiplicity recovers 14 from the ten proper targets.
The fold preserves total rank and intertwines complement:

    pi(1-a,1-b,1-c,1-d)=(2-a-b,1-c,1-d).

The generating polynomials show the same operation:

    (1+t)^4  ->  (1+t+t^2)(1+t)^2.

The missing coefficient in the first factor is exactly the forgotten
distinction between (1,0) and (0,1). One may also group the proper target
addresses by exponent of p; that gives another 3–4–3 partition, again
different from the original prime/composite partition.

**Connection to payment.** In one two-point fibre, row weights (+1,-1)
collapse to zero, while (+1,+1) collapse to two. Their original rank-one
Gram norms are both4; their new norms are0 and4. The quotient cardinality,
fibre size and complement diagram are identical. The sign/weight sidecar
decides the entire payment. This is exactly the elementary fibre mechanism
in (A6). No map from these divisor addresses to all Collatz trajectories
has been produced; that further connection remains an analogy.

A second, separate algebraic reading is useful for joint covariance:
for two two-dimensional spaces U,V,

    Sym^2(U direct_sum V) = Sym^2 U + (U tensor V) + Sym^2 V

has dimensions3–4–3, while End(U direct_sum V) has block dimensions4–8–4.
The proper exterior algebra of a four-dimensional space has dimensions
4–6–4. These identities do not identify the three carriers. Their concrete
lesson for the H/Delta Gram is to retain its mixed block, not to infer a
bound from the matching dimension10.

## 4. Changing generators: a useful Pell control, with its cost visible

Reuse, rather than re-prove as new, the matrices from
[the quadratic-cycle companion](quadratic_cycle_phase_20261015b.md):

    B=[[0,-1],[1,1]], U_u=[[1,u],[0,1]], C_u=U_u^(-1)B^(-1)U_u,
    B^3=C_u^3=-I,  B C_u=[[1,u],[u,u^2+1]].                  (F2)

B preserves K=[[2,1],[1,2]] and C_u preserves U_u^T K U_u.
For u!=0 the product has determinant1 and trace u^2+2>2. At u=2 it is
conjugate by [[1,-1],[0,1]] to [[3,4],[2,3]], the exact Q–Pell matrix
recovered alongside [THM-4630](../../01-canon/theorems/THM-4630-hofstadter-q-jitter-is-the-silver-ratio-arma-fixed-point-rho1-equals-1-minus-sqrt2.md).
Its Q white-sampling assumptions remain hypotheses.

Individual invariant metrics therefore do not provide a common nonexpanding
metric for arbitrary alternating words: the hyperbolic product would
contradict one. Transporting a metric along the whole word can make a formal
isometry, but its comparison with the original metric carries the growing
Pell factor. That factor cannot be erased from the target invoice.

The relevant positive control for terminal controllers is different. For
any real matrix A,

    ||A^T G A||_entrywise1 <= ||A||_infinity^2 ||G||_entrywise1, (F3)

where ||A||_infinity is the largest absolute row sum. This follows by
expanding entries and summing; the constant is sharp already for a rank-one
diagonal G supported at a row attaining that row sum. Native terminal
identifications have one signed entry per row (or zero), so this factor is
at most1, even with changing path lengths and changing terminal names.
Their extra improvement is precisely the weighted overlap in (A6).

An arbitrary basis change does not have that certificate semantics.
For example r=(1,2) and A=[[1,-2],[0,1]] give rA=(1,0). The original cube
maximum of (rz)^2 is9; the apparent fresh-cube maximum after changing basis
is1. But the correct new coordinate is A^(-1)z=(z1+2z2,z2), whose first
coordinate can be3. The energy is still9. The lost datum is the transformed
terminal cube (and possibly the marked ROOT coordinate), not a minor
normalization. This is why choosing an abstract basis cannot pay the gap.

## 5. What transfers, and the stopping boundary

| source -> target | map / preserved predicate | lost information and required sidecar | decisive test |
|---|---|---|---|
| coloured directed receipts -> graph presentation | encode typed ports; recover each original native arrow | auxiliary encoding paths are not native steps; keep clocks and weights | replay original arrows after decoding |
| directed observation union -> finite acceptance | F1 selects finite source relations | no uniform witness-depth bound | arbitrarily delayed merge followed by a single norm drop |
| four address bits -> divisor box | pi preserves rank and complement | first-pair order and fibre multiplicity | same fibre, Gram norms0 versus4 |
| polynomial/kernel viewpoint -> source relations | record witnessed e_s-epsilon e_t | ordinary linear span does not preserve bounded-label geometry | the unimodular9-to1 false-cube example |
| generator-local metrics -> common controller | retain the full product and coordinate conversion | switching toll | hyperbolic B C_2 versus individual finite orders |

The new research move is to combine **finite detectability** with a
**quantitative source-specific invoice**. F1 closes the issue of whether
an ideal signed quotient needs an actually infinite certificate. It does
not establish that its best Gram norm is below the required budget.
The annular experiments begin to pay that second obligation using the
actual weights; (A7) proves that the missing all-scale ingredient must go
beyond the halving submodel. This is a usable division of the remaining
work, with both the successful mechanism and its failure boundary retained.

No new META-PATTERNS card is promoted: the cross-domain examples support
existing cards, but the new uniform arithmetic mechanism is still missing.
