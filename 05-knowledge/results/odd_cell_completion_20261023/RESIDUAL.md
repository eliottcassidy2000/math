# Reuse residual potentials only along authenticated chart symmetries

23 October 2026. **PROVED elementary scoped transport rule.** This extends
the incoming [branch-conditioned test](../odd_cell_gluing_20261023/README.md),
which requires both residual-cycle constancy within each compatible chart
branch and agreement of the branches' projected baselines.

## Equivariant certificate reuse

Let S be the retained compatible chart branches. For each s, the transport
model has actual row and column sets, a multiset of labelled channels,
integer bounds, prescribed margins, and the SAME target observation P.
Do not replace these models by the union of their channels.

Suppose a finite group K acts on S. For each action s->g s retain actual
bijections of rows, columns and channels carrying the first model to the
second, preserving incidence, bounds and margins. Require, in addition,

    p_(g e) = p_e

for every channel's vector contribution to P. These conditions are checked
on channel labels, not inferred from equal anonymous signatures. Selected
source owners, path order, converse channels and folded black-port weights
must be among the retained conditions whenever the consumer uses them.

**Theorem.** It suffices to run the inherited residual-cycle potential test
on one branch per K-orbit. If x is its feasible table and h its potential on
each residual strongly connected component, then

    (g x)_(g e)=x_e,   (g h)(g v)=h(v)

give a feasible table and a valid potential certificate on every transported
branch. Their projected baseline is the same. Between distinct K-orbits,
baseline comparison is still required.

Proof. The retained bijections give a bijection of the entire feasible
integer-table sets. They carry positive-capacity residual arcs with their
signs, capacities and P-costs to the corresponding arcs. Thus SCCs, cycle
costs and the potential equations transport exactly. Also
P(g x)=sum_e p_(g e)x_e=P(x). Consequently constancy and its value are the
same throughout an orbit. Different orbits have no specified transport and
still require the incoming cross-baseline check. QED.

If the observation coordinates are permuted rather than fixed, the conclusion
is P(g x)=sigma_g P(x); the baseline must be proved fixed by these coordinate
permutations before calling it invariant. Actual-owner projections generally
need this stronger check. A merged class/colour projection can be invariant
under genuine tournament automorphisms, but only when its bucket and channel
construction is also equivariant.

## How this reduces the odd-cell problem

The comparison isomorphisms between two fixed attachment systems, when
nonempty, form a torsor under the source's cell/mark-preserving automorphism
group. The new small-cell kernel rule bounds root charts and computes these
actual maps. If the port model is constructed covariantly from these maps,
they supply the action above, and all chart branches lie in one orbit. Then
one potential certificate and its explicit transported copies replace a
separate residual solve for every root chart.

This condition is substantive. If a branch chooses a new untransported port,
uses a different clock, discards a source mark, or changes a colour weight,
the claimed action may fail. The theorem is a verification rule for the
supplied model, not a proof that every candidate wiring model has this
covariance. Nor does PASS prove that every relaxed table is native: retain
one authenticated actual table to connect the invariant projection to the
real tournament wiring.

## Odd-phase hostile to premature union

Take two three-by-three blocks, each with row and column margins one.
Branch s in F3 permits only the matching j=i+s in BOTH blocks. There are
three legitimate branches. Let P be the number of off-diagonal matches in
block one minus that number in block two. Each branch has P=0 and a unique
table. The same cyclic phase is essential.

Unioning the channels first permits two independent permutation matrices:
36 tables, instead of the three legitimate ones. Some have P=3 or -3.
They are valid tables in the enlarged relaxation but not compatible shared
chart choices. This is the odd-phase version of the incoming two-bit hostile,
not a newly discovered failure mechanism.

The checker also includes a different failure: two individually unique
branches can have different projected baselines. Their residual-cycle tests
both pass, but no P-preserving branch transport exists. Thus neither
per-branch constancy nor an unverified symmetry claim replaces the baseline
condition. Converse and modulo-eight guards remain conjuncts; reducing P
modulo eight only certifies its residue, not its integer value.

## Late incoming extension: shared linear rows require coupled moves

Commit af1ba9bc0e supplies the
[bounded coupled-move theorem and finite-domain compiler](../session41_owned_attachments_fibonacci_20261023/residual-cycles/README.md).
When the feasible set has extra shared-owner equations, testing just the
ordinary residual cycles is no longer an iff criterion. Several switches
can be individually forbidden and jointly feasible. The incoming theorem
instead tests every feasible conformally primitive integer kernel move
within the declared finite widths. Its completeness obligation is explicit.

Thus the potential test above remains exact for its stated network model,
and sufficient for a more constrained subset, but it must not be used as
a necessary test after extra linear rows have been imposed. The same
equivariant reuse principle applies to an authenticated permutation of a
complete coupled-move bank: a coordinate permutation preserving the bounded
integer model preserves its kernel, conformal order and width bounds, hence
carries its primitive moves and their observed costs. Completeness still
has to come from the incoming theorem or a verified generator; a few
sampled moves cannot supply it.

The incoming [owner-indexed enumerator](../session41_owned_attachments_fibonacci_20261023/residual-cycles/OWNER_INDEXED_ENUMERATORS.md)
also gives exact contraction on supplied forests and tree decompositions.
Our alignment quotient explains an additional input requirement: a tree
of local BLOCK ORBITS is not already a tree of whole compatible chart
states. Two rigid blocks have a relative phase at their shared middle cell.
Supply those states before applying the forest contraction. A forest has
no cycle holonomy, but it can still lose middle-owner data under an earlier
quotient.

## Relation to the Fibonacci and seven-term side lane

The [overlap-base theorem](fibonacci/REPORT.md) supplies the preceding,
nonabelian consistency test. An actual loop must be checked in the chart
group before its compatible branches enter this additive transport model.
The commutator in the seven-cell group gives a concrete reason that scalar
or abelian potentials cannot substitute for that step. Fibonacci interval
products likewise remain conditioned on the SAME shared endpoints; an
interior port introduced later requires splitting the compressed interval.

The reusable procedure is therefore: authenticate joint charts; retain or
prove their action on the actual port model; test one residual representative
per authenticated orbit; transport the certificate; compare baselines across
the remaining orbits. None of these steps assumes universal metagraph wiring.
