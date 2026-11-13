# Which questions now deserve the work?

**PROVED scoped / CITED inheritance / FINITE-EXACT / OPEN, separately below.**
The useful outside connection is a transfer of a proof obligation, not
an assertion that two counting problems are equivalent.

## 1. Euclidean descent and electrical networks: preserve less

**Standard comparison, with elementary proof.** A connected two-terminal
multigraph has the pair (t,f): spanning trees, and spanning two-forests
separating its terminals. Adding a parallel terminal edge gives
(t+f,f): trees using that edge correspond to the old separating forests.
Appending a series edge and moving the terminal gives (t,t+f): a
separating forest either isolates the new terminal or includes the new
edge and separates the old terminals. Start with one edge, pair (1,1).

The two operations are the unipotent matrices

    L=[[1,1],[0,1]], R=[[1,0],[1,1]].

Every positive coprime pair is reached: subtract the smaller coordinate
from the larger until (1,1); positivity and gcd are preserved and the
coordinate sum strictly decreases. Reversing the subtractions gives a
network. This is the same elementary descent underlying the
[Calkin--Wilf rational tree](https://cecas.clemson.edu/~calkin/Papers/recount.pdf).
The counting realizations are different. No novelty is claimed for this
standard comparison. `transfer_probe.py` checks all 255 words through
depth7 directly by spanning-tree and separating-forest enumeration.

Source: networks and coprime pairs. Target: a tournament selector with a
free auxiliary boundary coordinate. Map: compare their integer response
matrices and reverse rules. Preserved predicate: exact response under
the specified composition on each side. Lost information: simplicity,
complete orientation, and tournament-realizable boundary states. Needed
sidecar: an actual tournament constructor for each proposed rule.

The decisive test is now [the factor-five theorem](GROWTH.md): no
unipotent return can inhabit our existing half-repair context class,
even after visiting other endpoint types. A count-preserving translation
of all networks is impossible anyway: the simple cycle C7 has seven
spanning trees, while no tournament has seven Hamiltonian paths.

**Changed question (OPEN):** can an existential construction of the
prime-lane ternary children discard the paired repair and rebuild a
different interface? We should not demand a slow universal return of d.
This is a weaker requirement than the grammar we have been trying to
close. The existing spanning-tree/block dictionary is already
[THM-4986, spanning-tree twin](../../../01-canon/theorems/THM-4986-the-7-21-omission-proof-has-an-exact-spanning-tree-twin-block-floor-crevasse-lemma-and-the-first-creator-contrast.md);
this session adds the explicit obstruction to the proposed shear transfer.

## 2. Zaremba: the right quantifier is existential in the second coordinate

**CITED comparison, already proposed in the repo.** In the continued-fraction
problem the denominator is prescribed and a suitable numerator may be
chosen. Bourgain--Kontorovich's 2014 theorem obtains density-one coverage
for a bounded alphabet; it does not prove every denominator is represented.
Their treatment separates local congruence conditions, orbit growth, and
representation multiplicity. See the [primary paper](https://annals.math.princeton.edu/wp-content/uploads/annals-v180-n1-p03-p.pdf),
Theorem1.2 and sections1.2--1.3. The repo already makes this comparison in
[Session25, B6](../kps_s25_euler_field_bridge_law_20261111/README.md).

Source: scalar projections of continued-fraction matrix orbits. Target:
the projection (H,E)->H of a tournament grammar. Map: a construction word
to its matrix product, with the literal tournament retained as a witness.
Preserved predicate: existence of some auxiliary coordinate and word.
Lost information: endpoints and entry if one keeps only H. Needed
sidecar: actual entry proof and a bound on representation multiplicity,
not just a list of congruences. Our matrices also have nonunit determinants;
the cited analytic theorem cannot simply be applied to them.

**FINITE-EXACT hostile control.** The nine inherited matrices, starting
at (5,9), have every parity-allowed H residue for all moduli 2 through100
and additionally 128,210,256,315,420. These are full finite orbits, not
finite-depth samples. At modulus45 the orbit contains all 2,025 pairs;
at modulus420 it contains 44,100 pairs and all 210 odd H residues.
The residues of 7,21,53,084,231 pass every tested modulus.
Yet this same grammar has the proved O(X^(3/4)) upper bound in
[LANGUAGE](../endpoint_context_recursion_20261113/LANGUAGE.md).

**Changed question (OPEN):** for a proposed enlarged grammar, can we
prove enough distinct actual representations in every required target
fiber, after both congruence restrictions and growth have been accounted
for? Before a large residue search, test whether the language already
fails the weighted branching bound. More local coverage is not progress
against a proved global sparsity bound. Conversely, rapid individual
growth alone does not rule out a rich enough scalar projection.

## 3. Multiplicative atoms: distinguish a grammar obstacle from a proof obstacle

**Inherited PROVED reduction.** The desired odd-count monoid has atoms

    odd primes p!=7;
    7p for primes p>=5, p!=7;
    49,63,343.

The finite atoms already have witnesses. See
[THM-4094, two-prime-lane completeness](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md)
and [Session193, Theorems1--2](../session193_triangulated_proof_20261105/prime_induction/PROOF.md).
The latter proves a sufficient induction using 3h+2 and 3h+4 only when
the child lies in an infinite atom lane. It does not require a common
parent witness, retained old arcs, or a common endpoint preparation.

Source: factorization in the monoid of allowed odd integers. Target:
strong tournament witnesses for the two infinite lanes. Map: ordered
strong-component joins implement multiplication. Preserved predicate:
actual count. Lost information: marked endpoints; it also does not
produce a strong witness for a nontrivial product. Needed sidecar: a
strong realization of each atom. The cheapest test is factorization
before treating a hostile count as a new global obligation.

**New arithmetic triage, FINITE-EXACT:**

    53,084,231 = 43 * 1,234,517,

and both factors are prime. Thus the defect-one obstruction there is not
itself an atom-lane obstruction. The count is already globally attained
by the inherited direct construction; this factorization does not claim
new witnesses for its factors. Likewise the induced-extension failures
13->27,25->51,29->57,35->69 have composite children whose counts can be
assembled from attained smaller factors. The child63 is an exceptional
atom, but it is already seeded. In contrast, the inherited two-apex
41->127 obstruction has a prime child and rules out a genuine proposed
implementation of the ternary induction within that hereditary class.

**Changed question (OPEN):** for each unseeded atom N, can at least one
method realize N, assuming smaller allowed counts are attained? Different
methods may handle different N. Do not insist that one method be closed
under every target transition. This is a proved sufficient objective,
not a new reduction and not a proof of its missing existence premise.

## 4. Partition functions: arithmetic freedom is not geometric freedom

**CITED exact representation:** H(T)=I(Omega(T),2), with vertices of
Omega the odd directed cycles, adjacent when they intersect. This is
the Grinberg--Stanley odd-cycle formula, also Irving--Omar Corollary20;
see the [core-paper route](../../reference/CORE-PAPERS.md).
For an abstract clique K_m, I(K_m,2)=1+2m. Thus abstract conflict graphs
make every odd target trivial; tournament realizability carries the work.

Source: hard-core partition functions. Target: tournament odd-cycle
packings. Map: T->Omega(T). Preserved predicate: the count at fugacity2.
Lost information: orientation and the consistency of overlapping cycles.
Needed sidecar: a realizable cycle system, not only an independence
polynomial. The cheapest hostile test is the no-disjoint-cycles class.

That class is already completely analyzed in
[intercyclic rays](../feedback_word_hamilton_spectrum_20261020/intercyclic_rays.md):
its spectrum is the one-apex spectrum plus two explicit dyadic rays,
and it misses29 at every order. So a proposed universal rank-one
partition-function construction stops at29, not just at7 and21.

**Changed question (OPEN):** can a local operation control interactions
between two or more disjoint cycle families while keeping their actual
orientation data? The elementary equation for H alone is insufficient.
This lane should be tested at29 before any large target search.

## 5. Incoming result: exact multiplication can spend a finite interface

During this session, origin/main added
[Session302](../session302_ternary_ray_birth_20261113/README.md).
Its scoped result replaces a uniform TT2 module with C3 when the exterior
allows at most two runs through that module. Here F_2 counts ordered
covers by two nonempty paths. The response pair (H,F_2)
scales from (1,2) to (3,6), so the surrounding strong tournament count
multiplies by three. It constructs 63->189->567->1701.

This does not contradict the factor-five theorem: it uses a module
replacement, with no preserved nonzero half-repair pair. Its own bound
F_2(C)>=(|C|-1)H(C) shows why a fixed F_2/H interface cannot renew at
unbounded order. The positive move is real; the renewal resource is finite.
Session302 also applies the S-unit theorem to fixed exponential charts;
that scoped obstruction is inherited, not re-proved here.

The shared strategic lesson is precise: distinguish a count identity
from a replenishable interface. A promising next experiment must name
what resource it renews, and show the new block satisfies the next rule's
entry conditions. For the Chan--Pak atom induction, it may instead be
better to avoid demanding renewal at all and rebuild for each target.

## Decision after this session

1. Prioritize the inherited existential prime-lane induction, permitting
   core replacement and target-specific witnesses. Its scalar premise is
   weaker than our current paired grammar and is sufficient for Chan--Pak.
2. Retain the half-repair grammar as a controlled construction laboratory,
   with the all-order factor-five cost visible. Do not expect bigger
   contexts or slope detours alone to provide a slow return.
3. Test new families against growth, local admissibility, and actual
   entry separately. Use29,127 and53,084,231 for different questions;
   their obstructions are not interchangeable.

No every-atom construction, positive-density theorem for a new family,
or Chan--Pak proof has been obtained in this session.
