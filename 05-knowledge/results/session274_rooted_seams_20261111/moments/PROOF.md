# Seam-sector spectroscopy: forests hide the charge, loops reveal it

**Status: PROVED IN NOTE** for the identities and all-order moment expansion;
**FINITE-EXACT** for the small-order distinctness controls. Chan–Pak remains OPEN.

## 1. Carrier and provenance

Fix a Hamiltonian path P=(1,...,n). The free tile graph is G=K_n minus the n−1
path edges. A bit t_ij=1 means j→i, for i<j−1; path arcs are always forward.
The geometry lane proves that rooted cut-seams generate exactly t+Cut(G),
and that cycles of G classify these orbits. We call these *seam sectors*.
This is not the quotient by tournament isomorphism, nor the metagraph whose
lines complement every tile. All averages below are uniform over distinct
framed tilings in one seam sector.

The first-moment forest argument was independently found by the geometry
lane and this lane during the session. The whole-profile version is proved
there. Here we extract an explicit recurrence and show exactly where a
higher-order observable first sees the hidden information.

Historical inspiration: THM-2810, `factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go`,
uses moments to recover information that a single projection loses. This is
an analogy of observable hierarchies, not a transfer of its factorial theorem.
THM-4878, `graded-spine-grammars-cut-state-automaton-hankel-and-seam-ranks-cutwidth-censuses`,
already supplies finite-state path counting for bounded backward cutwidth.
Our sector average is different: backward cutwidth is unrestricted, but
averaging kills forest dependence exactly.

## 2. Explicit universal sector mean

Let mu_n be the mean of H in any seam sector. Set mu_1=1. Then

  mu_n = 2^(1−n) sum_S 2^r(S) (n−|S|)!,                    (1)

where S ranges over subsets of the n−1 path edges whose maximal consecutive
runs all have even length, and r(S) is the number of nonempty runs.

**Proof.** The geometry theorem permits taking independent fair tile bits.
For any vertex permutation, its chance to be Hamiltonian is the product of
edge weights: 1 on forward path arcs, 0 on backward path arcs, and 1/2 on
chords. Factor out 2^(1−n). The remaining product gives factor 2 to each
forward neighboring pair, factor 0 to each backward neighboring pair, and
factor 1 otherwise. Expand these as 1+1 and 1−1 respectively.
A selected set of required neighboring pairs can occur in a permutation
only if each of its maximal consecutive runs is consistently directed:
a mixed direction at an internal vertex would require two immediate
predecessors or two immediate successors. Contracting the runs leaves
(n−|S|)! permutations. The two possible directions of a length-l run
contribute 1+(−1)^l. Thus odd runs cancel and each even run contributes2.
The empty subset contributes n!. This proves (1). ∎

This is a weighted run-enumeration algorithm, not enumeration of tournament
classes. One may compute the bivariate polynomial

  F_m(z)=sum_S 2^r(S) z^|S|,   m=n−1,

by a three-state automaton (outside a run, odd run, positive even run):
start (A,O,E)=(1,0,0), and append one bit by

  (A,O,E) <- (A+E, 2z A+z E, z O).

The accepting polynomial is A+E. Then apply z^k ↦ (n−k)! and divide by2^(n−1).
This is an exact finite transfer description for the universal first moment,
not an instantaneous A000568 formula or a full Hamiltonian spectrum theorem.

Initial means for n=1,...,6 are 1,1,2,4,79/8,29. Nonintegrality causes no
problem: sector sums, rather than their means, are integers.

## 3. All moments are loop-consistency sums

Fix a representative tiling t and integer q≥1. Range over ordered q-tuples
of vertex permutations. Discard any tuple that uses a backward path arc.
Every chord used by a tuple prescribes the orientation needed for that
permutation. Under a cut switch x∈F2^V, it imposes

  x_i+x_j = t_ij + epsilon_ij,

where epsilon is the required tile bit. If repeated chord uses demand
opposite directions, discard the tuple. Otherwise let F be its chord graph.
The equations are consistent iff on EVERY cycle C of F,

  sum_C t_ij = sum_C epsilon_ij       in F2.              (2)

For a consistent tuple the probability is 2^(-rank(F)), where
rank(F)=|V(F)|−components(F), counting only vertices incident to F.
Consequently

  E_sector[H^q] = sum_(consistent tuples) 2^(-rank(F)).  (3)

**Proof.** Expand H^q as the sum of products of path indicators. Vertex
switch labels may be sampled uniformly even when G is disconnected: each
cut vector has the same number of preimages. The incidence matrix of F
has rank |V(F)|−components(F), and its image is exactly the edge labels
with zero cycle sums. Therefore each compatible system has the stated
probability. Sum indicators. ∎

For q=1, F is a subgraph of one Hamiltonian path, hence a forest: condition
(2) is vacuous and all sectors have the same first moment. For q≥2, unions
can contain loops and retain precisely the original sector's cycle charges.
Thus the hierarchy has a concrete meaning: *a single path cannot probe
holonomy; correlated paths can*. No continuum gauge-theory theorem is being
invoked.

## 4. Small decisive controls, not range evidence

`python3 moments/check.py` from the packet directory enumerates every framed
tiling at orders2..6, partitions it by actual cut-space translations, and
computes H by subset dynamic programming. This is a structural control,
not an attempt to enlarge Chan–Pak numerical coverage.

At order5 there are four sectors, each of size16. Every mean is79/8, but
the second moments are217/2,109,221/2,227/2, all different. Thus variance
alone distinguishes every sector in this initial charged case. At order6
there are32 sectors but only19 distinct second moments. This is a precise
hostile example to generalizing the order5 decoder. No claim that all
moments recover every labelled sector is made: framing symmetries can force
identical distributions.

## 5. Constructive reading

The geometry lane's rooted birth reaches every charge extension of an old
sector, yet mean H does not say which target a selected extension realizes.
Equation(3) explains what was missing from a scalar account: correlations
between path choices retain the new cycle coordinates. A next compiler can
retain a few physically selected mixed-path observables rather than only H
or its average. Whether such a bounded list suffices for exact odd-target
selection is OPEN; the order6 collision warns against replacing the full
interface by variance alone.
