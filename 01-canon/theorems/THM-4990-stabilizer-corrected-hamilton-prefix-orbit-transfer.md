---
id: THM-4990
title: "Integral path-orbit transfer requires pointwise stabilizers and becomes uncorrected at the last two unused vertices"
status: >
  PROVED, independently audited in Session232. FINITE-EXACT controls have
  their stated universes; inherited all-order theorem certificates remain
  explicit dependencies, not newly rerun claims. Chan-Pak completeness OPEN.
session: Session232-automorphisms-20261109
source: 05-knowledge/results/session232_automorphisms_20261109/orbit_transfer/PROOF.md
audit: 05-knowledge/results/session232_automorphisms_20261109/subgroup_transfer_audit/AUDIT.md
---

# THM-4990

[Full proof](../../05-knowledge/results/session232_automorphisms_20261109/orbit_transfer/PROOF.md) · [Independent audit](../../05-knowledge/results/session232_automorphisms_20261109/subgroup_transfer_audit/AUDIT.md)

Let G<=Aut(T), where T has n>=1 vertices. A directed prefix is an ordered path
on distinct vertices. Its state is (S,u), the used subset and final vertex.
Choose a representative for each G-orbit b of states, and define

    k_b=|{g in G: g fixes S pointwise}|,
    c_b=number of G-orbits of prefixes of state type b,
    N(b,a)=number of legal next vertices taking the representative into type a.

Then the exact integral transition is

    W(b,a)=N(b,a)*k_a/k_b,
    c_a=sum_b c_b*W(b,a).

Initialize c_b=1 on singleton vertex orbits. At the full layer the sum of c_a
is H(T)/|G|. This statement also holds for arbitrary nonempty finite digraphs.

## Proof and the two distinct measures

An ordered prefix's stabilizer is exactly the pointwise fixer of its used set.
That fixer acts on the allowed next vertices of target type a; every such
vertex stabilizer has size k_a. Thus all its orbits have size k_b/k_a and
their number W is integral. It does not depend on earlier ordering beyond
the state. Deleting the last vertex gives a unique prefix orbit, proving the
recurrence. Full prefixes have trivial stabilizer.

For groupoid masses m_b=c_b/k_b the recurrence instead uses N directly:
m_a=sum_b N(b,a)m_b, with singleton initial value1/k_b. N is correct for
mass, W for coarse orbits. Mixing initial coarse counts with mass transitions
is incorrect. A source followed by C3 has N=3,k_b=3,k_a=1 at the decisive
step, so W=1; the naive mixed recurrence incorrectly gives3 full orbits.
Replacing the pointwise fixer by the setwise endpoint-state stabilizer also
fails, as the separate C3=>K1 audit control demonstrates.

## Sharp terminal freeness

For tournaments k_b=1 whenever |S|>=n-2: a nonidentity automorphism supported
on at most two vertices could only transpose the endpoints of an arc, which
is impossible. Thus every transition with such a source has W=N. The bound
is sharp for n>=4: a transitive chain of n-3 vertices followed by C3 has a
prefix with stabilizer C3 at length n-3. More generally minimum group motion
mu gives freeness for |S|>n-mu.

## Connection with existing moving-core states

For a marked tournament (C,M), a cut-preserving subgroup acts freely on ordered
spanning one- and two-path covers. It consequently divides each of the seven
Session231 coordinates(h,q,A,G,T,S,E). Their bilinear ordinal composition is
compatible with division by the explicit product action. If two star-closed singleton-padded constructions use the same marked core
(C,M), every marked-core automorphism extends by the identity on the added
vertices. If their Hamilton counts differ by delta, then 2*|Aut(C,M)| divides
delta. A unit difference±2 requires a rigid marked core.

This is weaker than the inherited gcd(h,A,G)=1 content filter, not a new
unit-gain construction. C3 with singleton cut is rigid as a marked object but
has gcd(h,A,G)=3, so the converse fails. The twelve-bit Session231 core is
independently verified rigid; this observation does not solve baseline placement.

Finite controls cover all1099 labelled tournaments through order5 and all74521
endpoint-state orbits;328 transfer edges require correction. Separate audit
code verifies the recurrence and hostile stabilizer choices. The source paper
comparison is methodological only. The recurrence computes the exact orbit
statistic but neither proves its inherited7/21 exclusion by itself nor proves
unrestricted Chan–Pak completeness.

Inherited routes: LEM-003-aut-acts-freely-on-ham-paths.md;
THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md;
THM-4960-twelve-bit-moving-core-consecutive-odd-intervals.md.
