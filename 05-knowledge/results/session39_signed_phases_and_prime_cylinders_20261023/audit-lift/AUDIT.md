# All-orders lift audit

Status: PROVED scoped transfer and lower bound; FINITE-EXACT independent
rooted and exact-incidence controls. Earlier frozen audit packets unchanged.

## Rooted profile transfer

The final all-order lift proof passes. For L(T)=C3[T,{a},{b}], old-root
subsets split according to whether a and b are included; each resulting
rooted type is a deterministic extension of an old rooted type. For a
new root, the analogous extension uses old unrooted counts. These are
recoverable by summing rooted counts over actual roots and dividing by
the old motif size. The empty old subset is a separate fixed object.
There is no assumption of vertex transitivity or of a common isomorphism
for all subsets.

A strong old tournament of order N has degrees1..N-2. After lifting its
old degrees are2..N-1, while a,b have degrees1,N. The unique extrema
intrinsically recover the old tournament, so nonisomorphism is preserved.
Strongness is immediate from the cyclic block arrangement. An old
anti-isomorphism extended by swapping a,b gives an anti-isomorphism of
the lift. The old phase-blind equitable partition gains two singletons
and remains equitable with phase-independent quotient, so ordinary
unindividualized refinement cannot restore the lost phase.

For n>=11, take r=0,2,1 at residues0,1,2 modulo3. Then g=(n-2r)/3>=3.
The transferred pair requires cutoff k(n)>=g+1, implying k(n)>=(n-1)/3
and liminf k(n)/n>=1/3. This is confined to the stated static singly-rooted
histogram observer plus ordinary refinement. It is not a reconstruction
counterexample, indexed hypomorphy, a complexity bound, a sufficient cutoff,
or a claim about the narrower prime-cyclic subclass.

The final indexed-data guardrail is correct: on the fixed actual vertices
(0,0),(0,1),(g-1,0), phase+ gives a transitive triple, while phase- gives
(0,0)->(0,1)->(g-1,0)->(0,0). Lifting preserves it. Equality of rooted
histograms therefore must not be upgraded to equality on indexed subsets.

## Exact-E transfer checked independently

For an old vertex v, cyclic triangle counts change as c'_v=c_v+1.
For old u,v, e'_uv=e_uv. Also e'_va=e'_vb=1, and
c'_a=c'_b=e'_ab=N. Hence the outgoing row of an old vertex is obtained
by adding1 to all head-c entries and adding the pair(N,1). The row at a
is {(N,N)}. The row at b is {(c_v+1,1):v old}.

The old global c_v multiset is itself determined by E, since for each
root u, c_u=sum_(u->v)e_uv: each cyclic triangle through u has exactly
one outgoing edge at u. Thus equal exact E gives equal lifted exact E,
independently confirming the inherited lift theorem cited by the producer.
Self-converse preservation gives the corresponding converse equality.

For every n>=17, the same residue representation has g>=5. Apply at most
two lifts to the three-phase exact-E collision family to obtain three
intrinsically distinct strong self-converse classes at order n. This is
not a certificate of alternative endpoint-capacitated actual port wiring.

## Checks and scope

Standalone check_lift.py imports no producer code. It reconstructs both
phase tournaments and enumerates full rooted canonical profiles pointwise
for g3,4 and zero/one/two lifts, verifying equality, anti-isomorphisms and
unique new degree extrema. Orders checked:9,11,13,12,14,16. Separately it
reconstructs all three g5 phases and checks exact E directly in both
orientations after zero/one/two lifts, at orders15,17,19. All pass.

The final producer replay also passes its exact frozen receipt, including
the indexed triple guardrail. Its recorded canonicalizer dependency is
the separately frozen rooted checker. Producer pins and SHA256SUMS freeze
the approved exact bytes. None of these finite controls replaces the
all-order transfer proof.
