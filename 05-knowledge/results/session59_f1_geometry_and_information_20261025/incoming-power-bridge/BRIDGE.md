# Incoming conic powers meet cyclotomic dynamics and the graph degeneration

2026-10-25. **PROVED elementary bridge and scoped obstruction.** The good
conic, nested-set counting construction, canonical augmentation and graph
family are inherited. The two comparisons below are self-contained deductions,
not external novelty claims. No new theorem about integer-square existence,
blue paths, black-core connectivity, RH, or uniform tournament enumeration is
asserted. This is a focused integration pass, not recertification of the corpus.

## 0. Source pins and inherited correction boundaries

The incoming source is commit
`c3c7ae07e4c16eb0d1837254e7a2a777092d882b` of
`eliottcassidy2000/math`, specifically:

- `05-knowledge/results/ports_owners_f1_20261025/f1_bridge.md`, §§1–2;
- `05-knowledge/results/ports_owners_f1_20261025/anchor.md`, §§2–5;
- the corresponding `AUDIT.md`, proof reviews and corrections.

Session59 inputs are `cyclotomic-support/CYCLOTOMIC_SUPPORT_DYNAMICS.md`,
§§1–3, `graph-magic-degeneration/RESULT.md`, §§1–4, and
`toric-shared-roots/TORIFICATION_OBSTRUCTION.md`, §§1–3, under
`05-knowledge/results/session59_f1_geometry_and_information_20261025/`.
The exact source bytes inspected are recorded in `SOURCE_PINS.json`. The
cyclotomic and graph sources are also pinned to verified first checkpoint
`32693c14f4c76da929fc38b51a75c9acf9afb122`; the geometric addendum was a
local Session59 proof packet pending publication when inspected.

The immediate incoming successor
`e0e957e6749b3185cd4da3eff7631cb9c5fab342` was also inspected: its ONLY
change is the nine-line `ports_owners_f1_20261025/verify.out` recording the
complete normal/optimized replay as PASS. It changes no proof, hypothesis,
or correction used here. These commits are the incoming freeze boundary.

Incoming boundaries retained here:

1. The augmentation owner G_(a+1) is not the ordinary owner L_(a+1).
   Iterated rank transfer requires exact G-to-L recanonicalization. The regular
   pentagon is the literal failed-projection witness.
2. The parent decreases tournament order and feedback-vertex number by one,
   preserving maximum transitive size alpha. It does not automatically
   preserve strongness or a marked Hamiltonian path.
3. In the nested-set construction d=|A|−|B| must grow under powers. The old
   Boolean-family identification requires d=1 AND the stated power-of-two
   condition on |B|; d=1 alone is insufficient.
4. The good integral conic is UV=W², not the bad characteristic-two model
   x²+y²=z². Primitive-content normalization is an explicit operation.
5. The ordinary swap on X² fixes q points. The Euclid first-leg set has qd
   points, so for d>1 it cannot be identified with those blue diagonal paths.

## 1. Exact conic–cyclotomic bridge

Let C=Proj Z[U,V,W]/(UV−W²). The incoming Veronese isomorphism is

    v:P¹_Z → C,       [a:b] ↦ [a²:b²:ab].

For r≥1 its power map psi_r[U:V:W]=[U^r:V^r:W^r] satisfies

    psi_r ∘ v = v ∘ ([a:b] ↦ [a^r:b^r]).                 (1)

This is an equality over Z, including reduction in characteristic two.
On the dense torus of P¹, the coordinate a/b is simply sent to (a/b)^r.
Thus incoming conic powers and Session59's k=2 projective power action are
the same operation after the following explicitly stated realization.

### Proposition 1 (finite-field dynamics of the repaired conic)

Let Q be a prime power and choose a group isomorphism
C_(Q−1) ≅ F_Q^×. Adding an absorbing zero and taking projective classes gives

    P_(Q−1)(2) ≅ P¹(F_Q) ≅ C(F_Q).

Under this bijection the cyclotomic power F_r is psi_r. In particular

    #Fix(psi_r on C(F_Q)) = 2+gcd(Q−1,r−1).              (2)

For t≥0 and n≥1, the number of points with preperiod at most t and eventual
period dividing n is

    2+gcd(Q−1,r^t(r^n−1)).                              (3)

Exact tail-period counts follow by the same difference in t and Möbius
inversion in n as in Session59. These formulas include noninvertible r and
Q=2. The group identification is a choice; the power-equivariance itself
is independent of that choice.

**Proof.** The two singleton supports are fixed. Full-support projective
points have the unique form [u:1], u∈F_Q^×, and the map is u↦u^r.
A fixed point solves u^(r−1)=1, giving gcd(Q−1,r−1) solutions. The equation
psi_r^(t+n)(x)=psi_r^t(x) on the full support is
u^(r^t(r^n−1))=1, giving the gcd in (3). For a finite functional graph that
equality is equivalent to the stated tail-period condition. Equation (1)
transports the result to C. ∎

For example over F_7 the square map has three fixed conic points, two
additional periodic points of exact period two, one point of tail one and
period one, and two points of tail one and period two; the total is eight.
The two boundary points are included, not discarded as degenerate triples.

### Where the actual nested sets enter

Take B⊂A with m=|A|>n=|B|>0, initially coprime of opposite parity. The map

    (A,B) ↦ [m:n] ↦ [m²:n²:mn]

intertwines Cartesian powers (A^r,B^r) with the two power maps in (1).
Since m,n are coprime this gives an integral projective section, whose
reduction is defined at every prime; reduction commutes with powers. The
resulting finite-field points lie in the prime-field subpoints. Proposition1
applies to ALL finite-field points, not just this integer-reduction image.

The original orbit counts remain

    qd, (q²−d²)/2, (q²+d²)/2,  q=m+n, d=m−n,

with d_r=m^r−n^r. If noncoprime cardinalities are allowed, divide by
h=gcd(m,n) before forming an integral projective section; this commutes with
powers because gcd(m^r,n^r)=h^r. Retain h if the original cardinalities or
unnormalized orbit counts are required.

Reduction to a finite field forgets sizes, positivity, and most of d_r.
It does not identify a finite-field orbit with blue-path incidence or with
an additive integer-gap relation. In particular there is no conflict between
periodic finite-field dynamics and strictly growing positive d_r over Z.

## 2. Power maps survive the flat family, but their Frobenius property does not

Write the Session59 family as

    G_t = Spec Z[t][z,(1+tz)^−1],
    u ⊕_t v = u+v+tuv.

Its multiplicative coordinate after inverting t is x=1+tz. This is a
deformation of the dense torus, not a compactification of the full conic.

### Proposition 2 (unique power extension and a precise specialization obstruction)

For every r≥1 the map x↦x^r extends uniquely to a base-preserving
endomorphism of G_t, namely

    [r]_t(z) = ((1+tz)^r−1)/t
             = sum_(j=1)^r binom(r,j)t^(j−1)z^j.        (4)

The fraction is an integer polynomial before specialization. These maps
are group endomorphisms, obey [r]_t∘[s]_t=[rs]_t, and at t=0 become z↦rz.
Their tangent action is also multiplication by r. Applied coordinatewise,
they preserve the graph-potential and magic subfamilies and commute with
forest reconstruction, the shared-component obstruction map, and gluing.

For prime p the reduction of (4) modulo p is

    [p]_t(z) = t^(p−1) z^p.                             (5)

Consequently its additive fiber at t=0 is the ZERO morphism, not the
nonzero additive Frobenius z↦z^p. For a positive-rank graph carrier the
same distinction holds in forest coordinates. There is therefore NO
base-preserving integral extension of the generic torus p-power map to
this family whose additive fiber modulo p is the usual Frobenius.

**Proof.** Equation (4) is polynomial by the binomial theorem, and
1+t[r]_t(z)=(1+tz)^r is a unit. It therefore defines the indicated scheme
map. Inverting t identifies it with torus powering, proving composition
and the group law on a dense locus; these polynomial identities hold
integrally because the coordinate ring has no t-torsion. Alternatively,
[r]_t is the r-fold group sum by induction. This interpretation proves
commutation with every group-word construction used in graph potentials,
forest propagation, seams, and magic completion. The specialization and
tangent formulas follow from the j=1 term. In characteristic p all middle
binomial coefficients vanish, proving (5). A regular extension is unique:
the coordinate ring injects into its localization at t, so two pullbacks
of z agreeing there agree before localization. This proves the no-go with
its base-preserving and generic-power hypotheses. ∎

This is not a statement that flat specialization cannot carry any
Frobenius, nor that all Frobenius structures on G_a are impossible.

### The surviving characteristic-p repair: retain the base Frobenius

Over F_p[t], relative Frobenius is the base-preserving map

    G_t → G_(t^p),     z ↦ z^p,                          (6)

where the target is the pullback of the family along t↦t^p. Indeed

    (u⊕_t v)^p = u^p + v^p + t^p u^p v^p
               = u^p ⊕_(t^p) v^p.

Equivalently absolute Frobenius on the total characteristic-p scheme sends
(t,z)↦(t^p,z^p); it is not a morphism over the identity of the t-base. At
t=0 it gives the ordinary additive Frobenius and at t=1 it gives torus
Frobenius. On generic multiplicative coordinates, (6) sends
1+t^p z^p=(1+tz)^p, exactly as required.

Thus the distinction is between the SAME-base multiplication-by-p
endomorphism (4) and the BASE-TWISTED relative Frobenius (6). The base map
is indispensable information. No characteristic-zero integral commuting
lambda/Frobenius lifts of (6) are asserted. Ordinary square-root equations
r_i²=z_i are still not group doubling, and none of these maps is asserted
to preserve distinct positive square-valued magic entries.

## 3. Why the two other tempting connections do not improve the mathematics

### Parent forests are not graph-potential constraints on tournament paths

The incoming alpha=l parent component is a finite rooted tree, of height
at most 2^l−1−l. Applying Session59's graph-potential functor to that tree
is valid but tautological: a tree with v class-vertices gives
C_tree(A)≅A^(v−1), with no cycle constraints. The alpha=3 example has nine
class-vertices, so it produces the rank-eight split torus. This does not
identify class-vertices with tournament vertices or transfer canonical
owners, attachment responses, clone disagreement sets, Hamilton paths,
or blue/black predicates. A path and a star with nine vertices already
give the same abstract torus. An additional labelled construction would
have to specify the missing predicates before claiming a useful bridge.
The rank-transfer automata are not constructed by this observation.

### A good conic does not torify the square-root surfaces

The conic C is P¹ and genuinely toric. The Session59 root surfaces with at
least six distinct labels, over C and under the no-four-concurrent-lines
hypothesis, have a smooth projective resolution with p_g≥1. Thus no
nonempty open of those surfaces has a finite split-torus stratification:
a dense two-dimensional torus would make the function field rational,
contradicting p_g≥1. This is the existing Session59 birational obstruction,
not a new consequence of the incoming counting model. Conic powers,
finite-field dynamics, and smooth/flat graph degeneration do not remove
it. The obstruction does not cover every framework called F1 and does
not invalidate blue presentations which retain additive relations.

## 4. Exact checks and stopping boundary

`verify_bridge.py` checks the two bridges independently on the declared
finite prime-field universes and exact polynomial identities. The 7,950 gates passed, and the saved normal/optimized receipts agree
byte-for-byte; all gates use explicit exceptions. General proofs are above, not extrapolations from the checks. The parent
independently reviewed the characteristic-two Veronese charts and
Proposition2's polynomial existence, unit condition, localization uniqueness,
mod-p specialization, and base-twisted relative-Frobenius identity: PASS.
This review is distinct from the worker's bounded computational replay.

The integration yields one exact shared map (conic/cyclotomic powers) and
one precise obstruction with repair (same-base power specialization versus
base-twisted relative Frobenius). The parent-forest analogy adds no useful
new predicate-preserving result; the torification boundary remains intact.
No complete old-source audit, all-prime root-surface reduction theorem,
characteristic-zero lambda lift of the relative repair, or external
novelty determination was attempted. No repository files were edited and
no Git mutation was performed by this integration worker.
