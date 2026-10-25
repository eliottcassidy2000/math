# Independent third-wave audit: root surfaces and finite torification

**Verdict: PASS**, in the stated characteristic-zero and finite-locally-closed-
torification scope, for `TORIFICATION_OBSTRUCTION.md`, TheoremT1 and
CorollaryT2, read at2026-10-25 15:43UTC. This audit is separate from both
the base first-checkpoint audit and the resolution-of-B audit. The general
surface theorem, six-root K3 threshold, hourglass and full-surface applications,
and exact labelled-subset census were checked; no broad F1 no-go is implied.

## 1. Integrality, flatness, complete intersection, and normality

Choose any three independent defining linear forms L_i on the parameter
vector space. Replacing each remaining L_j by its expression in that basis
gives one monic quadratic for its root coordinate. Therefore the homogeneous
root ring is free of rank2^(n−3) over the polynomial ring in the three basis
roots; the equations are a regular sequence. The projective root surface is
a complete intersection of n−3 quadrics, dimension2 and degree2^(n−3).

Separately, the projective coordinate-squaring map is finite locally free
of rank2^(n−1): a target affine coordinate chart has n−1 independent monic
square equations in normalized root coordinates. Its pullback to the
parameter plane retains that rank. These two different degrees are correct:
the map pulls back O_P²(1) to O_X(2), multiplying self-intersection by4.

In k(P²), the n−1 square classes L_j/L_0 are independent: for any selected
numerator line its valuation is odd in precisely that selected factor.
Distinctness of lines is essential here. The generic root algebra is thus a
field of degree2^(n−1). Each affine finite-flat algebra over a base chart
injects into that field by torsion-freeness. It follows that the surface is
integral, not merely generically irreducible or connected. The construction
and argument are geometric over an algebraically closed characteristic-zero
field. The smooth/A1 local analysis below proves normality; equivalently,
complete-intersection S2 and isolated singularities imply Serre R1+S2.

## 2. All singularities and the root-sign count

At any parameter point choose a nonzero line value and normalize its root.
All other unit radicands split into étale local sign choices. Zero radicands
are exactly the branch lines through that point. With zero, one, or two such
lines, local root coordinates give a smooth ring: at a double point one has
x=u²,y=v² for transverse parameters x,y.

For three distinct concurrent lines, their forms have a unique relation
with all coefficients nonzero; otherwise two lines would be proportional.
The completed local root ring is

    k[[u,v,w]]/(w²−αu²−βv²),  αβ≠0.

This is a nondegenerate ordinary double point in characteristic zero.
There are n−3 nonzero root coordinates at the point. Their2^(n−3) sign
choices have one free simultaneous sign identification in projective space,
leaving2^(n−4) distinct points, each an A1 node. The no-fourfold-intersection
hypothesis excludes all remaining base-point types. No branch or singular
points are omitted by normalizing a nonzero coordinate, since a projective
root point always has one.

An independent linear-algebra check also classifies all crossings. If C is
the matrix of relations among the L_i, then the root Jacobian is2C diag(r_i).
After deleting columns where r_i=0, its rank defect is the dimension of the
space of relations supported on the vanishing forms. This is0 for at most
two distinct lines, and1 at a triple point. Thus the arrangement calculation
really locates every singularity, not just a chosen family of nodes.

## 3. Canonical divisor, crepant resolution, and K3 threshold

Adjunction for the complete intersection givesω_X=O_X(n−6). Resolving a
node is crepant, not merely a resolution of some kind: blow up the ambient
smooth threefold at the node, whose canonical divisor gains2E; the strict
transform of the multiplicity-two surface loses2E. Adjunction cancels the
two terms. Consequently the smooth minimal resolution has

    K_Y=π*((n−6)H),   K_Y²=(n−6)²·2^(n−3).

When n≥7 this is the pullback of an ample divisor, hence nef and big. It is
not generally ample on Y, because the exceptional curves have intersection0
with it. The draft correctly makes the big/general-type claim instead.

Normality givesπ_*O_Y=O_X; the projection formula therefore identifies the
canonical sections. The complete-intersection Koszul resolution gives the
Hilbert series(1−z²)^(n−3)/(1−z)^n and the required equality between graded
pieces and global sections. Its degree n−6 coefficient is p_g. In particular
n=6,7,8,9 give1,7,31,111; n=4,5 have p_g=0. For the hourglass, seven linear
sections give p_g=7; for the full surface,165 cubics minus6·9 quadratic-
linear products give111, with first Koszul syzygies only in degree4.

At n=6, the canonical bundle on Y is trivial. Koszul vanishing gives
H¹(X,O_X)=0. Locally at an A1 resolution the exceptional ideal's successive
quotients are O_P¹(2m), with vanishing H¹. Formal functions thus gives
R¹π_*O_Y=0, so H¹(Y,O_Y)=0. The connected smooth projective surface Y is
therefore a K3 surface in the standard algebraic definition. This step is
necessary; trivial canonical bundle alone would not prove K3.

## 4. Actual magic arrangements and every coordinate subset

The nine forms E,E±x,E±z,E±(x+z),E±(x−z) have eight triple points:
four at E=0, and four at[1:±1:0],[1:0:±1]. The remaining twelve intersection
points are double; none has multiplicity4. For the horizontal hourglass's
seven forms there are exactly three triples at infinity and twelve doubles.
A positive finite count of nodes is not being extrapolated to this claim:
the rational cross products of every pair of form vectors were computed,
normalized, and all incident lines checked exactly.

Any subset with at least four lines still has no fourfold intersection. If
its forms had rank at most2, all those lines would be concurrent, contrary
to that property. Thus every4≤n≤9 coordinate subset satisfies every
hypothesis of TheoremT1. The resulting invariants include:

| Surface | Triple points | Nodes | K_Y² | p_g |
|---|---:|---:|---:|---:|
| shared five-root base | 2 | 4 | 4 | 0 |
| specified six-root intermediate | 2 | 8 | 0 | 1 |
| either seven-root hourglass | 3 | 24 | 16 | 7 |
| full nine-root surface | 8 | 256 | 576 | 111 |

The number of nodes is not constant among all subsets of a fixed size. An
independent enumeration of every labelled subset reproduces the manuscript's
full n=5,…,9 triple-point table exactly. Those are labelled-subset counts,
not isomorphism classes or symmetry orbits. All six-root subsets are K3,
regardless of their varying number of nodes. The theorem does not assert
that every five-root subset is toric; the known shared base merely supplies
a positive five-root example.

As an additional independent finite check, at F17 all needed±1 and±2
radicands split. The checker constructs24 distinct projective singular root
points for H and256 for M, verifies every root equation, and verifies
Jacobian corank1 at each. It also checks Jacobian ranks symbolically for all
ordinary, single-line, double and triple strata. Finite-field checks do not
replace the characteristic-zero incidence hypothesis: small characteristics
can change the line arrangement.

## 5. Finite torification implies rationality; positive p_g contradicts it

For an irreducible surface, a finite union of locally closed subsets of
dimension≤1 cannot cover all geometric points. Hence a finite torification
must contain a two-dimensional torus stratum. Its closure is the whole
surface, and a dense locally closed subset is open. The surface therefore
has function field k(u,v) and is rational. The same argument applies to any
nonempty open subset, which has the same function field.

Every n≥6 root surface here has a smooth projective resolution with p_g≥1.
Geometric genus is invariant under birational maps of smooth projective
surfaces and is0 on rational surfaces. This contradicts the necessary
rationality of a finite torification. Consequently neither the surface nor
any nonempty open subset can have a finite decomposition into locally closed
split tori over C. In particular removing branch/zero/collision divisors
does not evade the obstruction.

This is stronger than the earlier *relative* branch obstruction, but it is
not a theorem against all structures called F1. It excludes finite geometric
torifications and, in particular, ordinary split toric varieties. The additive
blue-presentation pushouts remain valid. No conclusion about positive distinct
integer-square existence follows from the complex birational obstruction.

## 6. Primary-source and reproducibility boundaries

The primary [Auel–Singer manuscript](https://arxiv.org/html/2609.09351v1)
was opened independently. The full nine-root surface's nodes and numerical
invariants are existing results there, not new discoveries of this session.
The uniform square-root-cover argument and its applications above were
checked directly rather than accepted from that attribution. The
[López Peña–Lorscheid torification source](https://arxiv.org/abs/0903.2173)
was also checked for the decomposition framework; the rationality implication
used here was proved explicitly above.

Run `python3 check_torification.py` and `python3 -O check_torification.py`.
The source imports no producer implementation. All checks use explicit
exceptions and both runs produce identical outputs. The exact gate count,
subset table, node data and invariants are in `torification_certificate.json`;
input and output hashes are in `torification_manifest.json`. A development
run had a Python type-normalization error before producing a certificate;
it was repaired, and only the successful current source/output pair is
approved. No mathematical discrepancy was found.

## Final amendment: Euler stratification (2026-10-25 15:46UTC)

The final manuscript's added Section5 is also approved. For n lines with
t_2 doubles and t_3 triples, the four parameter strata have Euler
characteristics3−2n+t_2+2t_3,2n−2t_2−3t_3,t_2,t_3. Their *reduced*
preimages have finite étale degrees2^(n−1),2^(n−2),2^(n−3),2^(n−4).
Using reduced preimages is necessary: the scheme-theoretic branch fibers
are nonreduced. Multiplying and adding, then adding one per resolved node,
yields2^(n−4)(n²−9n+24). This is80 for the hourglass resolution and768
for the full resolution, agreeing with Noether's formula and the independently
computed p_g,K². The singular hourglass has Euler characteristic56, not80.
No finite-field point-count polynomial follows from this complex calculation.
The successful final checker includes these identities and n=4,…,9 controls.
