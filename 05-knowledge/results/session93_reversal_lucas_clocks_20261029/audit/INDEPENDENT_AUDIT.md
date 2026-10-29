# Session93 independent audit: reflections, holonomy, and exact clocks

Source pin: `59cb22d0e93b6c364c4ec5d2954d9fe14ed468d8`.
Date: 2026-10-29. This lane independently wrote its checker; no producer
code is imported. Status distinguishes elementary proof review, imported
classification, and finite exact checks. No novelty or major-open-problem claim.

## Source and inheritance status

The initial source packet contains Session90's square-root transport, centralizer
recurrence and nonflat embedded W2 witness, but **not** the all-order full
centralizer generation theorem. The new generator formula arrived directly
from the current holonomy lane and was independently checked here. Current
black/blue definitions were read from `THM-643-blue-black-line-parity-structure.md`,
including its corrected pure-blue and fibre-cap warnings. Current guardrails
and bounded startup routes were read. Sparse missing files are not treated as
missing mathematics. No classification result was inferred from historical
numerology or an unverified group order.

Anchor: full centralizer generation. Niche: exact fixed-point exception scope.
Wildcard: fixed-clock trace observation versus an unmarked conjugate trace.
Hostile controls are Paley7 (centralizer nontrivial but holonomy trivial),
dihedral prime-degree actions (one fixed point), and modulus5 (ramified memory).

## 1. Prime-degree gap: accepted with strict endpoint and cited dependency

The correct candidate says: if G is transitive of odd prime degree p and an
involution has `1 < f < sqrt(p)` fixed points, then G contains A_p, except
possibly p=11 with G=PSL2(11) or M11 in their eleven-point actions.

The classification input was checked directly in Gareth Jones,
[*Primitive permutation groups containing a cycle*, Theorem1.2(k=0)](https://arxiv.org/pdf/1209.5169).
It is CFSG-dependent. A transitive prime-degree action is primitive and contains
a p-cycle by Cauchy. The listed alternatives are affine, projective semilinear,
PSL2(11), M11, M23. No extra full-cycle hypothesis is necessary.

Independent elimination review:

* An affine involution in odd characteristic is x↦−x+b and fixes exactly one
  point. Thus the strict lower endpoint is essential. Dihedral D_p disproves
  any inclusive endpoint claim at arbitrarily large p.
* Prime projective degree `(q^d−1)/(q−1)` forces d prime. For odd d, a linear
  involution has a fixed eigenspace of dimension at least `(d+1)/2`: in odd
  characteristic its scalar square is a square because d is odd; in even
  characteristic its nilpotent part has square zero. The resulting number of
  fixed projective points exceeds sqrt(p).
* For d=2, prime degree forces q even; a nonidentity linear involution fixes
  one projective point.
* A genuinely semilinear involution has q=s² and, after scalar normalization
  and descent, fixes the subgeometry PG(d−1,s). Its point count
  `(s^d−1)/(s−1)` exceeds the square root of `(s^(2d)−1)/(s²−1)`.
* The ATLAS natural-action table shows the sole involution class of
  [M23 fixes7 points](https://brauer.maths.qmul.ac.uk/Atlas/v3/permrep/M23G1-p23B0),
  exceeding sqrt23. The sole involution class of
  [M11 fixes3 points](https://brauer.maths.qmul.ac.uk/Atlas/v3/permrep/M11G1-p11B0),
  inside the gap. The eleven-point exception therefore cannot be deleted.

This is a consequence of the cited classification, not a new elementary
classification theorem. Finite Lucas-word computations do not replace that
input. The older explicit three-cycle certificate for degree29 remains a
separate, less classification-heavy proof via Jordan's theorem.

## 2. Recursive triangles: accepted all-order generation mechanism

Let A_k=Aut(W_k), alpha_k the digit-negating reversal and C_k its centralizer.
Use the inherited intrinsic-block proof of A_k=A_(k−1)^3⋊C3. Then
`C_k ≅ C_(k−1) × A_(k−1)`, with coordinates `(c,a,beta a beta)`.

The new lane supplies a literal triangle generator, not merely a normal
closure assertion. Let x rotate the three top blocks; let v be a zero-prefix
digit rotation in A_(k−1), of order3, inverted by beta; and y=(v,1,1).
With s=alpha, t=x²alpha, u=y²alpha, the square-root triangle is

    h(u,t,s)=(1,v,v^−1).

Independently checked coordinates and composition order through k4. Indeed
z=y²x has order9, so its unique odd square root is z^5; substitution in
`y^−1 sqrt(y²x) x` yields the displayed tuple. These v generate A_(k−1).
Inherited triangles in child0 generate C_(k−1), yielding all C_k by induction.
No inference from a single embedded nonflat witness is being made.

There are exactly k(k−1)/2 displayed generators, and this is minimal:
`(B³⋊C3)_ab ≅ B_ab×C3`, since abelianization identifies the three coordinate
copies with no further relation. Hence A_j has elementary abelian
abelianization of rank j, and C_k has elementary abelian abelianization of rank
`1+...+(k−1)`. Any generating set must have at least this many elements.

Finite certificates independently generate C2 of order3 from one triangle and
C3 of order243 from three triangles. The k4 checks validate formulas but do
not enumerate its full group. The proof, not finite extrapolation, gives all k.

## 3. Blue quotient: accepted, but named paths retain the error

For a mark s, write X_s for paths fixed by reverse∘s and C_s for its
centralizer. Two transporters s→t differ by a target centralizer element;
therefore the induced map `X_s/C_s → X_t/C_t` is independent of transporter
and composes strictly. Every square-root triangle error vanishes on this
quotient. A nonidentity automorphism fixes no named full ordered path, so the
same error is visible before quotienting whenever X_s is nonempty.

The native partner obtained by reversing the non-path arcs is equivariant:
`K_(gP)(T)=g K_P(T)`. Thus holonomy preserves its **isomorphism class**, not
necessarily its named adjacency matrix. Blue colour and class-level native
partners need less information than named port observations.

For any subgroup L_s generated by fixed-base triangle errors, each fibre of
`X_s/L_s → X_s/C_s` has exactly `[C_s:L_s]` elements, by the free C_s action.
Normality is unnecessary. Full generation is special to W_k, not universal:
Paley7 has C3 but trivial holonomy.

## 4. Arithmetic: accepted, with single trace/history distinction

In `(Z/m)[phi]/(phi²−phi−1)` let a=phi². The fixed-clock pair observer is

    (u,v) ↦ (2u+v, 3u+4v), determinant5.

Its Smith normal form is diag(1,5). Every nonempty fibre has gcd(m,5)
elements; future observations satisfy y_(k+2)=3y_(k+1)−y_k. Thus this pair
classifies entire future fixed-clock traces exactly.

Over an unramified odd prime, a **single** norm-one trace has inversion fibres
{x,x^−1}. Do not transfer that statement to fixed-clock trace histories:
the two-coordinate observer is injective there. At modulus5, writing
eta=phi−3 gives eta²=0 and a=−1+eta; norm-one states ±(1+b eta) have traces
`±2(−1)^k` independent of b. That fivefold lost coordinate is real.
The order of a divides p−(5/p), but divisibility is not exact order and
neither gives a prime-distribution theorem.

### Additional received extensions, independently reviewed

The multi-clock observer with common initial trace and strides d_i has kernel
size `gcd(m,5F_(2d_1),...,5F_(2d_r))`: imposing initial trace zero sets v=−2u,
and each remaining equation becomes `5F_(2d_i)u=0`. This gives a direct
independent derivation without invoking a large transition matrix. Standard
Fibonacci strong divisibility compresses the gcd further. Direct iteration
checks strides (2),(3),(2,3),(2,4),(3,5) for every modulus2..20. At modulus6,
strides2 and3 individually lose factors3 and2 but jointly recover all36 states.

The final Lucas lane's explicit semilinear descent formulas, Fibonacci gap
identity `L_(2m+1)−F_(m+2)^2=2F_m F_(m−1)`, and parity split were read and
accepted. Its new Jordan certificates were independently replayed by direct legal-word
enumeration, reconstruction of the labelled reflection, permutation composition
and binary powering: degrees29,199,521 yield the claimed isolated prime cycles
of lengths13,23,7. This second route certifies those finite group conclusions
without invoking the full classification. The universal conclusion still uses
the independently reviewed classification argument.

The final holonomy lane's port orbit count2k+1 follows from first-nonzero-digit
position/value: each corresponding suffix is transitive under the relevant
factor. Three singleton orbits give precisely three globally fixed vertices
for k≥1. This does not collapse joint port information.

The metagraph torsor extension `(Y×D)/C ≅ D` after choosing a path base is
accepted: the invariant is c^−1d when P=cP0. Separate quotients instead lose
relative alignment. Likewise relabelling an entire legal native path-word
preserves its lengths, colours and target properties. A centralizer element
fixes the unique reversal pole, so a prescribed-pole black-core certificate
keeps that named pole. This transports existing certificates; it creates no
new connectivity or geodesic theorem. A nonempty path-word action is free
because its first full path already has trivial stabilizer.

The Lucas lane's final full-group decoder bound is also accepted: for an
equivariant deterministic map from a finite G-set of decoder states to
tournaments on the named G-domain, a state stabilizer is a subgroup of the
odd-order tournament automorphism group. Every state orbit therefore has size
divisible by the full2-part of |G|. This gives lower bounds2^24,2^194,2^517
for G=A29,S199,A521 respectively. It is not a bound on merely guarded
carriers, short certificates, or non-equivariant chosen decoders.

### Final bounded port and conic review

Parent's `ports/PORT_ORBITS_AND_LCA.md` is accepted at its exact recursive
scope. Immutable first-nonzero tags identify independent A_j factors. Folding
the side2 suffix by alpha turns the paired action into the same A_j action
on both sides. For a tuple inside one factor, pairwise longest-common-prefix
lengths recover the occupied prefix partitions, while directed first-branch
arcs determine their partial C3 rotations. One occupied child gives a unique
rotation to its target; two or three are compatible exactly when pair
orientations agree. Recursing proves sufficiency, including repeated ports.

The resulting one-/ordered-two-port orbit counts are2k+1 and8k²+1. A base
must hit every independent bottom C3 support of every A_j factor. One vertex
in each bottom triple suffices, yielding `(3^(k−1)−1)/2` anchors. This is not
the number of group generators or an arbitrary-element encoding length.
The proposed compact transporter needs at most rj occupied-prefix rotations;
O(r²j) digit comparisons construct/check the pair invariant. Emitting a full
permutation still incurs its output size; the stated bound does not hide that.

Independent finite audit directly enumerates C3's243 permutations and all
one-, two-, and three-port tuples on27 vertices. Complete tuple signatures
coincide exactly with actual orbits. Every candidate base of size0..3 has a
nonidentity stabilizer; `(3,9,12,15)` has trivial stabilizer, certifying minimum4.
The hostile `(3,6)` versus `(3,8)` at order9 has different folded signatures,
despite equal individual tags and equal original inter-block orientation.

The arithmetic lane's final conic identity is accepted by polynomial expansion:

    t²+z²−L_(2d)tz = −5F_(2d)² N(x),
    t=Tr(x), z=Tr(a^d x).

Our checker independently verifies every state for moduli1..24 and strides1..7.
At stride1 with5 invertible, the pair-observer matrix is invertible; therefore
norm-one states correspond bijectively to the conic `t²+z²−3tz=−5`.
The recurrence preserves it. Modulo5 that equation degenerates to
`(t+z)²=0`; the earlier fivefold phase loss is not repaired by calling it a
conic. No prime, density, or tournament-centralizer claim is imported.

### Final editorial scope check: exceptional decoder minimum

The Lucas lane's final degree11 sharpening is accepted: for
G=PSL2(11), of order660, an equivariant tournament-decoder state stabilizer is
odd and has order dividing165. Its only possible order greater than55 is165.
Such a subgroup would have index4; the nontrivial coset action of the simple
group G would be faithful, embedding a group of order660 into S4, impossible.
Thus every decoder-state orbit has at least12 elements. The actual orbit of
the Paley11 tournament under the inherited degree11 G action has stabilizer55
and size12, attaining the bound. That final orbit enumeration is producer
checked, not rerun by this audit; the lower-bound proof was independently
reviewed. This exact minimum concerns nonempty finite full-G equivariant
carriers decoding tournaments on the named eleven-point domain, not arbitrary
guarded carriers or compressed certificates. No new independent tests were
added for this editorial check; the count remains62,872.

## 5. Independent executable evidence and exclusions

`check_independent.py` uses only the Python standard library, with checks
active under optimization. It gives **62,872 checks**, identical output in
normal and `python -O` runs:

| Object | H | fixed paths | centralizer | triangle subgroup | blue presentations |
|---|---:|---:|---:|---:|---:|
| Paley7 | 189 | 9 | 3 | 1 | 3 |
| W2 | 3159 | 33 | 3 | 3 | 11 |

It independently tests every arithmetic state for moduli2 through40, including
composite and ramified moduli; literal triangle identities through k4; exact
small generated subgroups; and all fixed-base triangles of Paley7/W2.
Finite Aut groups here use known explicit structural descriptions, not a fresh
all-S_n automorphism exhaustion. Source Session90 supplies that separate check.

No access or permission blocker occurred. No source/remote edits were made by
this lane. Parent publication remains separate. No CI, full prime-distribution,
all-tournament holonomy, or major-conjecture closure is claimed.
