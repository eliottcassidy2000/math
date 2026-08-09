# One quadratic refinement, two radial gauges, and two graph scales

**Status:** research synthesis and frontier map.  The proved source is
[THM-3333](../01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md);
this reflection is not an additional proof source.  `LRC(14)` remains
**OPEN**, and none of the graph ranks below is asserted to be a physical
current.

## Outcome first

The long pass found a common object rather than another analogy:

```text
primitive lattice/torus slope [u]
        |
        | Gaussian square Phi(u)
        v
ordered unreduced Pythagorean null point P_u
        |                         |
        | Lorentz pairing         | triangular value
        v                         v
2 * intersection(u,v)^2          inradius / signed exradii
        |
        | radial rescaling P_u -> Q_u=P_u/||u||^4
        v
constant THM-2056 Kelvin/Lorentz cap.                         (1)
```

At integer scale, pairing `2` is exactly Farey adjacency.  At Kelvin-square
scale, cap `2/91^2` is exactly the existing sufficient LRC determinant gate.
The same triangular potential has edge polarization equal to the Lorentz
pairing and null-vertex value equal to the right-triangle radius.  For an
unscaled Euclid lift, graph topology appears at parameter rank `r`; every
positive integral right triple has side/hypotenuse rank `2r^2`.  Clique
filling kills both.

The deepest correction is equally simple: a sum-of-two-squares **value** is
not the carrier.  A represented value, with its ordered Euclid legs and raw
parity scale, is.

## Anchor / Niche / Wildcard portfolio

| Lane | Inherited target | What moved | Honest boundary |
|---|---|---|---|
| **Anchor — LRC(14)** | THM-2056 Farey determinant gate; THM-2596 Gram-owner state | Kelvin-square Lorentz cap and fixed-owner scalar state `(C,F,owner)` | owner switches, non-hull labels, phase, clocks, global exit remain absent; no row closed |
| **Niche — sums, triangles, triples** | sums of two squares, triangular identities, Pythagorean/Berggren parameter line | one triangular radius potential has Farey polarization; four signs give all in/exradii; polygonal hierarchy proved | norm values and unordered triples are lossy quotients |
| **Wildcard — graph topology** | THM-2975 graph/partial-cube boundary; THM-3145/3157 incidence versus Johnson carriers | torus intersection becomes Lorentz pairing; two-clique surgery gives ranks `r,2r^2` | flag completion is contractible, so no Hodge/LRC obstruction follows |

The niche and wildcard did not close the anchor, but they changed its state
representation: the LRC-facing output is no longer “Pythagorean analogy.”  It
is the exact cap (THM-3333 (44)) and exact mediant update (36).

## Inheritance pass

| Required inheritance item | Closest object found |
|---|---|
| proved mechanism | THM-2056: Kelvin denominator, owner-cone defect, finite acute Farey certificate |
| second proved mechanism | THM-2596: passive `GL_2(Z)` Gram-owner covariance and the endpoint-defect hostile |
| canonical hostile | THM-2975: a literal triangle-rich graph need not be the partial cube seen after controlled forgetting |
| corrected near miss | HYP-2941/MISTAKE-222: a shared triangular/Pascal ambient array is not a bridge between predicates |
| least-used relevant sidecar | THM-3145's subdivided-`K4` incidence partial cube, contrasted with THM-3157's Johnson octahedron |
| sign/origin sidecar | THM-2606: in/excenter labels form an affine `V_4` sign set, but Feuerbach chooses an origin rather than a group motion |

The convention mismatch was load-bearing.  THM-2596 uses `0<m<n` and
`(n^2-m^2,2mn,n^2+m^2)`; the new chart uses `m>n`.  The parameter swap has
determinant `-1`, so it reverses face orientation and conjugates branch
matrices.  The final theorem records this instead of silently identifying the
two formulas.

## Live concept board

The board was kept to six objects:

```text
1. primitive Gaussian/torus slope [u]
2. rank-one Lorentz null point Phi(u)
3. triangular/polygonal radius potential
4. Farey edge with its two incident faces
5. fixed-owner LRC state (C,F,owner)
6. graph pair versus clique-complex pair.                  (2)
```

The operation grid was:

| Object | Representation | Invariant | Native operation | Quotient/scale hostile |
|---|---|---|---|---|
| slope | primitive column / torus curve | `|det|`, intersection | Farey mediant, basis move | sign quotient loses oriented determinant |
| null point | ordered `(A,B,C)` / rank-one symmetric form | Lorentz form, content | ambient addition, `GL_2` symmetric square | primitive normalization changes edge weights |
| radius potential | `T(A)+T(B)-T(C)` | null value, Lorentz coboundary | sign flip, polarization | common polarization forgets in/excenter origin |
| face | `u,v,u+v,u-v` | determinant `4`, two face norms | Vieta mutation | unlabeled face pair retains only `|u dot v|` |
| LRC state | `(C,F_w,owner)` | exact fixed-owner defect | mediant update | owner change and phase absent |
| graph pair | missing biclique | relative graph `H_1` | clique filling | every relative class dies after filling |

Every positive result came from changing either the representation or the
operation, not from comparing two scalar sequences in isolation.

## Pull 1: the symmetric square explains the whole arithmetic carrier

For `X=(A,B,C)`, the matrix

```text
M(X)=[C+A  B; B  C-A]
```

has determinant `C^2-A^2-B^2`, while

```text
M(Phi(u))=2u u^T.                                         (3)
```

So Gaussian squaring, Euclid's Pythagorean parameterization, the Lorentz
cone, and the modular `Sym^2` representation are one object.  Polarizing the
determinant of the sum of two rank-one forms gives

```text
<Phi(u),Phi(v)>_L=2 det(u,v)^2.                            (4)
```

This answers why the square occurs.  The determinant is the algebraic
intersection of two torus slopes; Lorentz pairing forgets its sign and keeps
twice its square.  Minimal positive intersection is exactly the Farey edge.

The topology is therefore inherited, not discovered anew.  Projectivization
is the double-angle homeomorphism `RP^1 -> S^1`; ideal chords give the usual
Farey triangulation of the open disk and the dual trivalent tree.  Straight
segments in the unprojectivized cone are not a new surface.

The first quotient probe was decisive:

```text
65=8^2+1^2=7^2+4^2,
```

but the two representations have Lorentz pairings `2` and `32` against the
same endpoint.  Norm value is too coarse.  Ordered legs repair it.

The parity probe sharpened the repair.  Every Farey face uses all three
nonzero mod-two spinor colors, so its sorted raw Pythagorean content triple is
`(1,1,2)`.  Primitive normalization changes its sorted edge weights from
`(2,2,2)` to `(1,1,2)`.
If leg order is also forgotten, the adjacent spinors `(2,1)` and `(3,1)`
both become the same unordered `3-4-5` triangle.  This is a genuine graph
quotient collapse, not merely inconvenient notation.

## Pull 2: triangular numbers are a quadratic refinement, not a shared list

The exact mechanism is

```text
T(x+y)-T(x)-T(y)=xy.                                      (5)
```

For

```text
R(A,B,C)=T(A)+T(B)-T(C),
```

equation (5) becomes

```text
R(X)+R(Y)-R(X+Y)=<X,Y>_L.                                (6)
```

On a right-triangle null point, the quadratic bulk cancels and only the
linear boundary remains:

```text
R(a,b,c)=(a+b-c)/2=r.                                     (7)
```

This is the precise bridge that MISTAKE-222 demanded: source, target, map,
preserved predicate, and quotient loss are explicit.  It is not evidence
that every triangular sequence in the repository carries an LRC or Jacobian
predicate.

The parameter-level identity

```text
r=n(m-n)=T(m)-T(n)-T(m-n)                                 (8)
```

shows a descent from the Pythagorean side lengths back to the Euclid split.
Combinatorially, `r` is exactly the missing `K_(n,m-n)` cut between two
complete graphs meeting at one vertex.

The four sign classes made the connection more rigid.  Because
`T(-k)=T(k-1)`, changing a side sign toggles an edge count to the cycle rank
of the same complete graph.  The four signed defects are

```text
(r,-r_a,-r_b,r_c).                                        (9)
```

Their polarizations are identical.  Thus the edge scalar cannot choose an
in/excenter origin; the `V_4` sign label is a necessary sidecar.  THM-2606's
Feuerbach origin is compatible with this bookkeeping, but no circle motion
or quartic transfer was inferred.

The polygonal extension explains why triangles are special.  For the
`s`-gonal polynomial `P_s`,

```text
polarization on a Farey edge = 2(s-2),
value on a right triple       = (4-s)r.                   (10)
```

Squares keep a nonzero edge form but have zero null-vertex value.  Triangles
are the unique usual polygonal family with unit Lorentz polarization and
positive radius value.  This is an exact uniqueness statement inside this
family, not a claim that triangular numbers are universally privileged.

## Pull 3: the face has a Vieta mutation and a scale choice

For a Farey edge, let `x,y` be endpoint hypotenuses and `z,z'` the two
incident face hypotenuses.  The vector diamond yields

```text
z+z'=2(x+y),
zz'=(x-y)^2+4,                                             (11)
```

and each face triple lies on

```text
x^2+y^2+z^2-2xy-2xz-2yz=-4.                              (12)
```

This supplies a local scalar mutation `z -> 2(x+y)-z`.  It also gives the
acute norm-only repair

```text
h=sqrt(xy-1).                                              (13)
```

On a fixed THM-2056 owner cone, `(C,F_w)` at the two endpoints recovers both
the Gram matrix and the owner coordinate, and updates both components by the
same `2h`.  This is strictly stronger than the endpoint-`F` data refuted by
THM-2596.

The radial-scale pull then connected it directly to LRC.  Integer null points
make graph adjacency constant.  Kelvin-square null points

```text
Q_u=Phi(u)/C(u)^2                                          (14)
```

make the determinant gate constant:

```text
D(u)<=C(u)/91
iff max_i <Q_u,Phi(c_i)>_L<=2/91^2.                       (15)
```

The projective point is unchanged, but the norm `C(u)` is load-bearing.
Graph and gate become simple in different radial gauges.

## Pull 4: two graph scales survive, their flag topology does not

For two vertex subsets covering one complete graph, the graph pair has one
relative generator for every missing cross edge.  If the exclusive part
sizes are `p,q`, then

```text
H_1(X,Y;Z)=Z^(pq).                                        (16)
```

At the Euclid-parameter split, `(p,q)=(n,m-n)` and `pq=r`.  At the
side/hypotenuse split,

```text
(p,q)=(c-b,c-a),              pq=2r^2.                    (17)
```

These are honest included graph pairs, improving on a bare coincidence of
cycle-rank formulas.  But their clique complexes are a simplex versus two
simplices glued along a nonempty simplex.  Both are contractible and every
relative homology group vanishes.

This killed the first wildcard hope cleanly.  The graph rank is free on
deliberately omitted edges; it is not yet a natural Hodge obstruction.  To
be useful at the LRC frontier it would need a lawful map from missing cross
edges to typed proof obligations or transitions.  No such map was found.

## Concept-board update after all pulls

| Board object | What the session changed |
|---|---|
| slope | gained an exact torus-intersection interpretation and a rank-one symmetric form |
| null point | became the faithful represented-sum-of-two-squares carrier; raw content is explicit |
| triangular potential | unified edge index, inradius, exradii, edge/cycle toggles, and polygonal comparison |
| Farey face | gained a determinant-volume law and Vieta mutation; orientation needs a spin/face gauge |
| LRC state | compressed from a matrix pair to `(C,F,owner)` on a fixed acute cone; phase/global data still absent |
| graph pair | split into parameter rank `r` and triangle rank `2r^2`; flag completion supplies a decisive no-go |

The strongest cross-lane movement is from the third and fourth objects back
to the fifth: triangular and Pythagorean structure did not close LRC, but the
face norm did compress the exact owner-cone recursion.

## Exact hostiles and stopping reasons

| Temptation | Cheapest hostile | Strongest survivor |
|---|---|---|
| a sum-of-two-squares value labels a Farey vertex | the two norm-`65` representations pair `2` versus `32` | represented norm plus ordered legs |
| primitive Pythagorean triples carry the same graph | odd/odd normalization changes pairing `2 -> 1` | raw content can be restored from ordered leg parity |
| unordered triangles are enough | adjacent `(2,1),(3,1)` collapse to one `3-4-5` | compatible raw ordered spinor images |
| any null triples with pairing two form a Farey edge | `(4,3,5),(6,8,10)` pair to two without two raw lifts | theorem restricted to the `Phi` image |
| triangular polarization is spinor addition | `(8,16,18)!=(16,30,34)` | ambient polarization and spinor mediant are separate lawful operations |
| Pythagorean graph means Berggren tree | `(2,1)->(4,1)` has determinant two and pairing eight | Berggren branches are Farey prefix words |
| the face scalar restores intrinsic orientation | changing a spin lift swaps the two faces | oriented face plus spin/face label |
| graph rank is a Hodge obstruction | clique filling sends both ranks to zero | one-skeleton relative cycle invoice only |
| Lorentz cap closes LRC | cap is the already-sufficient THM-2056 gate | exact coordinate and fixed-owner recurrence |

## Frontier questions produced by the synthesis

### A. Fixed-owner scalar replay — highest-value next test

**OPEN.**  Rebuild the finite acute THM-2056 fan and replay every cone using
only

```text
(C(u),F_w(u)), (C(v),F_w(v)), owner-id.                    (18)
```

Positive control: every mediant inside one fixed owner cone must agree with
the full Gram-owner computation by (13).  Hostile control: deliberately cross
every owner boundary and tie; the compressed state must refuse transport
without a new owner label.  Success would make the arithmetic state smaller
and more scalar, but would still not attach phase or global exit.

### B. Norm-mutation dynamics with typed owners

**OPEN.**  The Vieta law (11) defines an exact norm mutation on each labeled
Farey edge.  Determine whether the THM-2056 finite residual is closed under a
finite set of such mutations after owner labels are attached.  The norm-`65`
hostile requires vertices to retain representation IDs; a norm-only Markov
tree is invalid.

### C. Natural meaning of the two graph ranks

**OPEN / currently negative signal.**  Search for a canonical action on the
missing biclique edges before assigning meaning to `r` or `2r^2`.  Candidate
vertices are proof obligations, boundary sections, or owner-switch events—not
generic cycle generators.  Any proposal must survive clique filling and
explain why one-skeleton cells, rather than filled triangles, are native.

### D. Parity channel versus radial scale

**OPEN.**  Every Farey face has one doubled raw Pythagorean vertex.  Test
whether THM-2632's selected mod-two owner channel predicts where the doubled
vertex lies relative to a signed THM-2056 owner cone.  Parity alone cannot
predict the defect sign; the decisive test must retain `(C,F,owner)`.

### E. Polygonal potentials as controls

**PROVED family, OPEN use.**  Equation (10) supplies a one-parameter control
family.  Because the square potential vanishes at every null vertex but not
on edges, it is a clean hostile against arguments that infer edge structure
from vertex defects.  Test whether this separates any proposed scalar
compiler before committing to triangular-specific language.

## Session status

- **PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED:** THM-3333's
  Gaussian/Lorentz/Farey graph isomorphism, two gauges, triangular and
  polygonal polarizations, in/exradius formulas, face mutation, fixed-owner
  scalar recursion, and two graph/clique scale laws.
- **REFUTED:** norm-value graph, normalized/unordered-triple graph,
  Berggren-edge identification, intrinsic orientation from the scalar,
  ambient-addition/mediant conflation, and flag/Hodge promotion.
- **OPEN:** LRC(14), owner-switch transport, phase/clock/global exit, a finite
  typed norm mutation, and any natural current carried by the graph ranks.

No meta-pattern was promoted from this session alone.  The reusable move is
already covered by the maintained patterns: type the shared object, expose
the quotient loss, restore the smallest sidecar, and apply a native operation
before compiling a scalar verdict.
