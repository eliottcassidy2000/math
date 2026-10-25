# Seven K3 deletions and the doubled cube inside the central K3

**Status: PROVED**, in the explicit characteristic-zero scope below;
[independent proof and replay audit: PASS](../audits/ROOT_DELETION_AUDIT.md).
Exact computations are
**FINITE-EXACT** supplements. The
elementary proofs make no historical-priority claim. No assertion about
distinct positive integer squares or LRC follows.

## 1. Inheritance and live board

The inherited theorem is T1 in
[Session59, TORIFICATION_OBSTRUCTION.md](../../session59_f1_geometry_and_information_20261025/toric-shared-roots/TORIFICATION_OBSTRUCTION.md),
with its independent
[TORIFICATION_AUDIT.md](../../session59_f1_geometry_and_information_20261025/positive-audit/TORIFICATION_AUDIT.md).
It proves that the seven-root surface has only A1 singularities, its minimal
resolution has `(p_g,K^2,e)=(7,16,80)`, and every six-root deletion resolves
to a K3. The new operation is **delete one specified root, then retain the
branch divisor and its incidence**, rather than compare the two invariant
lists alone.

The hostile is a branch node: the smooth resolved seven-root surface is
usually **not a finite double cover of the minimal K3**. The correct finite
quotient is a blowup of that K3. The near miss is applying a smooth-double-cover
formula directly to the singular branch curve. The least-used sidecar is the
set of paired restrictions on the deleted line; it records both component
splitting and the cube directions.

**Incoming recovery (live `origin/main`, commit `f7041582a8`).** The concurrent
[Session60 deck/K3 packet](../../session60_six_seven_geometry_20261025/deck-k3/DECK_CHARACTERS_AND_K3_QUOTIENTS.md),
especially its sections 2 and 5, independently establishes the seven deletion
maps, canonical sign twist, and the distinction between the finite quotient
and the minimal K3. These are recovered/shared mechanisms here. Its
[Kummer/Gale packet](../../session60_six_seven_geometry_20261025/deck-k3/KUMMER_GALE_BRIDGE.md)
also identifies the central K3 as a classical Kummer surface. The present
addition resolves the finite quotient explicitly and computes its reducible
branch curves, their doubled-cube intersection lattice, and the seven
nonconstant Walsh modes. No novelty claim is attached to the shared results.

| Lane / live object | Question and decisive control |
|---|---|
| Anchor: root deletion | Which surface is the actual smooth quotient? Resolve `w^2=st` locally. |
| Branch normalization | Does genus 17 conceal reducible curves? Compute squareclass rank on the deleted line. |
| Niche: intersection graph | Are the eight central components just labels? Compute their actual intersections. |
| Canonical forms | Does `7=1+6` arise from a morphism? Pull back K3 forms and calculate the cover eigenspaces. |
| Wildcard: seven supports | Which Boolean characters survive in the K3 divisor lattice? Diagonalize the intersection matrix, keeping null relations. |

## 2. General seven-to-six theorem

Work over an algebraically closed field of characteristic zero. Let seven
distinct lines `L_0,...,L_6` in `P^2` have no four concurrent. Let `X` be the
surface obtained by taking their projective square roots, and let `Y` be its
minimal resolution. Fix `i`, delete root `i`, and write `X_i` for the resulting
six-root surface, `V_i -> X_i` for its minimal K3 resolution, and `H` for the
pullback of its hyperplane divisor. In particular `H^2=8`, `h^0(H)=6`.
Let `t=t_i` count the triple points incident with the deleted line; `0<=t<=3`.

**Proposition D1 (the actual smooth quotient).**

1. Forgetting root `i` defines a finite flat double cover `X -> X_i`,
   with algebra `O_(X_i) + O_(X_i)(-H)` and branch divisor `B_i in |2H|`.
2. On `V_i`, that branch divisor is disjoint from the exceptional curves
   over the nodes of `X_i`. It has precisely `r=8t` ordinary nodes and no
   other singularities.
3. Blow up those `r` branch nodes: `sigma:S_i -> V_i`. If `E_1,...,E_r`
   are its exceptional curves, the strict transform `B'_i` is smooth
   (possibly disconnected), and

       B'_i ~ 2M,   M = sigma^*H - sum E_a.

   Its smooth double cover `f:Y -> S_i` is the minimal resolution of `X`.
   Thus `S_i` is the quotient by the lifted root-sign involution `tau_i`.

**Proof.** Any six of the lines span the three-dimensional linear-form
space: otherwise they would all be concurrent. Hence `L_i` is a linear
combination of the other forms. In the six root coordinates its pullback is
a section `q_i` of `O(2)`, and the seventh root satisfies `x_i^2=q_i`.
The homogeneous rank-two extension is exactly the stated double-cover
algebra. There is no base point for the coordinate projection: all six
retained roots zero would force all base coordinates, hence the omitted
root, to vanish.

A node of `X_i` lies over three retained lines. The deleted line misses that
base point by the no-fourfold hypothesis. Therefore the branch misses every
such node, and it also misses its exceptional curve on `V_i`. The double
cover is locally etale there, producing two copies of the same A1 resolution.

At an intersection of the deleted line with just one retained line, the
branch is smooth in the root chart. At a triple consisting of the deleted
line and two retained lines, the six-root chart is smooth with coordinates
`u,v`, and its branch is `alpha u^2+beta v^2=0`, with `alpha beta !=0`.
This is an ordinary node. The other four root coordinates are nonzero;
their signs modulo simultaneous sign give `2^4/2=8` points. Each incident
triple therefore creates eight branch nodes, and these exhaust the possibilities.

At a branch node use coordinates making the equation `w^2=st`. Blow up
the base, putting `t=s a` in one chart. The normalized pullback has
`w=s b`, `b^2=a`, hence is smooth with coordinates `(s,b)`. Its map to the
original cover is `(s,b) -> (s,s b^2,s b)`. This is the usual A1 resolution;
the other chart completes the exceptional `P^1`. Since the branch has
multiplicity two at the blown-up point, its strict transform subtracts
`2E_a`, giving `B'_i~2M`. The local computations cover every singularity,
so the smooth cover is exactly `Y`. QED.

The exceptional curve of this local A1 resolution is **not fixed pointwise**
by `tau_i`: in the displayed chart `b -> -b`. Its two fixed points are
already on the strict ramification curves. Thus no isolated fixed points
or extra fixed exceptional components must be added to the fixed locus.

## 3. The genus-17 branch curve splits into a cube

**Proposition D2 (components, genera, and intersections).** The branch
divisor on `V_i` has `2^t` smooth irreducible components, indexed by a
Boolean `t`-cube. Each has genus

    g_t = 1 + 2^(3-t)(2-t).

Two different components meet if and only if their labels differ in exactly
one bit. In that case they meet transversely at `2^(4-t)` distinct points.
No three components meet. The cases are:

| `t` | Components | Genus of each | Nodes of the union | Dual multigraph |
|---:|---:|---:|---:|---|
| 0 | 1 | 17 | 0 | one vertex |
| 1 | 2 | 5 | 8 | one edge, multiplicity 8 |
| 2 | 4 | 1 | 16 | square, edge multiplicity 4 |
| 3 | 8 | 0 | 24 | cube, edge multiplicity 2 |

**Proof.** Restrict the other six linear forms to the deleted line
`L_i=P^1`. Each incident triple gives a pair of proportional restrictions.
No restriction vanishes identically and no three have the same zero.
There are therefore `d=6-t` distinct zeros, with exactly `t` paired zeros.

Over the algebraically closed field, each proportional pair of squared
root coordinates factors into two possible linear relations. These `t`
independent choices label `2^t` components. More intrinsically, in
`k(P^1)^*/k(P^1)^{*2}`, the five root ratios generate a space of rank `d-1`:
their divisor parity vectors are the differences between the `d` zero
labels and one reference label. Distinct valuation coordinates prove
independence. Each component therefore has map degree `2^(d-1)` to `P^1`.
At each of the `d` branch points the inertia has order two. Riemann--Hurwitz
gives

    2g_t-2 = -2*2^(d-1) + d*2^(d-2)
            = 2^(d-2)(d-4).

This is the displayed formula. Locally at a paired zero the two linear
branches differ only in that pair's sign choice, so there are no
self-intersections of a component. Away from the paired zeros the curve
was already smooth. Thus the components themselves are smooth.

At a paired zero, fixing the other `t-1` sign choices selects one edge of
the cube. The `8` nodes above that zero distribute among `2^(t-1)` such
edges, equally by independent sign changes. Hence each edge contains
`2^(4-t)` nodes. The local ordinary-node equation proves transversality
and forbids threefold component intersections. QED.

The arithmetic genus of the whole divisor is always `17`, since
`B_i~2H`, `K_(V_i)=0`, and `H^2=8`. The independent component calculation
gives the same answer:

    p_a(B_i)=2^t*g_t + 8t - 2^t + 1 = 17.

Here `8t-2^t+1` is the first Betti number of the connected dual multigraph.
Thus the lost geometric genus is accounted for exactly by additional
components and their gluing. There is no hidden smoothing assumption.

## 4. Intersection lattice and its seven nonconstant modes

Let `C_epsilon` denote the components. On the K3, adjunction gives
`C_epsilon^2=2g_t-2`. Consequently their intersection matrix is

    M_t = 2^(4-t) [(2-t) I + A(Q_t)].

For a subset `S` of the `t` cube directions, the Walsh vector
`v_S(epsilon)=(-1)^(sum_(j in S) epsilon_j)` is an eigenvector with

    M_t v_S = 2^(5-t)(1-|S|) v_S.                 (D3)

This follows by summing the `t` coordinate flips: the cube adjacency
eigenvalue is `t-2|S|`. The `t` zero modes are actual divisor relations,
not merely numerical accidents. For the paired restricted forms indexed
by direction `j`, the equations `x_a-c_j x_b=0` and `x_a+c_j x_b=0`
are hyperplane sections of `X_i`. Their divisors are respectively the two
coordinate half-cubes. Hence

    sum_(epsilon_j=0) C_epsilon ~ H
      ~ sum_(epsilon_j=1) C_epsilon.              (D4)

Each component occurs with multiplicity one: near a generic branch point
the other linear factor is nonzero and the branch is reduced. These are
`t` linearly independent relations. The intersection matrix has nullity
exactly `t`, so the component classes span a rank `2^t-t` rational
subspace of `NS(V_i) tensor Q`, with signature `(1,2^t-t-1)`. No claim
about integral saturation or the complete Picard lattice is made.

At `t=3`, the eight rational curves have matrix `-2I+2A(Q_3)` and spectrum

    4 (multiplicity 1), 0 (3), -4 (3), -8 (1).

The seven **nonconstant Walsh labels** are precisely the seven nonempty
subsets of three directions. The three singleton labels select the three
actual hyperplane relations (D4); the three pair labels and full-support
label select the four negative directions. This is a typed bridge to the
seven-support model of
[reflection_bridge/seven.md](../../reflection_bridge_20261025/seven.md):
the support labels index a spectral basis for an actual K3 intersection
matrix. They are not seven curves, seven rank-zero points, or seven
independent divisor classes. The exact three singleton supports that are
rank-zero points in the pointed-monoid `P^2` convention here become the
null relations. Restoring the eight component labels and the marked three
cube directions restores the geometric incidence that the seven labels
alone do not contain.

The positive mode is the empty Walsh support, corresponding to the sum
`B_i=2H`. It cannot be discarded if one wants the whole intersection form.
Contracting the two opposite component labels into a single label would
again give the earlier Boolean seven-set, but would not preserve the curve
intersection graph. We make no geometric contraction-identification claim.

## 5. The explicit hourglass: one central deletion and six peripheral ones

Use the seven forms

    E, E+x, E-x, E+z, E-z, E+x+z, E-x-z.

There are three triple points, each on `E=0`; each peripheral line lies on
exactly one triple. Therefore the central deletion has `t=3`, and all six
peripheral deletions have `t=1`.

For the central K3 write the six roots as `A_+,A_-,B_+,B_-,C_+,C_-`.
On the branch line `E=0`, its eight components have explicit equations

    A_- = i (-1)^e1 A_+,
    B_- = i (-1)^e2 B_+,
    C_- = i (-1)^e3 C_+,
    C_+^2 = A_+^2 + B_+^2,             e in {0,1}^3.

Each is a smooth plane conic. Changing just one sign makes the corresponding
coordinate zero, leaving exactly two intersection points; changing two
signs forces two coordinates, hence all three, to vanish and gives no
projective point. This directly recovers the doubled cube independently
of the general cover calculation.

| Deleted root | Minimal K3's nodes before resolution | Smooth branch components | Smooth quotient of `Y` | Fixed Euler characteristic |
|---|---:|---|---|---:|
| `E` | 0 | eight rational curves | `Bl_24(K3)` | 16 |
| any of the six peripheral roots | 8 | two genus-5 curves | `Bl_8(K3)` | -16 |

The number `8` in the first branch row is intrinsic: it is the three
independent pair-sign choices. The number `7` counts the seven different
coordinate deletions. Their K3 quotients are indexed, not asserted to be
seven pairwise nonisomorphic abstract surfaces.

## 6. Why every deletion gives the same seven canonical sections

On `S_i`, `K_(S_i)=sum E_a` and `M=sigma^*H-sum E_a`, so

    K_(S_i)+M=sigma^*H,
    K_Y=f^*(sigma^*H),       K_Y^2=2H^2=16.

The cover algebra, or its local differential `ds wedge dt / w`, gives

    f_*omega_Y = omega_(S_i) + (omega_(S_i) tensor M).

The invariant summand is `omega_(S_i)`; it consists of the pullback of the
unique K3 two-form. The anti-invariant summand contributes `h^0(H)=6`.
Thus

    H^0(Y,K_Y) = invariant dimension 1 + anti-invariant dimension 6.

This is a morphism-level explanation of `7=1+6`. It is compatible with the
global canonical-character description in the parent packet: the invariant
canonical line for deletion `i` is the line indexed by the complementary
six-element support. The action on an adjunction coordinate section has
the residue determinant twist; it is not just the naive sign of that
coordinate.

The branch count is a further check. Its normalization has Euler characteristic

    e(Fix(tau_i))=2^t(2-2g_t)=16t-32.

The finite smooth quotient therefore satisfies

    e(Y)=2 e(S_i)-e(Fix(tau_i))
        =2(24+8t)-(16t-32)=80.

The holomorphic Euler characteristic can also be computed without a smooth
branch shortcut:

    chi(O_Y)=2 chi(O_(S_i)) + M(M+K_(S_i))/2
             =4+H^2/2=8.

This agrees with `q=0`, `p_g=7`, and Noether's identity. The branch incidence
changes the quotient's blowup count and the fixed curves, while the canonical
invariants remain constant because the exceptional contributions cancel.

## 7. Connection contract, failure boundary, and reproduction

| Source -> target | Map / correspondence | Preserved | Loss and sidecar |
|---|---|---|---|
| Seven-root surface -> six-root K3 | forget a specified root, then resolve | function-field degree two; K3 two-form as invariant canonical line | finite quotient requires branch-node blowups |
| Deleted line -> branch components | proportional restriction pairs -> independent sign choices | squareclass rank, genera, component number | individual branch locations and coordinate scalars need the restricted forms |
| Branch curves -> cube | component sign labels, intersections as edges | complete intersection multigraph with multiplicities | embedding and divisor classes need `H` and the equations |
| Cube -> seven supports at `t=3` | nonconstant Walsh characters | exact spectral eigenvalues and null relations | empty mode, component basis, and marked directions must be retained |

The cheapest hostile is one triple on the deleted line: the branch then has
eight nodes, so the smooth quotient is a blowup, and its branch has two
genus-5 components rather than a smooth genus-17 curve. The central hourglass
provides the strongest hostile: eight rational components. The strongest
surviving uniform statement is arithmetic genus 17 plus the exact cube formula.

The no-fourfold hypothesis is essential. Four concurrent lines allow three
restrictions to coincide at one point; independent paired-zero splitting
and the node resolution no longer apply. Characteristic two destroys both
the sign-cover and the local quadratic calculations. Over nonclosed fields,
components can be permuted by Galois and need not be individually rational.

Run from the repository root:

```text
python -B 05-knowledge/results/six_seven_geometry_20261025/root_deletion/check_deletion.py
python -B -O 05-knowledge/results/six_seven_geometry_20261025/root_deletion/check_deletion.py
```

The checker uses every one of the 36 seven-line subsets of the nine magic
lines and each of its seven deletions (252 labelled cases), comparing exact
rational incidence against independently computed restriction squareclasses.
It checks the whole cube spectra and component arithmetic for `t=0,1,2,3`,
and explicitly intersects the eight central conics over `F_17`, where `i`
and all the two-point intersection fibers split. The finite-field check is
only a control of the displayed equations, not a replacement for the
characteristic-zero proof. Explicit exceptions remain active under `-O`.
