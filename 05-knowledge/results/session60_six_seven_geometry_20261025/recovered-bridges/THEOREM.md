# Six-branch genus two is the exact bridge: root surfaces, pants, and Fano pairings

25 October 2026. **PROVED, with classical inputs explicitly cited; finite controls separate.**
This is a scoped extension of the root-surface hierarchy inherited from Session59,
not a new classification of Humbert–Edge curves, genus-two Jacobians, or Richelot
isogenies. The original square-point and plane six/seven questions remain open.
All geometric statements below are over **C**, unless another base is explicitly
specified. Branch labels and the chosen transverse line are retained throughout.

## 0. What really connects

There are three maps, not a match of numerals:

1. A six/seven-root surface has a transverse curve cover with exactly **one/seven
   six-branch genus-two quotients**. Their indexing characters are exactly those
   supporting the surface's one/seven holomorphic two-form eigenspaces.
2. The genus-two quotient turns the old pants cuff/seam interface into an
   explicit codec on its sixteen two-torsion classes. A branch pairing gives
   the precise four-element subgroup that cuff restriction forgets.
3. The old Fano/octonionic matching `{13,26,45}` is such a branch pairing. It
   therefore supplies a maximal isotropic subgroup and a **Richelot (2,2)-isogeny**
   on each corresponding genus-two Jacobian. A rational slice of the actual
   magic hourglass produces seven nonsingular explicit Richelot targets.

The repair is to transport the **marked pairing**, not the whole Paley or
PSL(2,7) action as an automorphism of a fixed complex surface. An exact
Möbius obstruction below demonstrates that distinction on the displayed slice.

## 1. Inheritance and entry conditions

The primary current input is Session59's
`toric-shared-roots/TORIFICATION_OBSTRUCTION.md`, Theorem T1. For distinct
spanning forms `L_0,...,L_(n−1)` in three variables, with no four concurrent
branch lines, the projective simultaneous-square-root surface `X_n` is an
integral normal complete intersection of `n−3` quadrics in `P^(n−1)`.
Its minimal resolution `Y_n` has `K=(n−6)H`, and for `n=6,7` has `p_g=1,7`.
The possible triple points give only A1 nodes on `X_n` and do not affect these
holomorphic-form counts. This packet uses that proved result; it does not
re-audit the entire surface hierarchy or drop the concurrency hypothesis.

The older exact mechanisms are:

- `sixes_and_sevens_20261006.md`, §3 and THM-4553: Paley labels `F_7`, the
  point-stripping ladder, and the three pairs `{13,26,45}` through vertex zero;
- `plane_cycles_pants_and_six_seven_20261011b.md` and
  `pants_boundary_holonomy_20261011b.md`, §3: a genus-two theta-pants marking
  has four compatible cuff vectors and four seam classes over each one;
- the current shape/ports packet's discipline: preserve the actual consumer
  and all marking data. No first-integral or quotient count substitutes for
  a scalar Jacobian, native arithmetic guard, or actual geometric lift.

Exact source paths, repository commit and content hashes are in `PROVENANCE.json`.
No assertion in this note promotes the old knight blocking conjecture, Lyapunov
order-six conjecture, S6 manuscript, plane colouring target, or Collatz target.

## 2. Restrict to a genuine curve, not a symbolic six-element set

Choose a line `ell ≅ P^1` in the base plane, avoiding every intersection of two
branch lines and not equal to a branch line. Write `l_i=L_i|ell`, whose zeroes
`b_i` are consequently distinct. Let `C_n` be the inverse image of `ell` in
`X_n`; it is smooth because the base line meets branch divisors transversely
at their smooth points and avoids surface nodes.

Its function field is

    C(P^1)(sqrt(l_1/l_0),...,sqrt(l_(n−1)/l_0)).

The square classes are independent: a nonempty product of these ratios has
odd valuation at a selected `b_i`, `i!=0`. Thus `C_n` is connected and its
deck group is

    G_n = F_2^n / <(1,...,1)>,  |G_n|=2^(n−1).

The inertia at `b_i` is the class `e_i`. Their sum is zero and no other
relation exists. Every branch index is two. Riemann–Hurwitz gives

    2g(C_n)−2 = 2^(n−1)(−2+n/2),
    g(C_n) = 1 + 2^(n−3)(n−4).                         (2.1)

In particular `g(C_6)=17` and `g(C_7)=49`. These are **not** genus-six and
genus-seven curves. Equivalently the curve is a complete intersection of
`n−2` diagonal quadrics in `P^(n−1)`, a classical Humbert–Edge curve of type
`n−1`. The indexing shift is explicit.

**Section-class guard.** A line in the **base** pulls back to `2H`, not `H`,
because `f^*O_(P^2)(1)=O_X(2)`. Our line avoids the nodes, so its smooth
inverse image on the resolution has class `2H`. Adjunction independently gives

    2g(C_n)−2 = (2H)·((n−6)H+2H) = 2(n−4)2^(n−3).

This is an extra diagonal **quadric** section, unlike a root-coordinate
hyperplane section of class `H`. In particular a general `H`-section of
the six-root K3 has genus five, while this base-line `2H`-section has genus
17; for the seven-root surface these two sectional genera are 17 and 49.

### Every even subset gives an actual quotient

Characters of `G_n` are even subsets `S` of the branch labels:

    chi_S(e_i) = (−1)^[i in S].

For nonempty even `S`, the quotient by `ker(chi_S)` is the double cover

    D_S: y^2 = product_(i in S) l_i,

where this equation uses the line bundle `O_(P^1)(|S|/2)`; a root at infinity
is retained, rather than silently discarded as a missing factor. Its genus is

    g(D_S)=|S|/2−1.                                    (2.2)

Pullback identifies `H^0(D_S,K)` with the `chi_S` eigenspace in
`H^0(C_n,K)`: differentials descending through `ker(chi_S)` are invariant
under that kernel; their invariant-under-all-`G_n` part is zero since
`P^1` has no holomorphic differentials. Conversely invariant holomorphic
differentials descend through a finite curve quotient. In a local chart
`t=z^e`, invariance makes the differential a multiple of `z^(e−1)dz`,
which is precisely pullback of a regular `dt`-differential.

These eigenspaces exhaust all differentials. Their dimensions sum to (2.1),
also directly by the even-binomial identity. Summing pullback maps on
Jacobians therefore gives an isogeny

    product_(S even, |S|>=4) Jac(D_S)  -->  Jac(C_n).    (2.3)

Indeed its differential is an isomorphism, so the image is an abelian
subvariety of the same dimension and its kernel is finite. This argument
claims an isogeny, **not an isomorphism or a product principal polarization**.
It does not claim the genus-two factors are simple or mutually nonisogenous.
The general decomposition is classical; see Auffarth–Lucchini Arteche–Rojas,
[Theorem 2.1 and Remark 1.2](https://arxiv.org/html/1905.12690v2).

| roots | curve genus | four-branch elliptic quotients | six-branch genus-two quotients |
|---:|---:|---:|---:|
| 6 | 17 | 15 | 1 |
| 7 | 49 | 35 | 7 |

Thus `17=15+2`, whereas `49=35+7*2`. This is a quotient construction,
not a claim about torus factors of the original surface.

### Why these are the same one/seven sectors as `p_g`

On a diagonal complete intersection, the residue numerator changes by
`product_i epsilon_i` under a coordinate-sign change; the defining quadrics
are unchanged. At `n=6`, the constant canonical numerator therefore has
character `chi_{all six}`. At `n=7`, the canonical numerators are the seven
coordinates `x_a`. The numerator `x_a` multiplies the residue character by
`epsilon_a`, giving

    character of x_a·Res = product_(i!=a) epsilon_i
                        = chi_(all seven except a).   (2.4)

Both formulas descend through simultaneous projective sign, since the
resulting supports have even size. Resolving nodes preserves these forms.
Consequently the size-six characters have two simultaneous geometric
readouts: a one-dimensional surface two-form sector and a two-dimensional
curve one-form sector arising from a genus-two double-cover quotient.
There is **no pullback of a two-form to a nonzero one-form** being asserted:
the bridge is the common character and the actual quotient maps.

## 3. Pants gives a geometric reconstruction, not only a dimension count

There are two useful levels of the pants bridge.

### 3.1 The whole root curve

Cut the punctured base sphere into `n−2` pants along a tree with `n` leaves.
The three cuffs of each pair of pants separate three nonempty disjoint sets
`A,B,C` partitioning the branch labels. Their inertia vectors are
`g_A=sum_(i in A)e_i`, and similarly for `B,C`.
These are distinct nonzero vectors, with `g_A+g_B+g_C=0`. They span a
four-element subgroup. Each connected lift of one pair of pants has degree
four, Euler characteristic `−4`, and six boundary circles (two over each
cuff), so it is a sphere with six holes.

There are `2^(n−3)` such pieces over each base pair of pants, and
`2^(n−2)` lifted circles over each internal cuff. Cap the lifts of the
original punctures using the ramified branch disks. Every resulting piece
still has genus zero. The connected gluing graph therefore has

    V=(n−2)2^(n−3), E=(n−3)2^(n−2),
    b_1=E−V+1=1+2^(n−3)(n−4).

At six roots this is `(V,E,b_1)=(32,48,17)`; at seven it is `(80,128,49)`.
The checker constructs the actual coset graph, including parallel edges,
and verifies connectedness, not merely the Euler subtraction. This is a
pants realization of the same curve genus; it does not make either curve
Euclidean. Both are already hyperbolic as compact Riemann surfaces.

### 3.2 A genus-two quotient and the exact lost two bits

Let `D` be a double cover with six labelled branch points. Its two-torsion
has the standard concrete symplectic model

    W = {v in F_2^6 : sum(v_i)=0}/<(1,...,1)>,
    B([A],[B])=|A intersect B| mod 2.                   (3.1)

Each nonzero class has a unique two-element representative `{i,j}` and
corresponds to the divisor class `[w_i−w_j]`. Complementary even subsets
give the same class. The pairing is well-defined, alternating and
nondegenerate. The divisor/Weil-pairing identification is recalled in
Castryck–Folsom–Hubrechts–Sutherland,
[§6, after Lemma 6](https://math.mit.edu/~drew/proclms104.pdf).

A perfect matching of the six branch points gives three disjoint branch
cuts on the sphere. Their lifts are the three cuffs in a theta-pants
decomposition of `D`. For example use pairs `12,34,56`; their mod-two
classes `a_1,a_2,a_3` satisfy `a_3=a_1+a_2`. Use the symplectic pairing to
regard `v in W` as the binary character `B(v,−)`. Its cuff vector is

    epsilon=(B(v,a_1), B(v,a_2), B(v,a_3)),
    epsilon_3=epsilon_1+epsilon_2.                     (3.2)

The kernel is `L={0,a_1,a_2,a_3}`: all four classes have zero cuff data.
A choice of dual seam classes `b_1,b_2` with `B(a_i,b_j)=delta_ij` and
`B(b_1,b_2)=0` gives the inverse

    v = tau_1 a_1 + epsilon_1 b_1
      + tau_2 a_2 + epsilon_2 b_2,
    tau_i=B(v,b_i).                                    (3.3)

This is exactly the inherited theta-pants four-bit codec, now with branch
points and Jacobian two-torsion as its input. The zero class and `a_1` have
identical cuff vectors but different global binary characters; the first
unramified double cover is disconnected, the second connected. Retaining
seam bits repairs that failure. The chosen seam transport is a marking,
not something recovered from cuff data alone.

Permuting branch labels acts faithfully on `W` and preserves `B`.
Since `|Sp_4(F_2)|=15*8*3*2=720=6!`, this realizes `S_6≅Sp_4(F_2)`.
This is a statement about marked two-torsion/monodromy, not 720 automorphisms
of a fixed generic genus-two curve. The checker exhausts all 720 actions,
all 15 matchings and all 16 classes for each codec.

## 4. The old Fano matching becomes a Richelot kernel

The recovered six-point set is `F_7^*`, with the old matching

    {1,3}, {2,6}, {4,5}.                               (4.1)

It comes from the three Fano lines through zero and is precisely the
pairing used by the old octonionic `J`; no integrability statement about
`S^6` is involved. In the model (3.1), the three pair classes are mutually
orthogonal, are nonzero, and sum to zero. Therefore

    L_Fano={0,[13],[26],[45]}

is a two-dimensional maximal isotropic subgroup of `Jac(D)[2]`.
The classical Richelot construction gives the principally polarized
quotient `Jac(D)/L_Fano`, by a degree-four `(2,2)`-isogeny.
Its target is either a genus-two Jacobian or a product of elliptic curves;
this distinction must be tested. A matching selects **one of fifteen**
possible kernels, not a unique unmarked choice. See Bruin–Doerksen,
[§4](https://arxiv.org/html/0902.3480v2).

For seven branch labels `F_7`, omit `a` to obtain `D_a`, and use

    M_a = {{a+q,a+3q}: q in {1,2,4}}.                  (4.2)

This gives seven specified genus-two Jacobians with seven specified
maximal isotropic kernels. The affine Borel group `x->u x+b`,
`u in {1,2,4}`, preserves this construction as labelled data:
`M_a -> M_(u a+b)`. At `a=0`, `x->2x` cycles the three nonzero kernel
classes. On the whole four-dimensional space `W` it fixes no nonzero
class: every nonzero class has a unique two-element representative, and
no two-element subset is a union of its two three-cycles.

**Typed conclusion.** The old Fano matching has a genuine image in abelian
geometry. The octonionic multiplication law, tournament Hamiltonian paths,
knight moves and Collatz affine carries are not transported by this map.
The Borel equivariance is a relabelling rule between marked configurations,
not a group of holomorphic symmetries of a fixed arrangement or slice.

## 5. An actual rational magic slice with seven nonsingular targets

Keep the old hourglass forms and choose this explicit label order:

    (L_0,...,L_6)=(E−x−z, E, E+x, E−x, E+z, E−z, E+x+z).

Omitting label zero is the inherited six-root intermediate. On
`[E:x:z]=[s:t:2s+5t]`, the restricted forms are

    (−s−6t, s, s+t, s−t, 3s+5t, −s−5t, 3s+6t).

Their seven distinct branch points in `u=t/s` are

    (−1/6, infinity, −1, 1, −3/5, −1/5, −1/2).

The line actually avoids all intersections of the nine magic branch lines;
the other two restricted roots are `−1/4` and `−3/4`, also distinct from
these seven. Thus this is a legal transverse slice of the real source
arrangement, not an arbitrary sextuple substituted for it.

For each omitted label `a`, write its three quadratic factors using (4.2):
`F_1,F_2,F_3`. A factor involving infinity is allowed to have affine degree
one, retaining its homogeneous degree two. Let `Delta` be the determinant
of their coefficient rows in basis `(1,u,u^2)`. The seven values are

| omitted a | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| Delta | −162 | 474 | −537 | −105 | −60 | −129 | −68 |

All are nonzero. Set `B_i=F_j'F_k−F_j F_k'` cyclically. The product
`B_1 B_2 B_3` is a squarefree sextic in every case, so over C it defines the
Richelot target. The checker records all factors and all seven target
polynomials in `normal.json`. For example at `a=0` the source factors are

    1−u,   6u^2+9u+3,   −25u^2−20u−3,

and the target factors can be taken as

    105u^2+114u+33,   25u^2−50u−23,   6u^2−12u−12.

The bracket formula is the classical construction, as also recalled in
Castryck–Decru–Smith, [§3.2](https://nutmic2019.imj-prg.fr/confpapers/HashRichelot.pdf).
Constant normalization is harmless over C. Over Q a quadratic twist can
matter: this packet does **not** assert that an arbitrarily rescaled target
model is the Q-defined Richelot quotient. The kernel and the source are
rational; an arithmetic implementation should use the normalized factors
`G_i=Delta^(−1)B_i` and the cited correspondence convention, retaining its
specified base-field twist. No rational square point is inferred.

### Two hostile controls and their strongest repair

**A pairing need not give a genus-two target.** The smooth genus-two curve
`y^2=(u^2−1)(u^2−4)(u^2−9)` has three independent-looking branch pairs, but
the coefficient determinant of this splitting is zero. The Richelot target
is a product, not a smooth genus-two Jacobian. Keep the determinant/target
smoothness gate. This shows why an abstract matching alone does not certify
the stronger target category.

**A label symmetry need not lift to the fixed curve.** In the `a=0` example,
the old label permutation `sigma:x->2x mod7` requires
`infinity -> −1 -> −3/5 -> infinity`. The unique Möbius transformation
with these first three images is

    M(u)=(-25u−19)/(25u+15).

But it sends branch point `u=1` (label 3) to `−11/10`, not `−1/2` (label 6),
which `sigma` requires. Therefore this particular label action does not
come from an automorphism of that genus-two curve. Every genus-two
automorphism descends through its canonical hyperelliptic map, so the
Möbius test is necessary. The repaired statement is the exact symplectic
label action and the equivariant family of kernels, not an automorphism
of the fixed slice. This does not assert that every other permutation fails.

## 6. Concept board and stopping boundary

| Source -> target | Actual map / preserved datum | Necessary sidecar | Hostile / boundary |
|---|---|---|---|
| n-root surface -> transverse root curve | base-line pullback; same deck characters | chosen transverse line and branch labels | line through branch intersection changes curve |
| size-six character -> genus-two quotient | quotient by character kernel | six actual branch points, not just support size | same support, varying cross-ratios |
| one/seven two-form sectors -> one/seven genus-two factors | common character indexing | no change of differential degree claimed | one-form pullback is not restriction of a two-form |
| six branch points -> genus-two pants codec | even subsets, Weil pairing, three branch cuts | seam transport and matching | zero cuffs conceal four global characters |
| Fano/octonionic matching -> Richelot quotient | three pair classes span a maximal isotropic subgroup | marking, polarization and determinant | smooth sextic can have split Richelot target |
| Paley Borel -> marked family action | affine relabelling of omitted point and pairing | full branch configuration | explicit order-three Möbius obstruction |
| curve quotients -> finite-field trace decomposition | direct character expansion over good odd primes | geometric curve, good reduction, base field | degree does not equal rational fiber cardinality |

The exact result is a new **repository bridge application** of classical
constructions. It gives no reduction from the integer magic-square problem
to a finite two-torsion calculation, no identification of the K3 surface
with a genus-two Jacobian or Kummer surface, and no automatic full
PSL(2,7), Fano-plane, or Paley action on the fixed complex source.
The old plane-colouring6/7 alternative is not decided by these covers.

## 7. Reproduction

Run, in this directory:

    python3 check.py > normal.json
    python3 -O check.py > optimized.json
    cmp normal.json optimized.json

The checker uses explicit exceptions, not optimization-sensitive assertions.
It checks character dimension identities for `4<=n<=9`, constructs the two
lifted pants graphs, exhausts the sixteen two-torsion classes and all720
label actions, checks all fifteen Lagrangians/codecs, checks Borel21
relabelling, and verifies every displayed Richelot factor, determinant,
squarefreeness condition and hostile. In addition, at all good primes
`5<=p<=101` for this slice it directly counts projective root fibers and
compares their trace with the sum of all even-subset double-cover traces.
These finite controls supplement the proofs; they are not a global point
search, an independent proof of the cited Richelot theorem, or a
verification of a base-field twist normalization.
