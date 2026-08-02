---
id: THM-3270
title: "Gamma(3) peripheral-cusp fixed-sheet transport and holonomy obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For G=PSL2(Z)=C2*C3, N=barGamma(3),
  and T=SR, the four N-cusps N\G/<T> are equivariantly identical to
  P1(F3), hence to the four sheets of the tetrahedral A4 action.  The cusp
  represented by g carries the N-conjugacy class of g<T^3>g^-1; in the
  Bass-Serre quotient this is the oriented boundary of the face opposite
  the fixed sheet g(infinity).  The four peripheral H1 classes sum to zero
  and any three form an integral basis.  Conjugation transports their
  N-classes independently of the chosen modular lift.  Therefore a
  peripheral class is an exact minimal combinatorial sidecar for the local
  fixed sheet used by THM-3230.  It gives a coherent section on a connected
  transport groupoid iff quotient holonomy fixes a sheet: iff the A4
  holonomy lies in one C3 point stabilizer (or the S4 holonomy lies in one
  S3 point stabilizer after peripheral orientation is forgotten).  A kernel
  loop never obstructs the mark, while a nonzero V4 loop or two distinct
  pure-C3 stabilizers do.  This constructs no geometric/rational section,
  cross-place norm gluing, Keller cofactor, C3/S4 exclusion, or JC(2)
  theorem.
audit: >
  An independent audit rederived the cusp double-coset bijection, parabolic
  normalizer gate, integral tetrahedral face basis, A4/S4 orientation rule,
  lift-independent conjugacy transport, subgroup fixed-point criterion,
  V4 and two-C3 hostiles, and twelve-flag completion.  Fresh normal and
  optimized runs byte-match the stored transcript and declared hashes.
source: jc-global-mark-gamma3-2026-08-02
depends_on:
  - THM-3141-quartic-v4-modular-congruence-shadow-and-gamma3-sidecar-boundary
  - THM-3145-bass-serre-two-three-tree-and-tetrahedral-congruence-quotient
  - THM-3072-a4-flag-three-c2-tomography-and-edge-cycle-cospan
  - THM-3230-marked-c3-trace-centered-norm-and-terminal-prefactor-recovery
related:
  - THM-3037-cusp-braid-s4-lift-dichotomy-and-common-sheet-owner-boundary
  - THM-3067-tetrahedral-modular-two-three-flag-quotient-and-origin-loss
  - THM-3095-marked-s3-affine-lift-half-face-and-oriented-tetrahedral-frame
script: 04-computation/jc_gamma3_peripheral_cusp_sheet_transport_thm3270.py
output: 05-knowledge/results/jc_gamma3_peripheral_cusp_sheet_transport_thm3270.out
script_sha256: 0db5f6def1d8cbe5917279d56a3b892f7f420f75728ca4865e3633a61c8834be
output_sha256: 561273d03e8eadc1eb4af08100b40a3927edaa71ed8284a6663a22112bb6deb6
hash_basis: LF-normalized bytes
---

# THM-3270 -- the lost level-three cycle carries a sheet exactly when it is peripheral

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Inheritance and exact answer

[THM-3141](THM-3141-quartic-v4-modular-congruence-shadow-and-gamma3-sidecar-boundary.md)
shows that the tetrahedral four-point quotient forgets the level-three word
`T^3`.  [THM-3145](THM-3145-bass-serre-two-three-tree-and-tetrahedral-congruence-quotient.md)
identifies the lost kernel with the rank-three fundamental group of the
subdivided tetrahedron.  [THM-3230](THM-3230-marked-c3-trace-centered-norm-and-terminal-prefactor-recovery.md)
needs one of the four locally fixed sheets before its centered norm can
recover the terminal cubeclass.  The question is whether the lost
`Gamma(3)` coordinate can transport that sheet.

The answer has two parts.

```text
arbitrary element or scalar invariant of Gamma(3):  no sheet is selected;
peripheral cusp class of Gamma(3):                   exactly one sheet;
global section:                                      iff holonomy fixes it.
                                                                    (1)
```

Thus the kernel does contain the right four-valued carrier, but only in its
distinguished peripheral subset.  It does not erase the global holonomy
obstruction.

## 2. Four cusps are the four sheets

Put

```text
G=PSL_2(Z)=<S,R | S^2=R^3=1>,
T=SR,
pi:G -> PSL_2(F_3)=A_4,
N=ker(pi)=barGamma(3),
X=P^1(F_3).                                                (2)
```

Use the compatible matrices of THM-3141, so `T` is the translation
`z |-> z+1`.  The stabilizer of `infinity` in `G` is `<T>`.  Hence the cusp
set of `N` is

```text
Cusps(N)=N\G/<T>.
```

Normality of `N` gives the exact double-coset identification

```text
N\G/<T>  ~=  A_4/<bar T>  ~=  P^1(F_3),
N g<T>   |-> pi(g)(infinity).                             (3)
```

The last map is well-defined because `N` reduces trivially and `<T>` fixes
`infinity`.  It is bijective because `bar T` has order three, so the right
side has `12/3=4` elements.  It is equivariant for the left `A_4` action.

The cusp represented by `g` has peripheral subgroup

```text
P_g=g<T^3>g^-1 <= N.                                     (4)
```

Replacing `g` by `n g T^k`, with `n in N`, changes `(4)` only by
`N`-conjugacy.  Conversely, the reduced-word centralizer/normalizer of the
primitive parabolic subgroup `<T^3>` in `PSL_2(Z)` is `<T>`.  Thus two
subgroups in `(4)` are `N`-conjugate exactly when their representatives lie
in the same double coset of `(3)`.  We may therefore write

```text
C_N={ [P_g]_N : g in G },       |C_N|=4.                 (5)
```

for the four **peripheral classes**.

The stabilizer in `A_4` of the class represented by `g` is

```text
pi(g)<bar T>pi(g)^-1.                                    (6)
```

It is the unique `C_3` subgroup fixing the sheet `pi(g)(infinity)`.  Hence
`(3)` gives an intrinsic equivariant bijection

```text
kappa:C_N -> X,
[g<T^3>g^-1]_N |-> the unique sheet fixed by
                         pi(g)<bar T>pi(g)^-1.            (7)
```

In particular, a pure-`C3` local inertia subgroup selects one and only one
peripheral class: the class whose stabilizer is that inertia subgroup.

This carrier is minimal among equivariant quotients.  An `A_4`-equivariant
quotient of `A_4/C_3` corresponds to an intermediate subgroup

```text
C_3 <= H <= A_4.                                         (8)
```

There is no such proper intermediate subgroup: `A_4` has no subgroup of
order six.  Therefore every equivariant quotient of the four peripheral
classes is either still four-valued or is a point.  Any invariant
identification of two sheet-cusps forgets all four.

## 3. The peripheral loops are the four opposite faces

Use the Bass--Serre tree from THM-3145.  Its `N` quotient is the once-
subdivided `K_4`; suppressing the six degree-two vertices preserves first
homology.  The word

```text
T^3=(SR)^3                                                 (9)
```

has six alternating Bass--Serre edges.  Modulo `N` it becomes one triangular
face of `K_4`.  Its quotient permutation `bar T` fixes one tetrahedral
vertex and cycles the other three, so that triangle is the face opposite
the fixed sheet.  Conjugating `(9)` gives the other three faces.

Choose the positive modular orientation and let

```text
c_x in H_1(N,Z) ~= H_1(K_4,Z),        x in X,             (10)
```

be the boundary of the oriented face opposite `x`.  The tetrahedral
orientation gives

```text
sum_(x in X) c_x=0,                                      (11)
```

and any three of the four `c_x` form an integral basis of the rank-three
cycle lattice.  One proof chooses a spanning tree of `K_4`: in the three
chord coordinates, every three displayed face boundaries have determinant
`+1` or `-1`.  Consequently the four vectors are pairwise distinct even up
to sign.  A peripheral subgroup class, equivalently the primitive ray
`{c_x,-c_x}`, recovers `x`; generic elements of `H_1(N,Z)` do not.

This explains both the rank-three and the four-valued data in THM-3145:
the four cusp loops have one relation rather than one cusp having been lost.

There is also an exact parity boundary.  Reduction

```text
PGL_2(Z) -> PGL_2(F_3)=S_4                              (12)
```

normalizes `N` and permutes its four **unoriented** peripheral subgroup
classes.  A determinant-minus-one reflection reverses every face
orientation, so on `(10)` it sends

```text
c_x |-> -c_(a x).                                        (13)
```

Thus the subgroup/ray transports under full `S_4`, while a chosen positive
peripheral generator transports only under `A_4` unless one also records an
orientation bit.

## 4. Lift-independent transport and the exact section obstruction

Let `B` be a connected path groupoid and let

```text
rho:B -> A_4                                               (14)
```

be a sheet-transport cocycle.  For a path `gamma`, choose any modular lift
`g_gamma in G` of `rho(gamma)`.  Transport a peripheral class by

```text
[P]_N |-> [g_gamma P g_gamma^-1]_N.                       (15)
```

This is independent of the lift.  Indeed another lift has the form
`n g_gamma`, with `n in N`, and its answer differs from `(15)` by the inner
`N`-conjugation by `n`.  Via `kappa`, equation `(15)` is exactly ordinary
sheet transport:

```text
kappa(gamma.[P])=rho(gamma) kappa([P]).                   (16)
```

Therefore a family of local pure-`C3` fixed sheets is coherent along a path
exactly when its peripheral classes are coherent.  No hidden choice of
modular lift remains.

Fix a base object and let `H<=A_4` be the holonomy image of its loops.  A
global coherent peripheral mark exists exactly when `H` fixes one class.
By `(7)`, this is equivalent to

```text
there exists x in X with Hx=x
  iff H <= Stab_(A4)(x) ~= C_3.                           (17)
```

For full unoriented `S_4` transport the identical argument gives

```text
there exists a coherent peripheral ray
  iff H <= Stab_(S4)(x) ~= S_3.                           (18)
```

The local prescribed marks must of course select the same class under each
path; `(17)--(18)` are the loop obstruction once one base mark is fixed.

Three sharp boundary cases follow.

1. A loop lying in `N` acts by inner conjugation on `(5)`, hence fixes all
   four peripheral classes.  Level-three winding by itself never obstructs
   a sheet mark.
2. A nonidentity `V_4` translation has no fixed point on `X`, hence one such
   quotient-holonomy loop forbids a coherent mark.
3. Two distinct pure-`C3` point stabilizers generate `A_4`.  Thus fixed
   sheets attached to two different stabilizers cannot be a single global
   section after both meridians are transported to one base point.

The third statement is immediate from the subgroup lattice: two distinct
Sylow-three subgroups cannot lie in a proper subgroup of `A_4`.  It is the
smallest global hostile to mistaking four valid local marks for one global
mark.

## 5. Completion to the flag carrier

Let `D` be the three nonzero `V_4` translations, equivalently the three
perfect matchings of the tetrahedron.  Combining the peripheral owner with
one matching direction gives

```text
C_N x D  ->  X x D,
([P],d)  |-> (kappa([P]),d).                              (19)
```

This is an `A_4`-equivariant bijection onto the twelve flags of
[THM-3072](THM-3072-a4-flag-three-c2-tomography-and-edge-cycle-cospan.md),
and the diagonal `A_4` action is regular.  Hence the level-three peripheral
class supplies exactly the four-valued point/owner coordinate, while the
resolvent direction supplies the independent three-valued coordinate.

The distinction between point labels and aggregated tomography remains.
Equation `(19)` reconstructs one flag from two **joint point labels**.  It
does not reconstruct a signal from two separate quotient averages and does
not remove THM-3072's three-dimensional blind character sector.  Nor does a
cusp class choose the direction `d` or the cyclic root of the `C3` block.

## 6. Consequence for the marked centered norm

At a pure-`C3` quartic completion, THM-3230 needs only the fixed sheet
`alpha`, not a cyclic ordering of the three moving roots.  Once an actual
quartic root realizes the peripheral class `(7)`, its operation

```text
divide by (T-alpha), trace-center the cubic root, take its norm             (20)
```

is invariant under cycling the moving orbit.  The peripheral cusp class is
therefore the exact **combinatorial** amount of marking needed for that
local operation.  A coherent section satisfying `(17)` transports which
root must be used, and THM-3230 then recovers the local classes `[K]` and
`Lambda` on each supplied completion.

This does not yet glue the resulting local units or divisors.  In
particular, the theorem does not supply:

1. a rational or polynomial section of the quartic graph order;
2. inclusion of the singleton idempotent in the affine Zariski-main open;
3. compatible uniformizers, residue-field maps, or signed divisor transport
   between different Jelonek components;
4. the sheetwise cofactor/inverse-different class required by THM-3064;
5. triviality or forbidden support of the recovered `Lambda`; or
6. an exclusion of pure `C3`, `A4`, `S4`, a degree-four Keller map,
   `JC(2)`, or `DC(2)`.

For a connected cover with full `A_4` or `S_4` holonomy, `(17)--(18)` in
fact prohibit a global sheet section on the whole base.  Any successful
Keller use must therefore live on a smaller boundary groupoid with
point-stabilizer holonomy, or replace a global section by genuinely
branchwise data plus a gluing law.  The level-three cusp sidecar identifies
that obligation; it does not solve the geometric realization problem.

## 7. Exact companion

Run

```bash
python3 04-computation/jc_gamma3_peripheral_cusp_sheet_transport_thm3270.py
python3 -O 04-computation/jc_gamma3_peripheral_cusp_sheet_transport_thm3270.py
```

Both modes byte-match the stored transcript.  The dependency-free companion
checks:

1. the `A_4` modular image and the nontrivial matrix
   `T^3=((1,3),(0,1))`, which is the identity modulo three;
2. all four double cosets and their matching `C3` sheet stabilizers;
3. the four opposite-face chains, relation `(11)`, and all four unimodular
   three-face determinants;
4. all `48` oriented `A_4` transports and all `96` unoriented `S_4`
   transports, including the odd orientation reversal;
5. the regular twelve-flag completion `(19)`;
6. all ten subgroups of `A_4` and the fixed-section criterion `(17)`; and
7. all six pairs of distinct `C3` stabilizers, the `V_4` hostile, and the
   kernel-loop positive control.

Every truth-bearing test uses an explicit exception and remains live under
optimized execution.

**QED.**
