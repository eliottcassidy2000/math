# Exact gradients vanish, but tetrahedral opposition creates genuine graph `H^1`

**Status: elementary proof plus VERIFIED-EXACT typed synthesis; not a new
LRC/JC bridge theorem.**  This note uses audited
[THM-3494](../01-canon/theorems/THM-3494-weighted-lift-primitive-coordinate-discriminant-atlas.md),
[THM-3496](../01-canon/theorems/THM-3496-marked-graph-kummer-degree-square-and-finite-coefficient-frobenius-flux-extinction.md),
[THM-3487](../01-canon/theorems/THM-3487-two-twenty-four-state-fibonacci-bundles-cycle-type-obstruction.md),
and newly promoted
[THM-3479](../01-canon/theorems/THM-3479-literal-half-twist-relation-current-two-transplant-certificate.md).
[THM-3497](../01-canon/theorems/THM-3497-berggren-t4-calibration-languages-and-harmonic-splitting.md)
was promoted during this synthesis after an independent hostile audit, so its
three calibration quantifiers, two minimal languages, and harmonic splitting
are now proved dependencies rather than reserved evidence.  The independent
standard-library companion reconstructs every finite claim below without
importing any source companion.

## 1. Verdict and inheritance

THM-3479's audit is decisive: its thirteen owner-oriented edge differences
are `delta v`, hence lie in `B^1`; all six absolute cycle pairings vanish even
though the bridge, both weighted `K4` tree factors, and their product are
nonzero in every one of 72 labelled charts.  THM-3494 has the identical
gradient verdict multiplicatively: `g_ij=v_j/v_i` is exact even when some
edge square classes are nonzero.

There is nevertheless a lawful nonzero graph-cohomology construction.  It is
not the gradient.  It uses the self-duality of a tetrahedron:

```text
vertex potential --opposite face--> oriented face boundary --edge pairing-->
graph H^1 of the K4 one-skeleton.                              (1)
```

On THM-3494's square classes, (1) canonically retains the product of four
vertex volumes and is provably nonzero for the audited cubic views
`(w,x,y,z)`.  On each odd-characteristic `K4` wing of THM-3479, it is a full
three-dimensional isomorphism; the two-wing map fires in all 72 endpoint
role charts.  On the complete role graph it loses exactly the leaf direction,
and the original hub-to-leaf bridge is the one scalar sidecar restoring it.

This does not produce a physical current or a map to the D5 Kummer line.  The
opposite-face construction pays for its escape from exactness with a
two-dimensional face/opposition structure, an orientation sign in odd
characteristic, and a chosen edge pairing.  No such sidecar is part of the
distinguished THM-3494 or THM-3479 edge gradient.

The inheritance board is:

| object | native datum | audited status | surviving nonzero object |
|---|---|---|---|
| THM-3494 primitive views | multiplicative `K4` gradient `v_j/v_i` | exact, graph class zero | opposite-face square-class image; nonzero for cubic `(w,x,y,z)` |
| THM-3496 D5 chart | marked `H^1_graph(C7;F13)` seam | nonzero | deck character, then marked Kummer line |
| THM-3487 Fibonacci repair | `H^1_graph(C6;V4)` seam | nonzero | full `V4` deck defect |
| THM-3497 Berggren fixed drift | `eta in H^1_grp(D4;F2)` plus the 192-state automaton | PROVED nonzero character; language stricter | one period bit, not full seam or acceptance language |
| THM-3479 endpoint roles | additive role-graph gradient | exact, graph class zero | oriented opposite-face class on both `K4` wings; bridge sidecar |
| U_full ancestry | two-marginal tensor observer | noncohomological kernel | common-base joint lift still missing |

The closest proved mechanism remains THM-3496's marked seam/deck-defect map.
The canonical hostile remains THM-2538's two-by-two checkerboard.  The
corrected near miss is the nonzero-chart/zero-class distinction in THM-3479.
The least-used sidecar is tetrahedral vertex-to-opposite-face duality.

## 2. Five different meanings of `H^1`-shaped notation

For a connected graph `G=(V,E)` and additive coefficient group `A`,

```text
delta:C^0(G;A) -> C^1(G;A),
B^1=im(delta),
H^1_graph(G;A)=C^1(G;A)/B^1.                         (2)
```

An exact edge cochain can be nonzero on every edge while representing zero
in (2).  The number

```text
beta_1(G)=|E|-|V|+1                                  (3)
```

is only the dimension of a cycle space over a field.  It is not a class.
The present profiles are

```text
                  dim B^1       beta_1
K4                    3             3
C7                    6             1
role graph            7             6.               (4)
```

Here `K4` means its graph one-skeleton.  Filling its four triangles gives the
tetrahedral two-sphere and changes the relevant complex.  The faces below are
an external opposition template; they are not added as graph 2-cells.

For a group `Gamma` acting trivially on `A`,

```text
H^1_grp(Gamma;A)=Hom(Gamma,A).                        (5)
```

This is group cohomology, not graph cohomology of a state graph.  A nonzero
group character can become a graph coboundary on a Cayley graph because the
two complexes allow different degree-zero cochains.

Finally,

```text
0 -> ker(e_L) tensor ker(e_R)
  -> U tensor V
  -> U x_k V -> 0                                    (6)
```

is a marginal-observer sequence.  Its left term is a tensor kernel, not
cocycles modulo coboundaries.  A checkerboard picture does not manufacture a
graph complex.

## 3. The nonzero map: tetrahedral opposition

Order the vertices of a `K4` as `(0,1,2,3)`.  Let `F_i` be the oriented
triangular face opposite vertex `i`, and put

```text
c_i=(-1)^i partial(F_i),
Omega(v_0,v_1,v_2,v_3)=sum_i v_i c_i.                 (7)
```

Regard the oriented edge chains in (7) as graph cochains through the standard
oriented-edge pairing.  Then

```text
partial c_i=0,                 sum_i c_i=0,
pi_* Omega = sign(pi) Omega pi_*       for pi in S4. (8)
```

Thus constants are killed.  The companion computes an integral split minor
of determinant `-16`.  Consequently, over every field of characteristic not
two,

```text
C^0(K4;k)/k*1 --[Omega]--> H^1_graph(K4;k)            (9)
```

is an isomorphism.  This is a discrete Hodge-star-like transform, not
functorial passage of `delta v` to cohomology.  Equation (8) also states the
price of odd-characteristic covariance: tetrahedral orientation is a sign
local system.

Characteristic two is different and exactly matches THM-3494's square-class
carrier.  The cut and cycle spaces cease to be complementary.  The image of
(7) in `H^1_graph(K4;F2)` has rank one; all four `[c_i]` become the same
nonzero class `kappa`.  For any `F2`-vector coefficient group `A`,

```text
[Omega(a_0,a_1,a_2,a_3)]
   =kappa tensor (a_0+a_1+a_2+a_3).                  (10)
```

Apply (10) to `A=K^*/K^{*2}` and vertex volumes `v_i`.  In multiplicative
notation,

```text
[Omega_sq(v)]=kappa tensor [v_0 v_1 v_2 v_3].         (11)
```

This is invariant under changing the common reference basis: multiplying
every `v_i` by `u` changes the product by `u^4`, a square.  Signs and
tetrahedral orientation also disappear in characteristic two.  Equation
(11) is therefore a canonical `S4`-equivariant one-bit slice of the four
volume square classes.  It loses the other two graph-`H^1` coordinates and
all individual index divisors.

### The audited cubic atlas makes (11) genuinely nonzero

THM-3494's degree-three row supplies four primitive views `(w,x,y,z)`.  Up to
irrelevant signs and squares, their change-of-power-basis volumes are

```text
v_w=1,
v_x=C^3 f,                 f=9P-27Q-2,
v_y=C^3(27Q),
v_z=C^6 I_z(P,Q).                                      (12)
```

The exact printed `I_z` core satisfies

```text
I_z(2/9,0)=32/2187 !=0,                               (13)
```

and `(2/9,0)` lies on the irreducible linear divisor `f=0`.  Hence `f` divides
neither `27Q` nor `I_z`; its valuation in

```text
v_w v_x v_y v_z=C^12 f(27Q)I_z                       (14)
```

is exactly one.  The product in (14) is not a square in `Q(P,Q,C)`.  Equations
(11)--(14) prove a genuinely nonzero class

```text
[Omega_sq(v_w,v_x,v_y,v_z)]
    !=0 in H^1_graph(K4;K^*/K^{*2}).                 (15)
```

This is the strongest lawful THM-3494 transport currently available.  It
does not contradict the theorem: the six ratio edges `v_j/v_i` still form an
exact multiplicative coboundary.  The new map starts from the four vertex
volumes and uses opposite faces; it does not reclassify the gradient.

## 4. THM-3479: two full `H^1` wings plus one lost bridge

The audited role graph consists of two `K4`s sharing the hub `u5`, together
with the leaf edge `u5--u7`.  Let `Omega_L` and `Omega_R` be (7) on the two
ordered wings and set

```text
Omega_role=Omega_L+Omega_R,
b(v)=v(u7)-v(u5).                                    (16)
```

The exact integral split minor is `-256`.  Therefore, over any field of
characteristic not two,

```text
C^0(G_role;k)/k*1
   --([Omega_role],b)--> H^1_graph(G_role;k) direct_sum k
```

is an isomorphism.  Equivalently,

```text
rank(Omega_role)=6,
ker(Omega_role)=<constants, leaf-only>,
rank(Omega_role,b)=7.                                 (17)
```

The leaf-only potential is the sharp hostile: its `H^1` image is zero while
its bridge scalar is nonzero.  Thus the bridge is not another cycle
coordinate; it is exactly the tree direction omitted by the two tetrahedral
wings.

Using THM-3479's five audited endpoint role values at the certified primitive
embedding, the independent probe finds

```text
both K4 H^1 wing classes nonzero: 72/72 charts,
bridge b nonzero:                   72/72 charts.      (18)
```

The cut/cycle split also has full rank at the certified prime, and nonzero
reduction proves exact cyclotomic nonvanishing.  The digest of all 72
`(Omega_role,b)` rows is

```text
4141229331e5ff9f9d1ca303c205d7a47a3edc478e0630ae4aa7988072886e8e. (19)
```

Again, the audited current is still `delta v` and still zero in absolute
graph `H^1`.  Equations (16)--(19) define a different algebraic observable.
It is neither THM-3479's matrix-tree determinant nor a physical edge current.

There is also a sharp canonicity obstruction.  The role graph has 72
automorphisms.  Exact calculation gives

```text
H^1_graph(G_role;F13)^Aut(G_role)=0.                  (20)
```

Because `13` does not divide `72`, Maschke semisimplicity gives the same zero
multiplicity for invariant linear functionals.  The D5 Kummer line carries
no supplied role-graph action, so there is no source-native
`Aut(G_role)`-equivariant scalar map from (20) to that line.  Selecting a wing
orientation, a cycle functional, or a sign-twisted target can produce a
scalar, but that selection is precisely a missing sidecar.

## 5. The separate lawful cross-complex map: cyclic transgression

Let `A` be killed by a prime `p`, orient `C_m`, and mark the degree-`p` cover

```text
pi:C_(mp)->C_m,              tau(j)=j+m.              (21)
```

Seam sum gives `H^1_graph(C_m;A)~=A`.  If `delta h=pi^*g`, then

```text
tau^*h-h=sum_(edges of C_m) g.                        (22)
```

This produces the marked isomorphism

```text
T_(m,p,A):H^1_graph(C_m;A) -> H^1_grp(C_p;A),
[g] |-> (tau |-> seam(g)).                            (23)
```

Degree `p` kills the seam and degree `p+1` restores it.  The exact instances
are

```text
(m,p,A)=(7,13,F13): THM-3496, kill 13 / restore 14;
(m,p,A)=(6,2,V4):   THM-3487, kill 2 / restore 3.     (24)
```

THM-3496 continues the first line, after its explicit markings and residue
gate, to the Kummer line.  The second line ends at a full `V4` deck defect
unless another character is chosen.  Equations (23)--(24) are a common
schema, not a coefficient map between the two rows.

## 6. Berggren characters and automata retain less than their seams

The affine fixed-drift group is `D4=<G,T>` of order eight, with

```text
|Hom(D4,F2)|=4,             dim H^1_grp(D4;F2)=2.     (25)
```

The period character has `eta(G)=eta(T)=1`; on branch words it is length
parity and creates the `96+96` state split.  But restriction of all four
characters to the translation `V4={0,p,q,r}` gives only

```text
(0,0,0), (0,0,0), (1,1,0), (1,1,0) on (p,q,r).       (26)
```

Every character kills the nonzero seam `r=G^2`, although THM-3487 proves
that `r` repairs the six-cycle just as `p` and `q` do.  The missing sidecar is
the full translation holonomy.

Nor are the variable/fixed word languages character kernels.  Independent
enumeration gives

```text
variable: accepted 4 < character kernel 6,
fixed:    accepted 34 < sign kernel 96.               (27)
```

The one-letter word `C` is the sharp hostile in both cases: it passes the
relevant scalar gate but fails the required cycle-type match.  Promoted
THM-3497 proves these strict containments and, separately, that `eta` controls
the alternating fixed-drift chambers with level densities `1/6` and `3/16`.
That harmonic role does not turn `eta` into the acceptance predicate.

## 7. Sharp transport no-go after the positive result

For graph maps `f` and additive coefficient maps `phi`,

```text
phi_* f^*(delta v)=delta(phi_*f^*v).                  (28)
```

Thus ordinary graph pullback, coefficient change, and passage to cohomology
send both distinguished gradients to zero.  The exact companion checks all
`4^7=16,384` maps from `C7` into the primitive-view `K4`; every pulled-back
seam telescopes.

The positive map `Omega` does not evade (28) by clever coefficients.  It
changes the operation: vertices are first sent to opposite 2-faces and only
their boundaries are retained on the 1-skeleton.  This is the minimal new
sidecar.

Three further obstructions prevent a forced D5/Berggren bridge:

1. THM-3494's coefficient group `K^*/K^{*2}` is 2-torsion, so every additive
   map from it to `F13` is zero.  A chosen square-class character can land in
   `F2`, but no source-native map identifies the resulting `K4` class with
   the Berggren period character or acceptance automaton.
2. THM-3479's certified split prime proves nonvanishing; it is not a semantic
   coefficient object.  Reducing the exact characteristic-zero endpoint
   values to `F13` requires an integral lattice, a prime/residue marking, and
   denominator control.  No such realization is present, and (20) still
   forbids an unmarked scalar functional.
3. `Hom_Add(F2,F13)=Hom_Add(F13,F2)=0`.  Hence the D5 and Berggren lines share
   degree-multiplication grammar, not a direct natural transformation.

The sharp no-go is therefore:

> No composite made only from graph maps, additive coefficient maps, and
> ordinary passage to cohomology sends the THM-3494 or THM-3479 distinguished
> exact gradient to the nonzero THM-3496 Kummer class or to a nonzero
> Berggren seam/character.

The sharp survivor is equally precise:

> Tetrahedral opposition sends the THM-3494 cubic four-view volume packet to
> a nonzero square-class graph `H^1`, and sends every audited THM-3479 endpoint
> role chart to two nonzero odd-characteristic `K4` classes.  These are new
> transformed observables, not transported gradients.

After arbitrarily choosing a role orientation, a cycle functional, a
coefficient realization in `F13`, and all THM-3496 markings, one could compose
linear arrows to its Kummer line.  That composite would encode the choices,
not a theorem relating the source mathematics.  A lawful physical bridge
still needs a source-defined phase/holonomy observable and a proof that the
endpoint current factors through it.

## 8. The ancestry boundary remains a tensor kernel

THM-3479's endpoint worker sums left and right atoms separately before
forming `phase*ax*by`; no shared ancestry key survives.  In the full
`13 by 13` categorical model of (6),

```text
joint dimension=169,
marginal image dimension=25,
mixed kernel dimension=144=(13-1)^2.                 (29)
```

The smallest hostile is

```text
[ 1 -1 ]
[-1  1 ],                                             (30)
```

with zero row and column marginals.  This is a model of the universal loss
type, not a dimension claim about U_full's unexposed atom spaces.  Neither
`Omega_role` nor the bridge scalar supplies a section of (6): both are formed
after the marginal endpoint values already exist.  The missing sidecar is a
joint common-base lift formed before marginalization.

## 9. Final map ledger

| source -> target | verdict | information lost | required sidecar |
|---|---|---|---|
| THM-3494 four volume classes -> `H^1_graph(K4;K^*/K^{*2})` | lawful via `Omega_sq`; nonzero for cubic `(w,x,y,z)` | two of three graph coordinates, individual index divisors | opposite-face structure; already canonical mod 2 |
| THM-3494 ratio gradient -> any graph `H^1` by ordinary functoriality | zero | all loop periods | a genuinely non-gradient operation/local system |
| THM-3479 endpoint potentials -> role graph `H^1` | lawful via oriented `Omega_role`; both wings nonzero in 72/72 | leaf bridge and physical semantics | two wing orientations and edge pairing; bridge restores algebraic data |
| THM-3479 edge gradient -> absolute graph `H^1` | zero | all loop periods | source-defined phase/holonomy |
| role `H^1` -> unmarked D5 Kummer line | no canonical nonzero map | automorphism/orientation and coefficient data | cycle functional, residue realization, D5 markings |
| `H^1_graph(C7;F13)` -> deck character -> Kummer line | lawful and marked (THM-3496) | physical word amplitude and ancestry | physical realization remains absent |
| `H^1_graph(C6;V4)` -> `H^1_grp(C2;V4)` | lawful marked instance | none of the full seam | no cross-domain claim |
| full `V4` seam -> fixed-drift `eta` | one-bit quotient | kills nonzero `r` | full translation holonomy |
| U_full marginals -> common-ancestry current | not determined | joint ancestry | common-base joint lift |
| D5 `F13` line <-> Berggren `F2` line | zero additive Hom | all nonzero coefficients | a new mixed-primary realization object |

## 10. Independent exact probe

Run

```bash
python -B 04-computation/typed_h1_coboundary_character_marginal_synthesis_probe_20260816.py
python -B -O 04-computation/typed_h1_coboundary_character_marginal_synthesis_probe_20260816.py
```

Both runs reproduce

```text
05-knowledge/results/typed_h1_coboundary_character_marginal_synthesis_probe_20260816.out
```

The companion verifies the graph profiles; the integral split indices
`16/256`; all 24 tetrahedral covariance identities; the characteristic-two
rank-one collapse; the THM-3494 divisor witness (13); all 72 THM-3479
opposite-face/bridge rows; the zero invariant role-`H^1` line; both cyclic
transgressions; every affine `D4` character; all 48 variable raw pairs and
192 fixed states; the full marginal rank; the cross-primary Hom obstruction;
and all 16,384 primitive-view-to-cycle gradient pullbacks.  Normal and
optimized outputs agree exactly.  LF SHA-256 values are recorded in the
matching results index entry after replay.
