# The hourglass has a Campedelli quotient: the Fano gauge kernel becomes a deck group

2026-10-25. **PROVED in the stated characteristic-zero scope**, with
[independent audit](../audits/FANO_AUDIT.md). The elementary
code identities and geometric deductions below have full arguments. The
universal 5,596-system census and the specific marking counts are
**FINITE-EXACT** with two enumeration paths. The general Campedelli
construction and its two triple-point singularity types are **CITED /
RECOVERED CLASSICAL**. No historical novelty claim, fundamental-group
calculation, or integer-square existence assertion is made.

## Inheritance and live board

The closest geometric mechanism is Session59's [Theorem T1, root surfaces](../../session59_f1_geometry_and_information_20261025/toric-shared-roots/TORIFICATION_OBSTRUCTION.md),
with its [independent audit](../../session59_f1_geometry_and_information_20261025/positive-audit/TORIFICATION_AUDIT.md):
seven distinct spanning forms, no four concurrent branch lines, and the
minimal resolution of the resulting diagonal-quadrics surface. The closest
code mechanism is the [Fano gauge/Hamming packet, section 1](../../arithmetic_braids2_20260917_fano_code.md):
the kernel of the seven-point-to-seven-line sign map is the simplex code.
These are the two inherited objects we actually identify.
Incoming commit `f7041582a8` independently develops the seven canonical
characters and K3 deletion/character correspondences in
`session60_six_seven_geometry_20261025/deck-k3/`; those are concurrent
shared mechanisms, not novelty claims of this code packet.

**Final incoming update:** checkpoints `c24c123d86` and `a144b255fc` arrived
before our close-out. Their independently audited
[seven-line theorem](../../session60_six_seven_geometry_20261025/seven-lines/RESULT.md)
also proves the16/12/2 census and goes further: the minimal models are
Campedelli, properly elliptic, and Enriques, respectively. The16 free choices
give exactly two complex isomorphism classes. These are now recovered results;
the general5596-system lemma and the shared regular-conic-action criterion
remain additional connections developed here.

The canonical hostile is a concurrent triple whose three binary labels sum
to zero. The corrected near miss is treating all quotient singularities as
crepant: the bad triple gives `1/4(1,1)`, not `A1`. The least-used sidecar is
the *choice of Fano marking on the seven branch coordinates*. It decides
whether the action is free. MISTAKE-642 in `01-canon/MISTAKES.md` warns that
an affine Hadamard gauge is not automatically a linear code; here the zero
word, XOR closure, incidence kernel, and weights are checked explicitly.

Anchor: root-cover quotient geometry. Niche: the recovered Fano sign kernel.
Wildcard: realize the supposedly troublesome Boolean reflection by moving
between two code markings. The six live concepts are:

| Concept | New operation or consequence | Retained sidecar |
|---|---|---|
| Seven canonical sections | restrict deck action to a rank-three subgroup | canonical determinant twist |
| Fano gauge kernel | use it as an actual projective sign group | point labels and linear origin |
| Branch arrangement | test freeness by triple incidence | actual rational/complex concurrence |
| Campedelli quotient | divide `K^2=16` by eight | free action, or discrepancies |
| Boolean `3+1+3` reflection | conjugate the code subgroup | two distinct Fano structures |
| Six-root K3 boundary | compare deletion with smaller deck quotients | quotient order and singularities |

## 1. One exact sequence joins the two packets

Work over `C`. Let `X` be the seven-root surface associated to distinct
spanning linear forms `L_1,...,L_7` on the plane, with no four concurrent
zero lines. Write `Y -> X` for its minimal resolution. Inherited T1 gives

```
X = intersection of four diagonal quadrics in P^6,
K_X=O_X(1),  p_g(Y)=7,  q(Y)=0,  K_Y^2=16,
G = F2^7 / <(1,...,1)>       (the full projective sign group).
```

Choose a bijection `i -> v_i` from the seven root coordinates to the seven
nonzero vectors of `V=F2^3`. Define

```
S = { (a dot v_i)_i : a in V }                 [7,3,4] simplex code,
C = { w in F2^7 : sum_i w_i v_i = 0 }          [7,4,3] Hamming code,
sigma(w) = sum_i w_i v_i.
```

Every nonzero simplex word has weight four, so `1` is not in `S` and `S`
embeds in `G`. Directly, `C=<1> direct-sum S`: both sides have dimension
four, and all pairwise coordinate sums of nonzero binary vectors vanish.
Consequently there is an exact sequence

```
0 -> S -> G -> V -> 0,       [w] -> sigma(w).                 (1)
```

The projective quotient is legitimate because `sum_i v_i=0`. The inertia
element flipping the i-th root maps to `v_i`. Thus `X/S -> P^2` is the
eight-sheeted cover with its seven branch lines labelled by the seven
nonzero residual deck elements. This is a *syndrome quotient of the actual
root-sign group*, not merely a common seven-element indexing set.

The earlier arithmetic packet defines

```
(Bw)(a) = sum_(a dot v_i=0) w_i,    a != 0.
```

Its formula `Bw = (sum_i w_i) + a dot sigma(w)` proves `ker B=S` and
`im B=C` in the chosen point/normal polarity. Therefore the sign changes
invisible to the inherited Fano line-gauge observer are exactly the group
by which we now quotient the geometric root surface. No octonion
multiplication or tournament path-count predicate is transported.

## 2. The geometric freeness test is precisely a Fano-incidence test

For `a!=0`, the diagonal action `(-1)^(a dot v_i)` has four minus signs
and three plus signs. A fixed projective point must have all its nonzero
coordinates in one sign block. On `X`, at most three coordinates vanish,
since zero coordinates correspond to branch lines through one plane point.
The four-coordinate vanishing option is excluded. The only possibility
is that exactly the three plus coordinates vanish. Their labels form the
Fano line

```
F_a = {v != 0 : a dot v=0}.
```

Hence:

**Freeness criterion.** `S` acts freely on `X`, and hence on `Y`, iff no
concurrent triple of branch lines is a Fano line under the chosen marking.

The lift to `Y` is canonical, since the minimal resolution of the nodes is
unique. A fixed point upstairs maps to a fixed point downstairs, so there
are no concealed exceptional fixed points in the free case. This checks
the actual resolution rather than only the smooth locus of `X`.

When the action is free, put `Z=Y/S`. It is smooth and the quotient is
etale of degree eight. Its canonical divisor pulls back to `K_Y`, giving
`K_Z^2=2`. The invariant part of `H^1(Y,O_Y)` is zero, so `q(Z)=0`.
Etale multiplicativity of the holomorphic Euler characteristic gives
`chi(O_Z)=(1+7)/8=1`, whence `p_g(Z)=0`. Moreover `K_Z` is nef and big
because its pullback is; thus `Z` is minimal of general type. It is a
numerical Campedelli surface. If the arrangement has `t` concurrent
triples, `X/S` has exactly `t` nodes, one orbit of eight original nodes
over each triple; `Z` resolves those nodes.

There is a second check on `p_g=0`. For arbitrary root signs `epsilon_i`,
the actual canonical action on the coordinate section `x_i` is
`(product_j epsilon_j) epsilon_i`, by the residue/adjunction determinant
factor. On `S` all words have even weight, so this becomes the seven
distinct nontrivial characters `a -> (-1)^(a dot v_i)`. None is invariant.
This twist matters for the full group `G`; it cannot be discarded there.

## 3. Bad triples measure exactly the failure of crepant quotienting

Let `r` be the number of actual concurrent triples which are Fano lines
in the marking. Each bad triple has eight nodes in `X`. Its node stabilizer
inside `S` is exactly the order-two subgroup attached to that Fano line.
The eight nodes make two orbits of four. Locally the node has an equation

```
u^2+v^2+w^2=0,
```

and the stabilizer changes all three signs. Write the node instead as
`A^2/{+1,-1}` using its quadratic generators. The sign change on all three
generators lifts to `(s,t)->(i s,i t)`, so the quotient singularity is
`A^2/mu_4`, with scalar action: type `1/4(1,1)`.

Its resolution has one exceptional curve `E` of square `-4`. Adjunction
gives `K.E=2`, so in `K_res=f^*K_sing+aE` one has `a=-1/2` and the
canonical-square correction is `a^2 E^2=-1`. The good triples remain
`A1` and contribute no discrepancy. Since the action on `X` has no fixed
divisor, `K_(X/S)^2=16/8=2`. The minimal resolution `Z` of `X/S` therefore
has

```
# A1 points of X/S               = t-r,
# 1/4(1,1) points of X/S         = 2r,
q(Z)=p_g(Z)=0,   K_Z^2=2-2r.                              (2)
```

For the cohomology statement one can use rationality of these quotient
singularities and invariant cohomology of `X`; equivalently, the lifted
action on `Y` fixes the relevant `(-2)`-curves pointwise, with quotient
normal bundle `O(-4)`. This identifies `Y/S` with the smooth resolution
locally. Formula (2) does *not* classify the minimal model or Kodaira
dimension when `r>0`. An ample canonical class of the singular quotient
does not settle those questions through noncanonical singularities. The
incoming theorem cited above supplies the missing geometry: an explicit
genus-one pencil for `r=1`, and four (-1)-curve contractions to an Enriques
surface for `r=3`. These classifications are inherited, not consequences
of formula (2) alone.

## 4. Every arrangement in the stated seven-root family has free choices

There are exactly `7!/168=30` distinct simplex subgroups on a fixed set
of seven coordinate labels. The factor 168 is `|GL(3,2)|`; relabelling
all `v_i` by one linear automorphism leaves the simplex subspace unchanged.
Conversely its seven weight-four words recover the complementary Fano
lines, so no additional identifications occur.

The actual concurrent triples form a partial Steiner triple system: two
different triples cannot contain the same pair of lines, since that pair
has a unique intersection and no four lines are concurrent. We enumerate
**all 5,596** partial Steiner triple systems on seven labelled points,
including combinatorial systems not realizable by complex line arrangements.
Each has at least eight Fano structures disjoint from it. The lemma is
**FINITE-EXACT**, independently replayed in
[the audit receipt](../audits/fano_independent_receipt.json). The resulting
**PROVED deduction using that finite lemma** is: every seven-root surface covered by T1 admits at
least eight distinct free simplex subgroups and degree-eight Campedelli
quotient maps. This counts subgroups/maps, not isomorphism classes of
the resulting surfaces.

The entire finite histogram is compact:

| Actual triples in enlarged universe | Free Fano structures : number of systems |
|---:|---|
| 0 | 30 : 1 |
| 1 | 24 : 35 |
| 2 | 18 : 70; 20 : 315 |
| 3 | 16 : 735; 17 : 840 |
| 4 | 14 : 2100; 15 : 210 |
| 5 | 12 : 1050 |
| 6 | 10 : 210 |
| 7 | 8 : 30 |

The bound eight is sharp in this enlarged combinatorial universe, attained
by a full Fano system. No characteristic-zero realization of that extreme
is asserted. The same conclusion for the actual magic family has the
stronger explicit controls in the receipt: all 36 labelled seven-subsets
have 12 to 20 free subgroups.

## 5. Hourglass census and an explicit free action

Order the hourglass forms as

```
(E, E+x, E-x, E+z, E-z, E+x+z, E-x-z).
```

Their three concurrent triples are `(0,1,2),(0,3,4),(0,5,6)`. Among the
30 Fano structures, the obstruction count is

| Bad triples r | Number of simplex subgroups | Singularities of X/S | K^2 of resolution |
|---:|---:|---|---:|
| 0 | 16 | 3 A1 | 2 |
| 1 | 12 | 2 A1 and 2 of type 1/4(1,1) | 0 |
| 3 | 2 | 6 of type 1/4(1,1) | -4 |

There is no `r=2`: two lines of a Fano pencil force the third pairing of
the remaining points. Here is an explicit free labelling, with integers
read as three binary bits:

```
(v_E,v_(E+x),v_(E-x),v_(E+z),v_(E-z),v_(E+x+z),v_(E-x-z))
       = (1,2,4,7,5,6,3).
```

The three actual triple sums are `7,3,4`, all nonzero. Thus every nontrivial
element of its displayed simplex code acts freely on the actual hourglass
resolution. The receipt supplies all eight codewords and every marking.

For a short independent count: a fixed triple belongs to six of the 30
Fano systems. The three prescribed star triples extend to exactly two
Fano systems (choose one of the two Latin-square completions on the three
opposite pairs). Each selected star triple occurs in four further systems,
which contain neither other star triple. This gives `r=3:2`, `r=1:12`,
and the remaining `r=0:16`.

## 6. The Boolean reflection becomes a genuine automorphism between markings

Take the incidence-adapted labelling

```
(7,1,6,2,5,4,3)
```

in the same form order. Label 7 is full Boolean support, and the three
opposite pairs are singleton/complementary-pair supports. The parameter
automorphism `(E,x,z)->(E,-x,-z)` fixes the central root label and swaps
all three opposite pairs. It lifts to the root surface by permuting root
coordinates. On these labels its permutation is exactly

```
J(7)=7,        J(v)=v+7 for v!=7,
```

the modified Boolean complement from [reflection_bridge/seven.md](../../reflection_bridge_20261025/seven.md).
It therefore *does* have a geometric realization on this seven-root
surface. This does not turn it into an automorphism of the monoid `P^2`
support poset, nor into a Fano collineation: those earlier obstructions
retain their original domains.

Instead, `J` conjugates one simplex subgroup to another. It exchanges
the two Fano completions of the actual three-star, and hence induces an
isomorphism between the two `r=3` quotient surfaces. They have the same
six noncanonical singularities; neither is a free Campedelli quotient.
More generally this reversal fixes none of the 30 simplex subgroups and
pairs them into `8+6+1` pairs in the `r=0,1,3` strata. An abstract reason
for no fixed subgroup is that an involution in `GL(3,2)` has fixed
subspace dimension two and therefore cycle type `1^3 2^2` on its seven
nonzero vectors. Our permutation has type `1^1 2^3`.

The strongest survivor of the old failed incidence bridge is thus a
*correspondence between markings*, not an automorphism of one marking.
The extra coordinate is the chosen Fano structure. The initially open
isomorphism question is now answered by the incoming theorem: the intrinsic
bicanonical map makes the full D12 orbit test complete in the free case.
Its two orbits have sizes4 and12 and give exactly two complex surface classes.
The corresponding converse for the properly elliptic quotients remains open.

## 6a. The same subgroup acts regularly on the eight central K3 conics

The [root-deletion packet](../root_deletion/DELETION_AND_CUBE.md),
independently checked in [the root audit](../audits/ROOT_DELETION_AUDIT.md),
exhibits eight conics on the central six-root K3. Their three sign bits
record the relative signs within the three opposite root pairs. For any
Fano marking put

```
c=v_0,  d_1=v_1+v_2,  d_2=v_3+v_4,  d_3=v_5+v_6.
```

An element `a` of the simplex group translates the conic cube by
`(a dot d_1,a dot d_2,a dot d_3)`. The seven nonzero labels have sum zero,
so `d_1+d_2+d_3=c`. Each `d_j` is nonzero. The j-th concurrent triple
is bad exactly when `d_j=c`, equivalently when the other two differences
are equal. Thus no bad triples means three distinct nonzero differences
with nonzero sum; over `F2` these form a basis. Conversely a basis gives
no bad triple. We obtain the actual geometric equivalence

```
S acts freely on the seven-root hourglass resolution
    iff S acts regularly on the eight central K3 conic components.      (3)
```

For `r=1` the translation rank is two and there are two orbits of four
components. For `r=3` it is one and there are four orbits of two. This
does not assert that the action on the *whole* K3 is free; only its action
on the set of components is regular in the first case.

For example choose root labels `(7,6,5,4,2,3,1)`. Then
`(d_1,d_2,d_3)=(3,6,2)` is a basis and the three triple sums are `4,1,5`.
The seven nonconstant Walsh characters of this component cube and the
seven canonical coordinate characters of `S` are now identified by the
specific map

```
nonempty J subset {1,2,3} -> sum_(j in J) d_j -> unique root label v_i.
```

It preserves the character of the same group action. The three singleton
Walsh labels have zero intersection eigenvalue and are actual divisor
relations, while the other four have negative eigenvalues. Therefore this
does not identify seven canonical forms with seven independent divisor
classes. The empty Walsh character is the positive class and remains a
necessary eighth coordinate for the whole intersection form. All 5,040
markings, including the nonfree ranks, are checked independently in the
root-audit receipt.

## 7. Connection contracts, source recovery, and stopping boundary

| Source -> target | Map and preserved predicate | Loss / sidecar / hostile |
|---|---|---|
| Arithmetic Fano gauge kernel -> root deck subgroup | identical linear code `S`, used as diagonal signs | arithmetic sizes, octonion product and tournament counts lost; keep labels |
| Root cover -> Campedelli cover | quotient by `S`, residual inertia is syndrome label | root signs modulo `S` lost; concurrency/Fano intersection tests freeness |
| Seven supports -> seven root sections | chosen labelling, canonical `S`-characters | support order not transported; determinant twist retained |
| Boolean reflection -> geometric quotient isomorphism | root permutation conjugates `S` to `J S J^-1` | does not preserve one Fano plane; two markings are necessary |
| Seven to six in literature | coordinate/reflection quotients and period constructions | generic smooth hypotheses matter; nodal hourglass not automatically Enriques |

**Primary sources checked 2026-10-25.** Alexeev–Pardini,
[*Explicit compactifications of moduli spaces of Campedelli and Burniat surfaces*,
v3, Definition 3.1 and Lemma 3.4](https://arxiv.org/html/0901.4431v3),
identify the classical seven-line cover and its independent/dependent
triple singularities. Oudompheng,
[*Periods of an arrangement of six lines and Campedelli surfaces*,
v1, sections 1.1–1.3](https://arxiv.org/html/1106.4846v1),
describes the diagonal-character construction and quotient reflections.
These general mechanisms are recovered, not discoveries here. Its
Proposition 1.9 assumes a smooth Campedelli surface for the Enriques
quotient conclusion; the following nodal warning prevents its automatic
application to our hourglass. The new packet's targeted contribution is
the explicit bridge to the inherited gauge kernel, the marking census,
and the reflection correspondence, with their exact domains.

The birational classification for `r>0` and the two free isomorphism classes
are supplied by the incoming audited theorem, as recorded above. The full
classification of the properly elliptic quotient isomorphism classes remains
open. The positive result already supplies an
actual geometric map from the seven-root hourglass to a well-understood
surface family; it gives no positive integral square point and no LRC(14)
or tournament conjecture implication.

## Reproduction

From the repository root:

```
python -B 05-knowledge/results/six_seven_geometry_20261025/fano_codes/check.py 05-knowledge/results/six_seven_geometry_20261025/fano_codes/receipt.json
python -B -O 05-knowledge/results/six_seven_geometry_20261025/fano_codes/check.py 05-knowledge/results/six_seven_geometry_20261025/fano_codes/receipt.optimized.json
```

Both runs are byte-identical. The script imports no producer/repository
implementation and uses explicit exceptions. It checks all 128 sign words,
all 5,040 label bijections and their 168-element fibers, all 30 subgroups,
every allowed root zero-set, every magic seven-subset, and all 5,596
partial Steiner systems. The latter universe is generated independently
by used-pair recursion and by clique recursion in the triple compatibility
graph. It also verifies the incidence-kernel/code-syndrome identities and
the Boolean reflection using the actual rational form coefficients. The
local quotient invariant-ring check is finite supporting evidence; its
geometric type and discrepancy are proved in section 3. The receipt
contains the source SHA-256 and full consequence data.
