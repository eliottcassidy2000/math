---
id: THM-2779
title: "Bockstein--symplectic decoder/frame torsor and Heisenberg root-degree gate"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For every prime p, the finite Heisenberg group H_p has center and
  commutator F_p, normalized symplectic frames encode its central generator,
  and its minimal faithful permutation degree is exactly p^2.  The modular
  decoder ambiguity is the affine socle torsor K_0+F_p N_p, formally twin to
  a fixed-vector frame shear but not canonically or physically identified
  with it.  At p=13 every literal 13-root action kills the commutator center;
  169 states are sharp.  This is a root-degree/typing obstruction, not a
  common-ancestry lift, row exclusion, or LRC(14).
source: heisenberg-root-gate/root-2026-07-28
depends_on:
  - THM-2771-joint-c7-c13-right-wing-mixed-spectrum-and-commuting-square-no-go
  - THM-2772-carrier-allocation-pullback-k4-segre-and-mixed-face-obstruction
related:
  - THM-2356-finite-field-chirp-gram-tomography-and-bockstein-pairing
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
  - THM-2697-filtered-affine-handoff-germ-category-and-base-signature-holotopy-boundary
  - THM-2782-semantic-arm-right-wing-local-unit-and-endpoint-deck-boundary
  - THM-2785-a12-minuscule-carry-fourier-character-and-common-atom-boundary
script: 04-computation/lrc14_heisenberg_decoder_frame_root_degree_thm2779.py
output: 05-knowledge/results/lrc14_heisenberg_decoder_frame_root_degree_thm2779.out
script_sha256: ef6e9f9bcb4f11152d291342a11ae215245d1d19b96c49940a01ba9ea850cbd9
output_sha256: 1feb463864015035ab8d7fcfcddf9cfe8b0ec0a3ed36481f2f66d6a9149182e6
hash_basis: LF-normalized bytes
---

# Merge source -- all-prime decoder/frame torsors and Heisenberg root-degree gate

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

**Scope warning.**  The group-theoretic theorem below is unconditional.  Its
THM-2771 specialization is coefficient-side, and its THM-2772 specialization
is endpoint-frame-side.  Their common central scalar does **not** identify
their physical carriers.  No common-ancestry lift, endpoint-dipole/root-deck
intertwiner, semantic arm, row exclusion, or LRC(14) conclusion is proved.

## 1. Audited inheritance and exact dependencies

The candidate uses exactly two current proved dependencies.

1. **THM-2771,
   `THM-2771-joint-c7-c13-right-wing-mixed-spectrum-and-commuting-square-no-go`.**
   Status: `PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED`.
   It supplies the intrinsic mod-13 Bockstein collapse `S_beta`, the printed
   decoder `K_beta`, the decoded chart column
   `C=(0,2,0,10,0,0,0)`, and the coefficient identity
   `C*N_7=-N_7`.  It explicitly leaves target convolution, common ancestry,
   and the endpoint-dipole/root-deck map open.
2. **THM-2772,
   `THM-2772-carrier-allocation-pullback-k4-segre-and-mixed-face-obstruction`.**
   Status: `PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED`.
   It supplies the endpoint plane `V=F_13^2`, the determinant mixed face
   `det(s,t)`, the count `2,184` of normalized ordered frames, and the formal
   correction `c_j=-a det(s_j,t_j)`.  It explicitly leaves the common
   endpoint atom and the physical realization of the translation square open.

THM-2697 is related vocabulary for a filtered Bockstein/carry handoff, but no
claim from it is used.  THM-2356 and THM-2633 occur only in the connection
ledger in Section 9 and are not proof dependencies.

The closest proved mechanism is THM-2772's determinant mixed face.  The
canonical hostile is THM-2771's commuting operation square: a nonzero mixed
coboundary is not a commutator class.  The corrected near miss is the idea
that the numeral `-a` identifies the two carriers.  It identifies only a
one-dimensional output after choices.  The least-used sidecar is the socle
of the modular target group algebra, which classifies all decoders rather
than just the printed one.

## 2. Uniform finite-Heisenberg theorem

Let `p` be any prime and put

```text
H_p=F_p^3
```

with multiplication

```text
(A,B,C)(A',B',C')
  =(A+A', B+B', C+C'+A B').                         (1)
```

Write `z=(0,0,1)`, let

```text
omega((A,B),(A',B'))=A B'-A' B,
```

and let `pi:H_p -> F_p^2` forget the central coordinate.

### Theorem A: commutator, center, and frames

For every prime, including `p=2`,

```text
[(A,B,C),(A',B',C')]=z^omega((A,B),(A',B')),         (2)

Z(H_p)=[H_p,H_p]=<z> isomorphic F_p,                  (3)

H_p/Z(H_p) isomorphic F_p^2.                          (4)
```

The set

```text
Fr_1(V)={(s,t) in V^2:det(s,t)=1}                     (5)
```

is a principal homogeneous `SL_2(F_p)`-set and has

```text
|Fr_1(V)|=p(p^2-1).                                   (6)
```

More generally there are `p(p^2-1)` ordered frames of every prescribed
nonzero determinant.  Projection `(s,t) -> s` maps `Fr_1(V)` onto
`V-{0}`; for fixed `s`, every fibre is the affine `F_p`-torsor

```text
t_0+F_p s.                                            (7)
```

For arbitrary lifts `x,y in H_p`,

```text
[x,y]=z^det(pi(x),pi(y)).                              (8)
```

Thus a normalized frame is exactly a projected pair whose commutator is the
chosen central generator.  Central changes of the two lifts do not affect
the commutator.  For a marker `a in F_p`,

```text
[x^(-a),y]=z^(-a det(pi(x),pi(y))).                    (9)
```

Equation (9) is the exact Heisenberg encoding of THM-2772's formal
`-a det(s,t)` correction.

### Theorem B: the modular decoder is a socle torsor

Let

```text
R_p=F_p[C_p]=F_p[u]/(u^p-1)
   =F_p[epsilon]/(epsilon^p),       epsilon=u-1.       (10)
```

Suppose

```text
S=epsilon V,                       V(0)!=0.            (11)
```

Then the decoder equation

```text
S K=epsilon                                             (12)
```

has exactly `p` solutions, and they form the affine torsor

```text
Dec(S)=V^(-1)+F_p epsilon^(p-1)
      =K_0+F_p N_p,                                    (13)

N_p=1+u+...+u^(p-1)=epsilon^(p-1).                     (14)
```

Equivalently, multiplication by `S` has rank `p-1` and kernel equal to the
one-dimensional socle `F_p N_p`.  Choosing any one target coefficient and
requiring it to be zero selects a unique decoder, but that section depends
on the chosen target label.

The fibre in (7) and the decoder set in (13) are both affine `F_p`-torsors.
After choosing a nonzero first frame vector and a basepoint in each fibre,
there are equivariant bijections

```text
K_0+lambda N_p  <-->  t_0+lambda s.                    (15)
```

There is no canonical bijection in (15): `N_p` is a target-group-algebra
socle class, while `s` is an endpoint translation vector.  The valid shared
statement is that both gauge actions preserve their central observable:
`S K` on the decoder side and `det(s,t)` on the frame side.

### Theorem C: sharp faithful permutation degree

The minimal faithful permutation degree of `H_p` is

```text
mu(H_p)=p^2                                             (16)
```

for every prime, including `p=2`.

Every action of `H_p` on fewer than `p^2` points kills the full commutator
center `Z(H_p)`.  Every faithful action of degree exactly `p^2` is transitive
and is a coset action on a noncentral order-`p` subgroup.  The bound is sharp:
on `F_p^2` define

```text
(A,B,C).(b,c)=(b+B, c+C+A b).                          (17)
```

This is a faithful transitive action of degree `p^2`.  It is the coset action
for the core-free stabilizer

```text
{(A,0,0):A in F_p}.                                    (18)
```

The center acts in (17) as `p` disjoint `p`-cycles, so the central height is
literally visible.

For odd `p`, the minimal transitive `H_p`-sets have `p+1` isomorphism types,
one for each projective line in `F_p^2`: central shifts of a fixed
noncentral order-`p` subgroup are conjugate, while inner conjugacy preserves
its projected line.  This refinement is not needed for (16).

## 3. Proof

### 3.1 Group law and commutator

The extra coordinate in (1) is the cocycle `f((A,B),(A',B'))=A B'`.
Associativity is the identity

```text
A B'+(A+A')B''=A'B''+A(B'+B'').                       (19)
```

The identity is `(0,0,0)` and

```text
(A,B,C)^(-1)=(-A,-B,-C+AB).                            (20)
```

Substitution into the commutator gives (2).  Nondegeneracy of the
determinant form says an element commuting with both `(1,0,0)` and
`(0,1,0)` has `A=B=0`, proving the center in (3).  Every central value occurs
as

```text
[(1,0,0),(0,d,0)]=z^d,
```

so the derived subgroup is the whole center.

For every nonzero `s`, the equation `det(s,t)=d` is a nonempty affine line
parallel to `s`, with `p` points.  This gives (6), (7), and the count for
every nonzero `d`.  Taking `s,t` as the two columns of a matrix identifies
`Fr_1(V)` simply transitively with `SL_2(F_p)`.  Equations (8) and (9) follow
from (2) and centrality of commutators.

### 3.2 Decoder classification

Since `V` is a unit, multiplication by `S=epsilon V` has the same kernel as
multiplication by `epsilon`.  In
`F_p[epsilon]/(epsilon^p)`,

```text
ann(epsilon)=F_p epsilon^(p-1).                         (21)
```

One solution of (12) is `V^(-1)`, and every other solution differs from it
by (21), proving (13).  The binomial congruence
`binom(p-1,j)=(-1)^j mod p` gives

```text
(u-1)^(p-1)=1+u+...+u^(p-1),
```

including at `p=2`, proving (14).

### 3.3 Lower bound and sharp action

Every orbit of a `p`-group has prime-power size.  If the total action degree
is less than `p^2`, each orbit has size `1` or `p`.  A one-point orbit kills
the whole group.  On a `p`-point orbit, the image is a `p`-subgroup of
`S_p`; because

```text
v_p(p!)=1,
```

that image has order at most `p` and is abelian.  Hence every such orbit
kills `[H_p,H_p]=Z(H_p)`.  The center therefore lies in the kernel of the
whole action, proving the lower bound.

If a faithful action has exactly `p^2` points, it cannot be a union of
orbits of sizes `1,p`, for the same reason.  It must be one orbit of size
`p^2`, so its stabilizer has order `p`.  A central stabilizer has nontrivial
core.  Conversely a noncentral order-`p` stabilizer is core-free: a nonzero
core would equal the stabilizer and make it normal; then its commutators
would lie in its trivial intersection with the center, forcing the
stabilizer itself to be central.

Finally, composing two transformations in (17) gives

```text
(b,c)
 -> (b+B'+B, c+C'+C+(A'+A)b+A B'),
```

which is exactly the product law (1).  The orbit of `(0,0)` is all
`F_p^2`.  If an element fixes every `(b,c)`, then `B=0`; evaluating at
`b=0` gives `C=0`, and evaluating at `b=1` gives `A=0`.  Thus (17) is
transitive and faithful.

## 4. Characteristic-two boundary

The proof above is uniform, but the internal structure at `p=2` must not be
silently treated as the odd-prime case.  One has

```text
H_2 isomorphic D_8,
```

with one identity, five elements of order two, and two elements of order
four.  Explicitly,

```text
(A,B,C)^2=z^(AB).                                       (22)
```

Thus the quotient `H_2/Z(H_2)=F_2^2` carries the quadratic refinement
`q(A,B)=AB`.  The full `SL_2(F_2)` frame torsor exists in the endpoint plane,
and every normalized frame still has commutator `z`, but not every frame
change lifts to a center-fixing automorphism of `D_8`: such an automorphism
must preserve `q`.  Consequently the normalized-frame torsor must not be
called an automorphism torsor uniformly in `p`.

The minimal faithful permutation degree is nevertheless

```text
mu(H_2)=4.
```

Formula (17) is the ordinary faithful action of `D_8` on a square.  The
noncentral order-two stabilizers split into two conjugacy classes, because
the anisotropic line `(1,1)` has only order-four lifts.  This is the exact
`p=2` exception to the odd-prime `p+1` minimal-action classification, not an
exception to the degree theorem.

## 5. Exact THM-2771 specialization

In the notation of THM-2771, the intrinsic Bockstein target collapse is

```text
S_beta=(0,0,8,0,0,0,0,0,9,9,9,9,8)
       in F_13[C_13].                                   (23)
```

In epsilon coordinates,

```text
S_beta=epsilon V_beta,                 V_beta(0)=12.     (24)
```

The printed normalized decoder is

```text
K_beta=(3,5,5,5,7,12,2,8,2,9,8,11,0).                 (25)
```

Theorem B sharpens the single computation to the complete classification

```text
Dec(S_beta)={K_beta+lambda N_13:lambda in F_13}.         (26)
```

The final coefficient of (25) is zero, so it is exactly the section selected
by the last-target-coefficient normalization.  In epsilon coordinates,

```text
V_beta K_beta=1+3 epsilon^12,                            (27)

V_beta^(-1)=K_beta+3N_13.                               (28)
```

Thus the `3 epsilon^12` in THM-2771 is not a defect: it is the chosen socle
gauge.  The full inverse exists and lies three steps away in the decoder
torsor.

The gauge is invisible only after the final uniform chart convolution.
If `C(lambda)` is the target-zero chart column obtained from
`K_beta+lambda N_13`, exact row augmentations give

```text
C(lambda)
 =(0,2,0,10,0,0,0)
  +lambda(0,3,3,7,0,0,0).                              (29)
```

All thirteen local chart columns are distinct, but

```text
sum_e C(lambda)_e=-1,
C(lambda)N_7=-N_7                                      (30)
```

for every `lambda`.  Therefore THM-2771's pointwise constant
coefficient correction is decoder-gauge invariant, while its local chart
representative is not.  This identifies exactly what uniform convolution
forgets.

## 6. Exact THM-2772 specialization

For `V=F_13^2`,

```text
|Fr_1(V)|=13(13^2-1)=2,184.                             (31)
```

For each fixed nonzero first endpoint move `s`, its thirteen normalized
completions are

```text
t_0+lambda s,                  lambda in F_13.           (32)
```

All have the same commutator `z`, and after multiplication by the marker
all produce

```text
z^(-a).                                                (33)
```

This recovers the exact scalar `c_j=-a` in THM-2772.  It does not identify
the `2,184` frames with THM-2772's `2,184` nonzero target sectors: equal
cardinality is not a map.  Nor does (32) identify endpoint frame shear with
the decoder gauge (26) without a new physical intertwiner.

## 7. The exact relation—and nonrelation—between THM-2771 and THM-2772

After choosing `z`, both constructions admit a scalar projection to the
same central copy of `F_13`:

```text
THM-2771 decoder gauge
  -> gauge-invariant coefficient c_j=-a
  -> z^(-a),

THM-2772 normalized endpoint frame
  -> -a det(s_j,t_j)=-a
  -> [x_(s_j)^(-a),x_(t_j)]=z^(-a).                    (34)
```

Across seven charts, the formal product is

```text
product_j z^(c_j)=z^(-7a),                              (35)
```

which cancels the formal obstruction `z^(7a)`.  Equation (35) is a clean
central-extension restatement of the finite Cech invoice.

The map in (34) preserves only the final scalar correction.  It destroys, on
the THM-2771 side, the target-convolution history and local chart
representative; on the THM-2772 side, it destroys endpoint origin, the
actual translation square, and the frame itself.  Missing are:

```text
source: one common physical ancestry carrying both structures;
map: an affine target-dipole/root-deck intertwiner compatible with endpoints;
preserved predicate: common current, semantic word, clock, endpoint order;
sidecar: root label plus nontrivial central/Bockstein height;
decisive test: exhibit z acting nontrivially on the same physical packet. (36)
```

The Heisenberg group packages the desired compatibility law; it does not
construct the compatibility.

## 8. The thirteen-root gate

Let

```text
rho:H_13 -> S_13                                      (37)
```

be any literal permutation action on thirteen root labels.  Since a
`13`-subgroup of `S_13` has order at most `13`, the image of `rho` is either
trivial or cyclic of order `13`.  Hence

```text
Z(H_13)=[H_13,H_13] <= ker(rho).                        (38)
```

In the nontrivial case, `rho` factors through one nonzero linear functional

```text
H_13/Z(H_13)=F_13^2 -> F_13,                            (39)
```

so its kernel has order `13^2=169`: it kills the center and one horizontal
line.  Thus a thirteen-root deck can see at most one abelian translation
coordinate.  It necessarily erases the determinant/commutator mixed face
that (34) is meant to retain.

The obstruction is sharp, not merely negative.  Formula (17) gives a
faithful transitive action on

```text
F_13^2,                         |F_13^2|=169.            (40)
```

Therefore `169`, not `13`, is the minimum size of a literal permutation
sidecar carrying the whole Heisenberg compatibility.  A 13-dimensional
phase representation is a different object.  If `zeta` is a primitive
13th root, `T_B f(x)=f(x+B)`, and
`M_(-A)f(x)=zeta^(-Ax)f(x)`, then

```text
rho(A,B,C)=zeta^C T_B M_(-A)                            (40a)
```

is faithful and obeys (1), because commuting `M_(-A)` past `T_(B')`
contributes `zeta^(AB')`.  Here the center acts by scalars rather than
permuting thirteen basis labels, and projectivization kills that scalar
center.  Thus (40a) does not contradict (38)--(40).

## 9. Connection ledger and cheapest tests

### LRC(14)

THM-2356 already retains target and first-Bockstein spaces of order `169`.
That numerical scale is compatible with (40), but the theorem's chosen
`F_169` identifications are not a Heisenberg action.  A meaningful transfer
would have to specify which coordinate is the horizontal root translation,
which is central height, and show the same physical packet transforms by
(17).  The cheapest decisive test is not another cardinality comparison:
construct the two root moves and check that their operation commutator acts
as the nontrivial `c -> c+1` center on all `169` states.  Failure means the
current `169`-object is an abelian shadow, not the needed sidecar.

### Planar Jacobian/quartic frontier

At `p=2`, the same extension is

```text
1 -> C_2 -> D_8 -> V_4 -> 1.                            (41)
```

THM-2633 excludes `D_4=D_8` as the **whole** transitive quartic Keller
monodromy, leaving `A_4,S_4`.  Thus the characteristic-two Heisenberg group
is presently a hostile for the live quartic monodromy, not a replacement for
the `V_4` Kummer torsor.  A separate `D_8` embedding problem above a
resolvent base would be differently typed and would need its own base map,
ramification, and Jelonek compatibility.  No JC exclusion follows here.

### Tournament lens

The intrinsic pairwise observable on endpoint moves is `det(s,t)`, with
dependent pairs tied and independent pairs carrying twelve possible
nonzero weights.  It is a weighted symplectic relation, not a tournament.
Reducing determinants to a binary quadratic-residue orientation would lose
the exact central exponent needed in (9).  The native carrier is therefore
the Heisenberg central extension (or the weighted symplectic graph), with
the full determinant as sidecar; forcing a tournament is information-losing.

## 10. Exact evidence

Scratch companion:

```text
04-computation/lrc14_heisenberg_decoder_frame_root_degree_thm2779.py
05-knowledge/results/lrc14_heisenberg_decoder_frame_root_degree_thm2779.out
```

LF byte hashes:

```text
script  53a0c0520588fbd79cab26ab1fd15b6a95aebe007cb12a5c30258244727ed579
output  cd0d8febf8883c2a9940a4c843df54ed9e9f4c94b94987fd6f64a213352943a2
```

Reproduce with:

```bash
python3 04-computation/lrc14_heisenberg_decoder_frame_root_degree_thm2779.py
python3 -O 04-computation/lrc14_heisenberg_decoder_frame_root_degree_thm2779.py
```

The companion uses explicit exception gates, not Python `assert`.  It checks:

1. the group law, inverse, cocycle identity, center, derived subgroup, and
   determinant commutator for `p=2,3,5,7,13`;
2. the complete faithful `p^2` action, core-free stabilizer, center cycle
   structure, and full frame/shear census for those primes;
3. the `p=2` order census and quadratic-refinement boundary;
4. THM-2771's printed primitive Bockstein rows, decoder, decoded rows,
   rank-12 multiplication map, all thirteen decoder gauges, epsilon-socle
   identity, full inverse repair, and local-versus-uniform chart boundary;
5. the standard thirteen-label quotient hostile and the sharp 169-point
   faithful control.

The lower-bound proof, classification of exact-degree actions, and the
noncanonical comparison in (34) are mathematical arguments rather than
enumeration.

The companion has `46` explicit exception gates and no Python `assert` node.
Normal, optimized, and stored transcripts agree byte-for-byte.

```text
script_sha256 = ef6e9f9bcb4f11152d291342a11ae215245d1d19b96c49940a01ba9ea850cbd9
output_sha256 = 1feb463864015035ab8d7fcfcddf9cfe8b0ec0a3ed36481f2f66d6a9149182e6
hash_basis    = LF-normalized bytes
```

## 11. Honest frontier

The result closes the abstract algebra asked of THM-2779:

```text
decoder ambiguity = a one-dimensional socle torsor;
fixed-vector frame ambiguity = a one-dimensional shear torsor;
determinant correction = a Heisenberg central commutator;
literal faithful root permutation degree = 169, sharply.                (42)
```

It also proves the exact stopping law:

```text
thirteen root labels kill the center,
so they cannot literally carry the determinant commutator.              (43)
```

What remains is not more finite group theory.  It is the physical
intertwiner in (36), or a proof that the existing LRC carrier cannot support
one.  Until that is constructed, (34) is a precise algebraic twin and
(43) is a root-degree obstruction, not an LRC(14) proof step.
