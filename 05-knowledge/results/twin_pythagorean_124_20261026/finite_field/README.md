# A paired finite-field design carried by the Pythagorean branches

**Status: PROVED elementary bridge; FINITE-EXACT full-group and field-32
controls.** The inherited point-stripping ladder and Paley/F8 dictionary are
credited below. The contribution of this packet is their unpunctured paired
design, its exact Berggren transport, and the explicit endpoint-loss and
larger-field boundaries. No novelty claim is made for the classical designs
or their abstract automorphism groups. No twin-prime, abc, or Hamilton-spectrum
conjecture is resolved.

## 1. Inheritance and the map contract

The closest mechanism is [THM-4553, point-stripping ladder](../../../../01-canon/theorems/THM-4553-the-point-stripping-ladder-psl27-borel-torus-and-the-collatz-frobenius.md):
`PSL2(7) -> C7⋊C3 -> C3` has action degrees `8 -> 7 -> 6`, after the
required choices of poles and markings. The inherited
[conic packet](../../pythagorean_form_conic_graphs_rp_margin_20261018.md)
already corrects the whole three-branch image at 7 to **PGL2(7)**, not PSL2(7).
The branch conventions and integer tree are fixed in the concurrent
[tree result](../tree/RESULT.md).

The incoming [Session69 tournament report](../../session69_prime_ports_power_aliquot_20261026/session69-tournament/REPORT.md),
read at `5e645b2a88`, already proves the trace rule on F8, the full order-21
Paley action, both Fano systems, and the action-module decomposition
`1+3+3`. Its field-32 hostile concerns the trace rule. Our field-32 hostile
below concerns preservation of affine planes; these are distinct predicates.
The [earlier marked 1-2-4 note](../../aliquot_power_124_20261026/tournaments/TOURNAMENTS_124.md)
already separates the dependent trace-zero Frobenius orbit from a normal
basis. We do not replace that distinction by a dimensional `1+2+4` claim.

The canonical hostile is the extra nonsquare determinant branch. The corrected
near miss is confusing one design's full group with the conic image group.
The least-used sidecar is the **ordered pair** of designs, rather than either
single Fano plane. Our live board is: conic points; affine planes over F2;
determinant square class; branch ancestry; pole-selected tournament; and
multiplicative norm circles.

| Connection | Source and target | Preserved predicate | Lost data / sidecar |
|---|---|---|---|
| Reduction and parametrization | integer Euclid pair to conic over F7 | quadratic equation and branch action | integer size, ancestry and primality; retain the integer pair for these |
| Exponential marking | eight conic points to additive F8 | paired design and the projective group action described below | circle multiplication; retain the marking and pole |
| Branch transport | branch word to point plus one bit | which member of the ordered design pair is carried to which | point alone loses this bit |
| Pole restriction | seven remaining points to Paley7 | exact arcs under the square-determinant stabilizer | full PGL action does not preserve the chosen orientation |

These are typed maps, not identifications of all arithmetic, geometric and
tournament predicates on the same-sized set.

## 2. The common slope convention and the eight-point marking

Write an integer primitive triple as

    Phi(m,n)=(m²−n²,2mn,m²+n²),  m>n>0.

The common slope is **r=n/m**. Introduce the explicitly reciprocal auxiliary
coordinate **t=m/n=1/r**, with 0 and infinity exchanged. The root `(m,n)=(2,1)`
has `r=4`, `t=2` over F7. We use t only to make the following marking and
group generators short. The conic parametrization in these two coordinates is

    C(r)=[1−r²:2r:1+r²],
    C(t)=[t²−1:2t:t²+1],       C(t=infinity)=[1:0:1].       (1)

Because −1 is not a square in F7, the conic has no point with z=0. It is the
eight-point affine circle `x²+y²=1` after normalizing z=1. The pole removed
later is `t=infinity`, equivalently `r=0`; it is not the integer tree root.

Let `F8=F2[alpha]/(alpha³+alpha+1)`. Mark the points by

    E(t)=alpha^t for t in F7,       E(infinity)=0.             (2)

In basis `(1,alpha,alpha²)`, the list for `t=0,...,6,infinity` is
`(1,2,4,3,6,7,5,0)`, using bit coordinates. In the common r gauge the same
marking is `E(1/r)`, with the reciprocal conventions just stated. This is a
chosen bijection, not a homomorphism of the circle group to a field group.

## 3. The paired-design theorem

On `X=P1(F7)`, define two collections of four-element subsets:

    D+ = { B : sum_(t in B) E(t) = 0 },
    D− = { B : sum_(t in B) E(t)^3 = 0 }.                     (3)

Sums here are in F8, equivalently coordinate XOR. Each collection is a
Steiner quadruple system `S(3,4,8)`: every triple belongs to exactly one block.

**Theorem F1.** The two designs have 14 blocks each and are disjoint. Every
square-determinant projectivity preserves each one. Every nonsquare-determinant
projectivity exchanges them. Thus their swap bit is exactly the well-defined
quotient character

    PGL2(F7) -> F7*/(F7*)² ≅ C2.                            (4)

*Proof.* Three distinct field elements x,y,z complete to the fourth
`w=x+y+z`, distinct from all three. This gives D+, and the bijective cube map
gives D−. Each system has `binom(8,3)/4=14` blocks. If x,y,z,w are distinct
with sum zero, then

    x³+y³+z³+w³=(x+y)(x+z)(y+z) != 0.                     (5)

So the two block systems are disjoint.

The transformations `T:t->t+1` and `S:t->−1/t` generate PSL2(F7). Under (2),
T becomes `x->alpha*x`, and S becomes `x->L(x)+1`, where the F2-linear map L
fixes 1 and exchanges alpha with alpha². These eight values are checked
directly from the power table. Affine maps preserve zero-sum four-subsets,
so PSL preserves D+. This is an affine representation: **L is not a field
automorphism**.

The projectivity `J:t->3t` has nonsquare determinant and becomes `x->x³`.
Thus D− is the inverse image of D+ under J. Since `J²:t->2t` becomes
Frobenius `x->x²`, J exchanges D+ and D−. Normality of PSL in PGL now makes
PSL preserve D− as well. Finally PGL has the two cosets PSL and J·PSL,
proving exactly the two cases of (4). Scalar representatives multiply a
determinant by a square, so (4) is independent of representative. QED.

**Full-group boundary, FINITE-EXACT.** Enumeration of all 8! permutations,
independent of the projective matrix enumeration, gives

    Aut(D+)                         order 1344,
    Aut(D+,D−)                      = PSL2(F7), order 168,
    Aut({D+,D−})                    = PGL2(F7), order 336.    (6)

The checker compares the actual permutation sets, not only their orders.
The first group is the affine group AGL3(F2): every affine map preserves the
planes; conversely, after translating f(0) to 0, preservation of
`{0,x,y,x+y}` forces `f(x+y)=f(x)+f(y)`. Thus the single design has extra
symmetries. The paired design is exactly what recovers the projective action
in this full finite oracle. The remaining 42 four-subsets are outside the two
designs; we do not need an assertion about their orbit to prove F1.

## 4. Exact Berggren transport and the first lost bit

Use the inherited branch names

    U(m,n)=(2m−n,m),   A(m,n)=(2m+n,m),   D(m,n)=(m+2n,n).

Their common-r actions are respectively

    r -> 1/(2−r),     r -> 1/(2+r),     r -> r/(1+2r).       (7)

In the auxiliary t coordinate they are `2−1/t`, `2+1/t`, `t+2`, with matrices

    U=[2 −1;1 0],    A=[2 1;1 0],    D=[1 2;0 1].          (8)

The determinants are +1,−1,+1. Since −1 is nonsquare modulo 7, U and D
preserve D+ and D− separately; A exchanges them. The two outer branches
generate PSL: `D^−1 U=S`, while D is `T²` and hence generates T. Adding A
gives PGL. The exact field-coordinate transitions, valid also at zero, are

    D: x -> alpha²*x,
    U: x -> alpha²*(L(x)+1),
    A: x -> alpha²*(L(x)+1)^6.                              (9)

Here exponent six realizes inversion on F8* and sends 0 to 0. Consequently
the transported design bit of a branch word is the parity of the number of
A letters. This is a word cocycle, not a function of the final conic point.

**Minimal endpoint hostile.** Acting successively from `(2,1)`,

    UU: (m,n)=(4,3),  Phi=(7,24,25),   bit 0;
     A: (m,n)=(5,2),  Phi=(21,20,29),  bit 1.               (10)

Both have `r=t=6` over F7 and the same projective conic point: their triple
residues `(0,3,4)` and `(0,6,1)` differ by scalar 2. The root and its three
one-step children have four distinct endpoints, so depth two is the first
possible such collision. All 16 states `(point,bit)` are reachable by the
three branches, as the exact finite-state replay verifies. The extra bit
repairs this specific loss, without recovering integer height, branch word,
or whether an odd leg is prime. The actual twin-prime graph map still needs
the integer factorization retained in the [tree lane](../tree/RESULT.md).

## 5. What survives after choosing a pole

Blocks through `t=infinity` give two disjoint seven-line Fano systems on F7.
For D+ and D− the base translates can be taken as `{0,1,3}` and `{0,1,5}`.
Their union is exactly the 14 directed triangles of the chosen Paley
tournament `i->j iff j−i in {1,2,4}`. These inherited facts are recovered
directly by the checker. With the pole fixed, the square-determinant
projectivities are `t->a t+b` with a square, an order-21 group; they are
exactly the full Paley automorphism group. Nonsquare affine maps reverse its
orientation and exchange the two Fano systems.

This explains the role of the marking: the full eight-point action has no
invariant tournament because PSL acts doubly transitively. Puncturing and
selecting the square orbital makes an intrinsic asymmetric binary relation
available. It does not turn the whole conic into a tournament. In particular,
the circle determinant is not this relation. In normalized t coordinates,

    det(C(t),C(u)) =
      2(t−u)(1+tu)/((1+t²)(1+u²)).                       (11)

The distinct points t=1 and t=6 are antipodal, so (11) is zero; Paley has a
strict arc between them. Three antipodal pairs remain after deleting the pole.
Those are honest ties, not arcs to orient arbitrarily.

The replay gives `H(Paley7)=189`, `|Aut|=21`, and `tc=9`. No claim about
global forbidden H values follows. Raw-H exclusions 7 and 21 remain the
separately cited [THM-1370, H-spectrum omits 7 and 21](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md);
normalized tc has its own hypotheses and current finite/universal boundary.
Neither the number 7 of punctured points nor the group order 21 supplies a
map preserving Hamilton counts for arbitrary tournaments.

## 6. Two form obstructions and the first tested field boundary

**Circle group versus field groups.** The eight-point F7 circle is the norm-one
group in `F49=F7[i]`, `i²=−1`. It is cyclic of order eight: `(2,2)` has
successive orders demonstrated by its square `(0,1)`, fourth power `(−1,0)`,
and eighth power `(1,0)`. The additive group of F8 is `C2³`, while F8* has
order seven. There is no group isomorphism from the circle to either, and
every homomorphism from its order-eight group into F8* is trivial. Deleting
one circle point does not repair multiplication: the remaining set is not
a subgroup. The affine-plane marking preserves incidence as proved above,
not the norm-circle multiplication.

**Characteristic-two form hostile.** Over F8, the naive equation
`x²+y²=z²` is `(x+y+z)²=0`, a double line with all partial derivatives zero.
Its projective point set has nine points, exactly the same count as the
smooth conic `xz+y²=0`; the latter has gradient `(z,0,x)`, never zero on a
projective solution. The common count does not identify their schemes or
quadratic forms. The direct enumeration checks both point sets and the
smoothness condition.

**Finite continuation test.** Replace F7/F8 in (2) by a Mersenne prime p and
the field of order p+1, and ask whether `S:t->−1/t` transports to an affine
map of that additive field. The complete tested universes are p=3,7,31,
all p−1 primitive-generator markings in the fixed displayed fields:

| p | field polynomial | primitive markings | affine transports of S |
|---:|---|---:|---:|
|3|x²+x+1|2|2|
|7|x³+x+1|6|6|
|31|x⁵+x²+1|30|0|

For p=31, use bit coordinates and alpha=x. The source affine plane
`{0,1,2,3}` is sent to `{1,0,18,14}`, whose XOR is 29 instead of zero.
This is a concrete failure of the preserved predicate, not merely a group
order mismatch. Every other primitive generator has an explicit witness in
the receipt. An affine-plane-preserving bijection must be affine by the
`{0,x,y,x+y}` argument in §3. Thus no tested marking at 31 repairs this map.
This statement is **FINITE-EXACT for 3,7,31 only**, not a classification of
all Mersenne primes or all possible, non-exponential markings.

## 7. Reproduction and stopping boundary

The checker is standard-library-only, uses explicit failure gates rather than
Python assertions, and prints its full JSON receipt without required arguments:

```sh
python3 -I -B check_finite_field.py > replay.json
python3 -I -B -O check_finite_field.py > replay.optimized.json
cmp replay.json receipt.json
cmp replay.json replay.optimized.json
```

Its universes are all 8! point permutations; all nonsingular 2-by-2 matrices
over F7 modulo their literal permutation action; all 70 four-subsets and
56 three-subsets; all conic points; all 7! Paley vertex permutations and
Hamilton words; all projective points over F8; and all primitive-generator
markings at the three displayed primes. `SOURCE_PINS.json` records the exact
inherited inputs. The paired design and word bit are sufficient finite-state
data for this action. Extending the incidence-preserving marking beyond the
tested fields, recovering the integer primality predicate from finite colour,
or transporting arbitrary Hamilton counts would require additional structure;
none is inferred here.
