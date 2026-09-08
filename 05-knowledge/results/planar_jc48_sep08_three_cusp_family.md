# Marked D4 quotients throughout the ordinary three-cusp (4,6) class

**Status: PROVED / FINITE-EXACT / INDEPENDENTLY AUDITED.**
Exact algebra and named singularity controls accompany the proof. The
actual representative braid supplier is the separate
[three-cusp braid bundle](planar_jc48_sep08_three_cusp_braid.md); the
family result and the representative supplier have both passed independent audit.

## 1. Class and conclusion

Consider an irreducible complex affine curve admitting a birational
polynomial normalization (U(t),V(t)) with deg U=4 and deg V=6. Assume
its affine singularities are exactly three ordinary cusps and ordinary
nodes, with no other singularities. Then its projective infinity branch
has type (2,7) or (2,9); the corresponding numbers of finite nodes are
four and three, respectively.

Every curve in this declared class has four positive curve meridians
generating its affine complement group and satisfying the marked Artin
D4 relations. The three cusp pairs and three distinct node pairs can
be chosen with the marking of the representative braid supplier.
Consequently, if such a curve were the whole irreducible nonproperness
set of a nonautomorphic polynomial Keller map, the
[three-cusp passport](planar_jc48_sep08_three_cusp_passport.md) gives
mapping degree at least eight.

This is a statement about a specified polynomial-normalization class;
no claim that all Keller curves admit these coordinates is made.
The complement is a quotient of Artin D4, with no assertion of equality.
No all-degree Keller exclusion or polynomial source realization follows.
The fixed-sheet involution classifier with exactly three nodes is not
used for the four-node stratum; its extra node contribution would need
its own supplier.

The closest proved mechanisms are the corrected proper-incidence and
marked resolution construction in the
[two-cusp family](planar_jc48_sep08_two_cusp_family.md), §4, and the
new representative's exact local-cluster braid certificate. The named
hostile is the even double-cover parameter (s,a)=(0,0). The corrected
near miss is deleting good curves because two cusp U-values co-project,
or transferring an infinity-nine statement across a degeneration by
an unproved global isotopy. The needed sidecars are ordinary cusp jets,
strict finite local braid certificates, the genus change, and proper
marked transport within each constant infinity stratum.

The five live concepts are the complete coefficient family, finite
cusp jets, infinity characteristic terms, proper good loci, and marked
meridian transport. The fixed finite braids survive the change of
infinity type; the global pair is trivialized only within a constant
infinity stratum. These are distinct operations in the proof.

## 2. Complete two-parameter normal form

The three cusps have distinct normalization preimages. Each preimage
is a common zero of U' and V'. Since U' has degree three, its three
roots are exactly those distinct preimages, each simple. There is no
additional repeated-projection-critical chart to recover in this class.
Normalize two preimages to -1,1 by an affine parameter change, and write
s for the third. Scale and translate the target coordinates so the
leading coefficients of U,V are one and U(0)=V(0)=0. Then s!=+-1 and

    U'=4(t^2-1)(t-s).

The sextic derivative is 6(t^2-1)(t-s)(t^2+a t+q) for some a,q.
Replacing V by V+(3/2)(as-q)U is a genuine invertible linear target
shear and makes q=as. Thus the complete representative family is

    U=t^4-(4s/3)t^3-2t^2+4st,
    V=t^6+(6/5)(a-s)t^5-(3/2)t^4
        +2(-as^2-a+s)t^3+6as^2t.                        (1)

It satisfies

    U'=4A(t),  V'=6A(t)(t^2+a t+as),
    A(t)=(t^2-1)(t-s).                                  (2)

The version with zero constant term in the final quadratic differs
from (1) by a multiple of U; it must not be confused with the literal
representative at (s,a)=(2,-4/3), which is exactly (1).

At any e in {-1,1,s},

    (U''V'''-V''U''')(e)=48 A'(e)^2(2e+a).               (3)

The ordinary-cusp condition is therefore precisely

    O(s,a)=(s^2-1)(a^2-4)(a+2s) != 0.                  (4)

This includes the condition that the three preimages are distinct.
Their individual jets are 192(a-2)(s+1)^2, 192(a+2)(s-1)^2,
and 48(a+2s)(s^2-1)^2. Their nonvanishing guarantees local ordinary
cusp branches, before excluding other preimages of those images.

Birationality holds throughout (4), not merely at a generic parameter.
Indeed [C(t):C(U,V)] divides both four and six, so it is one or two.
In the degree-two case its field involution fixes the unique pole of U;
as an automorphism of P1 it is an affine involution t -> beta-t. Both
U and V would be even polynomials in t-beta/2. At its fixed point
both derivatives vanish and both third derivatives vanish, giving a
zero determinant in (3). The common critical parameters are exactly
the three roots of A, all with nonzero jets, a contradiction.
This proof uses the jets, since a degree-three derivative gcd by itself
would not exclude the quadratic cover.

Finiteness over the image follows from the monic equation for t over
C[U,V]. The projective normalization is P1, its image is a sextic,
and the single infinite parameter maps to [0:1:0]. The explicit
(s,a)=(0,0) even pair is the excluded quadratic-cover control: the
jet at t=0 vanishes, exactly where the argument requires a unit.

## 3. Only two infinity strata

At z=1/t, put X=U/V and Z=1/V. The leading terms are X~z^2 and
Z~z^6, with contact six along the original infinity line. Exact cleared
numerators give

    [z^7](Z-X^3)=h7=4(3a+2s)/5.                         (5)

If h7 is nonzero the infinity branch is (2,7). On the zero line
 a=-2s/3, put k4=3-4s^2/3. After the lower even removal,

    [z^9](Z-X^3-k4 X^4)=h9=16s(s^2-9)/27.              (6)

All lower coefficients in the indicated expressions vanish. A zero
of (6) has s=0,3 or-3. On this line each of those parameters makes
one of the ordinary cusp factors (4) zero. Thus no infinity type
at least eleven occurs in the declared ordinary three-cusp class.

The starting spaces are therefore the irreducible open sets

    T7={(s,a) in A2: O(s,a)(3a+2s)!=0},
    T9={(s,-2s/3): s not in {0,1,-1,3,-3}}.              (7)

The vanishing of k4 at s=+-3/2 does not affect (6) or the infinity
type, and those points are not deleted for that reason. On a good
curve the genus identity is

    10=3+(m_infinity-1)/2+N,

so N=4 on T7 and N=3 on T9.

Projection collisions are also not bad-curve conditions. The cusp
U-values are

    U(1)=-1+8s/3, U(-1)=-1-8s/3, U(s)=2s^2-s^4/3.

Their differences include

    U(s)-U(1)=-(s-1)^3(s+3)/3,
    U(s)-U(-1)=-(s+1)^3(s-3)/3.                          (8)

Thus s=0 or+-3 can co-project cusps on T7 while their target points
remain distinct. These parameters are retained whenever the intrinsic
curve is good. They are already excluded on T9 for the independent
ordinary-cusp reason in (4), not to simplify a projection argument.

## 4. Explicit intrinsic good loci and their properness

For m=7,9, define G_m in T_m by the following construction. First
remove extra preimages of each prescribed cusp target. This condition
is explicit in the original coefficient coordinates. Define

    Rminus=45a^2s^2-162a^2s-27a^2+80as^3-216as^2-216as
           +128s^3-432s^2,
    Rplus =45a^2s^2+162a^2s-27a^2+80as^3+216as^2-216as
           -128s^3-432s^2,
    Rthird=36a^2s^2-216a^2+28as^3-108as+11s^4-126s^2+675.

Require Rminus*Rplus*Rthird!=0. In fact the three resultants

    Res_t((U(t)-U(e))/(t-e)^2, (V(t)-V(e))/(t-e)^2)

are, in order e=-1,1,s,

    -4(s+1)^2 Rminus/675,
    -4(s-1)^2 Rplus/675,
    (s^2-1)^2 Rthird/675.                               (9)

The first divided polynomial is monic quadratic in t. At t=e it is
U''(e)/2, which is nonzero by (4). Thus these incidences have no
unremoved diagonal point and are finite over the parameter space.
Their images are closed, and (9) excludes precisely an extra preimage
of each cusp image. It also excludes coincident cusp images. It does
not exclude two distinct cusp images that merely share a U-value.

Over this open base, resolve the three cusp sections and the infinity
branch in the proper P2 family. This is a fixed, simultaneous marked
resolution on each T_m. For finite cusps the pullback of the initial
centre ideal is (t-e)^2 times a unit ideal; U''(e) is a unit. The
ordinary jet (3) pays the subsequent fixed-order centres. Excluding
extra preimages before making these blowups is essential: otherwise
the prescribed normalization would not automatically lift through the
chosen centre sequence.

At infinity the first three charts are Z=X Z1, Z1=X Z2, Z2=X Z3.
The curve meets the third exceptional line at Z3=1. The strict infinity
line is at Z3=0 and the other marked boundary intersection is at
infinity in that chart. With Y=Z3-1, the remaining curve has orders
(2,m-6) after removal of lower even powers. Those changes replace Y
by Y minus a polynomial in X, fix the exceptional axis X=0 and its
marked intersections, and divide by no even coefficient. The final
odd coefficient h7 or h9 is a unit on its stratum. The remaining fixed
point-blowup sequence resolves the curve together with the infinity
boundary, including at k4=0. This is a resolution on each constant-type
stratum, not an asserted resolution trivialization across their boundary.

The lifted normalization P1 times the base is proper and is now
unramified, with nonzero differential on each branch: its only finite critical parameters were the three resolved
cusps, and infinity is resolved. Its diagonal in the fibre product with
itself over the resolved surface is open by unramifiedness and closed
by separatedness. Its off-diagonal double incidence is therefore
proper over the base. The pairwise-distinct triple incidence is proper
for the same reason, by excluding the three open-and-closed diagonals.
On the double incidence, zero tangent determinant is a closed condition.
Remove the closed parameter images of nontransverse double pairs and
triple points. The remaining open locus is, by definition, G_m.

This construction is an explicit incidence definition of G_m; expansion
of its last elimination polynomials is not needed. All incidence maps
whose images are removed have just been proved proper. It is exactly
the locus where the three cusp targets have no extra preimage and all
other affine singularities are ordinary nodes. At infinity no extra
branch can occur because the polynomial normalization has a unique
infinite preimage and its branch is already resolved. Thus G_m is a
Zariski open subset of T_m, with no unexplained properness assumption
on raw off-diagonal affine pairs.

## 5. Nonemptiness and persistence across the infinity boundary

The exact braid representative (s,a)=(2,-4/3) belongs to G9: its complete
geometry, three nodes and local cluster disks are certified separately.
Every root-tube, endpoint-order, local-cluster and critical-boundary
inequality in that supplier is strict and concerns finitely many fixed
rational paths or polynomial coefficient bounds on fixed compact disks.
The actual resultant coefficients depend polynomially on (s,a).
Therefore all those strict inequalities persist in some Euclidean
parameter neighborhood of the representative, using the same rational
paths, centres and radii. This supplies the same actual marked braid
relations, not just the same endpoint permutations.

For completeness this neighborhood meets G7, rather than merely T7.
The three ordinary cusp jets and no-extra-preimage resultants remain
nonzero near the representative. Each of its three ordinary node pairs
persists uniquely by the implicit function theorem applied to
U(t1)-U(t2)=V(t1)-V(t2)=0: the transverse tangent determinant is nonzero.
Their disjoint normalization neighborhoods remain distinct from all
cusp sections. Every nearby parameter off the h7=0 line has infinity
(2,7). Birationality has already been proved throughout the ordinary
locus, so its rational sextic genus budget is ten. The three cusps,
the three persisting nodes, and infinity now account for nine units
of delta invariant. Exactly one unit remains at finite distance.

No additional unibranch singularity is possible: every finite critical
normalization parameter is one of the three already prescribed cusps.
An additional nontransverse double point or triple point has delta at
least two. The remaining unit therefore gives exactly one further
ordinary node. It cannot merge into a persisting singular neighborhood,
where the established local ordinary model and lack of extra cusp
preimages are stable. Consequently sufficiently close parameters off
the infinity-nine line lie in G7. The extra node goes out toward
infinity in the limit; no global cross-stratum isotopy is asserted.

The six fixed local lassos still surround exactly their three persisting
cusps and three persisting nodes. The discriminant is nonzero on each
fixed local boundary, and its zero count is stable there; the
prescribed singularities already supply that complete local count.
Thus the actual local-pair identification, not just the global words,
persists. The fourth node is not among those six marked local disks.
This proves both nonemptiness of G7 and the existence of a G7 supplier
with the marked D4 quotient. It uses strict certified finite paths plus
the genus argument, not deformation invariance across unequal infinity
resolution graphs.

Independent exact controls make the retained exceptional loci concrete.
At (s,a)=(0,1), the discriminant is

    -65536/244140625 * u^3(u+1)^6(9u+5)^2
      * (625u^3+699u^2+663u+289)^2.

This is a good G7 curve with four nodes and two distinct cusps over
u=-1. At (s,a)=(3,1), it is

    -5308416/244140625 * (u-7)^3(u+9)^6(9u+17)^2
      * (625u^3-10881u^2+49707u+171549)^2,

again a good G7 curve with four nodes and a cusp coprojection.
At (s,a)=(3/2,-1), where k4=0, it is

    -81/256 * (u-3)^3(u+5)^3(16u-45)^3
      * (16u^3+36u^2+27u-765)^2,

a good G9 curve with three nodes. In each control the cusp fibre gcds
are exactly (t-e)^2, all jets are nonzero, and the residual node factor
is squarefree and avoids every cusp projection and every critical value
of U. Its discriminant squares therefore give ordinary nodes by the
analytic-root argument of the braid supplier. The co-projected cusps
account for exponent six; they are not a triple image or a missing node.

## 6. Connectedness and marked proper-SNC transport

Each nonempty G_m is a Zariski open subset of an irreducible affine
plane or line, and is path connected over C. One can join two points
on their complex affine line while avoiding its finitely many bad
points; the line cannot be contained in the bad locus because its
endpoints are good. Thus the good projection-collision and k4=0
controls lie in the same good component as the respective suppliers.

Over G_m the ordinary nodes form a finite étale multisection. Blow up
that entire multisection, without choosing global node labels, in the
resolved proper surface family of §4. Together with the cusp and
infinity blowups this produces a smooth proper family with relative
simple-normal-crossing reduced total divisor formed by the curve,
the infinity line and all exceptional components. Local labels, if
needed, are supplied étale locally; the construction and its complement
are global. All blowup centres lie in the deleted total divisor, so
its complement is the original affine curve complement in each fibre.

Along a compact parameter path, choose local lifting vector fields in
the normal-crossing divisor charts tangent to each divisor stratum.
Patch these fields and use properness to obtain the complete flow.
It trivializes the pair along the path and identifies the affine
complements. This is the marked proper-SNC mechanism already proved
in the predecessor family, here with three cusp sections and constant
node count on each G_m. The complex normal orientation preserves
positive meridians. The cusp neighborhoods and node neighborhoods are
also carried to their counterparts, so their local meridian pairs and
simultaneous access transport are retained. Individual node labels may
permute along a loop, but along a chosen path the three named node
neighborhoods can be followed to three distinct nodes at its endpoint.

Apply this within G9 to the literal certified representative, and within
G7 to a nearby supplier from §5. The four generating positive meridians
and the six actual local pairs are transported. At the destination they
need not be the standard generators of a particular chosen U fibre;
they are four actual positive curve meridians with the required group
and local marking. This is sufficient for the passport consumer.

No hypothetical Keller map or source is transported through parameter
space. At a destination curve, any hypothetical Keller monodromy is
applied to these actual transported loops, and the local retained-page
supplier is applied to its own cusp and node neighborhoods. The
isomorphism of marked target complements pays the group/access data;
actual retention and Euler specialization are supplied by the
hypothetical map at that destination. This distinction prevents a
topological group quotient from being mistaken for a polynomial source.

Finally every curve in the declared class has a representative in one
of the G_m by §2 and the infinity exhaustion. The parameter and target
normalizations are affine isomorphisms, so they preserve the relevant
marked complement predicates. This proves the family statement.

## 7. Exact verification and accepted audit

The [standalone source](../../04-computation/planar_jc48_sep08_three_cusp_family.py)
verifies the complete derivative/shear laws, three jets, no-extra-cusp
resultants, harmless projection coincidences, lower and leading infinity
coefficients, exclusion of higher infinity by the ordinary factors,
genus inventories, the quadratic-cover hostile, and the four named good
curve controls. Its universe is those explicit symbolic identities and
four fixed rational parameter points. No bounded parameter census is
used to infer openness, connectedness or marked transport.

```sh
python3 -B 04-computation/planar_jc48_sep08_three_cusp_family.py
python3 -B -O 04-computation/planar_jc48_sep08_three_cusp_family.py
```

Normal and optimized replays pass **102 always-active gates** and
reproduce the frozen output byte for byte. Frozen pins are:


- `planar_jc48_sep08_three_cusp_family.py`: 5,571 bytes; SHA256 `4f09e896a714213266ec3483bd595f6aae24c5a51f65b201a2b574adb8162fbb`.
- `planar_jc48_sep08_three_cusp_family.out`: 1,459 bytes; SHA256 `2631e28ff306b26f45a306c037956bdf0015cc66ef531202288e7dcb3fbb72c2`.

The [independent audit](planar_jc48_sep08_three_cusp_family_audit.md) accepts the complete analytic proof and all source/replay controls, with a separate23-check original-integral and local-series reconstruction. The representative braid is a proved supplier; its source, output and witnesses are unchanged. Shared navigation is owned by root.
