# The next infinity branch: a genuine marked A6 survivor

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The boundary A6 representation remains valid. Its required extension to the
actual affine complement is now **EXCLUDED** by the independently audited
[two-meridian theorem](planar_jc48_sep07_next_braid.md). The geometry and frozen
computation below are unchanged; this is not a Keller map or counterexample
to JC(2).

## 1. Procedural move and inherited stopping boundaries

The [five-node curve](planar_jc48_sep06_global_curve.md) was excluded
by Nori's strict inequality. Its [four-node successor](planar_jc48_sep06_boundary_plumbing.md)
was at equality, but its two tetrahedral infinity pieces could not move
six labels when the curve meridian was a single three-cycle. Preserve
the finite (2,5) cusp and raise the infinity cusp once more. This changes
the actual marked boundary, rather than merely retrying the old group.

The operation has a small symbolic supplier. In the family

```text
U=t^4+a t^3+b t^2,
V=c t^6+(3ca/2)t^5+d a t^3+d b t^2,
```

the t^3 coefficient remaining after subtracting the formal expression
c U^(3/2)+h U, with h=-3a²c/8-3bc/2, is
a(7a²c+12bc+16d)/16. Cancel it by choosing a=b=1,c=16,d=-19.
The formal expansion is a parameter-selection device only. Exact rational
charts below verify the resulting geometry independently of that device.

The resulting polynomial curve is

```text
U=t^4+t^3+t^2,
V=16t^6+24t^5-19t^3-19t^2.                      (1)
```

It has one affine (2,5) cusp, three ordinary nodes and infinity (2,11).
Nori's margin is now -2. More significantly, its **actual** infinity group
has a surjective A6 representation with the required single-three-cycle
meridian. Thus neither of the two established infinity obstructions
excludes (1). Extension of that representation over the full affine
complement is still a separate necessary condition.

## 2. Complete geometry

The derivative gcd is t and gcd(U,V)=t². At zero the exact local shear is

```text
V+19U-19U²=-14t^5-41t^6-38t^7-19t^8.
```

Together with U of order two, this gives the unique finite cusp of type
(2,5), and its image has no other normalization preimage.

For an unordered off-diagonal pair write p=s+t and q=st. The complete
divided-difference equations reduce to

```text
(2p+1)q=p(p²+p+1),
M=-p² h(p),           h=4p³+10p²+15p+14.
```

The value of p(p²+p+1) at p=-1/2 is -3/8, so the denominator is lawful.
The discriminant of h is -20972 and

```text
Res(h,p(2p²+3p+4)(2p+1))=220864.
```

The pair discriminant is -p(2p²+3p+4)/(2p+1). Thus the three roots of h
give six distinct ordered pairs. The remaining p=0 pair is just the cusp
diagonal. The ideal of the two divided differences and the tangent
determinant has Groebner basis (s+t³+t²+t,t⁴); every off-diagonal pair
is transverse. To exclude a shared three-point fibre, divide V-b by U-a:

```text
R=-3t³+(16a+5)t²+8at-24a-b,
9 rem(U-a,R)=(256a²+232a+49)t²
             +(128a²-8a-3b)t-384a²-16ab-201a-8b.
```

The three coefficient polynomials generate the unit ideal over Q[a,b].
Since R always has degree three, any three-point fibre would force this
impossible divisibility. These arguments also show generic injectivity,
so the proper polynomial normalization is birational and projective
degree six. All affine singularities have been exhausted.

At infinity set z=1/t, X=U/V and Z=1/V, and define

```text
rho=Z/X³-256,
tau=rho/X+15360,
eta=tau/X-753664.
```

The first odd remainder is exactly of order eleven:

```text
Z-256X³+15360X⁴-753664X⁵ = -(17/1024)z^11+O(z^12).
```

The successive charts (X,Z), (X,Z/X), (X,Z/X²), (X,rho), (X,tau),
(X,eta), (X/eta,eta), (X/eta²-1/4848615424,eta) have orders
(2,6),(2,4),(2,2),(2,2),(2,2),(2,1),(1,1),(1,1).
Here eta=-17408z+O(z²), so the final branch is transverse and away
from divisor intersections. The original infinity line exits at the
third blowup, rho=-256. The finite cusp multiplicities are (2,2,1,1)
and the infinity multiplicities are (2,2,2,2,2,1,1). Therefore

```text
delta_finite+delta_infinity+nodes = 2+5+3=10,
sum m_j=18,       D²=4,        D²-2N=-2.
```

The genus equality independently checks the complete singularity inventory.

## 3. Actual boundary and explicit representation

The seven infinity blowup centres are
L; L∩E1; L∩E2; E3; E4; E5; E5∩E6.
The middle three single-component centres are free. The graph is

```text
L--E3--E2--E1
    |
    E4--E5--E7--E6
            |
          curve
```

In the order L,E1,E2,E3,E4,E5,E6,E7 its weights are
(-2,-2,-2,-2,-2,-3,-2,-1). The negative intersection matrix has
determinant -1, and its marked abelianization gives
(-6,-4,-8,-12,-10,-8,-7,-14) times the curve meridian.
Use the same positive-fibre plumbing convention as the audited predecessor.
With fibres l,a,b,c,d,g,e,f, the vertex relations are

```text
l²=c, a²=b, b²=ac, c²=lbd,
d²=cg, g³=df, e²=f, f=ge mu,
```

together with commutation of adjacent fibres and [f,mu]=1. These are
the relations of the actual marked infinity exterior, obtained from the
displayed curve charts; no full affine presentation is asserted.

Here is a raw zero-based permutation witness:

```text
l=(0,1,3,2,5,4), a=(1,2,0,4,5,3), b=a²,
c=f=identity,
d=(1,2,4,0,3,5), g=(2,4,3,1,0,5),
e=(0,4,3,2,1,5), mu=(1,2,0,3,4,5).
```

Composition is (pq)(i)=p(q(i)). The exact source checks **every original
vertex relation and crossing commutator**, even parity, and the order of
the generated group: it is 360, hence A6. The arrow is a single
three-cycle. Thus this is a genuine surjective marked representation of
the infinity group. It need not factor through the global complement.

The failure mechanism of the previous proof is precise. The two groups
generated by (l,a) and (g,e) both have order 60. The former acts
transitively on six letters, while the latter fixes exactly label five.
They share a five-cycle, whereas the preceding obstruction shared a
three-cycle between groups of order at most twelve. An order-60 subgroup
can have a faithful six-point action. Sharing five moved labels therefore
does not force a common fixed label. This hostile was also recovered
independently by the certificate agent before seeing this producer.

## 4. Next consumer and concept-board update

The inherited odd-cusp theorem still forces any whole-support Keller
realization of (1) to have degree six, generic retained/deleted counts
three, one actual cusp point, all three nodes omitted, and A6 monodromy.
The new boundary witness pays **only** the infinity relations. The next
cheap decisive test is an exact global braid constraint, using the
already audited rational-tube method on the equality predecessor.
Do not claim a finite cover or polynomial source from this representation.

This changes each live lane as follows. Infinity now has a concrete
survivor; the source lane still owes an actual finite map and polynomial
coordinates. The unchanged first-coordinate pair equation supplies a
useful exact collision transport while its residual polynomial drops
from degree four to three. The explicit affine surface and its finite
ramified control distinguish boundary separation from etaleness. The
torsion connection has no newly proved map into this group problem;
it remains an orthogonal source-descent lane.

The [standalone source](../../04-computation/planar_jc48_sep06_next_infinity.py)
and [frozen output](planar_jc48_sep06_next_infinity.out) declare one curve
and one labelled witness, with no census claim. Reproduce by running
the source normally and with `python3 -B -O`; its checks remain active
under optimization. The manifest records the final hashes and gate count.

The [independent analytic and source audit](planar_jc48_sep06_next_infinity_audit.md)
passes, including an independent resultant/discriminant path and a separate
inverse-inclusive group closure. Both source modes match all63 frozen gates.
