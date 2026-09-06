# An exact (4,6) curve family: nodes, cusps, triples and tangencies

**Status: PROVED + independently audited exact algebraic controls.**
The [root audit](planar_jc48_sep06_curve_probe_audit.md) includes a second
independent analytic review. This classifies the explicit one-parameter
family below, not of all intrinsic `(4,6)` curves and not a construction of
a Keller map. JC(2) remains OPEN.

**Current boundary update:** the later independently audited
[cusp passport](planar_jc48_sep06_cusp_passport.md) excludes the three cusp
parameters left open by the argument below. Together the two notes exclude
every parameter in this explicit family as a sole nonproperness support.
The three-parameter OPEN statements below describe this note's original
intermediate boundary, which the linked result now supersedes.

## 1. Target, inheritance and the complete collision algebra

Let `lambda` range over `C`, and let `C_lambda` be the image of

```text
nu_lambda(t)=(U(t),V_lambda(t))
            =(t^4+t, t^6+t^2+lambda t).                 (1)
```

The source of the question is the recovered nodal nonproperness exclusion
in [planar_jc48_sep06_infinity.md](planar_jc48_sep06_infinity.md), together
with the intrinsic normalization degrees in
[THM-4122](../../01-canon/theorems/THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction.md).
The least-used datum restored here is **which collision pairs share a
parameter and a target image**. A discriminant of a projected resultant
alone loses that information.

The working board is: finite normalization; divided differences; unordered
pair sums; projected collision multiplicity; tangent/critical jets; actual
target fibres. The initial positive control is `lambda=1`. The hostile is
`lambda=-1`: the collision scheme is reduced, while its projection has
repeated roots caused by a triple image. The nearby nonnodal boundary
`lambda=2` instead has a nonreduced collision scheme and an ordinary
tangency. All coordinates and coefficient laws in (1) are retained.

Define the divided differences

```text
N(s,t)=(U(s)-U(t))/(s-t)
      =s^3+s^2t+st^2+t^3+1,
M_lambda(s,t)=(V_lambda(s)-V_lambda(t))/(s-t)
      =s^5+s^4t+s^3t^2+s^2t^3+st^4+t^5+s+t+lambda.
```

For `p=s+t`, `q=st`, these become

```text
N=p^3-2pq+1,
M=p^5-4p^3q+3pq^2+p+lambda.
```

On the collision scheme `p` is a unit, `q=(p^3+1)/(2p)`, and `4pM=-H`,
where

```text
H_lambda(p)=p^6+2p^3-4p^2-4lambda p-3.
```

Consequently the full divided-difference algebra is exactly

```text
C[s,t]/(N,M)
  = C[p,p^(-1),s]/(H_lambda(p),
                    s^2-ps+(p^3+1)/(2p)).               (2)
```

The inverse is `t=p-s`; thus this is an isomorphism, not just an eliminant.
Since `H(0)=-3`, (2) has complex dimension twelve for **every** `lambda`.
The ordered-pair swap fixes `p` and exchanges the two quadratic roots. Its
fixed points are diagonal pairs, detected by

```text
Delta_pair(p)=p^2-4q=-(p^3+2)/p.                        (3)
```

The image in (1) is an irreducible closed affine curve. The map is finite:
`t` satisfies the monic equation `t^4+t-U=0`, so `C[t]` is finite over
`C[U,V]`. The finite collision scheme shows that a generic parameter has
no distinct partner; hence the finite map is birational and is the
normalization. Its intrinsic pole pair is exactly `(4,6)` for every
`lambda`, including the exceptional values below.

## 2. Three disjoint exceptional divisors on the parameter line

Put

```text
C(lambda)=128lambda^3-288lambda-283,
E(lambda)=(lambda+1)(lambda^3-lambda^2+3lambda+1),
T(lambda)=(lambda-2) B(lambda),
B(lambda)=3125lambda^5+6250lambda^4-6625lambda^3
          -24604lambda^2-18669lambda-4371.              (4)
```

These polynomials are squarefree and pairwise coprime. Thus they have,
respectively, three, four, and six distinct complex roots, with no overlaps.
The source verifies the exact polynomial identities and Euclidean gcds
used in this section over `Q`, without numerical root classification.

**Critical parameters.** Direct elimination gives

```text
Res_t(U',V_lambda')=-8 C(lambda),
Res_p(H_lambda,p^3+2)=C(lambda).                         (5)
```

If `U'(t)=V'(t)=0`, then

```text
4t^3+1=0,             lambda=(3/2)t^2-2t.
```

The three critical parameters have distinct values of `lambda`: equality
for two different roots would force their sum to be `4/3`, making the
third root `-4/3`, which does not satisfy `4t^3+1=0`. Thus exactly one
critical parameter lies above each root of `C`. At such a point `U''!=0`,
and the determinant of the second and third derivative vectors is, modulo
`4t^3+1`,

```text
U'' V'''-V'' U''' = -12t(15t+4) != 0.                  (6)
```

The nonvanishing follows from its gcd one with `4t^3+1`. The branch is an
ordinary cusp, since its first nonzero independent jets have orders two
and three.

**Triple images.** This step retains actual fibre coordinates, following
root's independent coefficient comparison. Suppose three distinct
parameters share a target `(a,b)`. Reducing `V_lambda-b` modulo `U-a` gives

```text
-t^3+(a+1)t^2+lambda t-b.
```

The monic cubic with those three roots divides the quartic `U-a`. Set
`A=a+1`. Comparison with

```text
(t+A)(t^3-A t^2-lambda t+b)=t^4+t-a
```

forces, and is forced by,

```text
lambda=-A^2,      b=1-A^3,      A^4-2A+1=0.             (7)
```

Elimination of `A` gives exactly `E(lambda)=0`; explicitly
`A=(lambda^2+1)/2`, so each such parameter has one possible triple target.
Conversely, the cubic in (7) has three distinct roots: its discriminant
is coprime to `A^4-2A+1`. The fourth quartic root `-A` is distinct because
the cubic evaluated there is `1-4A^3`, also coprime to that quartic in `A`.
It follows that (7) supplies exactly three preimages of `(A-1,1-A^3)`.
No target fibre can have four distinct parameters, because the displayed
nonzero cubic remainder cannot vanish at four distinct roots of `U-a`.

**Tangency parameters.** The pair-sum polynomial satisfies

```text
Disc_p H_lambda=4096 T(lambda).                         (8)
```

At a root of `T`, `H` has exactly one double root and four simple roots.
Here is a certificate of that assertion without relying on a generic
discriminant heuristic. The next-to-last nonzero subresultant of `H,H'` is

```text
2560 [ a(lambda) p+b(lambda) ],
a=375lambda^4-1134lambda^2-813lambda-113,
b=315lambda^3-144lambda^2-907lambda-405,
```

and `gcd(a,T)=1`; the last subresultant is `-4096 T`. Therefore the gcd of
`H,H'` has degree exactly one at every root of `T`, which is precisely one
double root and no other multiple root. Since `C` is nonzero there, none
of the quadratic roots in (2) coincide. The double pair-sum root produces
two swapped collision points, each of local length two; all others have
length one.

## 3. Complete affine singularity classification for this family

At a collision with `s!=t`, the Jacobian determinant of the two divided
differences, restricted to their common zero set, is

```text
det d(N,M) = -D(s,t)/(s-t)^2,
D(s,t)=U'(s)V'(t)-U'(t)V'(s).                           (9)
```

For immersed branches, local collision length one is exactly transverse
intersection; length two is ordinary tangency. To see the latter assertion,
use a target coordinate noncritical on either branch, write them locally
as graphs, and eliminate one of the two source parameters. The local
collision length is the vanishing order of the difference of the graphs.
Length two therefore gives distinct quadratic terms after their common
tangent, the ordinary tacnode. This argument uses immersion; (6) treats
the separate critical case.

If `H` is squarefree and (3) is nonzero at its roots, (2) is a product of
twelve copies of `C`. Thus every ordered collision is reduced, and (9)
gives distinct tangents for every pair. The triple criterion (7) determines
whether three unordered pairs represent one common image. With the
critical points separately counted, all affine singularities are therefore:

| Parameter condition | Complete affine singularity inventory |
|---|---|
| `C E T != 0` | Six ordinary nodes |
| `C=0` | One ordinary cusp and five ordinary nodes |
| `E=0` | One ordinary triple point and three ordinary nodes |
| `T=0` | One ordinary tacnode and four ordinary nodes |

For the cusp row, one simple root of `H` has a double quadratic root in
(2), yielding a diagonal length-two point; the other five pair sums give
the five nodes. It is the unique critical parameter from (5). The cusp
cannot share its image with another branch: the fibre polynomial would
then have a double root at the cusp and one distinct root, so the cubic
remainder again divides the quartic and forces (7). Its cubic has distinct
roots, a contradiction. Equivalently `C` and `E` are disjoint.

For the triple row all six unordered collisions are transverse. Three of
them occur at the unique triple target from (7), and the others are three
distinct double targets. For the tangency row, `C E != 0` gives immersed
branches and no triples; (8) gives exactly one length-two unordered
collision and four simple ones. Finally a unique immersed parameter maps
to a smooth curve point by the analytic constant-rank theorem. Since the
normalization is finite, this accounts for every affine singularity; none
can be hidden outside the collision/critical equations.

The first row is excluded as a **sole** Keller nonproperness curve by the
inherited nodal theorem. Section 5 supplies a separate Euler argument that
also excludes the triple and tangency rows, leaving only the three cusp
parameters as OPEN controls. No curve in this family is asserted to occur
as a Keller nonproperness curve.

## 4. Named rational controls and the projection-discriminant hostile

The exact projection resultant is

```text
R_lambda(s)=Res_t(N,M_lambda)
 =3s^12-2lambda s^10+6s^9+7s^8+(4-4lambda)s^6
  +(10-2lambda^2)s^5+(5+4lambda)s^4+2s^3
  +(3-lambda^2-2lambda)s^2+(4-2lambda^2+2lambda)s
  +2+3lambda-lambda^3.
```

An independent reconstruction through (2) is

```text
Res_p(H_lambda,2ps^2-2p^2s+p^3+1)=-64 R_lambda.
```

Its discriminant retains all three mechanisms, with different powers:

```text
Disc_s R_lambda = -12288 C(lambda) E(lambda)^6 T(lambda)^2. (10)
```

At `lambda=1`, the collision ideal has lexicographic basis

```text
4t+3s^9+3s^8-2s^7+s^6+10s^5+3s^4-s^3+2s^2+6s+4,
R_1=3s^12-2s^10+6s^9+7s^8+8s^5+9s^4+2s^3+4s+4.
```

The second polynomial is squarefree. Adding either `s-t` or `D` to the
ideal gives the unit ideal. Thus it has twelve distinct off-diagonal
ordered collisions, one partner per parameter, and six node images.

At `lambda=2`,

```text
H_2=(p+1)^2(p^2-p-1)(p^2-p+3),
R_2=s^2(s+1)^2(s^4-s^3+2s-1)
                   (3s^4-3s^3+2s^2-2s+5).
```

The collision ideal has a `t`-linear lex basis, and adjoining `D` gives
`(s+t+1,s^2+s)`. Thus the only tangent collision is the swapped pair
`t=0,-1`, with common image `(0,0)` and tangent slope two.
The exact identity `V_2-2U=t^2(1-t^2)^2` gives quadratic graph coefficients
`1` and `4/9` on those two branches, so the contact is exactly two.

At `lambda=-1`, the common fibre `(0,0)` is exactly
`gcd(U,V_{-1})=t(t^2-t+1)`. Its three roots are distinct. At `t=0` the
tangent vector is `(1,-1)`; at a root `z` of `z^2-z+1` it is
`(-3,5-4z)`. The determinants are `2-4z` and `12(z_2-z_1)`, all nonzero.
The projected resultant factors as

```text
R_{-1}=s^2(s^2-s+1)^2
       (3s^6+6s^5+5s^4+4s^3+9s^2+10s+4).
```

Yet adding `D` to `(N,M_{-1})` gives the unit ideal, so the zero-dimensional
collision scheme is reduced. The repeated projection roots occur because
each of the three parameters has two different partners. This is an exact
hostile to interpreting every zero of `Disc R` as a tangency or a
nonreduced collision. The pair algebra and actual fibre sidecar separate
the mechanisms.

## 5. One arbitrary nonnodal point: retaining its actual fibre

Here is a direct conditional consequence of the page/Euler mechanisms that
does not assume a normal-crossing specialization rule at the exceptional
point. Suppose the **whole** nonproperness set of a hypothetical planar
Keller map is irreducible with normalization `A1`. Assume all singularities
are ordinary nodes except one point `z`, at which there are `r` local
branches of any type. Write `N` for the number of ordinary nodes,
`d>1` for mapping degree, `a` for the generic actual affine degree, and
`delta=d-a`, so `1<=delta<=d-1`. Retain the actual count
`n_z=#F^(-1)(z)` and the ordinary-node missing overlaps `omega_p>=0`.

The normalization has `2N+r` marked preimages, so

```text
chi_c(S\({z} union nodes))=1-2N-r,
chi_c(A2\S)=N+r-1.
```

At each ordinary node, exact specialization gives actual count
`d-2delta+omega_p`. Constructible Euler integration now yields

```text
1=d(N+r-1)+a(1-2N-r)
     +sum_nodes(d-2delta+omega_p)+n_z
 =(r-1)delta+n_z+sum_nodes omega_p.                     (11)
```

No actual-versus-formal intersection equality at `z` was inserted. Its
full fibre count remains visible in (11).

If `r>=3`, its right-hand side is at least two, impossible. If `r=2`,
nonnegativity forces

```text
delta=1,          n_z=0,          every omega_p=0.       (12)
```

The generic deleted-boundary length over the unique component is one.
Each boundary prime contributes a positive integer `e_D f_D`, so all
boundary ramification indices equal one. This is a statement about
weighted generic lengths, valid for every `d`, not a count of special
boundary points. The finite envelope is unramified in codimension one,
hence finite etale by branch purity, and connectedness over simply
connected `C^2` forces degree one. This contradicts `d>1`.

> Thus a whole irreducible nonproperness curve with only ordinary nodes
> and one additional singularity having at least two local branches is
> impossible for a nonautomorphic planar Keller map.

The argument is an explicit specialization of the inherited Euler/page
framework, with a new use here of the exceptional point's actual fibre.
Paper II Section 5 already has the general conductor/overlap ledger; we
make no external priority claim. The stronger single-point conclusion
above does not follow merely by setting its conductor term to zero, nor
does it assert that an arbitrary nonproperness curve has this singularity
inventory. Root independently checked (11), including the degree-independent
weighted-length purity step.

For (1), the ordinary triple has `r=3` and the tacnode has `r=2`. These
exclude all four roots of `E` and six roots of `T` as sole supports. Together
with the six-node exclusion, the only remaining values in this explicit
one-parameter family are

```text
C(lambda)=128lambda^3-288lambda-283=0.                  (13)
```

These three cases have one cusp and five nodes and remain OPEN here. At
the cusp `r=1`, so (11) gives the necessary ledger

```text
n_cusp+sum_(five nodes) omega_p=1.                      (14)
```

Consequently the cusp has zero or one actual affine preimage. If it has
one, every node overlap vanishes; if it has zero, exactly one node has
overlap one. Since each node satisfies
`omega_p>=max(0,d-2a)`, (14) also forces `d<=2a`. None of these necessary
conditions establishes existence or excludes the remaining cusp cases.

## 6. Verification and scope

The source checks universal polynomial identities, the explicit exceptional
polynomials and their gcds, the linear subresultant, the critical jets and
triple fibre factorization. Independently of the pair-sum presentation it
computes the ordinary divided-difference Groebner ideals at the three named
rational parameters, including tangent and diagonal ideals. Every gate is
always active under optimized Python. No numerical root search is used.

```sh
python3 04-computation/planar_jc48_sep06_curve_probe.py
python3 -O 04-computation/planar_jc48_sep06_curve_probe.py
```

Normal and optimized replay are byte-identical: **214 always-active gates**.
The source retains all universal polynomial laws and the three named
direct-ideal controls; its bounded extra passport universe is
`2<=d<=12`, `1<=delta<=d-1`, `2<=r<=5`, and `0<=n_z,sum omega<=2`.
The five-node inequality is additionally controlled on `2<=d<=25` and
`1<=a<=d-1`. Neither finite range proves the all-degree Euler argument.

```text
source SHA256 97c28adbf35451a8e521486680a2b79267073430fdc8fe69a96e9435ffa923a1
output SHA256 9d227fc32a01f24d13105ea4938cd0485c00ffe95225c723df4e6e293c4f786f
```

The saved output is [planar_jc48_sep06_curve_probe.out](planar_jc48_sep06_curve_probe.out).
The family exclusion leaves exactly the three OPEN parameters (13) and
concerns only (1) as the whole support. It supplies neither an entry theorem
for arbitrary Keller pairs nor an exclusion of all intrinsic `(4,6)` curves.
