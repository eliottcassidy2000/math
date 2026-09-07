# Changing the DG carrier: the complete source-linear layer has no regular conjugate

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
Work over `k=C`. On the inherited explicit surface W, every nonconstant global regular function of degree at most one in the original source variable t fails to have a global regular conjugate with nonzero constant source Jacobian. This includes the entire affine-linear span of the five native global generators `t,b,u,j,v`, not only functions of t. The theorem imposes no degree bound on the proposed conjugate.

The excluded class contains boundary-separating functions with no critical point anywhere on W. For these, either a rational residue already obstructs a primitive, or a rational primitive exists and fails regularity on the two components of a special fibre. This is a restriction on a specified class of coordinates on one explicit surface, not a general Keller exclusion.

## 1. Inheritance, operation, and scope

The immediate proved suppliers are the incoming [DG surface](planar_jc48_sep06_dg_surface.md), [finite quartic map and relative-primitive obstruction](planar_jc48_sep06_dg_finite_map.md), and [discrete carrier rigidity](planar_jc48_sep06_discrete_carrier.md), read at the integrated baseline `89aca0881928`.

The exact object and both complete charts are

```
W = (P1_s x P1_z) \ {z=s^2},
U0 = A2_(x,t),       Uinfty = A2_(r,b),
x = 1/r,            t = -r^2-r^4b,
b = -x^2-x^4t,      D = {r=0},      W\D=U0,
omega = dx wedge dt = r^2 dr wedge db.
```

The global ring is exactly the ring of functions regular in both charts. In particular its restriction to the dense U0 is contained in `k[x,t]`. The incoming finite map `(t,b)` and its full integral basis prove that the surface and its boundary separators are genuine global objects; the calculation below uses the two charts directly and requires no assumed Keller realization of W.

The old carrier `k[u,p,y]` collapses D. The new operation is to change the first-coordinate carrier to global functions that include the separator b. We retain

```
u=x^2t,       j=xt,       v=x(1+x^2t),       b=-x^2-x^4t.
```

The closest inherited negative result excludes first coordinates in `k[t]`. The inherited hostile is that `omega=d(x dt)` globally although x itself is not globally regular; exactness of a global form does not give a prescribed global relative primitive. The corrected near miss is to infer a Keller pair from boundary separation, a rational conjugate, or a critical-free first coordinate.

The least-used sidecars now restored are the full source-t degree, the two components of a special fibre, and residues of the actual relative one-form. [THM-2740, polynomial-coordinate first-target triangularity](../../01-canon/theorems/THM-2740-polynomial-coordinate-first-target-triangularity-and-mixed-resolvent-shear-closure.md) is a related proved result: it assumes the first function is already a polynomial coordinate. Our descent lemma below does not make that assumption and shows why the new critical-free hostile is not such a coordinate. No literature-priority claim is made.

## 2. Exact classification of the entire global source-linear layer

**Theorem 1.** Put `L={F in O(W): deg_t(F|U0)<=1}`. Then

```
L = span_k(1,t,j,u,v,b).                                      (1)
```

More explicitly, all its elements, uniquely, are

```
F=A(x)t+B(x),
A(x)=a0+a1*x+a2*x^2+a3*x^3+a4*x^4,
B(x)=a4*x^2+a3*x+B0.                                        (2)
```

Its full second-chart expression is

```
F = B0-a2-a1*r-a0*r^2
    -b*(a4+a3*r+a2*r^2+a1*r^3+a0*r^4).                      (3)
```

**Proof.** Begin with arbitrary polynomials A,B, without a degree restriction on x. On the second chart the coefficient of b is `-r^4 A(1/r)`. Its regularity forces `deg A<=4`, since no term from B can cancel a b-dependent pole. The b-independent expression is `B(1/r)-r^2 A(1/r)`. Its negative powers vanish precisely when `B=a4*x^2+a3*x+B0`. Conversely these conditions give the polynomial expression (3). This proves the whole intersection, rather than just exhibiting six regular functions. Linear independence follows from the distinct coefficients in A and B0. The formula in generators is

```
F=a0*t+a1*j+a2*u+a3*v-a4*b+B0.
```

In particular, a member with A identically zero is constant. Its actual boundary restriction is

```
F|D=B0-a2-a4*b.                                             (4)
```

Thus the part `a4!=0` really separates D and escapes every boundary-collapsing algebra. That does not settle the conjugate question.

## 3. An all-degree source descent lemma and the global conclusion

Use `J(F,g)=F_x g_t-F_t g_x`. The following elementary lemma is valid for arbitrary `A!=0,B in k[x]`, independently of global regularity on W.

**Lemma 2.** If `F=A(x)t+B(x)` and `g in k[x,t]`, then

```
J(F,g) in k[x]  iff  g=h(F)+q(x),    h in k[T], q in k[x].    (5)
```

Consequently

```
J(F,k[x,t]) intersect k[x] = A(x) k[x],
ker J(F,-) = k[F].                                          (6)
```

Here the first intersection refers to the image of the linear Jacobian operator. It is an equality, not merely an ideal containment.

**Proof.** Write the highest positive t-degree of g as n, with coefficient `g_n(x)`. The coefficient of `t^n` in its Jacobian is

```
n A'(x) g_n(x) - A(x) g_n'(x).
```

If the Jacobian has no t-dependence, this is zero. Equivalently the derivative of `g_n/A^n` in `k(x)` is zero, so `g_n=c_n A^n` for a scalar `c_n`. Subtracting the full polynomial `c_n F^n` preserves the Jacobian and lowers the t-degree. Repetition terminates and gives (5), with no bound on the initial n or the x-degrees. Conversely `J(F,h(F)+q)=-A q'`. Differentiation is surjective on `k[x]` in characteristic zero; this proves the image equality and the kernel statement.

**Theorem 3.** For every nonconstant `F in L`, there are no `g in O(W)` and `lambda in k*` such that

```
dF wedge dg = lambda omega.                                (7)
```

**Proof.** Restriction to U0 gives `J(F,g)=lambda` with `g in k[x,t]`. Equations (5)-(6) force A to be a nonzero constant, since a nonconstant polynomial A cannot divide a nonzero constant. In the global classification (2), this also forces B to be constant, hence `F=a0*t+B0`. All its affine-plane polynomial mates have the form

```
g=-(lambda/a0)x+h(F).
```

On Uinfty the first term has a nonzero simple pole `-(lambda/a0)/r`, while the second is a polynomial in the global function F. The pole cannot cancel. This is incompatible with `g in O(W)`.

For comparison, `(a0*t+B(x),-(lambda/a0)x)` is an entirely valid polynomial constant-Jacobian pair on the affine plane for arbitrary B. Global regularity of both coordinates on W is essential to the theorem. Nor does the argument assert that a general polynomial endomorphism is an automorphism.

Any nonconstant outer polynomial `Phi(F)` is excluded as a first coordinate as well: the equation `Phi'(F)J(F,g)=lambda` in `k[x,t]` forces `Phi'` constant. Pullback of this excluded family under an independently certified omega-preserving automorphism of W is also excluded, simply by applying its inverse to a proposed mate. This last statement requires an actual global automorphism; no formal or rational replacement is implicit.

## 4. A complete critical-point sidecar in this layer

The simple mixed pencil already fails: for `F=t+lambda*b`, `lambda!=0`,

```
F=(1-lambda*x^4)t-lambda*x^2
```

has four affine critical points `x^4=1/lambda`, `t=-1/(2x^2)`. But this elementary test does not cover the whole layer.

**Proposition 4.** The nonconstant members of L with no critical point anywhere on W are precisely

```
F=a(x-h)^2(x-k)^2 t + a[x^2-2(h+k)x] + B0,
a*h*k != 0,                                                (8)
```

where h=k is permitted.

**Proof.** On U0, criticality requires `A(x)=0` and `A'(x)t+B'(x)=0`. Every simple root of A supplies a critical point. At a multiple root it is avoided precisely when B' is nonzero there. If `deg A=4`, all roots must therefore be multiple, so its factorization is `a(x-h)^2(x-k)^2`. The coefficient constraint (2) gives the stated B. Its derivative at h is `-2ak`, and at k is `-2ah`; these are nonzero exactly when h and k are nonzero. On D, (3) has `F_b=-a4=-a`, so there is no boundary critical point.

For completeness, the other degree cases cannot be globally critical-free. A nonzero constant A gives F linear in t, and dF vanishes on all D. Degrees one and two force an affine critical point, either at a simple root or because B'=0 at a repeated root. In degree three, absence of affine critical points requires `A=a3(x-h)^3`, with `B'=a3!=0`. On D, however, `F_b=0` and `F_r=-a1-a3*b`, which vanishes at `b=-3h^2`. These exhaust all possibilities over C.

Thus (8) gives a full positive submersion family inside the class excluded by Theorem 3. Absence of critical points is genuinely insufficient here.

## 5. Distinct roots: a rational residue obstruction

For any `A!=0`, the function field is `k(x,t)=k(F)(x)`. On a generic fibre, the relative equation `J(F,g)=lambda` becomes

```
partial_x g |F = -lambda/A(x).                             (9)
```

For (8) with `h!=k`, the residue at x=h of `-1/A(x)` is

```
2 / [a(h-k)^3] != 0.                                      (10)
```

A derivative of a rational function has zero simple-pole coefficient in every Laurent expansion, including over the coefficient field `k(F)`. Therefore these globally critical-free, boundary-separating members have no rational conjugate at all. This is an exact period obstruction; polynomial degree tests are unnecessary.

## 6. Fourfold root: a rational primitive and its full special-fibre failure

When h=k!=0, write

```
F=a(x-h)^4 t + a(x^2-4hx) + B0,
f0=B0-3ah^2.
```

There is a genuine rational conjugate

```
g0=1/[3a(x-h)^3],       J(F,g0)=1.                        (11)
```

Moreover every rational solution with this bracket is `g0+H(F)` for `H in k(T)`, since the constants of `partial_x` in `k(F)(x)` are exactly `k(F)`. A generic fibre of F is a projective x-line with the single point x=h removed: the infinity point belongs to D. Thus the generic fibre is A1, so the previous high-genus obstructions do not apply.

The decisive retained special fibre is

```
F-f0=a(x-h)[(x-h)^3 t+x-3h].                              (12)
```

It has two disjoint irreducible components on W. The first is `E1={x=h}`, an affine line, contained in U0. The second is the closure E2 of `(x-h)^3 t+x-3h=0`; they are disjoint because the second factor at x=h equals `-2h!=0`.

Here is the complete E2, including its boundary point. With `z=1/(x-h)` on U0,

```
t=2h*z^3-z^2,
x=h+1/z,
r=z/(1+hz),
b=-2h^5*z^3-7h^4*z^2-8h^3*z-3h^2.                       (13)
```

The original and boundary charts cover this parametrization; z=0 is its D point `b=-3h^2`. This identifies E2 with A1, rather than retaining only its punctured source-chart part. Both components occur with multiplicity one.

Along E1, `F-f0` has order one and g0 has order minus three. To make `g0+H(F)` regular there, H must have a pole of order three at f0. But along E2 the same `F-f0` also has order one, whereas g0 is regular: it equals `z^3/(3a)` there. The required H therefore introduces an uncancelled pole on E2. No rational fibre-constant correction can repair both components. This is an independent geometric explanation of Theorem 3 in the generic-A1 case.

The smallest clean control is

```
F=(x-1)^4t+x^2-4x=t-4j+6u-4v-b,
F|Uinfty=-6+4r-r^2-b(1-r)^4.
```

It separates D by `F|D=-6-b` and is globally critical-free: at its only zero of F_t, x=1, the derivative F_x is -2, while on D its derivative F_b is -1. Nevertheless its rational conjugate `1/[3(x-1)^3]` cannot be made global. Its special value is -3. This is a sharp hostile to replacing global regularity by smoothness of the first-coordinate map, rational relative exactness, generic A1 geometry, or boundary separation.

## 7. Concept board and stopping point

| Lane | Actual transfer and preserved object | Lost data restored / next question |
|---|---|---|
| Anchor: changed global carrier | The two-chart intersection classifies every source-t-linear global first coordinate; top-degree descent excludes every polynomial mate with no degree bound | A higher-source-t-degree global coordinate remains outside the theorem |
| Niche: rational relative primitive | The function-field map `t=(F-B)/A` converts the bracket equation to the exact one-form `-dx/A` | Distinct-root residues retain the period that smoothness forgets |
| Wildcard: generic A1 escape | Fourfold A admits an exact rational primitive and a globally critical-free F | The two separate special-fibre pole orders obstruct global regularity |
| Incoming surface atlas | The global separator b removes the old carrier's boundary collapse | Separation is insufficient for a constant-Jacobian pair; the original finite quartic remains a positive finite-map control |

The new boundary is exact: a possible global constant-Jacobian pair on this specific W must have both coordinates of source-t degree at least two in the fixed inherited chart. The theorem does not classify that remaining layer, finite maps of arbitrary coordinates, actual six-sheet envelopes, or JC(2). It does not infer non-rationality of a flow from non-local-nilpotence.

## 8. Exact verification and frozen reproduction

The standalone [source](../../04-computation/continuing10_20260907_dg_linear_carrier.py) imports no previous producer. It checks the symbolic two-chart classification and full generator span; independent Laurent coefficient ranks through x-degree eight; the actual highest-t identity through degree six; fifteen unpruned exact polynomial coefficient systems with positive constant-A controls; both whole-parameter submersion derivatives and the nonzero rational residue; the rational conjugate; the complete special-fibre factorization; both charts of its A1 component and actual D point; and the literal h=1 hostile. The finite checks validate mechanisms and hostiles. The all-degree claims use the proofs above.

```
python3 -B 04-computation/continuing10_20260907_dg_linear_carrier.py
python3 -B -O 04-computation/continuing10_20260907_dg_linear_carrier.py
```

Both runs must give identical actual LF output. The parent owns independent audit, status promotion, repository filing and Git; this report changes no maintained file or scarce identifier.

Both frozen runs pass **88 always-active exact gates**, with byte-identical actual LF stdout.

- Source SHA256: `0f803f5ea98256bb8d4dd2685732885b60b90ed4514b3a5deb8d02a46f85241e`.
- Normal and optimized stdout SHA256: `3bb9dc940479599197a0fea02269c6485d80e80b1be26282c298e72452e45a42`.

The [independent audit](continuing10_20260907_dg_linear_carrier_audit.md) accepts the scoped proof.
Repository filing changes only routing and status prose; frozen source,
output and certificate bytes match the independently reviewed originals.
