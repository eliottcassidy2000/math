# Independent audit: the complete source-linear DG carrier has no global mate

**Status: PASS; complete analytic proof accepted with independent exact
controls. No producer repair requested.** This audits
`continuing10_20260907_dg_linear_carrier.md` on the specified surface over C.
It does not promote a theorem about arbitrary global first coordinates or a
general Keller/JC exclusion.

## 1. Object, inherited hypotheses and main proof

I read the frozen producer report and source, the integrated
`planar_jc48_sep06_dg_surface.md` and `planar_jc48_sep06_dg_finite_map.md`, and
the relevant discrete-carrier scope. The actual two charts cover the whole
surface:

```
U0=Spec C[x,t], Uinfty=Spec C[r,b],
x=1/r, t=-r²-r⁴b, b=-x²-x⁴t, D={r=0}.
```

Thus a global regular function restricts to a genuine polynomial on U0.
Restriction is injective because U0 is dense. Regularity is equivalently
polynomiality in both charts; no chosen subalgebra is substituted for the
full global ring. The inherited classical name of the surface is not a
dependency of any new calculation in this audit.

For `F=A(x)t+B(x)`, the b-dependent part of its second expression is
`-b r⁴A(1/r)`. No b-independent term can cancel its poles, so deg A<=4.
The remaining negative Laurent coefficients force exactly
`B=a4 x²+a3 x+B0`. This yields the six-dimensional layer and the stated
boundary restriction. In particular A=0 gives only constants. This is an
all-x-degree argument, not a conclusion extrapolated from a rank bank.

The all-degree descent lemma is sound for arbitrary nonzero A and arbitrary
B, even outside this global layer. If n>0 is the t-degree of a proposed g,
the coefficient of t^n in J(F,g) is `n A' g_n-A g_n'`. Its vanishing says
`(g_n/A^n)'=0` in C(x). A rational function with zero derivative in
characteristic zero is constant, so `g_n=c A^n`. Subtracting the **full**
polynomial cF^n cancels its top t-degree and preserves J. Induction gives

```
J(F,g) in C[x] iff g=h(F)+q(x),
J(F,h(F)+q)=-Aq'.
```

Because polynomial differentiation is onto C[x], the image intersection is
exactly A C[x], and the kernel is C[F]. A nonzero constant Jacobian forces A
constant. Within the classified global layer B is then constant too, and
every polynomial mate contains a nonzero multiple of x. On Uinfty that term
has a simple pole; adding a polynomial in the global F cannot cancel it.
This proves the full global conclusion without a bound on mate degree.

The outer-polynomial extension follows from the ordinary chain rule in
C[x,t]; a nonconstant Phi'(F) cannot divide a nonzero constant. The stated
transport by an actual omega-preserving global automorphism is also valid.
It does not license transport by a formal automorphism or a rational chart
change.

An important type boundary is retained: `omega=dx wedge dt` extends as
`r² dr wedge db`, so it vanishes to order two on D. A constant **source**
Jacobian equation is not a claim that the corresponding global map is
everywhere etale. The theorem excludes the stated mate equation directly.

## 2. Critical-free classification and rational residue

At a root xi of A, an affine critical point exists if A'(xi) is nonzero, or
if both A'(xi) and B'(xi) vanish. Thus all roots must be multiple and B' must
be nonzero at each if the affine critical locus is empty.

In degree four this leaves exactly `a(x-h)^2(x-k)^2`; the two values of B'
are -2ak and -2ah. The condition ahk!=0 is both necessary and sufficient,
including h=k. On D the derivative with respect to b is -a, so there is no
additional critical point. Degrees zero through three are correctly
discharged: the only affine-critical-free cubic possibility has a triple
root, but its first normal derivative vanishes at b=-3h² on D. The lower
degrees either have an affine critical point or, for constant A, have dF=0
on all D. This proves the entire globally critical-free sidecar over C.

For h!=k, the field identification C(x,t)=C(F)(x) is exact because A!=0.
The mate equation becomes the rational derivative equation
`partial_x g|F=-lambda/A`. The coefficient of `(x-h)^-1` in `-1/A` is the
derivative at h of `-1/[a(x-k)^2]`, namely `2/[a(h-k)^3]`. It is nonzero;
derivatives of rational functions have no simple Laurent-pole coefficient.
Consequently there is no rational mate in the distinct-root case. This is
an independent obstruction, not needed to make the polynomial descent work.

## 3. Fourfold root: complete special fibre and pole obstruction

For `F=a(x-h)^4t+a(x²-4hx)+c`, let f0=c-3ah² and assume ah!=0. Directly,
`g0=1/[3a(x-h)^3]` satisfies J(F,g0)=1. Every rational solution differs by
H(F), because the constants of partial_x on C(F)(x) are exactly C(F).

The source factorization is

```
F-f0=a(x-h)[(x-h)^3t+x-3h].
```

The second factor at x=h is -2h, so the components are disjoint on U0 and
both occur simply. Independently, put w=1-hr. The full second-chart equation
is

```
F-f0=-a w Q,
Q=h²(3-hr)+b(1-hr)^3.
```

At w=0 the second factor is 2h², nonzero. This proves disjointness on the
other chart as well. At r=0, Q=3h²+b, so there is exactly the stated simple
E2 boundary point, b=-3h². No divisor supported on the boundary was omitted.

The coordinate z=1/(x-h) on E2 in U0 agrees with r/(1-hr) in Uinfty. The
denominator is a unit on the latter component: the explicit identity

```
w * (-h²-bw²)/(2h²) = 1-Q/(2h²)
```

holds. Eliminating z independently from

```
r(1+hz)-z=0,
b-(-2h⁵z³-7h⁴z²-8h³z-3h²)=0
```

gives exactly Q. The two pieces are the z-line localized at z and at 1+hz.
These opens cover A1 because `(1+hz)-hz=1`. In particular z=0 is the retained
boundary point; z=-1/h belongs to U0 with x=0. This verifies the whole E2,
rather than its punctured first-chart part.

The generic fibre is also correctly A1. At generic level Y the same
parameter has

```
x=h+1/z,
t=(Y-f0)z⁴/a+2hz³-z².
```

Its b-expression is polynomial in z, and at z=0 takes the boundary value
`(c-6ah²-Y)/a`. The one missing point of the projective x-line is x=h.
Thus no high-genus argument is being silently used for this family.

The mate g0 has pole order three on E1 but is regular on E2, where its value
is z³/(3a). The function F-f0 has order one on both components. Regularity
of g0+H(F) on E1 would force H to have pole order three at f0; its necessary
leading coefficient is `8a²h³/3`, nonzero. On E2 that same H(F) then has
order minus three and cannot cancel the regular g0. This proves the global
failure of every rational correction. No assertion that cancellation on E1
alone is already sufficient is needed.

## 4. Independent finite controls and reproduction

The independent source imports no producer. Its paths are all new audit
files. It uses native sparse Q[x,t] operations and standard-library rational
Gaussian elimination for:

- complete negative-Laurent coefficient matrices through arbitrary-control
  x-degree 16, recovering dimension six from degree four onward;
- 24 unpruned polynomial-mate coefficient systems on six fresh A,B cases,
  with t-degree bounds 0,1,2,4 and x-degree eight;
- exact full-power image identities through order seven;
- all 35 multisets of roots from {-2,0,3} in degrees zero through four,
  including repeated zero roots and the missing-boundary cubic hostile.

The universal parameter checks use direct local coefficient differentiation
for the residue, a separate resultant for E2, its coordinate-ring inverse and
two-open cover, the generic fibre parametrization, and both divisor orders.
The affine-plane positive pair and the source volume's order-two vanishing
are also checked. The finite controls corroborate the analytic proofs; they
do not replace the all-degree descent or the whole-parameter classification.

```
python 04-computation/continuing10_20260907_dg_linear_carrier_audit.py
python -O 04-computation/continuing10_20260907_dg_linear_carrier_audit.py
```

Both modes pass **138 always-active exact gates**. Their actual LF outputs
are byte-identical, 1083 bytes. SymPy is used only for the stated symbolic
parameter checks; the sparse polynomial and rational-rank engines are
independent standard-library implementations.

Producer pins at review:

```
source 0f803f5ea98256bb8d4dd2685732885b60b90ed4514b3a5deb8d02a46f85241e
report fc9a081173fe48167b7f16de784c88e8ed014bc72c8779275f45a5234cb73448
output 3bb9dc940479599197a0fea02269c6485d80e80b1be26282c298e72452e45a42
```

Independent final pins:

```
source 55a5f0b24ff364b9e96767bac9bd57e63c319315375b5d0565526234089ca850
output 1d145d652dee4c15f56f990afd98e9d025d9a50d53937b9a9e5fc0f06d0d6d06
```

All producer bytes remain unchanged. The parent owns filing and promotion.
