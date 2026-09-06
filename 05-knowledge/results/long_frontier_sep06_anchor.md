# Residue moments close the smallest original phase uniformly

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
In the two-anchor q=5,h=4 model, the smallest positive zero s of P(-s)
satisfies `1/110<s<1/90` and `sQ(-s)<-160`. Only the C interlacer is needed.
Other positive zeros of the original first row and general actual Laurent
noncancellation remain **OPEN**.

## 1. Inheritance and the remaining target

The incoming [effective-anchor theorem](continuing2_20260906_effective_anchor.md)
and [independent referee](continuing2_20260906_effective_anchor_audit.md)
were read first. They prove e4>1/100 when both prescribed C,D interlace B,
and prove negativity at every original phase s>=118163898523. This already
makes the earlier [qualitative tail](open_frontier_sep06_quadratic_anchor.md)
effective; no part of that tail argument is reproved here.

The canonical hostile is the [linearly anchored escaping family](open_frontier_sep06_laurent.md).
The corrected near miss is promoting a bounded rational shape probe to
all-root negativity. The underused sidecar is the positive residue measure
of C/B, before replacing interlacing by scalar coefficient bounds. The five
live concepts are original-root elimination, the two moment anchors, C/D
residue geometry, the full carried response, and the coupled-window cone.
The source-to-target map below turns one interlacer into explicit coefficient
inequalities, then uses the original zero to cancel its third coefficient.
It forgets the detailed root locations; those are still needed for the other
first-row phases. No fixed-shape scan supplies a quantifier in this theorem.

## 2. Positive residues give a new coefficient envelope

Let a_i>=0, i=0,...,4, and suppose their elementary coefficients satisfy
`e1=13,e2=55`. Thus `sum a_i^2=59`. Define

```text
B(v)=product(v-a_i)=v^5-13v^4+55v^3-e3 v^2+e4 v-e5,
C(v)=v^4-12v^3+45v^2-(2/3)e3 v+(3/7)e4.
```

Assume C weakly interlaces B. After cancelling common roots, C/B has simple
poles on the nonnegative axis and positive residues. The sign follows from
the alternating root order; its behavior 1/v at infinity makes the residue
sum one. Equivalently

```text
C(v)/B(v)=sum_i w_i/(v-a_i),   w_i>=0, sum_i w_i=1,
mu_j=sum_i w_i a_i^j.
```

Cancelled nodes may have zero weight. This includes repeated and zero B roots.
Comparing the first five coefficients at infinity gives

```text
mu_0=1, mu_1=1, mu_2=3,
mu_3=e3/3-16,
mu_4=16e3/3-373-(4/7)e4.
```

Cauchy on the functions sqrt(v),v^(3/2) gives `mu_1 mu_3>=mu_2^2`, hence
e3>=75. The Gram matrix of 1,v,v^2 is positive semidefinite, so

```text
0<=det(mu_(i+j))_(i,j=0)^2
  =(e3-75)(135-e3)/9-(8/7)e4.
```

Because e4>=0, these imply the uniform envelope

```text
75<=e3<=135,
0<=e4<=(7/72)(e3-75)(135-e3)<=175/2.                   (1)
```

Also AM-GM gives `0<=e5<=(13/5)^5<119`. The lower e3 bound is attained in
the C-only weak class by `(0,0,3,5,5)`: cancellation leaves
`C/B=(v-2)/(v(v-3))`, with weights 2/3 at0 and 1/3 at3. The upper bounds
are not asserted sharp. The additional D interlacer is unnecessary for (1).

## 3. The smallest phase lies in a fixed short interval

Keep all original Laurent shifts and the induced contiguous coefficients:

```text
G_B=(1,13,55,e3,e4,e5),
G_C=(1,12,45,(2/3)e3,(3/7)e4),
G_D=(1,11,36,(5/12)e3,(1/7)e4),
beta=z^-1 G_B, C_raw=z^-1 G_C, D_raw=z^-1 G_D,
O=sum_j binom(14,2j+1)z^j, E=sum_j binom(14,2j)z^j,
P=O star beta,
Q=(O^2+z^-1 E^2) star(beta^2+2z C_raw D_raw).
```

Tuple notation means ordinary coefficients. No real-rootedness of D is
needed for this theorem, although the original two-interlacer class has it.
The original phase polynomial is

```text
P(-s)=182-20020s+2002e3 s^2-3432e4 s^3+2002e5 s^4.
```

For 0<s<=1/110, its first two terms have nonnegative sum and

```text
e3-(12/7)e4 s+e5 s^2 >= 75-15/11 >0.
```

Thus P(-s)>0 there. At s=1/90, dropping the negative e4 term gives

```text
P(-1/90)<=182-20020/90+2002*135/90^2+2002*119/90^4<0.
```

Moreover the derivative of P(-s) on 0<=s<=1/90 is at most
`-20020+4004*135/90+8008*119/90^3<0`. Consequently there is exactly one
root in this interval; it is simple and is the smallest positive root.
Write `I=[1/110,1/90]`. We next prove Q<0 at every original root in I.

## 4. Original-root elimination proves a strict full-response margin

Write b=e4,f=e5. At an original root, solve for the third coefficient:

```text
e3=10/s-1/(11s^2)+(12/7)b s-f s^2.                    (2)
```

The complete raw Q has exponent range -1,...,8, with coefficient 28 at
exponent -1; upper coefficients may vanish on the weak boundary.
Substitute (2) into sQ(-s) only after forming those complete coefficients.
The exact resulting polynomial is

```text
sQ(-s)=A_2(s)b^2+A_11(s)bf+A_1(s)b
        +A_02(s)f^2+A_01(s)f+A_0(s),
A_2 =-39330s^7(19601s+31920)/49,
A_11=26220s^8(9605s+24708)/7,
A_1 =30s^3(600168371s^3+1898076564s^2-166142340s+1753752)/77,
A_02=-2185s^9(7735s+24708),
A_01=-5s^4(1724814091s^3+5260283632s^2-716499174s+13585572)/11,
A_0=-7H(s)/121,
H(s)=9335990950s^4+19125419885s^3-1420269695s^2
     +19755450s-63866.                                (3)
```

For s in I, A_2,A_02 are negative. The cubic inside A_01 is at least
`13585572-716499174/90>0`, so A_01 is negative as well. Discarding the
negative monomial in A_1 and using s<=1/90 gives A_1<2; direct monotonicity
of its positive factors gives A_11<1/1000.

The quartic H is decreasing on I: H''<0 throughout I follows by bounding
both positive terms of H'' at1/90, and H'(1/110)<0. At the right endpoint,

```text
H(1/90)=4379140411/656100 >6000.
```

Using b<=175/2, f<119 and b,f>=0 in (3), we therefore obtain

```text
sQ(-s)<-42000/121+175+(175/2)*119/1000<-160.            (4)
```

This proves the announced uniform strict negativity at the smallest
original phase. The numerical interval and margin are convenient bounds,
not claimed optimal. The original zero is retained in every step; no root
of a contiguous replacement is substituted.

## 5. What is closed and what remains

The initial phase branch now has a uniform sign certificate throughout a
class larger than the two-interlacer model. The C-only repeated boundary
`(0,0,3,5,5)`, previously useful for exposing the role of D in tail bounds,
is an accepted positive control here. The linearly anchored geometric hostile
fails the second coefficient anchor and does not contradict this theorem.

The attempted alternative of eliminating e5 leaves an indefinite quadratic
form in e3,e4 at every positive phase. Its Hessian determinant is

```text
-1031232600s^12(10966105s^2+72692884s+144097056)/49 <0.
```

Thus unconstrained quadratic-form negativity cannot settle the other roots.
The residue envelope is the surviving extra information. Its extension to
higher moments or to the D residue measure is the next precise test.
This note does not close the other original phase branches between1/90 and
the inherited effective tail, and does not prove general two-rung separation.

## 6. Exact certificate and reproduction

The standalone [producer](../../04-computation/long_frontier_sep06_anchor.py)
and [frozen output](long_frontier_sep06_anchor.out) pass **106 explicit exact
gates**. The universe is the displayed formal polynomial identities and
rational interval inequalities, plus two named root-shape controls and one
excluded prior hostile. It is not a shape census. All arithmetic uses the
standard-library `Fraction` type, and checks remain active under `python3 -O`.

The producer obtains the moments from formal division of C by B, forms the
moment determinant, builds the complete negative-support response from the
odd and even binomial rows, and verifies the full six-term identity (3).
The strict control has roots `(1,3,9,22,30)/5`; the weak C-only boundary has
roots `(0,0,3,5,5)`, with `D(5)=-25/4` explicitly checked. A separate ordinary
carrier-convolution path checks the original first row and full response
for both controls, including the raw exponent -1 contribution. The old
geometric R=512 hostile is checked to fail e2=55. The proof uses exact
identities and uniform inequalities, rather than inferring its quantifiers
from the controls.

Reproduce from the worktree root:

```bash
python3 04-computation/long_frontier_sep06_anchor.py > /tmp/anchor-normal.out
python3 -O 04-computation/long_frontier_sep06_anchor.py > /tmp/anchor-optimized.out
cmp /tmp/anchor-normal.out /tmp/anchor-optimized.out
cmp /tmp/anchor-normal.out 05-knowledge/results/long_frontier_sep06_anchor.out
```

Normal and optimized outputs are byte-identical. SHA-256 manifest:

```text
source  8b1fb8a11a86df1f82a1e54552913b35d2d7bd7b9b4a38cb3c843531dc6db8ff
output  a851496e1fdb59ffa0fd38376c44140a32477b4638b84d0a1eee8c35ad0968ca
semantic 5f16b842b02c65e18f1c8cd341129029ce30ef986f16dc374d8fbdfd3b75a912
```

The [independent analytic/source audit](long_frontier_sep06_anchor_audit.md)
passes. Its separate ordinary-carrier symbolic reconstruction verifies
the entire response and all interval bounds in28 always-active gates;
normal and optimized outputs agree. No literature-priority claim, new
canonical ID, or closure of the remaining phases is made.
