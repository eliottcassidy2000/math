# A signed polynomial dual for one injection atom

2026-10-05. **PROVED:** the exact signed minorant, uniform geometric error,
monotone recovery, coefficient-conditioning formula, and conditional
moment-to-atom certificate. **FINITE-EXACT:** the declared controls below.
**OPEN:** independent moment lower bounds that certify presently unresolved
Collatz sources, and universal positivity. The program checks arithmetic
consequences of supplied moment intervals; it does not certify arbitrary
intervals as truthful data about the Collatz measure.

[Program](../../04-computation/experiments/collatz_atom_polynomial_dual_20261005.py)
and [output](collatz_atom_polynomial_dual_20261005.out).

## 1. Inheritance and the new source coordinate

Three nearby proved mechanisms guide this construction:

* [THM-2237, truncated Boolean moments and the top atom](../../01-canon/theorems/THM-2237-truncated-boolean-moment-interval-and-parity-top-atom-majorants.md)
  gives exact atom intervals and pointwise polynomial duals on a finite
  integer support. The transferable ingredient is the inequality direction
  of a minorant, not its particular Boolean support.
* [THM-2842, ordered positive-cone multiplier observability](../../01-canon/theorems/THM-2842-ordered-positive-cone-vandermonde-multiplier-observability.md)
  makes cardinal selectors and their external readouts explicit. Its finite
  quotient and unit promise are not imported. Here the input is a specified
  infinite discrete support, and the multiplier degree is allowed to grow.
* [THM-2117, scalar Fourier insufficiency and full Toeplitz separation](../../01-canon/theorems/THM-2117-scalar-fourier-insufficiency-and-full-toeplitz-separation.md)
  distinguishes a few scalar positivity tests from the actual separating
  functional. It motivates retaining the entire tested polynomial; it does
  not by itself give a Collatz atom bound.

The current injection measure is inherited from
[adaptive mixture flow, P5](collatz_adaptive_mixture_flow_20261005.md):

    p_j=lambda(6j+3),  j>=0;    p_j>=0;    sum p_j=1.

Its generating function F(z)=sum p_j z^j and the interpretation of positive
atoms are in [Fourier atom positivity](collatz_fourier_atom_positivity_20261005.md).
A positive p_j is equivalent to a finite ROOT route for the particular
integer6j+3. That equivalence is inherited, not assumed for every integer.

Ordinary positive size moments of lambda diverge, as proved in
[inverse predecessor sections, section 4](inverse_predecessor_sections_20261005.md).
The compact source encoding below avoids those divergent observables:

    x_j=2^(-j),
    eta=sum_(j>=0) p_j delta_(x_j),
    M_l=int x^l d eta(x)=F(2^(-l)),  l>=0.              (1)

In particular M_0=F(1)=1. Every M_l is finite and belongs to[0,1].
This is a Hausdorff moment problem for the **source-index distribution**.
Its variable x is not the price r in the incoming
[Pascal boundary note](collatz_pascal_boundary_leaf_section_20261005.md).
Changing the price prior there does not produce the source moment bounds
required here.

The board is **selected source / compact embedding / signed minorant /
moment oracle / precision cost / ROOT extraction**. The anchor is a lower
bound on p_m. The niche is an explicit infinite-support moment dual. The
wildcard is composing that floor with a finite weight-threshold compiler.
The hostile is an infinite-support measure with one missing atom despite
strict positivity of every finite Hankel matrix. The least-used sidecar is
the signed coefficient error budget.

## 2. An explicit minorant with an all-tail error bound

Fix a target m>=0 and a damping degree d>=0. Define

    q_(m,d)(x) = (x/x_m)^d
                 (x-x_(m+1))/(x_m-x_(m+1))
                 product_(0<=j<m) (x-x_j)/(x_m-x_j).     (2)

The empty product is1. This is a rational polynomial of degree m+d+1.
It has exactly the pointwise properties needed for a lower dual:

    q_(m,d)(x_m)=1;
    q_(m,d)(x_j)=0        for j<m or j=m+1;
    q_(m,d)(x_j)<=0       for j>=m+2.                  (3)

To prove the last line, all the factors indexed j<m have negative numerator
and denominator at a point below x_(m+1), so are positive. The factor with
root x_(m+1) is nonpositive, and the damping factor is nonnegative. Hence
q_(m,d)<=1_(target m) on the complete support, with no omitted tail nodes.

Put

    C_m=product_(l=1)^m (1-2^(-l))^(-1),  C_0=1.

For x=x_j, j>=m+2, one has 0<x<=x_m/4. The damping is at most4^(-d).
The absolute value of the x_(m+1) factor is at most1. Each head factor is
bounded above by x_i/(x_i-x_m). Therefore

    -C_m4^(-d) <= q_(m,d)(x_j) <=0.                  (4)

The constant is uniformly bounded:

    C_m <=32/9.                                      (5)

For m=0,1 this is immediate. The first two reciprocal factors give8/3.
For the remaining factors, induction on the finite product gives

    product_(l=3)^m (1-2^(-l))
       >=1-sum_(l=3)^m 2^(-l) >=3/4.

Combining the two estimates proves(5). No infinite-product theorem or
numerical limit is needed.

**Atom-dual theorem.** Let

    A_(m,d)=int q_(m,d) d eta.

Then, for every probability on the specified support,

    A_(m,d) <= p_m <= A_(m,d)+C_m4^(-d).              (6)

This follows by integrating(3)–(4); the total off-target mass is at most1.
Moreover A_(m,d) increases to p_m as d increases: on every tail node the
negative value is multiplied by x_j/x_m<=1/4, and all cardinal values
are unchanged. Thus

    p_m>0  iff A_(m,d)>0 for some finite d.           (7)

The implication does not require positivity at other atoms. The result
is not an assertion that every p_m is positive.

For m=0 the selector is particularly simple:

    q_(0,d)(x)=x^d(2x-1),
    A_(0,d)=2M_(d+1)-M_d.                            (8)

The general formula eliminates the finite head above the chosen source
coordinate and keeps the rest of the infinite support on the negative side.

## 3. Exact moment input and the conditioning bill

Expand q_(m,d)(x)=sum_(l=0)^D c_l x^l, D=m+d+1. Formula(1) gives

    A_(m,d)=sum_l c_l M_l=sum_l c_l F(2^(-l)).        (9)

Suppose independently certified rational intervals satisfy

    M_l in [a_l,b_l],  0<=l<=D,  a_0=b_0=1.

They must all concern the same actual injection measure. Define

    L=sum_(c_l>=0)c_l a_l + sum_(c_l<0)c_l b_l,
    U=sum_(c_l>=0)c_l b_l + sum_(c_l<0)c_l a_l.

Then the rigorous certificate is

    L <= p_m <= U+C_m4^(-d).                         (10)

In particular L>0 certifies the selected atom. An approximate positive
central value without the signed interval calculation is insufficient.
The bounds may be intersected with[0,1], but the implementation preserves
the unclipped values so the raw dual remains inspectable.

The coefficient norm has an exact product formula:

    B_(m,d)=sum_l |c_l|
      =2^(md)(2^(m+1)+1)
         product_(j=0)^(m-1) (2^j+1)/(1-2^(j-m)).     (11)

All nonzero polynomial roots are positive; its coefficients alternate
in sign, apart from initial zeros. Therefore sum|c_l|=|q_(m,d)(-1)|,
which evaluates to(11). If each moment has absolute error at most epsilon,
the readout error is at most epsilon B_(m,d). The exact product can grow
very large with m and d. Uniform small *tail error* in(6) is therefore not
uniform cheap *numerical access* to the moments.

The exported `moment_interval_bounds` validates exact input types, interval
order, elementary probability bounds, and normalization. It does not check
that an arbitrary interval packet is realizable, or that it describes
lambda. Truthful joint moment bounds are a proof premise. No inferred ROOT
certificate is returned merely because an unsupported packet yields L>0.

## 4. Two hostile controls and the strongest surviving conclusion

**Hankel positivity points in the wrong direction by itself.** For any
polynomial P with P(x_m)=1, nonnegative weights give

    p_m <= int P(x)^2 d eta(x).

Optimizing such squares gives atom upper bounds, not the lower certificate
needed here. Strictly positive moment matrices do not reverse that inequality.
For example, start with p_j=2^(-j-1), remove atom1, and renormalize. Every
finite Hankel matrix is still strictly positive definite: a nonzero polynomial
has only finitely many roots, while the measure retains infinitely many
distinct positive support points. Nevertheless p_1=0. Every valid signed
readout A_(1,d) is nonpositive, as required by(6).

For the original geometric law the moments are explicitly

    M_l=2^l/(2^(l+1)-1).

Removing atom m changes these to

    [M_l-2^(-m-1)2^(-ml)]/[1-2^(-m-1)].              (12)

These exact formulas provide an independent oracle for positive and missing
target tests. They do not use the selector implementation to define the law.

**The same finite bank cannot invent an omitted atom.** Suppose the only
measure information consists of finitely many retained atom values and total
mass1. For any unretained target, another probability measure keeps all those
values and places the remaining mass at a different unretained node. Its
target atom is zero. Every lower bound logically valid for that information
must consequently be nonpositive. Adding moment intervals deduced solely
from that same information does not remove this completion.

This last statement is relative to those data. Independent functional
equations, certified new moment estimates, or other arithmetic restrictions
can eliminate the completion. The dual theorem specifies exactly how a
successful new estimate would be converted into an atom floor. It does not
claim to have obtained that independent estimate for every Collatz source.

## 5. Conditional composition with a finite ROOT compiler

The [weight-threshold receipt compiler](collatz_weight_threshold_receipts_20261005.md)
provides `threshold_receipts(epsilon)`, the exact finite set of sources whose
W-weight is at least a supplied rational epsilon>0, with strict actual ROOT
words. Its proof and input guards belong to that companion package.

Consequently an independently proved moment packet with L>0 in(10) composes
as follows:

1. Retain the exact target n=6m+3 and the moment premises.
2. Set epsilon=L. Since lambda(n)=W(n), the target belongs to the proved
   W-superlevel set.
3. Run the finite threshold compiler and retrieve its actual first-hit word
   at n. Replay that word against n.

If the target is absent, the claimed positive moment lower bound was not a
valid premise about lambda; absence does not refute Collatz convergence.
This is a finite certificate extractor once a quantitative lower bound is
proved. It is not a source of the lower bound itself.

There is also a direct source-specific deadline, proved in
[floor transport and deadlines](collatz_floor_transport_deadlines_20261005.md).
Write epsilon=p/q in lowest positive terms, put M=floor(2q/p), and define

    B=(isqrt(1+4M)-3)//2.

For a nonroot source with W>=epsilon, its hidden counter sum N=L+K obeys
N<=B, and its actual odd ROOT rank is at most B: every nonroot receipt has
K>=1, so its rank L+1 is at most L+K=N. The counter bound follows from
W<=2/[(N+1)(N+2)]. Once the moment floor is proved,
one may therefore replay the selected source for at most B odd steps,
rather than enumerate every other superlevel source. Both procedures depend
on the same independent positive premise; neither creates it.

The finite tests here use already completed injection sources3,9,15,21,27.
From a separately replayed head of sixteen atoms and probability normalization,
they obtain sufficient moment intervals and strictly positive dual floors at
d=1,3,1,1,11 respectively. The exact p_4=lambda(27) is1/11274451650.
These are controls of the representation and error calculation; they add no
new ROOT source. The same head gives no positive floor at its omitted next
index16 in any of the thirteen tested damping degrees.

## 6. Reproduction and next obligation

```text
python -B 04-computation/experiments/collatz_atom_polynomial_dual_20261005.py
python -B -O 04-computation/experiments/collatz_atom_polynomial_dual_20261005.py
```

The exact universe is m,d=0,...,12 and nodes j=0,...,64; closed-form geometric
and missing-atom moment oracles at every required order; five positive Hankel
determinants for the missing-atom law; sixteen actual injection sources3,9,...,93
with literal ROOT cap512; and eleven malformed exact-input controls. The
all-tail signs, product bound, convergence rate, and conditioning formula
are proved above rather than inferred from this finite universe.

The next admissible research input is a rigorous bound for the signed
combination(9), or its individual moment values, whose proof supplies
information beyond the same finite set of already rooted atoms. The
construction gives a direct target for such a bound and a finite receipt
extraction route if it becomes positive. Universal positivity remains open.
