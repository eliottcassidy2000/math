# Independent referee: no source-linear logarithmic witness for a W2 submersion

**Status: PASS, independent analytic audit and independent exact controls.**
The [source-linear witness theorem](continuing14_20260908_logarithmic_linear_witness.md)
is accepted with its stated all-degree quantifiers. Its degree restriction
is on the polynomial witness `H=A(x)t+B(x)`, not on the first function.
There is no inference about witnesses of t-degree at least two or about
general impossibility of a nonzero unit of order one.

## 1. Independent analytic reconstruction

Write `F=sum f_j(x)t^j`. The literal equation `J(F,H)=F` has coefficient
recurrences

    A f_j' - j A' f_j - (j+1) B' f_(j+1) = f_j.              (1)

If `A=0`, its largest nonzero coefficient is impossible. If `A!=0`, use
the rational polynomial coordinate `v=H` over `C(x)`. The equation is
`A (partial_x F at fixed v)=F`; its nonzero coefficient solutions differ
by constants. Hence `F=Y(x)P(v)`, with `Y'/Y=1/A` and `Y` rational.
This change of coordinate is used only over `C(x)`, and the original
vertical divisors are restored before drawing a source conclusion.

A rational logarithmic derivative has only simple poles with nonzero
integer residues and no polynomial part. Therefore `A` is nonconstant
and squarefree. Its root residues are exactly the orders of `Y`.
Every root of `P` is simple, since a repeated root would give a repeated
component of `F=0` on the nonempty locus `A!=0`, contradicting submersion.
At a root alpha of A the order of `H-B(alpha)` along the ORIGINAL vertical
divisor is exactly one, with first coefficient `A'(alpha)t+B'(alpha)`.
Thus the order contributed by `P(H)` is either zero or one.

Polynomiality and source submersion now force each residue to be `+1`
or `-1`. A positive residue requires `P(B(alpha))!=0`; a negative one
requires a simple root there to cancel its simple vertical pole. Larger
positive orders make a repeated zero divisor, and larger negative orders
cannot be canceled by the squarefree P. These implications are necessary,
not merely sufficient tests on a selected family.

For `N=deg A>=2`, expansion at infinity annihilates the signed root
power sums through degree `N-2`. The zeroth sum makes `N=2n`, with n
roots of each sign. If n>=2, these power sums include degrees one through
n. Newton identities identify the two monic root polynomials, contradicting
disjointness of the signs. Thus `N=2`. For `N=1`, the single residue is
again one of the two signs. This independently pays the unbounded-degree
step; a finite search over A is not used.

## 2. The exact three rows and their boundaries

The only possible rational factors Y, up to a constant absorbed into P,
are `x-alpha`, `1/(x-beta)`, and `(x-alpha)/(x-beta)` with alpha!=beta.
They give exactly the producer's positive, negative and mixed rows.
The conditions on `P(B(alpha))` and `P(B(beta))` ensure polynomiality.
They also prove source submersion: a critical point with `F!=0` is already
excluded by `J(F,H)=F`, while every zero component is simple and disjoint
from the other zero components, with a nonzero displayed transverse
derivative. Constant-P and constant-R cases are retained, including units
that vanish rather than have order one.

In the actual second chart `x=1/r,t=-r^2-r^4 b`, the positive row always
has a pole. When B is nonconstant its leading polynomial part dominates;
when B is constant the source condition makes the leading `x P(B)`
nonzero. In the mixed row the two required values of B differ, so B is
nonconstant; `Y` tends to one and `A t` to a constant, leaving the
nonzero polynomial leading term of `P(B(x))` as a pole.

For the negative row use the polynomial expression

    F=(-t+C(x)) R(b0+(x-beta)(-t+C(x))).                    (2)

If C has positive degree k and R degree M, its leading pole degree is
`k+M(k+1)>0`. If C is a nonzero constant and M>=1, the leading pole
degree is M. The only globally regular cases are C=0 with any permitted
R, and C constant with R constant. For C=0, the second-chart expansion
starts at `r^2 R(b0)` with `R(b0)!=0`; for constant R it is a nonzero
multiple of `C+r^2+r^4 b`. Both derivatives vanish along every point of
the added divisor. Thus none of the complete source-submersive rows
extends to an everywhere submersive function on W2.

## 3. The order-one hostile and the sign repair

For `F=t(1+xt), H=-xt`, the rational mate is `G=H/F=-x/(1+xt)`.
The two source zero components have H-values zero and PLUS ONE.
The labelled primitive principal parts are therefore `(0,+1/F)`.
A preliminary prose version wrote the second sign incorrectly; the
producer's exact labelled-component gate and this audit both require +1.
The final frozen primary pins below contain the repair. The obstruction
to simultaneous cancellation is unchanged by that repaired sign.

The explicit inverse `t=F/(1-H), x=-H(1-H)/F` gives the rational field
`C(x,t)=C(F)(H)` and the rational constants `C(F)`. Every rational mate
with Jacobian one differs by a rational function of F. A scalar pole
cannot cancel the nonzero part on one zero component without introducing
it on the other. The unit is consequently nonzero of exact order one.
The function F is globally regular, but its entire added divisor is
critical, as required by the classification.

The coordinate family `F=-t+C(x)` instead has polynomial mate x and
zero unit. The degree-three A hostile with residues `(+1,-2,+1)` has
an actual repeated factor and fails source submersion. Together these
controls prevent both a false converse and an overly strong residue claim.

## 4. Independent finite universe and accepted bytes

The [independent engine](../../04-computation/continuing14_20260908_logarithmic_linear_audit.py),
[output](continuing14_20260908_logarithmic_linear_audit.out), and
[certificate](continuing14_20260908_logarithmic_linear_audit_certificate.json)
use forty explicitly declared normal-form pairs: twelve positive, twelve
negative and sixteen mixed. Each is tested by the separate coefficient
recurrence (1), its complete affine critical ideal via Groebner reduction,
and its full actual-chart Laurent expansion. Exactly six are globally
regular, and each has the entire added divisor critical. Separate controls
check Newton reconstruction, the corrected signed primitive, the zero-unit
coordinate case and the higher-degree residue hostile. No producer engine
is imported. There are 234 always-active exact gates per mode.

Frozen primary packet accepted by this referee:

    report c1f140294186bf6384b538b0bf929a3005d02fd3ab1b3524d8450c4f44caa27c
    source f2712ea491ea0a12b3ec148f8d4623ee460770b09cd2fd40da7e2a3e187141bd
    output de52d7a4edf315162d1867ccf49d58f8e4d7ed48a03ed206acd227e107fa375f
    cert   6fca12419f46fcaa2311e3373288afef01a0a6b4f4b6fd37e0000c778eafec49

The producer reports 467 always-active gates. Those and the independent
234 gates supplement the separate analytic proofs. Reproduce the audit
in the repository with normal and optimized Python on the linked source;
raw LF stdout and regenerated certificate must match in both modes.

Accepted consequence: any polynomial repair witness for a nonzero
order-one source unit of a global W2 submersion must have genuine t-degree
at least two. This is a necessary witness-degree barrier, not existence,
not a general order-one exclusion, and not a planar Jacobian theorem.
