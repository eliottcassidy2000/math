# Independent audit of the squarefree quadratic-witness global obstruction

**Verdict: ACCEPTED analytically, with independent exact controls.** No
mathematical repair is requested. The audited report is
`continuing15_20260909_squarefree_quadratic_witness.md`, SHA256
`2b78aedb18e49e396c4477f71f6cd861e203f882dc58094e2fe85841962208f3`;
its source SHA256 is
`162360786a05ed8d592f635c61ee00094522055cd2e8ac2f90401ad63b02b875`.
The verifier below imports no primary engine and does not replay it.

The conclusion is exact: no everywhere submersive global W2 first
function F can have a genuinely quadratic polynomial repair witness H
whose generic discriminant E+4Av is squarefree in the original x.
There is no degree bound on F or on the coefficients of H. The repeated
discriminant corridor remains outside this theorem.

## 1. The generic pole reduction is exhaustive

On H=v, y=2At+B and dF/F=dx/y. Squarefree D gives no finite poles.
For degree zero or one the infinity pole is double, impossible for a
logarithm. For degree at least three there are no infinity poles, also
impossible for the logarithmic derivative of a nonconstant rational
function on the proper normalization. In degree two the infinity
residues are +/-1/sqrt(d2(v)); their integral values force d2 to be a
constant reciprocal integer square. Hence deg A<=1 and E has degree
two. The proof does not confuse a squarefree affine equation with its
completion, or positive genus with the only possible obstruction.

## 2. Both complete source-submersive classifications pass

For constant A=a, the displayed coordinates U,V form a polynomial
automorphism: their difference recovers x and their sum then recovers t.
The Jacobian 4as is nonzero and H=h0+UV/(4a). Comparing monomials gives
every polynomial eigenfunction, rather than only a rational ansatz:
after an orientation choice it is U^n P(H), ns=1. Source submersion
forces n=1, P squarefree and P(h0)!=0. These conditions are also
sufficient, including constant P. Thus leading E=1 in this case.

For linear A=a(x-alpha), polynomial division by 2A gives the exact
shear tau=t+L(x), with constant remainder b. The transformed witness is
H=aX tau^2+b tau+cX+d, c!=0. Applying the already proved source-linear
classification to -H in the ORDERED coordinates (tau,X) has the correct
sign. At the positive root rho, -2a rho=1. The complete form is

    F=(tau-rho)/(tau+rho) P(H),
    P(h_minus)=0 simple, P(h_plus)!=0.

The expression is polynomial by the factor H-h_minus=(tau+rho)
[aX(tau-rho)+b]. The conditions force b!=0 and include the reversed
residue choice by relabelling rho. They also force leading E=1. No
polynomial denominator cancellation has been assumed without payment.

## 3. Every original boundary case is retained

In the constant-A form, B of degree at least two makes F have a pole.
For affine B, a nonzero x coefficient in U again gives a pole: if UV
has no pole, the only exceptional case has V proportional to t and
H tends to h0, where P is nonzero. The remaining U=2at+u0 gives a
global function only when P is constant or u0=0. In both cases the
entire added divisor is critical; in the latter case F has exact
normal order two. The controls verify the actual substitution
x=1/R,t=-R^2-R^4B, not a substituted polynomial chart mistaken for W2.

In the linear-A form, a nonconstant L or a constant L outside the two
roots gives a pole. At L=rho the cancelled polynomial expression is
regular with exact normal order two and the whole added divisor
critical. At L=-rho, the ratio has order -2 and P(H) has order +1,
leaving a genuine pole of order one. The nonzero coefficient 2a rho
prevents a hidden cancellation. These cases exhaust the shears.

## 4. Hostiles and scope

The repeated-pencil source submersion F=x(xt+1),
H=xt+lambda F^2 lies outside the hypothesis. Its admissible fixed
double root is not removable before taking dx/y. The genuine minimum
quadratic-witness family F=(x+2t)P(t^2-x^2/4), P squarefree/P0!=0,
likewise shows that affine submersion and exact source unit order one
do not prove globality. Conversely the globally regular hostile
t(1+xt) has its added divisor critical. The primary keeps all three
distinctions correctly.

There is no claim about arbitrary nonsquarefree quadratic pencils,
higher witness degree, the entire response ring on W2, or planar JC.
The original response module and original t-degree are preserved.

## 5. Independent controls

The independent engine uses direct original brackets, cancelled
polynomial formulas, and actual Laurent boundary expressions over a
declared bank of 56 source-submersive pairs. Four separate ideals of
both source partial derivatives provide additional exact checks. It
also verifies the symbolic coordinate automorphism and a forbidden
repeated-factor source example. These finite controls support the
all-coefficient analytic argument; they do not supply its quantifiers.

There are 217 always-active gates. All fourteen globally regular sample
survivors have the entire added divisor critical.

Run `python continuing15_20260909_squarefree_quadratic_audit.py` and
the same command with `python -O`. Output and certificate are frozen
byte-identically with LF newlines. The certificate path relocates from
04-computation to 05-knowledge/results when filed in the repository.
