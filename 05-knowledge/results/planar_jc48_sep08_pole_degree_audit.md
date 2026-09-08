# Independent audit of the rational primitive pole-degree gate

**Status: INDEPENDENT ANALYTIC + EXACT SOURCE + REPLAY AUDIT PASS.**
The [primary](planar_jc48_sep08_pole_degree.md) proves the stated
componentwise rational-mate obstruction and the literal quartic's
rational nonintegrability. Its affine critical point already excludes
polynomial mates; the additional rational conclusion has a separate
proof and is not credited to that criticality.

## 1. The exact generic component and its local degree

I read the complete proof and source. Nonconstant F|D guarantees a
generic transverse point on D. Removing finitely many critical values
of F and F|D makes the fibre smooth on W and F_b nonzero at these
points. For a proposed rational G, also removing vertical pole values
ensures its meromorphic restriction is defined on every compact
normalized component. A horizontal pole does not prevent restriction;
the differential equation then forbids its pole on the smooth fibre.
No global regularity hypothesis on ambient G is required.

In the literal chart the form is -r² dr/F_b. Since r is a uniformizer
on that fibre, a primitive has local expansion G(P)+u r³+O(r⁴), u!=0.
It is therefore nonconstant on the component containing P and has
local map degree three. This is an exact order, not a lower bound
inferred from the divisor class of omega.

The differential is regular at all points of the generic fibre in W.
A pole of a meromorphic function has a derivative pole one order
higher. Hence the primitive can have poles only on the projective
boundary, with total degree exactly the displayed sum of (pole order
of eta minus one). Simple differential poles are impossible, and all
residues vanish under exactness. The degree of a nonconstant map to P1
equals its pole divisor degree, which cannot be smaller than one local
fibre multiplicity. This proves the necessary bound three separately
on each component meeting D. Degrees at different D points are not
summed without a shared target value.

## 2. All branches of the residue-free example are included

The octic section has precisely its order-six zero at x=0 and its
order-two zero at infinity. The quartic section is nonzero at both.
At zero, Ebar(s,0) has order three in s. Its initial factor after
s=y²h is -(h-4/9)²(h-1/9). The simple h=1/9 branch is regular for
eta. The two h=4/9 branches split by the displayed nonzero generic
coefficient e. Their simple implicit roots produce two actual power
series in y, and together with the first branch exhaust Weierstrass
degree three. Each lifts to one smooth original branch parametrized
by x, with y=x²; no hidden double cover is counted.

The derivative order five in y gives exactly order-two eta poles in
x. Their coefficients are even Laurent series in x, so their residues
vanish identically. I checked the corrected e² value and the fixed
orientation: the sign differs from the predecessor's form displayed
up to sign, while all pole orders and residues agree.

At infinity the full numerator r² is retained. The two nonzero
implicit roots a²=729/512 exhaust the order-two Weierstrass equation
there. The actual eta orders are three on each branch. Consequently
the complete boundary order list is (-2,-2,0,3,3), with total primitive
pole budget two across all compact components. Every component through
a transverse D point has at most that budget and needs local degree
three, a contradiction. No genus or generic irreducibility claim is
required for this argument.

## 3. Hostiles and the exact strengthened conclusion

The supplied affine point x²=54/97, t=7081/1296 has both derivatives
zero, H=85/36 and F=-77/36. It explains the elementary polynomial
exclusion. The rational pair (x²,t/(2x)) demonstrates why it cannot
explain rational nonintegrability. The primary makes this distinction.

The existing boundary-linear rational mate has generic pole degree
eight and satisfies the new bound. The pair (t,-x) retains the constant
D exception. The rational map u+1/u has the same two residue-free
double poles as the new example, but lacks the order-two differential
zero at D; it confirms that those poles alone do not obstruct exactness.
The u³ control makes the degree-three local bound sharp as a compact
curve principle, without asserting a DG realization of that degree.

The rational strengthening of the safe quartic criterion also passes:
its already proved local calculation makes eta holomorphic on every
compact generic component. Meromorphic restriction of a rational mate
would then be holomorphic and constant, contradicting eta!=0. This
uses no extra ambient regularity and leaves the frozen predecessor
unchanged. The general degree bound and zero residues are necessary
conditions, never an existence theorem for a primitive.

## 4. Reproduction and pins

Normal, optimized, and frozen outputs agree on all634 bytes and68
always-active gates. The source directly verifies both full charts,
every finite and infinite branch derivative, the implicit-root tests,
the complete boundary order list, criticality, and positive/hostile
rational controls. The analytic argument supplies the unbounded
degree and componentwise conclusions.

    python3 -B 04-computation/planar_jc48_sep08_pole_degree.py
    python3 -B -O 04-computation/planar_jc48_sep08_pole_degree.py

Source SHA256: `58743ec60c599028800c560482dae7a7928bd7b97a9818e531e94d47c5bab151`.

Output SHA256: `3ffdc88c37dac410055cff384b84c7c516dbd044201656ecc3a4fb9932ef13bf`.

The frozen source/output retain their pre-audit RESERVED production
label. This independent audit and the promoted primary govern status.
JC(2) remains OPEN.
