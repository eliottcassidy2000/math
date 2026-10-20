# Direct exclusion of monic terminal cubic/quintic Keller pairs

20 October 2026. **PROVED; independently audited PASS in the stated scope.**
The full coefficient identities and polynomial-degree argument were
independently reconstructed and checked. This addendum has a separate
[independent audit](review/CUBIC_QUINTIC_AUDIT.md).

## Exact theorem

There are no P,Q in C[x,y] with constant nonzero Jacobian and respective
y-degrees three and five, both having constant nonzero leading coefficients.
In particular, an actual monic Keller source whose generic degree semigroup
is <3,5> is excluded: the separately proved characteristic classification
and inherited constant-leading reduction force precisely this terminal pair.

The direct theorem is stronger than a fixed-semigroup theorem. The proof
uses neither a classification of elliptic polynomial maps nor a Mason
inequality. It is an exact coefficient reduction followed by degrees of
polynomials in x. It makes no claim about unbounded terminal pairs or JC(2).

## 1. Lawful normalization and all positive-degree equations

Constant rescaling makes both coordinates monic. A source shear
t=y+s(x), determined by one third of the cubic quadratic coefficient,
centers the cubic and preserves the constant Jacobian. Write the actual pair

    p=t^3+a t+b,
    q=t^5+c t^4+d t^3+e t^2+f t+g,
    a,b,c,d,e,f,g in C[x].

With prime denoting x-derivative at fixed t, the coefficients of the ENTIRE
Jacobian p_x q_t-p_t q_x, from t^6 down to t^0, are

    -3c',
    5a'-3d',
    -ac'+4ca'+5b'-3e',
    -ad'+3da'+4cb'-3f',
    -ae'+2ea'+3db'-3g',
    -af'+fa'+2eb',
    -ag'+fb'.

Thus c is constant, and d=5a/3+D with D in C. The constant target shear
q -> q-Dp removes D and preserves the nonzero constant Jacobian.
The renamed remaining coefficients satisfy d=5a/3. Integrating the
t^4,t^3,t^2 equations successively gives

    e=4ca/3+5b/3+E,
    f=5a^2/9+4cb/3+F,
    g=2ca^2/9+10ab/9+2Ea/3+G,
    E,F,G in C.                                    (1)

No x-dependent target shear was used. In particular the actual scalar
terms b,g are retained; dropping them would invalidate the equations.

## 2. The exact first integral and residual Jacobian

Put

    U=-5a^2/9+4cb/3+F,
    V=4ca/3+10b/3+2E,
    H=-5a^3/27+4cab/3+Fa+5b^2/3+2Eb.

Substitution of (1) into the whole Jacobian gives the identity

    J(p,q)=(Ua'+Vb')t+Ub'-(aV/3)a'.               (2)

Furthermore H'=Ua'+Vb'. Since J(p,q)=j in C*, the coefficient of t
vanishes, so

    H is constant,    j=Ub'-(aV/3)a'.             (3)

These are polynomial identities in x, not merely formal asymptotic
conditions. The checker validates (1), (2), and H' universally.

## 3. Polynomial degrees force a nonconstant Jacobian or zero

If a is constant and b nonconstant, the term 5b^2/3 is the unique
highest-degree term of H, contradiction. If b is constant and a
nonconstant, the term -5a^3/27 uniquely dominates, contradiction.
Thus either both are constant, giving J=0 in (2), or both are nonconstant.

In the latter case let r=deg a>0 and s=deg b>0. The mixed term ab has
degree r+s, strictly below max(3r,2s). To verify this strict inequality:
if s<2r, then r+s<3r; if s>=2r, then r+s<2s because r<s.
The other terms in H have still smaller degree.
Therefore constant H requires cancellation between its a^3 and b^2 terms:

    3r=2s, hence r=2k, s=3k for an integer k>=1.

Write A,B for their nonzero leading coefficients. The leading terms of U,V
are respectively -(5/9)A^2 x^(4k) and (10/3)B x^(3k).
Consequently (3) has leading term

    -(35k/9) A^2 B x^(7k-1).                     (4)

Indeed Ub' contributes -(15k/9)A^2B and -(aV/3)a' contributes
-(20k/9)A^2B. They reinforce instead of cancel. The coefficient is
nonzero in characteristic zero, and 7k-1>=6. Thus j is nonconstant,
contradicting j in C*. This proves the theorem.

The first-integral leading relation B^2=A^3/9 is compatible with this
calculation, but is not needed except to establish the degree tie.
It does not annihilate A^2B.

## 4. Controls and scope

The checker `verify_cubic_quintic.py` independently expands the raw
Jacobian, all integrations, the full reduced identity, and the leading
coefficient for each symbolic positive integer k through an exact general
coefficient formula.

The family a=u^2, b=u^3/3, c=E=F=G=0 gives H=0 identically.
It is an essential hostile control: H constant does NOT force a,b constant.
For u=x^k, its whole Jacobian is exactly

    -(35k/27) x^(7k-1).

Thus the proof must use the final constant-Jacobian demand, not just the
first integral. These polynomial pairs have monic y-degrees (3,5) and
are primitive over C(x), so they are nonvacuous hostile curves of delta four.

The proof is direct and uses no external theorem beyond elementary
characteristic-zero differentiation. Applying it to an arbitrary original
monic <3,5> source additionally uses the separately audited
constant-leading and characteristic terminal-reduction inputs.

Reproduction:

    python3 verify_cubic_quintic.py > cubic-quintic-results.json
    python3 -O verify_cubic_quintic.py > cubic-quintic-results-optimized.json
    cmp cubic-quintic-results.json cubic-quintic-results-optimized.json
