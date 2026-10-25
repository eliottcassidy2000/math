# The genus-five branch curve has five elliptic quotient factors

**Status: PROVED, independently audited** in [AUDIT_BRANCH_GEOMETRY](AUDIT_BRANCH_GEOMETRY.md). This is a specialization of the classical character decomposition for an elementary abelian cover of P¹, not a historical-novelty claim. It supplements §2a of `BRANCH_CUBE_AND_DOUBLE_COVER.md`.

Let C/Q be the smooth genus-five curve with function field

    Q(t)(sqrt(t),sqrt(1+t),sqrt(1−t),sqrt(2−t)).

The four square classes are independent by valuations at the four finite zeros, so its deck group is (Z/2)^4. Its five branch places are {−1,0,1,2,∞}. A nontrivial character determines a double cover of P¹, branched at the places with odd valuation of its product of radicands. Exactly five characters have four branch places: the four products of three radicands, and the product of all four. Each such quotient is a genus-one curve with a rational point at infinity or a rational branch point. The remaining ten nontrivial characters have two branch places and genus zero.

The five quotient equations, with their actual rational coefficients retained, are:

| Omitted branch place | Quotient equation | Rational isomorphism class |
|---|---|---|
| −1 | y²=t(1−t)(2−t) | E₀: y²=x³−x, via x=t−1 |
| 2 | y²=t(1+t)(1−t) | E₀, via x=−t |
| 0 | y²=(1+t)(1−t)(2−t) | E₁: y²=(x+1)(x−1)(x−2), via x=t |
| 1 | y²=t(1+t)(2−t) | E₁, via x=1−t |
| ∞ | y²=t(t+1)(t−1)(t−2) | E₂: Y²=(x−2)(x+2)(x−1), via x=2/t, Y=2y/t² |

The birational transformations between smooth projective genus-one curves extend to isomorphisms. Consequently

    Jac(C) is Q-isogenous to E₀² × E₁² × E₂,

where

    j(E₀)=1728,   j(E₁)=21952/9,   j(E₂)=35152/9.

**Proof of the isogeny, including its direction.** Each quotient map C→E_j induces pullback Jac(E_j)→Jac(C). The sum of these five homomorphisms defines a map from the five-dimensional product to Jac(C). The five spaces of pulled-back regular differentials have distinct nontrivial deck characters, and each is nonzero since the maps are separable in characteristic zero. They are therefore linearly independent; because g(C)=5, they span all regular differentials. Equivalently, the dual norm map from Jac(C) to the product has an isomorphism on cotangent spaces. That map, and hence its dual, is an isogeny. All quotient maps and equations are over Q, so this conclusion is over Q. No claim about the isogeny kernel or integral polarization follows.

The equality of j-invariants alone is not used to identify factors over Q; the displayed coordinate changes preserve the actual twists. The genus-five components of the hourglass branch become copies of C over Q(√2); the decomposition of Jac(C) is already defined over Q.

The cheap hostile checks are the minus sign in t(1+t)(1−t), and the factor4 after the reciprocal change for the quartic. `check_geometry.py` checks both rational identities and computes the three j-values from their root differences. The classical generalized Fermat description and differentials are discussed in Hidalgo, *Holomorphic differentials of Generalized Fermat curves*, §4.3: https://arxiv.org/html/1710.01349 . The short quotient and isogeny proof above is self-contained.
