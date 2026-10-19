# The natural elliptic translation route meets an exact odd-divisor obstruction

**2026-10-19. PROVED elementary square-class obstruction; a separate finiteness corollary uses Faltings.** This investigates a mechanism outside the finite rational-map catalogue excluded by the preceding theorem. It does not exclude every elliptic contact construction.

## 1. Why a fixed elliptic shift would matter

By the standard elliptic canonical-height estimate, a non-torsion arithmetic progression Q_0+nP on one fixed elliptic curve has rational-coordinate height O(n²). If a rational parameter t(Q) converted its consecutive points into distinct compatible Pythagorean petals, a common integer scaling of the first k petals could have O(k³) label bits: each primitive hypotenuse would have O(k²) bits and the logarithm of their lcm is at most the sum. This is a proposed route around the proved exponential obstruction for finite global rational catalogues on P¹.

Existence, distinctness, positivity, and compatible ports are substantial additional conditions. In particular the simplest normalized x-coordinate construction fails. The obstruction below retains the finite square classes that can become squares after an isogeny; it is stronger than checking only an obvious chord formula.

## 2. The Pythagorean spread is reciprocal elliptic duplication

Let d be a nonzero rational number and

    E_d: Y²=X³-d²X,
    f(t)=4t(t²-1)/(t²+1)².

For a point Q where the formulas are defined, the duplication identity is

    X(2Q)=(X(Q)²+d²)²/[4X(Q)(X(Q)²-d²)].

Hence

    f(X(Q)/d)=d/X(2Q).                                  (1)

Let T=(0,0), a rational two-torsion point. Addition by T gives X(R+T)=-d²/X(R), so equivalently

    f(X(Q)/d)=-X(2Q+T)/d.                               (2)

This is an exact connection between the rational Pythagorean spread and elliptic duplication. It also explains why one must consider isogeny pullbacks when studying square contacts: functions that are nonsquares downstairs can sometimes become squares after duplication.

## 3. A general lemma for affine combinations of translated x-coordinates

Work over a characteristic-zero field, and allow algebraic closure when discussing zeros and their multiplicities. Let

    E: y²=x³+a x+b

be smooth, and let P=(p,q) with q≠0. Thus P is neither the identity nor a point of order two. Put

    x_0=x(Q), x_1=x(Q+P),
    g=(y(Q)-q)/(x(Q)-p), w=x_0-x_1.

The chord formula gives

    x_0+x_1=g²-p,
    w²=Delta(g)=g⁴-6p g²+8qg-3p²-4a.                    (3)

For example, substitute y=q+g(x-p) into the curve equation and divide out the root x=p. The remaining quadratic has sum g²-p and discriminant Delta(g), proving both formulas.

The function g has simple poles at O and -P, so it defines a degree-two map E→P¹. Infinity is unramified, because those two poles are distinct and simple. Its four finite branch values are therefore the four distinct roots of Delta. The involution exchanges Q and -P-Q, and sends w to -w.

### Lemma 1: equal coefficients

For constants A,B with B≠0, the divisor of

    F_+=A+B(x_0+x_1)

is even if and only if A=Bp. In that case F_+=B g².

**Proof.** By (3), F_+=B(g²-r), with r=p-A/B. If r=0 the divisor is even. Otherwise its two zero values on the g-line are the distinct opposite numbers s,-s with s²=r. To have an even zero divisor on E, each must be a branch value of the degree-two g-map. But

    Delta(s)-Delta(-s)=16qs,

which is nonzero because q,s≠0. Thus they cannot both be branch values. At least one unramified preimage is a simple zero. ∎

### Lemma 2: arbitrary unequal nonzero coefficients

For constants A,B,C with BC≠0 and B≠C, the divisor of

    F=A+B x_0+C x_1

is never even. In particular this includes opposite coefficients.

**Proof.** Put

    beta=(B+C)/2, gamma=(B-C)/2≠0, alpha=A-beta p.

Equation (3) gives F=alpha+beta g²+gamma w. At a branch point of the degree-two g-map, w is a local parameter of order one while g-g_0 has order two. If F vanishes there, its gamma w term makes that zero simple. Thus an even zero divisor could have no zeros at branch points.

At an unramified zero of F, its conjugate under w↦-w is nonzero: simultaneous vanishing would force w=0. Consequently the multiplicity of that zero equals the multiplicity of the corresponding root of the norm polynomial

    N(g)=(alpha+beta g²)²-gamma² Delta(g).

Its leading coefficient is beta²-gamma²=BC≠0, its cubic coefficient is zero, and its linear coefficient is -8gamma²q≠0. If every zero of F were even, this quartic would be a constant times a square quadratic over the algebraic closure. Its zero cubic coefficient forces that quadratic's linear coefficient to vanish, which also forces the square's linear coefficient to vanish. This contradicts -8gamma²q≠0. ∎

Combining the two lemmas gives the exact classification: for BC≠0, the divisor of A+B x_0+C x_1 is even if and only if B=C and A=Bp. The exceptional function is B g². Over Q, an actual square identity additionally requires B to be a rational square. Passing to a geometrically connected isogeny cover defined over Q does not change that constant-field requirement.

The poles in both lemmas are O and -P, each of order two, because the two x-coefficients are nonzero. Thus the exhibited odd zero multiplicity is exactly the obstruction to an even principal divisor; it is not merely a failure of one attempted square-root formula.

### Corollary 3: the obstruction survives every isogeny pullback

A characteristic-zero isogeny of smooth elliptic curves is étale. Pulling back a rational function therefore preserves the multiplicity at each point above a zero or pole. A function with an odd zero cannot become a square after any isogeny pullback, nor after composing that pullback with a translation.

In particular the exceptional possibility that an even-divisor nonsquare is killed by a two-isogeny does not rescue the functions excluded by Lemmas 1–2. The same result also rules out an attempted two-function repair using unequal affine coefficients of translated x-coordinates. Arbitrary higher-degree parameter functions remain outside this classification.

## 4. No normalized-x contact identity for a non-torsion shift

Let P∈E_d(Q) be non-torsion. For any signs epsilon,eta consider the function

    H(Q)=1+epsilon f(X(Q)/d)+eta f(X(Q+P)/d).            (4)

Set R=2Q+T and A=2P. Using (2), H is the pullback of

    1-(epsilon/d)X(R)-(eta/d)X(R+A).                    (5)

The point A is non-torsion, so its y-coordinate is nonzero.

If epsilon=-eta, (5) has opposite coefficients and Lemma 2 gives an odd zero. If epsilon=eta, Lemma 1 says its divisor could be even only when

    1=-(epsilon/d)X(A), or X(A)=-epsilon d.

On E_d, the points with x-coordinate ±d have y=0 and are two-torsion, contradicting non-torsion of A. Thus (5) has an odd zero in every sign case. Corollary 3 transfers that obstruction through Q→2Q+T.

**Theorem 4.** For every non-torsion rational shift P on E_d and every sign pair, the normalized-x contact function (4) has an odd divisor. In particular it is not a square in Q(E_d), even after any additional étale isogeny pullback.

This rules out the natural fixed-elliptic-shift identity, including its finite two-torsion square-class escape. It does not rule out a different rational parameter t:E→P¹ or a different invariant genus-one pencil on the contact surface.

## 5. Pointwise finiteness and a general necessary condition

A separate deep input strengthens the conclusion. For a nonconstant rational function H on an elliptic curve with at least one odd zero or pole, the smooth projective normalization of Z²=H is a geometrically connected double cover branched at a positive even number r of points. Riemann–Hurwitz gives

    genus = 1+r/2 >= 2.

Faltings' theorem then gives finitely many rational points on that cover. Hence only finitely many rational base points Q, away from the finite exceptional zero/pole set, can have H(Q) a rational square.

Applied to (4), this proves: for fixed E_d and fixed non-torsion P, only finitely many rational Q give a square contact between the petals represented by X(Q)/d and X(Q+P)/d. The union over all four sign pairs is still finite. Consequently this normalized-x arithmetic progression cannot produce an infinite rational contact chain even by choosing signs pointwise. In fact the degree-four duplication pullback has four points above each original odd-divisor point, so its double-cover genus is at least five. No effective bound or enumeration follows from this use of Faltings.

More generally, any proposed infinite rational fixed-shift construction on an elliptic curve must make its contact function have an **even divisor**. A square identity is sufficient, but a nonsquare with even divisor may instead define an unramified genus-one double cover. That finite square-class possibility remains relevant for other parameter maps. It is precisely what the elementary proof excludes for the normalized-x candidate.

This supplies a concrete criterion for the next search: seek an invariant genus-one correspondence with even contact divisors, then separately prove rational liftability, distinctness, positivity, and compatible port choices. Merely replacing a rational map by an elliptic orbit does not ensure those properties.

## 6. Status and references

The duplication, chord-discriminant, and odd-divisor arguments are elementary identities and divisor arguments proved above. The rational-point finiteness in Section 5 is explicitly an application of Faltings, not a new proof of that theorem. Primary reference: G. Faltings, “Endlichkeitssätze für abelsche Varietäten über Zahlkörpern,” Inventiones Mathematicae 73 (1983), 349–366, Satz 7 on p. 365. [Publisher DOI](https://doi.org/10.1007/BF01388432); [original paper scan](https://math.uchicago.edu/~drinfeld/Deligne%27s_conjecture_Manin_conf/Faltings_argument/Faltings.pdf).

From the repository root run `python3 05-knowledge/experiments/verify_elliptic_chord_obstruction_20261019.py`. The standard-library verifier checks the general chord quadratic and discriminant, the opposite-branch obstruction, the general norm leading/cubic/linear coefficients for arbitrary nonzero translated-x coefficients, the square-quartic obstruction, duplication and two-torsion translation, and independent normalized norm quartics for all four sign choices. Its published receipt is `05-knowledge/results/elliptic_chord_obstruction_controls.json`.
