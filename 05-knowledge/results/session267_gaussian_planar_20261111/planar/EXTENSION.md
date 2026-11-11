# Positive all-row closure in a rational boundary carrier

**PROVED:** finite Newton-support criterion equivalent to infinitely many rigid coefficient debts, for a specified rational generating family. **OPEN:** realization of the actual nonlinear Keller source discrepancy in this carrier. This is a positive extension of PROOF.md, not a claimed reconstruction of the inherited G_m family.

## 1. An all-row test with no row cutoff

For a nonzero polynomial A(t,x), define

w(A)=max{2j−r : [t^r x^j]A≠0}.

Set w(0)=−infinity. Let P,Q∈K[t,x] be coprime, Q(0,x)=1, and expand the rational function F=P/Q uniquely in K[x][[t]]:

F(t,x)=Σ_(r≥0) F_r(x)t^r.

**Theorem.** For any integer B, the following are equivalent:

(a) Every monomial t^r x^j of F satisfies 2j−r≤B.

(b) w(Q)≤0 and w(P)≤B.

Thus an infinite coefficient-cone restriction is decidable by two finite support tests on the reduced rational expression. Reduction to coprime P,Q is necessary.

**Proof of sufficiency.** Q=1+Q_+, where Q_+ is divisible by t. Every monomial of Q_+ has weight≤0. The t-adically convergent geometric series 1/Q=Σ_(k≥0)(−Q_+)^k also has only weights≤0. Multiplication by P proves (a).

**Proof of necessity.** The case P=0 is immediate after reducing Q=1. Suppose P≠0. Put V=w(P), W=w(Q); since Q contains 1, W≥0. Substitute x=u/t² and form

P~=t^V P(t,u/t²),    Q~=t^W Q(t,u/t²).

These are polynomials in t,u, neither divisible by t. They are coprime: after inverting t the substitution is a ring isomorphism, so any common factor could only be t, which is absent. Condition (a) says

S=t^B F(t,u/t²)=t^(B−V+W) P~/Q~ ∈ K[[t,u]].        (1)

To see that the formal substitution is lawful, the coefficient at a given t^a u^j is uniquely the original coefficient at t^(a−B+2j)x^j; all exponents a,j are nonnegative under (a). Multiplication by the finite denominator verifies (1) in the formal ring.

We use the standard local divisibility fact

K(t,u) ∩ K[[t,u]] = K[t,u]_(t,u).                  (2)

For completeness, if a/b lies in the formal ring, then a∈b K[[t,u]]. The power-series completion is faithfully flat over the polynomial local ring A=K[t,u]_(t,u), so b K[[t,u]]∩A=bA; hence a/b∈A. Equivalently, contraction of a principal ideal from the completion is unchanged. Thus a reduced denominator of a rational formal series must be a unit at the origin.

If B−V+W<0, the reduced denominator in (1) contains t, impossible by (2). Otherwise Q~ itself must be a unit. But Q~(0,0)≠0 exactly when W=0: the sole original constant-in-t term of Q is 1, which becomes t^W, and a monomial with positive x exponent cannot contribute a constant in u. Therefore W=0 and B−V≥0, proving (b). QED.

The proof works over any field. The local-algebra step concerns rational-function regularity, not Gaussian positivity or a dimension count.

## 2. Application to the inherited row debts

For m≥5 the inherited terminal degree is d_m=floor((m−1)/2)+2. Vanishing of every rigid bracket debt is exactly

deg F_m≤d_m+1=floor((m+5)/2),

or 2 deg F_m−m≤5. Let the tail generating function

R(t,x)=Σ_(m≥5)F_m(x)t^m

be rational, with reduced expression P/Q normalized by Q(0,x)=1. Then **every rigid bracket debt at every row m≥5 vanishes if and only if**

w(Q)≤0 and w(P)≤5.                                (3)

There is no finite-row extrapolation in (3). The hypotheses require the complete tail generating function, not a rational approximant to a finite prefix. If one starts with a rational full series, subtracting its first five terms gives this rational tail before applying the theorem.

A useful positive family is R=t^5 A(t,t²x)/Q(t,t²x), where A,Q are polynomials, Q(0,0)=1. This has w(R)≤−5 and automatically satisfies (3), provided the expression is interpreted in its t-adic expansion. More general allowed numerators may have weight up to 5.

Hostile: F=1/(1−tx) has F_m=x^m. Its denominator contains a weight-one monomial, so no finite B can bound all coefficient weights. Any fixed number of passing early rows can be produced by delaying that pole, e.g. 1/(1−t^N x^N) for N larger than the checked range. These are rational generating families, not asserted actual Keller sources.

## 3. Recursive closure and honest complexity

Let A_B consist of rational functions admitting reduced denominators normalized at t=0 with w(Q)≤0 and numerators of weight≤B. Then

A_B + A_C ⊆ A_max(B,C),    A_B A_C ⊆ A_(B+C).

The proof follows either the support theorem or common denominators. Cancellation cannot increase the resulting support budget. Differentiation gives

∂_x A_B ⊆ A_(B−2),     ∂_t A_B ⊆ A_(B+1).

These are exact graded closure rules. They let a recursive expression tree propagate an integer budget and a finite rational expression. Inversion is not unrestricted: a numerator can introduce a forbidden pole. A rational factor with numerator and denominator both normalized at t=0 and both of weight≤0 has an inverse in A_0.

This is not a bounded-size finite-state theorem: the rational numerator and denominator can grow. It is a finite **certificate per expression**, with a scalar weight budget, and it is consistent with the lower bound in PROOF.md §4. The source family must stay inside the lawful operation carrier; arbitrary multiplication can spend the budget, and inversion is the crucial exit gate.

The precise next research conjecture is therefore not 'all debts disappear by Gaussianity'. It is: actual normalized source discrepancy formation admits a rational-expression carrier whose allowed operations preserve A_5 (after the necessary tail subtraction). Proving that assertion would discharge all rigid bracket debts at once by (3), leaving scalar payment, inhomogeneous depth compatibility, chart entry, and polynomial termination separately open.

## 4. Row-shift intertwiners and their depth obstruction

There are genuine cross-row identities behind the Gaussian bridge. Write a=x²+6. Then

L_(m+1)(a θ) = (m/(m+1)) a L_m θ,                 (4)
Λ_(m+1)(a F) = [12(m+1)/(2m+1)] Λ_m(F),          (5)

where (5) applies for deg F≤2m. Equation (4) follows by differentiation; (5) follows either from the Student density or the Gaussian radial factor. Under homogenization, multiplying by a means multiplying by 6(G²+H²), which commutes with rotation.

The Gaussian Laplacian also commutes with rotation. Define

A_m F=(x²+6)F''−2(2m−1)xF'+2m(2m−1)F.

Then ΔT_mF=T_(m−1)(A_mF) for deg F≤2m, and for m≥2,

A_m L_m = ((m−1)/m) L_(m−1) A_m.                 (6)

The top two coefficients of A_m cancel as needed to land in degree≤2m−2.

These identities are positive transport laws for full Gaussian classes, but they do not automatically propagate the moving source-depth flag. For example, m=16 has d_16=9 and d_17=10; θ=x^9 is admissible at row16, but aθ has degree11 and is inadmissible at row17. Similarly A_m preserves the degree of low-degree monomials rather than lowering it, so a shrinking flag can reject its output. This cheap hostile prevents claiming a holonomic all-row repair merely from (4)–(6).

The rational cone carrier instead addresses exactly the moving slope-1/2 flag. Its failure mode is concrete: a positive-weight pole in the reduced generating denominator. Testing whether actual source recursion creates such a pole is a sharper next target than an unrestricted analogy with Gaussian invariance.

## Reproduction

From the session packet root, run `python3 planar/check_extension.py`. Tests check new identities and representative positive/hostile rational carriers; the all-parameter theorem follows from the proof above, not the finite tests.
