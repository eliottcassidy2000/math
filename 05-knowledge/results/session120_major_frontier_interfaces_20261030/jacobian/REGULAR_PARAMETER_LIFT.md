# Complete regular-lift classification for a monomial parameter degeneration

**PROVED scoped parameter theorem, independently audited separately.** This
extends the attributed fixed-fibre image classification, not the Jacobian
conjecture. No literature-priority claim is made. The parameter t is not a
Poisson variable; the bracket is in x,y only, with {P,Q}=P_x Q_y-P_y Q_x.

Let K be any characteristic-zero field, m>=2, d=m-1, and B in K[y]. Put

    P_t=t y^m x+B(y).

Define the ideal of fixed marked responses admitting a regular lift by

    J={R(T) in K[T]: exists Q in K[t,x,y] with {P_t,Q}=R(P_t)}.

It is an ideal since multiplying Q by H(P_t) multiplies its bracket by H(P_t).
The target R is independent of t. Allowing arbitrary t-dependent responses is
a different problem.

## Theorem

The complete classification is

1. If B=beta is constant, J=(T-beta).
2. If B=beta+c y^s is nonconstant, c!=0 and s divides d, then
   J=(T-beta)^(d/s).
3. In every other case J=0.

In case2 a polynomial companion for the ideal generator is

    Q=[c^e-(c+t y^(m-s)x)^e]/(d t),   e=d/s.

The numerator is divisible by t. A companion for any multiple of the generator
is obtained by multiplication by its polynomial quotient evaluated at P_t.
In case1, Q=-yx/d is a companion for T-beta.

## Necessity for nonconstant B

Assume R!=0 has such a companion. Over K(t), the fixed-fibre classification
and coefficient comparison in P_t give

    Q= -R(P_t)/(d t y^d)+H_t(P_t),   H_t(T) in K(t)[T].       (L1)

This formula follows directly from (d/dy)(-1/(d t y^d))=1/(t y^m).
For precision, write Q as a polynomial in P_t over K(t,y); then each coefficient
has derivative fixed by the corresponding coefficient of R. Constants of y
derivation are K(t), so the integration constant H_t is polynomial in P_t,
not an arbitrary rational function of it.

Generic polynomiality in x,y also forces R(B(y)) divisible by y^d: all positive
x-degree terms of R(B+t y^m x)/y^d have nonnegative y-degree, since m>d; its
constant-in-x term is R(B)/y^d. The H_t(B) term is polynomial in y and cannot
cancel a negative y power. Put F(y)=R(B(y))/y^d in K[y].

Expand the finitely many coefficients of H_t as Laurent series at t=0. If its
smallest order is -r<0 with leading polynomial h(T)!=0, substitution P_t=B+O(t)
has leading coefficient h(B(y))!=0, because B is nonconstant. Thus substitution
cannot raise the t-order. In (L1) the first term has a nonzero t^-1 coefficient
-F(y)/d. It follows that r=1 and cancellation requires

    F(y)/d = h(B(y)),     hence R(B(y))/y^d belongs to K[B(y)].   (L2)

Since R(B)!=0, equation(L2) implies y^d belongs to the field K(B). Write
 y^d=U(B)/V(B) with coprime U,V in K[T]. Coprimality remains after substitution
(Bézout's identity), so V(B) divides the polynomial U(B) only if V is constant.
Consequently y^d=F_0(B) for a polynomial F_0.

Differentiate:

    d y^(d-1)=F_0'(B) B'.

Every root of B' over an algebraic closure is therefore zero. In characteristic
zero this forces B=beta+c y^s (including s=1 when B' is constant). Degrees in
 y^d=F_0(B) then force s|d. From (L2), with e=d/s,

    R(B)/y^d = c^e R(B)/(B-beta)^e

is polynomial in B exactly when (T-beta)^e divides R. This proves necessity.
All reasoning descends to K, since beta=B(0) and c are coefficients of B.

## Sufficiency and the constant boundary

For B=beta+c y^s and e=d/s, P_t-beta=y^s(c+t y^(m-s)x). Thus

    Q= [c^e-(P_t-beta)^e/y^d]/(d t)

is the displayed polynomial. Its second term is the particular companion in
(L1) for R=(T-beta)^e; its first term is constant in x,y. Hence its bracket is
exactly R(P_t). Multiplying by H(P_t) proves sufficiency for the whole ideal.

If B=beta is constant, the generic image ideal requires R(beta)=0. Conversely
for R=T-beta the polynomial -yx/d satisfies

    {t y^m x+beta, -yx/d}=t y^m x=P_t-beta.

Multiplication gives every response divisible by T-beta. This proves case1.
The Laurent injectivity argument for nonconstant B was not used in this case.
Indeed the uniform kernel here is larger than K[t,P_t]: (P_t-beta)/t=y^m x
is a polynomial first integral. No false uniform-kernel identification is needed.

## What the theorem resolves, and what it does not

For m=5 and B=y^3, s=3 does not divide d=4, so J=0. Nevertheless R=T^2 has a
polynomial companion on every individual numerical fibre, as the companion
paper proves. For B=y^2, s=2 divides4 and J=(T^2), with the explicit regular
companion -y^3x/2-t y^6x^2/4. Thus the distinction is a complete divisibility
classification in this degeneration family, not merely one failed example.

This is a regularity and marked-response theorem. It neither supplies a
nonconstant-Jacobian counterexample nor gives a general chart-entry reduction
for JC(2). The inherited nonzero-A fixed-fibre sector remains attributed to
THM-3326/3348; the present theorem asks the stricter joint-parameter question.
