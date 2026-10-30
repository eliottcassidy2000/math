# Independent theorem audit: regular parameter lifts

**Verdict: VALID**, over every characteristic-zero field K. This is a
restricted polynomial-derivation theorem, not a general Jacobian-conjecture
claim. Let m>=2, d=m-1, P_t=t y^m x+B(y), and let

    J_B={R(T) in K[T]: there is Q in K[t,x,y] with {P_t,Q}=R(P_t)}.

Then, with beta=B(0), the complete classification is:

1. B constant: J_B=(T-beta).
2. B=beta+c y^s, c nonzero, s positive and s divides d:
   J_B=(T-beta)^(d/s).
3. Every other B: J_B=0.

The bracket differentiates x and y only; t is a parameter. This statement
concerns fixed target R in K[T], not arbitrary t-dependent targets.

## Necessity, with all cancellation issues explicit

Work first over the characteristic-zero field L=K(t). The independently
audited affine-in-x theorem gives every possible companion in L[x,y] as

    Q=-R(P_t)/(d t y^d)+H_t(P_t),    H_t in K(t)[T].          (1)

For nonconstant B, expand the finitely many coefficients of H_t in Laurent
series at t=0. If their maximum pole order r is greater than one, the leading
coefficient after substitution is H_{-r}(B(y)). It is nonzero: substitution
K[T]->K[y], T->B(y), is injective. The other summand in (1) has no pole of
order greater than one, so such a leading pole cannot cancel. Hence r<=1.
The coefficient of 1/t now requires

    R(B(y))/y^d = d H_{-1}(B(y)) in K[B(y)].                (2)

This treats arbitrary t-dependent coefficients in the kernel term; assuming
that H_t itself was polynomial in t would miss the relevant cancellation.
No control at other t-poles is needed for necessity: regularity at zero is
already necessary for a polynomial companion. Possible failures of generic
polynomiality likewise cannot invalidate this necessary implication.

Assume R nonzero. Then R(B) nonzero. Equation (2) gives y^d in K(B). Write
this as y^d=F(B) for F in K(T). In fact F is polynomial. For otherwise, after
passing to an algebraic closure, a finite pole alpha of the reduced rational
function F has a preimage under the nonconstant polynomial B, giving a
finite pole of F(B). This contradicts F(B)=y^d. Coprimality prevents any
cancellation at a preimage. Polynomiality descends to K because F was in
K(T).

Now y^d=F(B). Every root alpha of nonconstant F must have all roots of
B-alpha equal to zero. Hence alpha=B(0)=beta and B-beta=c y^s for some
positive s and c in K*. Also F has its unique root beta, so

    F(T)=c^(-e)(T-beta)^e,      d=se.

This proof needs no assumption that K is algebraically closed: the root
argument is made after extension, and beta=B(0), c the leading coefficient,
and all exponents descend explicitly. Finally (2), together with injectivity
of substitution, says R(T)=F(T)S(T) for S in K[T]. Therefore R is divisible
by (T-beta)^e, exactly as claimed.

If B is constant, the preceding injectivity argument is intentionally NOT
used. Instead the generic affine-in-x classification over K(t) already
requires R(beta)=0. This is enough for necessity in that case.

## Sufficiency and explicit generators

For B=beta+c y^s with s|d and e=d/s, put

    Qgen=[c^e-(c+t y^(m-s)x)^e]/(d t).                    (3)

The numerator is divisible by t, so this is polynomial in K[t,x,y].
It equals - (P_t-beta)^e/(d t y^d)+c^e/(d t), hence its bracket with
P_t is exactly (P_t-beta)^e. Every permitted R=(T-beta)^e S(T) has
companion Qgen S(P_t), by the derivation product rule.

For constant B=beta, the generator has companion Qgen=-xy/d, because

    {t y^m x+beta,-xy/d}=t y^m x=P_t-beta.

Again multiply by S(P_t) for any R=(T-beta)S(T). Thus generic polynomiality,
absence of all other parameter poles, and global existence are established
constructively, not left as unverified assumptions.

## Boundary and kernel cautions

* A nonconstant B=beta+c y^s with s not dividing d has zero regular-lift
  ideal even when each fixed fibre admits the target under discussion.
* The zero target always has Q=0; the nonzero-R hypothesis was only used to
  infer y^d in K(B).
* No assertion that the uniform kernel is K[t,P_t] is made. For constant B,
  y^m x=(P_t-beta)/t is already a polynomial kernel element, demonstrating
  why that assertion would be false.
* The construction specializes harmlessly at t=0, unlike a merely generic
  primitive. Characteristic zero is used in d, derivatives and coefficient
  constants throughout.

`check_regular_lift.py` independently checks the explicit generators and
multiples using exact three-variable polynomial arithmetic. It also tests
membership y^d in K[B] by leading-term elimination, independently of the
claimed shape classification, over all 240 B of degree 1–4 with coefficients
in {-1,0,1} and nonzero leading coefficient, for d=1–8. This is a finite
regression universe, not the proof of the theorem. Normal and optimized
Python runs are compared byte-for-byte.
