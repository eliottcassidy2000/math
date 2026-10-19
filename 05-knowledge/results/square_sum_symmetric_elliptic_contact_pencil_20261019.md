# A symmetric elliptic pencil of rational square-sum contacts

**Status.** A genuine non-torsion elliptic translation with a global contact identity. This clears the existence question for a fixed elliptic contact shift, but does not yet give an admissible long AP-paired chain: the exhibited contacts use the same port at every interior petal, and the full real orbit is not confined to the admissible spread region.

Put
\[
f(t)=\frac{4t(t^2-1)}{(t^2+1)^2},\qquad
s=t+u,\quad p=tu,\quad D=s^2+(p-1)^2=(1+t^2)(1+u^2),
\quad J=p^2+6p+1-s^2.
\]
The elementary identity
\[
D^2(1+f(t)+f(u))=D^2+4s(p-1)J
\tag{1}
\]
gives a pencil of symmetric biquadratic curves
\[
H_\lambda(t,u)=4J-2\lambda D-\lambda^2s(p-1)=0.
\tag{2}
\]
Indeed, with \(Y=D+\lambda s(p-1)\),
\[
D^2(1+f(t)+f(u))-Y^2=s(p-1)H_\lambda.
\tag{3}
\]
Every rational point on (2) therefore gives a rational `++` contact. There are rational boundary points
\[
(t,u)=\left(0,\frac{2}{\lambda+2}\right),\qquad
(t,u)=\left(0,\frac{\lambda-2}{2}\right)
\]
when the displayed affine expressions are defined.

## One completely explicit smooth member

For \(\lambda=2\), divide (2) by \(-4\) to obtain
\[
C:\quad (t+2)u^2+(t^2-4t-1)u+2t^2-t=0.
\tag{4}
\]
It is symmetric in \(t,u\). With
\[
v=2(t+2)u+t^2-4t-1,
\]
its smooth projective model is the genus-one quartic
\[
E:\quad v^2=t^4-16t^3+2t^2+16t+1.
\tag{5}
\]
The quartic is squarefree, already modulo 7. The point \(O=(t,u)=(1,1)\), or \((t,v)=(1,2)\), supplies a rational origin.

The vertical root exchange followed by coordinate exchange gives the QRT map
\[
\tau(t,u)=\left(u,\frac{u(2u-1)}{t(u+2)}\right).
\tag{6}
\]
As usual, the rational formula extends across its finitely many apparent indeterminacies to an automorphism of the smooth projective curve. The vertical root exchange is the involution of a degree-two map to the line; the coordinate exchange is another involution with a fixed point \(O\). Each is an elliptic reflection, so their composition is translation by a rational point \(P\). With this choice of origin the coordinate exchange is negation. Consequently the first-coordinate function \(t:E\to\mathbb P^1\) satisfies
\[
t(Q+P)=u(Q),\qquad
1+f(t(Q))+f(t(Q+P))
=\left(1+\frac{2(t(Q)+u(Q))(t(Q)u(Q)-1)}{(1+t(Q)^2)(1+u(Q)^2)}\right)^2.
\tag{7}
\]
This is a function-field identity, not merely a list of sample contacts.

## The translation is non-torsion

Exact point counts for (5), including its two rational points at infinity, are
\[
\#E(\mathbb F_7)=12,\qquad \#E(\mathbb F_{13})=20.
\]
The quartic is squarefree at both primes, so these are good reductions. Prime-to-residue-characteristic torsion injects into good reduction. Combining the two primes shows that every rational torsion point has order dividing 4: the 7-primary part is excluded by the count at 13, the 13-primary part by the count at 7, and every other primary part divides both counts.

However, starting at \(O=(1,1)\), (6) gives the five distinct points
\[
(1,1),\quad (1,1/3),\quad (1/3,-1/21),\quad
(-1/21,23/287),\quad (23/287,5543/8159).
\]
Thus \(P\) is non-torsion.

For a fixed rational initial point, standard elliptic height theory now gives
\(h(t(Q+nP))=O(n^2)\). Rational petal/contact data along this orbit have quadratic bit height. A common multiple of all rational denominators in a prefix of length \(k\) has at most \(O(k^3)\) bits. This arithmetic observation does not imply that these data form a positive AP-paired chain.

## The two remaining failures are explicit

1. Every edge in (7) uses `++`. At an interior petal the two contacts therefore use the same port. Replacing a parameter by another encoding of the same petal does not change this physical incidence. The opposite-port condition needed for the inherited Hamiltonian square-sum cycle is not supplied.
2. The full real orbit cannot simply be declared admissible. For example, \(t=1/3\) has \(f(t)=-24/25\), whose absolute value exceeds the admissibility limit \(1/2\). The displayed orbit already reaches this value at its first non-boundary step.

The rational map (6) acts on an elliptic curve, so the previously proved degree-three obstruction for global rational maps of \(\mathbb P^1\) does not apply. The normalized-x congruent-curve obstruction also does not apply: the quartic coordinate in (5) is a different degree-two elliptic parameter map.

## The second-sign cover has genus nine

The entire nonzero rational pencil has a uniform obstruction to making each edge available with either same-sign port choice. Set
\[
a=\frac{t+u}{tu-1},\qquad z=\frac{tu+1}{tu-1},\qquad v=\frac{t-u}{tu-1}.
\]
On \(H_\lambda\),
\[
8z^2=L_+(a):=2(\lambda+2)(1+a^2)+\lambda^2a,
\qquad
8v^2=L_-(a):=2(2-\lambda)(1+a^2)-\lambda^2a.
\tag{8}
\]
Conversely,
\[
t=\frac{a+v}{z-1},\qquad u=\frac{a-v}{z-1}.
\]
Thus (8) is a birational model. For every nonzero rational \(\lambda\), the two binary quadratics \(L_+,L_-\) each have two distinct geometric roots, and have no common root. Their discriminants factor as
\[
\operatorname{disc}L_+=(\lambda^2-4\lambda-8)(\lambda^2+4\lambda+8),
\quad
\operatorname{disc}L_-=(\lambda^2+4\lambda-8)(\lambda^2-4\lambda+8),
\]
and their resultant is \(64\lambda^4\). None of the discriminant factors has a rational zero. At \(\lambda=\pm2\), one quadratic has a simple root at infinity and a simple root at zero. Therefore the map \(E_\lambda\to\mathbb P^1_a\) is a geometrically connected degree-four biquadratic cover, branched at exactly four values, with two simple ramification points over each. Riemann–Hurwitz gives genus one. The displayed boundary point supplies a rational point.

The first contact square root is
\[
y_+=1+\frac{\lambda a}{1+a^2}.
\]
The additional `--` contact requires a square root of
\[
2-y_+^2=\frac{B_\lambda(a)}{(1+a^2)^2},\quad
B_\lambda(a)=a^4-2\lambda a^3+(2-\lambda^2)a^2-2\lambda a+1.
\tag{9}
\]
For nonzero rational \(\lambda\), this quartic has four distinct roots. Indeed, putting \(w=a+a^{-1}\) gives
\(B_\lambda/a^2=w^2-2\lambda w-\lambda^2\); its two \(w\)-roots are distinct, and neither is \(\pm2\) because that would force \(\lambda\) to be irrational.

None of these four roots is a branch value of (8). Direct elimination gives
\[
\operatorname{Res}_a(B_\lambda,L_+)
=\lambda^4(\lambda^2-8\lambda-16)^2,
\quad
\operatorname{Res}_a(B_\lambda,L_-)
=\lambda^4(\lambda^2+8\lambda-16)^2,
\tag{10}
\]
which are nonzero for rational \(\lambda\ne0\). Equivalently, a short polynomial identity proves the exclusion without invoking a resultant algorithm: for \(L=c(1+a^2)+da\),
\[
c^2B_\lambda-L\bigl(L-2(d+c\lambda)a\bigr)
=(d^2+2cd\lambda-c^2\lambda^2)a^2.
\]
For \(L_+\), the right coefficient is \(\lambda^2(\lambda^2-8\lambda-16)\), and for \(L_-\) it is \(\lambda^2(\lambda^2+8\lambda-16)\). If \(c=0\), its zero/infinity roots are visibly not roots of the binary quartic \(B_\lambda\).

Each of the four zeros in (9) consequently has four distinct unramified lifts to \(E_\lambda\). These are sixteen simple zeros. Its poles are all even: the two possible pole values \(a=\pm i\) have denominator order two and numerator \(\lambda^2\ne0\). There is no pole at infinity. Hence the smooth double cover adjoining the second contact square root has exactly sixteen branch points and genus
\[
2g-2=16,\qquad g=9.
\]
In particular, the second square class has odd divisor and cannot become a square under any isogeny pullback. By Faltings's theorem, for each fixed nonzero rational \(\lambda\) there are only finitely many rational edges on this pencil member supporting both `++` and `--` contacts. This is a noneffective finiteness assertion. It blocks the proposed simultaneous-second-square repair of this pencil; it does not classify all possible alternating constructions using different curves or different parameter functions.

The external theorem used for the last sentence is G. Faltings, *Endlichkeitssätze für abelsche Varietäten über Zahlkörpern*, Invent. Math. **73** (1983), 349–366, Satz 7, p. 365, DOI [10.1007/BF01388432](https://doi.org/10.1007/BF01388432).

## Exact replay

From the repository root run `python 05-knowledge/experiments/verify_symmetric_elliptic_contact_pencil_20261019.py`. It uses only the standard library and the colocated `rational_contact_certificate_core.py`; it checks (1)–(6), the five-point orbit, the contact squares, both good reductions, and both point counts, and the second-sign branch identities. The normal and `python -O` outputs agree byte for byte.

The receipt is [symmetric_elliptic_contact_pencil_controls.json](symmetric_elliptic_contact_pencil_controls.json).
