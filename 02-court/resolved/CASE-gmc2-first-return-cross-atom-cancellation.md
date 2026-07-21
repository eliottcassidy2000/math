# CASE — first-return primitiveness does not imply atomwise vanishing

**Filed and resolved:** codex-2026-07-21-S83

**Against:** THM-1770(B)--(D), the pair-only/star closure inferred from them,
and any universal `2 × primitive return` cutoff

**Verdict:** UPHELD by exact counterexample; THM-1770(A) survives, (B)--(D)
are retracted

**Correction:** MISTAKE-211; THM-1770 amended in place

## The invalid inference

At the least return level (m_*), every balanced word is indeed a primitive
balanced charge multiset. Distinct atoms use distinct coefficient monomials.
THM-1770 then inferred

\[
\sum_A C_A\,x^A=0\quad\Longrightarrow\quad C_Ax^A=0
\text{ for every atom }A.
\]

That implication would be valid if the left side were asserted to be the zero
polynomial for all coefficient values. The nullcone hypothesis supplies only a
scalar equality at one coefficient point, where distinct monomials can cancel.

## A star-support counterexample at the first return

Let

\[
P=aZ^6+bW^2+cW^{18},\qquad a,b,c\ne0.
\]

The charges are (+6,-2,-18). The least balanced return has length four,
and there are exactly two primitive atoms there:

\[
(+6)+3(-2)=0,
\qquad
3(+6)+(-18)=0.
\]

Using \(\mathbb E[Z^rW^s]=\delta_{rs}r!\),

\[
\mathbb E[P^4]
=4\cdot6!\,ab^3+4\cdot18!\,a^3c.
\]

Set

\[
c=-\frac{6!}{18!}\frac{b^3}{a^2}.
\]

Then every coefficient is nonzero but the two first-return atom monomials
cancel. This is already a star support: one positive charge and two negative
charges. It directly refutes both atomwise separation and the claimed star
closure from the first equation. Later moments may still close this support;
the point is that they are indispensable.

## A separate cutoff failure caused by radial channels

The same issue persists even when all charges are only \(\pm1\). Put

\[
P=W(a+cU)+Z(b+dU),\qquad U=ZW.
\]

Then

\[
\mathbb E[P^{2j}]=\binom{2j}{j}L(H^j),
\qquad H=U(a+cU)(b+dU).
\]

There are nonzero choices with (L(H)=L(H^2)=0) but (L(H^3)\ne0).
Thus the primitive return is (R=2), moments through (2R=4) vanish, and
the first obstruction is at (3R=6). A second exact support with charges
((-1,-4,-1,+2)) has (R=3), vanishes through (2R=6), and dies at
(3R=9). These examples refute the cutoff extrapolation without refuting
the local span-six certificate in commit `888c763ea`.

## What remains valid

1. THM-1770(A): least-return balanced words are primitive, so there are no
   *composite* balanced words at that level.
2. The exact two-character/binomial case: there is only one balanced tuple,
   so no cross-atom cancellation is possible.
3. Named multilevel resultant certificates, including the
   \(\{\pm1,\pm3\}\) constant-coefficient span-six calculation.

What does **not** follow is that pair monomials enter the moment ideal one at a
time. The correct target is

\[
(c_pc_n)^N\in
\sqrt{\langle \mathbb E[P],\mathbb E[P^2],\ldots\rangle}
\quad(p>0>n),
\]

proved by a genuine multilevel elimination. HYP-8765 records the candidate
radial-channel tower.

## Quotient audit

For a monomial (Z^aW^b), charge (q=a-b) is not enough. A balanced word
also has radial height (h=a+b), and after charge balance its Gaussian weight
is

\[
\left(\frac{\sum h_i r_i}{2}\right)!.
\]

The faithful combinatorial state is therefore at least `(charge, height)`,
together with the coefficient monomial. Passing to charge alone preserves
which words survive angular averaging but destroys their factorial weights
and their possible phase cancellation.
