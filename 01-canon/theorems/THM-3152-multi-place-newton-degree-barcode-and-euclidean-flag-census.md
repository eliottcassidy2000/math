---
id: THM-3152
title: "Multi-place Newton degree barcodes and Euclidean-flag closure through r=1098"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  A common rational factor of several polynomial rows has degree in a finite
  denominator-weighted Newton-capacity barcode at every p-adic place; an
  empty multi-place intersection proves coprimality.  For the resonant exact
  quadratic factorial pair, adjoining the exact first Euclidean remainder
  closes all 56 three-exit residual rows with 1001<=d<=1100.  Consequently
  every exact three-term quadratic factorial window beginning at 1<=r<=1098
  contains a nonzero moment.  This is not an all-height theorem or FC(3).
audit: >
  The first hostile pass found that finite Newton slopes omit a common
  coordinate factor; the statement and both referees now retain the exact
  zero-root capacity.  Independent hull extraction, polynomial construction,
  and degree-set dynamic programming rederived the Euclidean quotient,
  reproduced both d=1001 two-place certificates, all 56 residual closures,
  and every narrowing-prime trace.  Planted factors v and v+1 survive the
  main observer.  Fresh normal, optimized, stored, and independent transcripts
  agree with their declared hashes.
source: root/frontier-synthesis/euclidean-newton-beach-2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3131-prime-resonance-newton-slope-separation
  - THM-3142-prime-power-predecessor-newton-separation-and-composite-window-census
  - THM-3143-two-step-prime-resonance-euclidean-newton-separation
related:
  - THM-3146-three-step-prime-resonance-full-euclidean-newton-separation
  - THM-3148-fixed-offset-frobenius-endpoint-resultant-classification
script: 04-computation/factorial_multiplace_newton_degree_barcode_thm3152.py
output: 05-knowledge/results/factorial_multiplace_newton_degree_barcode_thm3152.out
script_sha256: f804d3996abe4df981dbf7db877af4aeca9218df64b0ac382af876a3cdca15a0
output_sha256: 6404e71fe13a79188612f55c25a6fef8ea24d65831ac05a9c774dff93ae81377
independent_script: 04-computation/factorial_multiplace_newton_degree_barcode_independent_audit_thm3152.py
independent_output: 05-knowledge/results/factorial_multiplace_newton_degree_barcode_independent_audit_thm3152.out
independent_script_sha256: 427ba0e8c7e5b25efdf56b00200ab201b534f150c6702d3d0ba4a1ce956f32de
independent_output_sha256: 8537f5dd5dcd6aa6a73ecbad78e582f75b857287eb9bcd8316da3d9ef86a907c
hash_basis: LF-normalized bytes
---

# THM-3152 -- multi-place Newton degree barcodes and Euclidean-flag closure through r=1098

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

Let

~~~text
L(t^m)=m!,
q(t)=a+bt+ct^2,                         abc!=0.              (1)
~~~

For every integer \(1\le r\le1098\), the three consecutive factorial
moments

~~~text
L(q^r),                 L(q^(r+1)),                 L(q^(r+2)) (2)
~~~

cannot all vanish.

The abstract mechanism is a local common-factor theorem.  At each rational
prime, retain not only the lower Newton slopes shared by several polynomial
rows, but also each slope denominator, its common horizontal capacity, and
the coordinate-root capacity.  The support of a finite product then contains
the degree of every common rational factor.  Empty support intersection
across places is a rigorous coprimality certificate.

The uniform lemma and the first Euclidean row below are proved for arbitrary
nonzero rational polynomial rows.  The range \(r\le1098\) is
**FINITE-EXACT**.  Nothing here asserts a fixed-prime-bank theorem, an
all-height exact-quadratic theorem, arbitrary support in one variable,
\(SFC(3)\), or the full three-variable Factorial Conjecture.

## 2. Inheritance and the resonant pair

THM-3124 proves that a bad window at \(r\) forces

~~~text
d=r+2,                        b/a=-1/d.                       (3)
~~~

After division by \(a\), put \(v=dc/a\) and

~~~text
A_k=A_k^(d)(v)=L((d-t+vt^2)^k),
n=d-2,
P=A_n,                         Q=A_(n+1).                    (4)
~~~

The window is bad exactly when \(P,Q\) have a common complex root.  Since
they lie in \(\mathbb Q[v]\), this is equivalent to a nonconstant common
rational factor.

The THM-3124 recurrence gives a constant-state exact constructor:

~~~text
A_0=1,                         A_1=d-1+2v,

A_(k+1)=d^k(d-k-1)
       +2(k+1)(2k+1)v A_k
       +k(k+1)(1-4dv)A_(k-1).                               (5)
~~~

Thus the polynomial sequence is generated in \(O(d)\) ring operations,
rather than by repeated multinomial expansion.  This recurrence is the
computational carrier for the finite census.

## 3. Local common-factor degree barcode

Let \(f_1,\ldots,f_s\) be nonzero polynomials in \(\mathbb Q[v]\), and fix a
rational prime \(p\).  For a finite lower Newton slope \(\sigma\) of \(f_i\)
over \(\mathbb Q_p\), let \(\ell_i(\sigma)\) be its horizontal length,
including root multiplicity.  For a slope common to all rows, write

~~~text
sigma=a_sigma/b_sigma in lowest terms,       b_sigma>0,
c_p(sigma)=min_i ell_i(sigma),
z_p=min_i ord_v(f_i).                                         (6)
~~~

Define the finite set

~~~text
D_p(f_1,...,f_s)
 ={m_0+sum_sigma m_sigma:
      0<=m_0<=z_p,
      m_sigma in {0,b_sigma,2b_sigma,...,c_p(sigma)}}.        (7)
~~~

Every \(c_p(\sigma)\) is a multiple of \(b_\sigma\), because the endpoints
of a Newton segment have integral coordinates.

> **Degree-barcode lemma.**  If \(F\in\mathbb Q[v]\) divides every \(f_i\),
> then
>
> ~~~text
> deg F in D_p(f_1,...,f_s).                                  (8)
> ~~~
>
> Consequently, for any finite prime bank \(\mathcal P\),
>
> ~~~text
> intersection_(p in P) D_p(f_1,...,f_s) intersection Z_(>0)
>     = empty                                                   (9)
> ~~~
>
> proves that the rows are coprime over \(\mathbb Q\).

To prove the lemma, first separate the coordinate factor \(v^{m_0}\).  Its
common multiplicity is at most \(z_p\).  Factor the remaining part of \(F\)
over \(\mathbb Q_p\).  All roots of one irreducible factor have the same
valuation.  If that valuation is \(-a_\sigma/b_\sigma\) in lowest terms,
then \(b_\sigma\) divides the ramification index and hence the irreducible
factor degree.  Newton's theorem bounds the total multiplicity at this
valuation in row \(i\) by \(\ell_i(\sigma)\).  The contribution of \(F\) at
\(\sigma\) is therefore a multiple of \(b_\sigma\) no larger than
\(c_p(\sigma)\).  Summing the coordinate and finite-slope contributions
proves (8), and intersecting (8) over primes proves (9).  QED.

The coordinate block in (7) is essential.  The hostile family
\(f_i=v g_i\) would otherwise lose a genuine common degree-one factor because
finite Newton slopes omit the root at zero.  A nonempty intersection remains
inconclusive: it constructs neither a common factor nor compatible residual
polynomials.

## 4. The first Euclidean row

At \(d=n+2\), the leading two terms of the quotient \(Q/P\) are

~~~text
2(n+1)(2n+1)v - 2d(n+1)/(2n-1).                             (10)
~~~

Indeed, the leading coefficient of \(P\) is \((2n)!\), and after cancelling
the leading term the coefficient of degree \(n\) is

~~~text
-4dn(n+1)(2n-2)!.                                           (11)
~~~

Clear the denominator in (10) and define

~~~text
R_1=(2n-1)(Q-2(n+1)(2n+1)vP)+2d(n+1)P.                     (12)
~~~

Then

~~~text
deg R_1<=n-1,
gcd_Q(P,Q)=gcd_Q(P,R_1).                                    (13)
~~~

The second identity holds because (12) is an invertible scalar combination
over \(\mathbb Q\).  Hence every common factor of \(P,Q\) divides all three
rows \(P,Q,R_1\), and the sharper ledger

~~~text
D_p(P,Q,R_1)                                                 (14)
~~~

is lawful.  This is the first **Euclidean--Newton flag**: Euclidean operations
preserve the target common factor, while the Newton barcode probes its degree
through a new valuation projection.  Further exact remainder rows can be
adjoined without changing the abstract lemma.

## 5. Exact first residual witness at d=1001

The first integer not covered by the three inherited exits is

~~~text
d=1001=7*11*13,
d-1=1000=2^3*5^3,
d-2=999=3^3*37,
deg(P,Q,R_1)=(999,1000,998).                                (15)
~~~

For the original pair, the common blocks are

~~~text
p=2:
 (sigma,c,b)=(63/32,32,32),(127/64,64,64),
             (255/128,128,128),(511/256,256,256),
             (1023/512,512,512),

D_2(P,Q)^+={32,64,96,...,992};

p=3:
 (sigma,c,b)=(26/27,27,27),(242/243,243,243),
             (728/729,729,729),

D_3(P,Q)^+={27,243,270,729,756,972,999}.                    (16)
~~~

The two positive degree sets are disjoint, so the barcode lemma alone proves

~~~text
gcd_Q(A_999^(1001),A_1000^(1001))=1.                        (17)
~~~

The Euclidean flag removes the two smallest binary blocks and the smallest
ternary block:

~~~text
D_2(P,Q,R_1)^+={128,256,384,512,640,768,896},
D_3(P,Q,R_1)^+={243,729,972}.                               (18)
~~~

These sets are again disjoint.  Equations (16)--(18) are exact individual
coprimality proofs, not statistics from the bounded census.

## 6. Finite-exact census and closure through r=1098

Call \(d\) a **three-exit residual** when

~~~text
d is composite,
d-2 is composite,
d-1 is not a prime power.                                   (19)
~~~

This is the complement of the uniform exits in THM-3131, THM-3142, and
THM-3143.  There are exactly 56 such integers in

~~~text
1001<=d<=1100.                                               (20)
~~~

For every one, the companion constructs \(P,Q,R_1\) over the integers and
intersects the exact barcodes (14) over the fixed bank

~~~text
{2,3,5,7,11,13,17,19,23,29,31,37,41,43,47}.                (21)
~~~

Every positive intersection is empty.  Thus all 56 pairs are coprime over
\(\mathbb Q\).

Eleven of these 56 integers have \(d-3\) prime:

~~~text
1012,1016,1024,1036,1042,1054,1066,1072,1090,1096,1100.    (22)
~~~

THM-3146 independently closes those rows, leaving 45 residuals relative to
all four currently proved exits.  The computation deliberately retains the
larger 56-row universe, so THM-3146 is related but is not a dependency.

Now combine the pieces.  THM-3142 already closes every \(d\le1000\).  For
\(1001\le d\le1100\), use THM-3131 if \(d\) is prime, THM-3142 if \(d-1\)
is a prime power, THM-3143 if \(d-2\) is prime, and the 56-row census
otherwise.  Hence every resonance \(3\le d\le1100\), equivalently every
window start \(1\le r\le1098\), is closed.  This proves the theorem.  QED.

As a hostile recovery test, the original-pair barcode also closes all 180
rows with \(4\le d\le250\) and \(d-1\) not a prime power, using primes
\(p\le2d\).  That reproduces already closed territory by a different exact
mechanism and is not needed for the new range.

## 7. Exact controls and independent replay

The main companion uses integer polynomial arithmetic, exact \(p\)-adic
valuations, rational slopes, and finite degree-set dynamic programming.  Its
universe and inherited filters are explicit.  It checks:

- recurrence (5) against direct multinomial coefficients for \(d\le15\);
- the quotient and remainder identities for \(3\le d\le30\);
- the discriminant-square positive lift
  \[
  (4d)^n A_n^{(d)}((x+1)/(4d))
    =L(((2d-t)^2+xt^2)^n),                                  \tag{23}
  \]
  whose coefficient of \(x^j\) is a positive exponential integral;
- every 180-row recovery case and all 56 rows in (20);
- the complete two-place data (16)--(18); and
- planted common factors \(v+1\) and \(v\), both of which retain degree one.

The independent companion does not import the main implementation.  It
reconstructs the moment recurrence, Euclidean row, determinant lower hull,
zero-root block, and degree dynamic program.  It reproduces all 56 closures,
every narrowing-prime trace, and both \(d=1001\) certificates.

Run

~~~text
python3 04-computation/factorial_multiplace_newton_degree_barcode_thm3152.py
python3 -O 04-computation/factorial_multiplace_newton_degree_barcode_thm3152.py
python3 04-computation/factorial_multiplace_newton_degree_barcode_independent_audit_thm3152.py
~~~

and compare byte-for-byte with the declared outputs.

## 8. Closed-form language and next frontier

For each place, (7) is exactly the support of the finite generating product

\[
B_p(x)=(1+x+\cdots+x^{z_p})
       \prod_\sigma
       (1+x^{b_\sigma}+x^{2b_\sigma}+\cdots+x^{c_p(\sigma)}). \tag{24}
\]

Thus the barcode is simultaneously a Newton invariant, a restricted
subset-sum language, and an efficiently computable sequence support.  A
bitset multiplication of (24) replaces factorization over \(\mathbb Q\);
multi-place coprimality is support disjointness.  This is the precise survivor
of the broader sequence/closed-form analogy.

The connection contract is

~~~text
source:     exact factorial pair and exact Euclidean remainder rows
target:     supports of the finite products B_p
map:        lower Newton segment -> (slope, capacity, denominator)
preserved:  existence and degree of every rational common factor
lost:       residual polynomials, root phases, cross-slope correlations
sidecar:    coordinate block, several rows, and several p-adic places
test:       empty positive support intersection
~~~

The next all-height problem is digital: derive the blocks in (24) from the
base-\(p\) digits of \(d-2,d-1\), then prove that two or three such digit
languages have empty positive intersection, routing any survivor to the next
Euclidean row.  THM-3148 and THM-3153's second-row lane address fixed
offsets; (24) addresses arbitrary composite residuals.  The two mechanisms
are complementary, but no all-height termination is asserted.

This also states the safe cross-thread lesson.  Bare projected slopes can
collide, just as bare tournament scores or unlabelled operation fibres can
collide.  Retaining the coloured fibre—denominator, capacity, and exact
Euclidean depth—restores a usable obstruction.  No tournament relation is
forced onto the factorial roots.

**End of proof.**
