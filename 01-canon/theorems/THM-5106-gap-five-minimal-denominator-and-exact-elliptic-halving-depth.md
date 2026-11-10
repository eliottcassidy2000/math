---
id: THM-5106
title: Gap-five minimal denominator and exact elliptic halving depth
status: PROVED elementary arithmetic, independently audited; finite divisor uniqueness certificate. Classical maps synthesized without a priority claim. No Chan–Pak consequence.
source: 05-knowledge/results/session242_congruent_descent_20261110/arithmetic/ARITHMETIC.md
session: Session242, 2026-11-10
depends_on: []
---

# THM-5106: why12,31,41,49 give5, and where halving stops

Let N>0 be squarefree. If positive rational squares
(p/d)^2,(q/d)^2,(r/d)^2 form a strictly increasing AP with common difference N,
with positive integer p,q,r,d, then

    (4 if N odd, else2) · (3 if 3 does not divide N, else1) divides d.

Square residues modulo16 force8|Nd²; modulo3 forces3|Nd². Hence for N=5,
12|d. This is sharp, since41²−31²=49²−41²=720=5·12². At denominator12,
(31,41,49) is the unique positive increasing numerator triple: exhausting the
positive divisor pairs (q−p)(q+p)=720 with q²+720 square proves uniqueness.
This is a least-denominator theorem, not a bound on all congruent witnesses.

For any such AP, define a=(r−p)/d, b=(r+p)/d, c=2q/d and

    P=(N(a+c)/b,2N²(a+c)/b²),
    Q=(q²/d²,pqr/d³) on E_N: y²=x³−N²x.

Then Q=−2P. Indeed the triangle identities yield
x(2P)=c²/4 and y(2P)=c(a²−b²)/8. Thus our packet has
P=(25/4,75/8), Q=(1681/144,62279/1728).

Furthermore Q∈2E_5(Q) but Q∉4E_5(Q). Duplication makes each of x,x−5,x+5
a rational square. A putative half of−P must therefore pass those square tests.
All four possible halves of Q differ from−P by rational two-torsion; their
x-coordinates are25/4,−4,45,−5/9. None passes the necessary tests, so none is a
double. No Mordell–Weil rank, generator claim or external2-descent criterion
is required.

The two geometric maps have different positivity chambers. The direct AP
point (X,Y) gives square-sum labels ((X−2N)/2,X/2,(X+2N)/2), positive iffX>2N.
Applying the inverse triangle map and then the quadratic clique construction
instead requires (X²+N²)²>8NX(X²−N²). Our Q passes the first and fails the
second. Neither test may silently replace the other.

The canonical Berggren parent of(9,40,41) is(7,24,25), changing area
squareclass5 to21 and losing positive clique labels. The integer-height descent
preserves the Pythagorean equation, not the fixed area class or chamber.
This is a scoped obstruction to that proposed route, not all possible descents.

[Full proof and provenance](../../05-knowledge/results/session242_congruent_descent_20261110/arithmetic/ARITHMETIC.md)
include exact formulas, normalization scales6/12/24, source paths and331
independent primitive-triangle identity checks. [Independent audit](../../05-knowledge/results/session242_congruent_descent_20261110/audit/REVIEW.md).
The classical congruent-number correspondence is credited, not claimed new.
No arithmetic-to-Hamiltonian-path-count map or Chan–Pak proof is asserted.
