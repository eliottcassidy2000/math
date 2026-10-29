# Minimal recurrences in port arity

**PROVED elementary below; FINITE-EXACT support and Hankel controls in
memory.py.** All sources in this note are uniform group elements. Port
arity is the variable in the recurrence; tree depth is held fixed.

## M1. Positive finite support gives the exact recurrence order

Suppose X has distinct support S and strictly positive rational masses w_a.
Set m_r=E X^r for r>=0, with 0^0=1. Then

    Q_S(E)m=0,             Q_S(t)=product_(a in S)(t-a),        (M1)

where E shifts r by one. This is the minimal constant-coefficient linear
recurrence valid from r=0, of order |S|. It is also the smallest dimension
of a linear realization m_r=u M^r v. If the sequence starts at r=1 instead,
remove the zero support point, if present; it contributes only at r=0.

Proof: each a^r is annihilated by E-a. For d=|S| the first d by d Hankel
matrix factors as V diag(w_a) V^T, where V_(i,a)=a^i. Its determinant is

    product_a w_a * product_(a<b)(b-a)^2 > 0.

Every linear realization has Hankel rank at most its dimension, and the
diagonal realization on S attains d. This is the finite positive version of
[THM-2810, factorial-Hankel faithfulness](../../../01-canon/theorems/THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go.md).
For the general weighted-automaton theorem see
[Kiefer, section 4](https://arxiv.org/html/2009.01217v1).
The result is about linear memory, not bit complexity, algebraic formulas,
or a lower bound for every conceivable computation.

## M2. The numerical semigroup hidden in Paley fixed points

For p>3 and the affine wreath actions in FOURIER, let X_k=fix(g), g uniform
in A_k. At the root the fixed count has support {0,1,p}, each with positive
probability. The branching rule is a sum of that many independent copies
of X_(k-1). At every depth k>=1,

    supp(X_k)={a+pb : a,b>=0, a+b<=p^(k-1)}.                  (M2a)

Proof: independent bottom-level choices with all upper root actions identity
give this whole set. Conversely root-fixed counts 0 and 1 are contained in
the support of a sum of p child counts by padding with zero; induction proves
that no additional counts occur. Positivity prevents support cancellation.

Set N=p^(k-1). Reflect the support at its maximum pN. Then

    pN-supp(X_k) = <p-1,p> intersect [0,pN].                 (M2b)

Indeed the left side consists of a(p-1)+bp with a+b<=N. Any nonnegative
representation by p-1 and p can replace p copies of p-1 by p-1 copies of p
until 0<=a<p. If a+b>N, its value p(a+b)-a exceeds pN. Thus the length
constraint adds no further restriction in this interval.

The semigroup <p-1,p> has exactly (p-1)(p-2)/2 gaps, all below p^2.
An elementary proof is to reduce modulo p-1: in residue r, 0<=r<p-1,
the first semigroup element is rp. There are r missing nonnegative values
in that residue. Summing r gives the genus; the last gap is p^2-3p+1.
This is the consecutive-generator case of the classical two-generator
formula, recorded in [Eliahou--Ramirez Alfonsin, theorem 2.1](https://imag.umontpellier.fr/~ramirez/semigroups.pdf).

Combining M1 and M2 proves the minimal arity recurrence orders

    d_(p,1)=3,
    d_(p,k)=p^k+1-(p-1)(p-2)/2,             k>=2.            (M2c)

| p | depth 2 | depth 3 | depth 4 |
|---|---:|---:|---:|
| 7 | 35 | 329 | 2387 |
| 11 | 77 | 1287 | 14597 |
| 19 | 209 | 6707 | 130169 |

These are proved orders for specified affine actions. memory.py checks
Boolean branching support against independent coin reachability at every
displayed depth. At p=7,11 and depth two it constructs the positive integer
histogram, verifies the characteristic recurrence on successive windows,
and certifies full Hankel rank modulo 1,000,000,007. The graph/Hankel upper
bound and this nonzero modular minor together certify the rational rank.

## M3. Ternary centralizers have a particularly explicit recurrence

For the recursive C3 family, with k>=1,

    supp(X_k)={0,3,6,...,3^k},
    supp(Y_k)={3,9,15,...,3^k},              Y_k=fix(c), c in C_k.

The first follows by threefold sums and the zero-root event. The second
follows from Y_k=Y_(k-1)+2X_(k-1), with the two factors independent and
Y_0=X_0=1. Therefore the minimal arity orders are respectively

    3^(k-1)+1,                 (3^(k-1)+1)/2.

At depth three, the 27-vertex centralizer port sequence b_r satisfies

    b_(r+5) -75 b_(r+4) +2070 b_(r+3) -25650 b_(r+2)
             +136809 b_(r+1) -229635 b_r =0,          r>=0.

Its characteristic roots are 3,9,15,21,27, all with positive masses; five
is minimal. This reverses the usual viewpoint in the repo: fixing arity
gives a polynomial in depth, but fixing depth and allowing every arity
forces exponentially growing linear memory. It does not contradict the
bounded-arity Fourier savings in F3.

## Boundary of the mechanism

[THM-2070, horizontal Wick embedding and dihedral cancellation](../../../01-canon/theorems/THM-2070-horizontal-wick-embedding-and-dihedral-return-cancellation.md)
proves why a feasible support need not survive a signed or complex-weighted
observation. Here all actual fixed-count masses are positive, so every
support point contributes to the Hankel rank. Replacing these masses by a
signed character sum requires the explicit cancellation analysis in FOURIER.

The semigroup is a support device, not a tournament, and its additive
generators do not encode joint attachment incidence. This is the exact
coordinate map: fixed count x -> deficit p^k-x; it preserves occurrence and
destroys group elements, path endpoints, and transporter multiplicities.
