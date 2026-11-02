# Weighted zeta first loses Hamiltonian stability at order five

Status: PROVED by the algebra below; exact rational witness independently auditable.
Baseline: `291b4d9757842c02aa451beb65acfade2848609b`.
No literature-priority claim. This is a recovered-mechanism extension of the
June 15 binary spectral/correlation distinction, not progress on RH.

## Definitions and theorem

A weighted tournament is a real matrix A with Aii=0, Aij>=0 and Aij+Aji=1
for i!=j. Set P_A(z)=det(I-zA), Z_A(z)=1/P_A(z), H(A)=sum over ordered
Hamiltonian paths of their edge-weight products, V=sum_{i<j}Aij Aji,
tau=sum of directed triangle weights, u=det(I+A), v=(-1)^n det A,
Delta=H-u-|v|. Cycles are counted up to cyclic rotation.

1. At every order n<=4, P_A determines H and Delta.
2. At order five, P_A together with the single incidence scalar

       M23(A)=sum_{edges e} Aij Aji * tau(A[e^c])

   determines H and Delta exactly, via

       H=1+2tau+2M23-2[z^5]P_A,                 (1)
       tau=-[z^3]P_A,  V=-[z^2]P_A.

3. P_A alone does not determine either H or Delta at any order n>=5,
   even with V and tau also supplied (they are already spectral).
   At n=5 there is a continuous two-parameter family of exact cospectral
   witnesses with the same skew spectrum and all spectral/walk moments.

The claim is about exact reconstruction. It does **not** rule out spectral
sufficient lower bounds, spectral ingredients in a stability proof, or a
spectral proof of a particular inequality.

## Proof of the reconstruction boundary

The weighted odd-cycle formula gives H=1+2tau at n=3,4 (and H=1 at n=1,2).
The determinant coefficient [z^3]P_A is -tau, so statement 1 follows.
At n=5 there cannot be two vertex-disjoint odd cycles, so
H=1+2tau+2C5. The only nonzero full cycle covers are a directed 5-cycle
or a disjoint 2-cycle and 3-cycle. Their determinant signs give

       [z^5]P_A=M23-C5.

This proves (1), hence statement 2. This is precisely the first order at
which a digon/triangle correlation can hide in the determinant. In binary
tournaments every digon weight is zero, so M23=0 and this obstruction is
absent; the old binary first-blindness order is six, not five.

## A continuous isospectral weighted family

Use indices 1,...,5. Let K(a,b) be skew with only

       K14=a,  K25=b,  K41=-a,  K52=-b

nonzero. Put A=(J-I+K)/2 and

       w=(1,1,-2,0,0)^T,
       Q=I-ww^T/3,
       B=QAQ^T=(J-I+QKQ^T)/2.

Because w^Tw=6 and w^T1=0, Q^T=Q, Q^2=I, Q1=1. Thus A and B
are orthogonally similar, with the same characteristic polynomial, P,
Z, u, v, skew characteristic polynomial, and all closed/open walk moments
1^T A^k1. For |a|,|b|<=1 both A and B are weighted tournaments: in
QKQ^T the only possible nonzero upper entries are

       (14,24,34)=(2a/3,-a/3,2a/3),
       (15,25,35)=(-b/3,2b/3,2b/3).

All these have absolute value <=2/3.

For a short hand-check put x=(2a,-a,2a)/3 and y=(-b,2b,2b)/3.
The nonconstant entries of QKQ^T are Ki4=xi and Ki5=yi. Each triangle
weight is one quarter of 1 plus its three cyclic K-products. Thus the
edge/complement sum is

    16 M23(B)=10-||x||^2-||y||^2-x.y
      -sum_(j<k)(xj xk+yj yk)
      +sum_i (xi^2 yj yk+yi^2 xj xk), {i,j,k}={1,2,3}.

Here x.y=0, both pair sums vanish, the norms are a^2,b^2, and the two
quartic sums are each 2a^2b^2/27. For A every triangle has weight 1/4.
Consequently

       M23(A)=(10-a^2-b^2)/16,
       M23(B)-M23(A)=a^2 b^2/108.

Equation (1) then gives the exact identity

       H(A)=15/2,
       H(B)-H(A)=Delta(B)-Delta(A)=a^2 b^2/54.       (2)

These finite polynomial identities are also checkable by independently
summing the 120 Hamiltonian path monomials. A reproducible exact witness
uses a=1/4,b=1/2 and gives the difference 1/3456.

## Explicit rational witness

A has all offdiagonal entries 1/2 except
A14=5/8,A41=3/8,A25=3/4,A52=1/4. B is

    [0,    1/2,  1/2,  7/12, 5/12]
    [1/2,  0,    1/2, 11/24, 2/3 ]
    [1/2,  1/2,  0,    7/12, 2/3 ]
    [5/12,13/24, 5/12, 0,    1/2 ]
    [7/12, 1/3,  1/3,  1/2,  0   ].

Their common polynomial is

    P(z)=1-(155/64)z^2-(5/2)z^3-(1019/1024)z^4-(37/256)z^5.

The exact invariants are

| Invariant | A | B |
|---|---:|---:|
| V | 155/64 | 155/64 |
| tau | 5/2 | 5/2 |
| u | 233/1024 | 233/1024 |
| v | -37/256 | -37/256 |
| M23 | 155/256 | 2093/3456 |
| H | 15/2 | 25921/3456 |
| Delta | 7299/1024 | 197081/27648 |

## All larger orders and all good residue rings

Adjoin r>=1 deterministic singleton blocks in ordinal sum. The weighted
Hamilton-path sum is unchanged (all paths must traverse the blocks in
order), P is unchanged (block triangular with a zero-diagonal singleton),
V and tau are unchanged, u is unchanged, and v becomes zero. Thus the same
H and Delta difference survives at every n>=6. This proves statement 3
at all orders, in particular throughout 10,...,1499. It does not violate
any proved lower bound at those orders.

The entries lie in R=Z[1/6]. For every integer modulus m>=2 coprime to 6,
there is a canonical reduction R -> Z/mZ. Orthogonal similarity and the
polynomial identities remain true after reduction; the difference
1/3456 is a unit, hence remains nonzero. Therefore exact polynomial
Hamilton-path evaluation is not determined by the characteristic polynomial
in any of those residue rings, including modulo 5 and modulo 11.
This is algebraic reduction, not an ordered/probabilistic tournament over a
finite ring: positivity, |v|, and real inequalities do not reduce modulo m.
At moduli divisible by 2 or 3 these rational entries have no canonical
reduction; use cleared integral data with its denominator sidecar instead.

## The repaired zeta tangent and a division-free all-modulus version

For a weighted tournament over a commutative ring, determinant expansion
always gives

       P_A(z)=1-Vz^2-tau z^3+O(z^4),
       Z_A(z)=1+Vz^2+tau z^3+O(z^4).              (3)

This identity uses no division by 2 or 3, and holds coefficientwise modulo
EVERY integer modulus for an integral edge-weight model. A real rational
weighted model must first localize at its denominators or retain a clearing
factor. For A over Z_p and x in pZ_p of valuation d>=1,

       v_p(Z_A(x)-1)>=2d,
       equality iff V is a p-adic unit.

Indeed P_A(x) is a unit, and Z_A(x)-1=x^2(V+x tau+O(x^2))/P_A(x).
If V=0 exactly, the deterministic THM-4093 triangle tangent returns:
v_p(Z_A(x)-1)>=3d, with equality iff tau is a p-adic unit.
This is an adjacency-determinant zeta statement, not a Kubota--Leopoldt
p-adic zeta or an irrationality/RH assertion.

## Connection contract to Session138

- Source: old closed-orbit zeta / OCF spectral-correlation split, particularly
  `07-reflections/the-zeta-function-and-the-ocf-read-complementary-halves.md`
  and THM-499/502/505 (June 15); THM-4093 (August 25).
- Target: the exact sharp order-five weighted Delta, V, tau region in
  `session138.../stability/SHARP_LOW_ORDER_STABILITY.md`.
- Map: A -> det(I-zA), enhanced by edge/complement incidence M23.
- Preserved: all determinant coefficients, V, tau, u, v; exact H and Delta
  after adding M23 at n=5.
- Lost: which digon and triangle occupy disjoint support. Full spectrum does
  not restore it; neither does its critical-line/skew symmetrization.
- Sidecar: M23 suffices at n=5; the sharper analytic nonnegative proof
  additionally organizes it with complementary triangle curls and C4,C5.
- Cheapest decisive test: the displayed rational pair, using exact
  permutation determinant and independent Held--Karp path DP.
- Not proved: a reduction of the finite stability gap 10..1499; an all-n
  fixed-dimensional correlation sidecar; or anything about zeros of the
  classical Riemann zeta function.
