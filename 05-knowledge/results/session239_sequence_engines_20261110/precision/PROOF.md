# Precision-first A000568: bounded cycle counts and two uniform congruences

**Status: PROVED by the elementary argument below; independently audited in
`../audit/`.** Exact numerical verification is separate. This refines the
inherited odd-order valuation in
[THM-305](../../../../01-canon/theorems/THM-305-two-adic-valuation-tournament-count.md).
No external novelty or Chan–Pak conclusion is claimed.

Let a(n)=A000568(n), the number of tournament isomorphism classes, n>=0.
For positive odd n define U_n=a(n)/2^((n-1)/2). THM-305 proves that U_n is an
odd integer; the proof below also recovers that fact while retaining more bits.
All congruences of fractions use their odd denominator inverse in Z/(2^b).

## Theorem 1: a precision-dependent cycle cutoff

For an odd partition lambda=(l_1,...,l_r) of odd n, let
G(lambda)=sum_{i<j} gcd(l_i,l_j), and m_j its part multiplicities. Put

    z_lambda = product_j j^(m_j) m_j!,
    z_lambda = 2^v z_odd,
    delta(lambda) = G(lambda) - (r-1)/2 - v.

Its contribution to U_n in Davis's formula is exactly

    2^delta(lambda) / z_odd.

For every such partition,

    delta(lambda) >= L(r) := (r-1)^2/2 - v2(r!).                 (1)

Consequently U_n modulo 2^b is the sum only over odd cycle counts r with
L(r)<b. L(1)=0, L(3)=1, L(5)=5, L(7)=14, L(9)=25. In particular the first
five bits require at most three cycles, independently of n. The cutoff is
O(sqrt(b)); there is no need to generate all odd partitions of n.

**Proof.** The Davis exponent is (n-r)/2+G. Subtracting the normalization
(n-1)/2 and removing the powers of two in z gives the displayed contribution.
All part lengths are odd, so v=sum_j v2(m_j!). The multinomial integer
r!/product_j m_j! shows v<=v2(r!). Since every gcd is at least one,
G>=r(r-1)/2, proving (1). L(r+2)-L(r)=2r-v2(r+1)>0 for positive odd r.
Also v2(r!)<=r-1, which proves quadratic growth and the stated cutoff.
For r=1 the contribution is 1/n. Every r>=3 has L(r)>=1, so U_n is
2-adically integral and odd. Since a(n) is an integer, this establishes the
inherited exact valuation as well. Terms with delta>=b vanish modulo 2^b.

**Complexity and limitations.** For fixed b let R be the largest admitted
odd r. There are O(n^(R-1)) nondecreasing r-part candidates, and each term
uses O(R^2) gcd operations and arithmetic on O(b+R log n) bits. The supplied
implementation never constructs a(n), whose binary length is quadratic in n.
This is polynomial in numeric n at fixed precision, not polynomial in the
binary input length log n, and it does not speed up every full-count query.
For b=1 the direct answer is 1; for b=2 Theorem2 gives a constant-size formula.
The new response-state modular engine in `../burnside/` is an alternative
that can be faster at moderate n and larger b.

## Theorem 2: the first two normalized bits

For every odd n>=3,

    a(n) / 2^((n-1)/2) == -n (mod 4).                         (2)

This is a uniform arithmetic statement about A000568, not a finite pattern
or a statement about the values H(T).

**Proof.** Modulo4, r>=5 vanish by (1). For r=3 with distinct parts,
delta=G-1>=2, hence these vanish. If lambda=(x,x,y), x!=y, then

    delta = x + 2 gcd(x,y) - 2.

This is less than2 exactly when x=1. For odd n>=5 the only extra term is
lambda=(1,1,n-2), contributing 2/(n-2). All-equal triples x>=3 have
 delta=3x-2>=7 and vanish. Thus U_n == 1/n+2/(n-2) (mod4), which equals
-n because n is odd. At n=3 the distinct-multiplicity formula is inapplicable:
lambda=(1,1,1) has odd denominator3 and contributes2/3. Directly U_3=1
agrees with (2). n=1 is genuinely excluded: U_1=1 but -1==3 (mod4).

## Theorem 3: the first three normalized bits and coprime triples

For odd n>=5, let C_n be the number of unordered triples of distinct positive
odd integers a<b<c, with a+b+c=n and gcd(a,b)=gcd(a,c)=gcd(b,c)=1. Then

    U_n == n^(-1) + 2(n-2)^(-1) + 4 C_n (mod 8).              (3)

**Proof.** Again r>=5 vanish. In a distinct triple, G is an odd integer >=3,
and delta=G-1<3 iff G=3, equivalently the three parts are pairwise coprime.
Each contributes 4/(abc)==4 modulo8. A repeated triple has
 delta=x+2gcd(x,y)-2, an odd integer; delta<3 requires delta=1, hence x=1.
Its contribution is2/(n-2). The all-equal case with x>=3 vanishes as above.
The one-cycle contribution is1/n. This exhausts all partitions.
The n=3 exception is essential: the naive right side would be5 modulo8,
whereas U_3=1. The code handles n=1,3 explicitly.

## Verification, hostile controls and proof leverage

Run from repository root:

    python3 05-knowledge/results/session239_sequence_engines_20261110/precision/precision.py
    python3 -O 05-knowledge/results/session239_sequence_engines_20261110/precision/precision.py

`results.json` and `results.optimized.json` record 312 comparisons against a
separate full integer Davis sum, all positive odd n<=51 and precisions1..12.
They check the inherited valuations, (2), (3), and larger exact mod8 values
at101,301,1001. The large residues are not full A000568 computations.
The audit supplies a separate divisor-square Davis path and checks.

This is the requested algorithm-to-proof move: choose the required precision,
prove which cycle strata cannot contribute, classify the surviving strata,
and obtain a uniform congruence. It is not legitimate to interpret a residue
of the total number of classes as support of a particular Hamiltonian-path
value. `../proof_bridges/` records the extra value information needed there.
