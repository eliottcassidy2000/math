# A 2-regular Hamilton-path bank and an algebraic irrational series

**Scope.** This note applies to the Session 158 strong transitive-deletion (one-apex) bank. It supplies a finite-state description, an explicit characteristic-two algebraic irrational series, exact dyadic moments, and an alternative single-block proof of a congruence theorem already established in Session 157. It does not prove exact integer-spectrum coverage, and no claim of novelty against the literature is made.

Sources: Session 158's [binary bank](../../session158_creative_20261103/binary_bank/NOTE.md) and Session 157's [strong one-apex modular universality](../../session157_gaps_20261103/universal/NOTE.md). In particular, Section 3 re-proves a result of S157 using a different block mechanism; it is not a new universality theorem.

## 1. Finite-state representation

For n>=1 take S158's A(1)=2 and

    A(2n)=A(n)+n,
    A(2n+1)=A(n)+1-(n mod 2).

Set A(0)=1 solely for arithmetic convenience. Both recurrences then hold also at n=0. Put F(n)=n+A(n). For every n>=1, F(n) is the Hamilton-path count of the strong one-apex tournament whose orientation word is the ordinary binary expansion of n followed by 0, as proved in S158. The arithmetic extension F(0)=1 is not used as a tournament witness.

For t=n mod 2, a bit b updates the column vector (A(n),n,t,1)^T by

    b=0:  [[1,1, 0,0], [0,2,0,0], [0,0,0,0], [0,0,0,1]],
    b=1:  [[1,0,-1,1], [0,2,0,1], [0,0,0,1], [0,0,0,1]].

Thus A and F are 2-regular over the integers: the four coordinate sequences generate a finitely generated Z-module closed under the maps x(n) -> x(2n), x(2n+1). For each q, reducing this state modulo q gives a finite automaton, so A(n) and F(n) modulo q are 2-automatic. In fact F(n) is odd for all n>=0: F(0)=1, and each of F(2n), F(2n+1) has the parity of F(n).

The generating function a(z)=sum_{n>=0} A(n) z^n has the formal Mahler equation

    a(z)=(1+z)a(z^2) + z^2/(1-z^2)^2 + z/(1-z^4).

This follows by summing the even and odd recurrences, using sum n x^n=x/(1-x)^2 and sum (1-(n mod 2))x^n=1/(1-x^2).

## 2. Exact dyadic means

For every integer k>=1, writing M=2^k,

    sum_{0<=n<M} A(n) = M^2/4+M,
    sum_{0<=n<M} F(n) = 3M^2/4+M/2.

Indeed, S_1=3, and splitting the next range into evens and odds gives S_k=2S_(k-1)+2^(2k-3) for k>=2. The F identity adds M(M-1)/2. The n=0 extension is included as stated; exclude its contribution 1 to get sums over actual strong witnesses 1<=n<M. These are means over indexed words, not averages over distinct counts.

## 3. A single repeating word covers every permitted residue

This is an alternative proof of S157's modular universality result. Fix q>=1 and write q=2^v r with r odd. Choose j>=max(1,v) such that 2^j=1 mod r (vacuously for r=1), put Q=2^j, and let B be the fixed bit word 0^j1^j. From any positive binary prefix n with accumulator A, appending B changes

    n -> Q^2 n+Q-1,
    A -> A+(Q-1)n+1.                                      (1)

The zeros add n+2n+...+(Q/2)n=(Q-1)n, and the first 1 adds one. Starting at prefix 1, write n_t for the value after t copies of B. After the first copy, n_t=-1 mod 2^v and n_t=1 mod r; both residues remain fixed thereafter. Therefore each subsequent block adds the same D to A modulo q, where

    D=2 mod 2^v,       D=1 mod r.

Thus gcd(D,q)=gcd(2,q). For t>=1, F(n_t) travels through all residue classes modulo q if q is odd, or all odd residue classes if q is even. Appending one final 0 converts every prefix into a strong tournament witness. The n_t are strictly increasing and F(n_t)>=n_t+2, so every eligible residue is attained by infinitely many distinct counts.

This gives a finite, explicit regular-language family 1 B^* 0 for each modulus, not a spectrum classification. The first q/gcd(2,q) values after t=1 already cover the eligible residues. The constructive bound on the required number of vertices is 3+2j(1+q/gcd(2,q)), where j can be chosen as a suitable multiple of ord_r(2), with ord_1(2) interpreted as 1. S157 instead uses a two-adic suffix followed by a separate odd-prime CRT/lifting argument.

For reference, starting n_0=1, the exact count is

    F(n_t) = 3 + ((Q+2)^2/(Q+1)^2)(Q^(2t)-1) + 2t/(Q+1).

This expression is an integer, as is clear from (1). Its geometric-plus-linear shape is a concrete bridge between fixed substitutions and the count spectrum.

## 4. A quadratic algebraic irrational in characteristic two

Define g(n)=((F(n)-1)/2) mod 2 for n>=0, and g(z)=sum g(n)z^n in F_2[[z]]. Since F is odd, g is well-defined. It is 2-automatic because it is an output map of F modulo 4. It is not eventually periodic: given any proposed period p>=1, use the block B in Section 3 for q=4p. Once one block has been appended, the prefixes n_t all agree modulo p, while each following block adds D=2 mod 4 to F(n_t). The n_t increase without bound, contradicting the proposed eventual period. Consequently g(z) is not rational over F_2(z), since a rational power series over a finite field has eventually periodic coefficients.

The integer recurrences imply, modulo 2,

    g(2n)=g(n)+(n mod 2),
    g(2n+1)=g(n)+floor(n/2)+1.

Splitting g(z) into even and odd terms and using g(z^2)=g(z)^2 in characteristic two yields

    g(z)=(1+z)g(z)^2 + z(z^3+z+1)/(1+z)^6.

Equivalently,

    (1+z)^7 g(z)^2 + (1+z)^6 g(z) + z(z^3+z+1)=0.

Because g is not rational, this polynomial is irreducible over F_2(z), and g has algebraic degree exactly two. This is an explicit algebraic irrational derived from Hamilton-path count residues; it is not a statement that the integer-valued counts cover all odd numbers.

## Reproducible controls

`python3 check.py` checks 127 actual tournaments by independent subset/endpoint dynamic programming, exact dyadic sums through k=12, and modular block calculations for 2<=q<=500. It checks finite cases only; the proofs above establish all stated universal claims.
