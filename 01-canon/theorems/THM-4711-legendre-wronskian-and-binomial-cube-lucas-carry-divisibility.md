---
id: THM-4711
title: "Lucas-carry Wronskians for Legendre and Zagier sequences and the binomial-cube radical"
status: "PROVED — all primes p >= 7 in the stated families; independently audited"
session: clock_duality_oct18l
related: [HYP-9330, HYP-9311]
---

# THM-4711 — Legendre Wronskian and binomial-cube Lucas carry divisibility

**PROVED, independently audited.** This proves the stated Legendre, six
Zagier, and binomial-cube cases. HYP-9330 remains OPEN in its remaining scope.

## 1. Exact scope and inheritance

For a prime p >= 7 and r in {2,3}, put

```
u_n = binom(2n,n)^r,
F_(r,p)(x) = sum_(k=0)^(p-1) u_k x^k,
V_(r,p)(x) = sum_(k=0)^(p-1) ((u_(p+k)-2^r u_k)/p) x^k  modulo p.
```

The divisions are integer divisions. Then:

1. **Legendre / HYP-9330(B):** F_(2,p) is squarefree and
   gcd(F_(2,p), V_(2,p)) = 1 in F_p[x].
2. **Binomial cube / HYP-9330(C), stronger form:** the monic gcd of
   F_(3,p) and V_(3,p) equals rad(F_(3,p))
   in F_p[x], including roots over the algebraic closure and the singular
   root 64x=1 when p=3 modulo 4.
3. **Six sporadic cases / HYP-9330(B):** for each of Zagier A–F,
   F_p is squarefree and gcd(F_p,V_p)=1 at every prime p>=7; see §7.

The proof gives stronger polynomial identities, not only tests at rational
CM parameters. Part (B) also holds for all six Zagier sporadic sequences by §7 below.
No assertion is made here about the split-prime slope law, Apéry zeta(3)'s carry
defect, or all weight-one/two families.

Inheritance: HYP-9311's elementary factorial-block shift identifies the
carry defect. The [earlier Deuring companion proof](../../05-knowledge/results/overnight9_20260906_jets_deuring.md)
§3 uses a divided reflection and a degree-bounded Wronskian. Its companion
is different from V; the present argument identifies the actual V first.
The corrected near miss is MISTAKE-638: a square factorization of F alone
says nothing about V. The hostile is characteristic p: D(z^p)=0, so zero
derivative does not mean constant. The retained sidecars are the parameter
tangent, the terminal degree, and the carry indices.

## 2. Integer shift and normalization

Write m=(p-1)/2 and H_n=sum_(j=1)^n 1/j, with H_0=0. For k<=m the
factorial-block calculation gives

```
binom(2p+2k,p+k) = 2 binom(2k,k) (1+2p(H_(2k)-H_k))  modulo p^2.
```

Indeed, split each factorial into complete p-blocks and a last partial
block. A complete block of units is (p-1)! modulo p^2 because
sum_(j=1)^(p-1)1/j=0 modulo p. The unit tails give 2pH_(2k)-2pH_k;
the multiples of p give the leading binomial(2,1)=2. Raising to r yields

```
(u_(p+k)-2^r u_k)/p = 2^r u_k * 2r(H_(2k)-H_k)  modulo p.       (1)
```

For m<k<p, both binom(2k,k) and binom(2p+2k,p+k) have p-valuation
exactly one: their factorial valuations are respectively 1-0-0 and
3-1-1 (3p+...<p^2 for p>=7). Their rth powers are divisible by p^r.
Consequently F's and V's coefficients at these indices vanish modulo p
for r>=2. This is a proof of the tail, not a division by a singular
harmonic sum modulo p.

Use z=4^r x and define, in F_p[z],

```
f_r(z) = sum_(k=0)^m ((1/2)_k/k!)^r z^k,
g_r(z) = sum_(k=0)^m ((1/2)_k/k!)^r * 2r(H_(2k)-H_k) z^k.
```

Then F_(r,p)(x)=f_r(4^r x), V_(r,p)(x)=2^r g_r(4^r x).
All denominators in these finite definitions are prime to p.

## 3. Weight one: a Wronskian with its Frobenius term retained

Set f=f_2, g=g_2, theta=zD, and

```
L = theta^2 - z(theta+1/2)^2.
```

Their coefficient recurrences imply the polynomial identities

```
Lf=0,    Lg = -2 theta f + z(2 theta+1)f.                      (2)
```

The top boundary vanishes because m+1/2=0 modulo p, including the
derivative of its squared factor. Equivalently, f log z+g is a formal
logarithmic solution of L. This notation can be eliminated using (2).

For W=fg'-f'g, the ordinary differential form of (2) gives

```
D((1-z)(f^2+zW))=0.                                           (3)
```

The polynomial on the left inside D has degree at most 2m+1=p.
Its constant coefficient is 1, and its value at z=1 is zero. Thus it is
1-z^p, not the constant 1. Dividing by 1-z in the polynomial ring yields

```
f(z)^2 + z(f(z)g'(z)-f'(z)g(z)) = (1-z)^(p-1).               (4)
```

Also f(0)=1 and f(1)=binom(2m,m)=(-1)^m modulo p, by Vandermonde,
since (1/2)_k/k!=(-1)^k binom(m,k) modulo p. At any root alpha of f
in the algebraic closure, alpha!=0,1 and

```
-alpha f'(alpha) g(alpha) = (1-alpha)^(p-1) != 0.              (5)
```

Both simplicity and coprimality follow. No external simplicity theorem
is needed. Identity (4) is also an explicit certificate that the ideal
(f,g) is the unit ideal: modulo (f,g), it contains (1-z)^(p-1), whereas
f(1)!=0, so these polynomials have no common irreducible factor.

## 4. Weight two: differentiate the correct Clausen deformation

Work first in characteristic zero with formal power series. For a=1/4
or 3/4 define

```
H_a(z) = sum_(k>=0) h_(a,k) z^k,       h_(a,k)=(a)_k^2/(k!)^2,
J_a(z) = sum_(k>=0) h_(a,k)*2*(sum_(j=0)^(k-1)1/(a+j)-H_k) z^k.
```

Let f_3^infty and g_3^infty denote the same series as in §2 before
truncation. The precise one-parameter identity to differentiate is

```
[ 2F1(1/4+rho/2, 1/4+rho/2; 1+rho; z) ]^2
 = 3F2(1/2+rho,1/2+rho,1/2+rho; 1+2rho,1+rho; z).            (6)
```

This is classical Clausen (DLMF16.12.2); no new square identity is
claimed. Coefficientwise differentiation at rho=0 gives

```
f_3^infty = H_(1/4)^2,       g_3^infty = H_(1/4) J_(1/4).    (7)
```

To check the tangent exactly: the logarithmic derivative on the right
is 3 sum 1/(1/2+j) - (2+1)H_k = 6(H_(2k)-H_k). On the left,
D_rho H_rho|0=J_(1/4)/2, so the square derivative is HJ.
The denominator slopes 2 and 1 are essential. Equation (6) is not a
claim that every higher parameter derivative equals an integer shift.

Euler's transformation (DLMF15.8.1), applied before differentiation,
has a rho-independent prefactor:

```
2F1(1/4+rho/2,1/4+rho/2;1+rho;z)
 = (1-z)^(1/2) 2F1(3/4+rho/2,3/4+rho/2;1+rho;z).
```

Therefore, exactly over Q[[z]],

```
f_3^infty = (1-z) H_(3/4)^2,
g_3^infty = (1-z) H_(3/4) J_(3/4).                           (8)
```

This derivative of a common identity, with its Euler factor retained,
is the requested carry-compatible deformation. An arbitrary deformation
of a square factorization would not suffice.

## 5. Integral reduction and complete tail accounting

Choose (a,d)=(1/4,(p-1)/4) when p=1 modulo4, and
(a,d)=(3/4,(p-3)/4) when p=3 modulo4. Put

```
B_p(z)=sum_(k=0)^d h_(a,k)z^k,
C_p(z)=sum_(k=0)^d h_(a,k)*2*(sum_(j<k)1/(a+j)-H_k)z^k,
S_p(z)=1 (p=1 mod4),   S_p(z)=1-z (p=3 mod4).
```

All displayed denominators are units at p. For d<k<p, (a)_k has
exactly one factor divisible by p and k! is a unit. Thus h_(a,k)
has valuation two, whereas the harmonic factor has valuation at least
minus one. Both H_a and J_a have coefficients divisible by p on that
whole omitted range. The coefficients of g_3^infty are also integral
through degree p-1: beyond m, the cubic numerator has valuation three
and its harmonic factor has valuation at least minus one.

Reduce (7) or (8) coefficientwise through degree p-1. The left sides
truncate to f_3,g_3; the right sides truncate to S_p B_p^2,S_p B_p C_p.
Their degrees are at most m<p, so these are identities of polynomials,
not merely congruences modulo z^p:

```
f_3 = S_p B_p^2,             g_3 = S_p B_p C_p.               (9)
```

Every root of f_3 is a root of S_p B_p, and hence of g_3. Thus
rad(f_3)|g_3 and part2 follows by the invertible scaling z=64x.
There is no omitted singular exception: for p=3 modulo4, the factor
1-z occurs in g_3 as well as f_3. The case B_p(1)=0, if it occurred,
would not affect this divisibility argument.

## 6. The cube gcd is exactly the radical

A single finite-degree lemma covers both Legendre and the quarter factors.
Suppose a=-d modulo p with 0<=2d<p, and B,C are the truncated equal-parameter
hypergeometric period and harmonic companion through degree d, as in §5.
Their coefficient equations, with L_a=theta^2-z(theta+a)^2, give

```
R=B^2+z(BC'-B'C),      (1-z)R'=2aR,      R(0)=1.
```

The top term of BC'-B'C cancels, so deg R<=2d<p. The coefficient
recurrence (n+1)R_(n+1)=(n+2a)R_n is uniquely solvable through that
range. Since 2a=-2d modulo p, it gives

```
B^2+z(BC'-B'C)=(1-z)^(2d).                              (10)
```

Furthermore B(0)=1 and B(1)=binom(2d,d)!=0 modulo p. Thus the same
root argument as (5) proves B squarefree and gcd(B,C)=1. When d=0,
B=1 and C=0 and the conclusions hold directly. Apply this lemma to
both quarter branches. S_p is squarefree and coprime to B_p, so (9) gives

```
monic gcd(f_3,g_3) = monic(S_p B_p) = rad(f_3).            (11)
```

No additional common root or multiplicity is hidden in the carry defect.
This refines HYP-9330(C), without promoting that hypothesis for other
families. The Legendre case is the same lemma with d=(p-1)/2; the
symmetric square changes how its period and companion enter F and V.

## 7. The same argument proves (B) for all six Zagier families

This extension uses the following precise general lemma. Let integers a,b,c
with p not dividing b or c define an integer sequence u_0=1,u_1=b and

```
n^2 u_n - (a(n-1)^2+a(n-1)+b)u_(n-1)
          +c(n-1)^2u_(n-2)=0,                    n>=2.
```

Assume u_(p+k)=b u_k modulo p for 0<=k<p and that
h(x)=1-ax+cx^2 has two distinct roots. Let F,V be the first block and
its carry defect as in §1, with u_1=b, and set G=V/b modulo p. Then

```
F^2+x(FG'-F'G)=h^(p-1),             gcd(F,h)=1.          (12)
```

In particular F is squarefree and gcd(F,V)=1.

**Proof, retaining the block ends.** Put

```
L=theta^2-x(a theta^2+a theta+b)+cx^2(theta+1)^2,
M=2theta-x(2a theta+a)+2cx^2(theta+1).
```

The recurrence gives LF=0 modulo p. Subtract b times the recurrence
at n=k from the recurrence at n=p+k and divide by p, using the first
block Lucas congruence. For 1<=k<p this gives LG=-MF. At k=1 the
out-of-block u_(p-1) term has coefficient cp^2 before division and
vanishes afterwards. At k=p use the recurrences at n=2p and n=p:
the leading 4p^2u_(2p) and p^2u_p terms also vanish after division.
At degree p+1 the only possible polynomial tail is cp^2G_(p-1),
and the forcing tail is 2cpF_(p-1); both vanish modulo p. Degree0 is
zero on both sides. Thus LG=-MF holds as a whole polynomial identity.

Writing W=FG'-F'G, direct subtraction of these two differential
identities gives

```
D{h(F^2+xW)}=0.
```

The degree is at most2p: degF,degG<=p-1, with the top term of W
cancelled. Hence h(F^2+xW)=Q(x)^p for a polynomial Q in F_p[x]
of degree<=2. It has constant1. At each of the two distinct roots
of h, Q vanishes, so Q=h. Dividing by h proves (12)'s identity.

To justify the singular points, the differential equation LF=0 has
ordinary form

```
(x h F')'+(cx-b)F=0.
```

At a simple root alpha of h, alpha!=0, a zero of F of order
1<=r<=degF<p would contribute leading coefficient
alpha*h'(alpha)*r^2 times its nonzero local leading coefficient.
That is impossible. Thus gcd(F,h)=1. F(0)=1; at every root of F,
(12) now makes -xF'G=h^(p-1) nonzero. This proves the lemma.

For the six Zagier families, the integer recurrence and first-block Lucas
property are known ([Malik–Straub Theorem3.1](https://arminstraub.com/downloads/pub/lucascongruences.pdf),
Tables1–2 and the base-p digit formula (9); also Straub Theorem1.2). The actual
parameters and discriminants are

| Family | (a,b,c) | a^2-4c |
|---|---|---|
| A, Franel | (7,2,-8) | 81 |
| B | (9,3,27) | -27 |
| C | (10,3,9) | 64 |
| D, Apéry zeta(2) | (11,3,-1) | 125 |
| E | (12,4,32) | 16 |
| F | (17,6,72) | 1 |

Every hypothesis of the lemma therefore holds at every prime p>=7.
This promotes exactly the six sporadic cases of HYP-9330(B). It does
not assume Straub's derivative interpolation or the conjectured detector
slope. At p=3 or5 some discriminants, b, or c vanish; this extension
makes no claim there. The omitted discriminant condition is necessary:
Zagier D at p=5 has h=(1+2x)^2 and F=(1+2x)^4, and (B) fails.

## 8. Obstructions to automatic extension

For Apéry zeta(3), CFVZ Theorem2 gives F_p=B_p^2 or
(t^2-34t+1)B_p^2, and Straub Theorem1.3 identifies V_p with a
truncated formal derivative series. These statements alone do not give
a parameter identity whose derivative is that precise series, nor its
complete truncation law. In particular, a differentiated expression
2B_p C_p+D_p F_p only forces vanishing at B_p's roots; it need not
vanish at simple roots of a separate singular prefactor. One must retain
the prefactor's tangent or prove its roots separately. Apéry zeta(3) remains OPEN
here, as do the remaining families in HYP-9330.

The cube proof works because the Euler prefactor in (6) is independent
of rho. That is an exact condition to seek in other families, not an
inference from their weight or a numerical match of singular polynomials.

## 9. Sources, verification, and audit

- [DLMF16.12.2](https://dlmf.nist.gov/16.12.E2): classical Clausen identity.
- [DLMF15.8.1](https://dlmf.nist.gov/15.8.E1): Euler transformation.
- [Straub, arXiv:2301.12248v1](https://arxiv.org/html/2301.12248v1),
  Theorem1.3: cited for the sporadic derivative route, not needed for (1).
- [Caruso–Fürnsinn–Vargas-Montoya–Zudilin, arXiv:2510.23298v1](https://arxiv.org/html/2510.23298v1),
  Theorems2–4: prior square-factor results for Apéry, alternating Domb,
  Almkvist–Zudilin; they do not by themselves identify V_p.

The all-prime conclusions rest on §§2–7. The
[independent audit](../../04-computation/cm_detector_sporadic/audit_bc_20261018l/REPORT.md)
accepts every stated all-prime conclusion, the carry and terminal-index
checks, and the explicit singular-prime exclusions. No new mathematical
error was found; the header was synchronized with the extended scope
before promotion.

Reproduce the parent checks:

```bash
python3 04-computation/cm_detector_sporadic/proof_bc_20261018l.py --bound 199 --out /tmp/proof_bc.json
```

The committed receipt passes 33,049 live gates: 40 exact rational
coefficients of the square and its parameter tangent, every prime7..199
for the full binomial polynomial identities and integer carry/tail
coefficients, and the unchanged existing `mechanism` on Legendre,
cube, Franel and Apéry zeta(3) controls. The latter remains finite evidence.
The parent run also passes under `python3 -O`, with identical output.

The independent standard-library polynomial engine checks 77 odd primes
3..397 for the binomial identities, and a separate exact recurrence
engine checks 450 cases (six Zagier families times 75 primes7..397).
Four sporadic generators are compared with direct binomial sums at every
n=0..100. The existing detector's `mechanism` is exercised on the actual
integer sequences, separately from the derived polynomial formulas.
The theorem deliberately retains p>=7 despite the extra small-prime
binomial checks. All source/receipt hashes are in the audit manifest.
