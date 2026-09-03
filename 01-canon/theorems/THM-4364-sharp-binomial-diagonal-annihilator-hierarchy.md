---
id: THM-4364
title: "Sharp binomial diagonal annihilator hierarchy"
status: >
  PROVED FINITE DEPTH-MODULE THEOREM RELATIVE TO THM-4308 + VERIFIED-EXACT
  + INDEPENDENTLY AUDITED. For ell>=2, m>=ceil(ell/2), and q,d>=0, the
  primitive alternating q-simplex functional on the diagonal r=2n-ell
  annihilates pi_m(P_d) exactly when q<ceil(ell/3), m+q>=ell, and
  d<=m+q-ell. Each first-excluded order, row, and depth boundary has an
  explicit unit hostile. The row-eight triangular and row-nine/ten
  tetrahedral functionals are specializations. No bracket, all-row
  membership, seam-entry, Keller-pair, JC(2), or DC(2) conclusion is asserted.
source: root + diagonal-hierarchy scout + independent boundary referee / next-sharp session, 2026-09-02
depends_on:
  - THM-4308-source-normal-bracket-hasse-truncation-through-row-eight
related:
  - THM-4358-source-normal-s4339-row-ten-delayed-depth-extinction
  - THM-4361-source-normal-beta-zero-row-ten-joint-depth-extinction
  - THM-4362-tetrahedral-diagonal-annihilator-and-sharp-depth-threshold
  - THM-3417-exchangeable-krawtchouk-shell-inversion-and-radial-high-minor-converse
mistake_firewall:
  - MISTAKE-522
primary_script: 04-computation/binomial_diagonal_annihilator_hierarchy_thm4364.py
primary_output: 05-knowledge/results/binomial_diagonal_annihilator_hierarchy_thm4364.out
primary_script_sha256: 67b6605bf1a41cafb0b301f18a277c3aba5109bed6848253604aadb52e90abf5
primary_output_sha256: ce6713e23325e350dedc2a0cae6c7c6f9c52163d66d4a39851a237cda90fa658
independent_referee_script: 04-computation/binomial_diagonal_annihilator_hierarchy_independent_referee_thm4364.py
independent_referee_output: 05-knowledge/results/binomial_diagonal_annihilator_hierarchy_independent_referee_thm4364.out
independent_referee_script_sha256: 145d4379255f1c363a978e21a43f3efe9ab7d5610491479545428c02cd5788cd
independent_referee_output_sha256: f5090ba77e1a836b7c7af43037679a968538a118c83b3d5a0dfd9ec4d1ddcecb
hash_basis: raw LF bytes
audit: >
  PASS WITH PARITY REPAIR. The primary performs 560,623 exact checks. An
  import-free boundary referee derives both closed evaluation formulas, the
  global iff, and all three unit boundary families independently, with
  7,528,886 grouped exact assertions. Normal/-O/-I/frozen LF streams match.
---

# THM-4364 -- Sharp binomial diagonal annihilator hierarchy

**PROVED FINITE DEPTH-MODULE THEOREM RELATIVE TO THM-4308 + VERIFIED-EXACT
+ INDEPENDENTLY AUDITED. THIS IS AN EXACT THEOREM ABOUT `pi_m(P_d)`. IT DOES
NOT ASSERT A BRACKET SOLUTION, AN ALL-ROW LIFT, SEAM ENTRY, A KELLER PAIR,
`JC(2)`, OR `DC(2)`.**

## 1. Statement and inheritance

Work in THM-4308's source-normal characteristic-zero setting. Write

```text
h_(n,r)=[x^r t^n]h,
pi_m = projection retaining exactly rows 0<=n<=m,
P_d=span{x^a u^b p^c y^e: a,b,c,e>=0, a+b<=d}.            (1)
```

For integers

```text
ell>=2,       s=ceil(ell/2),       m>=s,       q,d>=0,
rho=ceil(ell/3),
```

define the primitive functional

```text
L_(m,ell,q)(h)
 =sum_(n=s)^m (-1)^(n-s) C(m+q-n,q) h_(n,2n-ell).         (2)
```

Every binomial coefficient in `(2)` is ordinary: its upper argument is at
least `q`. The final coefficient is `+/-C(q,q)=+/-1`. The exact theorem is

```text
L_(m,ell,q) annihilates pi_m(P_d)
 iff
q<rho,              m+q>=ell,              d<=m+q-ell.   (3)
```

In particular,

```text
L_(m,ell,q) annihilates pi_m(P_0)
 iff q<ceil(ell/3) and m+q>=ell.                           (4)
```

The inheritance pass was:

- closest proved mechanism: THM-4308's exact monomial spanning and expansion
  laws, specialized by THM-4362 at `(ell,q)=(10,3)`;
- canonical hostile: a finite left-kernel row need not persist at the next
  depth, order, or endpoint;
- corrected near miss: the general hostile parity is the terminal-row parity
  `(ell+d) mod 2`, not depth parity alone;
- least-used sidecar: the interval of all depth-zero packet degrees, rather
  than only its smallest member.

The live concept board was

```text
diagonal intercept | packet degree | starting row | terminal overflow
simplex order | zero collar | Euclidean endpoint | unit hostile.           (5)
```

## 2. One monomial and two exact evaluations

THM-4308 gives the expansion

```text
g=x^a u^b p^c y^e
 =x^(a+2b+e)t^(b+c+2e)(1+x^2t)^(c+e).                    (6)
```

Put

```text
N=c+e,       n0=b+c+2e,       r0=a+2b+e,       n1=n0+N.
```

Its `k`-th term is

```text
C(N,k)x^(r0+2k)t^(n0+k),                    0<=k<=N.      (7)
```

The packet meets the diagonal in `(2)` if and only if

```text
r0+2k=2(n0+k)-ell  iff  a=2c+3e-ell.                     (8)
```

This condition is independent of `k`, so the whole packet lies on the
diagonal or all of it misses. Under `(8)`,

```text
N>=rho,              n0>=s,              n1=a+b+ell.     (9)
```

The first two bounds follow from
`3(c+e)>=2c+3e>=ell` and `2(c+2e)>=2c+3e>=ell`. The last
identity says that terminal row is exactly depth plus diagonal intercept.

There are two useful closed forms. A packet missing `(8)` or beginning after
row `m` contributes zero. Otherwise, if `n1<=m+q`, coefficient extraction
over the complete packet gives

```text
L(g)=0,                                                if N>q,
L(g)=(-1)^(n0-s) C(m+q-n1,q-N),                       if N<=q.  (10)
```

Indeed, with `R=m+q-n0>=N`, the unsigned identity is

```text
sum_(k=0)^N (-1)^k C(N,k)C(R-k,q)
 =[z^q]z^N(1+z)^(R-N)=C(R-N,q-N).                       (11)
```

No negative-top convention occurs. If instead `n1>m+q`, put

```text
delta=n1-(m+q)>0,             K=m-n0=N-q-delta>=0.
```

The retained partial packet has the exact value

```text
L(g)=(-1)^(n0-s+K) C(N-q-1,delta-1).                    (12)
```

Reversing the finite sum proves `(12)`:

```text
sum_(k=0)^K (-1)^k C(N,k)C(q+K-k,q)
 =(-1)^K[x^K](1+x)^(N-q-1)
 =(-1)^K C(N-q-1,delta-1).                              (13)
```

These formulas expose separately the three pieces of `(3)`: packet order,
retained endpoint, and allowed depth.

## 3. Why both parameter bounds are necessary

For `P_0`, equations `(8)` and `a=b=0` give `2c+3e=ell`. Packet degrees are
in bijection with the full integer interval

```text
rho<=N<=floor(ell/2),
c=3N-ell,                  e=ell-2N,                    (14)
```

and every such packet satisfies

```text
n0=ell-N,                  n1=ell.                       (15)
```

If `q<rho` and `m+q>=ell`, every packet has `N>q` and `(10)` vanishes.

Conversely, suppose first that `m+q<ell`. Choose

```text
N=max(rho,ell-m).
```

Because `m>=s`, this belongs to `(14)` and begins no later than row `m`.
Here `N>q`, `delta=ell-(m+q)>0`, and `(12)` has positive magnitude

```text
C(N-q-1,delta-1).                                       (16)
```

If instead `m+q>=ell` but `q>=rho`, use the same `N`. Then `N<=q`, and
`(10)` has positive magnitude

```text
C(m+q-ell,q-N).                                         (17)
```

Thus failure of either parameter bound is already visible in `P_0`. This
proves `(4)` and the necessity half of `(3)` for those two gates.

## 4. The depth threshold and corrected parity

Assume the two parameter inequalities in `(3)` and let
`a+b<=d<=m+q-ell`. Every diagonal packet then satisfies

```text
n1=a+b+ell<=d+ell<=m+q,                N>=rho>q.          (18)
```

Equation `(10)` kills it. Linearity and `(1)` prove the positive direction.

At the first excluded depth set

```text
d0=m+q-ell+1,                 S=ell+d0=m+q+1,
epsilon=S mod 2,
a=d0,       b=0,       e=epsilon,       c=(S-3epsilon)/2. (19)
```

All exponents are nonnegative, `(8)` holds, `a+b=d0`, and

```text
n1=S=m+q+1.
```

This is the one-row-overrun case `delta=1`; `(12)` gives the unit

```text
L_(m,ell,q)(g)=(-1)^(m-s).                               (20)
```

The same monomial belongs to every later `P_d`, so nesting proves failure for
all `d>=d0`. This completes `(3)`. The parity in `(19)` is essential: using
`d0 mod 2` works on THM-4362's even intercept `ell=10` but fails for odd
`ell`.

The other two first-excluded boundaries also have units. For
`q=rho,m=ell-rho`, use `(14)` with `N=rho`. For admissible `q<rho` but
`m=ell-q-1`, use the same minimal packet. In both cases the value is
`(-1)^(m-s)`. Hence none of the three inequalities in `(3)` can be relaxed
at its first boundary.

## 5. Existing consumers and figurate continuation

The specialization

```text
(ell,q,m,d)=(8,2,8,2)
```

gives the triangular row-eight `P_2` functional from THM-4308:

```text
15h_(4,0)-10h_(5,2)+6h_(6,4)-3h_(7,6)+h_(8,8).          (21)
```

Its first-depth hostile is `x^3p^4y in P_3`, with value `+1`. The two
specializations

```text
(ell,q,m,d)=(10,3,9,2),        (10,3,10,3)
```

give respectively

```text
35h_(5,0)-20h_(6,2)+10h_(7,4)-4h_(8,6)+h_(9,8),         (22)

56h_(5,0)-35h_(6,2)+20h_(7,4)-10h_(8,6)
 +4h_(9,8)-h_(10,10).                                   (23)
```

These are THM-4358's `H_A` and THM-4361's `H_C`; their first-depth hostiles
are `x^3p^5y` with value `+1` and `x^4p^7` with value `-1`. The hierarchy
proves only that these are universal depth-module hyperplanes. The cited
theorems separately evaluate bracket-selected jets against them.

For `q=2`, polynomial binomial continuation is exactly the integer extension
of triangular numbers:

```text
C(z,2)=T(z-1),             T(-N)=T(N-1).                 (24)
```

The two zero collar rows after `m` have weights

```text
C(1,2)=T(0)=0,             C(0,2)=T(-1)=0,
```

while the first row beyond them has

```text
C(-1,2)=T(-2)=T(1)=1.                                   (25)
```

Thus reflection of the negative triangular address is exactly the unit that
makes the threshold sharp. For general `q`, the same mechanism is the
`q`-row zero collar followed by `C(-1,q)=(-1)^q`.

The consumer-complete natural address for one packet is not a bare ordinal
but at least

```text
(ell,N,a+b,n0): intercept, packet order, depth, start row.                (26)
```

The endpoint follows as `a+b+ell`. Forgetting `N` destroys the order test;
forgetting depth destroys the endpoint test. The resemblance to signed
Krawtchouk transforms in THM-3417 stops at finite differences: no tournament
or LRC statement follows.

## 6. Audit and scope

The primary performs 560,623 exact checks of the finite-difference identity,
partial boundaries, visible source generators, `P_0` phase diagram, unit
hostiles, and the three named stencils. The import-free referee independently
checks both closed forms, 11,031 parameter triples, 5,022,057 direct generator
evaluations, 2,118,784 admissible projected columns, and all three unit
boundaries through `ell=120`, for 7,528,886 grouped exact assertions. Normal,
optimized, isolated, and frozen LF streams agree.

The proof after the spanning law is integral, and its sharp hostiles are
units. The identification of `(1)` with the geometric depth module is imported
only in THM-4308's characteristic-zero setting. This theorem does not provide
row existence, a bracket solution, infinite-depth membership, polynomial
termination, seam entry, a Keller pair, `JC(2)`, or `DC(2)`.

Reproduce from the repository root:

```text
python3 -B 04-computation/binomial_diagonal_annihilator_hierarchy_thm4364.py
python3 -B -O 04-computation/binomial_diagonal_annihilator_hierarchy_thm4364.py
python3 -B 04-computation/binomial_diagonal_annihilator_hierarchy_independent_referee_thm4364.py
python3 -B -O 04-computation/binomial_diagonal_annihilator_hierarchy_independent_referee_thm4364.py
```
