---
id: THM-3730
title: "Positive distinct two-cube support exponent, abscissa, and critical divergence"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Positive distinct
  two-cube representations have counting asymptotic
  kappa X^(2/3)+O(X^(1/3)).  Both the indexed and deduplicated Dirichlet
  series have abscissa 2/3 and diverge at the boundary.  Every split of every
  prime p=2 mod 3 into two positive unequal summands gives a singleton
  representation fibre; this supplies the deduplicated critical mass.
  No support asymptotic, support residue, or general fibre classification is
  claimed.
source: root + lrc14-cover-defect-bridge / 2026-08-22
audit: >
  PASS.  The support probe checks complete THM-463 divisor fibres for all
  6,825 represented values m<=2,000,000 (6,889 indexed pairs).  The boundary
  probe checks valuations and split injectivity for all 387,587 constructed
  values through p<=5000, and complete divisor fibres for all 308 targets
  through p<=101.  The independent probe checks complete coordinate fibres
  for all 12,494 targets through p<=800 and cross-prime injectivity for all
  387,587 values.  It also checks primitive reduction and the non-cube-free
  hostile 71*7^4=16^3+55^3.  The earlier cube-free sieve is verified but is
  not load-bearing.  Normal and optimized streams byte-match every frozen
  transcript.
depends_on:
  - THM-463-two-cube-representations-are-a-divisor-property-on-the-split-axis
related:
  - THM-2000-support-harmonic-abel-dini-figurate-surface
script: 04-computation/two_cube_support_abscissa_thm3730.py
output: 05-knowledge/results/two_cube_support_abscissa_thm3730.out
script_sha256: 04fb56c3b5145d640eb9af7b91610eeae6d5d2c9489b555416d484954524555c
output_sha256: e5e886bbb40dc51e0d1b38922391180f82d598dbc940dd33aa9ade5f26deb12b
boundary_script: 04-computation/two_cube_inert_prime_boundary_thm3730.py
boundary_output: 05-knowledge/results/two_cube_inert_prime_boundary_thm3730.out
boundary_script_sha256: b79a44eaa336e56b6f5482fa58395c67709450627cabc4dab08804fa41e80831
boundary_output_sha256: cf73df015dcefc8c7848cabe3bf793e80c6e47c922646a8dc47650ca8a4108ea
independent_audit_script: 04-computation/two_cube_inert_prime_boundary_independent_audit_thm3730.py
independent_audit_output: 05-knowledge/results/two_cube_inert_prime_boundary_independent_audit_thm3730.out
independent_audit_script_sha256: 533cce2de4b32d4da664cef62455a5f6e8ab0c5cd28c5cc6cfa309d39da36886
independent_audit_output_sha256: 238f6a20a3bd2fc864e1a5111978d89fb04e93c3074edc42fb6381d1983c2dc5
hash_basis: raw LF bytes
---

# THM-3730 -- positive distinct two-cube support diverges at its critical boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

Define

```text
r_+(m)=#{(x,y) in Z_(>0)^2: x<y, x^3+y^3=m},
C_+={m:r_+(m)>0},
R(X)=sum_(m<=X)r_+(m),
C(X)=#(C_+ intersect [1,X]).
```

Put

```text
kappa=Gamma(4/3)^2/(2 Gamma(5/3)).
```

Then

```text
R(X)=kappa X^(2/3)+O(X^(1/3)).                         (1)
```

The indexed representation series

```text
F(s)=sum_(0<x<y)(x^3+y^3)^(-s)
    =sum_m r_+(m)m^(-s)
```

has abscissa of absolute convergence `2/3`, diverges at `s=2/3`, and

```text
lim_(sigma downarrow 2/3) (sigma-2/3)F(sigma)
 =Gamma(4/3)^2/(3 Gamma(5/3)).                         (2)
```

The deduplicated support series

```text
D(s)=sum_(m in C_+)m^(-s)
```

also has abscissa `2/3` and

```text
D(2/3)=infinity.                                       (3)
```

Consequently `C_+` has counting exponent `2/3`, natural and logarithmic
density zero, and finite harmonic mass `D(1)`.

## 2. Lattice asymptotic and indexed series

Scale the bounded region

```text
Omega={(u,v):0<u<v, u^3+v^3<=1}.
```

Its boundary is rectifiable, so the unit-square boundary comparison gives

```text
#(t Omega intersect Z^2)=area(Omega)t^2+O(t).          (4)
```

Coordinate exchange halves the first-quadrant superellipse away from its
measure-zero diagonal.  With `z=u^3`,

```text
area(Omega)
 =1/2 integral_0^1 (1-u^3)^(1/3)du
 =Gamma(4/3)^2/(2 Gamma(5/3))
 =kappa.                                               (5)
```

Taking `t=X^(1/3)` in (4) proves (1).  Abel--Stieltjes partial summation
then gives the indexed abscissa and boundary divergence.  The leading term
contributes `kappa*(2/3)` to the real-axis critical coefficient, while the
`O(X^(1/3))` error is integrable after multiplication by
`sigma-2/3`.  This proves (2).

## 3. Support exponent and abscissa

THM-463 sends every representation injectively to its good divisor
`d=x+y`.  Therefore

```text
r_+(m)<=tau(m).                                        (6)
```

For every `epsilon>0`, the elementary divisor bound gives a constant
`A_epsilon` such that `tau(m)<=A_epsilon m^epsilon`.  From (1) and (6),

```text
C(X)>=R(X)/(A_epsilon X^epsilon)
     >>_epsilon X^(2/3-epsilon),
C(X)<=R(X)<<X^(2/3).                                   (7)
```

Hence

```text
lim_(X->infinity) log C(X)/log X=2/3.                  (8)
```

For `sigma>2/3`, `D(sigma)<=F(sigma)<infinity`.  If
`0<=sigma<2/3`, choose `epsilon<2/3-sigma`; the support partial sum
through `X` is at least `X^(-sigma)C(X)` and is unbounded by (7).
For `sigma<0`, each term is at least one.  Thus the support abscissa is
`2/3`.  The divisor bound alone does not decide its endpoint.

## 4. Every inert-prime split is a singleton fibre

Let `p>=5` be prime with `p=2 mod 3`, choose `1<=x<p/2`, set
`y=p-x`, and write

```text
q_p(x)=x^2-xy+y^2=3x^2-3px+p^2,
m_p(x)=p q_p(x)=x^3+y^3.                              (9)
```

Modulo `p`,

```text
q_p(x)=3x^2 !=0,             so v_p(m_p(x))=1.         (10)
```

Consider any competing positive distinct representation
`m_p(x)=u^3+v^3`.  Put `g=gcd(u,v)`, `U=u/g`, and `V=v/g`.
Since `g^3|m_p(x)`, (10) implies `p` does not divide `g`; the primitive
sum `U^3+V^3=m_p(x)/g^3` still contains `p`.  THM-463's split lemma says
that an inert prime cannot divide the Eisenstein cofactor of a primitive
representation.  It follows that

```text
p|(U+V),                  hence p|(u+v).               (11)
```

Let `d=u+v` and `e=v-u`.  Positivity and distinctness give
`0<e<d`, so

```text
4m_p(x)=d(d^2+3e^2)>d^3.                              (12)
```

The original cofactor in (9) satisfies `q_p(x)<p^2`, hence
`m_p(x)<p^3`.  Equations (11)--(12) yield

```text
0<d<(4m_p(x))^(1/3)<4^(1/3)p<2p.
```

The only positive multiple of `p` in this interval is `d=p`.
THM-463's good-divisor injection now forces

```text
{u,v}={x,p-x},                 r_+(m_p(x))=1.          (13)
```

No cube-free hypothesis is used.  The first hostile is

```text
p=71, x=16, y=55, q_p(x)=2401=7^4,
```

and (13) still holds.

For fixed `p`, `q_p(x)` is strictly decreasing on `1<=x<p/2`, so its
`(p-1)/2` values are distinct.  Values from two different inert primes are
also distinct: equality would give two representations with different
pair-sums, contradicting (13).

## 5. Critical support divergence

The reciprocal sum of primes `p=2 mod 3` diverges by an elementary Euler
product.  Let `chi` be the nonprincipal real character modulo three.
The grouped series

```text
L(1,chi)=sum_(k>=0)(1/(3k+1)-1/(3k+2))
```

is finite and strictly positive, and bounded partial sums of `chi` give
`L(s,chi)->L(1,chi)` for real `s` decreasing to one.  For `s>1`,

```text
zeta(s)/L(s,chi)
 =(1-3^(-s))^(-1)
   product_(p=2 mod 3)(1+p^(-s))/(1-p^(-s)).           (14)
```

If `sum_(p=2 mod 3)1/p` converged, then
`log((1+z)/(1-z))<=4z` for `0<=z<=1/2` would uniformly bound the product
in (14), contradicting `zeta(s)->infinity` and
`L(s,chi)->L(1,chi)>0`.

Every value in (13) satisfies `m_p(x)<p^3`, so disjointness and positivity
give

```text
D(2/3)
 >=sum_(p=2 mod 3,p>=5) ((p-1)/2)p^(-2)
 >=(2/5)sum_(p=2 mod 3,p>=5)1/p
 =infinity.                                            (15)
```

This proves (3).  More explicitly, the left side truncated at `X` is at
least the first sum in (15) over `p<=X^(1/3)`.

## 6. Collision tax, density, and scope

For `Re(s)>2/3`, THM-2000's exact support/multiplicity split is

```text
F(s)=D(s)+sum_m (r_+(m)-1)_+ m^(-s).                  (16)
```

The first nonzero tax is
`1729=1^3+12^3=9^3+10^3`.  The theorem determines the indexed asymptotic,
both abscissae, and both boundary divergences.  It does **not** determine an
asymptotic for `C(X)`, a pole or critical coefficient for `D`, the
collision-tax asymptotic, or general taxicab fibres.

The bound `C(X)=O(X^(2/3))` gives natural density zero.
Since `D(1)<infinity`, its logarithmic density is also zero.

## 7. Reproduction

```bash
python3 -B 04-computation/two_cube_support_abscissa_thm3730.py
python3 -B -O 04-computation/two_cube_support_abscissa_thm3730.py
python3 -B 04-computation/two_cube_inert_prime_boundary_thm3730.py
python3 -B -O 04-computation/two_cube_inert_prime_boundary_thm3730.py
python3 -B 04-computation/two_cube_inert_prime_boundary_independent_audit_thm3730.py
python3 -B -O 04-computation/two_cube_inert_prime_boundary_independent_audit_thm3730.py
```

Each normal/optimized pair must agree byte for byte with its frozen
transcript.  **QED.**
