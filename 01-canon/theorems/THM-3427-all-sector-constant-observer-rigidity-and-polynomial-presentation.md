---
id: THM-3427
title: "All-sector constant-observer rigidity and polynomial presentation"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT / INDEPENDENT
  AUDIT REQUIRED.  Let K have characteristic zero, d>=2, a!=0, and
  P=ax+b+g(x)z^d with nonconstant g.  Every generic fiber-character sector
  has an explicit polynomial differential presentation of dimension
  deg(rad(g)); the wrap sector requires a separate evaluation-at-t splitting
  and punctured-line model.  For 1<=sigma<=d, the constant observer
  [z^(sigma-1)] is generically zero exactly when
  g=c(x-alpha)^e with e>1 and d divides sigma(e-1).  Its exact integral
  annihilator is then ((P-(a alpha+b))^(sigma(e-1)/d)), and is zero otherwise.
  THM-3424 is the sigma=1 specialization.  No full multiroot integral module,
  polynomial mate, new Keller case, or JC(2) conclusion is claimed.  This
  candidate is not a proved dependency before independent audit and explicit
  status promotion.
source: root-2608-jc-all-sector-observer-packet-2026-08-15
depends_on:
  - THM-3419-generic-kummer-response-regular-sector-rank
  - THM-3422-one-root-nonlinear-integral-hamiltonian-response
  - THM-3348-linear-z-generic-puncture-response-and-one-root-valuation
related:
  - THM-3424-nonlinear-monomial-fiber-unit-observer-rigidity
  - THM-3418-one-monomial-nonlinear-fiber-keller-classification
script: 04-computation/jc_all_sector_constant_observer_thm3427.py
output: 05-knowledge/results/jc_all_sector_constant_observer_thm3427.out
script_sha256: 849756abbcab179974864ca7338f53a5f510e8744f4667a6827939cd70857e0f
output_sha256: ebfeb771585beda8ab24af0b90e6909fc1f5f1405d78311527b51043f54d5eea
hash_basis: LF-normalized bytes
---

# THM-3427 -- all-sector constant-observer rigidity and polynomial presentation

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT / INDEPENDENT
AUDIT REQUIRED.**  Nothing in this file is a proved dependency before a
different agent audits it and the status is explicitly promoted.

## 1. Statement and connection contract

Let `K` be a field of characteristic zero, let `d>=2`, and put

```text
P=ax+b+g(x)z^d,                  a in K*,
D=D_P=Jac(P,-),                  C_P=K[x,z]/D(K[x,z]).    (1)
```

Assume `g` is nonconstant.  For `1<=sigma<=d`, let `C_(sigma-1)` be the
target fiber-weight `sigma-1 mod d` sector and define its constant observer

```text
theta_sigma=[z^(sigma-1)] in C_(sigma-1).                 (2)
```

Then

```text
theta_sigma tensor_(K[P]) K(P)=0
iff g=c(x-alpha)^e, e>1, and d divides sigma(e-1).         (3)
```

In the exceptional case put

```text
q_sigma=sigma(e-1)/d,                 beta=a alpha+b.     (4)
```

The exact integral annihilator is

```text
Ann_(K[P])(theta_sigma)=
  ((P-beta)^q_sigma)       in the case (3),
  (0)                      otherwise.                    (5)
```

For `sigma=d`, the divisibility in `(3)` is automatic; its generic observer
vanishes exactly for a single repeated root.  For `sigma=1`, `(3)--(5)` are
exactly the unit-observer theorem THM-3424.  Thus THM-3424 is the first
character of this packet, rather than a separate mechanism.

There is also an explicit presentation of every generic sector.  Normalize
as in Section 2 and set

```text
F=K(t),                  B=F[x,g^(-1)],
S=rad(g),                N=deg(S),

L_sigma(q)=sigma q+(t-x)(sigma(g'/g)q-dq').              (6)
```

Then, including the separately justified wrap character,

```text
C_(sigma-1) tensor_(K[P]) F
 ~=B/L_sigma(B)
 ~=F[x]/L_sigma(S F[x]),                                 (7)

dim_F F[x]/L_sigma(S F[x])=N.                            (8)
```

The wrap sector has the additional punctured-line presentation

```text
B/L_d(B) ~= B/(g partial_x B) ~= H^1_dR(Spec B/F),        (9)
```

up to nonzero scalar factors.

| item | exact content |
|---|---|
| source | the generic THM-3419 sector and its constant monomial observer |
| target | a polynomial differential quotient, observer annihilator, and degree normal form |
| map | `q z^sigma -> L_sigma(q) z^(sigma-1)` |
| preserved | character, generic rank, and polynomial exactness |
| destroyed by localization | the positive annihilator depth in `(5)` |
| required sidecars | root principal parts, the single leading-degree resonance, and evaluation at `x=t` for `sigma=d` |
| sharp hostile | forcing the wrap input to be `qz^d` without splitting off `p(t)` loses weight-zero inputs |

This classifies the displayed observers, not all torsion in a multiroot
integral sector.

## 2. Generic operator for the nonwrap sectors

Replace `P` by `(P-b)/a` and `g` by `g/a`.  This scales `D` by a unit and
changes the `K[P]` generator affinely, so it preserves `(3)--(8)`.  Assume

```text
P=x+g(x)z^d.                                               (10)
```

On the generic fiber `P=t`, THM-3419 gives

```text
A=F[x,g^(-1),z]/(z^d-(t-x)/g)
 =direct_sum_(j=0)^(d-1) B z^j.                           (11)
```

For `1<=sigma<d`, every weight-`sigma` input is uniquely `qz^sigma`, with
`q in B`.  Direct differentiation at fixed `t` gives

```text
D(qz^sigma)
 =z^(sigma-1)[sigma q+z^d(sigma g'q-dgq')]
 =z^(sigma-1)L_sigma(q).                                  (12)
```

Therefore the first isomorphism in `(7)` is immediate for every nonwrap
sector.

## 3. Polynomial representatives and exact intersection

The rational function `g'/g` has only simple poles at the roots of `g`, so

```text
L_sigma(S F[x]) is contained in F[x].                     (13)
```

We claim the natural map

```text
F[x]/L_sigma(SF[x]) -> B/L_sigma(B)                       (14)
```

is an isomorphism.  It suffices to check this after a faithful splitting
extension.  Write locally at a root `alpha` of multiplicity `e`

```text
y=x-alpha,                    g'/g=e/y+(regular).          (15)
```

If `q=c y^(-p)+...`, `p>0`, then the leading term of `L_sigma(q)` is

```text
(t-alpha)(sigma e+dp)c y^(-p-1).                          (16)
```

The coefficient is nonzero.  Hence every pole of order at least two in a
target can be removed successively.  The remaining simple-pole vector is
removed by a polynomial `q` with prescribed values at the distinct roots;
ordinary interpolation supplies it because the simple-pole coefficient is
`(t-alpha)sigma e q(alpha)`.  Thus every class in the right side of `(14)`
has a polynomial representative.

Conversely, if `L_sigma(q)` is polynomial, `(16)` rules out every pole of
`q`.  A nonzero value `q(alpha)` would create an uncancelled simple pole, so
`q` vanishes at every root.  Therefore

```text
L_sigma(B) intersect F[x]=L_sigma(SF[x]).                 (17)
```

This proves `(14)` and the second isomorphism in `(7)`.  Equation `(8)` now
follows from THM-3419's rank-`N` theorem.  Notice that root multiplicities
enter the operator but cancel from its cokernel dimension.

## 4. The wrap input and the evaluation-at-t repair

For `sigma=d`, the target has weight `d-1`, but the actual input has weight
zero and is an arbitrary `p in B`; it is not literally an arbitrary
`qz^d`.  Directly,

```text
D(p)=-dg p' z^(d-1).                                     (18)
```

Evaluation `p |-> p(t)` is defined on `B`, because `g(t)` is a nonzero
element of `F`.  Its kernel is `(x-t)B`.  Hence every `p in B` has a unique
decomposition modulo its constant value

```text
p=p(t)+qz^d,
q=-g[p-p(t)]/(x-t) in B.                                 (19)
```

Since `D(p(t))=0`, equations `(12)` and `(19)` give

```text
D(B)=D(z^dB)=z^(d-1)L_d(B).                              (20)
```

This proves `(7)` for the wrap sector, after which Section 3 applies
unchanged.  It also proves

```text
L_d(B)=-dg partial_x B                                   (21)
```

as subspaces of the target coefficient ring.  The map

```text
[h] |-> [-h dx/(dg)]                                     (22)
```

now gives `(9)`.  In particular, the polynomial classes

```text
[g x^j/S],                       0<=j<N,                  (23)
```

form an `F`-basis of the wrap sector, up to nonzero scalars.  Equations
`(19)--(21)` are the load-bearing wrap sidecar; omitting them is the generic
analogue of omitting the integral wrap quotient in THM-3422.

## 5. One resonance and one basis swap

Put `r=deg(g)`.  If `q` is a polynomial of degree `m` with leading
coefficient `c`, then

```text
lc_degree_m(L_sigma(q))
 =[dm+sigma(1-r)]c.                                      (24)
```

Thus the polynomial matrix of `L_sigma` is degree-lower-triangular and has
at most one zero diagonal, at

```text
m_sigma=sigma(r-1)/d.                                    (25)
```

The domain `SF[x]` begins in degree `N`.  If `m_sigma` is not an integer or
`m_sigma<N`, every diagonal from degree `N` onward is nonzero, and

```text
[1],[x],...,[x^(N-1)]                                    (26)
```

is a basis of the sector.

Suppose `m_sigma` is an integer at least `N`.  Use the nonzero diagonal
columns below `m_sigma` to reduce

```text
L_sigma(Sx^(m_sigma-N))                                  (27)
```

to a polynomial `R_sigma` of degree less than `N`.  It is nonzero: otherwise
the triangular matrix would have cokernel dimension `N+1`, contradicting
`(8)`.  If `rho_sigma=deg(R_sigma)`, a basis is obtained from `(26)` by
replacing `[x^rho_sigma]` with `[x^m_sigma]`.  The scalar normalization of
`R_sigma` depends on choices, but the one-hole basis-swap mechanism does not.

The constant observer vanishes precisely when this resonant defect removes
the constant row.  Sections 6--7 prove that this happens only for a one-root
`g`.  For a multiroot `g`, an accessible degree resonance can still occur,
but it replaces a different low-degree combination; equal generic rank does
not mean equal observer position.

## 6. Nonwrap constant-observer rigidity

Fix `1<=sigma<d`.  By `(7)`, generic vanishing of `theta_sigma` is equivalent
to

```text
L_sigma(q)=1                    for q in SF[x].            (28)
```

The local calculation in Section 3 sharpens: at every root `alpha_i` of
multiplicity `e_i`, `q` must have a simple zero.  Its leading coefficient is
nonzero only if `sigma e_i!=d`; higher zero order cannot supply the constant
right side.  Thus, with `N=deg(rad(g))`,

```text
m=deg(q)>=N.                                              (29)
```

The infinity coefficient `(24)` forces

```text
dm=sigma(r-1).                                            (30)
```

Define

```text
B_sigma(q)=sigma(g'/g)q-dq'.                             (31)
```

Under `(30)`, `B_sigma` has no nonzero polynomial kernel.  Indeed,
`B_sigma(h)=0` gives `h^d=Cg^sigma`, hence `d|sigma r`; but `(30)` gives
`sigma r=sigma mod d`, impossible for `1<=sigma<d`.

Expand `q` at `t=infinity` and write its leading term `t^Lq_0(x)`.  In

```text
sigma q+(t-x)B_sigma(q)=1,                               (32)
```

kernel-freeness forces

```text
L=-1,                         B_sigma(q_0)=1.             (33)
```

Every coefficient remains divisible by `S`; put `n=deg(q_0)>=N`.  After
multiplication by `g`, `(33)` is

```text
sigma g'q_0-dgq_0'=g.                                    (34)
```

If `n>1`, its top coefficient forces `dn=sigma r`, contradicting the same
congruence.  Therefore `n=1`, so

```text
N=1.                                                      (35)
```

The unique root descends to `K`; write `g=c(x-alpha)^e`.  Equation `(30)`
is exactly

```text
e>1,                   d divides sigma(e-1),
m=q_sigma=sigma(e-1)/d.                                  (36)
```

Conversely, assume `(36)`.  Put `u=(x-alpha)/(t-alpha)` and `m=q_sigma`.
Then an explicit solution of `(28)` is

```text
q=sum_(n=1)^m A_n u^n,

A_n=
 d^(n-1)(m-1)! /
 [(m-n)! product_(j=m-n)^(m-1)(dj+sigma)].               (37)
```

The coefficient recurrence is

```text
A_(n+1)[sigma e-d(n+1)]
 +A_n[dn-sigma(e-1)]=0,                                  (38)
```

with `A_1(sigma e-d)=1`; it terminates exactly at `n=m`.
This proves `(3)` for every nonwrap character.  Setting `sigma=1` recovers
the generic theorem in THM-3424.

## 7. Wrap residues and the exact integral annihilator

Under `(22)`, the wrap constant observer is the differential class of
`dx/g`, up to a nonzero scalar.  THM-3348 proves

```text
[dx/g]=0 in H^1_dR(Spec B/F)
iff g=c(x-alpha)^e with e>1.                              (39)
```

This is precisely `(3)` for `sigma=d`.  The separate residue argument is
essential: in `(31)`, `B_d` can have the polynomial kernel generated by `g`,
so the nonwrap two-infinity proof must not be copied across the wrap seam.

It remains to restore integral depth.  In the one-root case THM-3422 gives a
weighted bilateral chain in sector `sigma-1`.  Its unique zero arrow occurs
after the orbit positions

```text
ell=0,1,...,q_sigma-1,                                   (40)
```

because the arrow numerator is

```text
d(ell+1)-sigma(e-1).                                     (41)
```

It is nonzero before the last position and zero there.  Thus `(5)` is the
exact minimal annihilator, not just an upper bound.  In every multiroot or
nonselected case the generic observer is nonzero, so no nonzero polynomial
in `P` can annihilate it integrally.  This proves `(5)`.

For a one-root multiplicity `e`, exactly `gcd(d,e-1)` of the `d` constant
observers are torsion.  Their annihilator exponents are the unequal
intercepts `q_sigma=sigma(e-1)/d` from THM-3422.  For every multiroot `g`, all
`d` displayed constant observers are nontorsion.  This last sentence does
not assert that the full multiroot module is torsion-free.

## 8. Exact replay, hostiles, and boundaries

The exact companion checks the operator coefficients, closed one-root
solutions `(37)`, the wrap evaluation splitting `(19)`, wrap residue vectors,
the two incompatible infinity degrees, minimal integral arrow counts, and
exact `Q(t)` linear systems on a declared low grid.  Its sharp multiroot
hostiles include

```text
(d,sigma;e_1,e_2)
 =(2,1;3,4), (4,2;3,4), (6,2;3,4), (6,3;3,4).            (42)
```

Each reaches the candidate polynomial degree and fails the exact system.
The nonwrap profile grid uses `2<=d<=6`, all nonwrap characters, one to three
roots at the canonical positions `(-1,1,3)`, and multiplicities one through
three.  The wrap residue grid uses one to four roots at `(-2,1,4,7)` and
multiplicities one through six.  These fixed positions scope the finite
evidence only; the proof is root-position independent.
The reproduction command is

```bash
python3 04-computation/jc_all_sector_constant_observer_thm3427.py
```

The boundaries are exact:

- if `g` is constant nonzero, THM-3419 gives `C_P=0` integrally; every
  observer is zero and its annihilator is the unit ideal, outside `(3)--(5)`;
- if `g=0`, the Kummer and logarithmic formulas are undefined, but again
  `C_P=0` because `D=a partial_z` in coordinates `(P,z)`;
- at `d=1`, only the wrap mechanism remains and the generic conclusion is
  THM-3348; the present packet statement assumes `d>=2`;
- characteristic zero is load bearing in the pole coefficients and residue
  primitives;
- generic vanishing in `(3)` is vertical torsion, not a polynomial mate.
  MISTAKE-374 forbids transferring that vanishing back as integral exactness.

No full multiroot integral decomposition, new polynomial mate, new Keller
stratum, or conclusion about the remaining cases of `JC(2)` follows.  QED
for the provisional candidate.
