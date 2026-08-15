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
script_sha256: 9358e4a63d51fbbb0cc074a905b2d5bae7e3a78608599b4b08999413f31cf194
output_sha256: 4055ad0afbcca0f8a0eb3a52700ffbc47a3326650364338e557839d39dfccd34
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

## 5. The resonant defect always removes the top low row

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
columns below `m_sigma` to choose a monic polynomial `q_res in SF[x]` of
degree `m_sigma` such that

```text
R_sigma=L_sigma(q_res),                  deg(R_sigma)<N. (27)
```

Define

```text
B_sigma(q)=sigma(g'/g)q-dq'.                             (28)
```

Normalize `q_res` to have leading `x`-coefficient one.  Expand its
coefficients at `t=infinity`, and let its maximal `t`-order be `kappa`:

```text
q_res=t^kappa q_0(x)+lower t-orders.                     (29)
```

The monic top coefficient has `t`-order zero, so `kappa>=0`.  Every Laurent
coefficient remains divisible by `S`.  In the coefficient of
`t^(kappa+1)` in `(27)`, only the `tB_sigma(q_res)` term occurs.  Hence

```text
leading_t(R_sigma)=B_sigma(q_0),               deg(B_sigma(q_0))<N. (30)
```

Put `j=deg(q_0)`, so `N<=j<=m_sigma`.  Its leading `x`-coefficient under
`B_sigma` is

```text
[sigma r-dj] lc(q_0) x^(j-1).                            (31)
```

This coefficient never vanishes.  For `sigma<d`, equality `dj=sigma r`
together with `dm_sigma=sigma(r-1)` would give
`d(j-m_sigma)=sigma`, impossible.  For `sigma=d`, one has
`m_sigma=r-1`, while a zero in `(31)` would require `j=r>m_sigma`.
Therefore

```text
deg(B_sigma(q_0))=j-1.                                   (32)
```

Combining `(30)` and `(32)` gives `j=N` and the sharper exact defect law

```text
deg(R_sigma)=N-1.                                        (33)
```

Consequently the resonant basis swap is not merely existential.  If `(25)`
is accessible, a basis is

```text
[1],[x],...,[x^(N-2)],[x^m_sigma].                       (34)
```

For `N=1`, the low list is empty and the constant row is removed.  For
`N>1`, the constant remains in `(34)`.  Together with `(26)`, this proves

```text
[1]=0 in F[x]/L_sigma(SF[x])
iff N=1 and m_sigma is an integer at least one.           (35)
```

The unique geometric root in `N=1` is Galois fixed and descends to `K`.
Writing `g=c(x-alpha)^e`, condition `(35)` is exactly `(3)`.  Thus the
all-character observer theorem follows directly from the one-hole normal
form, including the wrap character.  An accessible multiroot resonance
always replaces `[x^(N-1)]`, never the constant.

## 6. Closed selected primitives and the independent wrap audit

The basis argument already proves generic exactness.  There are useful
closed witnesses.  For `1<=sigma<d`, assume `(3)`, put

```text
u=(x-alpha)/(t-alpha),                m=q_sigma.
```

Then `L_sigma(q)=1` for

```text
q=sum_(n=1)^m A_n u^n,

A_n=
 d^(n-1)(m-1)! /
 [(m-n)! product_(j=m-n)^(m-1)(dj+sigma)].               (36)
```

The coefficient recurrence is

```text
A_(n+1)[sigma e-d(n+1)]
 +A_n[dn-sigma(e-1)]=0,                                  (37)
```

with `A_1(sigma e-d)=1`; it terminates exactly at `n=m`.
For `sigma=d`, write the normalized one-root coefficient as
`g=gamma(x-alpha)^e`.  The actual weight-zero input

```text
p=(x-alpha)^(1-e)/[d gamma(e-1)]                          (38)
```

satisfies `D(p)=z^(d-1)`.  Independently, under `(22)` the wrap observer is
the differential class of `dx/g`, up to a nonzero scalar.  THM-3348's
residue theorem gives

```text
[dx/g]=0
iff g=c(x-alpha)^e with e>1,                              (39)
```

which audits `(35)` at the wrap seam.  The separate wrap construction is
essential: `B_d` can have the polynomial kernel generated by `g`, even
though the degree bound `j<=m_d=r-1` keeps that kernel out of `(31)`.
Setting `sigma=1` in `(35)--(37)` recovers THM-3424 as the first-character
specialization, without copying its proof.

## 7. Exact integral annihilator

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
the two incompatible infinity degrees, `74` exact resonant defect reductions,
minimal integral arrow counts, and exact `Q(t)` linear systems on a declared
low grid.  Its sharp multiroot hostiles include

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
