---
id: THM-4340
title: "U-zero repeated-cusp planar Jacobian extinction"
status: >
  PROVED RELATIVE TO THM-4230/4327 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. In the inherited reduced (2,3), exact-weight-twelve seam, subject
  to K=2848/45-(7/6)Delta, the entire U=0, WZ!=0 stratum is extinct. The two
  residual repeated faces have exact completed local equations
  x^2=z^5-psi(tz) and x^2=X^3-psi(tX). Every possible critical-value order
  gives either a positive-genus tail of strictly positive Keller-form order,
  a rational tail, or a horizontally persistent cusp with regular
  simultaneous normalization. The Lambda=0 A23 point is a distinct owner
  address and remains covered by THM-4327. The walls U=0,WZ=0, seam entry,
  JC(2), and DC(2) remain open.
source: root + u-zero conductor/referee agents / planar-Jacobian next-sharp session, 2026-09-02
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
related:
  - THM-4289-a23-blowdown-observer-kahler-dualizing-quotient
  - THM-4297-general-lambda-zero-central-and-tail-keller-extinction
  - THM-4339-clean-interior-cubic-edge-planar-jacobian-extinction
mistake_firewall:
  - MISTAKE-522
  - MISTAKE-531
  - MISTAKE-540
primary_script: 04-computation/jc2_m12_u0_repeated_cusp_extinction_thm4340.py
primary_output: 05-knowledge/results/jc2_m12_u0_repeated_cusp_extinction_thm4340.out
primary_script_sha256: 9f0ddd473b70ac835f1c95ed1ccdb6bb864e94aa7e952b4b2506030119ea357e
primary_output_sha256: 4cf0e26bb59bb92cda68e0223d0c2e5bdce5ff66be06f42ea26bdc2a4f159953
independent_audit_script: 04-computation/jc2_m12_u0_repeated_cusp_extinction_independent_audit_thm4340.py
independent_audit_output: 05-knowledge/results/jc2_m12_u0_repeated_cusp_extinction_independent_audit_thm4340.out
independent_audit_script_sha256: 8beb9aa68b9a2f55658ad19ac5b02189a0f3402186eb6b02c821ae9479d6b228
independent_audit_output_sha256: 938625e0264d1d50c29767190c5dab661f510096e7f1b5bd605f24eb1d634d24
hash_basis: raw LF bytes
audit: >
  PASS AFTER FOUR PRESENTATION REPAIRS. The primary certificate reconstructs
  every active source monomial, the exact cusp factorizations, Hasse
  hostiles, conductor/Jacobian quotients, tail scales, differential orders,
  normalized genus ledgers, and a cancellation counterexample. The
  import-independent SymPy audit reconstructs the full sixteen-term source,
  differentiates at fixed original coordinates, solves the hostile critical
  sections, checks honest DVR base changes, and repeats every tail ledger.
  The referee required the normalized-genus correction, the inherited
  K/Delta relation, and explicit separation from the Lambda=0 address.
  Normal and optimized streams byte-match both frozen outputs.
---

# THM-4340 -- U-zero repeated-cusp planar Jacobian extinction

**PROVED RELATIVE TO THM-4230/4327 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED. IN THE INHERITED EXACT-WEIGHT-TWELVE SEAM, `U=0,WZ!=0` IS
EXTINCT. THE WALL `WZ=0`, SEAM ENTRY, `JC(2)`, AND `DC(2)` REMAIN OPEN.**

## 1. Statement and inheritance

Work over an algebraically closed field of characteristic zero in the
inherited reduced `(2,3)`, exact-weight-twelve seam. Retain the complete
sixteen-term residual source of THM-4230, including

```text
e=-1376/135,                 K=2848/45-(7/6)Delta.        (1)
```

Then no nonautomorphic planar Keller pair lies on

```text
U=0,                         WZ!=0.                       (2)
```

All other lower coefficients and `Lambda=W+Z` are arbitrary subject to the
inherited seam hypotheses and `(1)`. This is a theorem inside that seam; it
does not prove entry into the seam or `JC(2)`.

THM-4327 already proves every reduced-face case of `(2)` except exactly:

```text
A: u!=0, E5=alpha^2-4Wu=0;                               (3)
E: u=alpha=Delta=0, E3=eta^2-4eW=0.                     (4)
```

The repeated roots are interior because

```text
S0=-alpha/(2W)!=0,                 V0=-eta/(2W)!=0.      (5)
```

The inheritance pass was:

- closest mechanism: THM-4327's positive good-differential order followed
  by proper-flat degree conservation;
- canonical hostile: THM-4289's warning that regular dualizing does not
  imply ambient Kahler extension;
- corrected near miss: the normalized genus after a persistent
  `A_(r-1)` is `g_Pick-floor(r/2)`, not the raw Pick genus;
- least-used sidecar: the monomial buffer multiplying the exact Keller
  residue at the cusp.

The live board was

```text
critical section | self-similar variable | conductor/Jacobian quotient
residue buffer | persistent cusp | stable tail | conserved degree.       (6)
```

## 2. The two exact one-series normal forms

Start from the inherited exact source equation

```text
F_Q=(s^2-p)(1-QH(p,sp))-Q*s^2/2.                         (7)
```

### 2.1 The `(2,5)` face

Put

```text
Q=sigma^60, s=S, p=sigma^-12 P, t=sigma^12,
z=P^-1, rho=tz, G=sigma^12 F_Q, Phi=z^6G.                (8)
```

For the coefficient of `p^i(sp)^j`, collect by `n=i+j` and define

```text
f_n(S)=sum_(i+j=n)c_(i,j)S^j,
A(S,rho)=sum_(n=1)^5 f_n(S)rho^(5-n).                    (9)
```

Direct expansion of all fifteen active `U=0` rows plus the base terms gives

```text
Phi=(1-S^2rho)(A(S,rho)-z^5)-S^2rho^6/2.                (10)
```

The first factor is a unit at `(S0,z,t)=(S0,0,0)`. After dividing by it,
the `S`-dependent term is

```text
B(S,rho)=A(S,rho)-S^2rho^6/[2(1-S^2rho)].                (11)
```

On `(3)`, `B(S,0)=W(S-S0)^2` and `B_SS(S0,0)=2W!=0`.
The parametrized formal Morse lemma over `k[[rho]]` supplies a unique
critical section `S_c(rho)`. With

```text
psi(rho)=B(S_c(rho),rho),                                (12)
```

a `rho`-preserving formal coordinate change gives the exact completed form

```text
x^2=z^5-psi(tz),                         psi(0)=0.        (13)
```

### 2.2 The `(2,3)` face

Use

```text
Q=sigma^12, s=sigma^4S, p=sigma^-4P, X=P^-1,
V=SP, t=sigma^4, rho=tX, G=sigma^4F_Q, Phi=X^4G.         (14)
```

Case E deletes precisely the rows with `p`-exponent greater than three.
The twelve surviving source rows give, directly,

```text
Phi=(1-V^2rho^3)(A(V,rho)-X^3)-V^2rho^6/2.              (15)
```

The `rho^3` coefficient retains the `K V^2`, `zeta_3 V^3`, and `Z V^4`
terms. On `(4)`, the divided transverse function satisfies
`B(V,0)=W(V-V0)^2`. The same Morse argument yields

```text
x^2=X^3-psi(tX),                         psi(0)=0.        (16)
```

These are completed-local automorphisms with unit Jacobian. They preserve
the base and root coordinate; no ramified root cover is being substituted
for an exceptional function field.

## 3. Cancellation-complete self-similar cusp lemma

Both residual germs have the form

```text
x^2=z^m-psi(tz),                         m in {5,3}.      (17)
```

If `psi=0`, or `r=ord_rho(psi)>=m`, the right side is `z^m` times a unit.
After extracting its square root, `z=w^2,x=w^m` is a finite regular
simultaneous normalization. The `A_(m-1)` cusp persists horizontally and no
positive-genus vertical component appears.

Suppose `1<=r<m` and `psi(rho)=rho^r u(rho)`, `u(0)!=0`. Then

```text
x^2=z^r[z^(m-r)-t^r u(tz)].                              (18)
```

An honest dominating DVR base change and integer weighted blowup are

```text
t=tau^[2(m-r)],       z=tau^[2r]Z,       x=tau^[mr]Y.    (19)
```

Write `r=2q+epsilon`, `epsilon in {0,1}`, and divide `Y` by `Z^q`.
The reduced normalized tail is

```text
Y_0^2=Z^epsilon[Z^(m-r)-u(0)].                           (20)
```

Its branch polynomial is squarefree of odd degree `m-r+epsilon`. The affine
axis is smooth, its projective closure has one smooth point at infinity,
and that point is its sole attachment. Hence it creates no graph cycle and

```text
g_tail=delta(A_(m-1))-delta(A_(r-1)),
delta(A_j)=floor((j+1)/2).                               (21)
```

The complete local list is

| cusp | `r` | normalized tail | genus | persistent delta |
|---|---:|---|---:|---:|
| `(2,5)` | 1 | `Y^2=Z(Z^4-c)` | 2 | 0 |
| `(2,5)` | 2 | `Y_1^2=Z^3-c` | 1 | 1 |
| `(2,5)` | 3 | `Y_1^2=Z(Z^2-c)` | 1 | 1 |
| `(2,5)` | 4 | `Y_2^2=Z-c` | 0 | 2 |
| `(2,3)` | 1 | `Y^2=Z(Z^2-c)` | 1 | 0 |
| `(2,3)` | 2 | `Y_1^2=Z-c` | 0 | 1 |

The valuation direction used here is

```text
v(z^m-psi(tz)) >= min(mv(z),r(v(t)+v(z))).               (22)
```

Equality can fail on the balanced face, so the reverse inequality is
invalid. For `m=3`, `psi(rho)=rho`, `t=pi^2`, and
`z=pi(1+pi)`, both summands have value three but their difference has value
four. Equation `(20)`, not a reversed inequality, resolves the entire
cancellation face; all of its nonzero roots are simple.

## 4. Exact Keller residues and tail orders

For a self-similar cusp carrying

```text
omega=sigma^k z^B dz/x,                    t=sigma^d,     (23)
```

the `r`-tail order in sigma units is

```text
ord(omega)=k+(B+1-m/2)d*r/(m-r).                         (24)
```

For odd `m`, the conductor threshold `B=(m-1)/2` already makes the second
coefficient `1/2>0`.

On H5 the exact inherited residue becomes, up to a unit,

```text
phi^*eta_0=sigma^50 z^4 dz/x.                            (25)
```

Indeed `G=z^-6Phi` gives `G_P=-z^-4Phi_z` modulo `Phi`, so
`dS/G_P=z^4 dz/Phi_S`. With `(m,d,k,B)=(5,12,50,4)`:

| `r` | `ord(z)` | `ord(x)` | form order | genus |
|---:|---:|---:|---:|---:|
| 1 | 3 | `15/2` | `115/2` | 2 |
| 2 | 8 | 20 | 70 | 1 |
| 3 | 18 | 45 | 95 | 1 |
| 4 | 48 | 120 | 170 | 0 |

The half-integral first row becomes order `115` after the quadratic base
change; further ramification only multiplies the positive order.

On V3W, differentiating at fixed original `S`, rather than fixed `V`, gives

```text
G_P=X^-3(V Phi_V-X Phi_X) mod Phi,
phi^*eta_0=sigma^14 X^3 dX/x                             (26)
```

up to a unit. With `(m,d,k,B)=(3,4,14,3)`:

| `r` | `ord(X)` | `ord(x)` | form order | genus |
|---:|---:|---:|---:|---:|
| 1 | 2 | 3 | 19 | 1 |
| 2 | 8 | 12 | 34 | 0 |

Every positive-genus tail therefore has positive Keller-form order. Its map
to the good elliptic target is constant because a nonconstant morphism of
smooth characteristic-zero curves has nonzero differential.

## 5. The conductor firewall and exact buffer

For the special cusp

```text
A_m=k[[z,x]]/(x^2-z^m)=k[[tau^2,tau^m]],       m odd,   (27)
```

the conductor is `tau^(m-1)k[[tau]]`, whereas the ambient Jacobian numerator
ideal is `(tau^m,tau^(2m-2))A_m`. Thus

```text
conductor/J_f:
m=3: {tau^2}, length 1;             m=5: {tau^4,tau^6}, length 2.   (28)
```

Regular dualizing extension alone therefore does not imply ambient Kahler
extension, exactly as warned by THM-4289. The actual Keller residues are
stronger: on the special faces their `z^4` and `X^3` numerators lie in the
Jacobian ideal, and `(24)` proves that this buffer pays every tail scale.

## 6. Hasse selectors, exact hostiles, and what survives

For H5, write

```text
f5=u+alpha S+WS^2,
f4=Delta+eta S+xi_10 S^2+beta_11 S^3+ZS^4,
f3=e+Phi S+Theta S^2+zeta_3S^3,
f2=8/3+KS^2.                                             (29)
```

At `S0`, the first three critical-value coefficients are

```text
h1=f4,
h2=f3-(f4')^2/(4W),
h3=f2-f4'f3'/(2W)+f4''(f4')^2/(8W^2).                   (30)
```

They can all vanish. An exact allowed hostile is

```text
W=u=1, alpha=2, S0=-1,
Delta=Z=5936/105, eta=4Delta, xi_10=6Delta, beta_11=4Delta,
K=-8/3=2848/45-(7/6)Delta,
Phi=e=-1376/135, Theta=zeta_3=0.                         (31)
```

It has `(h1,h2,h3)=(0,0,0)` but

```text
psi(rho)=(-528019/18225)rho^4+(11008/405)rho^5+...,      (32)
```

so its tail is rational. For V3W take

```text
e=W=-1376/135, eta=-2e, V0=1,
Phi=-16/3, xi_10=8/3, Delta=0, K=2848/45, Z=1,          (33)
```

with `Theta=beta_11=zeta_3=0`. Its first Hasse selector vanishes but
`psi(rho)=-3rho^2+...`, again leaving a rational tail.

Hence “some displayed Hasse coefficient must be nonzero” is **REFUTED**.
The repaired, exhaustive selector is `ord_rho(psi)`: values `1,2,3` for H5
and value `1` for V3W detect every positive-genus tail; later values give
only rational or horizontally persistent geometry.

## 7. Corrected global genus and proper-flat extinction

The toric Pick genera are seventeen on H5 and sixteen on V3W. If an
`A_(r-1)` persists generically after the split, the normalized genus is

```text
g_norm=g_Pick-delta(A_(r-1))=g_Pick-floor(r/2).          (34)
```

The tail restores exactly the difference in `(21)` and attaches once. The
global checksum is therefore:

| face | `r` | `g_Pick` | persistent delta | `g_norm` | tail genus |
|---|---:|---:|---:|---:|---:|
| H5 | 1 | 17 | 0 | 17 | 2 |
| H5 | 2,3 | 17 | 1 | 16 | 1 |
| H5 | 4 | 17 | 2 | 15 | 0 |
| H5 | `>=5` or `psi=0` | 17 | 2 | 15 | 0 |
| V3W | 1 | 16 | 0 | 16 | 1 |
| V3W | 2 | 16 | 1 | 15 | 0 |
| V3W | `>=3` or `psi=0` | 16 | 1 | 15 | 0 |

In both cases the unchanged main face has genus four and the inherited graph
has `b1=11`; thus `4+11+g_tail=g_norm`. All other local pieces and toric
subdivision chains are rational.

The cusp owner addresses are

```text
H5:      p=rho^-1, s=S0 in G_m;
V3W:     p=rho^-1, s=V0*rho -> 0, sp=V0 in G_m;
Lambda:  main-face top infinity, s -> infinity, p/s^2 -> 1.   (35)
```

They are disjoint completed neighborhoods. Equations `(10)` and `(15)`
retain `Z`, so they remain valid when `Lambda=W+Z=0`; conversely,
THM-4327 Section 6's `A_23` calculation requires only `W!=0` and remains
valid on both cusp discriminants.

After finite base change, normalization, and regularization, every component
is therefore an unchanged THM-4327 constant component, a rational component,
or a positive-genus tail made constant by Sections 3--4. Retaining the actual
positive component multiplicities, proper-flat intersection gives

```text
deg(L_generic)=sum_Gamma mu_Gamma deg(L|Gamma)=0.         (36)
```

The generic response degree is positive for a nonautomorphic Keller pair,
which contradicts `(36)` and proves `(2)`. **QED.**

## 8. Unexpected exact organizers and generated tasks

The tail scale `r/(m-r)` gives a typed Stern--Brocot organizer: reduce the
pair `(r,m-r)` and use the corresponding positive rational address. It
preserves the exceptional ray and minimal base-change denominator. It loses
the parity `epsilon`, leading coefficient `u(0)`, residue buffer `(k,B)`,
and owner address, so those are mandatory sidecars. Ranking the address by a
natural number is harmless only with that typed data retained.

There is also a precise tournament warning. Orienting two candidate monomial
terms by which has smaller valuation records the generic initial monomial,
but equality is a tie, and the tie is exactly where the tail curve `(20)`
appears. Arbitrarily directing the tied edge destroys the coefficient and
genus data. The lawful object is a valuation preorder plus the initial
polynomial on every tie class, not an ordinary tournament.

The reusable mechanism suggested by `(24)` is broader than this theorem:
for any odd cusp `x^2=z^(2g+1)-psi(tz)`, a residue numerator with exponent
`B>=g` has positive contribution on every split tail. The formula is proved
formally here; applying it to a new global seam still requires its own exact
owner atlas and proper-flat component ledger.

The sharp next problems are:

1. intersect `U=0` with `W=0` or `Z=0`, where the transverse Hessian or
   interior-root condition can fail;
2. combine THM-4339's `K=0` zero-exit filtration with simultaneous
   normalization;
3. merge its `U+W=0` root collision with the separate `A_23` ladder; and
4. test the odd-cusp buffer lemma on higher-weight faces before promoting it
   beyond this exact local setting.

These are generated tasks, not extensions of theorem `(2)`.

## 9. Reproduction and scope

Run both independent certificates in normal and optimized modes and
byte-match their frozen outputs:

```bash
python3 -B 04-computation/jc2_m12_u0_repeated_cusp_extinction_thm4340.py
python3 -B -O 04-computation/jc2_m12_u0_repeated_cusp_extinction_thm4340.py
python3 -B 04-computation/jc2_m12_u0_repeated_cusp_extinction_independent_audit_thm4340.py
python3 -B -O 04-computation/jc2_m12_u0_repeated_cusp_extinction_independent_audit_thm4340.py
```

The primary certificate is standard-library only. The independent path uses
SymPy but imports neither the primary script nor its identities. What is
proved is exactly `(1)--(2)`, relative to the inherited seam and THM-4327's
proper-flat interface. The organizers in Section 8 do not enlarge the gate.
