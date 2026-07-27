---
id: THM-2505
title: "First-collision Abel germ and geometric-depth covariance"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For two nonzero
  nonnegative rational step
  packets with disjoint support, the complete p-adic covariance shells
  form an absolutely convergent Abel germ. On the real Abel interval the
  germ is strictly negative and is exactly minus the expected collision
  mass at a geometrically sampled depth. Its order at zero is one less
  than the first positive Perron collision depth, and its leading
  derivative is the negative first-collision mass. For prime p, every
  nonzero root-character germ has the same exact order and a nonzero
  leading coefficient, with the signed and energy ledgers of THM-2471.
  All the germs are rational functions, and every colour remains nonzero
  at every rational Abel radius in (0,1]. One fixed radius, even with the
  full colour vector, does not determine collision depth.
  At a fully minimal squarefree multi-prime corner, THM-2474's Ramanujan
  ledger is exactly the full mixed Gregory difference; Pareto collision
  corners are therefore the Newton frontier of the positive collision
  series. The germ does not retain an ordinary landing frequency, target
  quotient, source/arrival time leg, or deep-root sheet. No scalar row is
  excluded and LRC(14) remains open.
source: codex-2026-07-27-first-collision-abel-germ
depends_on:
  - THM-2306-owner-normalized-disjoint-support-and-first-collision-shell
  - THM-2323-primitive-fixed-colour-cross-correlation-and-same-gauge-word-alignment
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
related:
  - THM-2313-biprime-pareto-collision-frontier-and-91-unit-current-shell
  - THM-2333-abel-target-fibre-sum-landing-and-zero-fibre-boundary
  - THM-2474-squarefree-first-collision-primitive-character-saturation
  - THM-2504-endpoint-tournament-no-go-and-root-chart-holonomy
script: 04-computation/lrc14_first_collision_abel_germ_thm2505.py
output: 05-knowledge/results/lrc14_first_collision_abel_germ_thm2505.out
script_sha256: 6cff4ab706fc377c12d38c9aec4e7d373ea8c38daeaecacad5ee5796f3c39090
output_sha256: b9e20fa20fc88269b4067a8c2243088a1da092be7be988a0244219a4c7eef773
hash_basis: working-tree bytes (LF)
---

# THM-2505 -- first collision is the zero of an Abel germ

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2306 selects the least scale at which two disjoint owner-normalized
packets overlap.  THM-2471 then opens the preceding shell into all twelve
root colours.  The selection is intrinsic: it is the order of vanishing of
one absolutely convergent analytic germ.  Equivalently, the germ samples the
collision mass at a geometric random depth.

This also identifies the discrete/continuous dictionary precisely.  Perron
depth is the discrete variable, its shell current is a Gregory difference,
and the Abel parameter packages all depths into a Maclaurin germ.  The
squarefree Ramanujan sign in THM-2474 is the corresponding mixed Gregory
difference at a Newton-frontier corner, not a separate cyclotomic accident.

## 1. The Abel covariance germ

Let `p` be prime and let `A,B:T->R_(>=0)` be nonzero rational-valued step
functions with rational breakpoints.  Assume

```text
A B=0                                      almost everywhere.       (1)
```

Write

```text
P h(y)=(1/p)sum_(u in F_p)h((y+u)/p),

A_s=P^s A,                 B_s=P^s B,

I_s=integral_T A_s B_s.                                      (2)
```

Then `I_0=0`, every `I_s>=0`, and

```text
I_s -> a b>0,

a=integral_T A,            b=integral_T B.                   (3)
```

Consequently the first collision depth

```text
r=min{s>=1:I_s>0}                                            (4)
```

exists.

For `k in F_p`, define the exact residue current

```text
J_s(k)
 =sum_(n congruent k mod p)
    A_hat(p^s n) conjugate(B_hat(p^s n)).                     (5)
```

For `k!=0`, put

```text
Gamma_k(z)=sum_(s>=0)J_s(k)z^s,

Gamma(z)=sum_(k!=0)Gamma_k(z).                               (6)
```

Every series in (6) is absolutely convergent for `|z|<=1`, and is analytic
for `|z|<1`.  The aggregate has the three exact descriptions

```text
Gamma(z)
 =sum_(s>=0)(I_s-I_(s+1))z^s

 =sum_(q!=0) z^(nu_p(q))
    A_hat(q)conjugate(B_hat(q))

 =-(1-z)sum_(s>=1)I_s z^(s-1),             |z|<1.             (7)
```

In particular,

```text
Gamma(z)<0,                              0<z<1,               (8)

Gamma(1)=-a b.                                                 (9)
```

Here (9) is the absolutely convergent boundary value.

There is also an exact operator form.  Let `E_s` be the Fourier projection
onto frequencies divisible by `p^s`, let `U h=h composed_with (x->px)`, and
let `W_z` kill the zero mode and multiply every nonzero frequency `q` by
`z^(nu_p(q))`.  Then

```text
W_z=sum_(s>=0)z^s(E_s-E_(s+1))

   =Id-(1-z)sum_(s>=1)z^(s-1)E_s,                           (10)

W_z U=z U W_z,                    P W_z=z W_z P,             (11)

Gamma(z)=<W_z A,B>.                                         (12)
```

The zero mode cancels in the second expression of (10).  Equations (11)
make the Abel parameter an exact clock-degree variable rather than a formal
weight.

### Proof

The Perron multiplier identity gives

```text
(P^s h)_hat(n)=h_hat(p^s n).
```

Parseval and Cauchy--Schwarz therefore give

```text
I_s=sum_n A_hat(p^s n)conjugate(B_hat(p^s n)),                (13)
```

with absolute convergence.  Separating `n` by its residue modulo `p` gives

```text
sum_k J_s(k)=I_s,

J_s(0)=I_(s+1),                                                (14)
```

and hence the first line of (7).  For `k!=0`, the map

```text
(s,n) -> q=p^s n,                 n congruent k mod p          (15)
```

is injective and its image consists precisely of the nonzero integers whose
`p`-free part is congruent to `k`.  Thus

```text
sum_(s,k!=0)|J_s(k)|
 <=sum_(q!=0)|A_hat(q)B_hat(q)|
 <=||A||_2||B||_2.                                             (16)
```

This proves absolute convergence and the second line of (7).  Since `I_0=0`,
ordinary summation by parts gives the third line.  Equation (8) follows from
`I_s>=0` and (3).  At `z=1`, (1), Parseval, and the zero Fourier mode give

```text
Gamma(1)
 =sum_(q!=0)A_hat(q)conjugate(B_hat(q))
 =integral A B-a b
 =-a b.                                                        (17)
```

Finally, (3) is the standard `L^2` Perron convergence used in THM-2306:
`P^s h` converges to `integral h`. QED.

## 2. Geometric depth and the exact first-collision derivative

For `0<=z<1`, let `S_z` be geometric on the positive integers:

```text
Pr(S_z=s)=(1-z)z^(s-1),                  s>=1.                 (18)
```

Equation (7) is exactly

```text
-Gamma(z)=E[I_(S_z)].                                         (19)
```

Thus strict negativity is not cancellation luck: it is the expectation of a
nonnegative collision observable.  More sharply,

```text
ord_(z=0) Gamma=r-1,

Gamma^(r-1)(0)=-(r-1)! I_r.                                  (20)
```

Indeed, minimality in (4) gives

```text
I_s-I_(s+1)=0,                    0<=s<r-1,

I_(r-1)-I_r=-I_r.                                        (21)
```

So the first spatial collision is recoverable from the germ without first
choosing a shell.  Near zero,

```text
-Gamma(z)=I_r z^(r-1)+O(z^r).                                (22)
```

This is the exact Abel--Dini boundary datum: the first nonzero Maclaurin
coefficient records both collision depth and collision mass.

## 3. Every root-colour germ has the same order

For each `s`, put

```text
U=A_s,                         V=B_s,

U_u(y)=U((y+u)/p),             V_u(y)=V((y+u)/p),

alpha_(s,k)(y)=(1/p)sum_u U_u(y)zeta_p^(-ku),

phi_(s,k)(y)=(1/p)sum_u V_u(y)zeta_p^(-ku).                  (23)
```

Fourier expansion of (20) gives the physical identity

```text
J_s(k)=integral_T alpha_(s,k)conjugate(phi_(s,k)).            (24)
```

For every `k!=0`,

```text
ord_(z=0) Gamma_k=r-1.                                       (25)
```

If

```text
c_k=[z^(r-1)]Gamma_k=J_(r-1)(k),                              (26)
```

then

```text
c_k!=0                                      for every k!=0,

sum_(k!=0)c_k=-I_r,

sum_(k!=0)|c_k|^2>=I_r^2/(p-1),

max_(k!=0)|c_k|>=I_r/(p-1),                                  (27)
```

and some `k` satisfies

```text
Re c_k<=-I_r/(p-1).                                          (28)
```

### Proof

If `s<r-1`, then `I_(s+1)=0`.  Nonnegativity implies pointwise that the
two root totals

```text
sum_u U_u(y),                     sum_u V_u(y)                 (29)
```

cannot both be positive.  Hence every cross-root product
`U_(u+t)(y)V_u(y)` vanishes, every root correlation vanishes, and

```text
J_s(k)=0                         for every k.                   (30)
```

At `s=r-1`, `UV=0` but `(PU)(PV)` has integral `I_r>0`.  The rational
nonnegative cyclic correlation

```text
C(t)=integral_T sum_u U_(u+t)(y)V_u(y)dy                       (31)
```

has `C(0)=0` and positive total.  If its transform vanished at one nonzero
`k`, its rational coefficient polynomial would be divisible by
`Phi_p=1+X+...+X^(p-1)`, forcing all entries of `C` equal.  This contradicts
`C(0)=0` and `sum C>0`.  Since

```text
sum_t C(t)zeta_p^(-kt)=p^2 J_(r-1)(k),                        (32)
```

every leading coefficient is nonzero.  Equations (11) and (18) give the
signed sum; Cauchy--Schwarz gives (24), and averaging real parts gives (25).
QED.

For `p=13`, (25)--(28) say that all twelve THM-2471 colours are synchronized
as analytic germs: no colour is born earlier or later than another.

## 4. Rationality and all-colour survival at one Abel radius

The step-function hypothesis gives more than unit-disc convergence.  If `h`
is a rational step function with rational breakpoints and mean `m_h`, then

```text
H_s^h=p^s(P^s h-m_h)                                        (33)
```

is eventually periodic as a rational step function.  To see this, first take
one rational interval.  The value of `p^sP^s 1_[alpha,beta)` is the number of
integers in an interval with endpoints `p^s alpha-y` and `p^s beta-y`.
After subtracting `p^s(beta-alpha)`, the resulting step function depends only
on the fractional parts of `p^s alpha,p^s beta`.  Those are eventually
periodic because the endpoints are rational.  Rational linearity proves
(33).

Since every `H_s^h` has mean zero,

```text
I_s=a b+p^(-2s)R_s,                                         (34)
```

where `R_s` is eventually periodic and rational.  Likewise, for every
`k!=0`,

```text
J_s(k)=p^(-2s)R_(s,k),                                      (35)
```

with `R_(s,k)` eventually periodic in `Q(zeta_p)`.  Therefore `Gamma` and
every `Gamma_k` are rational functions of `z`, with Taylor radius at least
`p^2`.

For every rational

```text
0<z<=1,                                                     (36)
```

all primitive colours remain nonzero:

```text
Gamma_k(z)!=0                         for every k!=0.         (37)
```

Indeed, define the centered rational correlation

```text
K_z(t)=sum_(s>=0)z^s
  (p integral_T A_s(x+t/p)B_s(x)dx-p a b).                  (38)
```

Equations (33)--(35) make every `K_z(t)` rational when `z` is rational,
including at `z=1`.  Moreover

```text
K_z(0)-(1/p)sum_t K_z(t)=p Gamma(z)<0,                       (39)
```

where (9) handles `z=1`.  Thus `K_z` is nonconstant.  Its nonzero Fourier
transforms are `p^2 Gamma_k(z)`.  If one vanished, rational divisibility by
`Phi_p` would force all `p` entries of `K_z` equal, contradicting (39).

One rational radius does **not** retain the depth.  Let `C_p(j;m)` denote the
`j`-th half-open cell of the `p^m` grid.  Compare

```text
X:
 A_X=t 1_(C_p(p;2)),             B_X=1_(C_p(0;2)),

Y:
 A_Y=1_(C_p(p^2+p;3)),           B_Y=1_(C_p(0;3)).           (40)
```

The first collision depths are respectively `1` and `2`, while direct root
disintegration gives, for every `k!=0`,

```text
Gamma_k^X(z)=t(p^(-3)zeta_p^(-k)+p^(-4)z),

Gamma_k^Y(z)=(z/p^2)(p^(-3)zeta_p^(-k)+p^(-4)z).             (41)
```

At `z=1/2`, choose `t=1/(2p^2)`.  Then the complete primitive colour vectors
are identical although the collision depths differ.  The *full germ*, its
zero order, or a separate depth/mass coordinate is therefore load-bearing.

## 5. Gregory differences are the squarefree Ramanujan ledger

Let

```text
N=p_1...p_d
```

be squarefree.  On one local Perron cube, write

```text
I_S=integral_T
  (P_(product_(i in S)p_i)U)
  (P_(product_(i in S)p_i)V),          S subset {1,...,d}.    (42)
```

Assume the top is a fully minimal collision:

```text
I_{1,...,d}=I>0,

I_({1,...,d} minus {i})=0              for every i.           (43)
```

Positivity is upward under further Perron quotienting, so every proper-face
value in (30) is zero.  With `E_i` denoting the forward cube shift, define the
signed Gregory operator

```text
nabla_i=1-E_i.                                               (44)
```

Then

```text
(nabla_1...nabla_d I)_emptyset
 =sum_S (-1)^|S| I_S
 =(-1)^d I
 =mu(N)I.                                                     (45)
```

The right side is exactly THM-2474's sum of all primitive root currents.
Thus its alternating Ramanujan sign is the full mixed discrete derivative of
the collision table.

More generally, in the multivariable positive collision series

```text
M(z_1,...,z_d)=sum_(a in N^d) I_a z_1^(a_1)...z_d^(a_d),     (46)
```

the minimal support monomials are precisely the Pareto collision frontier.
The number of active axes at such a monomial is THM-2474's collision spectral
rank, while (33) is its local mixed Gregory coefficient.  This is the exact
Maclaurin/Gregory--Newton duality:

```text
Maclaurin support       = where a collision is born;

Gregory difference     = its primitive signed current;

cyclotomic transform   = the individual directional colours.             (47)
```

## 6. Exact information boundary

The Abel germ retains:

- the first collision depth and mass;
- the complete `p`-adic valuation-shell covariance;
- every chart-oriented nonzero root colour at the first collision; and
- at squarefree corners, the Newton-frontier rank and mixed Gregory sign.

It deliberately forgets:

- which ordinary frequency inside one residue progression lands;
- the THM-2305 source word versus its arrival-time copy;
- the THM-2365 target quotient and deepest-probe root; and
- the ancestry sheet needed to compare collision and deep roots.

The finite-cylinder hostile in the companion makes the first loss literal.
Its two interval atoms have a nonzero ordinary `q=1` Fourier product already
at depth zero, while every depth-zero root-residue aggregate is zero; the
individual terms cancel inside the residue class.  Therefore neither a
single ordinary term nor the scalar germ can replace the coloured vector and
its endpoint-Prony sidecar.

Equation (41) supplies a second sharp loss: even the complete colour vector at
one rational radius does not recover `r`.  Only the germ order or an attached
depth coordinate does.

This valuation Abel parameter is also distinct from THM-2333's base-index
Poisson parameter `rho^(||u||_1)`.  The two regularize different coordinates
and may be used bivariately, but neither identifies a target quotient with a
predecessor root.

Consequently THM-2505 changes the canonical description of the open object,
but it does not provide the off-diagonal temporal endpoint intertwiner left by
THM-2471, identify an owner/deep chart, remove one of the `165` scalar rows, or
prove LRC(14).

## 7. Exact companion

Run

```text
python3 04-computation/lrc14_first_collision_abel_germ_thm2505.py
python3 -O 04-computation/lrc14_first_collision_abel_germ_thm2505.py
```

The dependency-free `Fraction` companion:

- realizes first collision at depth three on exact `p^4` cylinders for
  `p=2,3,5,7,13`;
- checks every collision mass and shell coefficient;
- reduces every primitive leading transform exactly modulo `Phi_p`;
- checks that all coloured germs have order two and the aggregate leading
  coefficient is the negative collision mass;
- verifies all colours at three rational radii and the exact valuation-clock
  operator identities on `1,000` signed frequencies per prime;
- realizes two different collision depths with identical full colour vectors
  at radius `1/2`;
- verifies the Abel summation identity and strict sign at several rational
  parameters, including its constant tail;
- verifies the boundary covariance at `z=1` and the first four mixed Gregory
  signs; and
- checks the ordinary-frequency-versus-residue-aggregate hostile.

Normal and optimized runs must reproduce

```text
05-knowledge/results/lrc14_first_collision_abel_germ_thm2505.out
```

byte-for-byte.

An independent audit rederived the Fourier normalization, Abel
summation-by-parts identity, valuation-multiplier covariance, zero-order
formula, eventually periodic rational-step discrepancy, rational-function
upgrade, all-colour rational-radius theorem, squarefree Gregory sign, and the
fixed-radius depth hostile.  It independently ran normal and optimized
companions with byte-identical `PASS` transcripts before promotion.

QED.
