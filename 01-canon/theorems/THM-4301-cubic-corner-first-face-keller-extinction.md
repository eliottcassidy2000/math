---
id: THM-4301
title: "Cubic-corner first-face Keller extinction"
status: >
  PROVED RELATIVE TO THM-4103 + THM-4230 + THM-4299 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. At the exact-M=12 corner D=Lambda=0 with UZ!=0,
  every q-separable first Newton-face component has good-differential order
  at least s+3beta>0 and is Keller-constant. Every q-horizontal first face is
  binomial/rational. Hence only q-repeated discriminant refinements remain
  locally. A literal five-term genus-four face attains order seven. This does not
  close the cubic corner, seam entry, JC(2), or DC(2).
source: root / alternating LRC14-JC2 session, 2026-09-01
depends_on:
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4299-d-zero-square-face-elliptic-splitting-and-off-corner-extinction
related:
  - THM-4292-lambda-zero-repeated-face-keller-extinction
  - THM-4297-general-lambda-zero-central-and-tail-keller-extinction
  - THM-4298-weighted-face-source-normal-unimodular-visibility-transform
  - THM-4300-lrc14-size-preserving-response-staircase-and-index-297-ideal
primary_script: 04-computation/jc23_cubic_corner_first_face_keller_extinction_thm4301.py
primary_output: 05-knowledge/results/jc23_cubic_corner_first_face_keller_extinction_thm4301.out
primary_script_sha256: b22a016dcb24e66c21eca6ae5b4db2c4928e0f82364189e00a216b43f8bddfeb
primary_output_sha256: 1e998e9d8bd8d8368d6b9eb46b293df7a84cd06972b978932518bd9f255cabd8
independent_audit_script: 04-computation/jc23_cubic_corner_first_face_keller_extinction_independent_audit_thm4301.py
independent_audit_output: 05-knowledge/results/jc23_cubic_corner_first_face_keller_extinction_independent_audit_thm4301.out
independent_audit_script_sha256: d001bb870e468fd898d4861c240f25ea376298873df3eedb7885540f7a77c19d
independent_audit_output_sha256: 84af8fb5ac4230fdff581a0490bc39923f8549cb8596970c812e4357cbe8fe7e
hash_basis: raw LF bytes
audit: >
  PASS. The primary SymPy path reconstructs the literal Hhat support, the
  q-linear response, every possible first response order, horizontal
  binomials, and a sharp literal five-term genus-four hostile. A dependency-free
  rational-arithmetic path independently checks the response orders,
  inequalities, binomial components, hostile discriminant, irreducibility,
  genus, and differential order. The theorem text supplies the exhaustive
  two-monomial valuation bound covering all q^2 and mixed faces. Normal and
  optimized streams byte-match.
---

# THM-4301 -- Cubic-corner first-face Keller extinction

**PROVED RELATIVE TO THM-4103 + THM-4230 + THM-4299 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THE REPEATED-FACE TOWER, CUBIC-CORNER EXTINCTION,
SEAM ENTRY, `JC(2)`, AND `DC(2)` REMAIN OPEN.**

## 1. Statement and inheritance

Retain the exact-weight-twelve reduced `(2,3)` seam, put

```text
D=W^2-4UZ,                         Lambda=U+W+Z,
```

and work over `C` on

```text
D=Lambda=0,                       UZ!=0.                    (1)
```

Then `(U,W,Z)=(U,-2U,U)` with `U!=0`. At the triple contact use

```text
z=1/S,               r=P/S^2,               q=r-1,
t=sigma*z.                                                   (2)
```

> **Theorem.** Let `v` be a divisorial valuation centered at
> `sigma=z=q=0`, and let `Gamma` be an irreducible component of its first
> weighted face.
>
> 1. If the face is generically separable in `q` along `Gamma`, then
>    the pullback of the good target differential has vertical order at
>    least `s+3beta>0`; hence `Gamma` is Keller-constant.
> 2. If the irreducible factor defining `Gamma` is `q`-independent after
>    torus saturation, it is a binomial with rational normalization and is
>    again constant.
>
> Consequently a first-face component can survive only when its `q`-
> derivative vanishes identically for a genuine repeated `q`-factor. Such
> repeated/discriminant refinements are not closed here.

The inheritance pass is:

- closest proved mechanism: THM-4297's valuation of the good differential;
- canonical hostile: THM-4299's literal cubic face at `(1,-2,1)`;
- corrected near miss: MISTAKE-531 forbids transporting a simple-root packet
  through a repeated face, and THM-4299 forbids division by `q=r-1` here;
- least-used sidecar: the complete `q`-linear coefficient of the literal
  `Hhat`, including its forced `t^8` term.

The live board was

```text
{cubic unit, constant response, q-linear row, separable derivative,
 repeated discriminant, physical consumer}.                            (3)
```

## 2. Exact local support

Specialize THM-4230's exact source `(4)` and integral model `(9)` at
`W=-2U`, `Z=U`, then use `(2)`. Direct expansion, without division by `D`
or `q`, gives

```text
F=q(Hhat(r,t)-z^12)-t^12/2.                              (4)
```

Writing the resulting scaled source polynomial as `Hhat` and expanding at
`r=1+q` gives

```text
F=-t^12/2+q(h(t)-z^12)+q^2 A(t)+q^3 C(q,t),
C(0,0)=U!=0,                                               (5)

h(t)=c_1 t+c_2 t^2+c_3 t^3+c_4 t^4+c_5 t^5+c t^6
     +(8/3)t^8-3t^10,                                     (6)

c_1=alpha_11+beta_11,          c_2=upsilon_5+xi_10,
c_3=eta+zeta_3,                c_4=Delta+Theta,
c_5=Phi,                       c=7168/135-(7/6)Delta.      (7)
```

Thus three pieces never disappear:

```text
q^3 * unit,                         -t^12/2,
q * ((8/3)t^8 plus lower or higher rows minus z^12).       (8)
```

If `ell` is the first nonzero `t`-order in `h`, then exactly

```text
ell in {1,2,3,4,5,6,8}.                                  (9)
```

The fixed coefficient `8/3` makes `(9)` exhaustive. This is a literal
support statement, not Weierstrass preparation and not a genericity claim.

## 3. Uniform extinction of every q-separable first face

Write

```text
s=v(sigma)>0,             beta=v(z)>0,
tau=s+beta=v(t),          gamma=v(q)>0.                   (10)
```

Let `N` be the minimum monomial weight defining the first face. The two
coefficient-independent terms in `(8)` imply

```text
N<=3gamma,                         N<=12tau.               (11)
```

Suppose the initial `q`-derivative is nonzero at the generic point of
`Gamma`. Differentiation lowers the face weight by exactly `gamma`, so

```text
v(F_q)=N-gamma
      <=min(2gamma,12tau-gamma)
      <=8tau.                                             (12)
```

The last maximum occurs only at `gamma=4tau`.

The good-model form also crosses the wall directly; importing THM-4297's
off-wall conclusion is unnecessary. The inherited THM-4103 residue identity
and THM-4230 model are

```text
phi^*(dA/(2 C_target))=Q ds/(F_Q)_p,
Q=sigma^12,       s=sigma^-1 S,       p=sigma^-2 P,
F_Q=sigma^-2 G,   A=sigma^-4 X,       C_target=sigma^-6 Y,
eta_0=dX/(2Y).
```

Taking relative differentials gives `(F_Q)_p=G_P` and therefore

```text
phi^*eta_0=sigma^9 dS/G_P.                               (13)
```

Here the local equation satisfies `F=z^14 G`. Since
`P=(1+q)z^-2` and `S=z^-1`, differentiation at fixed `z,t` gives

```text
F_q=z^12 G_P,                         dS=-z^-2 dz.
```

Therefore, up to the global orientation sign,

```text
phi^*eta_0=-sigma^9 z^10 dz/F_q.                          (14)
```

Since `v(d_rel z)>=beta`, equations `(10)--(12)` and `(14)` give

```text
v(phi^*eta_0)>=9s+11beta-8(s+beta)=s+3beta>0.             (15)
```

The restriction of the good target differential to `Gamma` is zero. A
nonconstant map of smooth characteristic-zero curves has nonzero
differential, so every positive-genus `q`-separable component is constant.
A rational component maps constantly to the good elliptic target as well.

There is a finer response bound from `(6)`: if the first `q`-linear row has
order `ell`, then `N<=gamma+ell*tau` and

```text
v(phi^*eta_0)>=(9-ell)s+(11-ell)beta.                      (16)
```

For the seven values in `(9)`, `(16)` is positive; its weakest row is again
`s+3beta` at `ell=8`.

## 4. Derivative-null coefficient faces are rational

Suppose the irreducible factor defining a first-face component is
`q`-independent after torus saturation. The constant response in `(5)`
cannot participate: it is the
only `q`-degree-zero monomial. For every `q`-degree above one, the first
coefficient row is a single monomial. The only possible cancellation is
therefore between

```text
q*b_ell*t^ell                         and -q*z^12,          (17)
```

where `b_ell!=0` and `ell` is as in `(9)`. Substituting `t=sigma*z` and
removing `q*z^ell` turns `(17)` into

```text
z^(12-ell)=b_ell*sigma^ell.                               (18)
```

At the factor level, a `q`-independent irreducible component must divide
every coefficient of the initial Laurent polynomial in `q`. The `q^0`
coefficient and the first coefficient at each `q`-degree at least two are
single torus monomials. The only nonmonomial first coefficient is `(17)`,
and the presence of any simultaneous higher-`q` monomial prevents its
binomial factor from being common. Thus `(18)` is exhaustive, not merely a
classification of faces that are globally independent of `q`.

Equation `(18)` splits into `gcd(12-ell,ell)` primitive monomial curves;
every normalization is rational. These coefficient-cancellation faces add
no positive-genus carrier. Any remaining derivative-null face depends on
`q`, so its component is a repeated `q`-factor and belongs to the unresolved
discriminant-refinement tower.

## 5. The lower bound is sharp on a genuine genus-four tail

The weak row `ell=8` is realizable in the literal source equation. Set

```text
U=1,
alpha_11=beta_11=upsilon_5=xi_10=eta=zeta_3=Phi=0,
Delta=2048/45,                  Theta=-Delta.              (19)
```

Then `c_1=...=c_5=c=0`. At primitive weights

```text
(s,beta,gamma)=(1,2,12),               N=36,              (20)
```

the complete first face, including the `q^2` translation sidecar, is

```text
q^3+(2048/45)q^2t^4+(8/3)qt^8-qz^12-t^12/2.              (21)
```

In the toric chart

```text
x=sigma^2/z,                         y=q/z^6,              (22)
```

division by `z^18` normalizes `(21)` to

```text
G=y^3+(2048/45)x^2y^2+((8/3)x^4-1)y-x^6/2.               (23)
```

As a cubic in `y`, its discriminant is

```text
(24553315427 x^12-1282042880 x^8
 +247770240 x^4+486000)/121500.                           (24)
```

It is squarefree of degree twelve. The monic cubic `(23)` is irreducible
over `C(x)`: a rational root would be a polynomial of degree at most two;
coefficient comparison for `ux^2+vx+w` forces `w in {0,+1,-1}`, then `v=0`,
and gives an immediate contradiction in both `w=0` and `w^2=1` cases.
At infinity, putting `x=1/u`, `y=Y/u^2` gives three distinct roots of

```text
Y^3+(2048/45)Y^2+(8/3)Y-1/2.                             (25)
```

Thus the degree-three map to the `x`-line has twelve simple finite branch
points and no ramification at infinity. Riemann--Hurwitz gives genus four.

On this divisor

```text
v(F_q)=N-gamma=24,
v(phi^*eta_0)=9+11*2-24=7.                               (26)
```

Hence a genuine positive-genus cubic tail exists and attains the lower bound,
but its Keller map is constant. This is the canonical hostile against
replacing the theorem by a rationality assertion.

## 6. Exact failure boundary and generated task

The theorem stops precisely when the first face has a genuine repeated
`q`-factor. There `F_q` has larger valuation than `N-gamma`, so `(12)` cannot
be imported. The next native object is not a scalar critical value. For the
top cubic `Uq^3`, infinitesimal right-coordinate changes generate `(q^2)`,
and

```text
T^1=C[[q]]/(q^2)=span{1,q}.                               (27)
```

The raw two coordinates are the constant and `q`-linear rows in `(5)`, but
the `q^2 A(t)` translation feeds the depressed-cubic pair nonlinearly.
Therefore the required sidecar is the full prepared pair plus discriminant;
the raw quotient `(27)` is not valuation-complete.

This also gives a strict firewall against a tempting LRC analogy. THM-4300
may delete a mask only after proving it inactive on every row of its consumer.
Here `q=r-1` lies in the maximal ideal and

```text
F mod q=-t^12/2!=0.                                       (28)
```

Dividing by `q` creates a pole. No map from LRC masks/bodies to the cubic
pair is supplied, and no coverage or Keller predicate is transferred.

## 7. Reproduction and scope

Run

```bash
python3 -B 04-computation/jc23_cubic_corner_first_face_keller_extinction_thm4301.py
python3 -B -O 04-computation/jc23_cubic_corner_first_face_keller_extinction_thm4301.py
python3 -B 04-computation/jc23_cubic_corner_first_face_keller_extinction_independent_audit_thm4301.py
python3 -B -O 04-computation/jc23_cubic_corner_first_face_keller_extinction_independent_audit_thm4301.py
```

The result concerns first faces over the inherited exact-`M=12` corner `(1)`.
It does not extinguish repeated refinements, complete the component inventory,
prove seam entry, cross `U=0` or `Z=0`, or prove `JC(2)` or `DC(2)`.

**QED.**
