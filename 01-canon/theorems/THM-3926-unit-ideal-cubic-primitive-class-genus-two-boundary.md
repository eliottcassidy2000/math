---
id: THM-3926
title: "Unit-ideal cubic has primitive ramification class but genus-two boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. The normal
  THM-3907 unit-ideal cubic surface is smooth, has scalar units and class
  group Z^3, and its unique ramification divisor represents a primitive
  basis vector. Two explicit unramified vertical primes complete it to a
  class-group basis; deleting those three divisors gives a smooth
  quasi-finite etale degree-three open with scalar units and trivial class
  group. Thus the full THM-3922 boundary-basis invoice passes. Nevertheless
  the ramification normalization has genus two, contradicting THM-3920's
  rational-boundary theorem for every affine-plane open. The natural
  three-divisor deletion also has compactly supported Euler characteristic
  13, not 1. This is a sharp proof that the boundary-basis condition is
  necessary but not sufficient, and it closes THM-3907 as a Keller target.
source: jc_zero_debt_lift / post-THM-3922 THM-3907 class-group lane, 2026-08-23
audit: >
  TWO INDEPENDENT HOSTILE AUDITS PASS (jc_degree6_one_place and
  incoming_thm3926_audit/root, 2026-08-23). They independently reconstructed
  the factorial chart and its complete
  unit group, the six height-one boundary primes and their valuations, the
  Nagata quotient and primitive ramification class, and the different
  divisor. It checked both the cyclic Kummer and hyperelliptic genus-two
  models, the smoothness saturations, the six infinity punctures, purity of
  the natural etale open, and the compactly supported Euler ledger. Normal
  and optimized executions agree in all 50 active gates and, after the
  Windows CRLF stream is normalized to LF, byte-match the frozen LF output;
  all raw hashes and documentation checks pass. No mathematical repair was
  needed; one stale related-theorem slug was corrected during promotion.
depends_on:
  - THM-3907-unit-ideal-nonmonogenic-cubic-six-place-boundary
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
  - THM-3922-affine-plane-open-boundary-basis-class-group-obstruction
related:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
  - THM-3924-decic-cubic-index-five-ramification-class-obstruction
script: 04-computation/jc2_unit_ideal_cubic_class_genus_boundary_thm3926.py
output: 05-knowledge/results/jc2_unit_ideal_cubic_class_genus_boundary_thm3926.out
script_sha256: 5ff37ae063072265f413cd0039c1e4e4d0d2c84ffc1cb937a47a29c42e27500b
output_sha256: cdb53643de8207a9e338c8a8f752a7603164e3a905e48180a89d700eede64540
semantic_sha256: 8deb06257def95cc087c0d2b0819cde103f9742fe1bc541f1e0258bd5e462e95
hash_basis: raw LF bytes
---

# THM-3926 -- the class invoice passes, but the boundary has genus two

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. Euler characteristics
are asserted after specialization to `k=C`.

Let `R=k[A,C]`. The THM-3907 binary cubic is

```text
Phi=A X^3+C X^2Y+(AC-1)XY^2+A Y^3.                         (1)
```

Its Delone--Faddeev algebra `S` is free over `R` with basis
`(1,omega,theta)` and multiplication

```text
omega^2=A-A^2C+C omega-A theta,
omega theta=-A^2,
theta^2=-AC+A omega+(1-AC)theta.                            (2)
```

THM-3907 proves that `S` is a normal finite-flat cubic domain, is globally
nonmonogenic, and has irreducible squarefree discriminant

```text
Delta=-4A^4C^3-27A^4+30A^3C^2+A^2C^4
      -30A^2C-6AC^3+4A+C^2.                                (3)
```

Put `X=Spec S`. This theorem determines `Cl(X)` and then shows exactly why
the favorable answer does not produce an affine-plane atlas.

## 1. A factorial chart

On `A!=0`, the endpoint `d=A` of `(1)` is a unit, so `theta` generates the
cubic algebra. Its monic equation is

```text
P(T)=T^3+(AC-1)T^2+AC T+A^3.                               (4)
```

Because `omega theta=-A^2`, both `omega` and `theta` are units after
inverting `A`. If `theta+1` is also inverted, `(4)` solves for the base
coordinate and `(2)` solves for `omega`:

```text
C=-(theta^3-theta^2+A^3)/(A theta(theta+1)),
omega=-A^2/theta.                                          (5)
```

Conversely these formulas satisfy all three relations in `(2)`. Hence

```text
S_(A(theta+1))
  =k[A,A^(-1),theta,theta^(-1),(theta+1)^(-1)],             (6)
```

a factorial ring whose units are exactly

```text
k* A^m theta^n(theta+1)^p,                 m,n,p in Z.      (7)
```

The exact Jacobian ideal of `(2)` together with all rank-two minors is the
unit ideal. Since `(2)` has codimension two, this also proves the useful new
fact

```text
X is smooth.                                                (8)
```

## 2. The six Nagata boundary primes

There are three reduced height-one primes over `A=0`:

```text
P_0=(A,omega,theta),
P_1=(A,omega,theta-1),
P_C=(A,theta,omega-C).                                     (9)
```

Their quotients are all `k[C]`. Equation `(3)` gives
`Delta(0,C)=C^2`, so the cubic map is generically etale over the line
`A=0`; all three multiplicities in `div(A)` are one.

Let `rho_1,rho_2,rho_3` be the three roots of `rho^3=2`. Since

```text
P(-1)=A^3-2,                                                (10)
```

the other three boundary primes are

```text
Q_i=(A-rho_i,theta+1,omega-rho_i^2),       i=1,2,3.        (11)
```

Each quotient is again `k[C]`, and the roots of `A^3-2` are simple. The
complete chart-boundary divisor identities are

```text
div(A)       =P_0+P_1+P_C,
div(theta)   =P_0+2P_C,
div(theta+1) =Q_1+Q_2+Q_3.                                 (12)
```

For the middle row, at `P_C` the element `omega` reduces to the generic unit
`C`, so `omega theta=-A^2` gives order two. At `P_0`, the first relation in
`(2)` has leading equation `A+C omega=0`; hence `omega` and `theta` both
have order one. These are all components because `(6)` is the complement of
the six primes in `(9),(11)`.

## 3. Exact units and class group

Nagata localization for `(6)` presents `Cl(S)` by the six prime classes in
`(9),(11)` modulo the three rows in `(12)`. The relation matrix is

```text
      P0 P1 PC Q1 Q2 Q3
A      1  1  1  0  0  0
theta  1  0  2  0  0  0
th+1   0  0  0  1  1  1.                                  (13)
```

It has Smith invariants `(1,1,1)`. Thus

```text
Cl(S)=Z g direct_sum Z q_1 direct_sum Z q_2,                (14)

[P_C]=[P_1]=g,       [P_0]=-2g,
[Q_1]=q_1,           [Q_2]=q_2,       [Q_3]=-q_1-q_2.      (15)
```

This also proves `S*=k*`, rather than assuming it. Any unit of `S` is a
unit `(7)` on the factorial chart. Its divisor coefficients in `(12)` all
vanish only when `m=n=p=0`.

## 4. The ramification class is primitive

On `A!=0`, the different is generated by

```text
D_theta=P'(theta)=3theta^2+2(AC-1)theta+AC.                 (16)
```

Let `E` be the unique ramification prime over the irreducible discriminant
`(3)`. Tameness and its squarefree exponent give `D_theta` order one on
`E`. There are no other zero divisors on `A!=0`.

At `P_1`, equation `(16)` reduces to `1`. To inspect the two zero branches,
put `theta=A u`. Dividing `(16)` by `A` and setting `A=0` gives

```text
(D_theta/A)|_(A=0)=C-2u.                                   (17)
```

The two branches of `(4)/A^2` have `u=C` at `P_0` and `u=0` at `P_C`, so
`(17)` is respectively `-C` and `C`. Both orders are exactly one. Hence

```text
div(D_theta)=E+P_0+P_C,
[E]=-[P_0]-[P_C]=g.                                        (18)
```

The ramification class is primitive. More strongly,

```text
([E],[Q_1],[Q_2])=(g,q_1,q_2)                              (19)
```

is a basis of `Cl(S)`. Thus the strongest immediate THM-3922 class-group
tests all pass: the completion has free class group, scalar units, and a
primitive ramification divisor that extends to an effective prime-boundary
basis.

## 5. The ramification boundary has genus two

The favorable class calculation is not sufficient. Eliminating `C,omega`
from `(2),(16)` on the factorial chart gives the Kummer model

```text
A^3(2theta+1)=theta^2(theta^2+2theta-1),                   (20)

C=theta(2-3theta)/(A(2theta+1)),
omega=-A^2/theta.                                          (21)
```

The rational function on the right of

```text
A^3=theta^2(theta^2+2theta-1)/(2theta+1)                   (22)
```

has valuations nonzero modulo three at the four distinct points

```text
theta=0,       theta=-1+sqrt(2),       theta=-1-sqrt(2),
theta=-1/2.                                                  (23)
```

The connected cyclic cubic is totally ramified at all four. Infinity has
valuation `-3` and is unramified. Riemann--Hurwitz gives

```text
2g(E^nu)-2=3(-2)+4(3-1)=2,
g(E^nu)=2.                                                   (24)
```

The plane Kummer model `(20)` is singular at its `theta=0` presentation
point, but this is only a projection artifact: the coordinate `C` in `(21)`
is a uniformizer there. Saturating `(2),(16)` by `A` and taking the closure
in `X`, the exact Jacobian singular ideal is the unit ideal. Thus `E` itself
is the smooth affine normalization used in the Euler calculation below.

There is an independent model directly from THM-3907's repeated-root
incidence. Completing its quadratic equation gives

```text
Y^2=t^6-12t^3+8.                                            (25)
```

The sextic is squarefree, again proving genus two. The incidence curve,
`E`, and the branch curve `(3)` have the same function field. In particular

```text
g(E^nu)=g(Gamma^nu)=2.                                      (26)
```

Any degree-three Keller affine-plane source must omit `E`. THM-3920 proves
that every irreducible boundary curve of a normal affine completion
containing `A2` has rational normalization. Equation `(26)` is therefore an
unconditional no-atlas obstruction for this completion. It closes the
THM-3907 unit-ideal candidate even though `(19)` passes the class-basis gate.

## 6. The natural boundary-basis deletion passes algebra but fails Euler

The contrast can be made completely explicit. Delete the three prime
divisors whose classes form `(19)`:

```text
U_12=X minus (E union Q_1 union Q_2).                       (27)
```

The divisor localization sequence and `(19)` give

```text
Gamma(U_12,O)^*=k*,                   Cl(U_12)=0.            (28)
```

The exact singular ideal of `X` is empty by `(8)`. Since the finite map is
generically separable and has the single ramification divisor `E`, purity
makes its restriction to `X minus E` etale. Thus `(27)` is a smooth
quasi-finite etale open of generic degree three over the target plane. No
affineness assertion about `U_12` is needed or made.

Over `C`, its Euler characteristic gives a second independent failure. The
factorial chart `(6)` is

```text
G_m times (A1 minus {0,-1}),                               (29)
```

so its compactly supported Euler characteristic is zero. The `A=0`
boundary consists of three affine lines, with `P_0` and `P_C` meeting once,
and has Euler characteristic `2`. The three disjoint `Q_i` contribute `3`.
Therefore

```text
chi_c(X)=5.                                                 (30)
```

The smooth affine ramification curve has projective genus two and six
normalization places at target infinity, so

```text
chi_c(E)=2-2(2)-6=-8.                                       (31)
```

Each `Q_i` meets `E` once, at `C=5/rho_i`. Hence
`Q_i minus E` is a copy of `G_m` and contributes zero. It follows that

```text
chi_c(U_12)=5-(-8)=13 !=1=chi_c(A2).                        (32)
```

Thus even the most natural boundary-basis deletion has the correct scalar
units and trivial class group but is not an affine plane.

## 7. Design consequence

THM-3907 remains a valuable hostile example. It proves that all of the
following can coexist:

```text
unit coefficient ideal;       no represented scalar unit;
normal smooth S3 cubic;        S*=k*;
Cl(S) free of rank three;      primitive ramification class;
effective boundary-class basis.                                (33)
```

What it cannot pay is boundary rationality. The next positive cubic design
must preserve `(33)` while replacing the genus-two ramification curve by a
rational curve whose affine singular fibres satisfy THM-3920's address cap.
This is stronger and more precise than asking only for a one-place
discriminant or nonmonogenicity. THM-3928 adds an orthogonal warning:
degenerating a torus coefficient conic to two affine lines forces a high
Cardano fold and gives only one intrinsic resolvent three-direction. Neither
statement substitutes for this theorem's actual-completion class and boundary
audit.

Reproduce the exact packet with

```bash
python3 04-computation/jc2_unit_ideal_cubic_class_genus_boundary_thm3926.py
python3 -O 04-computation/jc2_unit_ideal_cubic_class_genus_boundary_thm3926.py
```

After platform newlines are normalized to LF, both streams must byte-match
the frozen LF output named in the metadata.
**QED.**
