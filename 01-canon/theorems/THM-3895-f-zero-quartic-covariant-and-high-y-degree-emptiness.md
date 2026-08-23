---
id: THM-3895
title: "F-zero quartic covariant and high-y-degree emptiness"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
  independent hostile audit.  Every polynomial square in the THM-3885 f=0
  residual has y-degree at most two.  A high-degree square canonically
  produces a degree-four covariant Z.  The candidate T must divide two
  coprime shifted fibres of Z and its remaining quadratic discriminant
  reduces to one of four normalized square ideals.  Their exact reduced
  bases force the required sidecar to vanish, contradicting its degree.
  Classification of the remaining quadratic-y coefficient functions, a
  polynomial-plane Keller atlas, and JC(2) remain OPEN.
source: jc_zero_debt_lift / post-THM-3888 elliptic denominator compression, 2026-08-23
audit: >
  SELF-AUDITED EXACT PROOF CANDIDATE.  The companion verifies the quartic
  covariant relation, its two-colour divisor identity, the normalized
  quadratic discriminant, the complete degree ladder, and the four exact
  square-root ideals over Q(f).  Their reduced grevlex bases are respectively
  [b0], [b0,b1], [b0,b1,b2], and [b0,b1,b2,b3].  It separately checks the
  linear-row terminal factors, the constant-row parity obstruction, and the
  B=0 positive square in 43 active gates.  Normal and optimized runs must
  byte-match the frozen output.  Independent audit must recheck the choice
  of equianharmonic sign, every degree comparison, the algebraic
  normalization A0=1, transcendence of f, and the exact Groebner bases.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
  - THM-3885-cusp-residual-f-zero-arm-dichotomy-and-quadratic-closure
related:
  - THM-3888-f-zero-equianharmonic-jacobian-and-two-section-integrality
  - THM-3893-cusp-residual-f-zero-two-arm-degree-six-closure
script: 04-computation/jc2_f_zero_quartic_covariant_high_y_emptiness_thm3895.py
output: 05-knowledge/results/jc2_f_zero_quartic_covariant_high_y_emptiness_thm3895.out
script_sha256: 64be7ed6b4c77b91a40e0861be7cd34bef16dbbd5e548613f581a5b8e13e66ed
output_sha256: 8e2aadf3d944ffc1dd3daa3add7f950efbf783e919c2d67060737e15a2386a56
semantic_sha256: b2c261b027357a48476f9c9d827d6068fcbc9ee483b633b3eeb27aa6a2c801df
hash_basis: raw LF bytes
---

# THM-3895 -- the f-zero residual has only three y-channels

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
independent hostile audit.**  Work over an algebraically closed field `k` of
characteristic zero.  Retain

```text
a=x+1,                 L=9x+4,
F=15x^2+15x+4,         K=y^2-F.                            (1)
```

Suppose `T,G in k[x,y]` satisfy the THM-3885 `f=0` square equation

```text
G^2=L^4-6aL^2T^2-8KT^3-3a^2T^4.                         (2)
```

Then

```text
deg_y(T)<=2.                                               (3)
```

Thus the apparent all-degree two-variable problem has the complete remaining
shape

```text
T=t_2(x)y^2+t_1(x)y+t_0(x).                               (4)
```

No bound on the three `x`-coefficient functions is inferred here.  In
particular, `(3)` is a dimension reduction, not a Keller counterexample or a
proof of `JC(2)`.

## 1. Inheritance and the missing object

THM-3888 identifies `(2)` with a `j=0` rational elliptic surface and closes
the integral Weierstrass-coordinate shell.  Its canonical hostile point has

```text
T_*=-2K/(3a^2),                                            (5)
```

which already has y-degree two and shows that a generic `k(x)[y]` lane is
not empty.  THM-3893 instead controls global descent through total degree six
when both labelled arms vanish.  What neither view exposes is that every
higher-y-degree square carries a **quartic polynomial sidecar**.  That
sidecar, rather than `T` itself, is the bounded object.

Pass to `C=overline{k(x)}` and regard `(2)` as an equation in `C[y]`.  It is
enough to prove `(3)` there.

## 2. The equianharmonic quartic covariant

Assume for contradiction that

```text
n=deg_y(T)>2.                                              (6)
```

Write `t_n` for the leading coefficient of `T`.  The unique highest term of
the right side of `(2)` is `-3a^2t_n^4y^(4n)`.  Hence
`deg_y(G)=2n` and, for

```text
sigma=lc_y(G)/(a t_n^2),
```

one has `sigma^2=-3`.  Choose

```text
d=-sigma,                 d^2=-3,                 d sigma=3. (7)
```

Define

```text
Z=(3aL^2+4KT+3a^2T^2-daG)/2.                             (8)
```

The two degree-`2n` terms in `(8)` cancel by `(7)`, so

```text
m=deg_y(Z)<2n.                                             (9)
```

Here `Z` is nonzero: if it vanished, the term `4K^2T^2` in `(10)` below
would have unique degree `2n+4`.

Direct substitution of `(2)` gives the exact quadratic relation

```text
Z^2-(4KT+3aL^2+3a^2T^2)Z
   +4K^2T^2+6aKL^2T+3a^2L^4=0.                          (10)
```

This identity is the bounded replacement for the original square root.

If `m>4`, the term `-3a^2T^2Z` has degree `2n+m`, strictly larger than
every other term of `(10)`: use `(9)` against `Z^2`, and `n>2,m>4` against
the remaining terms.  If `m<4`, the term `4K^2T^2` has unique degree
`2n+4`.  Both cases contradict `(10)`.  Therefore

```text
deg_y(Z)=4,                 lc_y(Z)=4/(3a^2).             (11)
```

The leading coefficient follows by comparing the two degree-`2n+4` terms in
`(10)`.

## 3. Two shifted quartic colours and the degree-six ladder

Put

```text
Phi(Z)=Z^2-3aL^2Z+3a^2L^4,
A=2K(2Z-3aL^2),
B=3a^2Z-4K^2.                                             (12)
```

Equation `(10)` is exactly

```text
Phi(Z)=T(A+BT).                                            (13)
```

In particular, `T|Phi(Z)` in `C[y]`.  The two intrinsic root colours are

```text
Phi(Z)=
 (Z-(3+d)aL^2/2)(Z-(3-d)aL^2/2).                         (14)
```

They are coprime because their difference is the nonzero unit `daL^2` in
`C[y]`.  Thus every root of `T` is assigned to one of two shifted fibres of
the quartic map `Z:P^1_y -> P^1`.  This is the exact binary observable; its
sidecar is the quotient in `(13)`, not an arbitrary tournament orientation.

By `(11)`,

```text
deg Phi=8,                  deg A=6,                  deg B<=3. (15)
```

The divisibility `T|Phi` first gives `n<=8`.  Let `E=Phi/T`, so
`deg E=8-n`.  If `n=7` or `8`, then a nonzero `BT` has degree at least seven
and cannot cancel `A`, while `B=0` leaves the degree-six polynomial `A`;
either conclusion contradicts `deg E<=1`.  Hence `n<=6`.

For `3<=n<=6`, the right side

```text
E=A+BT                                                    (16)
```

has degree at most five.  Its degree-six part can cancel if and only if

```text
deg B=6-n.                                                 (17)
```

This gives the complete high-degree ladder

```text
n:          3   4   5   6
deg B:      3   2   1   0.                                (18)
```

## 4. One normalized square problem closes all four rows

Treat `(13)` as a quadratic equation

```text
BT^2+AT-Phi=0.                                             (19)
```

Because `(17)` makes `B` nonzero, its discriminant must be the square
`(2BT+A)^2`.  Exact simplification gives

```text
A^2+4B Phi
 =4/(9a^4) * N_B,

N_B=(4K^2+B-3A_0)^3+27A_0^2(A_0-K^2),
A_0=a^3L^2.                                                (20)
```

The prefactor is already a square.  Since `A_0!=0` in `C`, choose
`c,e in C` with

```text
c^2=A_0,                    e^2=c,
Y=y/e,                      f=F/c,
b(Y)=B(eY)/A_0.                                            (21)
```

Then `K/c=Y^2-f`, and after removing the square factor `A_0^3=c^6`, the
remaining square test is

```text
N_b=(4(Y^2-f)^2+b(Y)-3)^3+27(1-(Y^2-f)^2).               (22)
```

The element `f` is transcendental over the constant field: already
`f^2=F^2/(a^3L^2)` is a nonconstant function of `x` (its numerator and
denominator have different degrees at infinity).  Hence exact computation
over `Q(f)` applies to the actual generic fibre.

For each `r=0,1,2,3`, take the complete universe

```text
b(Y)=b_0+b_1Y+...+b_rY^r.                                 (23)
```

The polynomial `(22)` has degree twelve and leading coefficient `64`.
After changing the sign of a putative square root, write its leading term as
`8Y^6`.  Coefficients in degrees eleven down through six determine the other
six square-root coefficients recursively, dividing only by `16`.  Let
`E_0,...,E_5` be the six residual low coefficients.  Exact grevlex reduction
over `Q(f)`, in variable order `(b_0,...,b_r)`, gives

```text
r=0:   GB(E_0,...,E_5)=[b_0],
r=1:   GB(E_0,...,E_5)=[b_0,b_1],
r=2:   GB(E_0,...,E_5)=[b_0,b_1,b_2],
r=3:   GB(E_0,...,E_5)=[b_0,b_1,b_2,b_3].                (24)
```

This is an exact symbolic universe, not a finite coefficient scan.  As a
transparent hostile edge check, in the linear row its last three residuals
are, up to nonzero scalar normalization,

```text
E_5=-6b_1(b_0-3),
E_4=-3(b_0^2-6b_0+b_1^2f),
E_3=b_1(96b_0f+b_1^2-288f).                              (25)
```

If `b_1!=0`, they force successively `b_0=3`, `b_1^2f=9`, and
`b_1^3=0`.  The contradiction is visible before Groebner reduction.

The common solution of every row in `(24)` is `b=0`, and it is a genuine
positive control:

```text
(4(Y^2-f)^2-3)^3+27(1-(Y^2-f)^2)
 =((Y^2-f)(8(Y^2-f)^2-9))^2.                             (26)
```

But `(17)` requires `b` to have degree exactly `6-n`, hence to be nonzero.
This contradiction closes all four rows in `(18)` and proves `(3)`.

## 5. Exact boundary and next computation

The proof does **not** say that every polynomial residual square is the base
point.  It says that every survivor lies in the three-channel universe `(4)`.
The generic hostile point `(5)` occupies its quadratic channel but fails
`k[x,y]` descent and the origin address.  A global survivor must additionally
satisfy THM-3885's `a=0` arm dichotomy, the `L=0` root polarization, and

```text
T(0,0)=0.                                                  (27)
```

The cheapest decisive successor is therefore the full quadratic-y
coefficient system over `k[x]`, split by the two arm addresses.  THM-3893
already kills its simultaneous-zero sector through total degree six; `(3)`
removes the purported high-y-degree escape but does not bound the remaining
`x` degrees.

Reproduce the exact covariant and square ideals with

```bash
python3 04-computation/jc2_f_zero_quartic_covariant_high_y_emptiness_thm3895.py
python3 -O 04-computation/jc2_f_zero_quartic_covariant_high_y_emptiness_thm3895.py
```

Both streams must byte-match
`05-knowledge/results/jc2_f_zero_quartic_covariant_high_y_emptiness_thm3895.out`.
