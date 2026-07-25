---
id: THM-2240
title: "The DC2 grade-response gauge is not a continuation state"
status: >
  PROVED + VERIFIED-EXACT. In the THM-2049 Ore boundary problem, equality
  of the current grade-six response is not even a one-step right-continuation
  state: the homogeneous syzygy C=1 leaves the grade-six cancellation
  unchanged but changes the normalized grade-seven residual by
  8q(4u^2-13u+13). More strongly, the q-integration-constant axes are
  infinite-dimensional and separately inject into the next residual, while
  the combined two-axis map is not injective. The complete pure grade-six
  affine correction fiber has an exact arbitrary-q response formula, and no
  member skips grade seven without a next-rung attachment. This is a kernel-
  memory and pure-rung theorem, not nontermination: THM-2049's
  associated-graded surjectivity and formal beta-adic lift remain unchanged.
source: codex-2026-07-25-dc2-continuation-gauge
depends_on:
  - THM-2049-the-DC2-Ore-boundary-correction-complex-is-acyclic
related:
  - THM-2044-explicit-rank-two-poisson-counterexample-by-symplectic-suspension
  - THM-2046-first-order-cotangent-pullbacks-cannot-cross-the-DC2-wall
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
script: 04-computation/dc2_grade_response_continuation_defect_thm2240.py
output: 05-knowledge/results/dc2_grade_response_continuation_defect_thm2240.out
script_sha256: f562de13b4c0d9c60b3fc52b5cc437b262290ceeb2dc3145999529f3ee445278
output_sha256: 4a8e51e87101aec7a1de7ecd2b03ed8191fcca21e7d3e7df8ec52347602109c9
hash_basis: working-tree bytes (LF)
---

# THM-2240 -- the grade response forgets its own continuation

Work in the coefficient-left Ore algebra

```text
O=Q[x,q][ell;delta],
delta=3x^2 partial_x+(2-6xq) partial_q.              (1)
```

This theorem uses the Weyl-ordered pair `(S_W,T_W)` of THM-2049. It does not
alter that theorem's formal lift. It determines exactly what its grade-six
response quotient forgets one grade later.

## 1. Filtration and symbol-lift convention

For a nonzero coefficient-left operator, put

```text
beta(sum_k a_k(x,q)ell^k)=min_k(v_x(a_k)-2k).        (2)
```

Let `u` be the commutative associated-graded variable represented by
`x^2 ell`. For `F(q,u)=sum_k f_k(q)u^k`, define the coefficient-left symbol
lift

```text
iota_p(F)=sum_k f_k(q)x^(p+2k)ell^k.                 (3)
```

Thus (3), not literal noncommutative evaluation at `(x^2 ell)^k`, is meant
whenever a correction is displayed using `u`.

For an operator `R` and an integer `g`, define its normalized grade-`g`
symbol by

```text
H_g(R)
 =sum_(k,j) [x^(g+2k)q^j ell^k]R q^j u^k.           (4)
```

Equivalently, in Rees notation,

```text
in_g(R)=x^g H_g(R).                                  (5)
```

Put

```text
D(u)=2u^2-10u+9,
K(u)=4u^2-13u+13.                                   (6)
```

The exact Weyl residual

```text
R_W=[S_W,T_W]-1                                     (7)
```

begins in grade six, with

```text
H_6(R_W)=2(u-1)(u-2).                                (8)
```

## 2. The current response and its complete pure-rung fiber

A pure grade-six correction is a pair

```text
(iota_5(Atilde), iota_6(Btilde)),
Atilde,Btilde in Q[q,u].                             (9)
```

THM-2049 computes its grade-six response as

```text
d_6(Atilde,Btilde)
 =(8/3)(u-2)partial_q Atilde
   +(D/9)partial_q Btilde.                           (10)
```

The standard section is

```text
Atilde_0=-(3/4)q(u-1),
Btilde_0=0.                                         (11)
```

It has response `-2(u-1)(u-2)` and therefore kills (8). Its exact residual
begins in grade seven:

```text
H_7^0=-(3/2)q(10u^2-36u+29).                        (12)
```

Every pure correction in the same grade-six cancellation fiber, and no
others, has the following form. Choose

```text
C(q,u) in Q[q,u],
J(q,u)=integral_0^q C(z,u) dz,
a(u),b(u) in Q[u],                                  (13)
```

and set

```text
Atilde=Atilde_0+D J+a,
Btilde=-24(u-2)J+b.                                 (14)
```

Indeed, differentiating (14) in `q` adds

```text
(D C,-24(u-2)C),                                    (15)
```

which is the full syzygy kernel in THM-2049:

```text
(8/3)(u-2)DC+(D/9)(-24(u-2)C)=0.                    (16)
```

Conversely, coprimality of `D` and `u-2` gives (15), and integrating in `q`
gives the two and only two integration constants `a(u),b(u)`. The lower
limit in (13) fixes the otherwise redundant constant part of `J`.

## 3. Exact grade-seven response of the full fiber

For the correction (14), the exact next symbol is

```text
H_7
 =H_7^0+L_J(J)+L_a(a)+L_b(b),                       (17)
```

where

```text
L_J(J)=4K(J+q partial_q J)+18 partial_u J,           (18)

L_a(a)=20(u-2)a
       +2(2u^2-6u+1)a',                             (19)

L_b(b)=D b+[u(u-2)(u-4)/3]b'.                       (20)
```

These are exact Ore identities, not a truncated commutative approximation.

### Derivation

The first two homogeneous layers of the Weyl pair are, in Rees-symbol
notation,

```text
T_(-1)=x^(-1)t(u),       t=(2/3)u(4-u),
T_0=qv(u),               v=(u-1)^2,                 (21)

S_(-2)=x^(-2)s(u),       s=u(u-3)(2u-9)/54,
S_(-1)=x^(-1)q r(u),     r=-u(u^2-3u-3)/18.         (22)
```

Split the coefficient derivation as

```text
delta=delta_0+delta_1,
delta_0=2 partial_q,
delta_1=3x^2 partial_x-6xq partial_q.                (23)
```

Through grade seven, only the `r=1` term of the Ore product contributes.
The `delta_0` bracket of the leading layers is the grade-six response;
the `delta_1` bracket of the leading layers and the `delta_0` bracket with
the next layers give grade seven. Terms with at least two applications of
`delta` begin in grade eight.

For a general grade-five `S` correction `iota_5(F)`, these grade-seven
pieces sum to

```text
M_S(F)
 =(2v-3t)F_u-15t'F
   +q(6t'-2v')F_q.                                  (24)
```

For a general grade-six `T` correction `iota_6(G)`, they sum to

```text
M_T(G)
 =18s'G+(6s-2r)G_u
   +q(2r'-6s')G_q.                                  (25)
```

Putting

```text
F=DJ,                    G=-24(u-2)J                 (26)
```

in (24)--(25) and simplifying gives (18). Setting `F=a,G=0` gives
(19), while `F=0,G=b` gives (20).

No nonlinear correction term can enter (17). The two correction components
have beta grades five and six, so their mutual commutator has beta at least

```text
5+6+2=13.                                           (27)
```

The same bound applies to the standard grade-five correction paired with a
grade-six gauge term.

## 4. An explicit splitter at the very next grade

Take the homogeneous syzygy

```text
C=1,                    J=q,                 a=b=0.  (28)
```

Relative to the standard section, this adds

```text
eta_S=iota_5(qD),
eta_T=iota_6(-24q(u-2)).                             (29)
```

The two grade-six responses are

```text
H_6([eta_S,T_W])=(8/3)(u-2)D,
H_6([S_W,eta_T])=-(8/3)(u-2)D,                      (30)
```

so they cancel. At grade seven, however,

```text
H_7([eta_S,T_W])
 =8q(2u-5)(2u^2-8u+5),                              (31)

H_7([S_W,eta_T])
 =-8q(u-2)(4u^2-22u+19).                            (32)
```

Their sum is the nonzero polynomial

```text
8qK=8q(4u^2-13u+13).                                (33)
```

Thus the gauged correction still begins in grade seven, but with

```text
H_7^1=q(34u^2-100u+121)/2,                          (34)

H_7^1-H_7^0=8qK.                                    (35)
```

In particular, if two pure grade-six corrections are identified when their
current responses under `d_6` agree, then the exact next-residual map is not
well defined on the quotient:

```text
ker(d_6) is not contained in ker(Delta H_7).         (36)
```

This is failure of a one-step, right-continuation state. It is not a claim
about a two-sided congruence.

The splitter is stable under a common next continuation. If the same
next-rung correction

```text
(iota_6(P),iota_7(Q))                               (37)
```

is added to both representatives, it adds the same `d_7(P,Q)` to both
grade-seven residuals. Cross-commutators between (9) and (37) begin in beta
grade fourteen, so their difference at grade seven remains (33). A
state-dependent next correction can respond differently only if the
forgotten sidecar has first been retained.

## 5. Infinite-dimensional memory, with the combined-kernel guardrail

The integration-constant axes alone make the lost memory
infinite-dimensional.

If `a` is nonzero of degree `n>=0` with leading coefficient `A`, then

```text
deg L_a(a)=n+1,
LC(L_a(a))=4(n+5)A!=0.                              (38)
```

For `n=0`, the derivative term vanishes and the same formula gives `20A`;
the zero polynomial maps to zero. Hence `L_a` is injective on `Q[u]`.

If `b` is nonzero of degree `n>=0` with leading coefficient `B`, then

```text
deg L_b(b)=n+2,
LC(L_b(b))=((n+6)/3)B!=0.                           (39)
```

Again the `n=0` derivative term vanishes and gives leading coefficient `2B`.
Thus `L_b` is separately injective.

The word **separately** is load-bearing. The combined map

```text
(a,b) -> L_a(a)+L_b(b)                              (40)
```

is not injective. The exact nonzero witness

```text
a=-146+140u-25u^2,
b=-680+300u                                         (41)
```

satisfies

```text
L_a(a)=-20(35u^3-248u^2+515u-306),
L_b(b)= 20(35u^3-248u^2+515u-306).                  (42)
```

No claim of combined `a,b` injectivity is made.

The arbitrary-`q` syzygy axis is also injective. Write

```text
C(q,u)=sum_(r>=0) q^r c_r(u),
J(q,u)=sum_(r>=0) q^(r+1)c_r(u)/(r+1).              (43)
```

Equation (18) has independent `q^(r+1)` layers

```text
[4(r+2)K c_r+18c_r']/(r+1).                         (44)
```

For nonzero `c_r` of degree `n`, the numerator in (44) has degree `n+2`
and leading coefficient

```text
16(r+2)LC(c_r)!=0.                                  (45)
```

Equivalently, the actual layer in (44) has leading coefficient
`16(r+2)LC(c_r)/(r+1)`. Hence every layer operator, and therefore the full
zero-constant `J` axis, is injective.

## 6. No pure grade-six correction skips grade seven

The full parameterization (13)--(14) makes a stronger scoped statement
possible:

> Every pure grade-six correction which kills (8) has exact residual beta
> equal to seven.

Suppose otherwise that (17) vanished. The terms `L_a(a)+L_b(b)` have
`q`-degree zero, so they cannot affect the positive `q` layers. For every
`r>=1`, the `q^(r+1)` layer comes only from (44). Equations (44)--(45) force

```text
c_r=0                    for every r>=1.             (46)
```

At `q^1`, a nonconstant `c_0` would make (44) have `u`-degree at least three,
which cannot cancel the quadratic coefficient of (12). If `c_0=lambda` is
constant, cancellation of the `u^2` coefficient would force

```text
lambda=15/32.                                        (47)
```

The `u` coefficient would then remain

```text
54-104(15/32)=21/4!=0.                              (48)
```

Thus the `q^1` layer cannot vanish. This proves the claim.

The word **pure** is essential. Section 6 excludes only a correction
consisting of the grade-five/grade-six pair (9). One may simultaneously add
a next-rung pair of beta grades six and seven as in (37), and THM-2049's
surjective `d_7` then kills the exposed grade-seven residual. The theorem
does not block that recursion.

## 7. Consequences and scope

THM-2049 proves that every current homogeneous residual has a correction and
therefore constructs a formal beta-adic lift. THM-2240 proves a different
fact: a current-response class does not determine the exact next residual.
The two statements are compatible because surjectivity neither removes nor
controls the kernel.

Consequently:

1. the particular correction ladder `6 -> ... -> 14` in THM-2049 is one
   section, not an invariant nontermination certificate;
2. a termination search must retain a representative-sensitive sidecar,
   such as the exact future residual profile or equivalent gauge data;
3. proving nontermination would require controlling every gauge sequence,
   while proving termination may use state-dependent higher-rung
   corrections; and
4. nothing here proves or refutes polynomial termination, supplies the
   missing `D`-column, proves `DC(2)`, or affects the exact associated-graded
   acyclicity.

The theorem is also deliberately separate from THM-2230. In the planar
Jacobian response fiber, the kernel is an exact target-shear orbit and the
response quotient is complete for that operation. Here, the Ore
associated-graded response has kernel elements which exact continuation
detects one grade later.

## 8. Exact companion

Run

```bash
python3 04-computation/dc2_grade_response_continuation_defect_thm2240.py
python3 -O 04-computation/dc2_grade_response_continuation_defect_thm2240.py
```

The standalone script independently rebuilds coefficient-left Ore
multiplication and the Weyl pair. It verifies (8), (12), the four component
identities (30)--(32), the arbitrary-`q` formula (18), the axiswise leading
terms (38)--(45), the combined-kernel witness (41), and the no-skip
calculation (47)--(48). A mixed-degree hostile polynomial is checked both by
the full Ore product and by (17), and a common next-rung correction preserves
the splitter. All arithmetic is exact, every gate uses an explicit raising
check rather than `assert`, and normal and optimized runs match the frozen
output.

QED.
