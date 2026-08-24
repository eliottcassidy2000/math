---
id: THM-3998
title: "Three-by-at-most-three source-weight support cannot repair the reduced 2:3 cell"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. After THM-3992's
  fixed node gauge, suppose A has exactly the forced source-grading weights
  2,0,-2 and C has weights 3,1 plus at most one arbitrary additional integer
  weight. Coefficient degrees in u are unrestricted. The bracket rows force
  a polynomial conic and scalar ODE whose degree invoices are incompatible;
  every possible collision of the optional C row is either absorbed or
  killed by an endpoint valuation. No such Darboux pair exists. This closes
  a three-by-two/three source-weight cell, not the full reduced 2:3 case.
source: root + bounded_residual + laurent_rows / planar Jacobian continuation, 2026-08-24
audit: >
  PASS (laurent_rows, 2026-08-24). The exact B2 graded pieces, forced
  weight-minus-two jet, core equations, all-degree divisibility contradiction,
  and complete optional-row collision list were independently reconstructed.
  Normal and optimized outputs are byte-identical. Independent finite direct
  Groebner replays provide hostile controls but are not the all-degree proof.
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
related:
  - THM-3979-two-color-formal-cusp-darboux-lifting
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
script: 04-computation/jc2_three_by_three_weight_support_thm3998.py
output: 05-knowledge/results/jc2_three_by_three_weight_support_thm3998.out
script_sha256: 06abe1c0845b282942d45a8a965a91a2b85277ac64fe029e4d114bf134adcc20
output_sha256: 8be9e9838a7e082a53bc64f5be6ee8a803213c9f57f64a594be9f8a30caacdb2
audit_script: 04-computation/jc2_three_by_three_weight_support_groebner_thm3998.py
audit_output: 05-knowledge/results/jc2_three_by_three_weight_support_groebner_thm3998.out
audit_script_sha256: bf561c4ae7712b413deae2c3315a36806103230692625d139cfc33e6f0aed822
audit_output_sha256: a7d633a98ad0dfc588d5c59b039518dcb19add320f1a480e0c4623dd35639c46
hash_basis: raw LF bytes
---

# THM-3998 -- a three-by-at-most-three source-weight obstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** This theorem
closes the stated finite source-grading support cell. It does not claim
`JC(2)`.

## Exact setting and gauge

Work over an algebraically closed field `k` of characteristic zero. In
`B_2` use the corrected canonical relations

```text
u=x^2 t,
p=t(1+u)=t+x^2t^2,
y=xtp,
B_2=k[x,u,p,y] subset k[x,t].
```

Use the THM-3973 source grading

```text
wt(x)=1,       wt(t)=-2,       wt(u)=0,
wt(p)=-2,      wt(y)=-3.
```

Its exact pieces are

```text
(B_2)_r=x^r k[u]                                      (r>=0),
(B_2)_(-h)=x^-h u^ceil(h/2)(1+u)^ceil(h/3) k[u]       (h>=1).
```

Begin **after** THM-3992's determinant-one target scaling, constant target
shear, `A`- and `C`-translations, and centering of the two node addresses.
Thus a hypothetical pair in this gauge has

```text
A(x,0)=gamma^2 x^2+a,
C(x,0)=gamma^3 x^3+(3a gamma/2)x,
A_t(0,0)=-2/(3 gamma a),                 a gamma !=0.
```

The remaining fifth-root action sends

```text
(A,C,a,gamma) -> (zeta^2 A,zeta^3 C,zeta^2 a,zeta gamma),
zeta^5=1.
```

The argument is uniform in `a,gamma`, so this residual action is harmless.
The support condition below is imposed **after** the gauge actions; support
need not be invariant under undoing them.

The nonzero value of `A_t(0,0)` forces a weight `-2` row. Indeed, among the
exact graded pieces, only

```text
x^-2 m(u),       m in u(1+u)k[u],
```

can have a nonzero `t` derivative at `(0,0)`. If
`m=u(1+u)M`, then `A_t(0,0)=M(0)=m'(0)`, so

```text
ord_u(m)=1,       m'(0)=-2/(3 gamma a).
```

## Frozen ansatz

Coefficient degrees in `u` are unrestricted. Freeze only the number of
homogeneous weights:

```text
A=x^2 f(u)+F(u)+x^-2 m(u),
C=x^3 g(u)+x k(u)+x^r ell(u),
```

where the last row may be absent and, if present, `r` is one arbitrary
integer. For `r=1` or `3`, it is already absorbed into `k` or `g`. For
`r>=0`, `ell in k[u]`; for `r=-h<0`, `ell` belongs to the coefficient ideal
in the exact negative-piece formula above. Boundary matching says

```text
f(0)=gamma^2,       F(0)=a,
g(0)=gamma^3,       k(0)=3a gamma/2.
```

Every genuinely new row has `ell(0)=0`: for negative weight this follows
from the module formula, and for nonnegative weight it follows because the
deleted-line polynomial has only weights `3` and `1` in `C`. A scalar
weight-zero constant was already removed by the `C`-translation.

Then no rows in this ansatz satisfy `J_(x,t)(A,C)=1`.

This is the first complete wider-support cell beyond THM-3973: that theorem
closes pairs with at most two weights in each entry, whereas the forced
THM-3992 row has three nonzero `A` weights and two or three `C` weights.

## Core cell: no optional row

For homogeneous rows, direct differentiation gives

```text
J(x^r alpha(u),x^s beta(u))
 =x^(r+s+1)(r alpha beta'-s alpha' beta).               (1)
```

The core bracket weights are `6,4,2,0`. The weight-six equation is

```text
2fg'-3f'g=0.                                            (2)
```

Hence `(g^2/f^3)'=0`. Its value at `u=0` is one, so `g^2=f^3`.
Put `q=g/f` in `k(u)`. Then `q^2=f`, `q^3=g`; since `q^2` is polynomial,
normality of `k[u]` gives `q in k[u]`, and `q(0)=gamma`. Thus

```text
f=q^2,                 g=q^3.                           (3)
```

The weight-four equation is

```text
2fk'-f'k-3F'g=0.
```

After `(3)`, it says `(k/q-3F/2)'=0`. The boundary values kill the constant:

```text
k=(3/2)qF.                                               (4)
```

The weight-two equation becomes

```text
F F'+4m q q'+2m'q^2=0,
```

so its boundary value integrates to

```text
F^2+4m q^2=a^2.                                         (5)
```

Finally the scalar row is

```text
-2mk'-m'k=1.                                            (6)
```

Differentiate `(5)`, substitute `(4)`, and eliminate `m,m'` from `(6)`:

```text
3F'(2F^2-a^2)=4q.                                       (7)
```

If `d=deg(F)>=1`, equation `(7)` gives

```text
deg(q)=3d-1.                                             (8)
```

But `(5)` and polynomiality of `m` give

```text
q^2 divides a^2-F^2,
```

which would require `2(3d-1)<=2d`. The gap `4d-2` is positive for every
`d>=1`. If `F` is constant, `(7)` has zero left side and nonzero right side.
Both cases are impossible.

## One optional `C` weight

The optional row contributes at weights

```text
r+3, r+1, r-1.
```

These can collide with a core weight in `{6,4,2,0}` only when

```text
r in {-3,-1,1,3,5,7}.                                  (9)
```

If there is no collision, all four contradictory core equations remain
unchanged. The values `r=1,3` merely enlarge the already arbitrary
polynomials `k,g`.

For `r=-3,-1`, the unique lowest row is

```text
-2m ell'-r m'ell=0.                                    (10)
```

Since `ord_u(m)=1`, a nonzero `ell` of integer valuation `v` would make
the leading coefficient proportional to `2v+r`. Equation `(10)` would
require `v=-r/2`, a half-integer in both cases. Hence `ell=0`.

For `r=5,7`, the unique highest row is

```text
2f ell'-r f'ell=0.                                     (11)
```

Here `f(0)=gamma^2!=0` and `ell(0)=0`. If `ord_u(ell)=v>=1`, the leading
term of `(11)` is `2 gamma^2 v` times `u^(v-1)`, nonzero in characteristic
zero. Hence `ell=0` again. This exhausts `(9)`, so one arbitrary extra
`C` weight never repairs the core obstruction.

## What “support” means here

This is a **source-grading weight-support** theorem, not a bounded
`tau`-Laurent-depth theorem. In THM-3989's variables

```text
s=xt,       tau=t,       x=s/tau,       u=s^2/tau.
```

Therefore

```text
x^r f(u)=(s/tau)^r f(s^2/tau).
```

If `deg(f)=D`, this row can occupy `D+1` Laurent powers and reach depth
growing with `D`. Since all coefficient degrees above are unbounded, the
Laurent support and cusp-pole depth are genuinely unbounded. The result
closes a finite **number of weights**, uniformly in polynomial degree and
Laurent depth.

The mixed residual

```text
C^2-A^3+(3a^2/4)A+a^3/4
 =gamma u+(3a/(2gamma))p+R(p,y),       R in (p^2,y),
```

is not truncated or specialized. The proof uses only its already proved
node gauge and forced `A_t` jet. Thus it applies for every `R`, but it does
not extract additional restrictions from the coefficients of `R`.

## Controls, reproduction, and boundary of the result

The certificate verifies:

1. `J(x,t)=1` as an ambient positive sign control;
2. THM-3973's rational pair
   `P=u(1+u)/(1+2u)^2`, `Q=(1+2u)^3/x` has `J(P,Q)=1`;
3. `p=x^-2u(1+u)=t(1+u)` passes the exact weight-minus-two module gate,
   while bare `t=x^-2u` does not;
4. the nodal boundary parametrization satisfies its cubic eliminant but has
   bracket zero;
5. `A=A_node-2p/(3gamma a), C=C_node` populates the frozen support cell and
   meets both the boundary and forced-jet filters, but is hostile rather than
   Keller.

Run from the isolated worktree root:

```text
python3 04-computation/jc2_three_by_three_weight_support_thm3998.py
python3 -O 04-computation/jc2_three_by_three_weight_support_thm3998.py
python3 04-computation/jc2_three_by_three_weight_support_groebner_thm3998.py
python3 -O 04-computation/jc2_three_by_three_weight_support_groebner_thm3998.py
```

The normal and optimized executions are byte-identical. The first certificate
ends with `ALL THM-3998 EXACT CHECKS PASSED`; the second gives four independent
direct-source Groebner hostiles and ends with
`ALL THM-3998 DIRECT GROEBNER REPLAYS PASSED`.

Open outside this result:

```text
* any additional weight in A;
* two or more genuinely additional weights in C;
* a bounded Laurent-support search independent of the weight support;
* vanishing or further classification of R;
* the full reduced (2,3) cell and JC(2).
```
