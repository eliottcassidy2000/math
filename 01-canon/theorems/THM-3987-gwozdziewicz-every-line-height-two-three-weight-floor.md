---
id: THM-3987
title: "Gwozdziewicz every-line gate and the height-two three-weight floor"
status: >
  PROVED FROM CITED INPUTS + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Any polynomial Darboux pair in the height-two completion B_2 would
  give a nonautomorphic planar Keller map, hence is immersed and noninjective
  on every affine source line. The two coordinate lines force, in each
  output, a surviving positive weight at least 2, a surviving even negative
  weight at most -4, and a distinct linear-jet weight 1 or -2. Across the
  pair one positive endpoint is at least 3 and one negative endpoint is at
  most -6; exact three-weight rows have opposite middle jets. Thus every
  output has at least three retained weights. In ordinary source degree, each
  output's degree-at-least-two tail and the pair's combined
  degree-at-least-three tail have no common projective direction zero.
  Together with THM-3974's 3x3 obstruction, all 2xm/mx2 cells disappear and
  the first live total-support cells are only 3x4 and 4x3. No Darboux pair or
  JC(2) counterexample is constructed or excluded in unrestricted support.
source: jc-zero-debt-lift / post-THM-3985 every-line support lane, 2026-08-24
audit: >
  PASS (root / jc-cohn3709, 2026-08-24). The audit independently checked
  the nonautomorphism-to-every-line contrapositive, the boundary valuation
  ord_D(t)=-2, the immersed-curve degree lemma, both axis-survival ledgers,
  odd-negative disappearance, the complete origin jet, the exact-three
  opposite-jet seam, and the projective directional-tail reformulation.
  Normal, optimized, and frozen outputs byte-match at CHECKS=511; hashes and
  documentation agree.
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
  - THM-3974-height-tower-few-weight-darboux-support-obstruction
related:
  - THM-3589-danielewski-central-arm-every-line-and-kummer-trace-darboux-gates
  - THM-3695-y0-collision-ring-danielewski-embedding-and-seven-piece-floor
  - THM-3983-coordinate-boundary-constancy-and-rational-place-budget
external:
  - "Gwozdziewicz, Injectivity on one line, arXiv:alg-geom/9305008, Theorem 1.1."
script: 04-computation/jc2_gwozdziewicz_every_line_weight_floor_thm3987.py
output: 05-knowledge/results/jc2_gwozdziewicz_every_line_weight_floor_thm3987.out
script_sha256: 14d43bcd6ded0a750ba2e471bc00656572ab714593d0793a0f8b765f82c1061c
output_sha256: 292c5de53e4ae164542b20edc8c65fb97bb190877ffbad82f0a6e4809e5dcb3f
semantic_sha256: a7d7d96e3874e1c811cc3d34a612900130d501bf11638a15f0969c6e1fe5726b
hash_basis: raw LF bytes
---

# THM-3987 -- every source line invoices both tails and a transverse jet

**PROVED FROM CITED INPUTS + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.** Work over an algebraically closed field `k` of
characteristic zero.
Inside `k[x,t]`, put

```text
u=x^2t,             z=1+u,
p=zt,               y=xzt^2,
B_2=k[x,z,p,y].                                             (1)
```

Give `B_2` the grading

```text
wt(x)=1,             wt(t)=-2,             wt(u)=wt(z)=0. (2)
```

Let `A,C in B_2` and suppose

```text
J_(x,t)(A,C)=1.                                            (3)
```

Delete only scalar weight-zero summands, which bracket to zero, and call all
remaining nonzero homogeneous pieces *retained*. Then every one of the two
outputs `P in {A,C}` has three distinct compulsory pieces:

```text
(positive tail)  x^r f_r(u),      r>=2,      f_r(0)!=0;
(negative tail)  weight -2j,      j>=2,      normalized coefficient at 0 !=0;
(linear jet)     weight 1 or -2,              nonzero linear coefficient. (4)
```

Moreover, if `r_A,r_C` and `j_A,j_C` denote the largest surviving degrees on
the two coordinate lines, then

```text
max(r_A,r_C)>=3,                 max(j_A,j_C)>=3.          (5)
```

Equivalently, one positive support endpoint across the pair is at least `3`,
and one negative support endpoint is at most `-6`. Consequently

```text
#supp_w(A)>=3,                    #supp_w(C)>=3.           (6)
```

If both supports have exactly three retained weights, they have the forms

```text
{r_A, epsilon_A, -2j_A},         {r_C, epsilon_C, -2j_C},
r_A,r_C>=2,   j_A,j_C>=2,        {epsilon_A,epsilon_C}={1,-2}. (7)
```

Finally, THM-3974 excludes the remaining `3x3` cell. Hence

```text
no 2xm or mx2 Darboux pair exists for any finite m;
#supp_w(A)+#supp_w(C)>=7;
the first live support cells are exactly 3x4 and 4x3.      (8)
```

Here *live* means not excluded by the present theorem and its stated
dependency; it does not assert existence.

## 1. Why the source Keller map cannot be an automorphism

The inclusion `(1)` makes

```text
H=(A,C):A2_(x,t) -> A2
```

a polynomial Keller map. It is not a polynomial automorphism. Indeed, near
the boundary `D=V(x,z)` of `X_2=Spec(B_2)`, the `D(z-1)` chart of THM-3973 is

```text
z(z-1)^2=x^3y.                                            (9)
```

At the generic point of `D`, `x` is a uniformizer, `z-1` is a unit, and

```text
t=(z-1)/x^2,                     ord_D(t)=-2.             (10)
```

Thus `t notin B_2`. If `H` were an automorphism, then

```text
k[x,t]=k[A,C] subset B_2 subset k[x,t],                  (11)
```

contradicting `(10)`. This is the whole nonautomorphism argument; finite
degree, properness, and an assumed collision are not being smuggled in.

Gwozdziewicz's cited Theorem 1.1 is stated over any algebraically closed
characteristic-zero field: a planar Keller map injective on one affine line
is a polynomial automorphism. Its contrapositive and `(11)` show that `H`
is noninjective on **every** affine source line. It is immersed on every such
line, because the invertible differential in `(3)` cannot kill a nonzero
line tangent.

## 2. The elementary immersed-curve invoice

We use the following sharp lemma. If

```text
gamma=(f,g):A1 -> A2
```

is an immersed noninjective polynomial curve, then

```text
deg(f)>=2,             deg(g)>=2,             max(deg(f),deg(g))>=3. (12)
```

A degree-one coordinate is itself injective. If a coordinate is constant,
immersion makes the derivative of the other coordinate nowhere zero; over
an algebraically closed characteristic-zero field that derivative is a
nonzero scalar, again making the curve injective. Thus both degrees are at
least two. If both are quadratic and `gamma(a)=gamma(b)` with `a!=b`, then
each quadratic derivative vanishes at `(a+b)/2`, contradicting immersion.
The nodal parameterization

```text
s |-> (s^2-1, s(s^2-1))                                  (13)
```

shows the degree pair `(2,3)` is sharp.

There is also a global directional form of `(12)`. Write the ordinary source
homogeneous expansions

```text
A=sum_(d>=0) A_d(x,t),             C=sum_(d>=0) C_d(x,t).
```

For every direction `[v] in P1`, apply `(12)` to the line `s |-> sv`. Then

```text
some A_d(v)!=0 with d>=2,          some C_d(v)!=0 with d>=2,
some A_d(v) or C_d(v)!=0 with d>=3.                         (12a)
```

Equivalently, each of the two finite form families

```text
{A_d:d>=2},              {C_d:d>=2}
```

has no common zero on `P1`, and the combined family

```text
{A_d,C_d:d>=3}                                             (12b)
```

has no common zero there. This is a necessary directional Newton gate, not a
sufficiency statement and not a claim that the top homogeneous form alone is
nowhere zero.

## 3. Exact survival on the two coordinate lines

THM-3973 gives the homogeneous modules

```text
(B_2)_r=x^r k[u]                                      for r>=0,
(B_2)_(-q)=x^(-q)u^ceil(q/2)(u+1)^ceil(q/3)k[u]       for q>=1. (14)
```

On the source line `t=0`, every negative piece vanishes, a weight-zero piece
restricts to a scalar, and

```text
[x^r f_r(u)]|_(t=0)=f_r(0)x^r.                          (15)
```

Apply `(12)` to `(A,C)|_(t=0)`. Each output therefore has a surviving
positive weight `r>=2`, and at least one of the two positive endpoints has
`r>=3`.

On `x=0`, every positive piece vanishes and a weight-zero piece again becomes
a scalar. Put `q=2j` or `2j-1` in the negative row of `(14)`. Substitution
`u=x^2t` gives

```text
[weight -2j]|_(x=0)=c*t^j       if its normalized coefficient has c!=0,
[weight -(2j-1)]|_(x=0)=0.                              (16)
```

Thus only even negative weights survive this line. Applying `(12)` to
`(A,C)|_(x=0)` forces, in each output, a surviving weight `-2j<=-4`, and at
least one output has `j>=3`, hence endpoint at most `-6`. This proves both
tails in `(4)` and all of `(5)`.

## 4. The origin jet supplies the third piece

At `(x,t)=(0,0)`, a homogeneous source polynomial can have a nonzero linear
term only in the two weights

```text
weight 1:       x f(u)=f(0)x+higher terms,
weight -2:      p g(u)=g(0)t+higher terms.               (17)
```

Every other weight has source order at least two after its scalar part is
removed. The Jacobian condition `(3)` says the two gradient rows are nonzero
and independent. Therefore each output contains at least one of the jet
pieces in `(17)`, and the pair collectively contains both jet directions.
Neither jet weight can equal the positive or negative tail already forced in
Section 3, proving `(6)`.

If an output has exactly three retained weights, the two tails use two slots,
so exactly one slot remains and it must be its unique jet weight. If both
outputs have size three, independence of their linear parts forces those
unique jets to be opposite, proving `(7)`.

## 5. Exact cusp log chart

There is a complementary denominator-depth form of `(3)`. In the function
field put

```text
s=y/p=xt,                         tau=H/p^2=t.             (18)
```

Then

```text
p=s^2+tau,        y=s(s^2+tau),        x=s/tau,
u=s^2/tau,        dx wedge dt=ds wedge dtau/tau.          (19)
```

Consequently every hypothetical pair satisfies the exact log-Darboux row

```text
tau(A_s C_tau-A_tau C_s)=1.                               (20)
```

This identity is recorded as a bridge to denominator-depth analysis. It is
not used to infer the weight floor, and no regularity at `tau=0` is claimed
for the individual expressions `s/tau` or `s^2/tau`.

## 6. Support closure, scope, and reproduction

The coordinatewise floor `(6)` closes every infinite `2xm` and `mx2` arm;
it is not inferred from a finite rectangle. This is exactly the logical
repair demanded by MISTAKE-442. THM-3974 independently closes `3x3` at height
two, so `(8)` follows.

The theorem concerns polynomial functions in `B_2` whose source bracket is
one. It does not extend the line restrictions to arbitrary rational
functions, does not claim every three-weight support is realizable, and does
not close `3x4` or larger cells. Gwozdziewicz is the sole cited input; the
nonautomorphism, immersion, line-degree, survival, and jet arguments are
proved here.

The companion verifies the homogeneous restriction formulas through hostile
weight ranges, the quadratic-midpoint and nodal controls, the directional
ordinary-tail implications, the complete linear jet matrix, the cusp log
identity, and the support-cell enumeration. These computations illustrate
the uniform proof rather than replace it.

Reproduce with

```bash
python3 04-computation/jc2_gwozdziewicz_every_line_weight_floor_thm3987.py
python3 -O 04-computation/jc2_gwozdziewicz_every_line_weight_floor_thm3987.py
sha256sum 04-computation/jc2_gwozdziewicz_every_line_weight_floor_thm3987.py \
  05-knowledge/results/jc2_gwozdziewicz_every_line_weight_floor_thm3987.out
python3 agents/check_docs.py
```

Independent hostile audit is requested especially for the contrapositive use
of Gwozdziewicz, the `t notin B_2` step, odd-negative disappearance on
`x=0`, and the exact-three opposite-jet conclusion.
