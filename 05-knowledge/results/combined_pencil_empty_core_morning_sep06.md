# A free relative scale defeats the combined midpoint abstraction

Status: **PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The model below retains both Euler identities, full nonnegative auxiliary
support, all-weight square-pencil compatibility, the common binomial
carrier, and its original joint first-row zero. Nevertheless a rational
instance has a positive doubled-model value, and a precisely specified
positive algebraic scale makes both model rows vanish. This is **not an
actual binomial-path or Laurent-support counterexample**. Actual
trinomial signed transport remains **OPEN**.

## 1. Inheritance and the combined target

The incoming [whole midpoint classes](nc2_midpoint_classes_overnight_hexagon_sep05.md)
and its direct [weighted-deletion](nc2_weighted_midpoint_overnight_hexagon_sep05.md)
and [Hadamard-transport](nc2_hadamard_transport_overnight_hexagon_sep05.md)
routes were read before this probe. For A2B3 they identify whole classes
`R_0=G1+G3` and `K_(N,1)=K_(N,2)=G2/2`. Alpha completion instead controls
`W=V+G1`; the two groupings must not be conflated.

The earlier [Euler-coupled hostile](laurent_transport_empty_core_next_sep06.md)
keeps the ordinary carrier and contiguous relations but fails weighted
compatibility. The incoming cubic Hadamard hostile keeps weighted
compatibility but lacks that carrier and those relations. This probe
keeps both collections. The least-used sidecar is the relative
normalization of the crossing factor `2t` against the squared hit term.
The live board is: full beta support; weighted square pencils; Euler
coupling; joint first zero; whole-class signs; relative phase scale.
No fixed-degree actual-return theorem is reproved and no external
priority is asserted.

## 2. An exact all-weight admissible model

Fix `q=5,h=4,x=1,g=q+2h+1=14,k=2h+1=9`. Put

```text
a_i=64^i, 0<=i<=4,
B(v)=product_i(v-a_i)=sum_(i=0)^5 b_i v^i,
c_i=3(i+1)b_(i+1)/(q+2+2i),
d_i=(2+3i)c_i/(q+1+2i),           0<=i<=4,
C(v)=sum c_i v^i,   D(v)=sum d_i v^i.
```

B,C,D are monic. Direct rational evaluation gives
`C(a_i)/B'(a_i)>0` and `D(a_i)/B'(a_i)>0` at every a_i; the full rational
ratios are frozen in the output. These signs force one root of each C,D
in every interval `(a_i,a_(i+1))`, exhausting their degree four. Thus both
strictly interlace B and have simple positive roots.

For **every** `w>0`, the polynomial

```text
L_w(v)=wB(v)^2-2vC(v)D(v)                              (1)
```

has ten simple positive roots. Its sign is positive at0, negative at each
B root, positive at each intervening C root, and positive for sufficiently
large v. These ten sign changes exhaust degree ten. At **w=0** its
degree drops to nine and its roots are0 and the roots of C,D, with any
multiplicities retained. No sampled-weight argument is used.

Define positive-coefficient polynomials by reversal,

```text
G_B(t)=product_i(1+a_i t),
G_C(t)=sum_(j=0)^4 (-1)^j c_(4-j)t^j,
G_D(t)=sum_(j=0)^4 (-1)^j d_(4-j)t^j.
```

They have constant coefficient one and negative roots. Replacing
`t=-1/v` shows that (1) is exactly the all-weight statement for
`wG_B^2+2tG_CG_D`, whose coefficients are nonnegative for `w>=0`.
The alpha factor below is the actual odd binomial row, so its analogous
all-weight compatibility is already covered by the incoming weighted
network theorem.

Use the complete Laurent factors and responses

```text
beta=t^-1 G_B,   C_raw=t^-1 G_C,   D_raw=t^-1 G_D,
O(t)=sum_j binom(14,2j+1)t^j,
E(t)=sum_j binom(14,2j)t^j,
A_double=O^2+t^-1 E^2,
P=O star beta,
V=O^2 star beta^2,
G1=(t^-1 E^2) star beta^2,
G2=2 O^2 star(t C_raw D_raw),
G3=2(t^-1 E^2) star(t C_raw D_raw),
W=V+G1,
R_skip=G2+G3,
Q_model=W+R_skip,
R_0=G1+G3,        K=G2/2.                               (2)
```

All lower indices are retained before the ordinary coefficientwise
products. In particular the doubled lower carry is present. The whole
class factorization and finite Hadamard root theorem apply to this
model's full PF inputs, but give no root-location comparison with P.
Explicit common shifts keep this invocation ordinary: `t^2 Q_model` is
the Hadamard product of `t(E^2+tO^2)` and `G_B^2+2tG_CG_D`;
`t^2 R_0` instead has first input `tE^2` and the same second input.
For the crossing class, `tK=(tO^2) star(G_CG_D)`.

At phase `t=-s<0`, the ordinary kernels and carriers are

```text
K_B=s^5 B(u^2/s),
K_C=s^4 C(u^2/s),   K_D=s^4 D(u^2/s),
H_i=(1+u)^14 K_i.
```

They satisfy the exact incoming identities

```text
K_B'=(2/3)u[(q+2)K_C+uK_C'],
(q+1)K_D+uK_D'=2K_C+(3/2)uK_C'.                        (3)
```

These follow coefficientwise from the definitions of c_i,d_i. The same
original zero is retained:

```text
[u^9]H_B=22s p(s),
p(s)=-10845877+73940768626816s
 -1665737635624910848s^2+182710607798008283136s^3
 -104915856919223074816s^4.                              (4)
```

## 3. Exact sign certificate and the scale operation

Let `lambda_d=2^-25` and

```text
I=(4938/1000,4939/1000),      J=lambda_d I.
```

Exact endpoint signs and a positive derivative on `[0,sup J]` prove
that p has a unique root s_* in J and that it is its smallest positive
root. Exact rational interval arithmetic on J proves, at `rho=-s_*`,

```text
V,W,G1,R_0,Q_model-V,Q_model <0,
G3,K,R_skip >0.                                        (5)
```

The original zero therefore coexists with positive whole-crossing and
combined-skip responses, even after both collections of constraints have
been retained. The stronger full response still survives at this scale.

Now introduce **any positive real** lambda:

```text
B_lambda(v)=lambda^q B(v/lambda),
C_lambda(v)=lambda^(q-1)C(v/lambda),
D_lambda(v)=lambda^(q-1)D(v/lambda).                     (6)
```

Both Euler identities and strict interlacing are preserved. Consequently
all-weight compatibility is preserved as well. At the retuned phase
`rho/lambda`, each ordinary K_i and H_i is literally unchanged. Thus the
same common binomial carrier and its original zero persist.

However the raw responses scale differently. Direct coefficient
comparison with all Laurent shifts retained gives

```text
P_lambda(rho/lambda)=lambda^x P(rho),
V_lambda(rho/lambda)=lambda^(2x)V(rho),
W_lambda(rho/lambda)=lambda^(2x)W(rho),
R_skip,lambda(rho/lambda)=lambda^(2x-1)R_skip(rho).       (7)
```

For an explicit ordinary-carrier derivation, `K_(i,lambda)(s,u)=K_i(lambda s,u)`
and likewise for H. The unscaled raw extraction identities are

```text
P(-s)=(-s)^(-x)[u^k]H_B,
W(-s)=s^(-2x)[u^(2k)]H_B^2,
R_skip(-s)=-2s^(1-2x)[u^(2k-2)]H_C H_D.
```

At the retuned s/lambda the three coefficient extractions are identical;
the displayed prefactors give every exponent in (7).

The extra factor `t` in a crossing class is responsible for the exponent
loss. Here x=1, so

```text
Q_model,lambda(rho/lambda)
 =lambda[lambda W(rho)+R_skip(rho)].                    (8)
```

The dyadic scale `lambda_d=2^-25` gives a rational-coefficient hostile:
its B-root tuple is

```text
(2^-25,2^-19,2^-13,2^-7,2^-1),
```

and its original zero has `-rho/lambda_d in I`. Direct exact interval
evaluation proves `Q_model,lambda_d>0` and `Q_model,lambda_d-V_lambda_d>0`.
This refutes the full signed-transport implication for the combined model.

There is also an exact double-cancellation instance. Define the positive
algebraic number

```text
lambda_0=-R_skip(rho)/W(rho)>0.                         (9)
```

It belongs to the real algebraic field `Q(s_*)`, of degree at most four;
the output records the
remainders of `sW(-s)` and `sR_skip(-s)` in that field. Equations (7)--(9)
give exactly

```text
P_lambda0(rho/lambda0)=0,
Q_model,lambda0(rho/lambda0)=0.                         (10)
```

Every retained model predicate above still holds at this scale. No
approximate algebraic number is used to infer the cancellation.

## 4. What fails, and what actual normalization retains

The first failed implication is that the two Euler identities and
weighted compatibility determine the **relative normalization** of a
hit square and a crossing contribution. They do not. The ordinary
carriers at the chosen zero cannot see this free scale, while the raw
crossing factor `2t` can.

This is not a variable gauge applied to an actual doubled path row.
Under a genuine replacement `t->lambda t`, that factor would also
become `2lambda t`. In the abstract model it remains `2t`.

For the genuine parameters `(A,B,h,x,z)=(2,3,4,1,1)`, the translated
beta row is the fixed composition polynomial `F_15(t)`, whose linear
coefficient is13. Our normalization fixes only the constant coefficients
of G_B,G_C,G_D to one. Restoring the genuine first coefficients fixes

```text
lambda_norm=13/(1+64+64^2+64^3+64^4)=13/17043521,
([t]G_B,[t]G_C,[t]G_D)=(13,12,11).                      (11)
```

The positive K and skip signs survive every positive scale. The bounded
normalization control (11) has a negative full Q_model; it is not the
double-cancellation scale (9). Thus fixing this scalar does not repair
whole-class or skip negativity, and no conclusion about the full response
of every scalar-normalized model is asserted. The remaining exact
binomial coefficient law is stronger still.

The next concrete test is therefore the full response on the combined
model **with an explicit relative coefficient anchor**, or a faithful
path argument that retains the crossing-edge normalization throughout.
No new parameter-family census, universal normalized-model theorem, or
actual Laurent noncancellation claim is made here.

## 5. Exact verification and audit status

The standalone verifier uses Fraction arithmetic, literal coefficient
convolutions, exact root brackets, and coefficientwise scaling identities.
Its final declared universe is this geometric-root model, its dyadic and
linear-coefficient-normalized controls, one actual q2/h1 control, and the
earlier Euler hostile as a weighted-compatibility rejection control.
The exploratory numerical locator is not a proof input or a minimality
claim. No repository mathematical producer is imported.

The [standalone source](../../04-computation/combined_pencil_empty_core_morning_sep06.py)
passes **156 exact gates**, and its [frozen output](combined_pencil_empty_core_morning_sep06.out)
is byte-identical under normal Python and `python3 -O`. All sign tests use
rational Horner interval arithmetic before outward rounding for display.
Representative outward integer bounds, for the local positive phase times
the indicated raw response, are:

| Model | sW | sR_skip | sQ_model |
|---|---:|---:|---:|
| Base | [-26466334502,-26447251665] | [1553,1603] | [-26466332948,-26447250063] |
| Dyadic scale | [-789,-788] | [1553,1603] | [765,815] |
| Genuine linear anchor | [-20188,-20172] | [1553,1603] | [-18634,-18570] |

The program also verifies both Euler relations coefficientwise, all
raw scaling identities, literal ordinary-carrier constructions at three
rational phases, and the inherited zero-preserving lowering identity.
It checks the actual `(-9,1,6)` control at phase -2, where `P=0,Q=-6305`,
and excludes the previous Euler-only hostile through `C(1)D(1)<0`.

Reproduce from the repository root:

```bash
python3 04-computation/combined_pencil_empty_core_morning_sep06.py > /tmp/combined_pencil_replay.out
python3 -O 04-computation/combined_pencil_empty_core_morning_sep06.py > /tmp/combined_pencil_optimized.out
cmp /tmp/combined_pencil_replay.out /tmp/combined_pencil_optimized.out
cmp /tmp/combined_pencil_replay.out 05-knowledge/results/combined_pencil_empty_core_morning_sep06.out
```

SHA256:

- Source: `2f6ce6480a1642a357c48ffde2185a36bec3c7b7d7be195489747447c67da605`.
- Output: `e02df806c315993d45c5a4a5c0748a24e2ad8a5867f29978acd866b8ad198904`.
- Semantic certificate: `119b853c921308c8e1a9c523fa50661ad435ef9bf6df93974e76294ac72fa3f4`.

Root independently audited the all-weight proof, the exact Euler and
ordinary-carrier identities, the carry exponents, the algebraic cancellation,
and the distinction from actual path normalization: **PASS**.
`three_ray_geometry` independently confirmed the structural proof, all ten
strict-interlacer residues, the displayed first polynomial, and linear-anchor
ratios: **PASS**.
Root then independently read the complete source and reproduced both
156-gate runs with frozen output agreement: **PASS**, including rational
interval bounds, the smallest-root isolation, the lower carry, and the
literal ordinary-carrier path. All three files are frozen.
