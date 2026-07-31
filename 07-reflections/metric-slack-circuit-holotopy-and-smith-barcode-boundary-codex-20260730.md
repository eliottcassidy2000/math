# Metric slack, circuit holotopy, and the Smith-barcode boundary

Status: `PROVED` for the abstract closed-slack, mod-`p` cokernel, and weak
ray-compactification lemmas below; `VERIFIED-EXACT` for the displayed finite
controls; `SYNTHESIS ONLY` for transfers between Gaussian moments, LRC, and
carry ledgers. This reflection proves no new GMC, LRC, or minimal-constant
result.

## Trigger and inheritance

Three live objects looked superficially identical:

1. a deterministic Farkas separator for an LRC status row;
2. a wall-stripped Gaussian resultant and its calibrated Sylvester map; and
3. a carry-ledger coefficient untouched by legal cells.

The comparison became useful only after retaining their differences.
[THM-2065](../01-canon/theorems/THM-2065-two-anchor-fejer-circuit-ray-collapse.md)
shows that a relation circuit can persist without making any grade
infeasible. [THM-2658](../01-canon/theorems/THM-2658-balanced-lift-helly-circular-arc-gain-nerve-and-wrap-boundary.md)
and [THM-2672](../01-canon/theorems/THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap.md)
show that coarse label nerves can manufacture virtual spheres which disappear
after physical-component/gain refinement.
[THM-2174](../01-canon/theorems/THM-2174-endpoint-phase-scale-obstruction.md)
shows that exact current can vary like `C/W` along a residue ray even when
finite phase data are fixed. [THM-2974](../01-canon/theorems/THM-2974-discriminant-cover-integral-order-smith-and-owner-boundary.md)
shows that a common normalized algebra may preserve inertia while its two
integral graph orders are nonnested and its owner is lost.
[THM-2975](../01-canon/theorems/THM-2975-modular-six-sheet-schreier-graphs-and-farey-partial-cube-boundary.md)
adds the graph-theoretic version: a lossy Bass--Serre incidence shadow can be
a partial cube while the literal six-state graph contains triangles. A clean
quotient topology is not evidence that the physical state graph shares it.

These are not annoyances. They identify the coordinates that every lawful
bridge must carry.

## 1. The closed metric-slack theorem

Let `P` be a nonempty compact convex set, `|J|>=2`, and

```text
d_j(x)=h_j-L_j(x).
```

Assume the closed family `{d_j<=0 : j in J}` is minimally infeasible: its full
intersection is empty and every proper subfamily intersects. Define

```text
sigma_* = min_(x in P) max_(j in J) d_j(x).               (1)
```

Then `sigma_*>0`, and the nerve of the relaxed family `{d_j<=sigma}` is

```text
boundary Delta^(|J|-1),   0<=sigma<sigma_*,
Delta^(|J|-1),            sigma>=sigma_*.                 (2)
```

Thus its sole reduced-homology bar is `[0,sigma_*)` in degree `|J|-2`.
Convex minimax gives

```text
sigma_*
 = max_(lambda in Delta_J)
     min_(x in P) sum_j lambda_j d_j(x),                  (3)
```

and every optimizer has full support. Rational polytopal data give rational
`sigma_*`.

The nerve vertices are `supp(lambda)`, the grade constraints. They are not
the cell columns `supp(lambda A)`. A support-minimal covector of a column
configuration is a cocircuit; a positive minimal dependency in the
homogenized affine constraint configuration is a circuit. Relative to `P`,
the latter must include any facets of `P` used by the certificate.

### Metric, not oriented-matroid, persistence

The infeasible pair `x>=1`, `x<=0` has death `1/2` under equal slack. Scaling
the first deficit by `100` changes death to `100/101` without changing the
sign pattern or minimal-infeasibility type. The sphere type is combinatorial;
its lifetime is row-normalization dependent.

For strict-open inequalities, zero margin may occur. The top simplex then
appears only for `sigma>sigma_*`; an endpoint at
`sigma=0=sigma_*` is decorated/ephemeral. Formula `(2)` is the closed compact
theorem and must not be silently reused.

## 2. Mod-p cokernels are not automatically Bocksteins

Let `C` be free abelian and `B<=C` the subgroup of legal integral repairs.
Put

```text
Q_p(C,B)=C/(B+pC).                                        (4)
```

If a defect has nonzero class in `(4)`, no legal integral repair kills it.
For a directed system preserving the repair subgroups, `(4)` is functorial,
so a class surviving every transition is a persistent obstruction.

This is a mod-`p` cokernel statement. Calling it a Bockstein requires a
specified short exact sequence and connecting homomorphism.

Nor is there generally a canonical primitive Farkas ray. For
`A=I_2,y=(-1,-1)`, both coordinate rays are primitive extreme separators.
Reproducible evidence must bind the primal instance and the full exact
extreme-ray/sign-support set, or specify a deterministic canonicalization.
The solver basis is not semantic data. The deterministic z312-to-z298 replay
repair at commit `effdad261` is the current LRC model for this discipline.

## 3. Weak projective-ray compactification

Fix finitely many low cores and positive residue rays

```text
z_(rho,m)=r_rho+mL,       delta_(rho,m)=a_rho/z_(rho,m).
```

For a fixed low core `c`, suppose the weak scalar obligation is

```text
sum_(u in c) delta(u)+delta_(rho,m) >= T,
R(c)=T-sum_(u in c)delta(u).                              (5)
```

Adjoin `m=infinity` and set `delta_(rho,infinity)=0`. The fibre is exactly

```text
R(c)<=0: full compactified ray;
R(c)>0:  m<=floor((a_rho/R(c)-r_rho)/L),                  (6)
```

where the second line is a possibly empty finite initial segment. A
height-invariant residue certificate therefore kills an entire infinite
cylinder without a label horizon. This is the lawful topology behind
THM-2970's one-high residue rays.

The strict boundary is different. At `R(c)=0`, every finite point satisfies
the strict obligation but infinity does not, so the fibre is not compact in
the one-point topology. THM-2174 gives the complementary hostile: finite
phase invariance alone does not make exact current height-invariant.

## 4. The lawful Gaussian Smith barcode

At one fixed width `M`, let `S` be the integer Sylvester map of calibrated
polynomials `U_M,V_M`, with nonzero Smith factors `s_i`. Define

```text
Q_(p,k)=Z^D/(im(S)+p^k Z^D).                              (7)
```

Then

```text
log_p |Q_(p,k)| = sum_i min(k,v_p(s_i)).                  (8)
```

Equation `(8)` is a genuine `p`-adic barcode refining
`p | Res(U_M,V_M)`. It is not yet a GMC transition object:

- the positive wall-stripped core `N_M` is not proved to be a class in `(7)`;
- no canonical map from width `M` to width `M+1` has been constructed; and
- THM-2974 warns that shared normalized structure need not identify an
  integral order or owner.

The next Gaussian test is therefore not another analogy. It is to construct
an actual width-transition map and ask whether the pointed Smith class and
the positive core are compatible with it.

## 5. Source-target contract

```text
source:      exact repair lattice or closed feasibility family;
target:      mod-p cokernel, Farkas support, or selected nerve;
map:         quotient by repairs / nonnegative separation / nerve;
preserved:   impossibility of the specified repair;
destroyed:   respectively positivity, integral torsion, physical owner,
             row normalization, and exact current height;
needed:      transition maps plus the discarded sign/owner/metric/current;
hostiles:    nonunique extreme rays, row scaling, THM-2065 relation circuit,
             THM-2658/2672 virtual sphere, THM-2174 C/W ray,
             THM-2974 nonnested integral orders.
```

The word **holotopy** is useful here only for the selected proof-space nerve:
it remembers how a minimal infeasible family becomes feasible. It is not a
claim about physical LRC Cech topology.

## 6. Exact controls and next tests

Run

```text
python 04-computation/circuit_holotopy_slack_persistence_probe.py
python -O 04-computation/circuit_holotopy_slack_persistence_probe.py
```

Both modes byte-match

```text
05-knowledge/results/circuit_holotopy_slack_persistence_probe.out.
```

LF-normalized SHA-256:

```text
script  0d750b0658f8dc7b3037905e44be32ea87e884057c084279712da8a5d9c040ae
output  bf6544f2dc71373e8ede515fb3d462f9e313a6000a088718a185cd60308e044f
```

The controls verify deaths `1/2` and `1/3`, the `100/101` row-scaling
hostile, a nonminimal contractible nerve, the exact finite ray
`m=0,...,4`, residue-phase invariance, and the strict zero-residual infinity
failure.

Cheapest honest successors:

1. for LRC, use full circular phase diameter on the projective ray class,
   then retain the physical packet labels that the residue quotient loses;
2. for Gaussian moments, build one explicit Smith-compatible width map before
   interpreting `(8)` persistently;
3. for carry ledgers, prove the quotient class has no later leakage rather
   than extrapolating from finitely many frozen deadlines.

That is the useful stopping boundary: topology organizes an already lawful
obstruction, but it cannot manufacture the missing transition map or owner.
