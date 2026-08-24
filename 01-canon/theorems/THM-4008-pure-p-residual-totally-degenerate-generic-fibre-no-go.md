---
id: THM-4008
title: "Pure-p residuals have totally degenerate source fibres and cannot occur in the reduced 2:3 cell"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY GEOMETRIC-AUDITED. In THM-3992's
  normalized reduced 2:3 cell, assume the
  entire residual is R=f(p), not merely that its first displayed terms are
  pure in p. Put H=lambda*p+f(p) and m=deg H. The source generic fibre is the
  smooth genus-m hyperelliptic curve V^2=p(q-H)(q+gamma-H). At q=infinity it
  has a semistable model whose special fibre has rational normalization and
  exactly m nodes; after regular resolution every special component is
  rational. The target generic fibre has good elliptic reduction on the same
  base change. Resolving the hypothetical generic fibre map would therefore
  give a map from a connected union of rational special components to a
  smooth elliptic curve, hence a constant special map, contradicting
  conservation of the positive generic pullback degree. Equivalently, the
  source Jacobian is purely toric of rank m and cannot have the good elliptic
  target as a quotient. Thus every pure-p residual is excluded. An exact
  maximum-weight-six p^3+y^2 face is a sharp hostile: it restores a j=0
  elliptic component matching the target good fibre. It is not the leading
  model when a higher-weight coefficient is nonzero. THM-4017 refutes the
  former application of this face to the sharp 5x5 data because p^4 is
  forced there. The full reduced 2:3 cell and JC(2) remain open.
source: root + generic_fibre_residual / planar Jacobian continuation, 2026-08-24
audit: >
  PASS (primary + independent SymPy-free geometric audit, 2026-08-24). The
  primary certificate checks the function-field reconstruction, all infinity
  exponents and node thicknesses through degree eight, the target good model,
  the live cubic, and the mixed y^2 failure boundary. The independent audit
  reconstructs the field identities on exact hostile samples and checks the
  branch/genus, resolution-graph, source/target scaling, degree-conservation,
  mixed elliptic, and Eisenstein-unit ledgers through degree 24. Normal and
  optimized streams match both frozen outputs. The pure-p no-go in Sections
  1--4 is unconditional. Section 6 is an exact maximum-weight-six boundary
  model; Section 6.1's torsion invoice requires explicit degree and clutch
  owners. Its former join to THM-4016's fixed sharp residual is refuted by
  THM-4017 before the stable-map gate.
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
related:
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-3999-companion-divisor-boundary-endpoint-and-class-ledger
  - THM-4007-live-two-three-third-normal-row-five-weight-floor
  - THM-4012-weighted-leading-face-good-elliptic-factor-observer
  - THM-4016-sharp-five-by-five-elliptic-attachment-nontorsion
  - THM-4017-sharp-weight-eight-specialization-obstruction-and-newton-ledger
script: 04-computation/jc2_pure_p_residual_reduction_thm4008.py
output: 05-knowledge/results/jc2_pure_p_residual_reduction_thm4008.out
script_sha256: 42dfd32bb5f0ae027d96fe738f5db4734878bdc46c739c028104077e4a5ece0e
output_sha256: 24accea9e998c0ec953ddc8c6f1aa213173ac863088a043141e19aa201ee5788
independent_audit_script: 04-computation/jc2_pure_p_residual_reduction_thm4008_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_pure_p_residual_reduction_thm4008_independent_audit.out
independent_audit_script_sha256: 3a1e0144770f77f9e0a2ca52cdce79816c646237a8c0652b5332f0cf38a7f63d
independent_audit_output_sha256: 508c278b4490fba1600f05ff9ed30ec9831a9934dca7cf6cb76b4fdeeb0ed608
independent_audit_semantic_sha256: 60e835b4dcf778e3734508cc1488812cc04c90a214236150c3a874818d78bb08
hash_basis: raw LF bytes
---

# THM-4008 -- pure-`p` residuals are totally degenerate

**PROVED + VERIFIED-EXACT + INDEPENDENTLY GEOMETRIC-AUDITED.**

Work over an algebraically closed field `k` of characteristic zero. Use the
normalized reduced `(2,3)` cell of THM-3992:

```text
s=xt,                 p=s^2+t,              y=sp,
u=s^2/t,

P(A,C)=C^2-A^3+(3a^2/4)A+a^3/4,
G=P(A,C)=gamma*u+lambda*p+R(p,y),

lambda=3a/(2gamma),   a*gamma!=0,            R in (p^2,y).       (1)
```

Assume throughout the theorem's closed lane that the **entire** residual is

```text
R(p,y)=f(p),                    f in p^2 k[p].             (2)
```

Put

```text
H(p)=lambda*p+f(p),             m=deg H.                   (3)
```

Then `m>=1`; `m=1` is exactly `f=0`, while every nonzero `f` gives `m>=2`.
The conclusion is

```text
no normalized reduced (2,3) Keller pair satisfies (2).     (4)
```

This is an all-coefficient statement about `(2)`. It does not follow from a
finite pure-`p` truncation, and it does not exclude a residual with any mixed
`y`-term.

## 1. Exact source and target generic fibres

Use the target value itself as parameter:

```text
q=G(A,C),                         K=k(q).                  (5)
```

The algebraic independence of `A,C` makes `q` transcendental. Since
`t=p-s^2`, equation `(1)` becomes

```text
s^2(q+gamma-H(p))=p(q-H(p)).                              (6)
```

With

```text
V=s(q+gamma-H(p)),                                       (7)
```

this gives the hyperelliptic model

```text
C_H: V^2=p(q-H(p))(q+gamma-H(p)).                         (8)
```

Conversely,

```text
s=V/(q+gamma-H(p)),       t=p-s^2,       x=s/t,           (9)
```

so `(8)` has exactly the source generic-fibre function field:

```text
K(C_H)=K(p,V)=k(x,t).                                    (10)
```

The polynomial on the right of `(8)` is squarefree over `K`. Indeed,
`p` is coprime to the other two factors because `H(0)=0` and `q` is
transcendental; the latter two factors differ by the nonzero constant
`gamma`; and a common zero of `H(p)-q` and `H'(p)` would make the
transcendental `q` equal to a critical value of a polynomial over `k`.
The same argument applies after replacing `q` by `q+gamma`. The degree is
`2m+1`, hence the smooth projective model of `(8)` has

```text
g(C_H)=m.                                                 (11)
```

The target generic fibre is the smooth elliptic curve

```text
E_T: C^2=A^3-(3a^2/4)A+q-a^3/4.                          (12)
```

The Keller map makes `k(x,t)/k(A,C)` finite. Restricting over `K=k(q)`
therefore gives a finite nonconstant morphism of smooth projective curves

```text
phi:C_H -> E_T.                                          (13)
```

This is the only use of the Keller condition in the generic-fibre argument:
it supplies the finite function-field inclusion and hence `(13)`.

## 2. The source has a rational `m`-node special fibre

Write

```text
H(p)=c_m p^m+...+c_1 p,                 c_m!=0.           (14)
```

At `q=infinity`, make the common finite base change

```text
q=rho^(-6m),               p=rho^(-6)z,
W=rho^(6m+3)V,             O=k[[rho]].                   (15)
```

Define the integral polynomial

```text
H_rho(z)=rho^(6m)H(rho^(-6)z)
        =c_m z^m+c_(m-1)rho^6 z^(m-1)+...+c_1 rho^(6m-6)z.
                                                                  (16)
```

Equation `(8)` becomes

```text
W^2=z(1-H_rho(z))(1+gamma*rho^(6m)-H_rho(z)).            (17)
```

The right side has constant odd degree `2m+1` and unit leading coefficient
`c_m^2`. Its standard hyperelliptic projective closure is flat over `O` and
has one smooth point at infinity. The special fibre is

```text
C_0: W^2=z(1-c_m z^m)^2.                                 (18)
```

It is irreducible, and its normalization is explicitly

```text
P1_r -> C_0,
z=r^2,                    W=r(1-c_m r^(2m)).             (19)
```

The `m` roots of `1-c_m z^m` are simple and nonzero. Each is an ordinary
node of `(18)`, and there are no other singularities. Thus `C_0` has rational
normalization and exactly `m` nodes, accounting for its arithmetic genus
`m`. For `m>=2` it is stable. For `m=1` it is the usual semistable rational
nodal cubic; stability is not needed below.

The thickness at every node is also exact. Near a root of `1-c_m z^m`, put

```text
B=1-H_rho(z),                 N=6m,
delta=gamma*rho^N/2.                                       (20)
```

Here `B` is a regular parameter, `z` is a unit, and the completed local ring
contains a unit `r_z` with `r_z^2=z`. With

```text
U=W-r_z(B+delta),             Z=W+r_z(B+delta),           (21)
```

equation `(17)` gives

```text
UZ=-z delta^2=unit*rho^(12m).                             (22)
```

Resolving this `A_(12m-1)` thickness singularity inserts only a chain of
rational curves. Doing so at every node produces a proper regular model
`mathcal C/O` whose strict special component is the `P1` in `(19)` and whose
other special components are rational chains. The special fibre is reduced
and connected. A further finite base change only increases the thickness in
`(22)` and again resolves by rational chains; it never creates a
positive-genus component.

## 3. The target has good elliptic reduction on the same base change

Under `(15)`, set

```text
A=rho^(-2m)X,                  C=rho^(-3m)Y.              (23)
```

Multiplying `(12)` by `rho^(6m)` gives

```text
Y^2=X^3+1-(3a^2/4)rho^(4m)X-(a^3/4)rho^(6m).             (24)
```

Its special fibre is

```text
E_0: Y^2=X^3+1,                                           (25)
```

which is smooth in characteristic zero. Hence `(24)` is a smooth proper
elliptic scheme `mathcal E/O`. This is actual good reduction, not merely a
nonnegative `j`-valuation assertion.

For `m=1`, `(18)` is multiplicative while `(25)` is good. This recovers the
reduction mismatch used independently for `R=0` in THM-3997. The present
argument keeps all `m` and does not require source genus one.

## 4. Direct curve-model contradiction

Base-change the hypothetical morphism `(13)` to `Frac(O)`. It is a rational
map

```text
mathcal C  -->>  mathcal E                                (26)
```

defined on the generic fibre. Since `mathcal C` is a regular surface and
`mathcal E` is projective, elimination of indeterminacy for surface maps
resolves `(26)` by a finite sequence of point blowups supported in the
special fibre:

```text
pi:mathcal C' -> mathcal C,       Phi:mathcal C' -> mathcal E.    (27)
```

Every new exceptional component is a `P1`; thus every irreducible component
of the special fibre of `mathcal C'` is rational. Any morphism from a smooth
projective rational curve to the elliptic curve `E_0` is constant, by
Riemann--Hurwitz. The constants agree at component intersections. Since the
special fibre is connected, `Phi` is constant on its whole reduced support.

Here fibre multiplicities cause no loss. The regular model obtained from
`(22)` has reduced special fibre. A later blowup at an intersection can give
an exceptional component the sum of the adjacent multiplicities, but every
component still maps constantly. For any relatively ample line bundle `L`
on `mathcal E`, flatness and constancy of relative degree (equivalently, the
intersection formula with the full special fibre) give

```text
deg(Phi_eta^* L)
 =sum_i e_i deg((Phi|C_i)^*(L|E_0))
 =0,                                                       (28)
```

where `sum_i e_i C_i` is the full special fibre, with its actual positive
multiplicities. On the other hand, `(13)` is finite and nonconstant, so

```text
deg(Phi_eta^* L)=deg(phi)*deg(L|E_T)>0.                   (29)
```

Equations `(28)--(29)` contradict each other. This proves `(4)`.

### 4.1 Equivalent toric-Jacobian sidecar

The same obstruction has an exact group-theoretic reading. The dual graph of
`C_0` has one rational vertex and `m` loop edges. The normalization sequence
for units gives

```text
1 -> G_m^m -> Pic^0(C_0) -> Pic^0(P1) -> 0,
Pic^0(C_0)=G_m^m.                                        (30)
```

Subdividing the loops by the rational chains from `(22)` leaves their first
Betti number `m`. Thus `Jac(C_H)` has purely toric semistable reduction of
toric rank `m=dim Jac(C_H)`. Any abelian quotient of a purely toric
semistable abelian variety is purely toric: by Poincare reducibility it has a
complement up to isogeny, toric rank is additive under products and invariant
under semistable isogeny, and equality of total toric rank with total
dimension forces equality on each factor. The elliptic quotient induced by
`(13)` would therefore be toric, contradicting the good reduction in `(25)`.

This sidecar is equivalent to the direct degree argument, not an additional
assumption in it.

## 5. The live forced cubic and the finite-truncation firewall

On THM-3997's live seam with `b=beta/gamma=0`, one has

```text
gamma=-a^3/2,                 lambda=-3/a^2,
[p^2]R=8/(3a^7),              [p^3]R=-1376/(135a^12).    (31)
```

If one additionally assumes that these are the **entire** residual, put

```text
z=p/a^5,                      Q=q/gamma,
h(z)=H(a^5z)/gamma
    =6z-(16/3)z^2+(2752/135)z^3.                         (32)
```

The source fibre is the genus-three curve

```text
v^2=z(Q-h(z))(Q+1-h(z)).                                 (33)
```

Its discriminant support in the `Q`-line is, up to a nonzero constant,

```text
Q^2(Q+1)^2
 (7396Q^2-7340Q+10935)
 (7396Q^2+7452Q+10991),                                  (34)
```

and its infinity fibre is the rational three-node curve

```text
v_0^2=Z(1-(2752/135)Z^3)^2.                              (35)
```

Thus the exact pure cubic truncation is excluded by the theorem. But
`b=d=0` does **not** say that `(31)` is the entire residual: terms such as
`y^2`, `p^2y`, and higher rows remain available. No use of `(31)--(35)` may
erase them.

## 6. Exact max-weight-six hostile: `y^2` restores the missing good factor

The first failure boundary is structural. Give `s,p,y` weights `1,2,3` and
suppose the **entire polynomial** `H` has maximum weighted degree six, with
top face

```text
H_6(p,y)=epsilon*p^3+kappa*y^2,
epsilon*kappa*(epsilon+kappa)!=0.                         (36)
```

The live coefficient in `(31)` makes `epsilon!=0`; `kappa=[y^2]R` is not
fixed by THM-3997. The absence of all terms of weight greater than six is
load-bearing. Under

```text
q=rho^-6,                 s=rho^-1 S,
p=rho^-2 P,               y=rho^-3 SP,                   (37)
```

the leading source equation factors as

```text
(S^2-P)(1-epsilon P^3-kappa S^2P^2)=0.                   (38)
```

The first component is rational. On the second put `T=SP`; it becomes

```text
E_kappa: kappa*T^2=1-epsilon*P^3.                         (39)
```

This is a smooth elliptic curve with `j=0`. Because `k` is algebraically
closed, choosing `u^3=-epsilon` and `v^2=kappa` sends `(39)` to
`Y^2=X^3+1`, exactly the target special fibre `(25)`.

The two components meet where

```text
P=S^2,       T=S^3,       (epsilon+kappa)S^6=1.           (40)
```

These are six transverse points. The dual graph has two vertices and six
edges, so its first Betti number is five; the reduction has toric rank five
and abelian part `(39)` of dimension one. The corresponding Newton polygon

```text
(0,1),(2,0),(4,2),(2,3),(0,4)
```

has six interior lattice points, matching total genus `1+5=6` in the
nondegenerate face.

A displayed lower row is not enough for this scaling. A monomial `p^i y^j`
of weight `2i+3j` enters the rescaled equation with rho-exponent
`6-(2i+3j)`. Any fixed nonzero higher-weight coefficient makes the proposed
family nonintegral unless a separate coefficient family gives it
compensating positive valuation. This is precisely the failure exposed by
THM-4017 at the sharp `5x5` point.

Therefore the proof of Sections 2--4 stops sharply at `kappa!=0`: the first
mixed term restores precisely the potential-good elliptic factor needed to
carry a map to the target. Reduction type alone cannot exclude this face.

### 6.1 Conditional torsion invoice on the six attachments

**CONDITIONAL.** There is nevertheless a precise next condition. Assume a
proper flat resolved model through `(36)--(40)` in which `E_kappa` owns
positive degree (for example, it is the unique positive-genus component) and
all six branches meet one connected contracted rational clutch. Every
rational component maps constantly; the degree owner makes the restriction
to `E_kappa` nonconstant, while the clutch owner gives all six points one
common image. After translation this restriction is an isogeny.

Let `zeta_6` be a primitive sixth root and `zeta_3=zeta_6^2`. The six points
in `(40)` form the orbit of

```text
sigma(P,T)=(zeta_3 P,-T),                                (41)
```

an order-six elliptic automorphism fixing the origin. In the Eisenstein
endomorphism ring, `sigma-1` is a unit. If `P_0` is one attachment point and
`psi:E_kappa->E_0` is the specialized isogeny, equality of the first two node
images gives

```text
psi((sigma-1)P_0)=0.                                     (42)
```

The kernel is finite and `sigma-1` is invertible, so

```text
P_0 is a torsion point of E_kappa.                        (43)
```

THM-4016 performs the all-order arithmetic test for the formal weight-six
truncation of THM-4007's sharp `5x5` coefficients: its six normalized points
are non-torsion. THM-4017 proves that the same sharp calculation forces a
nonzero `p^4` coefficient, which has rho-exponent `-2` here. Therefore those
points are not attachment points of this weight-six model for the fixed sharp
residual. The former direct join is **REFUTED**.

For comparison, THM-4012's actual exact maximum-weight-six branch sets
`[p^4]R=0` and has different attachment ratios
`(43/224,267/224)`. Its face-stability and exclusion are proved separately.

The resonance `epsilon+kappa=0`, higher weighted faces, and multiple leading
monomials require separate models.

## 7. Exact certificate and honest scope

The companion script verifies:

1. the function-field identity `(6)--(10)`;
2. degrees, squarefree controls, infinity scaling, normalization, node count,
   and local thickness through `m=8`;
3. the target good model on every corresponding base change;
4. the live cubic coefficients and discriminant support `(31)--(35)`; and
5. the mixed factorization, elliptic `j=0` model, transversality, and Newton
   genus ledger `(36)--(40)`.

Reproduction:

```bash
python3 04-computation/jc2_pure_p_residual_reduction_thm4008.py
python3 -O 04-computation/jc2_pure_p_residual_reduction_thm4008.py
python3 -B 04-computation/jc2_pure_p_residual_reduction_thm4008_independent_audit.py
python3 -B -O 04-computation/jc2_pure_p_residual_reduction_thm4008_independent_audit.py
```

Both outputs byte-match
`05-knowledge/results/jc2_pure_p_residual_reduction_thm4008.out` and end with
`ALL THM-4008 EXACT CHECKS PASSED`.

**PROVED:** every normalized reduced `(2,3)` pair whose entire residual lies
in `k[p]` is impossible. The independent streams match their frozen
`36`-gate transcript and have semantic SHA-256
`60e835b4dcf778e3734508cc1488812cc04c90a214236150c3a874818d78bb08`.

**Not proved:**

1. that a finite truncation with only displayed pure-`p` coefficients has no
   later mixed terms;
2. that `[y^2]R` vanishes, is nonzero, or is forced to a particular value;
3. that the mixed torsion invoice `(43)` applies to a residual with any
   fixed nonzero higher-weight coefficient;
4. exclusion of multiple-term higher weighted faces or the full reduced
   `(2,3)` cell;
5. `JC(2)`.
