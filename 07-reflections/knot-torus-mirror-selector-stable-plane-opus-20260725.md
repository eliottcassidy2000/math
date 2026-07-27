---
source: codex-2026-07-25-knot-torus-mirror-plane
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For every K_g=T(2,2g+1), g>=2, two switched cross-frequency inertia
  gauges give the exact real-Blanchfield/common-context floor
  P+Q+(g-1)|P-Q|.  The complete bank of all two-frequency inertia
  gauges is dominated by the switched pair and one odd-axis gauge, so
  it certifies exactly 2 per balanced mirror pair and no more.  On the
  abstract stable mirror plane this sharpens the defect interval to
  0<=delta<=2(g-1), gives two exact curvature budgets, forces rigidity
  at both endpoints, and proves genuine continuum ambiguity at every
  interior delta.  No ordinary unknotting equality or endpoint
  realizability is asserted.
depends_on:
  - THM-2191-catalytic-localization-of-the-gordian-metric
related:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2308-mirror-double-nakanishi-floor-and-sharp-stable-mixture-profile
  - THM-2346-global-allocation-anova-normal-form
  - THM-2348-prime-type-rectangularity-and-target-token-conditioning
  - knot-blanchfield-mirror-mixture-selector-opus-20260725
  - knot-torus-mirror-curvature-extreme-atlas-opus-20260725
  - knot-curvature-mixed-cohabitation-kernel-opus-20260725
external:
  - "Maciej Borodzik and Stefan Friedl, The unknotting number and classical invariants I, arXiv:1203.3225."
  - "Maciej Borodzik and Stefan Friedl, The unknotting number and classical invariants II, arXiv:1207.2413."
  - "Charles Livingston, Signature invariants related to the unknotting number, arXiv:1710.10530."
script: 04-computation/knot_torus_mirror_selector_stable_plane_probe.py
output: 05-knowledge/results/knot_torus_mirror_selector_stable_plane_probe.out
script_sha256: 788a39c429a3dd357fc9d68fb281bbe0fa13894202e7078f6a3cf926b106a8f5
output_sha256: 4f134d5e8460178415f4fc49f8bf30c4a9b370e214df18fa629ca7987d9dea4c
hash_basis: LF-normalized working-tree bytes
---

# The torus--mirror family has two rigid stable endpoints

## 1. Inheritance and the new question

The earlier `T(2,7)` calculation found one useful switched pair of
Levine--Tristram inertia gauges.  The sharper question is not merely
whether that pair generalizes, but whether another choice of two
frequencies can do better on the balanced mirror ray.

For

```text
K_g=T(2,2g+1),                 g>=2,

a=[K_g],                       b=[mirror(K_g)],

e=a+b,                         o=a-b,               (1)
```

the answer is exact:

```text
the switched pair and the odd axis dominate the entire
two-frequency inertia bank.                          (2)
```

This both proves the family lower bound and exposes its stopping
boundary.  The faithful object is a centrally symmetric dual polygon,
not a tournament.

## 2. Signature and nullity cells

Choose the chirality for which the simple Alexander roots of `K_g` have
Levine--Tristram signature and nullity

```text
sigma=1,3,...,2g-1,...,3,1,

eta=1.                                               (3)
```

On the complementary arcs the signatures are the even values

```text
0,2,...,2g,...,2,0,

eta=0.                                               (4)
```

The mirror reverses `sigma` and preserves `eta`.  These are the standard
`A_(2g)` inertia cells; the tridiagonal eigenvalue calculation is the
same one recorded for `g=3` in
`07-reflections/knot-blanchfield-mirror-mixture-selector-opus-20260725.md`.

For two frequencies `alpha,beta in S^1\{1}`, with the
Borodzik--Friedl convention `eta(1)=0` retained at the omitted endpoint,
define the additive gauge

```text
Phi_(alpha,beta)(J)
 =(1/2)[eta_J(alpha)+sigma_J(alpha)
       +eta_J(beta)-sigma_J(beta)].                 (5)
```

The common rank-one semidefinite crossing-change model shows

```text
|Phi_(alpha,beta)(J)-Phi_(alpha,beta)(J')|
 <=d_G(J,J').                                       (6)
```

Thus every such gauge survives addition of an arbitrary common context.

## 3. The switched determinant and exact mixture floor

Choose one root with signature `2g-1` and one with signature `1`.
Equation (5), together with its mirror switch, gives the coefficient
matrix on `(a,b)`

```text
          [ g    2-g ]
          [2-g    g  ],                              (7)
```

whose determinant is

```text
4(g-1).                                             (8)
```

For `P,Q>=0`,

```text
max(|gP+(2-g)Q|, |(2-g)P+gQ|)
 =P+Q+(g-1)|P-Q|.                                  (9)
```

Because the gauges are additive, for every common context `C`,

```text
d_G((#^P K_g)#(#^Q mirror(K_g))#C,C)
 >=P+Q+(g-1)|P-Q|.                                 (10)
```

The same expression is the exact real Blanchfield matrix size.  Indeed,
interchange `K_g` and its mirror if necessary so that `P>=Q`, and write

```text
N=P+Q,                         d=|P-Q|.
```

At the extreme and smallest roots, respectively,

```text
max(eta+sigma)=N+(2g-1)d,

max(eta-sigma)=N-d.                                 (11)
```

The largest off-root signature `2gd` does not exceed the first term
because `N>=d`.  Therefore the two maxima give

```text
mu_BF=(1/2)[N+(2g-1)d+N-d]
     =N+(g-1)d.                                    (11a)
```

All roots of
`Delta_(T(2,2g+1))` are simple and lie on `S^1`, so the remaining
Borodzik--Friedl input is

```text
eta_BF=max_(C^*) eta=N.                             (11b)
```

Thus `mu_BF>=eta_BF`.  The
Borodzik--Friedl formula therefore gives

```text
n_R((#^P K_g)#(#^Q mirror(K_g)))
 =N+(g-1)d.                                        (12)
```

Equations (10)--(12) are lower certificates.  They do not compute the
ordinary unknotting number.

## 4. Classification of every two-frequency inertia row

For an arbitrary pair `alpha,beta`, put

```text
epsilon=eta(alpha)+eta(beta),

d=sigma(alpha)-sigma(beta).                         (13)
```

On `se+to`, the gauge (5) is the row

```text
epsilon*s+d*t.                                     (14)
```

The parity cells (3)--(4) give the sharp possibilities

```text
epsilon in {0,1,2},

|d|<=2g-epsilon,

d=epsilon (mod 2).                                 (15)
```

Conversely, every parity-compatible integer in the displayed interval
occurs, so (15) is the exact row classification rather than only its
convex envelope.

Let

```text
L=2(g-1),                    X=|s|,       Y=|t|.
```

Every row in (14) is bounded by

```text
epsilon X+|d|Y
 <=max(2X+LY, 2gY).                                (16)
```

For `epsilon=0` and `2` this is immediate.  For `epsilon=1`,

```text
X+(2g-1)Y
 <=2X+(2g-2)Y              if X>=Y,

X+(2g-1)Y
 <=2gY                     if X<=Y.                (17)
```

Both extremal supports in (16) actually occur:

```text
the switched root pair gives 2X+LY;

a maximal-signature arc paired with a zero-signature arc gives 2gY.
                                                               (18)
```

Consequently the dual seminorm generated by **all** gauges (5) is
exactly

```text
q_2freq(s,t)=max(2|s|+L|t|, 2g|t|).                (19)
```

In particular,

```text
q_2freq(e)=2.                                      (20)
```

No selection, supremum, or convex combination of these two-frequency
inertia rows can certify more than two units per balanced mirror pair.
This is the precise classical-selector stopping boundary; it is not a
statement that the stable Gordian norm itself equals two there.

## 5. The abstract stable plane

Let `p=u_hash` be the stable Gordian norm supplied by THM-2191.  The
axis and odd-ray calibrations are

```text
p(a)=p(b)=g,

p(o)=2g.                                           (21)
```

The lower equality follows from the torus-knot signature and the
standard `g`-crossing upper route.  The second follows from the odd-axis
row in (18) and the triangle inequality.

Write

```text
p(e)=2g-delta.                                     (22)
```

The triangle inequality and (20) imply

```text
0<=delta<=L.                                       (23)
```

For positive mixtures define

```text
D(P,Q)=g(P+Q)-p(Pa+Qb),

D(P,Q)=P h(Q/P),                         P>=Q>=0.   (24)
```

Convexity of `p`, homogeneity, and mirror symmetry make `D` homogeneous,
symmetric, and concave.  Thus

```text
h(z)>=delta z.                                     (25)
```

Symmetry gives `h'_-(1)>=delta/2`; the supporting line at one gives

```text
h(z)<=delta(1+z)/2.                                (26)
```

Finally (9) gives

```text
h(z)<=Lz.                                          (27)
```

Hence the exact inherited defect region is

```text
delta z
 <=h(z)
 <=min(Lz,delta(1+z)/2).                           (28)
```

## 6. Sharp norm envelope

For `X=|s|` and `Y=|t|`, equations (25)--(28) give

```text
max((2g-delta)X, 2X+LY, 2gY)
 <=p(se+to)
 <=(2g-delta)max(X,Y)+delta Y.                     (29)
```

When `X>=Y`, this follows by writing the positive coefficients of
`a,b` as `X+Y` and `X-Y`.  When `Y>=X`, equality in the odd-axis triangle
forces

```text
p(se+to)=2gY.                                      (30)
```

The two sides of (29) are themselves unconditional norms with all data
(21)--(23):

```text
p_low=max((2g-delta)X,2X+LY,2gY),

p_high=(2g-delta)max(X,Y)+delta Y.                  (31)
```

They correspond respectively to

```text
h_high(z)=min(Lz,delta(1+z)/2),

h_low(z)=delta z.                                  (32)
```

For every interior `0<delta<L` and every `X>Y>0`, both comparisons with
the maximum in `p_low` are strict, so

```text
p_low(s,t)<p_high(s,t).                            (33)
```

Thus the inherited data do not determine an interior stable plane.

## 7. Curvature measure and two endpoint rigidities

Let

```text
kappa=-D^2h
```

be the canonical nonnegative curvature measure on `(0,1)`, and put

```text
A=integral_(0,1)(1-t)dkappa(t),

B=integral_(0,1)t dkappa(t).                       (34)
```

The Green-kernel representation of (28) is

```text
h(z)=delta z+integral G(z,t)dkappa(t),

G(z,t)=min(z,t)(1-max(z,t)).                        (35)
```

The two upper supports in (28) are exactly

```text
A<=L-delta,

B<=delta/2.                                        (36)
```

At `delta=0`, the second budget forces `kappa=0`; at `delta=L`, the first
budget forces `kappa=0`.  Therefore both endpoints are rigid:

```text
delta=0:
  h=0,
  p(se+to)=2g max(X,Y);

delta=L:
  h=Lz,
  p(se+to)=
    2X+LY  if X>=Y,
    2gY    if Y>=X.                                (37)
```

For `0<delta<L`, the upper-defect extremal is the single atom

```text
t_0=delta/(2L-delta),

m_0=(2L-delta)/2,

kappa=m_0 delta_(t_0).                             (38)
```

It saturates both budgets:

```text
m_0(1-t_0)=L-delta,

m_0t_0=delta/2.                                    (39)
```

The zero measure and (38), together with their convex interpolants,
already give continuum many distinct abstract profiles at each interior
`delta`.

## 8. Boundary and interpretation

For `g=3`, the new formulas recover the `T(2,7)` bounds

```text
L=4,

delta z<=h(z)<=min(4z,delta(1+z)/2).                (40)
```

The explicit Brittenham--Hermiller construction separately narrows the
actual `T(2,7)` range to `1<=delta<=4`; that family-specific upper route
is not assumed here.

At `g=1` the determinant (8) vanishes and `L=0`.  The switched gauges
collapse to one line.  This is an exact algebraic boundary, but it does
not by itself explain any trefoil symbiont phenomenon.

The theorem refines connected-sum unknotting behavior in an
operation-ready way:

```text
source:
  two-frequency signature/nullity inertia;

target:
  a dual polygon and a curvature measure on the stable mirror plane;

preserved:
  additive common-context lower bounds and every positive mixture ray;

destroyed:
  finite unknotting paths, catalysts, and actual realization of an
  abstract interior norm;

needed sidecar:
  a non-inertia invariant or explicit stable upper construction to
  decide the balanced value delta.                  (41)
```

## 9. Exact audit

Run

```text
C:\Users\Eliott\.cache\codex-runtimes\codex-primary-runtime\
dependencies\python\python.exe
  04-computation/knot_torus_mirror_selector_stable_plane_probe.py
```

and repeat with `-O`.  The script checks `g=2,...,30`, the determinant,
the exact mixture floor, every row allowed by (15), both norm envelopes,
both endpoint rigidities, and the atomic moment identities using exact
rational arithmetic.  These loops are hostile controls for the symbolic
all-`g` proof above, not a finite substitute for it.

Stored transcript:

```text
05-knowledge/results/knot_torus_mirror_selector_stable_plane_probe.out
```

Normal, optimized, and stored transcripts match.  Raw-byte SHA-256:

```text
script  788a39c429a3dd357fc9d68fb281bbe0fa13894202e7078f6a3cf926b106a8f5
output  4f134d5e8460178415f4fc49f8bf30c4a9b370e214df18fa629ca7987d9dea4c
```
