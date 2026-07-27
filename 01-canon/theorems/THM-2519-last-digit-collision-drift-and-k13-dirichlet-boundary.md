---
id: THM-2519
title: "Last-digit collision drift and the K13 Dirichlet boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Average the two-arm inverse-branch
  cospans with a fixed last base-13 difference u.  Their future-owner-weighted
  same-last-digit excess is exactly a conditional variance, equivalently the
  complete-graph K_13 Dirichlet energy of the thirteen predecessor-root
  masses.  Its twelve nontrivial u-Fourier coefficients are weighted squared
  norms of the corresponding high-frequency ladders.  For rational step data
  they vanish all together or are all strictly positive.  The construction
  is antipodally even and retains old deep data only on its ancestry sheets.
  A Boolean future-measurable delta-plus-replicas family has complete bulk
  square spectrum but zero last-digit drift, so bulk moment charge alone does
  not force a future collision colour.
source: codex-2026-07-27-last-digit-collision-drift
depends_on: []
related:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2513-anchored-first-or-second-moment-spectrum-and-pair-space-boundary
  - THM-2514-cyclic-k14-factor-chart-and-six-phase-ordinary-degree-reconstruction
  - THM-2518-perron-inverse-branch-owner-word-cospan-recovery
script: 04-computation/lrc14_perron_collision_cospan_thm2518_2519_referee.py
output: 05-knowledge/results/lrc14_perron_collision_cospan_thm2518_2519_referee.out
script_sha256: b09f9e8940914729583057cf5c92833a0f8dd85673b47bf1b2d5d71afe895e39
output_sha256: 2818c04f4ef56f352e0365504d7b2f8bf194bf8bdb0b04e4268a5340acba7898
hash_basis: working-tree bytes (LF)
---

# THM-2519 -- last-digit collision drift is conditional variance

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The depth-`L` inverse tree of `x -> 13x` has two distinct coordinates near
its root:

```text
last digit u in F_13,             higher address e mod 13^(L-1).       (1)
```

The last digit is the predecessor-root colour.  The higher address is the
transport stalk.  Keeping the first and averaging the second produces an
exact positive observable: the conditional variance of the thirteen
predecessor masses over a marked future owner.

This identifies a precise functional form for one missing drift.  It also
identifies its sharp nullspace.  A nonzero bulk first or square moment need
not retain any last digit; a response that is already measurable at the
future base is an exact hostile.

## 1. Perron and inverse-branch notation

Work on `T=R/Z` with normalized Haar measure.  For an integer `m>=1`, use
the normalized Perron operator

```text
(P_m F)(y)=1/m sum_(r=0)^(m-1) F((y+r)/m).                    (2)
```

Fix `L>=1` and put

```text
M=13^(L-1),                         N=13M=13^L.                (3)
```

Let

```text
F:T -> R,                 F bounded,
G:T -> R_(>=0),           G bounded,                           (4)
```

and define

```text
h=P_M F,                         W(y)=G(13y).                  (5)
```

All translations below are modulo one.  For `u in F_13`, define the
fixed-last-digit collision average

```text
B_u
 =1/M sum_(e=0)^(M-1) integral_T
    G(Nx) F(x) F(x+(u+13e)/N) dx.                             (6)
```

The hypotheses in (4) are enough for every identity below.  Nonnegativity
of `G` is used for positivity and the equality classification.  Rational
step hypotheses enter only in Section 5.

## 2. Collision semantics and the exact reduction

For

```text
d=u+13e,
x'=x+d/N,                                                       (7)
```

the two endpoints have the same depth-`L` future:

```text
13^L x'=13^Lx+d=13^Lx                 mod 1.                  (8)
```

If `u!=0`, then `d` is a `13`-adic unit.  The endpoints are still distinct
at time `L-1`, with oriented predecessor difference `u/13`, and first
coalesce exactly at time `L`.  If `u=0`, they already coalesce by time
`L-1`.  Thus `u=0` is the same-last-digit class, while `u!=0` records a
genuine last-step collision colour.  The address `e` records all higher
digits and is deliberately averaged, not discarded pointwise.

There is an exact one-level reduction of (6):

```text
B_u=integral_T W(y)h(y)h(y+u/13)dy.                           (9)
```

Indeed, split the `x`-circle into the `M` inverse branches

```text
x=(y+r)/M,                       r=0,...,M-1.                 (10)
```

Then `dx=dy/M`, `G(Nx)=G(13y)`, and

```text
x+(u+13e)/N=(y+r+e+u/13)/M.                                  (11)
```

For each fixed `r`, the index `r+e mod M` runs through all `M` branches.
The `e`-average of the second factor is `h(y+u/13)`; the `r`-average of the
first factor is `h(y)`.  This proves (9).

Equation (9) is also the exact toothpick recursion.  A depth-`L` pair has
been split into its last digit `u` and a depth-`L-1` stalk, and the latter is
precisely the Perron average `h=P_MF`.

## 3. Conditional variance and the `K_13` Dirichlet form

Define the last-digit conditional expectation

```text
(E_13 h)(y)
 =1/13 sum_(v in F_13)h(y+v/13)
 =(P_13h)(13y).                                               (12)
```

Both `W` and `E_13h` are invariant under `y -> y+v/13`.  Consequently
`E_13` is the orthogonal projection onto the last-digit-invariant subspace
of `L^2(W(y)dy)`: averaging over the thirteen fibre translations gives

```text
integral_T W(E_13h)(h-E_13h)=0.                               (13)
```

Using (9),

```text
B_0=integral_T Wh^2,

1/13 sum_u B_u
 =integral_T Wh(E_13h)
 =integral_T W(E_13h)^2.                                     (14)
```

Therefore the same-last-digit excess is exactly

```text
D_13(F;G,L)
 :=B_0-1/13 sum_(u in F_13)B_u

 =integral_T W(y)(h(y)-E_13h(y))^2dy
 >=0.                                                         (15)
```

This is a Pythagorean identity, not an asymptotic inequality.

To expose the graph object, put

```text
q_r(z)=h((z+r)/13),
qbar(z)=1/13 sum_r q_r(z).                                    (16)
```

Changing variables through the thirteen inverse branches of `z=13y`
turns (15) into

```text
D_13(F;G,L)
 =integral_T G(z)
    [1/13 sum_r (q_r(z)-qbar(z))^2] dz

 =integral_T G(z)
    [1/(2*13^2) sum_(r,s)(q_r(z)-q_s(z))^2] dz.               (17)
```

Thus the drift is the future-owner-weighted complete-graph `K_13`
Dirichlet energy of the thirteen predecessor-root masses.

Each individual collision residue has its own positive toothpick energy.
Fibre invariance of `W` lets us symmetrize (9):

```text
B_0-B_u
 =1/2 integral_T W(y)(h(y)-h(y+u/13))^2dy

 =1/26 integral_T G(z)
    sum_(r in F_13)(q_r(z)-q_(r+u)(z))^2dz
 >=0.                                                        (17a)
```

For `u!=0`, the edges `r--(r+u)` form one Hamiltonian `13`-cycle.
Because every nonzero `u` generates `F_13`, equality in (17a) for even one
`u!=0` forces all thirteen `q_r` to agree on `{G>0}`.  Hence

```text
D_13>0  implies  B_0>B_u for every u!=0,                      (17b)
```

without rationality.  The six antipodal pairs `{u,-u}` give six Hamiltonian
toothpick cycles which partition the edges of `K_13`; averaging their cycle
energies is exactly (17).  This is the graph-level self-similarity behind the
last-digit recursion.

The equality boundary is exact.  The following are equivalent:

```text
D_13(F;G,L)=0;                                                (18a)

h(y)=E_13h(y) for almost every y with G(13y)>0;               (18b)

for almost every z with G(z)>0,
  q_0(z)=q_1(z)=...=q_12(z);                                  (18c)

B_u=B_0 for every u in F_13.                                 (18d)
```

The implication from (18a) to (18b) uses `G>=0`.  Fibre invariance of the
weight upgrades (18b) to equality of all thirteen values in (18c), and
(9) then gives (18d).  Conversely (18d) makes the left side of (15) zero.
The zero branch therefore means exactly that the marked future owner sees
no last inverse digit.

## 4. Root-colour norms and the high-frequency ladder

Let `zeta=exp(2 pi i/13)`.  For `a in F_13`, define the fibre-character
projection

```text
H_a(y)=1/13 sum_(u in F_13)zeta^(-au)h(y+u/13)                (19)
```

and the `u`-Fourier transform

```text
Bhat(a)=1/13 sum_(u in F_13)B_u zeta^(-au).                   (20)
```

The projections satisfy

```text
H_a(y+v/13)=zeta^(av)H_a(y),
h=sum_a H_a.                                                  (21)
```

Because `W` is fibre-invariant, different fibre characters are orthogonal.
For real `F`, `H_(-a)=conjugate(H_a)`, and (9) gives

```text
Bhat(a)
 =integral_T W(y)H_(-a)(y)H_a(y)dy
 =integral_T W(y)|H_a(y)|^2dy
 >=0.                                                         (22)
```

In particular,

```text
D_13(F;G,L)=sum_(a in F_13^*)Bhat(a).                         (23)
```

With the Fourier convention

```text
Fhat(n)=integral_T F(x)exp(-2 pi i n x)dx,                    (24)
```

the Perron identity is

```text
hhat(k)=Fhat(kM).                                             (25)
```

Consequently (19) has the exact ladder expansion, in `L^2`,

```text
H_a(y)
 =sum_(k congruent a mod 13)
    Fhat(kM)exp(2 pi i k y).                                  (26)
```

The future collision colour `a` is therefore not carried by the bulk mean.
It is the weighted squared norm of the high-frequency ladder whose quotient
frequency has last base-`13` digit `a`.  For BV step data these are boundary
terms of size at most the corresponding high-frequency scale; their
nonvanishing is additional information, not a consequence of the bulk
moment.

The construction remains antipodally even.  Substitution
`y -> y-u/13` in (9) gives

```text
B_(-u)=B_u,
Bhat(-a)=Bhat(a).                                             (27)
```

Thus the twelve nontrivial colours occur in six equal converse pairs.  The
label `u` orients the two predecessor branches before contraction, but the
self-product in (6) forgets that orientation.  No tournament edge or
owner-loop direction may be inferred from (22).

## 5. Rational all-or-all law

Assume now that `F` and `G` are rational step functions.  Every integral in
(6) is rational, so

```text
B=(B_0,...,B_12) in Q^13.                                    (28)
```

If `Bhat(a_0)=0` for one `a_0!=0`, then the rational polynomial

```text
P(X)=sum_(u=0)^12 B_uX^u                                     (29)
```

vanishes at a primitive thirteenth root.  Since `Phi_13` is irreducible over
`Q` and both polynomials have degree at most twelve, `P` is a scalar multiple
of `Phi_13`.  Hence all coefficients `B_u` are equal.  Conversely a constant
vector has every nontrivial Fourier coefficient zero.

Combining this with (15), (18), and the norm formula (22) gives the sharp
dichotomy

```text
D_13(F;G,L)=0
  iff B_u is constant in u
  iff Bhat(a)=0 for every a!=0;

D_13(F;G,L)>0
  iff Bhat(a)>0 for every a!=0.                               (30)
```

Thus one positive last-digit drift simultaneously supplies all twelve
future root colours on one common rational `u`-indexed collision stalk.  It
does not choose twelve unrelated higher addresses.  Equation (27) remains
load-bearing: simultaneous colour survival is not oriented colour survival.

## 6. Sharp future-measurable replica hostile

Bulk degree-two charge does not force (30).  The failure already occurs for
a Boolean table in the exact hard branch of THM-2513.

Index a table by `(ell,s) in F_7 x F_13`.  Put

```text
a=1/2,
w_1=1/4,                         w_s=0 for s!=1,

A_(0,s)=a+w_s,
A_(ell,s)=w_s                    for ell!=0.                  (31)
```

Thus the target-zero column is the positive delta anchor

```text
A_(ell,0)=1/2 1_(ell=0),                                        (32)
```

the owner row is nonflat, and `A` is delta plus six replicas, so its first
ANOVA interaction is zero.  Its entrywise square has, for every `ell!=0`,
the nonzero anchored rectangle at `s=1`

```text
(3/4)^2-(1/4)^2-(1/2)^2+0=1/4.                               (33)
```

Hence `A^(circ 2)` has nonzero interaction and, by rational Galois
transitivity, all `72` primitive bulk mixed colours.

For each cell choose the rational Boolean interval function

```text
q_(ell,s)=1_[0,A_(ell,s)),
F_(ell,s)(x)=q_(ell,s)(Nx).                                   (34)
```

The mean of `F_(ell,s)` is exactly `A_(ell,s)`.  But

```text
(P_MF_(ell,s))(y)=q_(ell,s)(13y),                             (35)
```

which is constant on every last-digit fibre.  Equivalently, for every
integer `d`,

```text
F_(ell,s)(x+d/N)=q_(ell,s)(Nx+d)=F_(ell,s)(x).                (36)
```

Therefore for any positive rational Boolean choice of `G`, every table
`B_u` is independent of `u` and every nontrivial future collision colour
vanishes.  This hostile is nonconstant, Boolean, positively anchored, and
has complete bulk square spectrum.  What it lacks is precisely last-digit
variation.  It is an abstract future-measurable control, not asserted to be
an actual THM-2449 cover row.

Future measurability of `F` is sufficient but not necessary for zero drift.
The complete nullspace is (18): only `h=P_MF` must be last-digit invariant
on the owner support.  Components of `F` killed by `P_M` may still vary
inside the higher-address sheets.  Thus the hostile displays the mechanism
cleanly without pretending to classify its full preimage upstairs.

## 7. Lawful packet realization and the deep-sheet boundary

At a fixed lawful THM-2449 clock, write one response density before the deep
sum as

```text
F_j(x)=sum_r f_(j,r)(x),                                      (37)
```

where the finitely many `f_(j,r)` are nonnegative rational Boolean packet
factors.  The sum `F_j` can take the values `0,1,2`; no false idempotence of
that sum is used.  Let `G` be a positive rational Boolean future
owner--word event.  Refine (6) before summing:

```text
B_(u;r,r')(j)
 =1/M sum_e integral_T
    G(Nx)f_(j,r)(x)
          f_(j,r')(x+(u+13e)/N)dx.                            (38)
```

Every summand in (38) is Boolean on a finite inverse-ancestry fibre product,
and

```text
B_u(j)=sum_(r,r')B_(u;r,r')(j).                               (39)
```

Thus the rational `u`-DFT in Section 5 is taken only after a lawful
Boolean-before-Fourier construction.  If the old packet has a pointwise
deep diagonal zero on its first leg, multiplying by the second leg and by
`G(Nx)` preserves that zero in (38).  The old deep label can therefore be
retained as a sidecar.

It cannot generally be rebased at the future owner.  For the standard LRC
danger-indicator phase bank

```text
{d(cx-q/p):q in F_p},                                        (40)
```

translation by a nonzero address `d mod N` permutes the labels if and only if

```text
pcd/N is an integer.                                         (41)
```

Choose the canonical representative of `d` and put
`v=nu_13(d)`, so `0<=v<L`.  The two inverse branches first coalesce at time
`L-v`.
For the septimal source bank, (41) forces

```text
L-v<=nu_13(c),                                                (42)
```

and for the old deep thirteenth bank it forces

```text
L-v<=nu_13(c)+1.                                              (43)
```

In particular, a unit last digit `u!=0` gives `v=0`: a genuinely late
first-collision pair cannot preserve either old bank once `L` exceeds its
fixed valuation horizon.  Making `d` highly divisible restores a bank only
by moving the actual first collision back to that bounded horizon.  This is
an exact scale/sheet obstruction.

For a full `F_7 x F_13` response bank, (30) applies entrywise and to every
nonnegative rational combination of cells.  A single drifting marked cell,
zero-extended to the table, has every primitive table character because its
cell label is a delta mark.  That is a chart-marked statement, not an
intrinsic mixed interaction: signed ANOVA combinations can still cancel,
and a positive bulk mixed moment does not force entrywise drift by Section
6.  Proving last-digit drift on an actual live table is therefore a new,
precise obligation rather than a consequence already supplied by
THM-2513.

## 8. Relation to the cyclic `K_14` chart and final nonclaims

Equation (17) presents a genuine graph object.  At each future owner point
there are thirteen predecessor-root vertices, and their scalar masses carry
the `K_13` Dirichlet energy.  Adjoining the marked future owner gives a
natural set of fourteen semantic vertices.

This is not the vertex identification of THM-2514.  In THM-2514 the thirteen
finite vertices are table target/root labels and `Omega` resolves diagonal
pairs in a signed static edge chart.  Here the thirteen vertices are
physical inverse branches immediately before one future collision, and the
marked point is that future owner.  No theorem identifies these two
`F_13`-torsors, their affine gauges, their deep sheets, or their edge
weights.  Numerical equality of labels or characters is not such a map.

The proved gain and remaining boundary are therefore:

- last-digit drift is an exact nonnegative conditional variance and
  `K_13` energy;
- for rational data, positive drift gives all twelve future root colours on
  one common Boolean-before-DFT collision stalk;
- the quadratic self-correlation identifies `u` with `-u` and supplies no
  oriented owner-loop current;
- old deep data survives on the two ancestry legs but cannot generally
  descend to the future owner without its sheet residues;
- bulk first/square moment charge can lie entirely in the zero-drift future
  factor; and
- no signed degree has been converted into one runner-labelled Boolean
  event, no live row is excluded, and LRC(14) remains open.

The next decisive test is whether the actual lawful owner-word response has
positive last-digit variance on its owner support, or whether its zero branch
can be shown to force a previously excluded future-factor form. **QED.**

## 9. Exact finite-tower referee

Run

```bash
python3 04-computation/lrc14_perron_collision_cospan_thm2518_2519_referee.py
python3 -O 04-computation/lrc14_perron_collision_cospan_thm2518_2519_referee.py
```

Both runs reproduce the stored transcript byte-for-byte.  On the exact
`C_845 -> C_5` tower, the referee checks all thirteen fixed-last-digit
averages, `B_u=B_(-u)`, first-collision timing, and the identical values

```text
B_0-average_u B_u
 =conditional variance
 =K_13 Dirichlet energy
 =2638/142805.
```

It also checks a future-measurable Boolean hostile for which all `169` deck
needles are exactly `1/5`, so the drift vanishes.

An independent line-by-line audit rederived the weighted projection identity,
the Fourier normalizations and frequency ladder, the rational all-or-all law,
the hostile family, and the packet/deep scope.  It also forced the phase-bank
iff above to be stated only for the standard LRC danger indicator and the
valuation formula to exclude the diagonal address `d=0`.
