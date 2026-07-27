---
id: THM-2533
title: "Owner-weighted collision phase and mixed gain Radon ladders"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Multiplication by a
  1/13-periodic future owner commutes with
  every last-digit projector.  Hence the positive owner-weighted collision
  energies of THM-2522 are exactly squared norms of linear owner-attached
  root-character pieces.  On every one of the 165 live THM-2349 rows, one
  sufficiently late Boolean owner block therefore leaves a single lawful
  Boolean event with all twelve nontrivial mod-13 Fourier ladders.  Every
  ladder has a positive integer Fourier lift bounded by its finite grouped
  jump count.  Reapplying THM-2349's abstract carrier lemma to that same
  owner-marked subevent gives all twelve shallow colours of 91-unit marked
  deepest-comb triangles at one common future owner.  Separately, a lawful
  THM-2449/2513 first-or-delayed-square
  mixed table channel has at least four depth-one projections including
  its mean, hence at least three distinct gains and at least 216
  Galois-saturated owner-attached phase ladders.  If exactly four
  projections occur, the channel is confined to the guard toothpick stalk
  C_H intersect T^(-1)E_H; any mass off that stalk forces a fourth gain.
  For each fixed mixed channel the phase ladders assemble into one physical
  13-root fibre; every chosen cyclic-tournament Cayley derivative is a
  nonzero owner-attached odd linear signal, and its sawtooth inverse recovers
  the centred fibre.  The three proved gain ratios therefore give three
  slope-labelled, nonzero aligned-Radon-difference currents on each of the
  72 mixed channels.  On the Boolean eventwise branch the same predecessor
  potential has all 72 pointwise mixed strip modes and all 5,184
  Hilbert-valued cut modes.  At the common owner time, either the full
  gain-energy vector is flat and
  all twelve gains give 864 phase ladders, or every cyclic-tournament Cayley
  chart has a nonzero oriented gain-energy contrast; this latter contrast is
  a separate quadratic gain-index diagnostic.
  Honest Boolean branch intervals show that the full diagonal energy bank
  cannot recover root address.  The theorem does not select a terminal
  component, relation-lattice address, common lift across colours, corolla
  gain, target/chronological interpretation of its root orientation, scalar
  row exclusion, or LRC(14).
source: codex-2026-07-27-owner-weighted-collision-phase
depends_on:
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2388-thirteen-root-multiplicity-reflection-and-blocker-caged-toothpick-law
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2513-anchored-first-or-second-moment-spectrum-and-pair-space-boundary
  - THM-2522-intrinsic-collision-depth-toothpick-descent-and-late-owner-decoupling
  - THM-2532-cyclic-tournament-cayley-algebra-and-chi7-even-quotient
related:
  - MISTAKE-279
  - THM-2286-endpoint-prony-lift-bank-and-sharp-owner-multiplier-landing
  - THM-2312-sparse-root-bispectrum-positive-word-current
  - THM-2323-primitive-fixed-colour-cross-correlation-and-same-gauge-word-alignment
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2457-complete-atom-root-cosupport-graph-and-semantic-word-hostile
  - THM-2474-squarefree-first-collision-primitive-character-saturation
  - THM-2524-translated-chi7-hamilton-polarization-inversion
  - THM-2525-unit-guard-collision-floor-and-word-owner-cross-cospan-collapse
  - THM-2526-affine-skew-orientation-gauge-boundary
  - THM-2527-owner-weighted-all-mode-odd-bank-and-boolean-cut-coordinate
  - THM-2528-intrinsic-four-arm-boolean-path-and-joint-autocorrelation-scalarization
  - THM-2531-prime-necklace-guard-boundary-selector
script: 04-computation/lrc14_owner_weighted_collision_phase_thm2533.py
output: 05-knowledge/results/lrc14_owner_weighted_collision_phase_thm2533.out
script_sha256: b811ed197e3e29af3755c22c08113f66fdb41dafa75baba3b7d7092a68080f41
output_sha256: add3a2d638a29d55c2ffbf8fd6ef2c778fb2feecb899818c5b50fa877a625096
hash_basis: working-tree bytes (LF)
---

# THM-2533 -- collision energy lifts to owner-attached phase

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2522 proves two facts at the intrinsic first collision of every live
THM-2349 owner--word event:

```text
all twelve last-digit energies are positive;

one future owner block may be placed arbitrarily later.                 (1)
```

Those are quadratic statements.  The elementary but load-bearing point of
this theorem is that the future block is constant on each depth-one root
fibre.  It can therefore be moved *inside* every last-digit projector.  A
positive collision energy then certifies a nonzero **linear** phase piece of
the actual owner-attached event.

There are two applications which must remain typed separately:

```text
one THM-2349 Boolean event:
  all twelve root colours on the same lawful event;

one THM-2449/2513 mixed table character:
  at least three collision/target gain ratios on a complex channel.     (2)
```

The first line has no target character and hence no gain.  The second line
is a character sum, not one Boolean event.  Conflating them would turn the
stronger eventwise colour count into a false all-gain table theorem.

## 1. The invariant-multiplier module lemma

Let `zeta=exp(2 pi i/13)`.  For `a in F_13`, define the last-digit
projector on `L^2(T)` by

```text
(E_a f)(x)
 =1/13 sum_(u in F_13) zeta^(-au)f(x+u/13).                  (3)
```

The projectors are pairwise orthogonal, sum to the identity, and satisfy

```text
(E_a f)(x+v/13)=zeta^(av)(E_a f)(x).                          (4)
```

Let `g` be any bounded function invariant under translation by `1/13`.
Then

```text
E_a(gf)=g E_a f                                  for every a. (5)
```

Indeed, `g(x+u/13)=g(x)` in every summand of (3).  Thus the invariant
functions form the scalar ring over which the thirteen character spaces
are modules.

For `g>=0`, put

```text
Gamma_a(g,f)=integral_T g(x)|E_a f(x)|^2dx.                   (6)
```

Equation (5) gives the exact implications

```text
Gamma_a(g,f)>0  ->  E_a(gf)=gE_af is not zero;                (7)

g in {0,1} almost everywhere
  -> Gamma_a(g,f)=||E_a(gf)||_2^2.                            (8)
```

The Boolean hypothesis is used only for the equality in (8).  For an
arbitrary nonnegative weight, positivity in (6) still implies (7).

If

```text
g_R(x)=G(13^R x),                         R>=1,                (9)
```

then `g_R` is `1/13`-periodic.  Hence every sufficiently delayed owner
weight in THM-2522 is an allowed scalar in (5).

There is a useful Fourier reformulation.  With

```text
fhat(n)=integral_T f(x)exp(-2 pi i n x)dx,                    (10)
```

equation (4) says

```text
(E_a f)^hat(n)=fhat(n),                 n=a mod 13,
              =0,                       otherwise.            (11)
```

Consequently (7) says that a positive owner-weighted collision colour is
not merely an energy: the owner-attached function `gf` has a nonzero
ordinary Fourier coefficient in that exact residue class.

## 2. A sharp grouped-jump Prony bound

The preceding conclusion has a finite positive-frequency lift.  Let `W`
be a complex-valued periodic step function with finitely many **true net
jumps**.  Artificial repeated endpoints with total jump zero are combined
before the following definition.

Fix `a in F_13^*`.  Put two jump locations in the same projected class when

```text
x~y  iff  13(x-y) is an integer.                              (12)
```

For a class `c`, define its charged grouped jump

```text
A_(a,c)
 =sum_(x in c intersect Jump(W))
    Delta W(x) exp(-2 pi i a x),                              (13)
```

and let

```text
L_a(W)=#{c:A_(a,c)!=0}.                                      (14)
```

Choose one representative `x_c` of every active class.  Distributional
integration by parts gives, for `n=a+13h`,

```text
2 pi i(a+13h) What(a+13h)
 =sum_c A_(a,c)
    [exp(-2 pi i 13x_c)]^h.                                  (15)
```

Changing the representative in (15) changes neither the node nor the
grouped coefficient because the compensating thirteenth-root phase has
already been included in (13).

Suppose `E_aW!=0`.  Then `L_a(W)>=1`.  Otherwise (15) would make every
Fourier coefficient of `W` in the residue class `a` vanish; since
`a!=0`, there is no zero-frequency exception, contradicting (11).

The nodes in (15) are distinct and nonzero.  A nonzero exponential sum on
`L=L_a(W)` nodes cannot vanish at `L` consecutive integers: on such a
block, remove the nonzero starting powers and invert the ordinary
Vandermonde matrix.  Therefore, for every integer `H_0`, some

```text
H_0<=h<=H_0+L_a(W)-1
```

satisfies

```text
What(a+13h)!=0.                                              (16)
```

Taking `H_0=0` yields an exact positive lift

```text
n_a=a+13h,

1<=n_a<=a+13[L_a(W)-1]
       <=13L_a(W)-1
       <=13J(W)-1,                                           (17)
```

where `J(W)` is the number of true net jumps.  The stronger block statement
in (16) says that the positive lifts are syndetic in the gauge index.

This is a one-function version of the endpoint-Prony mechanism of
THM-2286.  If one also retains a nonzero inner product

```text
integral_T W conjugate(f)!=0,                                 (18)
```

Parseval gives a common frequency of `W` and `f`.  Multiplying their two
endpoint-current formulas gives the THM-2323 bound with at most
`J(W)J(f)` endpoint-difference nodes.  We will use only the sharper
one-function bound (17) for the owner-attached phase.

## 3. Eventwise consequence on all 165 live rows

Fix any one of the `165` live first-depth-one rows.  THM-2349 supplies a
shallow owner `j`, a nonempty terminal word `sigma`, and a finite clock `k`
such that the rational Boolean event

```text
F(x)=1_(E_j)(x)1_(Q_(j,sigma))(13^k x)                       (19)
```

has positive measure.  THM-2522 proves uniformly that its intrinsic first
digit is level zero:

```text
D_0=(I-E_0)F!=0,                L_*=1.                        (20)
```

Here `E_0` in `(I-E_0)` is the invariant projector from (3), and THM-2522's
nonzero last-character pieces are precisely

```text
D_(0,a)=E_aF,                         a in F_13^*.             (21)
```

Because `F` is rational-valued, cyclotomic irreducibility gives the exact
all-or-all law

```text
E_aF!=0                      for every a!=0.                  (22)
```

Now let `G:T->{0,1}` be any rational Boolean BV future owner or complete
owner--word block of positive mean

```text
rho=integral_T G>0.                                          (23)
```

Place it `R` steps after the depth-one collision and put

```text
g_R(x)=G(13^(R+1)x),
W_R(x)=F(x)g_R(x).                                           (24)
```

The exponent `R+1` is THM-2522's exact convention: the collision itself is
at time one and `R` is a later wait.  One common sufficiently large `R`
makes every owner-weighted colour positive.  More explicitly, with

```text
e_0=min_(a!=0)||E_aF||_2^2>0,
B=||F||_infinity,
V=Var(F),
W=Var(G),                                                    (25)
```

THM-2522 equation (55) permits any `R` satisfying

```text
13^(R+1)>BVW/(6rho e_0).                                     (26)
```

At such a delay, equations (8), (21), and THM-2522's weighted collision
identity give, simultaneously for all `a!=0`,

```text
Bhat^[0,R](a)
 =integral_T g_R|E_aF|^2
 =||E_aW_R||_2^2
 >0.                                                        (27)
```

Thus one and the same lawful Boolean event `W_R` has

```text
E_aW_R!=0                         for all twelve a!=0.        (28)
```

Apply Section 2 to `W_R`.  For every `a!=0`, there is a positive integer

```text
n_a=a mod 13,
1<=n_a<=13L_a(W_R)-1<=13J(W_R)-1,                            (29)

(W_R)^hat(n_a)!=0.                                           (30)
```

The twelve integers are distinct because their residues are distinct.
There is no assertion that their gauge quotients `(n_a-a)/13` agree.

This is the direct all-colour Boolean phase conclusion.  The event in (24)
retains literally:

- the THM-2349 source owner `j`;
- its selected terminal word `sigma` and clock `k`;
- the one common later owner--word block `G`;
- the physical first collision `L=1`; and
- the old shallow septimal and old deep `c_3` banks, each permuted by the
  unit depth-one translation as in THM-2522.

The new residue `a`, the shallow septimal label, and the old deep label are
three different typed coordinates.  Equation (30) does not identify them.

### The same future owner marks all twelve 91-unit triangles

The Boolean event in (24) has a second consequence which uses its carrier,
not merely its last-digit projections.  Retain the THM-2349 shallow blocker
`c_j` and deepest blocker `c_3`.  Up to the usual null endpoints, its
inherited geometry is

```text
nu_13(c_j)=1<nu_13(c_3),
support F subset E_j subset D_(c_j),
support F intersection D_(c_3)=empty.                       (30a)
```

Moreover `0<=W_R<=F`, and (28) makes `W_R` nonzero.  Both functions are
rational step functions.  Thus the abstract shallow-carrier triangle of
THM-2349 Section 1 applies anew with

```text
e=F,                    f=W_R.                               (30b)
```

For every prescribed `kappa in F_13^*`, it supplies integers
`X_kappa,Y_kappa,m_kappa` such that

```text
Y_kappa=X_kappa+m_kappa c_3,
gcd(m_kappa,91)=1,

nu_13(X_kappa)=nu_13(Y_kappa)=1,
X_kappa/13=Y_kappa/13=kappa mod 13,                          (30c)

(W_R)^hat(X_kappa)Fhat(X_kappa)Fhat(Y_kappa)!=0,             (30d)

(W_R)^hat(X_kappa)
 (1_(D_(c_3)))^hat(m_kappa c_3)
 conjugate(Fhat(Y_kappa))!=0.                               (30e)
```

The inherited finite invoice is

```text
|m_kappa|<=13J(F)J(W_R)+7J(F)-13.                            (30f)
```

Consequently the **same** future owner block marks twelve 91-unit
deepest-comb edge/triangle incidences, one in every nonzero shallow root
colour.  At `X_kappa`, the coefficient of `W_R` records both the original
terminal word and the future owner, while the accompanying coefficient of
`F` records the original terminal word before that owner was imposed.

This incidence is ordered by the owner-marked endpoint `X_kappa` versus the
bare endpoint `Y_kappa`; no assertion says that `(W_R)^hat(Y_kappa)=0`.
It is still target-neutral: it neither chooses
a terminal component nor identifies a relation-lattice phase, and it does
not orient `m_kappa` against `-m_kappa` in a quotient which forgets which
endpoint carried the future owner.

There is also a general analytic version.  At any active rational THM-2522
digit `m`, take `f=F_m=P_(13^m)F` and choose `R` sufficiently large for
that digit, with the same-form weight `g_R(y)=G(13^(R+1)y)`.  Equations (5)
and (50) of THM-2522 give twelve nonzero owner-attached phase pieces on
`g_RF_m`.  Only `m=0` in (19) is literally one Boolean event; at `m>0`,
`F_m` is a rational averaged response.

## 4. The lawful mixed-table channel

The eventwise conclusion has no target colour and therefore no gain ratio.
We now construct a different object which retains source and target
characters.

Choose one fixed lawful THM-2449 clock after the THM-2513
eventual-overlap threshold, within one fixed word/parity class, and write
the rational response densities

```text
F_(ell,s)(x)=sum_r h_(ell,r,s,0)(x),
A_(ell,s)=integral_T F_(ell,s),                               (31)
```

for `(ell,s) in F_7 x F_13`.  Every density is a finite sum of lawful
Boolean deep-labelled events.  The target covector moves one ordinary unit
and its blocker; the source skew moves the source blocker and its transported
word.  The guard-safe factor is common to every cell.  Put

```text
E_H={x:||Hx||<1/7},
C_H=T minus E_H,                     13 does not divide H.    (32)
```

Up to null endpoints,

```text
F_(ell,s)=F_(ell,s)1_(C_H)              for every ell,s.      (33)
```

There is also at least one stationary ordinary unit-safe factor.  Choose an
ordinary unit `q_k` unaffected by the source and target moves.  With

```text
D_q={x:||qx||<1/14},
g_safe(x)=1-1_(D_1)(x),
```

one has

```text
13 does not divide q_k,
F_(ell,s)(x)=F_(ell,s)(x)g_safe(q_kx)   for every ell,s.      (34)
```

Fix one mixed character `(kappa_0,b_0) in F_7^* x F_13^*`.
THM-2513 proves that at least one of the first and entrywise-square table
coefficients

```text
Ahat(kappa_0,b_0),
(A^(circ 2))hat(kappa_0,b_0)                              (35)
```

is nonzero.  We realize the successful branch by one finite rational step
channel.

If the first coefficient is nonzero, put

```text
S_(ell,s)(x)=F_(ell,s)(x),                    Q=0.             (36)
```

Otherwise choose a finite even `Q>=2` and put

```text
S_(ell,s)(x)
 =F_(ell,s)(x)F_(ell,s)(13^Qx).                             (37)
```

Such a common finite `Q` exists.  The two-BV mixing estimate gives

```text
|integral F_(ell,s)(x)F_(ell,s)(13^Qx)dx-A_(ell,s)^2|
 <=Var(F_(ell,s))^2/(12*13^Q).                               (38)
```

After the normalized `91`-entry character transform, the total error is at
most

```text
[sum_(ell,s)Var(F_(ell,s))^2]/(1092*13^Q).                  (39)
```

Choose the next even `Q` for which (39) is smaller than half the modulus of
the nonzero second value in (35).  Evenness retains the lawful source
subsequence `13^Q=1 mod 7`.  Therefore the mixed channel

```text
R_(kappa,b)(x)
 =1/91 sum_(ell,s)
   S_(ell,s)(x)zeta_7^(kappa ell)zeta_13^(bs)                 (40)
```

satisfies

```text
integral_T R_(kappa_0,b_0)!=0.                               (41)
```

Every channel in (40) still obeys the two common support laws (33)--(34).
In the square branch it also obeys

```text
support R_(kappa,b) subset C_H intersection (T^Q)^(-1)(C_H), (42)
```

where `T^Qx=13^Qx` and the superscript `-1` denotes set-theoretic
preimage, not a negative iterate.

## 5. Four projections, three gains, and the toothpick stalk

For the selected channel abbreviate

```text
R=R_(kappa_0,b_0),
K_a=E_aR,
mathcal A={a in F_13:K_a!=0},
k=#mathcal A.                                                 (43)
```

Fourier reconstruction on one depth-one fibre is

```text
R(x+u/13)=sum_(a in mathcal A)zeta_13^(au)K_a(x).             (44)
```

For generic `x`, the guard-forbidden root set is

```text
Z_x={u:x+u/13 in E_H}.                                       (45)
```

THM-2388 equation (7), applied over the parent `13x`, gives

```text
#Z_x=4-1_(E_H)(13x) in {3,4}.                                (46)
```

After the unit relabelling `v=Hu`, the elements of `Z_x` are cyclically
consecutive.

Suppose `k<=3`.  At every generic `x`, choose `k` consecutive forbidden
roots.  Equation (33) makes the left side of (44) zero there.  The resulting
square coefficient matrix is, up to nonzero column factors,

```text
(zeta_13^((aH^(-1))i))_(0<=i<k, a in mathcal A),              (47)
```

an ordinary Vandermonde on distinct thirteenth roots.  It is invertible,
so every `K_a(x)` vanishes almost everywhere, contradicting (41).  Hence

```text
k>=4.                                                        (48)
```

The mean projection is active:

```text
integral_T K_0=integral_T R!=0.                              (49)
```

Thus at least three distinct nonzero indices occur.  Choose

```text
a_1,a_2,a_3 in mathcal A intersection F_13^*,
lambda_i=a_i/b_0 in F_13^*,                    i=1,2,3.       (50)
```

The three gains are distinct.

The equality case in (48) has a sharper geometric description.  Assume
`k=4`.  At a generic `x` with `13x notin E_H`, equation (46) supplies four
consecutive forbidden roots.  The `4 by 4` version of (47) forces every
`K_a(x)` to vanish.  Reconstruction (44) then makes the entire root fibre
zero.  Together with (33), this proves the **toothpick-stalk confinement**

```text
support R
 subset C_H intersection T^(-1)(E_H)              almost everywhere.
                                                                    (51)
```

In the square branch, combine (42) and (51):

```text
support R
 subset C_H intersection T^(-1)(E_H)
              intersection (T^Q)^(-1)(C_H).                   (52)
```

Equivalently, if

```text
measure({x:R(x)!=0 and 13x in C_H})>0,                        (53)
```

then `k>=5`, so at least **four** distinct nonzero gains survive.  The
minimal three-gain case can occur only on a stalk which is guard-safe now
but guard-dangerous after one toothpick step.  Equation (51) does not prove
that this stalk is empty for a lawful table.

The lower bound is sharp for the support information used here.  Choose a
small rational base interval `I` about zero so that, for `y in I`,

```text
guard-danger roots:       Hu in {-1,0,1};
every stationary ordinary danger root: u=0.                  (54)
```

Let

```text
P(X)=(X-zeta_13^(-1))(X-1)(X-zeta_13),                        (55)
```

and define a cyclotomic BV step function on the root chart by

```text
R_*((y+u)/13)=1_I(y)P(zeta_13^(Hu)).                          (56)
```

It vanishes on every guard-danger root and on the common ordinary-danger
root.  Hence it obeys all the support constraints in (33)--(34), and its
support lies in the stalk (51).  Expanding (55) shows that all four
coefficients are nonzero, so its active projection set is exactly

```text
{0,H,2H,3H}.                                                  (57)
```

Its integral is the constant coefficient of `P` times `measure(I)`, namely
`-measure(I)!=0`.  Thus support geometry alone cannot improve (48).  The
control (56) is cyclotomic and signed; it is not asserted to be a lawful
nonnegative anchored table channel.

The sharpness survives rationality and nonnegativity at the **table** level.
Write the values in (56), in the standard cyclotomic basis, as

```text
R_*(u)=sum_(s=0)^11 c_(u,s)zeta_13^s,
c_(u,12)=0.                                                  (57a)
```

Every `c_(u,s)` is rational, and all are zero on the three forbidden roots.
On the other ten roots choose one rational constant `B` larger than the
absolute value of every negative `91c_(u,s)`.  Define a nonnegative rational
root-by-source-by-target table by

```text
A_(ell,s)(u)
 =B 1_(u safe)+91 1_(ell=0)c_(u,s).                          (57b)
```

The uniform background cancels in every nonzero source character.  For
`kappa,b!=0`, the normalized mixed transform of (57b) is exactly

```text
1/91 sum_(ell,s) A_(ell,s)(u)
  zeta_7^(kappa ell)zeta_13^(bs)
 =sigma_b(R_*(u)).                                           (57c)
```

It therefore has precisely the four collision modes

```text
{0,bH,2bH,3bH},                                             (57d)
```

and precisely the three gain ratios `{H,2H,3H}`.  Across all nonzero
`(kappa,b)` this is exactly `216` mixed modes.  The table vanishes on the
common guard-danger block, including the common ordinary-danger root, and
is nonnegative rational everywhere.  It is still not a lawful THM-2449
anchored table: the delta anchor, diagonal packet zeros, and continuation
semantics have not been imposed.  Hence any improvement beyond three gains
must use those lawful structural coordinates, not rationality,
nonnegativity, or the common support factors alone.

## 6. One late owner gives at least 216 phase ladders

Let `G:T->{0,1}` again be a rational Boolean BV owner--word block of
positive mean.  Choose one `R_own>Q` sufficiently large that

```text
Gamma(kappa,b,a)
 :=integral_T G(13^(R_own)x)|E_aR_(kappa,b)(x)|^2dx,

Gamma(kappa_0,b_0,a)>0,
                    a in {0,a_1,a_2,a_3}.                     (58)
```

One common delay exists by the BV mixing estimate, since there are only
four positive norms.  Put

```text
g_own(x)=G(13^(R_own)x),
W_(kappa,b)=g_own R_(kappa,b).                                (59)
```

The factor `g_own` is `1/13`-periodic.  The module lemma gives

```text
E_(a_i)W_(kappa_0,b_0)=g_own K_(a_i)!=0.                     (60)
```

Expanding (58) before the character sums shows that every scalar integral
lies in the real part of the relevant cyclotomic field.  It need not be
rational: even two equal rational branch cells can contribute a factor
`|1+zeta_13|^2`.  For independent units `v_7 mod 7` and `v_13 mod 13`,
Galois covariance is nevertheless exact:

```text
sigma_(v_7,v_13)Gamma(kappa,b,a)
 =Gamma(v_7kappa,v_13b,v_13a).                               (61)
```

Here the automorphism is not being asked to preserve an order on a
cyclotomic field.  The positive seed in (58) is nonzero, so every Galois
conjugate is nonzero.  The right side of (61) is independently an integral
of an absolute square against the nonnegative owner and is therefore
nonnegative in the distinguished real embedding.  It is consequently
strictly positive.  Nonzero preservation plus independent nonnegativity,
not rationality, is the load-bearing propagation step.

The ratio `a/b` is fixed.  Equations (50), (58), and (61) prove

```text
Gamma(kappa,b,lambda_i b)>0
   for every kappa in F_7^*, b in F_13^*, i=1,2,3.            (62)
```

These are three disjoint `72`-element character-diagonal gain slices.  By
(5), each positive entry linearizes:

```text
E_(lambda_i b)W_(kappa,b)!=0
   for every kappa,b!=0 and i=1,2,3.                          (63)
```

Thus at least

```text
3*6*12=216                                                   (64)
```

distinct owner-attached mixed-channel/phase-ladder incidences occur at one
common `Q` and one common owner time.  Section 2 supplies, separately for
each triple, a positive integer

```text
n_(kappa,b,i)=lambda_i b mod 13,

1<=n_(kappa,b,i)
 <=13L_(lambda_i b)(W_(kappa,b))-1                            (65)
```

with

```text
(W_(kappa,b))hat(n_(kappa,b,i))!=0.                           (66)
```

No common gauge quotient `h` across the `216` channels is asserted.

The function in (59) is a complex character sum, not one event.  Still,
(66) has a literal Boolean consequence.  In the first branch, expand

```text
S_(ell,s)=sum_r h_(ell,r,s,0).
```

In the square branch, expand (37) into the finite products

```text
h_(ell,r,s,0)(x)h_(ell,r',s,0)(13^Qx).                       (67)
```

Multiplication by `g_own` makes each summand a Boolean owner-attached event.
If the finite character-weighted sum (66) is nonzero, at least one such
Boolean deep-labelled summand has a nonzero Fourier coefficient at the
**same** integer frequency.  This selects existentially:

- one table cell `(ell,s)`;
- one deep label `r`, or two labels `(r,r')` in the square branch; and
- the same future owner block and positive frequency.

The selection is not canonical.  It need not preserve the other `215`
channels, one prescribed complete atom, or the Galois orbit after the term
has been chosen.

The square branch stays physically shallow.  Since `Q>=1`, translation by
`u/13` leaves the far factor fixed:

```text
F_(ell,s)(13^Q(x+u/13))=F_(ell,s)(13^Qx).                    (68)
```

The late owner is fixed for the same reason.  Hence the collision in
(58)--(63) remains at depth one.  The near lawful packet retains its
existing shallow and deep sidecars under the same label-permutation
criterion as THM-2522; neither the far square factor nor the owner creates a
new shallow/deep identification.

### The phase pieces form one lossless oriented root derivative

The linearized pieces in (63) are not merely unrelated scalar witnesses.
For any periodic function `f`, let root translation by `tau` be

```text
(P_tau f)(x)=f(x+tau/13),
```

and apply THM-2532's cyclic-tournament polynomial on this physical root
fibre:

```text
mathcal C_tau f
 =sum_(d=1)^12 (-1)^(d+1)f(x+d tau/13).                     (68a)
```

On the `a`-th projector piece,

```text
E_a(mathcal C_tau f)
 =[(zeta_13^(a tau)-1)/(zeta_13^(a tau)+1)]E_af.             (68b)
```

The multiplier is nonzero for every `a,tau!=0`.  Moreover every
`1/13`-periodic owner factor commutes with (68a):

```text
mathcal C_tau(gf)=g mathcal C_tau f.                         (68c)
```

Consequently:

- for each eventwise Boolean `W_R` in (24),
  `mathcal C_tau W_R!=0` for every `tau!=0`; the alternating root sum keeps
  the same terminal-word factor and future owner on every term, while the
  present shallow-owner root is the varying signal;
- for every fixed mixed `(kappa,b)`, the thirteen translates of
  `W_(kappa,b)` in (59) form one common root-fibre vector, and
  `mathcal C_tau W_(kappa,b)!=0` for every `tau!=0`; in particular this holds
  on all `72` channels in the three gain slices (62)--(63).

Let `mathcal B_tau` be the sawtooth polynomial of THM-2532.  Its inverse law
on each fibre is

```text
f-E_0f=mathcal B_tau mathcal C_tau f.                        (68d)
```

Thus the pair `(E_0f,mathcal C_tau f)` loses no root-fibre information.
In fact this root fibre is **exactly** the THM-2521 predecessor potential,
not merely an analogous representation.  Write, away from endpoints,

```text
q_r^f(z)=f((z+r)/13),
mu_f(z)=1/13 sum_r q_r^f(z),
p_r^f(z)=q_r^f(z)-mu_f(z).                                  (68e)
```

Then directly from the definitions,

```text
(f-E_0f)((z+r)/13)=p_r^f(z),

(mathcal C_tau f)((z+r)/13)=(C_tau p^f(z))_r.               (68f)
```

The first line is THM-2521 equation (16), and its signless map `S(p^f)` is
the corresponding `K_14` potential edge weighting.  The second line is
THM-2532's aligned-Radon difference on that same profile.  Thus (68d) is
literally its sawtooth inversion, and the phase lift has supplied the
previously missing common linear input.

For the Boolean event (24), THM-2349 has `k>=2`, and the root coordinates
factor more sharply:

```text
q_r^(W_R)(z)
 =1_(E_j)((z+r)/13)
   Q_(j,sigma)(13^(k-1)z)G(13^Rz).                          (68g)
```

The integer `r` disappears from the last two factors.  Hence all thirteen
predecessor coordinates share the literal terminal word and future owner;
only the present shallow-owner root varies.  This gives a real Boolean
owner-attached input before centring.  For the mixed channel, (68g) is a
complex character-summed predecessor profile whose individual summands are
Boolean as in (67).

The gain labels now give an exact slope-labelled Radon lift, although the
chart is not yet physically canonical.  For every `kappa in F_7^*`,
`b in F_13^*`, and `i in {1,2,3}`, use the incidence
`a=lambda_i b` from (62)--(63), fix `f=W_(kappa,b)`, and choose the affine
chart `tau=lambda_i`.  Equations (63) and (68b) give

```text
E_(lambda_i b)(mathcal C_(lambda_i)W_(kappa,b))
 =[(zeta_13^(lambda_i^2 b)-1)
    /(zeta_13^(lambda_i^2 b)+1)]
   E_(lambda_i b)W_(kappa,b)
 !=0.                                                       (68h)
```

By (68f), the function in (68h), read on its thirteen roots, is exactly

```text
(mathcal R_(lambda_i)-mathcal R_(-lambda_i))p^(W_(kappa,b))
 =C_(lambda_i)p^(W_(kappa,b)).                              (68i)
```

Thus every mixed `(kappa,b)` channel carries three owner-attached,
slope-labelled, nonzero THM-2521 aligned-Radon-difference currents, for
`3*6*12=216` proved slope--channel incidences.  This is an exact linear
predecessor statement, not merely the index coincidence `a=b tau`.  It
still does not identify the affine root chart with THM-2508's ordered
cut/carry sheet or turn the signed current into a Boolean terminal row.

There is a simultaneous eventwise consequence.  For `f=W_R`, the vector
`q^f(z)` is rational and Boolean.  Equation (28) makes `p^f` nonzero on a
set of positive measure, so THM-2521 Section 4 applies in every chosen
affine `K_14` chart, on each nonuniform root fibre.  Its full signless-
potential/cut bundle therefore has

```text
12*6=72       pointwise nonzero mixed strip modes,
12*6*12*6=5,184 nonzero Hilbert-valued cut modes.            (68j)
```

All of them retain the common terminal-word and future-owner factors in
(68g).  This upgrades the old abstract potential bridge to an actual
Boolean input before centring, but not to a Boolean output after the
signed potential and cut transforms.

The precise scalar seam is already visible in the root derivative.  Its
twelve coefficients sum to zero, so translation invariance of Haar measure
gives, for every integrable `f`,

```text
integral_T mathcal C_tau f=0.                               (68k)
```

Hence untwisted scalar averaging erases each `mathcal C_tau f`.  It does
not erase its positive-frequency ladders: if `n=a mod 13` with `a!=0`, then

```text
(mathcal C_tau f)^hat(n)
 =[(zeta_13^(a tau)-1)/(zeta_13^(a tau)+1)]fhat(n),          (68k')
```

so every Prony witness from (29)--(30) remains nonzero.  A closing scalar
consumer must pair this charged mode with a nonconstant target, arrival, or
terminal twist before integration; no claim is made that every function in
the larger 72/5,184 cut bundle has zero integral.

The resulting `C_tau p^f` and `S(p^f)` are nevertheless signed/complex
linear combinations, not one Boolean cospan or scalar row.  For the live
eventwise branch, THM-2526's inherited guard sheet supplies the physical
root sign `tau_H`; for the three mixed gain slopes no theorem identifies
`lambda_i` with that live sign.  In neither case does the root sign yet map
this potential to a target-active THM-2334 terminal relation current or a
chronological owner loop.

### All gains or an oriented gain-energy signal

At the same common owner time, retain the entire nonnegative vector

```text
gamma_a=Gamma(kappa_0,b_0,a),                  a in F_13.     (68l)
```

Equation (58) gives `gamma_0>0`.  Apply the cyclic-tournament Cayley
derivative `C_tau` of THM-2532 to the gain-index coordinate.  Its kernel on
`F_13` is exactly the constants.  Therefore exactly one of the following
two alternatives holds:

```text
gamma is constant:
  gamma_a=gamma_0>0 for every a;
  all twelve nonzero gains survive;
  Galois propagation gives 12*6*12=864 phase ladders;

gamma is nonconstant:
  C_tau gamma!=0 for every tau in F_13^*.                    (68m)
```

The second line is an oriented gain-energy signal: replacing `tau` by
`-tau` reverses its sign.  It does not orient a physical owner loop.
The entries of `gamma` are quadratic energies on collision-character
channels.  Equations (68f)--(68m) now identify their underlying linear root
fibre with a genuine THM-2521 predecessor potential.  What remains absent is
an identification of the **gain-index** action on `gamma` with the
predecessor-coordinate Radon action on `p^f`, and then with a Boolean
target/terminal current.  Thus (68m) is an exact gain-index diagnostic, not
that semantic transport theorem.

## 7. Honest phase-retrieval hostile

The passage from diagonal energy to the charged function in (5) is
essential.  Even the complete diagonal bank does not determine root
address.

Fix a generic rational base point `y_0`.  On its thirteen-root fibre, the
guard-safe set has `9` or `10` roots and the stationary-unit-safe set in
(34) has `11` or `12`.  Their intersection therefore contains at least
seven roots.  Choose two distinct roots `t,t'` in that intersection.  By
openness away from finitely many endpoints, there is a small rational
interval `I` about `y_0` for which both selected branches remain guard-safe
and `q_k`-safe.

Put

```text
B_t={(y+t)/13:y in I},
R_t=1_(B_t),                                                  (69)
```

and define `R_(t')` similarly.  These are honest rational Boolean BV step
functions satisfying the same common support factors.

Write a point of the circle uniquely, away from endpoints, as

```text
x=(y+r)/13,                 y in [0,1), r in F_13.             (70)
```

Directly from (3),

```text
E_aR_t((y+r)/13)
 =1/13 1_I(y)zeta_13^(-a(t-r)).                              (71)
```

Therefore every one of the thirteen projections is nonzero and

```text
|E_aR_t(x)|^2=1/169 1_I(13x),                                (72)
```

independent of both `a` and `t`.  For every common delayed owner
`g_M(x)=G(13^M x)`, `M>=1`,

```text
integral_T g_M|E_aR_t|^2
 =1/169 integral_I G(13^(M-1)y)dy,                           (73)
```

again independent of `a` and `t`.  Thus `R_t` and `R_(t')` have different
root addresses but identical complete owner-weighted diagonal energy
banks.  Their charged projections differ by the visible phase

```text
E_aR_(t')/E_aR_t=zeta_13^(-a(t'-t))                          (74)
```

on corresponding branches.

This is a BV realization, not a formal delta-vector.  It is not claimed to
be a full lawful THM-2449 anchored table.  It proves that the common guard,
the common stationary unit, and the complete diagonal energy bank cannot by
themselves invert root address.  Any such inversion must retain charged
phase or use additional table structure.

## 8. Exact connections and stopping boundary

The theorem sharpens several nearby mechanisms while preserving their type
boundaries.

### THM-2527--2531: Boolean consumers versus uncontracted phase

THM-2527 couples the late owner to a positive Boolean full-mask layer;
THM-2528 opens the four-arm signed path and scalarizes the 72/5,184 norm
bank; THM-2529 and THM-2530 give deep-comb and anchored-Gram consumers; and
THM-2531 selects an actual occupied-to-empty predecessor boundary.  Those
results make the phrase “no Boolean root consumer” obsolete.

The new coordinate here is different: equations (5), (29)--(30), and
(68e)--(68k') retain the **uncontracted charged phase** of the actual
owner-marked event, while (62)--(68i) put three mixed gains on literal
slope-labelled predecessor/Radon currents.  None of the existing Boolean
consumers is proved to intertwine that direct phase with a semantic target
or arrival character.  Conversely, this theorem does not replace their
positive Boolean factorizations or THM-2531's canonical boundary selector.

### THM-2524: autocorrelation versus predecessor phase

THM-2524's translated `chi_7` bank reconstructs the entire centred
collision autocorrelation.  The hostile (69)--(73) has identical diagonal
autocorrelation data at different addresses.  The present theorem does not
invert THM-2524.  It uses the retained event `F` *before* autocorrelation
and the module identity (5) to keep one charged phase ladder.

### THM-2312: cubic closure versus charged phase

The sparse-root bispectrum is invariant under a common root translation;
on a one-root fibre it is pure mass cubed.  Equation (74) is deliberately
not invariant.  It retains ordinary root phase, but it is neither a
gauge-invariant cubic terminal current nor a terminal-component phase.

### THM-2334: ordinary frequency versus relation address

The integer `n` in (29) or (65) is one physical Fourier frequency.
THM-2334's address is a relation-lattice orbit of many endpoint
decompositions.  A nonzero ordinary coefficient does not prove that a
prescribed relation-address orbit sum is nonzero.  The root-character and
target-twist actions remain different representations.  The pair
`(X_kappa,Y_kappa)` in (30c) does retain THM-2349's marked `c_3` edge, but
it does not select a prescribed THM-2334 target-twisted relation-orbit term.

### THM-2457 and THM-2474: different positivity inputs

THM-2457 starts with two disjoint semantic packets in one supplied oriented
chart and may coarsen complete atoms to obtain positive co-support.
THM-2474 starts at a squarefree collision corner and obtains every primitive
colour for two nonnegative packets.  Here the eventwise line of (2) has one
literal Boolean event and all twelve prime colours; the table line has no
disjoint-packet corner and therefore guarantees only three gain ratios.
The endpoint-Prony algebra is shared, but the positive input is different.

### Radon slope versus marked gain corolla

At the level of indices, the equation

```text
a=b tau                                                        (75)
```

identifies `lambda=a/b` with the candidate THM-2507/2508 toothpick slope
`tau`.  The phase lift (68e)--(68i) now realizes that label as the actual
THM-2521 aligned-Radon difference on a common predecessor profile: each
proved gain incidence supplies a nonzero `C_lambda p`.  What is still
unproved is an equality with the THM-2508 Radon transform of the lawful
THM-2449 ANOVA defect.  That representation additionally remembers an
ordered cut, carry cocycle, and target semantics.  Three slopes do not
reconstruct that general THM-2449/2508 strip defect; by contrast, any one
oriented `C_lambda` already reconstructs the newly built centred predecessor
profile by THM-2532.

The gain in THM-2315 is a different object: a primal target--target point
`[1:g]` on a marked projective line.  Its target swap acts by

```text
g -> g^(-1),
```

with fixed gains `+-1`.  Collision converse acts on (75) by

```text
lambda -> -lambda,
```

which has no fixed nonzero gain.  The involutions have different cycle
types on `F_13^*`, so no direct equivariant identification exists without a
new target/collision intertwiner.

### Final type ledger

The proved maps and their first losses are:

| source | map | exact target | first information not supplied |
|---|---|---|---|
| THM-2522 event collision energy | invariant-multiplier identity (5) | twelve phase ladders on one Boolean event | terminal component and relation address |
| the same late-owner event inside its THM-2349 carrier | abstract carrier triangle (30a)--(30f) | twelve owner-marked 91-unit deepest-comb incidences | target character, terminal component, unoriented quotient edge |
| mixed lawful table channel | guard Vandermonde + Galois + phase-lifted Cayley derivative | at least three gains and 216 phase ladders; three slope-labelled nonzero aligned-Radon differences per channel; flat energy gives 864, nonflat gives a separate odd energy contrast | Boolean target/terminal intertwiner, THM-2508 cut/carry identification, prescribed gain, canonical atom |
| eventwise Boolean predecessor | root-fibre map (68e)--(68g) + THM-2521 | rational nonuniform owner/terminal-word-attached `K_14` potential with all 72 pointwise mixed and 5,184 Hilbert cut modes | noncancelling scalar twist and a Boolean output proved to intertwine this direct phase with semantic target/arrival data |
| one phase ladder | grouped endpoint current (15) | bounded positive integer lift | common lift across colours or gains |
| complete diagonal energy bank | squared magnitude | collision power data | root address, by (69)--(74) |

In particular, this theorem proves none of the following:

1. a common gauge quotient `h` for the twelve event colours or the `216`
   mixed channels;
2. a terminal component, exact relation-lattice address, or prescribed
   THM-2457 complete atom;
3. a canonical target/chronological meaning for `tau`, an identification of
   a mixed gain `lambda_i` with the live guard sign, a Boolean oriented
   cospan, or an owner-loop current;
4. an identification of collision gain with the THM-2315 target corolla;
5. an old-deep rebase at the future owner;
6. a scalar row exclusion or LRC(14).

## 9. Exact referee and independent audit

Run

```bash
python3 04-computation/lrc14_owner_weighted_collision_phase_thm2533.py
python3 -O 04-computation/lrc14_owner_weighted_collision_phase_thm2533.py
```

Both executions must byte-match

```text
05-knowledge/results/lrc14_owner_weighted_collision_phase_thm2533.out
```

The dependency-free referee uses exact rational, cyclotomic, and finite-
field arithmetic.  It verifies `1,404` projector-module identities, `104`
Boolean-owner identities, all twelve eventwise colours, `882` grouped-jump
Prony minors and `240` frequency bounds, all twelve inherited marked
`91`-unit deepest-comb residues, `14,196` consecutive guard Fourier minors,
`9,295` four-mode stalk ranks, all `220` gain triples, the sharp three-gain
cyclotomic fibre, and `1,872` entries of its nonnegative rational table
completion.  The last completion is explicitly a support-sharp hostile,
not a lawful anchored THM-2449 table.

An independent line audit rederived the projector and owner identities,
the grouped-jump quantifiers, the lawful-clock placement, the repaired
Galois nonzero/nonnegativity step, and the sharp hostile boundary.  A second
audit checked the phase-lifted predecessor normalization (68e)--(68g), all `216`
slope--channel incidences, the pointwise `72` and Hilbert-valued `5,184`
mode consequences, and the zero-frequency boundary (68k)--(68k').  Normal
and optimized transcripts byte-matched the stored output and recorded
hashes.

The advance is exact but narrower: at the uniform physical first collision,
positive owner-weighted energy is no longer trapped in a magnitude-only
pair cospan.  It has a bounded, charged, owner-attached ordinary Fourier
phase lift; the same Boolean event feeds twelve owner-marked deepest-comb
triangles and a full pointwise `K_14` potential spectrum, while the mixed
table feeds three Galois-saturated gain slices whose phase lifts are genuine
slope-labelled Radon-difference currents.  Their untwisted scalar averages
vanish, locating the next obligation at a target/terminal twist rather than
at construction of a common predecessor profile.
**QED.**
