---
id: THM-2517
title: "Cubic anchored spectrum, FLT(3), and marked future Boolean lift"
status: >
  PROVED + VERIFIED-EXACT + LEAN-VERIFIED ARITHMETIC KERNEL +
  INDEPENDENTLY AUDITED.  For every nonnegative rational C_7-by-C_13 table
  with positive delta anchor and nonflat owner row, the entrywise cube has
  all 72 primitive mixed coefficients nonzero.  One nonflat target column
  already supplies six nonzero cubic rectangles: a zero would be a rational
  solution x^3=y^3+a^3, forbidden by FLT(3) except at (x,y)=(a,0).
  Hence the cube has all 5,184 affine-cut coefficients and is the least
  universal single pure-power moment; degrees one and two fail sharply on
  replicas and rational Pythagorean columns.  THM-2516 now unconditionally
  yields an owner-supported three-arm rational simplex.  A weighted one-copy
  lift, a three-time Boolean same-circle lift, and a four-time lift with one
  common future owner/word factor preserve the old deep diagonal.  A sharp
  zero-or-norm dichotomy for the seven lawful phase shifts gives either a
  direct covariant null-row detector or, in the full-support branch, a
  genuinely source-neutral common Boolean gate via the cyclic K_(7,7)
  one-factorization.  Every norm-gate cell contains the actual phase-zero
  owner once.  Its owner epoch is a seven-point torsor, not one fixed clock,
  so owner-loop emission, row exclusion, and LRC(14) remain open.
source: codex-2026-07-27-anchored-cubic-spectrum
depends_on:
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2466-delayed-word-simultaneous-drift-service-retention
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
  - THM-2512-lawful-interaction-cut-bundle-transplant-and-replica-dichotomy
  - THM-2516-anchored-moment-simplex-disintegration-and-owner-star-recovery
related:
  - THM-2456-two-root-replica-uniform-offset-boundary
  - THM-2460-idempotent-semantic-word-copy-and-word-block-cosupport-boundary
  - THM-2513-anchored-first-or-second-moment-spectrum-and-pair-space-boundary
  - THM-2515-haar-self-correlation-disintegration-and-rational-shift-recovery
script: 04-computation/lrc14_anchored_cubic_spectrum_thm2517.py
output: 05-knowledge/results/lrc14_anchored_cubic_spectrum_thm2517.out
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCAnchoredCube.lean
lean_root: 04-computation/lean/TournamentH7/TournamentH7.lean
script_sha256: c994979188df7812937787edc8aabc1f9735f48801df4e860bdabf89b1fdb5d7
output_sha256: ceab24add134f3048f32c16a634007dc9fa3e2ed10bac2b2c3c73a1b53222826
lean_sha256: 904c8e3153599be1fbf750ab6534638a75ce0f91b40848956624bcaba619fec0
hash_basis: working-tree bytes (LF)
---

# THM-2517 -- the cubic moment is universally charged

**PROVED + VERIFIED-EXACT + LEAN-VERIFIED ARITHMETIC KERNEL +
INDEPENDENTLY AUDITED.**

THM-2513 showed that the pair `(A,A o A)` cannot vanish spectrally, but its
two channels appeared essential: replicas kill the first and rational
Pythagorean columns can kill the second.  The cube changes the arithmetic
geometry.  A zero anchored cubic rectangle would be a rational point on the
Fermat cubic.  Exponent three is the first pure moment at which that locus
has no nontrivial nonnegative rational point.

The consequence is branch-free and scalar:

```text
lawful nonnegative rational anchored table A
  -> one table C=A o A o A
  -> all 72 mixed colours
  -> all 5,184 primitive affine-cut modes.                     (1)
```

The same cubic signal has three complementary realizations: a one-copy
rational weighting, an exact finite owner-supported simplex, and delayed
Boolean products on one circle.  The residual is no longer moment
cancellation or pair-space descent.  It is the source typing of the marked
future owner.

## 1. Local Fermat-cubic separation

Let

```text
A=(A_(ell,s))_(ell in F_7,s in F_13) in Q_(>=0)^(7 x 13)     (2)
```

obey

```text
A_(ell,0)=a 1_(ell=0),                    a>0,                 (3)
```

and assume `s -> A_(0,s)` is nonconstant.  Put

```text
C=A^(o 3),                         C_(ell,s)=A_(ell,s)^3.      (4)
```

Choose one target `s_*` with

```text
A_(0,s_*)!=a.                                                 (5)
```

For `ell!=0`, write

```text
x=A_(0,s_*),                         y=A_(ell,s_*).
```

The anchored rectangle of `C` is

```text
Gamma_(ell,s_*)=x^3-y^3-a^3.                                 (6)
```

The exact equality boundary is

```text
Gamma_(ell,s)=0                    iff (x,y)=(a,0).            (7)
```

For the forward implication, rewrite (6) as

```text
y^3+a^3=x^3.                                                 (8)
```

If `y=0`, injectivity of the odd cube gives `x=a`.  If `y!=0`, then
`a,y>0`; equation (8) also forces `x>0`.  Fermat's theorem for exponent
three over `Q` forbids (8).  The reverse implication is immediate.

Equation (5) and (7) show that

```text
Gamma_(ell,s_*)!=0                   for every ell!=0.        (9)
```

Thus one target column supplies all six local witnesses before ANOVA,
Fourier transform, Galois propagation, or cut construction.  In particular,

```text
I(C)!=0.                                                        (10)
```

There is also an exact global classification:

```text
I(C)=0

iff A_(0,s)=a and A_(ell,s)=0 for every s and every ell!=0.   (11)
```

Indeed, zero interaction makes every anchored rectangle zero, so (7)
applies column by column and source leg by source leg.  Conversely the pure
baseline in (11) has no interaction.  Owner-row nonflatness excludes exactly
this boundary.

## 2. Why the jump from two to three is structural

The first three pure moment loci have different geometry.

### Degree one: an affine replica line

The first moment vanishes precisely on

```text
x-y=a.                                                        (12)
```

Thus every nonconstant nonnegative replica profile of THM-2449 is a sharp
degree-one hostile.  On a nontrivial replica column, write `x=a+w`, `y=w`;
the square rectangle has the oriented sign

```text
x^2-y^2-a^2=2aw>0.                                           (13)
```

### Degree two: a rational Pythagorean hyperbola

The entire nonnegative rational `I(A^(o 2))=0` branch consists of columns

```text
x_s^2-y_s^2=a^2,                                             (14)
```

with the same `y_s` on all six nonowner rows.  Every nontrivial rational
column is parametrized by

```text
t_s=y_s/(x_s+a) in Q intersect [0,1),

x_s=a(1+t_s^2)/(1-t_s^2),
y_s=2a t_s/(1-t_s^2).                                       (15)
```

Its first rectangle has the opposite oriented sign

```text
x_s-y_s-a=-2a t_s/(1+t_s)<0                 when t_s>0.       (16)
```

The smallest familiar column is `(a,x,y)=(3,5,4)`.  Taking the anchor
column `(3,0,...,0)` and every other column `(5,4,...,4)` gives a
nonnegative rational anchored nonflat table whose square is exactly anchor
`9` plus replicas `16`.  Thus a single quadratic moment cannot suffice.

### Degree three: the Fermat curve

The cubic zero locus is

```text
x^3-y^3=a^3.                                                 (17)
```

Unlike the genus-zero conic (14), the Fermat cubic has no nontrivial positive
rational points.  Therefore degree three is the least universal **single
pure-power moment**.  This does not supersede THM-2513's lower-degree vector
observer `(A,A^(o 2))` or its ternary scalar tag.

The hypotheses are sharp.

- Without rationality, take `a=1`, `y=1`, `x=cuberoot(2)` in every
  nonanchor column; the real nonnegative cube is a replica.
- Without nonnegativity, take `a=1`, owner entry `x=0`, and all six
  nonowner entries `y=-1` away from the anchor column; the signed rational
  cube is a replica.
- Without `a>0`, take the zero anchor column and all-one remaining columns.
- Without nonflatness, (11) is the exact both-zero boundary.

## 3. Complete spectrum, cut bundle, and a denominator invoice

Every entry of `C` is rational.  The coprime Galois group

```text
Gal(Q(zeta_91)/Q)
 =(Z/7Z)^* x (Z/13Z)^*                                      (18)
```

acts transitively on the `72` mixed characters.  Equation (10) and
THM-2449's ANOVA/Parseval dictionary therefore give

```text
Chat(kappa,b)!=0                    for all kappa,b!=0.        (19)
```

Transpose `I(C)` to the row-zero defect `d_C(s,ell)=I(C)_(ell,s)`.
THM-2508/2512 gives

```text
Psi^C_(tau,a_0)(alpha,beta)
 =91 K_(alpha tau,beta)Chat(beta a_0,-alpha)!=0              (20)
```

for `tau,alpha in F_13^*` and `a_0,beta in F_7^*`.  Hence all

```text
12*12*6*6=5,184                                               (21)
```

primitive cut coefficients survive.  The same row-zero argument gives at
least `294` nonzero toothpick components and at least `3,528` nonzero raw
root-colour entries.

There is a rational quantitative floor.  If one integer `D>=1` clears every
entry of `A`, then (9) gives

```text
Gamma_(ell,s_*) in D^(-3)Z nonzero,
|Gamma_(ell,s_*)|>=D^(-3).                                  (22)
```

An anchored rectangle is a four-entry functional of `I(C)` of Euclidean
norm `2`.  Thus

```text
||I(C)||_2>=1/(2D^3),

sum_(kappa,b!=0)|Chat(kappa,b)|^2
 =||I(C)||_2^2/91 >=1/(364D^6).                              (23)
```

No individual conjugate magnitude floor follows without a cyclotomic height
bound.

## 4. Exact weighted one-copy deep lift

For the lawful THM-2449 table, write

```text
H_ell(r,s,t)=integral_T h_(ell,r,s,t)(x) dx,
A_(ell,s)=sum_r H_ell(r,s,0),                                (24)
```

where every `h` is a nonnegative rational Boolean step product and

```text
h_(ell,t,s,t)=0 almost everywhere.                           (25)
```

Define

```text
H^[3]_ell(r,s,t)=A_(ell,s)^2 H_ell(r,s,t).                   (26)
```

Then

```text
sum_r H^[3]_ell(r,s,0)=A_(ell,s)^3=C_(ell,s),
H^[3]_ell(t,s,t)=0.                                          (27)
```

The THM-2449 deep transform applies verbatim:

```text
J^[3](kappa,0,b)=Chat(kappa,b)/13,
sum_alpha J^[3](kappa,alpha,b)=0.                            (28)
```

For each of the `72` mixed colours, (19) therefore forces some
`alpha!=0` and `tau` with a nonzero old-sheet deep coefficient.  This is an
exact one-copy nonnegative rational weighting.  The scalar
`A_(ell,s)^2` is global data, not a local Boolean ancestry factor.

## 5. Exact owner-supported simplex: THM-2516 becomes unconditional

Apply THM-2516 at moment order `m=3`.  Its unanchored identity disintegrates
the cube into exact three-point/two-difference fibres:

```text
A_j^3
 =average_(u,v) integral_T F_j(x)F_j(x+u)F_j(x+v) dx.        (29)
```

Because (19) supplies a nonzero cube coefficient unconditionally, one
rational fibre carries that coefficient with no magnitude loss, and its
rational table has all `72` mixed and all `5,184` cut coefficients.

More strongly, let `G_owner` be any positive marked owner density of mean
`rho>0`.  THM-2516's anchored identity is

```text
rho A_j^3
 =average_(u_1,u_2,u_3)
   integral_T G_owner(x) product_(q=1)^3 F_j(x+u_q) dx.       (30)
```

Thus one rational **three-arm/four-point** simplex fibre has the complete
spectrum, every entry contains the same positive owner factor, and the
selected coefficient has modulus at least `rho` times the cube coefficient.
This closes the moment-existence premise left conditional in THM-2516.

The shifts in (29)--(30) are exact rational difference coordinates.  They
need not lie on the lawful phase-bank/word grid; owner support is present,
but ancestry lawfulness is not automatic.

## 6. Three-time Boolean lift on one circle

There is a canonical lawful-scale alternative to arbitrary translations.
For one response density `F=F_(ell,s)`, put

```text
A=integral F,               M=||F||_infinity,
V=Var(F),                   N=13^L,

Q^L_(ell,s)
 =integral_T F(x)F(Nx)F(N^2x) dx.                            (31)
```

Two applications of the exact BV covariance estimate give

```text
|Q^L_(ell,s)-A_(ell,s)^3|
 <=V^2/(12N) [A+M(1+1/N)]
 <=M V^2/(4N).                                                (32)
```

Indeed, first

```text
|integral F(x)F(Nx)dx-A^2|<=V^2/(12N),                       (33)
```

then use

```text
Var(F(x)F(Nx))<=M V(N+1)                                    (34)
```

against the last factor at frequency `N^2`.  The convergence is simultaneous
over the finite `7 x 13` bank.  Every table `Q^L` is nonnegative and rational.
For all sufficiently large `L`, its interaction is nonzero; rational Galois
then gives all `72` mixed and all `5,184` cut modes exactly.  If `D` is as in
(22) and `epsilon_L` bounds the maximum entry error in (32), the explicit
condition

```text
4 epsilon_L < D^(-3)                                        (35)
```

already preserves the selected cubic rectangle.

The literal Boolean refinement is

```text
mathcal H^[L,3]_ell(r_0,r_1,r_2,s,t)
 =integral_T
   h_(ell,r_0,s,t)(x)
   h_(ell,r_1,s,0)(N x)
   h_(ell,r_2,s,0)(N^2 x) dx.                                (36)
```

Every entry is a Boolean three-epoch intersection before sums or DFT.
Summing `r_0,r_1,r_2` at `t=0` gives (31), while

```text
mathcal H^[L,3]_ell(t,r_1,r_2,s,t)=0                         (37)
```

retains the old deep diagonal.  After summing the two future roots, (36)
converges entrywise to `H^[3]` in (26).  More explicitly, if
`W=Var(h_(ell,r_0,s,t))`, then the error is at most

```text
[V^2+A W V]/(12N) + M W V/(12N^2).                           (38)
```

Hence one common sufficiently large `L` retains a chosen old-sheet deep
coefficient for every mixed colour.  Taking `L` even makes
`N=N^2=1 mod 7`, so all three printed source labels have the same septimal
orientation.  The two future roots remain distinct ancestry sheets.

## 7. Adding one actual future owner/word factor

Let

```text
G(z)=E(z)Q(T^K z),                       rho=integral G>0      (39)
```

be a supplied positive Boolean future owner-to-word block of THM-2478, and
write `W_G=Var(G)`.  Delay it beyond all three response epochs:

```text
K^L_(ell,s)
 =integral_T F(x)F(Nx)F(N^2x)G(N^3x) dx.                    (40)
```

The same proof, now using

```text
Var(F(x)F(Nx)F(N^2x))<=M^2 V(1+N+N^2),                      (41)
```

gives the explicit invoice

```text
|K^L_(ell,s)-rho A_(ell,s)^3|
 <=1/(12N) {
      rho V^2[A+M(1+1/N)]
      +M^2 V W_G(1+1/N+1/N^2)
   }.                                                        (42)
```

Thus for every sufficiently large even `L`, `K^L` is a nonnegative rational
table with positive delta anchor and all `72`/`5,184` primitive modes.  Its
refined integrand is (36) multiplied by the same literal Boolean factor
`G(N^3x)`.  It still has (37), converges to `rho H^[3]`, and places every
entry on one common actual future owner/word support.  The delayed factor is
root-constant and factor-wise target-neutral relative to each earlier local
root chart once `L` exceeds the finite packet horizons; no rebase is used.

One type boundary is load-bearing.  A fixed phase-zero `G` is **not
septimally source-neutral**: powers of `13` are `+-1`, never `0`, modulo
seven.  Therefore the mixed `kappa` in (40) is the old algebraic table label,
not automatically the total physical source residue of the four-factor atom.
Adding a redundant translated danger detector is the fake second harmonic
forbidden by THM-2460: idempotence collapses it before DFT, and shifting only
that duplicate omits the mandatory word skew.  Shifting the entire future
owner/word bank is lawful but makes its mean phase-dependent and destroys the
simple `rho A^3` limit.

There is an exact tradeoff.  The Boolean union over the full lawfully shifted
orbit of `G` is source-neutral, but its semantic mark becomes "some shifted
owner/word phase," not the actual phase-zero owner.  Section 8 repairs even
that loss by scheduling the whole orbit at distinct epochs.

## 8. Zero-or-norm phase dichotomy and the Latin owner scheduler

The tradeoff at the end of Section 7 can be repaired without a fake harmonic.
Retain the full lawfully shifted owner/word blocks

```text
G_gamma,                         q_gamma=integral G_gamma>=0,

q_0>0,                           gamma in F_7,                 (O1)
```

where shifting includes the mandatory terminal `R gamma` word skew of
THM-2442/2460.  There are two exact branches.

### A zero phase creates a covariant rectangle

Suppose `q_lambda=0`.  Since `G_lambda` is Boolean and nonnegative, it is
zero almost everywhere.  Gate row `ell` of the three-response product by the
matching full block `G_ell` at a later even epoch.  This row family is
source-covariant.  Row `lambda` is exactly zero, while the owner row and
anchor converge to

```text
q_0 x^3,                              q_0 a^3.                 (O2)
```

At the nonflat target `s_*`, the owner-versus-zero rectangle tends

```text
q_0(x^3-a^3)!=0,                                             (O3)
```

by nonnegative odd-cube injectivity.  Hence all `72`/`5,184` modes survive
for every sufficiently large even delay.  A missing phase is itself a
source-covariant detector.

### No zero phase gives a cyclic norm

Now suppose every `q_gamma>0`, and put

```text
P=product_(gamma in F_7) q_gamma>0.                           (O4)
```

Choose seven separated future epochs

```text
D_d=(3+rep(d))L,                       d in F_7,               (O5)
```

with `L` even and beyond every finite packet horizon.  For each
`c in F_7`, define the assignment cell

```text
P_c^L(x)=product_(d in F_7) G_(c+d)(T^(D_d)x).                (O6)
```

Each cell uses every source phase exactly once.  Distinct cells are disjoint:
at every fixed epoch `d` they use distinct mutually exclusive source-danger
phases.  Therefore

```text
U^L=union_(c in F_7) P_c^L=sum_c P_c^L                       (O7)
```

is one Boolean gate.  Source translation permutes the seven `c`-cells, so
`U^L` is exactly source-neutral.  Higher-order mixing gives

```text
integral U^L -> 7P>0.                                        (O8)
```

Gate every row of the response cube by this same `U^L`.  The resulting
nonnegative rational table converges to

```text
7P A^(o 3).                                                   (O9)
```

Here is the needed higher-order mixing argument.  Approximate the finite
family of bounded rational step functions in `L^1` by trigonometric
polynomials of a common bandwidth `M`; the product error is controlled by
the usual telescoping sum because all other factors are uniformly bounded.
For `N=13^L>M+1`, a Fourier term in any one of the seven-factor products
in (O6), or in a ten-factor response-times-cell product, can contribute to
the integral only if

```text
sum_(j=0)^9 n_j N^j=0,                         |n_j|<=M.      (O9a)
```

If `J` is the largest index with `n_J!=0`, then

```text
N^J <= |n_J|N^J
    <= M sum_(j<J)N^j
    < N^J,                                                     (O9b)
```

a contradiction.  Thus only the constant Fourier term survives.  Passing
through the `L^1` approximations proves both (O8) and (O9), uniformly over
the finite phase/row/target bank.  This also records why separated epochs,
not a same-time product, are load-bearing.

Equations (9), (19), and rational Galois propagation now give all primitive
table and cut modes at every sufficiently large even `L`.  Fully root-refine
all three response sums as in (36): the first carries the old `(r,t)` chart
and the two future roots remain sidecars to be summed afterward.  This
retains the old deep diagonal, and each resulting ten-factor cell is a
literal Boolean intersection before DFT.  There are three response packet
slots and seven fully shifted composite handoff blocks.  Expanding the seven
`G_gamma=E_gamma Q_gamma o T^K` blocks adds seven terminal-word times, while
the response packets retain their pre-existing internal clocks.

This is the correct combinatorial carrier.  The table

```text
phase gamma = c+d                    versus future slot d     (O10)
```

is the cyclic one-factorization of `K_(7,7)`: the seven `c`-rows are its
parallel perfect matchings.  The intrinsic observable is assignment, not
orientation, so a bipartite factor chart is more faithful than a forced
tournament.

Every cell `P_c^L` contains the actual phase-zero owner/word block exactly
once, at slot

```text
d=-c.                                                         (O11)
```

Thus `U^L` guarantees an honest phase-zero future owner, not merely some
shifted owner.  But its epoch is the free `C_7` torsor coordinate `c`.
Summing all `c` restores source neutrality and forgets the clock; fixing one
`c` fixes the clock and breaks source neutrality.

This seven-slot cost is minimal at this interface.  If a seed conjunction
uses a phase set `S subset F_7` and every source translate `S+a` must contain
phase zero, then `-a in S` for every `a`, hence `S=F_7`.  Mutually exclusive
phases cannot share an epoch, so at least seven future slots are necessary.
A source-neutral deterministic choice of one common owner slot would be a
fixed point of the free `C_7` assignment torsor and is impossible.  Repairing
it requires exactly the missing anti-diagonal ancestry cocycle/sheet, not
another colour census.

The zero-or-norm law is therefore sharp:

```text
some q_gamma=0  -> the null covariant row creates a rectangle;

all q_gamma>0   -> the cyclic norm product gives one neutral Boolean gate.
                                                                  (O12)
```

What remains is branch-specific.  In the zero branch, row zero already has
an actual phase-zero owner at one fixed common epoch, but no
emission/source-arrival intertwiner.  In the full-support norm branch, the
actual owner lies at a free `C_7`-torsor epoch and needs either a gauge-safe
section or a clock-covariant intertwiner.  Moment cancellation, source
covariance, Booleanity, owner presence, and deep-sheet survival are no longer
the missing coordinates.

## 9. Why same-time tournament stars cannot replace the delays

At one time, the seven source-phase indicators `e_u` are mutually exclusive:

```text
e_u^2=e_u,                       e_u e_v=0 for u!=v.           (43)
```

Their Boolean algebra has only the eight atoms

```text
e_u,                    1-sum_u e_u.                          (44)
```

In the actual exhaustive seven-phase danger partition, the last atom is
null; allowing a "none" atom only strengthens the general no-go.

Every same-time Boolean polynomial therefore reduces to a constant plus
one-body terms.  More invariantly, any phase-translation-equivariant Boolean
output family `K_v` is determined by one cyclic kernel `phi`.  If `p_u`
denotes the physical phase-mass bank (with its target coordinate retained),
then

```text
K_v(x)=phi(v-L(x))                                      when L(x) exists,

khat(kappa,b)=phihat(kappa) phat(kappa,b)              for kappa!=0.  (45)
```

It cannot revive a vanished physical mixed source character.

A fixed-center star can appear to do so only by breaking the gauge.  On
uniform `C_7`, let `e_u=1_(x=u)` and put

```text
K_0=e_0,                       K_ell=e_0 union e_ell.          (46)
```

The masses are `(1/7,2/7,...,2/7)`, so every naive nontrivial **raw,
unnormalized** leaf Fourier sum is `-1/7` (normalized value `-1/49`), although
all `42` distinct-phase intersections are empty.  The signal is center
main-effect leakage; on the full ordered-pair
carrier, two-way ANOVA removes it.  A second epoch is structurally necessary
to make genuine edges `e_u(x)e_v(Nx)` nonempty.  The delayed lifts above are
not decorative tournament encodings: they create the missing adjacency
coordinate.

## 10. Exact companion and Lean kernel

Run

```bash
python3 04-computation/lrc14_anchored_cubic_spectrum_thm2517.py
python3 -O 04-computation/lrc14_anchored_cubic_spectrum_thm2517.py

cd 04-computation/lean/TournamentH7
lake env lean TournamentH7/LRCAnchoredCube.lean
lake build TournamentH7.LRCAnchoredCube
```

The exact companion:

- checks `7,500` nonnegative integer and `6,734` bounded-denominator rational
  cubic equality/floor controls;
- checks `60,000` one-target/six-source cubic rectangles;
- verifies the degree-one replica, degree-two `3-4-5`, signed cubic, and
  rational Pythagorean hostiles;
- checks `58` exact three-time initial-arc BV invoices and `182` delayed table
  entries at both septimal parities; and
- reproduces the one-hot `42`-empty-edge / raw `-1/7` marked-star gauge
  hostile, the `7 x 7` Latin scheduler, its unique invariant phase set, and
  all `64` zero-or-norm mean profiles.

The Lean module imports Mathlib's proved
`fermatLastTheoremThree : FermatLastTheoremFor 3`, transports it to `Q` via
`fermatLastTheoremFor_iff_rat`, and proves:

```text
cube_rectangle_eq_zero_iff,
cube_rectangles_nonzero_at_nonflat_column,
exists_column_with_all_cube_rectangles,
all_cube_rectangles_eq_zero_iff_pure_baseline.                (47)
```

The six-leg statements explicitly require `ell!=0`; the owner leg's full
four-cell rectangle is zero and is not miscounted.  The module is imported by
the project root, contains no `sorry`, and kernel-checks the global boundary
(11) as well as the local witness.  The simplified expressions in the Lean
boundary become full four-cell rectangles using hypothesis (3).  Its axiom
audit has only Mathlib's standard logical axioms.

Independent audits rederived the FLT equality boundary, the Pythagorean
degree-two classification, every BV constant in (32), (38), and (42), the
old deep-diagonal inheritance, even-parity source gauge, the fixed-`G`
source-neutrality defect, the zero-or-norm Latin scheduler and its minimality,
and the one-hot equivariance no-go.  No live scalar row is removed.  LRC(14)
remains open, now at the sharply typed question:

```text
zero branch: intertwine its fixed owner clock with emission/source-arrival;

full-support branch: either select one clock from the neutral C_7 torsor
without breaking gauge, or build a clock-covariant intertwiner;

otherwise prove the live rows forbid the corresponding structure.           (48)
```

**QED.**
