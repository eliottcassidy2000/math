---
id: THM-2478
title: "Delayed owner-handoff graft and deep-sheet rebase boundary"
status: >
  CLAIMED + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING. A fixed
  positive owner-to-word block may be delayed and inserted into any
  already nonzero endpoint-drift/common-root-service packet. On the old
  oriented root base it is one Boolean root-constant, target-neutral
  filter, so the old drift tends to the handoff mass times the old drift
  and the service tends to the same positive fraction, with explicit BV
  error O(13^(-L)). Thus no prior owner support is needed to install an
  honest future owner-to-word coupling. More strongly, the twelve complex
  densities of a fixed THM-2471 collision stalk may be delayed on their
  physical collision base; the same old nonzero target/deep drift survives
  in every future collision colour simultaneously. All products exist on
  the Boolean ancestry stalk before the finite root DFT. This is a lawful
  future-owner graft on the old endpoint/root base, not a rebase theorem:
  after transfer to the future owner, a past deep probe of 13-adic depth
  lambda becomes ancestry-sheet dependent for L>lambda. Source and arrival
  atom representations also remain inequivalent, and aggregate H-drift
  still need not force owner-loop drift. No temporal root identification,
  scalar-row exclusion, or LRC(14) closure is claimed.
source: temporal-cocycle-probe-2026-07-27
depends_on:
  - THM-2306-owner-normalized-disjoint-support-and-first-collision-shell
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2459-four-atom-drift-and-root-service-coarsening
  - THM-2466-delayed-word-simultaneous-drift-service-retention
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
related:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2445-twenty-four-cell-graft-owner-conditioning
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2474-squarefree-first-collision-primitive-character-saturation
script: 04-computation/lrc14_delayed_owner_handoff_graft_thm2478.py
output: 05-knowledge/results/lrc14_delayed_owner_handoff_graft_thm2478.out
script_sha256: a611026c7920247b8c5949b38fd1906664d735972b23ff1f330c942f49ec2228
output_sha256: b722d88df0c819ee066949e039ce8cce03d7fae4211cac80f5eed12e7fb17889
hash_basis: working-tree bytes (LF)
---

# THM-2478 -- delay the whole owner handoff, not only its terminal word

**CLAIMED + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-2466 delays one positive terminal word after fixing endpoint drift and
root service.  Its semantic application assumed that both fixed packets were
already supported on the corresponding owner.  THM-2471 shows why that
assumption is not automatic: source, arrival, and deep endpoint atoms are
different temporal copies.

The order can instead be reversed.  First form the entire positive handoff

```text
G(z)=1_(E_j)(z) 1_(Q)(T^K z),
```

and only then delay `G`.  This installs both ends of the semantic edge at
once.  On the old oriented root fibre it is one common Boolean base filter:

```text
old endpoint/root packet at x
  + G(T^L x)
  = old packet at x
    + future owner E_j at T^L x
    + its terminal word Q at T^(L+K)x.                 (1)
```

The construction is deliberately forward-looking.  It preserves the old
deep/root table and adds a future owner edge; it does not move the old deep
probe to the new owner time.

## 1. Delayed positive-handoff lemma

Put

```text
T(x)=13x mod 1.
```

Let `E,Q` be rational Boolean circle-step functions, let `K>=1`, and assume

```text
G(z)=E(z)Q(T^Kz),

rho=integral_T G>0.                                  (2)
```

Let `H` be a finite-dimensional real or complex Hilbert space.  Fix on one
circle base

```text
b:T->H,                         a Hilbert-valued BV function,

s:T->R_(>=0),                  a scalar BV function,

q=integral_T b,                D=||q||^2>0,

M=integral_T s>0.                                      (3)
```

For `L>=1`, put

```text
N_L=13^(L-1),

q_L=integral_T G(N_L y)b(y)dy,

M_L=integral_T G(N_L y)s(y)dy.                        (4)
```

Then

```text
q_L -> rho q,

M_L -> rho M.                                         (5)
```

Write

```text
C_G=min(rho(1-rho),Var(G)/12).
```

The exact BV invoices are

```text
||q_L-rho q||
 <=C_G Var_H(b)/N_L,

|M_L-rho M|
 <=C_G Var(s)/N_L.                                    (6)
```

Consequently, if

```text
N_L>=2 C_G Var_H(b)/(rho sqrt(D)),

N_L>=2 C_G Var(s)/(rho M),                            (7)
```

then the same delayed handoff obeys

```text
||q_L||^2>=rho^2D/4,

M_L>=rho M/2.                                        (8)
```

The conclusion holds on either prescribed parity class of `L`, or on any
cofinal subsequence.

### Proof

The function `G` is itself a rational Boolean step function of mean `rho`.
Apply the two-density dilation proof of THM-2466 to `G`, rather than only to
`Q`.  Weak-star convergence of `G(N_L y)` gives (5).  The periodic primitive
of `G-rho` has sup norm at most `rho(1-rho)`, while the two-BV Fourier
covariance estimate gives `Var(G)Var(f)/(12N_L)`.  Taking the better bound
gives (6), and the reverse triangle inequality gives (7)--(8). QED.

## 2. Why the graft is lawful on the old root/deep table

Use the old oriented root chart

```text
phi(y,u)=(y+u)/13,                  u in F_13.        (9)
```

The physical delayed handoff is exactly root-constant:

```text
G(T^L phi(y,u))
 =G(13^(L-1)y+13^(L-1)u)
 =G(N_L y).                                           (10)
```

Thus the same Boolean factor may be inserted into both semantic root packets.
Its square is itself, so same-root disjointness and every pre-existing
positive service table remain lawful before any smoothing or Fourier
expansion.

It is also target-neutral.  A THM-2365 target co-shift changes an atomic
coordinate `w_i x` by an integer representative `theta_i/13`.  Inside an
event evaluated at `T^Lx`, the change is

```text
13^L theta_i/13=13^(L-1)theta_i in Z.                (11)
```

Every base factor is one-periodic.  The later terminal word gains the still
larger integer factor `13^(L+K-1)theta_i`.  Hence the whole `G(T^Lx)` is fixed
under the lawful target action.  Multiplication by it therefore preserves:

- the diagonal zero `H(t,s,t)=0`;
- the target quotient `(b,a+h)` and its circulant projection;
- every already installed old deep probe `Delta_r(x)`; and
- the orientation of the old common root fibre.

Disintegrating the old pointwise endpoint-table density `B(x)` through (9)
gives the function `b(y)` in (3).  Equation (10) turns the physical gated
drift into exactly the first integral in (4).  The same identity on the old
root-count product gives the second integral.

For septimal source phases, choose one parity of `L` before delaying, because
`13^L=(-1)^L mod 7`.  The cofinal-parity clause in Section 1 then retains the
same phase convention.  It does not preserve an independently preselected
septimal Fourier coefficient without the shifted-word construction of
THM-2442.

## 3. LRC owner-to-word consequence

On each of the `165` rows, THM-2349 supplies a universal depth-one owner
`E_j`, a finite row-dependent `K`, and a terminal word `Q_(j,sigma)` such
that

```text
rho_(j,sigma,K)
 =measure(E_j intersection T^(-K)Q_(j,sigma))
 >0.                                                  (12)
```

Use this pair in (2).  Suppose a fixed old packet has already retained
endpoint drift `D>0` and common-root service `M>0`; no support of that packet
on `E_j` is assumed.  Sections 1--2 give, for every sufficiently large `L`,
one physical packet with

```text
old endpoint/deep drift at x,

old root service at x,

T^Lx in E_j,

T^(L+K)x in Q_(j,sigma).                              (13)
```

Thus the extra filter is not merely a terminal stratum named after `j`.  It
contains the actual source and terminal events of a positive semantic
owner-to-word coupling at the future epoch `T^Lx`.

If the old packet comes from the THM-2459 bounds

```text
D>=D_0/63001,

M>=M_0/16384,
```

then (8) gives

```text
D_L>=rho^2 D_0/252004,

M_L>=rho M_0/32768.                                  (14)
```

These are the THM-2466 denominators with the terminal-word mean replaced by
the mass of the whole owner handoff.  The conclusion bypasses THM-2466's
prior-owner-support hypothesis by installing a later owner, not by proving
that the old packet was emitted by that owner.

## 4. Every future collision colour can carry the same old drift

There is a stronger fibre-product form.  Fix the THM-2471 first-collision
packet belonging to any one positive owner handoff.  In its notation, on the
collision base put

```text
c_s(y)=sum_(u in F_13)A(y,u+s)F(y,u)>=0,

g_k(y)
 =alpha_k(y)conjugate(phi_k(y))
 =(1/169)sum_s c_s(y)zeta_13^(-ks).                  (15)
```

The functions `c_s` are rational BV step densities.  They are marginals of
the finite Boolean ancestry stalk of THM-2471, and

```text
c_0=0,

J(k):=integral_T g_k!=0             for every k!=0. (16)
```

Let `d` be THM-2471's owner-normalization/penultimate-collision multiplier,
let `R=13^K` be that packet's fixed terminal clock, and choose `L>=K+1`.
Its physical collision base is `13d x`.  For the old root-base drift `b(y)`
from (3), define

```text
q_(k,L)
 =integral_T g_k(d13^L y)b(y)dy.                    (17)
```

The old physical point is `x=(y+u)/13`.  A collision packet delayed by `L`
has base

```text
13d T^Lx=d13^L(y+u)=d13^L y mod 1,                  (18)
```

so (17) is root-constant on the old fibre.  The same multiplier kills every
old target co-shift:

```text
13d13^L(theta/13)=d13^Ltheta in Z.                  (19)
```

The base calculation alone is not enough; the complete Boolean stalk must
also be target-neutral.  Substitute the delayed collision base

```text
y_col=13d13^Lx
```

into THM-2471's stalk coordinates.  Up to sheet constants, their leading
old-time terms are

```text
X_(u,a)=13^Lx+sheet,              current/word leg,

Z_(u,a,b)=13^(L-K)x+sheet,        owner ancestry leg,

Y_(v,e)=13^Lx+sheet.              bare source leg            (19a)
```

The smallest multiplier is `13^(L-K)`, which is divisible by thirteen
because `L>=K+1`.  Hence every old `/13` target co-shift becomes an integer
in every atomic factor of the stalk.  The whole Boolean fibre product, not
only its collision base, is target-neutral.

For complex BV functions the same Fourier covariance proof gives

```text
||q_(k,L)-J(k)q||
 <=Var(g_k)Var_H(b)/(12d13^L).                       (20)
```

There are only twelve nonzero colours.  Equations (16) and (20) therefore
give one common sufficiently large `L` for which, simultaneously for every
`k!=0`,

```text
||q_(k,L)||>=|J(k)|sqrt(D)/2>0.                     (21)
```

This does not multiply a physical packet by an unexplained complex observer.
First form the nonnegative products with each `c_s(d13^Ly)`.  Each such
product is the marginal of the explicit Boolean ancestry fibre product.
Only after that physical construction take the finite transform in `s` in
(15).  Hence (21) says that one old nonzero target/deep drift survives in
every future collision colour on one Boolean-before-DFT stalk.

There is an exact character-label corollary.  THM-2365 supplies an old
nonzero finite component

```text
chi=(a,b,h),                    a!=0,

(b,a+h)!=(0,0).                                     (21a)
```

Take `b(y)` in (17) to be the scalar root-base density of this fixed
component.  Its integral is the nonzero coefficient `B(a,b,h)`.  Equations
(20)--(21), with the error threshold chosen for this scalar coefficient,
retain that **same** old component in every future collision colour.  In
particular choose the future colour

```text
k=a.                                                 (21b)
```

The two-scale joint tensor then has diagonal equality of the printed
character labels.  This is not an identification of the two root fibres and
does not retain the same ordinary frequency `X`, deep multiplier `m`, or
relation address.  The old `a` and future `k` live at the two scales displayed
in (18); only their finite character labels have been aligned.

The direct graft of Section 3 and the collision graft here have different
semantic strengths.  Section 3 puts the literal old orbit point `T^Lx` in
`E_j`.  Section 4 attaches a canonical future collision ancestry stalk whose
base is (18).  Its word-current branches occur at old time `L` and their
owner-ancestry branches at old time `L-K`, but the averaging does not assert
that the literal orbit point `T^Lx` is a distinguished source sheet.

## 5. Sharp rebase boundary: the past deep probe becomes a sheet weight

The preceding sections deliberately keep the old endpoint and deep probe at
the old physical point `x`.  If instead one rebases at

```text
z=T^Lx,
```

transfer adjointness gives

```text
integral_T B(x)G(T^Lx)dx
 =integral_T (P^L B)(z)G(z)dz.                      (22)
```

The old packet has become the ancestry average `P^L B`.  For an old deepest
probe

```text
Delta_(C,r)(x)=d_1(Cx-r/13),

C=13^lambda u,                 13 does not divide u,
```

the exact transferred expression is

```text
P^L(B Delta_(C,r))(z)
 =1/13^L sum_(a mod 13^L)
   B((z+a)/13^L)
   d_1(C(z+a)/13^L-r/13).                            (23)
```

If `L<=lambda`, the probe is independent of `a` and factors as

```text
d_1((C/13^L)z-r/13) P^L B(z).                       (24)
```

Even here its coefficient has lost `L` levels of thirteen-adic depth.

If `L>lambda`, the phase in (23) depends on

```text
a mod 13^(L-lambda).                                (25)
```

This residue is essential.  Take `z=0`, `a_0=0`, and

```text
a_1=13^(L-lambda-1)u^(-1) mod 13^(L-lambda).
```

The two phases are respectively

```text
0,                    1/13 mod 1.                  (26)
```

Here take the probe root `r=0`; for general `r` the two displayed probe
arguments differ by `1/13`.  At the strict `1/14` danger radius with `r=0`,
the first is dangerous and the second is safe.  Hence no deep-root map
depending only on `z` can replace (23).  A
future-owner rebase must retain the ancestry residue (25), exactly the kind
of sheet coordinate isolated in THM-2471 Section 7.

The mixing horizon in (7) can be arbitrarily larger than the finite
`lambda`.  Therefore no uniform choice of `L` simultaneously gives mixing
and a sheet-free new-owner deep probe.

### The time/charge cocycle

The same natural-extension coordinate explains both the gain and the loss.
Write

```text
x=(z+a)/13^L,                  a mod 13^L.           (27)
```

An old lawful target shift by `theta/13` sends

```text
x -> x-theta/13,

z -> z,

a -> a-13^(L-1)theta                  mod 13^L.      (28)
```

Thus the future base `z=T^Lx` is fixed and only the ancestry digit moves.
This is exactly why a late handoff depending on `z` is target-neutral and
mixes cleanly with the old charged table.

A genuine target shift on the future owner leg would instead require

```text
z -> z-theta/13,

x -> x-theta/13^(L+1).                              (29)
```

The shift in (29) is not the old THM-2365 target action.  Hence late renewal
cannot transfer the old target charge into the future owner leg: it retains
the old charged vector while adding a neutral future coupling.  When
`L>lambda`, the same moving digit in (28) is read by the old deep comb through
`a mod 13^(L-lambda)`, as in (25).  Target neutrality and the deep-sheet loss
are therefore two faces of one cocycle, not unrelated caveats.

## 6. Crossed-product and owner-loop obstruction

Let `mathcal L=P_13` be the Perron operator and let `U h=h composed_with T`
be Koopman pullback.  For every bounded mask `P_omega`,

```text
M_(P_omega) mathcal L^K
 =mathcal L^K M_(U^K P_omega).                      (30)
```

Thus arrival refinement and source refinement differ by the exact temporal
commutator

```text
M_(P_omega)mathcal L^K M_e
 -mathcal L^K M_(P_omega)M_e

=mathcal L^K M_e M_(U^K P_omega-P_omega).           (31)
```

Let

```text
d lambda(x)
 =rho^(-1)e(x)Q(T^Kx)dx
```

be the probability measure on the positive handoff source, and work on
`L^2(lambda)`.  Suppose `eP_O=e` and `eP_omega=0` for `omega!=O`.  Define
the two representations of the finite atom algebra by

```text
pi_src(P_omega)=M_(P_omega),

pi_arr(P_omega)=M_(U^K P_omega).                    (32)
```

The source representation is scalar: `pi_src(P_O)=I` and all its other atom
projections vanish.  The arrival projections can split the same source into
several positive temporal cells.  There exists a bounded **onto** operator
`S:L^2(lambda)->L^2(lambda)` satisfying

```text
S pi_src(a)=pi_arr(a)S
```

for every atom-algebra element `a` if and only if

```text
P_O(T^Kx)=1                         lambda-almost everywhere. (33)
```

Indeed, for `omega!=O` the source projection is zero.  Intertwining would
make the arrival projection annihilate the range; surjectivity makes that
arrival projection zero.  Summing the atom projections gives (33).
Conversely, if (33) holds, both representations are the same scalar one and
`S=I` works.  THM-2461's positive terminal words lie in different truth
states, so the nontrivial handoff is precisely on the non-intertwining side.

Aggregate `H`-drift also does not force drift on the owner loop from the
finite tensor identities alone.  On `F_13^3`, put

```text
H_O(r,s,t)=(1/4)1_(r-t=1),

H_N(r,s,t)=(1/4)1_((r,s,t)=(1,1,0)).                (34)
```

Both are nonnegative and their sum vanishes on `r=t`.  The owner tensor is
circulant, so

```text
D(H_O)=0.
```

The nonowner atom has no untwisted mass:

```text
sum_r H_N(r,0,0)=0,
```

but its orbit projection is the uniform value `1/(4*169)` on the `169` cells
with `r-t=1`.  Hence

```text
D(H_O+H_N)
 =(1/16)(168/169)/13^3
 =21/742586
>0.                                                  (35)
```

This is an exact hostile at the current nonnegative-tensor interface, not a
physical scalar-cover row.  It proves that the direct graft and future
collision graft do not by themselves establish `q_O!=0`.

## 7. Exact gain and stopping boundary

The proved candidate chain is

```text
fixed old nonzero endpoint drift + old root service
  -> delay one whole positive owner handoff
  -> retain both on the old root base
  -> install a literal future owner-to-word edge;

fixed old nonzero endpoint drift
  -> delay one canonical first-collision stalk
  -> retain that drift in all twelve future collision colours.    (36)
```

The construction removes prior owner support as a prerequisite for adding a
future semantic edge.  It does **not** prove any of the following:

- the old endpoint packet was emitted by the future owner;
- the old collision root equals a future owner/deep root;
- the old deep probe descends after rebasing at the owner;
- the owner-loop drift `q_O` is nonzero;
- a prescribed first-expiration clock or old exact relation address survives;
- any scalar row is excluded.

The remaining cospan is narrower.  One must either retain the ancestry sheet
in (23)--(25), use the forward graft without rebasing, or add a theorem that
couples the old target/deep current to the future collision stalk on more than
their common base.

## 8. Exact companion

Run

```text
python3 04-computation/lrc14_delayed_owner_handoff_graft_thm2478.py
python3 -O 04-computation/lrc14_delayed_owner_handoff_graft_thm2478.py
```

The dependency-free `Fraction` companion:

- integrates a positive owner-to-word block and fixed drift/service densities
  exactly on rational step grids through four delays;
- checks both BV invoices and simultaneous nonvanishing;
- constructs five Boolean root-shift cells, verifies all twelve exact
  cyclotomic collision colours, and couples each to the same old drift;
- verifies old-root constancy and target neutrality on the physical
  collision-base multiplier;
- tests the sheet-free and sheet-essential sides of (23)--(26); and
- reproduces the exact `21/742586` aggregate-drift/owner-loop hostile.

Both runs must reproduce

```text
05-knowledge/results/lrc14_delayed_owner_handoff_graft_thm2478.out
```

byte-for-byte.  Promotion to proved canon remains contingent on an independent
hostile audit of the temporal orientation and collision-base normalization.
