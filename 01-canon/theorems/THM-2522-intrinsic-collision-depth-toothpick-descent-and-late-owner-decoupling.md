---
id: THM-2522
title: "Intrinsic collision-depth toothpick descent and late-owner decoupling"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For a rational step response F, the m-th
  base-thirteen Perron digit is the replica defect
  (I-U P_13)P_(13^m)F.  Its spectrum consists exactly of the nonzero
  Fourier modes of F having 13-adic valuation m, and its jump current is
  exactly the true level-m endpoint current minus the replication of the
  next endpoint current.  These digits give a unique orthogonal
  reconstruction of F.  The digit at m is also exactly the K_13
  toothpick drift of histories which first collide at depth L=m+1; when
  it is nonzero, rationality makes all twelve last-digit colours nonzero.
  A later positive owner may be placed arbitrarily far after this fixed
  collision and sees every colour beyond an explicit BV threshold.
  For every THM-2349 shallow owner-word event on the 165 live rows,
  support in the unit guard carrier forces the level-zero digit, while
  support in a depth-one danger comb also forces the level-one digit.
  Thus the first collision is uniformly L=1, all twelve colours fire,
  and both the old shallow septimal bank and the old deep c_3 bank are
  retained.  This does
  not orient the antipodal toothpick, couple its positive norm to the
  signed target/deep current, prove owner-loop drift, remove a row, or
  prove LRC(14).
source: codex-2026-07-27-intrinsic-collision-depth
depends_on:
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2519-last-digit-collision-drift-and-k13-dirichlet-boundary
  - THM-2520-rational-jump-crt-dichotomy-and-delayed-owner-forcing
related:
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
  - THM-2518-perron-inverse-branch-owner-word-cospan-recovery
script: 04-computation/lrc14_intrinsic_collision_depth_thm2522.py
output: 05-knowledge/results/lrc14_intrinsic_collision_depth_thm2522.out
script_sha256: f26609a5fbe46fa334453fec113b78e9dd5a3a60ddefff55ff9a3965d46ef941
output_sha256: 79832e1a28117bbd6bee377b5f0558a7780933e6bf5490f4c6d9acaf5700b21a
hash_basis: working-tree bytes (LF)
---

# THM-2522 -- the collision clock is a Perron digit, not an owner delay

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The Perron tower has two indices which must not be conflated:

```text
m = the scale of a nonconstant base-13 replica defect;
L=m+1 = the first-collision time of its unit toothpicks;
R = an arbitrary later delay used only to sample a common owner.       (1)
```

The phrase “first nonconstant Perron level” is otherwise ambiguous:
`P_(13^m)F` itself can be nonconstant at every level.  The intrinsic
object is its **digit**, the part which fails to replicate from the next
level.  This theorem identifies that digit in four equivalent ways:

```text
replica-defect function
  = unreplicated endpoint current
  = exact 13-adic Fourier stratum
  = positive depth-(m+1) K_13 toothpick energy.                (2)
```

## 1. The exact Perron--replication recursion

Work on `T=R/Z` with normalized Haar measure and put

```text
T(x)=13x mod 1,

(Uf)(x)=f(Tx),

(Pf)(y)=1/13 sum_(r=0)^12 f((y+r)/13).                         (3)
```

Thus `U` is Koopman pullback, `P=P_13` is its Perron adjoint, and

```text
PU=I,
P^m=P_(13^m),
E:=UP,
(Ef)(y)=1/13 sum_(u in F_13)f(y+u/13).                         (4)
```

The operator `E` is the orthogonal projection onto the functions invariant
under `y -> y+1/13`.  For a real periodic step function `F`, define

```text
F_m=P^mF,
D_m=(I-E)F_m=F_m-UF_(m+1).                                    (5)
```

Then `P D_m=0`, `integral D_m=0`, and the exact function recursion is

```text
F_m=D_m+UF_(m+1).                                             (6)
```

There is an equally intrinsic endpoint-current recursion.  Let

```text
mu_m=dF_m                                                     (7)
```

be the distributional derivative formed from the **true net jumps** of
`F_m`.  If `T_*` denotes pushforward and

```text
Rep(sum_z a_z delta_z)
 =sum_z a_z sum_(Tx=z)delta_x                                 (8)
```

is unnormalized inverse-image replication, then

```text
mu_(m+1)=1/13 T_*mu_m,

mu_m=1/13^m (T^m)_*mu_0,                                     (9)

Gamma_m:=mu_m-Rep(mu_(m+1))=dD_m.                            (10)
```

Indeed, a branch in `PF` crosses a jump at `x` exactly when `y=13x`, and
the chain rule for `UF_(m+1)` replicates each child jump over all thirteen
preimages.  Thus

```text
mu_m=Gamma_m+Rep(mu_(m+1))                                   (11)
```

is the requested replicated-endpoint-current descent.  Combining repeated
endpoint presentations before (7) is load-bearing: artificial zero jumps
are not currents.

Because `D_m` has mean zero, the following equivalence is exact:

```text
Gamma_m=0  iff  D_m=0.                                       (12)
```

## 2. Fourier valuation, first depth, and unique reconstruction

Use

```text
Fhat(n)=integral_T F(x) exp(-2 pi i n x)dx.                   (13)
```

Equations (3)--(5) give

```text
(F_m)^hat(k)=Fhat(13^m k),                                   (14)

(D_m)^hat(k)
 =Fhat(13^m k),                    13 does not divide k,
 =0,                               13 divides k.              (15)
```

Thus `D_m` contains exactly the original Fourier modes whose nonzero
frequency has valuation `m`, after division by `13^m`.  For nonconstant
`F`, put

```text
m_*(F)=min{nu_13(n):n in Z minus {0}, Fhat(n)!=0}.             (16)
```

Fourier uniqueness makes the set in (16) nonempty.  Equations (10), (12),
and (15) prove

```text
m_*(F)
 =min{m>=0:D_m!=0}
 =min{m>=0:Gamma_m!=0}.                                      (17)
```

For a constant function set `m_*=infinity`; every digit and current defect
then vanishes.  The minimum in (17) is unique, but later active digits need
not form an interval.

There is also a coordinate-free characterization.  For every `r>=0`,

```text
D_0=...=D_(r-1)=0

iff F=U^rF_r

iff F is invariant under translation by 1/13^r.              (18)
```

Consequently, when `m_*=r`, `F` is an `r`-fold replica but not an
`(r+1)`-fold replica.  This is the physical meaning of the first digit.

Iterating (6) gives, for every `M>=1`,

```text
F-integral F
 =sum_(m=0)^(M-1) U^mD_m
   +U^M(F_M-integral F).                                     (19)
```

The summands `U^mD_m` are pairwise orthogonal: if `i<j`, adjointness and
`PD_i=0` give

```text
<U^iD_i,U^jD_j>=<P^(j-i)D_i,D_j>=0.                          (20)
```

Moreover,

```text
||F_M-integral F||_2^2
 =sum_(k!=0)|Fhat(13^M k)|^2 -> 0.                            (21)
```

Hence (19) converges in `L^2` to the unique orthogonal expansion

```text
F-integral F=sum_(m>=0)U^mD_m,

||F-integral F||_2^2=sum_(m>=0)||D_m||_2^2.                  (22)
```

Every original nonzero frequency belongs to exactly one depth in (22).
This is a base-thirteen Wold decomposition, not an asymptotic heuristic.

## 3. The digit is exactly one first-collision toothpick bank

Fix `m>=0` and put

```text
M=13^m,                    N=13M,                    L=m+1.   (23)
```

For `u in F_13`, define the unweighted fixed-last-digit table

```text
B^(m)_u
 =1/M sum_(e=0)^(M-1) integral_T
    F(x)F(x+(u+13e)/N)dx.                                    (24)
```

The exact THM-2519 branch reduction is

```text
B^(m)_u=integral_T F_m(y)F_m(y+u/13)dy.                      (25)
```

If `u!=0`, every address `d=u+13e` is a thirteen-adic unit.  The two
histories in (24) remain distinct through time `m` and first coalesce at

```text
L=m+1.                                                       (26)
```

This proves the off-by-one in (1).

Let `zeta=exp(2 pi i/13)` and split the digit into last-character pieces:

```text
D_(m,a)(y)
 =1/13 sum_(u in F_13)zeta^(-au)F_m(y+u/13),
                                      a in F_13^*.            (27)
```

Its Fourier series consists exactly of the terms

```text
k=a mod 13,                  Fhat(13^m k).                    (28)
```

The twelve pieces are pairwise orthogonal and

```text
D_m=sum_(a!=0)D_(m,a).                                       (29)
```

Taking the normalized `u`-transform of (25) gives

```text
Bhat^(m)(a)
 :=1/13 sum_u B^(m)_u zeta^(-au)
 =||D_(m,a)||_2^2
 >=0.                                                        (30)
```

In particular,

```text
B^(m)_0-1/13 sum_u B^(m)_u
 =||D_m||_2^2
 =sum_(a!=0)Bhat^(m)(a).                                     (31)
```

Equations (17), (26), and (31) identify the first positive unweighted
collision depth:

```text
L_*(F)=m_*(F)+1.                                             (32)
```

The collision energy is the `K_13` Dirichlet energy of the thirteen
predecessor values.  Each nonzero `u` is one Hamiltonian toothpick cycle;
`u` and `-u` have the same energy.  The bank is therefore positive but
intrinsically unoriented.

### Rational all-or-all phase law

Assume now that `F` is rational-valued.  Every `F_m(y+u/13)` is rational
away from finitely many jumps.  If one `D_(m,a)` vanishes for `a!=0`, then

```text
sum_(u=0)^12 F_m(y+u/13)X^u                                  (33)
```

vanishes at a primitive thirteenth root for almost every `y`.  Irreducibility
of `Phi_13` forces all thirteen coefficients in (33) to be equal.  Therefore

```text
D_m=0
 iff D_(m,a)=0 for one a!=0
 iff D_(m,a)=0 for every a!=0;                               (34)

D_m!=0
 iff Bhat^(m)(a)>0 for every a!=0.                           (35)
```

Thus every active rational digit fires the complete twelve-colour bank on
one common physical collision table.  Rationality is unnecessary for the
sum (31), but necessary for (34)--(35).

## 4. The finite endpoint horizon

Let the true jumps of `F` lie on a rational grid

```text
x_j=b_j/D,
D=13^nu d,                         gcd(d,13)=1,                (36)
```

and form the prime-to-thirteen current

```text
C_t=sum_(j:b_j=t mod d)Delta_j,             t in C_d.         (37)
```

For `m>=nu`, put

```text
U_m=13^(m-nu) mod d.                                           (38)
```

THM-2520's exact jump formula becomes

```text
mu_m
 =1/13^m sum_(v in C_d)C_(U_m^(-1)v) delta_(v/d).             (39)
```

Therefore the tower has a finite phase-bank horizon:

```text
C=0
  iff F_m=integral F for every m>=nu
  iff D_m=0 for every m>=nu;

C!=0
  iff D_(m,a)!=0 for every m>=nu and every a!=0.              (40)
```

Before `nu`, digits may vanish and reappear.  After `nu`, the bank is
permanently off or permanently on.  Its labelled current is periodic in
`m` with period dividing

```text
ord_d(13),                                                    (41)
```

while its amplitude scales by `13^(-m)`.  In particular, a nonconstant
rational step function always has

```text
m_*(F)<=nu.                                                   (42)
```

To find the first depth using only endpoints, compute the finite defects
`Gamma_0,...,Gamma_(nu-1)` by (9)--(11), then use `C` for the horizon
branch.  If all those defects and `C` vanish, (18) and (40) imply that `F`
is constant.

There is also an exact quantitative invoice.  Put `Q=13d`,

```text
Chat(b)=sum_t C_t exp(-2 pi i bt/d),
S_C=sum_t C_t^2,                                               (43)
```

and let `r_(m,a,b)` be the unique representative in `{1,...,Q-1}` with

```text
r_(m,a,b)=a mod 13,
U_m r_(m,a,b)=b mod d.                                        (44)
```

For `a!=0`,

```text
||D_(m,a)||_2^2
 =1/(4*13^(2m)*Q^2) sum_(b in C_d)
    |Chat(b)|^2 csc^2(pi r_(m,a,b)/Q)                         (45)

 >=S_C/(4*13^2*13^(2m)*d).                                   (46)
```

This is the exact cost of retaining a last digit after `m` Perron levels.

## 5. Owner time is independent of collision time

Fix one active digit `D_m`, hence the physical collision depth `L=m+1`.
Let

```text
0<=G<=1,          rho=integral G>0,          W=Var(G),         (47)
```

be any positive BV future owner or whole owner--word block.  Put the owner
`R>=0` steps after the collision:

```text
B^[m,R]_u
 =1/M sum_(e=0)^(M-1) integral_T
    G(13^(L+R)x)
    F(x)F(x+(u+13e)/13^L)dx.                                 (48)
```

The same inverse-branch calculation as (25) gives

```text
B^[m,R]_u
 =integral_T G(13^(R+1)y)F_m(y)F_m(y+u/13)dy,                (49)

Bhat^[m,R](a)
 =integral_T G(13^(R+1)y)|D_(m,a)(y)|^2dy
 >=0.                                                        (50)
```

The factor in (48) is common to both histories because they already agree
at time `L`.  Moving it later changes neither their address nor their first
collision time.

Let

```text
B=||F||_infinity,             V=Var(F),
K=13^(R+1).                                                   (51)
```

Perron contraction and the phase average give

```text
||D_(m,a)||_infinity<=B,
Var(D_(m,a))<=V/13^m,
Var(|D_(m,a)|^2)<=2BV/13^m.                                  (52)
```

The periodic BV covariance bound therefore yields

```text
|Bhat^[m,R](a)-rho||D_(m,a)||_2^2|
 <=BVW/(6*13^m*K).                                           (53)
```

For an active rational digit define

```text
e_m=min_(a!=0)||D_(m,a)||_2^2>0.                              (54)
```

All twelve owner-weighted colours are strictly positive whenever

```text
K>BVW/(6 rho 13^m e_m).                                      (55)
```

Doubling the right side gives the uniform lower bound

```text
Bhat^[m,R](a)>=rho e_m/2.                                    (56)
```

At the endpoint horizon, (46) recovers THM-2520's explicit sufficient
condition

```text
K>2*13^2*d*13^m*B*V*W/(3 rho S_C).                           (57)
```

When `W=0`, (50) is already `rho||D_(m,a)||_2^2`.  Because every parity
class of `R` is cofinal, a septimal convention can be fixed before choosing
the delay.  Conversely, if `D_(m,a)=0`, no owner and no delay can create
that colour.  Late mixing samples an existing collision digit; it does not
change or manufacture the collision depth.

If `G=1_E(Q composed T^K_0)` is a positive owner--word block, shift the
whole block in (48).  Its source occurs at `L+R` and its terminal word at
`L+R+K_0`.  This is the forward, no-rebase interpretation of THM-2478.

## 6. Carrier geometry bounds and then fixes the collision depth

For the LRC application, write

```text
D_c={x in T:||cx||<1/14}.                                    (58)
```

First let

```text
c=13^s a,                      s>=0,             13 does not divide a,

0<=F not identically 0,
support(F) subset D_c                                            (59)
```

up to null endpoints.  Put `h=P^sF`.  Nonnegativity and mean preservation
give `h>=0` and `h!=0`.  Moreover,

```text
h(y)>0  ->  F((y+r)/13^s)>0 for some r
          -> ||a(y+r)||<1/14
          -> y in D_a.                                      (60)
```

Thus `support(h) subset D_a`.

Suppose its level-`s` digit vanished.  Then `h=Eh`, so `h` would be
invariant under every `y -> y+u/13`.  A nonzero value would force its whole
translation orbit into `D_a`.  But

```text
intersection_(u in F_13)(D_a-u/13)=emptyset.                 (61)
```

Indeed, multiplication by the unit `a` permutes the thirteen equally
spaced phases `u/13`, and they cannot all lie inside the strict interval
of radius `1/14` about an integer.  Equations (60)--(61) contradict
`h!=0`.  Therefore

```text
D_s=(I-E)P^sF!=0,

B^(s)_0-1/13 sum_u B^(s)_u>0,

m_*(F)<=s.                                                    (62)
```

If `F` is rational step, all twelve depth-`s+1` colours are positive by
(35), and every sufficiently late positive owner sees all twelve by (55).

There is a sharp first-collision dichotomy.  If `m_*<s`, then

```text
L_*=m_*+1<=s,                                                 (63)
```

so the septimal phase bank belonging to the coefficient `c=13^s a` is
retained.  If `m_*=s`, equations (17)--(18) give

```text
F=U^s(P^sF) exactly,

L_*=s+1,                                                      (64)
```

and the same septimal bank lies exactly one level beyond its retention
horizon.  Thus source-sheet loss is confined to the exact `s`-fold replica
equality branch.  The full danger comb

```text
F=1_(D_(13^s a))=U^s1_(D_a)                                  (65)
```

is the sharp equality control.

For the live shallow case `s=1`, (62)--(64) reduce to two branches:

```text
D_0!=0:
  first collision L_*=1; shallow septimal and old deep banks retained;

D_0=0:
  F=U(PF), first collision L_*=2;
  old deep bank retained, shallow septimal bank not retained. (66)
```

The second line is an exact once-replicated equality, not a generic loss.

### All 165 live rows

THM-2349 supplies, on every one of the `165` first-depth-one rows, a shallow
owner `j`, a nonempty terminal word `sigma`, and a finite clock `k` such that

```text
F=1_(E_(j,sigma,k))
 =1_(E_j)(x)1_(Q_(j,sigma))(13^k x)                           (67)
```

is a nonzero rational Boolean event with

```text
support(F) subset E_j subset D_(c_j),
nu_13(c_j)=1.                                                 (68)
```

This already makes the depth-two digit nonzero by (62).  The same event has
a stronger carrier.  In THM-2349's notation,

```text
E_j subset A_0 subset C_H,

C_H={x:||Hx||>1/7},                  13 does not divide H.    (69)
```

The unit guard carrier contains no complete root orbit:

```text
intersection_(u in F_13)(C_H-u/13)=emptyset.                 (70)
```

Indeed, multiplication by `H` permutes the thirteen-grid, and one of its
phases is at distance at most `1/26<1/7` from an integer.  If `D_0=0`, then
`F=EF` would be root-orbit invariant, contradicting (69)--(70) and `F!=0`.
Therefore, uniformly on all `165` rows,

```text
D_0!=0,                 m_*=0,                 L_*=1.         (71)
```

Rationality and (35) now give all twelve depth-one collision colours.  The
shallow support argument also gives `D_1!=0`, so the event actually has two
consecutive full banks at `L=1,2`; the first one is the stronger semantic
choice.  Since THM-2349 has `k>=2`, its word factor is literally common on
each depth-one pair:

```text
Q(13^k(x+u/13))=Q(13^kx).                                   (71a)
```

Any positive BV future owner block can be added at an arbitrarily late
common time by Section 5.  Moving that block does not move the collision
from `L=1`.

The bank ledger at the first collision is exact.  Every address with
`u!=0` is a unit, so THM-2519's phase-bank criterion says:

```text
septimal bank {d(c_jx-ell/7)} retained only if
  L<=nu_13(c_j);

deep bank {d(c_3x-r/13)} retained if
  L<=nu_13(c_3)+1.                                           (72)
```

Here `L=1`, `nu_13(c_j)=1`, and `nu_13(c_3)>1`.  Therefore

```text
old deep c_3 bank:             retained by label permutation;
old shallow septimal bank:     retained by label permutation;
new collision F_13 bank:       all twelve colours nonzero.   (73)
```

The new collision colour and old deep label are different typed
coordinates.  The shallow septimal character is also retained as a bank,
but no theorem here identifies it with the new collision colour or proves a
signed mixed target/deep current.

## 7. Sharp controls and stopping boundary

Every qualification above has an exact hostile.

1. **Digits are not monotone.**  If

   ```text
   q=1_[0,1/13),
   F=q+2(q composed T^2),                                    (74)
   ```

   then only `D_0` and `D_2` are nonzero.  A digit can vanish and reappear.

2. **The off-by-one is sharp.**  For `F=q composed T^r`,

   ```text
   m_*=r,                         L_*=r+1.                    (75)
   ```

3. **Nonnegativity is load-bearing in Section 6.**  On two inverse branches
   over `D_a`, put equal positive and negative copies.  The resulting
   nonzero signed function is supported in `D_(13a)` but has `PF=0`.

4. **The unit hypothesis is load-bearing.**  If `13|a`, then `D_a` itself
   is `1/13`-invariant.  Taking `F=1_(D_(13a))` makes the proposed level-one
   digit zero.

5. **The unit guard and its radius are load-bearing.**  If `13|H`, then
   `C_H` itself is `1/13`-invariant, so its indicator has `D_0=0`.  If the
   excluded central radius is reduced below `1/26`, the half-grid orbit

   ```text
   {(2u+1)/26:u in F_13}
   ```

   lies wholly in the weakened carrier.  Its finite translation-invariant
   neighbourhood can be chosen inside that carrier, giving a nonzero step
   hostile with `D_0=0`.  Neither hostile is available at the live unit
   guard radius `1/7`.

6. **Rationality is load-bearing for all twelve colours.**  On the thirteen
   inverse branches over one interval, the real values

   ```text
   2+2 cos(2 pi a_0u/13)                                    (76)
   ```

   retain only the two colours `+-a_0`.

7. **Owner delay is load-bearing but cannot change depth.**  With

   ```text
   F=1_[0,1/14),                 G=1_[13/14,1),              (77)
   ```

   the depth-one digit is nonzero, while the owner at `R=0` misses every
   predecessor.  Later dilates see the digit.  If the digit were zero,
   every delay would remain zero.

8. **Orientation is absent.**  The square table always satisfies
   `B_u=B_(-u)`.  Positivity in (30) is not an oriented owner-loop current.

The all-row gain (71)--(73) therefore closes the zero-drift alternative at
the intrinsic first collision of the positive shallow owner--word event.  It
does not assert that THM-2520's eventual coprime current `C` is nonzero, and
it does **not**
show that a signed ANOVA/cut/target coefficient is nonzero on that event,
that the old packet was emitted by the late owner, that source and arrival
atoms agree, or that a scalar cover row is impossible.  No row is removed;
LRC(14) remains open.

## 8. Exact companion

Run

```bash
python3 04-computation/lrc14_intrinsic_collision_depth_thm2522.py
python3 -O 04-computation/lrc14_intrinsic_collision_depth_thm2522.py
```

Both runs reproduce the stored transcript.  The dependency-free exact
companion checks:

- the digit pattern `nonzero, zero, nonzero, zero` for the gapped tower;
- the sharp pattern `zero, zero, nonzero` and first collision `L=3` for a
  twice-replicated tile;
- the function and true endpoint-current recursions, orthogonality, and
  exact energy telescope;
- direct versus reduced collision tables at all three levels;
- `P_13 1_(D_(13a))=1_(D_a)`, the empty thirteen-translate intersection,
  depth-two energy `144/1183`, and all sharp sign/unit controls;
- the empty unit-guard orbit, consecutive guarded energies
  `174/5915,8724/142805` at depths one and two, and the nonunit/small-radius
  guard hostiles;
- the undelayed owner miss and the first positive delayed bank;
- `40` late-current recursion cells, `540` CRT phase addresses, and the
  exact `ord_5(13)=4` label cycle at a `13^2*5` conductor.

The exact referee proves the identities on finite rational step towers and
prints the theorem's consequence objects, not only its hypotheses.

Two independent audits rederived the operator recursion, valuation shells,
off-by-one collision clock, BV constants, endpoint horizon, and the
THM-2349 application.  The second audit separately attacked the strengthened
unit-guard step and confirmed the empty root-orbit argument and the `L=1`
phase-bank ledger. **QED.**
