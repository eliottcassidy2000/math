---
id: THM-2286
title: "Endpoint-Prony lift banks and sharp owner-multiplier landing"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. If an indicator of a
  union of R circular intervals has positive Fourier energy in a nonzero
  residue class modulo m, then that class contains an actual nonzero atom
  of absolute frequency at most mR-1. Applied to THM-2269, every one of
  the 165 live scalar profiles has a selected marked-owner set with an
  exact bounded Fourier lift in every nonzero mod-13 residue. If the
  selected owner is c_1 and its THM-2266 normalized pair has owner
  multiplier d<=6, a 7-unit multiple of that exact pair frequency lands
  in the marked spectrum; d=7 is the sharp uniform boundary. Tiny generic
  phase separation can make any prescribed bounded relation visible in
  safe-interval approximants, but supplies no quantitative return. Exact
  homometric supports show that even the full power spectrum cannot name
  the lift. No scalar profile is excluded, and LRC(14) remains open.
source: codex-2026-07-25-endpoint-prony-owner-landing
depends_on:
  - THM-2199-effective-positive-subspace-rank-lift
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2266-depth-one-deep-pair-centered-signed-dual-and-relation-atlas
  - THM-2269-marked-expiration-root-spectrum-and-branch-state-no-go
related:
  - THM-2218-integral-fourier-hasse-guard-hole-transform-and-family-knapsack
  - THM-2278-two-shallow-proper-root-spectrum-and-gap-ancestry-activation
  - THM-2280-centered-polynomial-grid-avoidance-and-bounded-generic-keller-fibre
  - THM-2284-thirteen-adic-anchored-rank-three-plucker-lift
script: 04-computation/lrc14_endpoint_prony_owner_landing_thm2286.py
output: 05-knowledge/results/lrc14_endpoint_prony_owner_landing_thm2286.out
script_sha256: 846a945221c4fe8342f2a40c4b3d2fed570ccfd60dd7a3e934f64e32ab55ade9
output_sha256: f2819a4726fbb9d2387ff9d5ee64708e6db50a17afbea5ddd8973af96dd0eec1
hash_basis: working-tree bytes (LF)
---

# THM-2286 -- endpoint-Prony lift banks and sharp owner landing

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2269 proves that every nonzero residue class modulo thirteen carries
positive Fourier energy for a selected exclusive-owner flow. It correctly
stops before claiming that one integer frequency in the class is nonzero.
For an arbitrary measurable set that stopping boundary is real.

The selected flow is not arbitrary: it is a finite union of circular
intervals. Its endpoints turn each arithmetic progression of Fourier
coefficients into a finite exponential sum. A Vandermonde argument then
recovers an actual bounded lift:

```text
positive energy in q=k mod 13
  + at most R circular components
  -> some exact q=k mod 13 with |q|<=13R-1.           (1)
```

This solves the existence and finiteness of an integer lift. It does not
align that lift with a relation, owner transition, or cut.

## 1. Endpoint-Prony lemma

Let `G subset R/Z` be, up to null endpoints, a union of `R>=1` disjoint
circular intervals, and put

```text
f=1_G,

f_hat(q)=integral_(R/Z) f(x)exp(-2pi i qx)dx.         (2)
```

Choose oriented endpoints `a_j,b_j`, `1<=j<=R`. Direct integration gives,
for every nonzero integer `q`,

```text
2pi i q f_hat(q)
 =sum_(j=1)^R[
    exp(-2pi i q a_j)-exp(-2pi i q b_j)
  ].                                                  (3)
```

Fix an integer modulus `m>=2` and a nonzero residue

```text
k in {1,...,m-1}.                                    (4)
```

For `h in Z`, set

```text
A_h=2pi i(k+mh)f_hat(k+mh).                          (5)
```

Equation (3) writes this sequence as

```text
A_h=sum_e gamma_e z_e^h,                             (6)

z_e=exp(-2pi i m x_e),
gamma_e=sigma_e exp(-2pi i k x_e),                   (7)
```

where the `x_e` are the at most `2R` endpoints and
`sigma_e in {+1,-1}`. Combine equal nodes `z_e` and discard zero combined
coefficients. The result has `L<=2R` distinct nonzero nodes.

We use the elementary consecutive-zero lemma:

> If
>
> ```text
> A_h=sum_(ell=1)^L gamma_ell z_ell^h
> ```
>
> has distinct nonzero `z_ell` and is not the zero sequence, then it cannot
> vanish at `L` consecutive integers.

Indeed, after shifting the first index to zero, the `L` equations have
coefficient matrix

```text
[z_ell^j]_(0<=j<L,1<=ell<=L).                        (8)
```

Its determinant is

```text
product_ell z_ell^(h_0)
product_(i<j)(z_j-z_i)!=0.                           (9)
```

Thus `L` consecutive zeros force every `gamma_ell=0`.

Suppose now that the residue class in (4) has positive energy:

```text
sum_(q congruent k mod m)|f_hat(q)|^2>0.             (10)
```

Since `k+mh` is never zero, equations (5)--(6) show that `A_h` is not the
zero sequence. It therefore cannot vanish throughout the `2R` consecutive
indices

```text
h=-R,-R+1,...,R-1.                                  (11)
```

For some index in (11),

```text
f_hat(q)!=0,            q=k+mh.                     (12)
```

The endpoints of this window satisfy

```text
|k-mR|<=mR-1,
|k+m(R-1)|<=mR-1,                                   (13)
```

so every member of the window obeys the same bound. We have proved:

> **Endpoint-Prony lift lemma.** Under (2), (4), and (10), there is a
> nonzero integer `q` such that
>
> ```text
> q=k mod m,
> |q|<=mR-1,
> f_hat(q)!=0.                                       (14)
> ```

Coincident endpoint images modulo `1/m` create fewer nodes and only improve
the argument.

## 2. Component complexity of comb Boolean sets

For a positive integer `s`, each comb

```text
D_s={x:||sx||<1/14}
```

has `s` circular components and at most `2s` boundary points. Its complement
and the guard set

```text
C_H={x:||Hx||>1/7}
```

have the same boundary count as their defining combs.

Consider any Boolean combination of

```text
C_H,
D_(q_1),...,D_(q_5),
D_(c_1),D_(c_2),D_(c_3).                            (15)
```

Put

```text
S=H+sum_(i=1)^5 q_i+sum_(j=1)^3 c_j.                (16)
```

The union of all defining boundaries has at most `2S` points. Between
consecutive boundary points every Boolean predicate in (15) is constant.
A proper selected set therefore has at most

```text
S                                                       (17)
```

circular components. This is sharp for arbitrary selections of `2S`
cyclic cells.

The image of a connected circular interval under multiplication by an
integer is connected. Hence multiplication maps cannot increase the number
of components of a finite union.

## 3. Every marked residue has a bounded exact lift

Use THM-2269's selected exclusive-owner stratum. In the notation of
THM-2255,

```text
A_0=C_H minus union_(i=1)^5 D_(q_i),

E_j=A_0 intersection D_(c_j)
          minus union_(ell!=j)D_(c_ell).             (18)
```

For the selected positive label, write

```text
c_j=13^lambda u,

F=T^lambda(E_j),              T(x)=13x mod 1.        (19)
```

Equations (15)--(18) show that `E_j` has at most `S` components.
Multiplication by `13^lambda` does not increase that number, so

```text
F has at most S circular components.                 (20)
```

THM-2269 proves, separately for every `k=1,...,12`,

```text
sum_(q congruent k mod 13)|1_F_hat(q)|^2>0.          (21)
```

Apply (14) with `m=13` and `R<=S`. For every nonzero residue there is an
actual integer atom

```text
q_k=k mod 13,
|q_k|<=13S-1,
1_F_hat(q_k)!=0.                                     (22)
```

This holds on every one of the `150` strict and `15` repeated-first live
profiles. It depends only on proved THM-2269, not on the stronger labelled
gap spectrum of THM-2278.

There is also one universal, although enormous, bound. Let `V_*` be the
explicit `197`-digit primitive speed ceiling in equation (23) of THM-2199.
THM-2203 identifies the selected original coordinates as

```text
(8H,16q_1,...,16q_5,16c_1,16c_2,16c_3).            (23)
```

Thus

```text
S<=V_*/8+8V_*/16=5V_*/8,                            (24)
```

and every atom in (22) may be chosen with

```text
|q_k|<=13(5V_*/8)-1.                                (25)
```

The companion prints this exact integer. Equations (22)--(25) are a finite
lift bank, not a selection rule matching a prescribed relation frequency.

## 4. Pinned low owner atoms

The owner geometry gives a more precise statement at its first few
multiples. THM-2269 has

```text
F subset D_u,
measure(F)>0.                                        (26)
```

For `1<=n<=3`, every `x in F` satisfies

```text
cos(2pi n u x)>=cos(n pi/7)>0                       (27)
```

after taking the centered representative of `ux`. Therefore

```text
Re 1_F_hat(nu)
 >=cos(n pi/7)measure(F)
 >0.                                                  (28)
```

By conjugation, (28) also holds as nonvanishing at the negative
frequencies. Thus the exact owner atoms

```text
q in {+/-u,+/-2u,+/-3u}                             (29)
```

always fire. The cutoff `3` in this pointwise cosine proof is sharp because
`cos(4pi/7)<0`.

## 5. Sharp owner-multiplier landing

The previous result extends qualitatively through multiplier six. For a
nonnegative integrable `f` supported in `D_u`, define the normalized
pushforward

```text
(P_a f)(y)
 =(1/a)sum_(r=0)^(a-1) f((y+r)/a),                  (30)
```

for positive integers `a`. It preserves mass and satisfies

```text
Fourier(P_a f,n)=f_hat(an).                          (31)
```

Put

```text
g=P_u f.                                             (32)
```

Then `g` is nonnegative, nonzero, and supported in

```text
I=(-1/14,1/14).                                      (33)
```

Fix `1<=d<=6` and put

```text
h=P_d g.                                             (34)
```

Its support lies in the proper circular arc

```text
dI,                 measure(dI)=d/7<1.              (35)
```

Hence `h` cannot be a nonzero constant. Some nonzero Fourier coefficient
must fire:

```text
h_hat(n)=g_hat(dn)=f_hat(ndu)!=0.                    (36)
```

The multiplier `n` can moreover be chosen prime to seven. Suppose otherwise
that every nonzero Fourier coefficient of `h` were supported on `7Z`.
Then `h` would be `1/7`-periodic. The complement of the arc in (35) has
length

```text
1-d/7>=1/7.                                          (37)
```

It contains a full period on which `h` vanishes, so periodicity would force
`h=0` almost everywhere, contradicting preservation of its positive mass.
Therefore:

> **Small owner-multiplier landing.** If `f>=0` is nonzero and supported
> in `D_u`, and `1<=d<=6`, then there is a nonzero integer `n`, with
>
> ```text
> 7 does not divide n,
> f_hat(ndu)!=0.                                     (38)
> ```

If `f=1_F` and `F` has at most `R` circular components, `n` can be bounded.
Choose a nonzero residue `k mod 7` in the Fourier support of `h`, and apply
the endpoint-Prony argument directly to

```text
f_hat(du(k+7j)).
```

There are at most `2R` endpoint nodes. A centered block of `2R` lifts gives

```text
0<|n|<=7R-1.                                         (39)
```

### Sharp failure from seven onward

For every `d>=7`, take

```text
g_d=1_[-1/(2d),1/(2d)).                              (40)
```

This interval lies in `D_1`. Multiplication by `d` maps it bijectively
almost everywhere onto the circle, so

```text
P_d g_d=1/d                                          (41)
```

almost everywhere. Equivalently,

```text
g_d_hat(nd)=0                    for every n!=0.      (42)
```

Thus the transition from `d=6` to `d=7` is the exact boundary:

```text
1<=d<=6:
  every positive owner-sector density has some exact d-multiple atom;

d>=7:
  one interval can kill every nonzero d-multiple atom.          (43)
```

No better residue-energy estimate can remove this obstruction, because
the hostile interval itself has positive mass inside the owner sector.

## 6. Conditional landing of the THM-2266 shallow pair

On each of the `120` interior profiles, THM-2266 gives a reduced coprime
pair which can be oriented as

```text
d u_1=e z,

gcd(d,e)=1,                  de<=757,                (44)
```

where `c_1=13u_1` and `z` is `H` or one of the `q_i`. Translating (44) to
the actual scalar row gives one of

```text
13e H-d c_1=0,

d c_1-13e q_i=0.                                    (45)
```

Suppose that the positive marked owner chosen by THM-2269 is this shallow
label `c_1`. Then its pre-expiration set satisfies

```text
F subset D_(u_1),             measure(F)>0.          (46)
```

If the oriented owner multiplier in (44) obeys

```text
d<=6,                                                (47)
```

equation (38) supplies a seven-unit `n` for which

```text
1_F_hat(ndu_1)
 =1_F_hat(ne z)
 !=0.                                                 (48)
```

Thus a multiple of the exact normalized pair frequency lands in the
selected marked spectrum. This is the first direct relation-to-marked-atom
interface in the scalar atlas.

There are exactly `3643` oriented coprime pairs with `de<=757`. Of these,

```text
1372
```

have `d<=6`, with counts

```text
d=1,2,3,4,5,6:
757,189,168,95,121,42.                              (49)
```

This is `37.66...%` of the oriented **shape atlas**, not of the `120` live
rows. Equation (48) also requires the independent owner-alignment
hypothesis (46). No row count or profile exclusion follows.

For centered length-`5/7` and length-`6/7` interval approximants, a Fourier
coefficient can vanish when its index is divisible by seven. The stricter
subatlas in which the other factor `e` is also a seven-unit has size

```text
1176,                                                (50)
```

with histogram

```text
649,162,144,81,104,36.                               (51)
```

The next lemma shows that this additional coefficient restriction is an
artifact of using unperturbed centered intervals when only spectral
visibility is required.

## 7. Generic phase separation for a prescribed relation

Let `w_1,...,w_s` be nonzero integer speeds and let

```text
m.w=0                                                (52)
```

be a prescribed nonzero integer relation. Partition the labels into blocks
`A,B` with nonzero partial frequency

```text
K=sum_(i in A)m_i w_i
  =-sum_(i in B)m_i w_i
  !=0.                                                (53)
```

Start from the relevant guard or ordinary safe interval `J_i`. Choose an
arbitrarily close inner interval `J_i'` whose length avoids

```text
m_i measure(J_i') in Z                              (54)
```

for every `m_i!=0`. There are only finitely many forbidden lengths. Smooth
`1_(J_i')` by a normalized squared-Fejer kernel whose degree is at least
`|m_i|`, obtaining

```text
0<=p_i<=1,
Fourier(p_i,m_i)!=0.                                 (55)
```

Translate the coordinate polynomial by an independent phase `theta_i`.
The Fourier coefficient at `K` of

```text
product_(i in A)p_i(w_i t-theta_i)                   (56)
```

is a finite Laurent polynomial in

```text
z_i=exp(-2pi i theta_i).                             (57)
```

The representation `(m_i)_(i in A)` contributes one distinct monomial,
with nonzero coefficient by (55). Hence the Laurent polynomial is not
identically zero. The same argument applies to block `B` at `-K`. The
product of the two coefficient polynomials is nonzero, so its zero set has
empty interior in the phase torus. The phases may therefore be chosen
arbitrarily close to zero so that both coefficients are nonzero.

Shrinking and translating the intervals changes their `L^1` errors by an
arbitrarily small amount. Consequently:

> Every prescribed finite relation can be made spectrally visible in both
> blocks of arbitrarily accurate positive safe-interval approximants.

This lemma supplies nonvanishing, not a uniform coefficient lower bound,
sign, marked correlation, or cover contradiction. It removes the
seven-unit restriction in (50) only at the level of visibility; it does not
promote the `1372` shapes to discharged rows.

## 8. Power spectra do not name a lift

The finite supports

```text
A={0,1,2,6,8,11},
B={0,1,6,7,9,11}                                    (58)
```

have the same oriented difference multiset. Therefore, for every real `y`,

```text
|sum_(h in A)exp(2pi i h y)|^2
 =|sum_(h in B)exp(2pi i h y)|^2.                    (59)
```

Yet `A` and `B` are neither translates nor reflected translates. Thus even
the full pointwise power spectrum can fail to identify which integer lifts
are present. A magnitude-only root or gap energy cannot canonically select
the atom in (22); it must retain complex phase, a pinned endpoint/address,
or a separate intertwiner.

The one-sheet interval from THM-2284 gives the complementary local no-go:
every nonzero root character can fire while one prescribed integer lift
vanishes. Equation (59) shows that supplying the entire magnitude function
still need not recover the lost support.

## 9. Connection and loss ledger

```text
source:
  THM-2269's every-residue marked energy and the finite endpoint
  complexity of its exclusive-owner set;

target:
  one bounded exact integer Fourier atom in every nonzero residue,
  plus a sharp conditional landing theorem for small owner multipliers;

map:
  differentiate interval indicators, read each residue progression as
  a finite exponential sum, and use a consecutive-zero Vandermonde gate;

preserved:
  residue, exact integer lift, selected owner set, endpoint rank,
  and an explicit finite bound;

destroyed:
  which endpoint supplies the atom, its complex phase/sign, relation
  alignment, current blocker service, and post-expiration ancestry;

cheapest hostile probes:
  coincident endpoint nodes, the d=7 interval (40), the homometric pair
  (58), and a THM-2266 shape whose marked owner is not c_1;

needed sidecar:
  align one bounded lift with a relation coefficient and the same
  ancestry/cut transition, or retain enough complex endpoint phase to
  control the signed whole residue class.                         (60)
```

## 10. Exact referee and scope

The companion uses exact integer and `Fraction` arithmetic. It:

1. verifies the `2R` centered lift windows and sharp `mR-1` arithmetic for
   moduli `13` and `7`;
2. reconstructs representative consecutive-power Vandermonde determinants
   with negative, zero, and positive starting exponents;
3. checks the complete oriented difference multiset in (58)--(59);
4. verifies the proper-arc/complement lengths for `d=1,...,6` and the sharp
   hostile containment from `d=7` onward;
5. replays the complete `3643/1822` THM-2266 atlas and the counts
   (49)--(51);
6. exhausts cyclic component maxima through `S=8`; and
7. reconstructs THM-2199's speed ceiling and the exact bounds (24)--(25).

The universal distributional, Vandermonde, support, and phase-separation
arguments are proved above rather than inferred from the finite checks.
Normal and optimized executions reproduce the stored output byte for byte.

The theorem proves exact lift existence and a sharp conditional pair
landing. It does not choose the same lift as THM-2276/2284, force the
selected owner to be `c_1`, price an owner switch, exclude any of the `165`
profiles, or prove LRC(14). QED.
