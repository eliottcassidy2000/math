---
id: THM-3224
title: "Complete LRC orbit Bernoulli gcd-carry and owner Hodge splitting"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The complete reflected-cell orbit has an exact Bernoulli-quadratic
  evaluation at every dilation, a nonzero residue-periodic g^-2 coefficient
  with an exact O(g^-3) remainder, and the sharp coefficient bound 9/49.
  Cellwise, this curvature is the harmonic part of a canonical cyclic Hodge
  split; the exact part necessarily retains residue and cell-owner data.
  This proves neither cellwise safety, physical entry, nor LRC(14).
source: root/frontier-sidecars-cont-2026-08-02
audit: >
  A self-contained exact companion checks 92,352 carry states through Q=25,
  the complete 18-lane Bernoulli owner table, four direct closed-orbit
  controls, a same-orbit arithmetic sign switch, two infinite cell-branch
  certificates, a zero-first-correction hostile, and a FINITE-EXACT
  168-owner scout.  The scout has no stabilization bound and is not a proof
  dependency.  Normal, optimized, and stored transcripts agree byte-for-byte
  and contain no assertion-dependent gate.  An independent audit rederived
  the Fourier normalization, carry periodicity, sharp range, nonvanishing,
  sign switch, and general Hodge identities.  It identified and repaired the
  finite-scout stabilization overclaim and added an explicit exact-remainder
  gate.
depends_on:
  - THM-3211-uniform-lrc-channel-limit-bernoulli-cubic-and-sharp-floor
  - THM-3200-fixed-lrc-channel-cleared-overlap-quasipolynomial-and-mass-recurrence-boundary
related:
  - THM-3171-global-high-channel-cell90-floor-and-all-width-uniform-two-star-law
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
script: 04-computation/lrc_second_owner_bernoulli_curvature_thm3224.py
output: 05-knowledge/results/lrc_second_owner_bernoulli_curvature_thm3224.out
script_sha256: 84cb3cd5d91b47e4d67918d7eaa49854ec525c9791694cfb679915e4947c7eaa
output_sha256: 413a487add6d53ea5f62734a13c40ccb966a01c64844de4d6fc03eff80a1fea7
hash_basis: LF-normalized bytes
---

# THM-3224 -- complete LRC orbit Bernoulli gcd-carry and owner Hodge splitting

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Universe

Let

```text
L0=168,  H=(1,2,3,4,6,12),
E={{1,2},{1,3},{1,4},{1,6},{1,12},
             {2,3},{2,4},{2,6},{2,12}}.
```

Fix coprime positive integers

```text
P<Q<=2P,  P+Q>=8,
```

and either orientation `(e,f)` of an edge in `E`.  Put

```text
d=gcd(L0,e,f),       T=L0/d,
a=e/d,               b=f/d,
K=|Qe-Pf|/d,
n_g=gPT-a,           m_g=gQT-b,
h_g=gcd(n_g,m_g),    A_g=n_g/h_g,  B_g=m_g/h_g.
```

The inherited nondegeneracy argument gives `K>0`.  For `g>=1`, let

```text
I_(g,j)=integral_[0,1]
 chi(gPu-e(j+u)/L0) chi(gQu-f(j+u)/L0) du,
J_g=sum_(j=0)^(T-1) I_(g,j),
```

where `chi(x)=1` for `||x||<=1/14` and is zero otherwise.

Write

```text
Bbar_2(x)={x}^2-{x}+1/6.
```

## 2. Exact complete-orbit law and second coefficient

For every `g>=1`, define

```text
Xi_g = Bbar_2((A_g-B_g)/14)-Bbar_2((A_g+B_g)/14),
Omega_g = h_g^2 Xi_g/(PQT).
```

Then the complete reflected orbit has the exact evaluation

```text
J_g = T/49 + T h_g^2 Xi_g/(n_g m_g).                 (1)
```

Moreover `h_g|K`, `gcd(h_g,14)=1`, and, if `u_g` is the inverse of
`h_g` modulo `14`, then

```text
Xi_g = Bbar_2((a-b)u_g/14)-Bbar_2((a+b)u_g/14).       (2)
```

Both `h_g` and `Omega_g` are periodic in `g` with period dividing `K`.
The exact second-order expansion, including its remainder, is

```text
J_g = T/49 + Omega_g/g^2 + R_g,                       (3)

R_g = h_g^2 Xi_g [T(Pb+Qa)g-ab]
      /[PQT g^2 n_g m_g].                             (4)
```

Thus `Omega_g` is the exact residue-periodic `g^-2` coefficient and
`R_g=O_(P,Q,e,f)(g^-3)`.

Across the stated eighteen ordered lanes, the complete owner table is

```text
Xi_g in {-8,-6,-5,-3,1,2,3,4,5,8,9}/49.              (5)
```

Consequently the old two-sided bound sharpens to

```text
|J_g-T/49| <= 9T h_g^2/(49n_gm_g)
            <= 9T K^2/(49n_gm_g).                     (6)
```

The coefficient `9/49` is sharp: `(P,Q;e,f;g)=(4,5;1,6;11)` has
`h_g=19`, `Xi_g=9/49`, and `Omega_g=1083/54880`.

There is no zero case in this universe.  In particular,

```text
J_g != T/49 and Omega_g != 0                           (7)
```

for every admissible channel, ordered lane, and `g>=1`.  Hence complete
cell-orbit cancellation stops sharply after the `g^-1` term.

The sign is not determined by the phase orbit.  On the single orbit
`(P,Q;e,f)=(3,5;1,12)`, one has

```text
g not congruent 4 (mod 31): h_g=1,  Xi_g=-5/49,
                             Omega_g=-1/24696;
g congruent 4 (mod 31):     h_g=31, Xi_g= 8/49,
                             Omega_g=961/15435.          (8)
```

Thus the same channel and phase path approaches the complete bulk from
opposite sides according to its arithmetic carry.

## 3. Cellwise coefficient and the endpoint-only obstruction

For a fixed cell owner `j`, let

```text
D_g=(L0gP-e)(L0gQ-f)=d_2g^2+d_1g+d_0,
d_2=L0^2PQ,  d_1=-L0(Pf+Qe),  d_0=ef,
N_(g,j)=D_g I_(g,j),
M=|Qe-Pf|.
```

THM-3200 and THM-3211 give, on each residue `r mod M` and for all
sufficiently large `g` in that residue,

```text
N_(g,j)=d_2L_jg^2+(d_2c_j+d_1L_j)g+kappa_(j,r),        (9)
```

where `L_j` is the cell bulk and `c_j` is the `Bbar_2` endpoint
coboundary.  Define

```text
q_(j,r)=[kappa_(j,r)-d_0L_j-d_1c_j]/d_2.               (10)
```

Then `q_(j,r)` is periodic in `r` with period dividing `M`, and the exact
tail identity is

```text
I_(g,j)=L_j+c_j/g+q_(j,r)/g^2
 -[(d_0c_j+d_1q_(j,r))g+d_0q_(j,r)]/(g^2D_g).          (11)
```

In particular, the last term is `O(g^-3)`.

There is no functional of only the two static phase endpoints that can
produce `q`.  The minimal-channel witness is the fixed cell
`(P,Q;e,f;j)=(3,5;6,1;90)`.  Its endpoints, `L=17/1680`, and
`c=-8213/1411200` are fixed, while

```text
g mod 3 = 1: q=-343771/395136000,
g mod 3 = 2: q=  48229/395136000,
g mod 3 = 0: q= -10057/131712000.                       (12)
```

Thus the missing coordinate at second order is a residue/carry owner, not
another phase-only Bernoulli endpoint.

Vanishing of the first correction does not remove this coordinate.  For
`(P,Q;e,f;j)=(5,6;1,12;56)`, one has `L=4/189`, `c=0`, and, for every
`g>=1`,

```text
N_(g,j)=17920g^2-(704/3)g+kappa_(g mod 9),              (13)
```

with `kappa` for residues `(1,2,3,4,5,6,7,8,0)` equal to

```text
(-64/3, 2188/3, 920, 1076/3, 1312/3,
 -44,    1820/3, 3460/3, 0).                            (14)
```

The corresponding `q` word is

```text
(-17/666792, 11483/13335840, 7243/6667920,
 1129/2667168, 1721/3333960, -697/13335840,
 9551/13335840, 18161/13335840, -1/3333960).            (15)
```

It is nonzero and has both signs.

## 4. Exact owner-dependent Hodge decomposition

Fix a residue `r mod M`.  Reflection gives

```text
q_(T-1-j,r)=q_(j,r).                                    (16)
```

Since summation of the cell expansions equals the complete-orbit expansion,

```text
sum_(j=0)^(T-1) q_(j,r)=Omega_r,                         (17)
```

where the `Xi` entering `Omega_r` is evaluated by `(2)` on that residue (and
depends only on `r mod K`).  Define the canonical owner potential

```text
D_(0,r)=0,
D_(j,r)=sum_(k=0)^(j-1) [q_(k,r)-Omega_r/T].             (18)
```

Then

```text
q_(j,r)=Omega_r/T + D_(j+1,r)-D_(j,r),                  (19)
D_(T,r)=0,               D_(T-j,r)=-D_(j,r).            (20)
```

Equation `(19)` is the repaired endpoint-owner law.  The first correction is
a pure phase coboundary.  The second correction has a nonzero harmonic part
`Omega_r/T`, namely the Bernoulli gcd-carry curvature, plus a residue- and
owner-dependent exact part.  The nonzero holonomy `(17)` is the unique cyclic
obstruction to a pure second potential.

For the control `(3,5;1,2)`, a **FINITE-EXACT scout** fitted from the five
dilations `80,81,82,83,160` has `156` positive entries, `12` negative
entries, and no zero entries, but sums to

```text
Omega=1/24696,                                           (21)
```

not zero.  No stabilization bound proves that those five samples already lie
on every eventual owner branch.  Thus `(21)` is a hostile finite control, not
an all-`g` sign-word theorem and not an input to `(17)--(20)`.

## 5. Proof

THM-3211 gives the exact concatenation

```text
J_g=T integral_[0,1] chi(n_gu)chi(m_gu)du.
```

After dividing by `h_g`, the map `u -> h_gu` preserves Haar measure.  Since
`gcd(A_g,B_g)=1`, the Fourier relation in the product is exactly
`(k,l)=(B_gt,-A_gt)`.  The resulting series is absolutely convergent:

```text
integral chi(Au)chi(Bu)du
=1/49 + 2 sum_(t>=1)
 [sin(pi Bt/7)/(pi Bt)] [sin(pi At/7)/(pi At)]
=1/49 + [Bbar_2((A-B)/14)-Bbar_2((A+B)/14)]/(AB).
```

This proves `(1)`.  The determinant identity

```text
Qn_g-Pm_g=Pb-Qa
```

shows `h_g|K`.  Also `gcd(T,a,b)=1`; hence any common divisor of `h_g` and
`T` would divide `T,a,b`, so `gcd(h_g,T)=1`.  Since `14|T`, `h_g` is a unit
modulo `14`, and division of `n_g=-a mod 14`, `m_g=-b mod 14` proves `(2)`.
Finally `h_g=gcd(n_g,m_g,K)`, so replacing `g` by `g+K` leaves it unchanged.
This proves periodicity.  Subtracting `Omega_g/g^2` from `(1)` gives `(4)`
by one denominator calculation.

For integers `A,B`, `Bbar_2((A-B)/14)=Bbar_2((A+B)/14)` exactly when the
folded residues `|A-B|_14` and `|A+B|_14` agree.  This is equivalent to
`A-B` being congruent to either sign of `A+B` modulo `14`, hence to
`7|A` or `7|B`.  Here `7` divides neither `a` nor `b`, and `h_g` is a unit,
so `(7)` follows.  The finite table of eighteen normalized lanes and six
units modulo `14` is exactly `(5)`; it proves the sharp bound `(6)`.

Equations `(9)--(11)` are direct division of the proved quadratic numerator
by its quadratic denominator.  Reflection `u -> 1-u` proves `(16)` exactly.
Summing the finite set of cell expansions and comparing with `(3)` proves
`(17)`, after which `(18)--(20)` are the canonical cyclic Hodge decomposition.

## 6. Exact audit and scope

Run

```text
python3 04-computation/lrc_second_owner_bernoulli_curvature_thm3224.py
python3 -O 04-computation/lrc_second_owner_bernoulli_curvature_thm3224.py
```

and compare LF-normalized bytes with the declared output.  The self-contained
companion uses exact integer and rational arithmetic.  It checks `92,352`
complete-orbit carry states through `Q<=25`, the complete lane table, four
direct closed-orbit controls, the arithmetic sign switch, both infinite cell
branch certificates, the zero-first-correction hostile, and one full
`168`-owner **FINITE-EXACT scout**.  The scout is explicitly excluded from
the proof graph because it has no branch-stabilization certificate.  Normal,
optimized, and stored transcripts agree byte-for-byte.

The theorem distinguishes four scopes:

- `q_(j,r)` is cellwise and retains the cell owner;
- the complete-orbit formula is exact for all finite `g` but forgets the
  owner distribution;
- an asymptotic coefficient does not replace a cell's pre-stabilization head;
- no statement supplies physical survivor entry, closes the rung, proves
  cellwise safety, or proves `LRC(14)`.

**QED.**
