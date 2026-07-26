---
id: THM-2367
title: "Septimal root averaging, graft drift, and cover alignment"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. For a
  lawful isolated blocker/graft tensor with 13|c and 13 not dividing q,
  the overlap is circulant exactly when
  nu_7(c)>nu_7(q); the statement holds for either an ordinary unit-safe
  graft or the guard-safe graft. This local drift is not monotone under
  multiplication: an exact 90/91-mass Boolean mask restores circulancy
  to the escaping pair (c,q)=(13,7), and a full nine-factor bare-owner
  packet is circulant on every valuation profile (1,b,c), though that
  packet has positive uncovered mass and is not a scalar cover.
  Conversely, the scalar cover supplies additive order-seven orbit
  rigidity. In a mixed guard/unit septimal layer, some blocker has
  nu_7(c_j)<=max(nu_7(H),nu_7(q_i)); if only c_3 is dominant, a sharp
  seven-bin inequality routes a low blocker into c_1 or c_2. In the
  hard W=k=1 lane the six lower factors have a mean-zero signed defect
  supported on only three absorber combs; its endpoint current gives an
  exact handoff congruence and a global event-capacity obstruction. This
  excludes the displayed local N=49 chamber from extending to a scalar
  cover. The mechanism is additive Phi_7 rigidity, not the
  multiplicative chi_7/Fano structure of THM-1156. It does not promote
  isolated drift through the full owner, land a target, exclude a
  scalar row, or prove LRC(14).
source: codex-2026-07-25-septimal-root-averaging
depends_on:
  - THM-2296-prescribed-expiration-return-or-bounded-ancestry-resonance
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
related:
  - THM-1156-tooth-seam-chi7-bipartition
  - THM-2362-thirteen-shift-successor-statistic-and-role-jet-floor
  - THM-2364-anchored-corner-forces-mixed-deep-blocker-colour
script: 04-computation/lrc14_septimal_root_graft_cover_thm2367.py
output: 05-knowledge/results/lrc14_septimal_root_graft_cover_thm2367.out
script_sha256: 269e87cd870b59c00a3889212c4520d4b568b47463c079029a7cce0b1501a58c
output_sha256: aa13c7bd1a4c25d6d2858c6639113129752fdc44f58092514bceee7c41ddc85b
hash_basis: working-tree bytes (LF)
---

# THM-2367 -- the H-drift obstruction has a septimal role

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2365 reduces target landing to departure of a nonnegative overlap
tensor from the circulant law `H(r,s,t)=G(r-t)`. This theorem identifies
the first arithmetic mechanism which creates or preserves that law and
then brings the scalar cover back into the picture.

The resulting chain is

```text
isolated lawful graft:
  drift iff the blocker is not strictly 7-deeper than the graft;

arbitrary extra mask:
  isolated drift can be cancelled exactly;

full scalar cover:
  some blocker nevertheless has an arithmetically escaping graft;

remaining debt:
  align that blocker with an eligible owner/target role and prevent
  the other owner factors from cancelling its drift.               (1)
```

## 1. The isolated lawful graft

Put

```text
p=13,

d(y)=1_(||y||<1/14),

g(y)=1-d(y),

u_L(y)=1_(||y||>=L/14),             L in {1,2}.     (2)
```

Thus `u_1=g`, while `u_2` is the guard-safe base factor. Boundary
conventions are immaterial below.

Let `c,q` be positive integers satisfying

```text
13|c,                       13 does not divide q.   (3)
```

For `r,t in F_13`, define

```text
H_(c,q,L)(r,t)
 =integral_T
    d(c x-r/13)
    g(c x-t/13)
    u_L(q x+t/13) dx.                              (4)
```

The signs in (4) are the lawful target-dual signs for

```text
ell=e_c-e_q:
```

the deep complement moves by `-t/13`, while the grafted safe factor
moves by `+t/13`.

Write

```text
h=gcd(c,q),              C=c/h,              Q=q/h. (5)
```

Then

```text
13|C,                    13 does not divide Q,
gcd(C,Q)=1.                                          (6)
```

The exact criterion is

```text
H_(c,q,L)(r,t)=Phi(r-t) for some Phi

 iff 7|C

 iff nu_7(c)>nu_7(q).                               (7)
```

It holds for both `L=1` and `L=2`.

## 2. Root averaging proves the criterion

Haar reduction by `x -> h x`, followed by

```text
z=Cx-t/13,
```

gives

```text
H_(c,q,L)(r,t)
 =integral_T
   f_(r-t)(z) A_(C,L)(Qz+Qt/13) dz,                (8)

f_a(z)=d(z-a/13)g(z),

A_(C,L)(Y)
 =1/C sum_(k=0)^(C-1)u_L((Y+k)/C).                 (9)
```

The apparently omitted term `Ct/13` is an integer. Multiplication by
`Q` permutes the `C` roots because of (6).

Root averaging selects the Fourier coefficients

```text
A_(C,L)^hat(n)=u_L^hat(Cn).
```

For `n!=0`,

```text
u_L^hat(Cn)
 =-sin(pi L C n/7)/(pi Cn).                        (10)
```

Since `L` is a seven-unit, all nonconstant coefficients in (10) vanish
exactly when `7|C`. In that case

```text
A_(C,L)(Y)=1-L/7
```

almost everywhere, and (8) depends only on `r-t`. With strict-open
intervals there is one aligned endpoint phase at which the finite root
count differs; this null exception does not change any integral.

Conversely suppose `7` does not divide `C`. Put

```text
s=L C mod 14,

k=min(s,14-s) in {1,...,6},

delta=k/14,

epsilon=+1  if s in {1,...,6},
epsilon=-1  if s in {8,...,13}.                    (11)
```

Exact Fourier inversion of (10) gives, almost everywhere,

```text
A_(C,L)(Y)
 =1-L/7
  -(epsilon/C)(d_delta(Y)-2delta),                 (12)

d_delta(Y)=1_(||Y||<delta).
```

On the slice `r=t+1`,

```text
f_1=1_I,

I=(1/14,1/14+1/13)                                 (13)
```

almost everywhere. Thus nonconstancy in `t` reduces to

```text
M_t
 =integral_(I+t/13)d_delta(Qy)dy.                  (14)
```

Write

```text
Q=13A+rho,                    1<=rho<=12.
```

Partition one turn of the centred `d_delta` arc by the offset
thirteen-grid into cell masses `x_0,...,x_12`. Then

```text
M_t
 =1/Q(
   2delta A+S_(rho t)
  ),

S_j=sum_(i=0)^(rho-1)x_(j+i).                      (15)
```

If all `S_j` were equal, then

```text
x_(j+rho)=x_j.
```

Since `rho` is a unit modulo thirteen, all `x_j` would be equal to
`2delta/13`. That is impossible: one of the arc or its complement has
length at least one half and contains a complete grid cell, giving some
cell mass `1/13` or `0`, while

```text
0<2delta/13<1/13.
```

Therefore (14) is nonconstant, proving the reverse direction of (7).
This one `r=t+1` toothpick slice is already a decisive drift witness.

## 3. Isolated drift is not monotone

The condition in (7) is not enough after multiplication by an arbitrary
fixed nonnegative mask. There is an exact Boolean hostile.

Take

```text
(c,q,L)=(13,7,1),

D=16562,
```

and let `W` be one on every half-open cell

```text
[j/D,(j+1)/D)
```

except the ten linear ranges

```text
[16555,16562), [0,13),
[1625,1651),   [2457,2463),
[3263,3289),   [4907,4927),
[7449,7455),   [9087,9113),
[10725,10751), [12363,12389).                     (16)
```

The first two ranges form one circular interval. Exactly `182` cells
are deleted, so

```text
mu(W)=16380/16562=90/91.                            (17)
```

Define `H_W` by inserting `W(x)` in the integrand (4). Direct exact cell
counting gives

```text
H_W(r,t)=G(r-t),

G(0)=0,

G(+1)=G(-1)=1078/16562=11/169,

G(a)=2002/16562=11/91       for a notin {0,+1,-1}. (18)
```

Thus a target-independent Boolean mask of mass `90/91` restores exact
circulancy to an isolated pair which drifts when `W=1`. Positivity,
large mask mass, and isolated-graft drift do not give a monotone route to
the full owner tensor.

## 4. Septimal peeling and a full-owner hostile

The following elementary lemma explains a larger hostile. If

```text
phi in {d,g,u_2},

7|N,
```

and `K` is integrable, then

```text
integral_T phi(y+alpha)K(Ny)dy
 =mu(phi)integral_T K.                              (19)
```

Indeed, every nonzero Fourier frequency of `K(Ny)` is a nonzero multiple
of seven, where all three displayed base factors have zero Fourier
coefficient.

Now put

```text
H=1,

(q_1,q_2,q_3,q_4,q_5)=(420,2940,2,3,5),

c_1=13*7*q_2,

c_2=c_1*7*13^(b-1),

c_3=c_2*7*13^(c-b),              1<=b<c.            (20)
```

These speeds have the blocker valuation profile

```text
(nu_13(c_1),nu_13(c_2),nu_13(c_3))=(1,b,c).        (21)
```

Use owner `j=1`, omit `H` in the THM-2309 packet, and use the lawful
target-dual basis

```text
eta=e_(c_2)-e_(q_1),

ell=e_(c_3)-e_(q_2).                               (22)
```

Let `H_E(r,s,t)` be THM-2365's full bare-owner tensor for this packet.
The low-frequency factor

```text
F_0(x)
 =u_2(x)g(2x)g(3x)g(5x)
```

is constant on the common `1/420` endpoint cells, and exact enumeration
gives

```text
mu(F_0)=182/420=13/30.                              (23)
```

After the `420`-root average, repeatedly apply (19) along

```text
q_2/q_1=7,

c_1/q_2=7*13,

c_2/c_1=7*13^(b-1),

c_3/c_2=7*13^(c-b).                                (24)
```

The factors at `q_1,q_2,c_1,c_2` peel with means

```text
6/7, 6/7, 1/7, 6/7.
```

Consequently

```text
H_E(r,s,t)=A J(r-t),

A=(13/30)(6/7)^3(1/7)=468/12005,                   (25)

J(0)=0,

J(+1)=J(-1)=1/13,

J(a)=1/7                    otherwise.              (26)
```

This is a full nine-factor, lawfully shifted, positive bare owner with
exactly zero H-drift for every `1<=b<c`, including the repeated-first
case `b=1`.

It is deliberately not a scalar-cover row. Its uncovered set has exact
mass

```text
(13/30)(6/7)^5
 =16848/84035
 >0.                                                (27)
```

Thus valuation type, packet rank, every owner factor, and lawful target
shifts still do not exclude the circulant branch. A valid exclusion must
use the global scalar cover or an equivalent nonlacunarity input.

## 5. The scalar cover forces top-layer alignment

Return to a genuine scalar cover. Put

```text
U={H,q_1,...,q_5},

M=max_(u in U)nu_7(u),

N=7^(M+1).                                         (28)
```

First suppose

```text
nu_7(c_j)>M             for j=1,2,3.                (29)
```

The blocker-safe set

```text
S=intersection_(j=1)^3 D_(c_j)^c
```

has measure at least `4/7` and is invariant under translation by `1/N`.
On almost every complete `N`-orbit in `S`, the scalar cover reduces to

```text
E_H union D_(q_1) union ... union D_(q_5),          (30)

E_H={x:||Hx||<1/7}.
```

Each such orbit has exactly `2N/7` guard incidences and `N/7` incidences
for each unit danger. Their total capacity is exactly `N`. Hence (30) is
an exact one-fold partition of every orbit point.

Take the six nontrivial orbit characters

```text
frequency=k*7^M,                   k=1,...,6.
```

A label `u` with `nu_7(u)<M` has zero transform there; a label in the
top layer `nu_7(u)=M` gives, after reduction to `Z/7`, one occupied
residue for a unit danger and two consecutive occupied residues for the
guard danger. Since the partition is constant one, the aggregate
top-layer mask has all six nonconstant transforms zero and is therefore
constant on `Z/7`.

Its total arc weight is

```text
W
 =2*1_(nu_7(H)=M)
  +#{i:nu_7(q_i)=M}.                               (31)
```

The top layer is nonempty, so `1<=W<=7`. Constancy forces `7|W`, hence

```text
W=7.
```

Therefore every one of the six labels in `U` has the same septimal
valuation `M`. In a primitive scalar row that common value is zero.

Contrapositively:

> If the guard/unit septimal valuations are not all equal, then
>
> ```text
> min_j nu_7(c_j)<=max_(u in U)nu_7(u).             (32)
> ```

This is the first global-cover input into the H-drift graft problem.

## 6. The only-c3-dominant seven-bin inequality

The live residual needs a more asymmetric statement. Assume

```text
nu_7(c_3)>M,

k=#{j in {1,2}:nu_7(c_j)<=M}.                      (33)
```

Restrict to the positive, `N`-orbit-invariant set on which every blocker
of valuation greater than `M` is safe. The orbit cover now uses the
guard, five unit dangers, and the `k` low blockers. Measured in
`1/7`-arc units, its total capacity is

```text
P=7+k.
```

Collapse the orbit into its seven top-layer residue bins. Let `W` be the
total top-layer arc weight among the guard, units, and low blockers, and
let `m_r` be the number of top-layer arcs occupying bin `r`. Sub-top
factors contribute uniformly

```text
(P-W)N/49
```

to each bin, while top factors contribute `m_r N/7`. Coverage of all
`N/7` sites in the bin gives

```text
7m_r>=W-k                       for every r.         (34)
```

If `W>k`, every `m_r` is positive and `W=sum_r m_r>=7`. Therefore

```text
W<=k        or        W>=7.                         (35)
```

There is a useful profile-free classification in the subcase

```text
M>0,                         nu_7(H)<M.             (36)
```

Let

```text
t=#{i:nu_7(q_i)=M},

b=#{j in {1,2}:nu_7(c_j)=M}.
```

Then `1<=t<=5`, `0<=b<=k`, and `W=t+b`. Equations (33)--(35) give the
complete list:

- `k=0` is impossible;
- if `k=1`, then `(t,b)=(1,0)`;
- if `k=2`, then

  ```text
  (t,b) in {(1,0),(1,1),(2,0),(5,2)}.              (37)
  ```

The last alternative is the exceptional saturated case `W=7`; the
other alternatives have `W<=k`. In particular, the `k=1` hard lane has
exactly one top unit and its unique low blocker lies strictly below
`M`. This conclusion no longer imports an unproved count of
seven-divisible `q_i` from a differently labelled branch.

The slack from low blockers is exactly why the perfect-tiling conclusion
of Section 5 cannot simply be reused when only `c_3` is dominant.

### The hard role lane is locally sharp

Positive owner mass on one orbit does not strengthen the `k=1`
conclusion. Take

```text
N=49,                  x_0=1/686,

H=1,

(q_1,q_2,q_3,q_4,q_5)=(7,148,171,172,243),

(c_1,c_2,c_3)=(195,16562,215306).                  (38)
```

The blocker thirteen-adic profile is `(1,2,3)`. Septimally, `q_1` has
valuation one, `c_2,c_3` have valuation two, and every other displayed
speed has valuation zero. Thus `M=1`, `c_3` is dominant, and the unique
low blocker is the present source `c_1`, with `W=k=1`.

On the orbit

```text
O_(x_0)={x_0+j/49:j in Z/49},
```

the strict masks are

```text
E_H:     {0,...,6,42,...,48},

D_(c_1):{11,...,17},

D_(c_2)=D_(c_3)=empty,

D_(q_1):{0,7,14,21,28,35,42},

D_(q_2):{35,...,41},

D_(q_3):{18,20,22,24,26,28,30},

D_(q_4):{19,21,23,25,27,29,31},

D_(q_5):{7,8,9,10,32,33,34}.                      (39)
```

The guard, `c_1`, and the four lower-unit masks are pairwise disjoint
and tile all `49` sites. The top unit adds seven double-covered sites.
The six sites

```text
{11,12,13,15,16,17}
```

are genuine exclusive `c_1`-owner sites. Every pattern in (39) persists
for

```text
|x-x_0|<1/3014284,                                 (40)
```

so this is a positive-measure chamber, not an endpoint coincidence.
It is not asserted to extend to a global scalar cover. It proves that
one-orbit capacity and positive exclusive-owner conditioning cannot
remove the strict-row source-role equality lane; cross-orbit or chamber
coherence is genuinely needed.

## 7. The signed event current is the cross-chamber invariant

The `k=W=1` lane has more global structure than one orbit sees. Retain
(36), call the unique top unit `q_*`, call the unique low blocker `c_*`,
and call the other two blockers `a,b`. Put

```text
A=D_(q_*) union D_a union D_b.
```

The seven-bin count used above applies on every orbit on which `a,b`
are safe. In each of the six bins outside `D_(q_*)`, the remaining
capacity equals the number of sites. Covering therefore forces an exact
one-fold partition there. Consequently the lower multiplicity

```text
L(x)
 =1_(E_H)(x)
  +sum_(q_i!=q_*)1_(D_(q_i))(x)
  +1_(D_(c_*))(x)
```

satisfies

```text
F=L-1=0                    almost everywhere on A^c. (41)
```

Thus the nonlinear nine-factor cover has become a nonzero-capable,
mean-zero integer step function on six lower frequency lattices whose
physical support lies in three absorber combs.

There is a convention-free endpoint form. Index the six lower factors
by `(v_i,a_i)`, where `(v_H,a_H)=(H,2)` and the five ordinary danger
factors have `a_i=1`. The distributional derivative of `F` is

```text
mu_L
 =sum_i sum_(n=0)^(v_i-1)
   (delta_((n-a_i/14)/v_i)-delta_((n+a_i/14)/v_i)).
                                                               (42)
```

Equation (41) says precisely

```text
mu_L restricted to T minus closure(A) is zero.       (43)
```

This repairs a tempting but false stronger claim: an individual lower
wall need not lie in an absorber, because an exit wall may hand off to
an entry wall at the same null point. What is forced is signed
cancellation of the entire event cluster.

The handoff arithmetic is exact. An exit of lower factor `i` coincides
with an entry of lower factor `j` if and only if

```text
14 gcd(v_i,v_j) divides a_i v_j+a_j v_i.             (44)
```

When (44) holds there are exactly `gcd(v_i,v_j)` such coincidences.
For two ordinary dangers this is

```text
14 gcd(v_i,v_j) divides v_i+v_j;
```

for a guard/danger pair it is the corresponding `H+2v` congruence.
These opposing-wall incidences are an intrinsic directed event
relation; no cosmetic tournament orientation is being imposed.

This gives a cheap global obstruction. For ordered `i!=j`, put

```text
h_(ij)
 =gcd(v_i,v_j) if (44) holds,
 =0             otherwise.
```

For `w in {q_*,a,b}`, put

```text
g_(iw)=gcd(v_i,w),              C_(iw)=v_i/g_(iw).
```

An entry branch of `v_i` samples a closed absorber arc of length `1/7`
on a `C_(iw)`-grid, with multiplicity `g_(iw)`. It therefore contributes
at most

```text
g_(iw)(floor(C_(iw)/7)+1)
```

absorbed entry atoms. Every remaining entry belongs to a balanced event
cluster and can be paired with an exit. The union bound yields the
necessary event-capacity inequality

```text
sum_i v_i
 <=sum_(i!=j)h_(ij)
   +sum_i sum_(w in {q_*,a,b})
      g_(iw)(floor(C_(iw)/7)+1).                    (45)
```

For the local chamber (38), the six lower speeds are

```text
(1,148,171,172,243,195).
```

Their total entry mass is `930`; no ordered handoff congruence (44)
holds; and the absorber term for `(7,16562,215306)` is only `432`.
Thus (45) fails:

```text
930>432.
```

The positive chamber in Section 6 cannot extend to a global scalar
cover. This is a genuine cross-chamber gain over the one-orbit count.

The event screen is not yet a row exclusion. For the arithmetically
lawful control

```text
H=1,

(q_*,q_2,q_3,q_4,q_5)=(7,2,3,4,5),

c_*=13,

a=98*13^2*60,                   b=13a,              (46)
```

the thirteen-adic blocker profile is `(1,2,3)`, the same hard septimal
roles hold, and (45) has left side `28`, no handoffs, and absorber
capacity `63`. Nevertheless it is not a scalar cover: `x=319/2000`
lies in `C_H` and in none of the eight strict danger combs. Hence
event capacity is necessary, not sufficient; the next invariant must
see values between the event atoms or their signed Fourier phases.

## 8. What is now missing

Sections 1, 5--7 separate arithmetic availability from analytic
transport:

```text
mixed U-layer + scalar cover
  -> some blocker c has nu_7(c)<=M
  -> some q in U gives a locally drifting lawful graft for c.        (47)
```

THM-2309 may use any of the six `U` labels as a graft; Section 1 includes
both the ordinary unit-safe and guard-safe base factors.

On a strict first-depth-one row, if `c_2` is the low blocker, then

```text
nu_13(c_1)=1<nu_13(c_2),
```

so THM-2349's shallow-carrier lemma and THM-2365 may be rerun with
`d=c_2`, rather than `c_3`. This aligns the arithmetic escape with an
eligible target role.

The exact unresolved role branches are sharper:

1. a strict row on which `c_1` is the only septimally low blocker, so
   the available escape sits at the present source rather than an
   excluded target; and
2. a repeated-first row, where `c_2` has the same thirteen-adic depth as
   `c_1` and only `c_3` satisfies the strict shallow/deep carrier
   hypothesis.

Even outside these branches, Sections 3--4 show that isolated drift is
not preserved by arbitrary extra owner masks. A conditioning,
deletion, source-shift covariance, or phase-cone lemma must still
transport the local witness into positive full `D_(H_E)`.

The character argument in Sections 5--6 is additive order-seven
(`Phi_7`) orbit rigidity. It is not the multiplicative Legendre
`chi_7` bipartition or Fano-line incidence of THM-1156. Those remain a
distinct possible refinement once a genuine incidence map is supplied.

No scalar profile is excluded. The ledger remains `165`, and LRC(14)
remains open.

## 9. Exact companion

The dependency-free exact companion:

- checks the root-average step formula for both `L=1,2` and all
  fourteen residues of `C`;
- checks all `288` nonconstant thirteen-window controls in the converse;
- verifies the three basic graft examples `(91,1),(91,7),(13,1)`;
- recounts all `16,562` Boolean-mask cells and all `169` shifted
  overlaps in (16)--(18);
- enumerates the `420` low-frequency cells in (23);
- checks the peeling ratios, constants, and every one of the `165`
  valuation profiles in the full-owner hostile;
- exhausts the top-layer weight and only-`c_3` seven-bin alternatives;
- reproduces the `49` masks and positive chamber in (38)--(40);
- verifies the signed event-current handoff and capacity arithmetic,
  including the `930>432` chamber exclusion and the lawful `28<=63`
  noncover control.

Run

```bash
python3 04-computation/lrc14_septimal_root_graft_cover_thm2367.py
python3 -O 04-computation/lrc14_septimal_root_graft_cover_thm2367.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_septimal_root_graft_cover_thm2367.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

Independent audit is pending. QED.
