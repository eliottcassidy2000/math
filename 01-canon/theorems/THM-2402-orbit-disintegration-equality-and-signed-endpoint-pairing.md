---
id: THM-2402
title: "Orbit-disintegration equality and signed endpoint pairing"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING. For
  T_q(x)=qx on the circle, if E lies over H and every H-fibre meets E,
  then q mu(E)-mu(H)=integral_H(N_E-1). Equality is equivalent to an
  almost-everywhere one-sheet measurable section and gives the exact
  Fourier transfer q Ehat(qn)=Hhat(n). For one-dimensional BV sets the
  distributional boundary obeys D N_E=(T_q)_*D1_E, so equality pairs
  every exposed endpoint above the interior of H with opposite signed
  endpoint mass in another root. The multiplicity support B_2 obeys
  X/(q-1)<=mu(B_2)<=X, sharply. Applied to THM-2396, equality in
  delta>=66/4459 is exactly a one-sheet forty-nine-root section; strict
  surplus quantitatively forces two-sheet base mass. Rational comb
  endpoints with the same T_49 image satisfy an exact congruence:
  if both speed labels are seven-units, their root displacement is
  divisible by seven. Thus a primitive interior paired event in the
  common-core packet must involve q_*, conditionally on such an event
  existing. No primitive event, canonical owner, terminal target, row
  exclusion, or LRC(14) proof is claimed.
source: codex-2026-07-26-orbit-equality-endpoint-pairing
depends_on:
  - THM-2396-common-core-forty-nine-orbit-word-incompatibility
related:
  - THM-2393-c3-safe-double-fibre-capacity-and-common-core-residual
  - THM-2399-physical-one-clean-forty-nine-orbit-sharpness
  - THM-2400-clean-parent-root-gauge-quotient-and-target-slope-boundary
script: 04-computation/lrc14_orbit_equality_endpoint_pairing_thm2402.py
output: 05-knowledge/results/lrc14_orbit_equality_endpoint_pairing_thm2402.out
script_sha256: 854520c367c075cfa52a5b56592f14c2ae189849cc94533f7f620259c659c975
output_sha256: f9c642cbf3236a63dfd2a926aa5e2e0ea239e072efd47deb8e79de469f90841f
hash_basis: working-tree bytes (LF)
---

# THM-2402 -- orbit equality is a one-sheet endpoint-paired section

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-2396 obtains the common-core clean-mass floor by integrating the
pointwise statement that every high-safe forty-nine-root orbit contains
at least one clean root. THM-2399 proves that the pointwise multiplicity
one is physically sharp on a single orbit, but its displayed scalar row
does not cover the circle and does not attain the integrated floor.

The present theorem identifies the actual equality object. Equality in
an orbit-disintegration floor is not merely "one clean point occurred
somewhere": it is a measurable one-sheet section on almost every base
fibre. In finite interval geometry, every change of selected sheet has
an exact opposite-sign endpoint partner. This exposes a new conditional
septimal route, while also isolating its missing input: no current
theorem forces a primitive interior sheet change.

## 1. The measurable cyclic-cover identity

Let

```text
T_q:R/Z -> R/Z,                 T_q(x)=qx,             q>=2,
```

and let `E,H` be measurable subsets of the circle. Define the root
multiplicity

```text
N_E(z)=sum_(j=0)^(q-1) 1_E((z+j)/q).                 (1)
```

Assume

```text
E subset T_q^(-1)(H)                    almost everywhere,

N_E(z)>=1                               for almost every z in H.  (2)
```

The first inclusion makes `N_E=0` almost everywhere off `H`. Branchwise
change of variables gives

```text
integral N_E(z) dz=q mu(E).                             (3)
```

Consequently

```text
X:=q mu(E)-mu(H)
  =integral_H (N_E(z)-1) dz
  >=0.                                                  (4)
```

This identity needs no interval, rationality, or regularity hypothesis.

## 2. Equality, sections, and Fourier transfer

For `j=0,...,q-1`, put

```text
H_j={z in H:(z+j)/q in E}.                             (5)
```

The following are equivalent:

```text
(i)   q mu(E)=mu(H);

(ii)  N_E=1_H almost everywhere;

(iii) the H_j partition H almost everywhere and
      E is, modulo a null set, the union of the selected branches
      {(z+j)/q:z in H_j}.                              (6)
```

Indeed, (4) has a nonnegative integer-valued integrand. It vanishes
exactly when there is one selected root over almost every point of `H`.
The branch representation is then unique almost everywhere.

Use the Fourier convention

```text
fhat(n)=integral f(x) exp(-2 pi i n x) dx.             (7)
```

Direct substitution on the `q` inverse branches yields, for every
integer `n`,

```text
Nhat_E(n)=q 1_Ehat(qn).                                (8)
```

Hence equality in (4) gives the complete base-frequency identity

```text
q 1_Ehat(qn)=1_Hhat(n)                  for every n in Z.  (9)
```

This determines only the Fourier modes at multiples of `q`; the
nonzero root characters still encode which sheet is selected.

## 3. Strict surplus and sharp multiplicity-support bounds

Let

```text
B_2={z in H:N_E(z)>=2}.                                (10)
```

On `B_2`,

```text
1<=N_E-1<=q-1,
```

while the defect vanishes elsewhere. Therefore

```text
X/(q-1)<=mu(B_2)<=X.                                  (11)
```

Both constants are sharp. For the upper equality, start with any
one-sheet section and add exactly one further sheet over a set of mass
`X`. For the lower equality, add all other `q-1` sheets over a set of
mass `X/(q-1)`.

The exact companion includes rational interval versions of both
boundaries. Thus no improvement of (11) is available from
disintegration alone.

## 4. BV boundary pushforward and signed pairing

Assume now that `1_E` has bounded variation on the circle. Its
distributional derivative is a finite signed measure `D1_E`. For every
smooth periodic test function `phi`,

```text
<D N_E,phi>
 =-integral N_E(z) phi'(z) dz
 =-q integral 1_E(x) phi'(T_q x) dx
 =-integral 1_E(x) (phi o T_q)'(x) dx
 =integral phi(T_q x) d(D1_E)(x).                     (12)
```

Thus

```text
D N_E=(T_q)_*(D1_E).                                  (13)
```

If equality holds in (4) and `1_H` is represented by its BV version,
then

```text
(T_q)_*(D1_E)=D1_H.                                   (14)
```

For a reduced finite union of circular intervals, `D1_E` is the signed
endpoint measure: `+1` at each entry and `-1` at each exit. Let `z_0`
be outside the boundary of `H`. Equation (14) says

```text
sum_(x in partial E, T_q x=z_0) jump_E(x)=0.           (15)
```

Therefore every exposed endpoint above `z_0` has opposite total signed
endpoint mass above another `q`-root of `z_0`. In the reduced interval
case all jumps are `+1` or `-1`, so at least one endpoint of the
opposite sign exists.

This is a **signed** assertion. An unsigned endpoint count, total
variation, or root energy cannot distinguish an entry from an exit.
Rationality is unnecessary for (13)--(15); it only turns the endpoint
bank into a finite exact rational object.

A sharp equality control is already nonconstant. For `q=5`, take

```text
H=[0,1/2),

E=[0,1/20) union [1/4,3/10).                          (16)
```

The selected branch is `j=0` over `[0,1/4)` and `j=1` over
`[1/4,1/2)`. The exit at `x=1/20` and entry at `x=1/4` both map to
`z=1/4` and cancel in (14). Thus equality does not force one globally
constant sheet.

## 5. The THM-2396 common-core equality face

Retain THM-2396's notation

```text
q_*=7u,                C_3=49V,                c_3=13C_3,

H_0=D_V^c intersection D_(13V)^c,

S={K=0} intersection D_(q_*)^c intersection D_(c_3)^c
    intersection
    (D_(C_1) union D_(C_2) union D_(C_3))^c,

delta=mu(S).                                           (17)
```

Because `C_3` and `c_3` are safe on `S`,

```text
S subset T_49^(-1)(H_0)                 almost everywhere.  (18)
```

THM-2396 proves

```text
N_S(z)>=1                              for almost every z in H_0,  (19)
```

and computes

```text
mu(H_0)=1-2/7+1/91=66/91.              (20)
```

Applying (4) gives the exact identity behind its floor:

```text
49delta-66/91
 =integral_(H_0)(N_S-1).                               (21)
```

Hence

```text
delta=66/4459
```

if and only if `S` is a one-sheet measurable section over `H_0`.
In that equality case,

```text
49 1_Shat(49n)=1_(H_0)hat(n)             for all n in Z,  (22)

(T_49)_*(D1_S)=D1_(H_0).                              (23)
```

All sets here are finite unions of rational comb intervals, so (23) is
a finite signed rational endpoint identity after choosing half-open
representatives.

If the floor is strict, put

```text
B_2={z in H_0:N_S(z)>=2},

X=49delta-66/91.                                      (24)
```

Then the sharp quantitative consequence is

```text
X/48<=mu(B_2)<=X.                                    (25)
```

THM-2399 does not instantiate the equality case. It constructs one
strict physical orbit with `N_S(z)=1`, while its scalar row has the
global safe point `3/14`; it is not a circle cover. Its correct role is
to show that (19) cannot be strengthened pointwise from the same local
orbit data.

## 6. Exact endpoint arithmetic at forty-nine roots

Every centered ordinary or guard comb boundary has a lift of the form

```text
x=(14a+epsilon L)/(14v),                epsilon in {+1,-1},

L=1  for an ordinary danger,
L=2  for a guard danger.                               (26)
```

Take a second endpoint

```text
x'=(14b+eta M)/(14w),                  eta in {+1,-1},  (27)
```

with `M in {1,2}`. If the two endpoints have the same `T_49` image,
choose the displayed rational lifts and write

```text
49(x-x')=m in Z.                                      (28)
```

Clearing denominators gives the exact identity

```text
49[w(14a+epsilon L)-v(14b+eta M)]
 =14m v w.                                            (29)
```

If both speed labels are septimal units,

```text
7 does not divide vw,
```

then division of (29) by seven shows

```text
7[w(14a+epsilon L)-v(14b+eta M)]=2m v w,

so                       7 divides m.                 (30)
```

The residue class of `m mod 7` is independent of changing either
endpoint lift by an integer, since that changes `m` by a multiple of
forty-nine.

The divisibility boundary is sharp. With

```text
v=7, w=1, L=M=1, a=b=0, epsilon=eta=-1,
```

the circle representatives are

```text
x=97/98,                  x'=13/14,

49(x-x')=3,                                             (31)
```

so a primitive displacement is possible as soon as one speed is
seven-divisible.

## 7. Conditional common-core label consequence

In the THM-2396 residual typing, the guard, the four lower `q` labels,
and the low common-core blocker labels are septimal units. The sole
remaining lower-layer seven-divisible label is `q_*`; the high pair
`C_3,c_3` has valuation at least two.

An endpoint of `S` above the interior of `H_0` cannot be contributed
solely by `C_3` or `c_3`, because those endpoints push to
`partial H_0`. Assign to each exposed interior endpoint any underlying
comb label whose boundary contains it. Combining (15) and (30) gives:

```text
If an opposite-sign paired interior endpoint event has
49-root displacement m with 7 not dividing m, then at least one
endpoint has q_* among its boundary labels.                         (32)
```

This is the exact surviving endpoint-typing mechanism. It is
conditional and does **not** prove that such a primitive event exists.
A one-sheet section can switch only through displacements divisible by
seven, or can have no interior switch at all. Even when (32) fires, it
does not by itself choose a canonical owner, align the endpoint with
THM-2400's target slope, or preserve a terminal complex phase.

The next sharp service is therefore:

```text
force at least one interior opposite-sign sheet switch with
primitive mod-seven displacement, then transport its q_*-labelled
jump into the canonical owner/target gauge.                         (33)
```

No row is excluded, the scalar ledger remains `165`, and LRC(14)
remains open.

## 8. Exact companion

The dependency-free exact companion:

- verifies a rational nonconstant one-sheet section and its signed
  interior endpoint cancellation;
- realizes both sharp boundaries in (11);
- checks the common-core arithmetic (20)--(21);
- exhausts `705600` endpoint pairs in a finite positive-control box,
  finding zero violations of (30);
- verifies the exact primitive hostile (31); and
- keeps the measure-equality, BV-pairing, and LRC-label scopes separate.

Run

```bash
python3 04-computation/lrc14_orbit_equality_endpoint_pairing_thm2402.py
python3 -O 04-computation/lrc14_orbit_equality_endpoint_pairing_thm2402.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_orbit_equality_endpoint_pairing_thm2402.out
```

after LF normalization. Every executable assertion remains active
under optimized Python.
