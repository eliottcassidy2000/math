---
id: THM-2503
title: "Sharp duplicate-probe active-dipole sweep"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For a
  positive weight supported where one ordinary
  or guard role fails and one helper role is safe, oppositely shifting
  those two probes through C_13 leaves both safe at 9--11 shifts in
  the ordinary case and 7--10 shifts in the guard case. The sharp
  profile censuses are 20,35,14 and 32,69,52,8. Total duplicate-probe
  service lies in [9rho,11rho] or [7rho,10rho], so at least 9 or 7
  shifts are active and one has mass at least the sharp 4rho/5 or
  2rho/3.
  With rational interval data every nonzero shift Fourier colour
  survives; some normalized coefficient has real part at most
  -T/156. The sharp nonzero energy is at least T^2/2028 for an
  ordinary failure and T^2/1690 for a guard failure. This sharpens
  the analytic survivor of THM-2384 but remains its type-incorrect
  duplicate-probe geometry, not a THM-2350 target action or canonical
  relation current. No scalar row is excluded.
source: mac-mini-2026-07-27-active-dipole-sweep
depends_on:
  - THM-2384-owner-pivot-primal-dual-obstruction-and-two-probe-repair-corner
related:
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2379-anchored-guard-unit-deletion-factor-repair-current
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
  - THM-2464-linked-blocker-depth-threshold-and-clock-cell-universality
  - THM-2502-endpoint-boolean-newton-carry-tournament-and-dipole-boundary
script: 04-computation/lrc14_active_dipole_sweep_thm2503.py
output: 05-knowledge/results/lrc14_active_dipole_sweep_thm2503.out
script_sha256: d7fd292162c6867280b9ca87fc5d6ad09e2299a9444e098e48fb1cb5c745195f
output_sha256: ed66e4aa4c525ee1e7efacdefc8c7e020d342aa6f57b2c0a45f2cb55c65c8957
hash_basis: working-tree bytes (LF)
---

# THM-2503 -- the sharp active two-probe sweep

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2384 already isolates the tempting but type-incorrect two-probe
repair: shifting a failed role and a helper in opposite directions
produces positive analytic mass, but the shift vector need not lie in
the owner-pivot dual and therefore is not a lawful target action.

The present theorem does not revisit that settled obstruction. It
classifies the underlying one-dimensional two-role sweep exactly. The
new information is its sharp profile atlas, its many-shift service
floor, and the full rational cyclotomic spectrum. These are useful
sidecars if a future theorem supplies the missing lawful target map;
they do not supply that map themselves.

## 1. Active duplicate-probe sweep

Work on `T=R/Z` and put

```text
p=13,

d_L(y)=1_(||y||<L/14),       u_L=1-d_L,       L in {1,2}. (1)
```

Let `v,a` be positive integers, let `beta,gamma in T`, and put

```text
V(x)=vx-beta,                    A(x)=ax-gamma.       (2)
```

Their strict-boundary preimages are finite, hence null for every
integrable density below. Let `w>=0` be integrable with

```text
rho=integral_T w(x)dx>0,

support(w) subset {d_L(V(x))=1} intersection {u_1(A(x))=1}. (3)
```

Thus the `V` role fails at the anchored shift and the helper `A` role
is safe there. For `ell in F_13`, define

```text
J_ell
 =integral_T w(x)
    u_L(V(x)-ell/13)
    u_1(A(x)+ell/13) dx.                              (4)
```

The opposite signs make `ell` an active balanced two-probe shift.
Because the failed `V` anchor remains inside `w`, this is duplicate-
probe service in the sense of THM-2384, not replacement of the unique
physical endpoint factors.

## 2. Pointwise identity and sharp profile atlas

At a fixed nonboundary point (hence for almost every `x`) write

```text
D_V={ell:d_L(V-ell/13)=1},

D_A={ell:d_1(A+ell/13)=1}.                           (5)
```

The number of shifts active in (4) is

```text
N(V,A)=13-|D_V union D_A|.                           (6)
```

The exact thirteen-shift count law from THM-2379/2384 gives

```text
|D_V|=2L-d_L(13V),

|D_A|=2-d_1(13A).                                   (7)
```

Inclusion--exclusion therefore gives the exact refinement

```text
N(V,A)
 =11-2L+d_L(13V)+d_1(13A)+|D_V intersect D_A|.      (8)
```

This identity shows `N>=11-2L`. It also exposes the complete profile
atlas. On the support (3), the possible `D_V` are:

```text
L=1:
 {0}, {-1,0}, {0,1};

L=2:
 {-2,-1,0}, {-1,0,1}, {0,1,2},
 {-3,-2,-1,0}, {-2,-1,0,1},
 {-1,0,1,2}, {0,1,2,3}.                             (9)
```

The helper profiles are the `12` nonzero singletons and the `11`
adjacent pairs which avoid zero:

```text
D_A={r},                         r!=0,

D_A={r,r+1},                    0 notin {r,r+1}.     (10)
```

Every listed profile occurs on a nonempty rational open cell. Counting
the Cartesian pairs in (9)--(10) by (6) gives

```text
L=1:   N=9,10,11       with counts 20,35,14;         (11)

L=2:   N=7,8,9,10      with counts 32,69,52,8.       (12)
```

The counts sum to `3*23=69` and `7*23=161`. Hence the ranges are sharp.

### Proof of the upper endpoints

The support condition says `0 in D_V` and `0 notin D_A`. For `L=1`,
if `D_V={0}` then the nonempty `D_A` adds a second residue; otherwise
`D_V` already has two. Thus `|D_V union D_A|>=2` and `N<=11`.
For `L=2`, every `D_V` has at least three residues, so `N<=10`.
Together with (8) and the realized profiles, this proves (11)--(12).

## 3. Integrated service consequences

Put

```text
T=sum_(ell in F_13)J_ell.
```

Fubini and (6) give

```text
T=integral_T w(x)N(V(x),A(x))dx.                    (13)
```

Therefore

```text
(11-2L)rho<=T<=(12-L)rho.                           (14)
```

Both endpoints are sharp at the profile level. Moreover

```text
J_0=0,                                             (15)
```

because `d_L(V)=1` on the support of `w`. Each `J_ell<=rho`, so (14)
forces at least `11-2L` nonzero shifts. Averaging over the twelve
nonzero shifts gives the preliminary bounds `3rho/4` and `7rho/12`.
The profile geometry improves them sharply.

For `L=1`, take

```text
E_1={2,4,6,8,10}.
```

Every `D_V` in (9) avoids `E_1`, while every singleton or adjacent
pair in (10) meets it at most once. Thus at least four shifts in
`E_1` are active at every point. For `L=2`, the same argument with

```text
E_2={4,6,8}
```

leaves at least two active shifts. Integrating and averaging on these
fixed hitting sets gives

```text
max_(ell!=0)J_ell
 >=4rho/5                         if L=1,

 >=2rho/3                         if L=2.            (16)
```

Both constants are sharp. For `L=1`, assign weight `rho/5` to the five
profile pairs

```text
({0},{4,5}), ({0},{6,7}), ({0},{10,11}),
({0,1},{2,3}), ({-1,0},{8,9}).
```

Then every nonzero `J_ell=4rho/5`. For `L=2`, assign weight `rho/3`
to

```text
({0,1,2,3},{8,9}), ({-3,-2,-1,0},{6,7}),
({-2,-1,0,1},{4,5}).
```

The largest entry is exactly `2rho/3`. These floors are substantially
larger than THM-2384's mixed-coefficient floor, but they remain mass
statements in the coarser auxiliary shift coordinate.

## 4. Fourier consequences

Let `zeta=exp(2*pi*i/13)` and normalize

```text
Jhat(k)=(1/13)sum_ell J_ell zeta^(k ell).            (17)
```

Equations (13)--(15) imply, without rationality,

```text
sum_(k!=0)Jhat(k)=-T/13.                             (18)
```

Hence some `k!=0` satisfies

```text
Re Jhat(k)<=-T/156.                                 (19)
```

Parseval and Cauchy--Schwarz on the twelve nonzero shifts give

```text
sum_(k!=0)|Jhat(k)|^2
 =(1/13)sum_ell J_ell^2-T^2/169
 >=T^2/156-T^2/169
 =T^2/2028.                                        (20)
```

For `L=1`, the five-profile equality packet in Section 3 has all
twelve nonzero entries equal, so (20) is sharp.

The guard profiles force a stronger inequality. Put

```text
h=(0,9,9,11,12,12,12,12,12,12,11,9,9)/130.         (21)
```

For each of the `161` guard profile pairs, let `f in {0,1}^13` be its
active-shift word. Direct comparison of the explicit lists (9)--(10)
gives

```text
<h,f>-(11/130)|f| in {0,1,2,3,4,6}/130,             (22)
```

with multiplicities `25,40,38,16,32,10`, respectively. Hence
`<h,J> >= (11/130)T`. Since `||h||^2=11/130`, Cauchy--Schwarz and
Parseval improve (20) to

```text
sum_(k!=0)|Jhat(k)|^2>=T^2/1690,          L=2.       (23)
```

This is sharp in the guard profile cone. Indeed `h` itself is the
following nonnegative combination of admissible `(D_V,D_A)` profiles:

```text
weight 1/130: ({-3,-2,-1,0},{4,5});
weight 1/65:  ({0,1,2},{4,5});
weight 1/65:  ({-2,-1,0},{6,7});
weight 3/130: ({0,1,2,3},{4,5});
weight 2/65:  ({0,1,2,3},{6,7});
weight 3/65:  ({-3,-2,-1,0},{8,9}).
```

Its total `T=1`, and its nonzero energy is exactly `1/1690`.

Using the service floors (14), the corresponding uniform consequences are

```text
L=1:  27rho^2/676,

L=2:  49rho^2/1690.                                (24)
```

All equality mixtures above occur in the stated affine class. Take
`v=1` and choose the helper speed `a` so large that `x->ax-gamma`
makes more than one full turn on every selected rational `V`-profile
cell. It then meets every required `A`-profile cell on a rational open
subinterval. Choose disjoint such subintervals and rational step
heights for `w` to realize the displayed masses exactly.

Now assume the phase packets and weight are rational finite-interval
step data, so every `J_ell in Q`. If `Jhat(k)=0` for one `k!=0`, the
rational polynomial

```text
P(X)=sum_(ell=0)^12 J_ell X^ell
```

vanishes at a primitive thirteenth root. Thus `Phi_13` divides `P`.
Both have degree at most twelve, so `P=c Phi_13`; all thirteen
coefficients are equal. Equation (15) forces `c=0`, contradicting
`T>0`. Consequently

```text
Jhat(k)!=0                         for every k!=0.   (25)
```

Rationality is load-bearing only for (25). Equations (13)--(24) hold
for arbitrary integrable `w` under (3).

## 5. Exact LRC type boundary

The map in (4) has a source and target only in the auxiliary probe
space:

```text
source: failed V role plus safe helper A inside w;

target: the thirteen opposite test translations ell;

preserved: positive mass, both anchored statuses, and rationality;

lost: endpoint replacement, owner-pivot duality, target phase,
      deepest harmonic, semantic word, and relation address.      (26)
```

THM-2384 already proves the decisive obstruction. A balanced relation
vector is a primal class in `K/L`; a lawful target co-shift is a
covector in `L^perp/<w>`. The two-factor shift suggested by the
omitted-unit/helper relation need not annihilate the owner-pivot row
space. Its nonzero colour can therefore coexist with zero canonical
target.

The genuine THM-2350 target of an atomic relation is instead

```text
pi(r_full)
 =(
  delta_a(u)-delta_a(v)+m 1_(d=a),
   delta_b(u)-delta_b(v)+m 1_(d=b)
  ).                                                   (27)
```

The auxiliary sweep becomes one THM-2350 axis response only after all
of the following sidecars are supplied:

1. one fixed owner pivot identifies the ordered roles as
   `(A,V)=(a,k_a)` or `(b,k_b)`, so `e_A-e_V` lies in the true dual
   `L^perp/<w>`;
2. the co-shift moves **every** occurrence of both labels at present
   and bare endpoints, rather than leaving the anchored copies inside
   `w` frozen; and
3. one common endpoint/deepest gauge retains the `m e_d` contribution.

Without these, a shift colour can cancel in the true target. If the
frozen packet contributes dipole residue `c` and (17) selects colour
`k`, the actual axis colour with the displayed positive-exponent DFT
convention is

```text
c-k+m 1_(d=axis).
```

Choosing `c=k-m 1_(d=axis)` gives canonical target zero
while `Jhat(k)` remains nonzero. Even after lawful axis typing, (17)
aggregates one fibre of the two-dimensional `F_13^2` target and need
not select the same target as a semantic word. THM-2380 still needs
endpoint-matched charged quadrature, while THM-2502's four-cell
cancellation remains when the old-status pieces are reassembled.

Thus the theorem completes the sharp analytic classification of the
duplicate-probe sweep but does **not** solve THM-2384's stated honest
target-side mass problem. No scalar row is removed; the ledger remains
`165`, and LRC(14) remains open.

## 6. Exact companion

Run

```text
python3 04-computation/lrc14_active_dipole_sweep_thm2503.py
python3 -O 04-computation/lrc14_active_dipole_sweep_thm2503.py
```

The dependency-free `Fraction` companion:

- constructs every rational open profile cell and verifies (7)--(8);
- reproduces all `69+161` profile pairs and both sharp censuses;
- checks the two universal hitting-set proofs and their sharp equality
  packets;
- verifies the guard dual gap census, its six-profile primal
  decomposition, and both sharp energy constants;
- checks the service and support floors on a rational mixture of every
  profile;
- reconstructs all twelve cyclotomic evaluation matrices over `Q`,
  each of rank twelve with constant kernel; and
- verifies every nonzero colour on both positive-control mixtures.

Normal and optimized runs must reproduce

```text
05-knowledge/results/lrc14_active_dipole_sweep_thm2503.out
```

byte for byte. Every truth-bearing executable check uses `require`;
there are no floating-point truth tests.

An independent audit reconstructed the full profile cone, found the
sharp hitting sets and stronger guard dual, checked affine equality
realization, and repaired the boundary quantifier and DFT target sign.
It verified normal, optimized, and stored transcripts and confirmed
that no auxiliary colour is promoted to a canonical target.

QED.
