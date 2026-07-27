---
id: THM-2578
title: "Pure-tooth Abel-normal resonance and component-weight obstruction"
status: >
  PROVED + VERIFIED-EXACT.  For the complementary tooth pair
  P_s=d_L(kx+s/13), 1-P_s with k>0 and L in {1,2}, THM-2573's one-sided
  logarithmic Abel normal at physical offset M vanishes unless k divides M.
  At resonance M=kh it is
  k cos(pi hL/7) zeta^(-hs)/pi^2; the cosine never vanishes, and the
  plus-sign target DFT lands exactly at q=h mod 13.  In the live deep bank
  M=m c_3, gcd(m,91)=1, one has 13|M and gcd(k,13)=1.  Hence off resonance
  the aggregate tooth components cancel, while resonance forces q=0.  A
  pure danger/safe gate therefore supplies no nonzero live target Abel
  normal.  For a gate-supported total-layer handoff with nonuniform positive
  component weights, the off-resonant normal is exactly their tooth-index
  Fourier mode r=M mod k.  Its monodromy is the oriented component local
  system of THM-2574; a connection makes the mode 13-periodic but choosing
  an integer lift of r moves one unit of target colour between connection
  and carry.  A fixed target-neutral half-circle filter gives a lawful
  positive sharpness control with a nonconstant M=13 normal.  Thus the next
  positive input is a lawful nonuniform total-layer boundary filter, not the
  bare gate or an unphysical component phase.  No such filter is constructed
  on THM-2569, and no THM-2334 current, row exclusion, or LRC(14) conclusion
  follows.
source: common-endpoint-seam-2026-07-28-pure-tooth-normal
depends_on:
  - THM-2573-logarithmic-abel-normal-and-common-endpoint-jump-pairing
related:
  - THM-2568-full-x-transition-annihilation-and-refined-pair-drift-boundary
  - THM-2569-stationary-diagonal-conditioned-paired-corner-and-frozen-future-role-boundary
  - THM-2574-oriented-tooth-component-holonomy-and-fixed-frequency-descent
script: 04-computation/lrc14_pure_tooth_abel_resonance_thm2578.py
output: 05-knowledge/results/lrc14_pure_tooth_abel_resonance_thm2578.out
script_sha256: 098d08d03862385e9eca05b8a32e05a95eae26eace8a3cb2e1d7a17de0612901
output_sha256: 9a4437bf345ba2e1688637388717ab135dbe149dd3d044f0a6b694ece0ac4663
hash_basis: working-tree bytes (LF)
---

# THM-2578 -- the bare tooth normal is trapped at zero target colour

**PROVED + VERIFIED-EXACT.**

THM-2573 turns the first singular physical normal of a disjoint endpoint
pair into a positive common-jump measure.  The most immediate LRC candidate
is the target-active `k_a` danger comb against its safe complement.  Its
boundary is highly symmetric: every one of its physical teeth has the same
weight.  That symmetry is fatal at the live deep offsets.

The exact mechanism is

```text
uniform tooth boundary
  -> physical-frequency resonance k|M
  -> one target colour q=M/k;

13|M and k a 13-unit
  -> every resonant colour is q=0.                            (1)
```

Off resonance, the normal does not disappear pointwise.  Its `k` tooth
components cancel.  This identifies nonuniform component weight as the
precise positive datum that an oriented-tooth reference would have to see.

## 1. The pure complementary tooth pair

Put `p=13`, `zeta=exp(2 pi i/p)`, and for `L in {1,2}` write

```text
d_L(y)=1_(||y||<L/14),                  u_L=1-d_L.            (2)
```

Fix a positive integer `k`.  For `s in F_13`, define

```text
P_s(x)=d_L(kx+s/13),

L_s=P_s,                         R_s=1-P_s.                   (3)
```

The `2k` boundary points are distinct and may be labelled

```text
x_(j,epsilon,s)
 =(j+epsilon L/14-s/13)/k mod 1,

j in {0,...,k-1},                epsilon in {+1,-1}.          (4)
```

At each point in (4), `Delta R_s=-Delta L_s` and
`(Delta L_s)^2=1`.  Thus THM-2573's positive handoff measure is

```text
nu_s=sum_(j,epsilon)delta_(x_(j,epsilon,s)).                  (5)
```

For a physical offset `M in Z`, its Fourier moment is

```text
integral exp(2 pi i Mx)dnu_s(x)

 =exp(-2 pi i Ms/(13k))
   [sum_(epsilon=+-1)exp(2 pi i M epsilon L/(14k))]
   [sum_(j=0)^(k-1)exp(2 pi i Mj/k)].                         (6)
```

The last bracket is the complete cyclic character sum.  Hence it is zero
unless `k|M`.  If

```text
M=kh,                                                        (7)
```

then (6) is

```text
2k cos(pi hL/7) zeta^(-hs).                                  (8)
```

The cosine in (8) never vanishes.  Such a zero would require

```text
2hL=7 mod 14,                                                (9)
```

which equates an even integer with an odd one.

By THM-2573's one-sided normalization, the logarithmic Abel normal is

```text
N_s(M)=0,                                      k does not divide M;

N_s(kh)=k cos(pi hL/7)zeta^(-hs)/pi^2.          k divides M. (10)
```

Using the canonical plus-sign transform

```text
Nhat(q;M)=1/13 sum_s N_s(M)zeta^(qs),                        (11)
```

equations (10)--(11) give the exact support law

```text
Nhat(q;M)!=0
 iff k|M and q=M/k mod 13.                                   (12)
```

This is an aggregate whole-boundary theorem.  No physical component has
been oriented or discarded.

## 2. The live deep bank forces the resonant colour to zero

In the LRC endpoint triangle the physical offset is

```text
M=m c_3,                         gcd(m,91)=1.                 (13)
```

The deepest speed has positive `13`-adic valuation, so `13|M`.  Every
guard/unit role `k`, including the target-active `k_a`, is a `13`-unit:

```text
gcd(k,13)=1.                                                 (14)
```

There are only two cases.

If `k` does not divide `M`, (10) gives zero before the target transform.  If
`k|M`, Euclid's lemma applied to `M=13(M/13)` and (14) gives

```text
13 divides M/k.                                             (15)
```

The unique colour in (12) is therefore `q=0`.  Consequently

```text
Nhat(q;m c_3)=0

for every q!=0, gcd(m,91)=1, and k in the live role bank.    (16)
```

The condition `7 does not divide m` ensures that the separate deepest-comb
Fourier coefficient is nonzero; multiplying by that coefficient cannot
change (16).  This is a target no-go for the **bare complementary tooth**,
not for a fully filtered endpoint.

## 3. Nonuniform component weights are the exact aggregate escape

Suppose a completed nonnegative endpoint has a gate-supported handoff
measure on the same points (4), but with total-layer weights

```text
w_(j,epsilon,s)>=0.                                         (17)
```

These weights must be computed after every common filter and every
coincident endpoint is inserted.  Define the tooth-index transform

```text
W_(r,epsilon)(s)
 =sum_(j=0)^(k-1)w_(j,epsilon,s)exp(2 pi i rj/k),             (18)
```

where `r=M mod k` is initially represented by `0<=r<k`, and put

```text
h=(M-r)/k.                                                   (19)
```

The weighted version of (6) is

```text
exp(-2 pi i Ms/(13k))
 sum_epsilon exp(2 pi i M epsilon L/(14k))
               W_(r,epsilon)(s).                            (20)
```

Thus an off-resonant gate-boundary normal can survive only if the positive
weight profile has a nonzero component mode `r!=0`.  Uniform weights have

```text
W_(0,epsilon)=k,                 W_(r,epsilon)=0 for r!=0,   (21)
```

which is exactly the cancellation in Section 1.  Conversely the positive
profile with weight `4` on one component and `1` on all others has every
nontrivial component coefficient equal to `3` after the uniform cyclic sum
vanishes.  Nonuniform positivity can therefore break the component
cancellation; uniform mass or total perimeter cannot.

Equation (20) is only the gate-supported part.  A general pair of total
layers may have additional common jumps away from (4); THM-2573 applies to
their complete aggregate separately.

## 4. Component holonomy and what THM-2574 changes

The component labels in (4) live on a cover.  Advancing the target lift by a
full turn gives

```text
x_(j,epsilon,s+13)=x_(j-1,epsilon,s).                        (22)
```

A covariant physical weight family therefore obeys

```text
w_(j,epsilon,s+13)=w_(j-1,epsilon,s),                        (23)
```

and (18) has holonomy

```text
W_(r,epsilon)(s+13)
 =exp(2 pi i r/k)W_(r,epsilon)(s).                           (24)
```

The compensated component field

```text
G_(r,epsilon)(s)
 =exp(-2 pi i rs/(13k))W_(r,epsilon)(s)                      (25)
```

is `13`-periodic.  Substituting `M=r+kh` into (20) gives

```text
zeta^(-hs)
 sum_epsilon exp(2 pi i M epsilon L/(14k))
               G_(r,epsilon)(s).                            (26)
```

This is the boundary-weight version of THM-2574's oriented tooth-component
local system.  It states exactly how that theorem changes category:

```text
positive aggregate boundary measure
  -> sums all tooth components and cancels off resonance;

complex oriented component + compensating connection
  -> retains one component character before aggregation.    (27)
```

The connection is not canonical target data.  Replacing the integer lift
`r` by `r+k` leaves `W_r` unchanged, multiplies `G_r` by `zeta^(-s)`, and
replaces `h` by `h-1`; the physical product (26) is unchanged.  Thus one
unit of apparent target colour moves between the connection and the carry.
Choosing `0<=r<k` is an explicit oriented trivialization, not a consequence
of a positive Boolean endpoint.

THM-2574 can diagnose and retain a component mode at fixed frequency.  It
does not prove that the positive THM-2569 total-layer weights have such a
mode, make the connection a lawful THM-2334 factor, or resolve the frozen
target-informed selector.

## 5. A lawful fixed filter proves sharpness

The pure-gate no-go (16) fails as soon as a genuine common filter makes the
total boundary weights nonuniform.  This can happen without an artificial
component phase.

Take `k=L=1`, `M=13`, and the fixed target-neutral filter

```text
H=1_[0,1/2),

L_s=H P_s,                         R_s=H(1-P_s).              (28)
```

The gate boundaries are

```text
x_(epsilon,s)=epsilon/14-s/13 mod 1.                         (29)
```

None coincides with `0` or `1/2`.  Hence the total-layer handoff measure
retains exactly those points in (29) which lie in the support of `H`; the
endpoints of `H` themselves are not common jumps.  Over all thirteen target
shifts, the exact retained pattern is

```text
s=0:             epsilon=+1;
s=1,...,5:       none;
s=6,7:           epsilon=-1;
s=8,...,12:      both signs.                                (30)
```

There are thirteen retained atoms in total.  At `M=13`, the phase of one
retained atom depends only on its sign:

```text
exp(2 pi i 13 x_(epsilon,s))
 =exp(2 pi i 13 epsilon/14).                                (31)
```

The moment sequence in (30)--(31) is nonconstant, so invertibility of the
finite `C_13` transform forces at least one `q!=0` Abel-normal coefficient.
This is a lawful common-target family: `H` is fixed and target-neutral, and
only the actual gate is co-shifted.

The example is an abstract endpoint control, not a realized THM-2569 packet.
It proves the no-go boundary is sharp and names the missing positive object:

```text
a lawful common filter whose total-layer traces weight the tooth boundary
nonuniformly while retaining word, owner, future root, and deep unit mode. (32)
```

## 6. Scope for the live seam

The updated chain is

```text
complete danger/safe pair + ordinary full-X sum
  -> zero in every target character;                         THM-2568

same pair + logarithmic Abel normal
  -> positive boundary handoff measure;                      THM-2573

bare k_a tooth + live deep offset
  -> off-resonant component cancellation or q=0 resonance;  THM-2578

nonuniform boundary mode + oriented connection
  -> exact possible fixed-component target carrier.          (33)
```

The last line is a carrier specification, not a current theorem on the live
packet.  The highest-value next test is to compute the actual aggregate jump
weights after a covariant completion of the THM-2569 selector and stationary
future factor, then ask whether any allowed `(m,r)` component mode in (18)
survives.  Freezing that selector would again violate MISTAKE-266.

No scalar row is excluded; the ledger remains `165`, and LRC(14) remains
open.

## 7. Exact companion

Run

```bash
python3 04-computation/lrc14_pure_tooth_abel_resonance_thm2578.py
python3 -O 04-computation/lrc14_pure_tooth_abel_resonance_thm2578.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_pure_tooth_abel_resonance_thm2578.out.
```

The dependency-free referee checks:

- `1,440` resonant and `50,560` off-resonant pure-tooth packets;
- `18,720` plus-sign target-transform identities and every cosine parity
  obstruction;
- `53,760` live `(k,c_3,m)` packets, split into `7,452` zero-colour
  resonances and `46,308` cancelling off-resonances;
- all `1,770` nontrivial uniform and positive-spiked component modes;
- `71,980` component-holonomy and `23,010` connection-descent identities;
  and
- the exact thirteen-atom fixed-filter control, with four distinct target
  moment profiles and no coincident filter/gate boundary.

The root-of-unity sum, cosine nonvanishing, Euclidean divisibility, and
finite-DFT implications are symbolic proofs above, not finite
extrapolations. **QED.**
