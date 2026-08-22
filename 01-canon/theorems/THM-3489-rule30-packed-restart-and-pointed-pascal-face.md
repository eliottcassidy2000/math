---
id: THM-3489
title: "Rule 30 packed restart, wrapped-center closure, and pointed Pascal face"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The finite-width packed
  orbit restarts after P_k steps with only its lift bit toggled.  Consequently
  every Rule 30 center at a depth k>=P_k is determined exactly by the period
  lift: c_k=epsilon_k for P_k<=k<2P_k and c_k=0 at k=2P_k.  At every proper
  terminal arc, the marked point is an explicit sparse high-order Pascal face
  which is disjoint from, and independently load-bearing beyond, the arc
  operator's low-order rank-loss sector.  Odd depths use only odd Hasse
  orders.  ANF degree and Walsh support still do not determine the calibrated
  point.  No Rule 30 prize is claimed.
source: root-rule30-next-targets-20260816
audit: >
  An independent hostile audit rederived the packed restart and wrapped-center
  closure, both pointed Pascal-face formulas, the exact transversality and
  coordinate-minimality statements, and the self-addressed ANF/Walsh boundary.
  It matched every declared verifier universe and confirmed that ordinary and
  optimized runs equal the stored transcript byte-for-byte: ACCEPT.
depends_on:
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
  - THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum
  - THM-3481-rule30-cyclic-arc-norm-rank-and-marked-innovation-spectrum
related:
  - THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary
  - THM-3468-rule30-radial-green-fold-innovation-discrepancy-and-fixed-seed-carrier-boundaries
  - THM-3476-rule30-depth-four-transverse-jet-barrier-and-slack-pascal-atlas
script: 04-computation/rule30_packed_restart_pointed_pascal_thm3489.py
output: 05-knowledge/results/rule30_packed_restart_pointed_pascal_thm3489.out
script_sha256: bd767efc9738b8c9c5b9afe95333fe7ee46948486ffd815d9aaf052ae76813f6
output_sha256: 9be3ad840be4ea88b5b4beacb08c05a571b2ab14b14c53caa366eb81a62499b6
hash_basis: raw bytes
---

# THM-3489 -- Rule 30 packed restart, wrapped-center closure, and pointed Pascal face

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-3481 computes exactly how a cyclic arc operator loses phase-current
information.  The present theorem asks the complementary pointed question:
which part of the surviving information is read by the physical endpoint?
There are two sharply different regimes.  A depth which has reached its base
period restarts exactly and its center closes from the lift bit.  Before that
wrap, the center lives on a sparse face of *high* Hasse orders, entirely above
the repeated-root sector killed by arc integration.

## 1. Inheritance, conventions, and live objects

Use THM-3458's inward right-edge packing

```text
R_0=1,
R_(t+1)=Phi(R_t),
Phi(x)=x xor ((2x) or (4x)),                           (1)
```

and let `Phi_w` be reduction modulo `2^w`.  If `x_i` denotes bit `i`, then

```text
(Phi x)_i=x_i xor (x_(i-1) or x_(i-2)),               (2)
```

with exterior negative bits zero.  For `k>=1`, let `P_k` be the seed-return
period modulo `2^k` and

```text
epsilon_k=bit_k(R_(P_k)),
P_(k+1)=2^(epsilon_k)P_k.                              (3)
```

Every `P_k` is a power of two, and THM-3458 proves

```text
ceil(k/2)<=P_k.                                       (4)
```

Retain THM-3471's phase readout and current, for `k>=2`,

```text
T_k(h)=b_k(h+k),
T_k(h+1)+T_k(h)=Q_k(h),
xor_(h mod P_k)Q_k(h)=epsilon_k,                      (5)
```

where `b_k(t)=bit_k(R_t)`.  The light-cone boundary gives

```text
T_k(-k)=0,
T_k(0)=c_k,
c_k=xor_(h=-k)^(-1)Q_k(h).                            (6)
```

For a power-of-two cycle `p`, write

```text
A_s q(h)=xor_(0<=i<s)q(h+i),
M_j(f)=xor_(0<=h<p) binom(h,j)f(h).                   (7)
```

The `M_j` are the phase Hasse moments at `X=1`.  They are distinct from
THM-3476's Hasse jets in the **transport-slack** marker.  Their shared Pascal
algebra is an analogy, not an identification of coordinates.

The inheritance pass is:

1. closest proved mechanism: THM-3481's repeated-root arc rank;
2. canonical hostile: phase rotations preserve every unpointed profile while
   moving the physical sample through both bit values;
3. corrected near miss: retain the light-cone-calibrated phase, not only the
   arc spectrum; and
4. least-used sidecar: the ordered terminal interval in the native phase
   coordinate.

The live board is the packed restart cocycle, the lift bit, the terminal
current, its pointed Pascal face, the innovation cube, and the physical
basepoint.

## 2. The packed restart identity

Put

```text
p=P_k,
e=epsilon_k.                                         (8)
```

### Theorem 2.1 (one-cycle restart)

For every `k>=1` and every `t>=0`,

```text
boxed:
R_(p+t)=R_t xor e*2^k          mod 2^(k+1).           (9)
```

More generally, if `q>=0` and `0<=s<p`, then

```text
R_(qp+s)=R_s xor (q mod 2)e*2^k mod 2^(k+1).         (10)
```

### Proof

The seed return modulo `2^k` and the definition of `e` say exactly

```text
R_p=1 xor e*2^k mod 2^(k+1).                          (11)
```

Equation (2) shows that toggling the top bit commutes with the truncated map:

```text
Phi_(k+1)(x xor 2^k)=Phi_(k+1)(x) xor 2^k.           (12)
```

Indeed, every lower output bit is unchanged, and the top output bit contains
the input top bit additively.  Iterate (12), start from (11), and use the
semigroup law:

```text
R_(p+t)
 =Phi_(k+1)^t(R_p)
 =Phi_(k+1)^t(1 xor e*2^k)
 =R_t xor e*2^k.                                     (13)
```

Repeated application proves (10).

This is stronger than a return-time statement: after one base cycle, the
entire lower `(k+1)`-bit orbit restarts, with precisely one possible sheet
toggle.

### Corollary 2.2 (all wrapped centers close)

If `k>=P_k`, then

```text
boxed:
c_k=(floor(k/P_k) mod 2) epsilon_k.                  (14)
```

To prove this, write `k=qp+s`, where `q>=1` and `0<=s<p`.  The highest set
bit of the untruncated row `R_s` is exactly `2s`.  It lies strictly below
position `k`: for `q=1`, `2s<p+s=k`; for `q>=2`, `2s<2p<=qp<=k`.
Taking bit `k` in (10) therefore leaves only the lift toggle.

The lower bound (4) now makes the complete physical boundary table

```text
k<p:          c_k remains the raw marked terminal arc;
k=p:          c_k=epsilon_k;
p<k<2p:       c_k=epsilon_k;
k=2p:         c_k=0.                                  (15)
```

At `k=2p`, the endpoint is necessarily innovative.  The highest set bit of
`R_p` is `2p=k`, so

```text
epsilon_k=bit_k(R_p)=1.                               (16)
```

The two toggles in (10) nevertheless cancel at time `2p`, giving `c_k=0`.
For `k>=2`, both endpoint equalities `k=p` and `k=2p` are even-depth cases,
because `p` is a power of two.  Thus any **odd innovative** strict-wrap depth
has

```text
p<k<2p, epsilon_k=1  ==>  c_k=1.                     (17)
```

No assertion is made that a strict wrap or either endpoint occurs infinitely
often.

### Corollary 2.3 (the residual open line vanishes on a strict wrap)

If `p<k<2p`, put `s=k-p`.  THM-3471's wrap reduction is

```text
c_k=epsilon_k xor xor_(h=-s)^(-1)Q_k(h).             (18)
```

Together with (15), this proves the new pointed boundary circuit

```text
boxed:
xor_(h=-s)^(-1)Q_k(h)=0.                              (19)
```

The residual profile need not vanish at other phases.  Equation (19) is a
calibrated endpoint law, not a phase-average cancellation.

## 3. Point evaluation is a Pascal face

The next result is abstract and applies to every function on a power-of-two
cycle.

### Lemma 3.1 (Pascal inversion)

For `p=2^d` and `0<=r<p`,

```text
f(r)=xor_(r subseteq_bits j<p) M_j(f).                (20)
```

Indeed, the Taylor map in (7) is the transpose of the finite Pascal matrix.
That matrix is self-inverse over `F_2`.  Equivalently,

```text
f(r)=xor_(0<=j<p)binom(j,r)M_j(f),                    (21)
```

and Lucas' theorem gives the bit-containment condition in (20).

### Theorem 3.2 (the profile face)

Let `0<s<p`, put `r=p-s`, and let `f` be any phase profile.  Then

```text
boxed:
f(p-s)=xor_(a subseteq_bits s-1) M_(p-s+a)(f).        (22)
```

There are exactly

```text
2^popcount(s-1)                                      (23)
```

moments in this face.

Within `d` bits,

```text
p-s=(p-1)-(s-1)                                      (24)
```

is the bitwise complement of `s-1`.  Therefore the supersets of `p-s` are
exactly

```text
(p-s)+a,       a subseteq_bits s-1,                  (25)
```

with no carries.  Lemma 3.1 proves (22), and (23) is the number of submasks.

### Theorem 3.3 (the current face)

For `q in V_p`, define its terminal interval

```text
L_s(q)=xor_(h=p-s)^(p-1)q(h).                         (26)
```

Then

```text
boxed:
L_s(q)=
 xor_(a proper-submask-of s)
 M_(1+((p-s-1) bit_or a))(q).                        (27)
```

Here a proper submask satisfies `a bit_and complement(s)=0` and `a!=s`.
The current face contains exactly

```text
2^popcount(s)-1                                      (28)
```

moments.

To prove the formula, insert Pascal inversion into (26).  The coefficient of
`M_j(q)` is zero for `j<p-s`; otherwise it is

```text
xor_(h=p-s)^j binom(j,h)=binom(j-1,p-s-1) mod 2.     (29)
```

The last equality follows immediately by summing Pascal's recurrence: every
interior term occurs twice and only the lower endpoint survives.  Now

```text
p-s-1=(p-1)-s                                       (30)
```

is the `d`-bit complement of `s`.  Lucas' theorem therefore writes

```text
j-1=(p-s-1) bit_or a,       a subseteq_bits s.       (31)
```

The choice `a=s` would give `j=p`, outside the phase table, and every proper
submask gives one valid order.  This proves (27)--(28).

### 3.4 Physical specialization and wrap boundaries

For Rule 30 put

```text
p=P_k,
w=floor(k/p),
s=k mod p.                                           (32)
```

THM-3471 gives

```text
c_k=(w mod 2)epsilon_k xor L_s(Q_k).                 (33)
```

Let the complete terminal profile be

```text
F_k=A_k Q_k.                                         (34)
```

When `s>0`, its physical phase has standard representative `p-s`, so (22)
gives

```text
boxed:
c_k=xor_(a subseteq_bits s-1)
        M_(p-s+a)(F_k).                              (35)
```

Equation (27) gives the equivalent current formula

```text
c_k=(w mod 2)epsilon_k xor
 xor_(a proper-submask-of s)
 M_(1+((p-s-1) bit_or a))(Q_k).                      (36)
```

The endpoint `s=0` must be kept separate.  At `k=p`, `F_k` is the constant
`epsilon_k` profile; its only nonzero Hasse coordinate is `M_(p-1)`.  At
`k=2p`, the profile is zero.

For a strict wrap, (19) turns (27) into the high-face cancellation

```text
xor_(a proper-submask-of s)
 M_(1+((p-s-1) bit_or a))(Q_k)=0.                    (37)
```

Equivalently, apply (22) to the residual profile `A_sQ_k`: its pointed face
has XOR zero.  This relation is imposed by the Rule 30 seed calibration; it
is not a universal consequence of arc rank.

Finally, every odd physical depth `k>=3` has even `p` and odd nonzero residue
`s`.  Then `p-s` is odd while every submask of `s-1` is even.  Therefore

```text
boxed:
every Hasse order in the physical profile face is odd. (38)
```

This is a parity localization of the point functional, not a prediction of
its value.

## 4. The pointed face is transverse to rank loss

Write

```text
s=2^a u,       u odd,
ell=2^a-1.                                           (39)
```

THM-3481 proves

```text
im(A_s)=(S+1)^ell V_p,                               (40)
```

or, in Hasse coordinates,

```text
M_0(f)=...=M_(ell-1)(f)=0,                           (41)
```

with every coordinate `M_j`, `j>=ell`, otherwise free.  Since `s<p` and
both `s` and `p` are divisible by `2^a`,

```text
p-s is a positive multiple of 2^a,
p-s>=2^a=ell+1.                                     (42)
```

Every pointed-face order in (22) is therefore strictly above the entire
forced-zero sector.  The same conclusion applies to `A_(p+s)`, because
`nu_2(p+s)=nu_2(s)`.

### Theorem 4.1 (independent load-bearing coordinates)

Every one of the `2^popcount(s-1)` coordinates in (22) is independently
load-bearing on `im(A_s)`.  More precisely, for each face order `j_0` there
are two profiles in `im(A_s)` which agree in every Hasse coordinate except
`M_(j_0)` and have opposite values at `p-s`.

### Proof

Identify a profile with its polynomial in

```text
F_2[X]/((X+1)^p),       Y=X+1.                       (43)
```

By (42), `Y^(j_0)` lies in the image ideal `(Y^ell)`.  Its Hasse vector has a
single one, at order `j_0`.  Since `j_0` contains `p-s` in bits, its
coefficient at `X^(p-s)` is

```text
binom(j_0,p-s)=1 mod 2.                              (44)
```

Adding this image profile toggles exactly the named Hasse coordinate and the
marked point.  This proves the claim.

Consequently no proper subset of the face coordinates universally determines
the endpoint on the arc image.  Arc integration can lose a large
low-repeated-root sector while leaving the marked functional fully exposed
on independent high orders.

There is an equally sharp current-side separation.  In the hard unwrapped
regime `k<p`,

```text
epsilon_k=M_0(Q_k),                                  (45)
```

whereas every index in (36) is at least `p-k>=1`.  On the ambient current
module, holonomy and the marked terminal functional are independent linear
functionals.  This does not say that the particular Rule 30 currents are
arbitrary; it says that no operator identity can infer the hard center from
the lift bit alone.

The restart theorem shows exactly where extra dynamics *does* add such an
identity: on a strict wrap, Rule 30 forces the otherwise independent
high-face relation (37).

## 5. ANF, Walsh, and the self-addressed physical recursion

Suppose `epsilon_k=1`.  Let

```text
Gamma_k:Z/P_k Z -> F_2^d,
d=log_2 P_k,                                          (46)
```

be the innovation-coordinate bijection below depth `k`, and let `tau_k`
represent phase addition by one.  Pull the terminal profile to the cube:

```text
mathcalF_k(Gamma_k(h))=F_k(h).                        (47)
```

The physical origin has an exact recursive address.  If
`kappa_1<...<kappa_(d)` are the earlier innovation depths, then

```text
Gamma_k(0)
 =(T_(kappa_i)(0))_i
 =(c_(kappa_i))_i.                                   (48)
```

The marked terminal phase is

```text
x_k^*=Gamma_k(-k)=tau_k^(-k)Gamma_k(0),              (49)
```

and therefore

```text
boxed:
c_k=mathcalF_k(
 tau_k^(-k)(c_(kappa_1),...,c_(kappa_d))).           (50)
```

At successive innovation depths this is a deterministic triangular,
self-addressed recursion on the preceding innovation-center prefix.  It is
not a probabilistic martingale: both `mathcalF_k` and `tau_k` change with
depth, and the selected prefix is the calibrated seed address rather than a
Haar sample.

### 5.1 What maximal ANF degree does at one point

Write the Boolean ANF of `mathcalF_k` as

```text
mathcalF_k(x)=xor_(B subseteq [d]) alpha_B product_(i in B)x_i. (51)
```

At the physical address,

```text
c_k=xor_(B subseteq support(x_k^*))alpha_B.          (52)
```

For an odd innovative `k`, THM-3481 gives `alpha_[d]=1`.  But this top
coefficient contributes to (52) only when `x_k^*` is the all-one address;
even then it can cancel with lower coefficients.  Maximal degree is therefore
an exact global obstruction and not a point-value rule.

### 5.2 Walsh inversion and the final signed carry

Put

```text
g_k(x)=(-1)^mathcalF_k(x),
What_k(B)=sum_x g_k(x)(-1)^(B dot x).                 (53)
```

Fourier inversion at the physical address is

```text
boxed:
2^d(-1)^c_k=
 sum_B What_k(B)(-1)^(B dot x_k^*).                  (54)
```

At an odd innovative depth with `d>=2`, THM-3481 proves

```text
What_k(B)=2u_B,       u_B odd,                        (55)
```

for every mode.  Dividing (54) by two gives

```text
sum_B u_B(-1)^(B dot x_k^*)
 =2^(d-1)(-1)^c_k.                                   (56)
```

Modulo `2^d`, the two possible right sides are both `2^(d-1)`.  Thus full
Walsh support and all mod-four residues are blind to the mark; the exact
signed alignment in (56) is the final binary carry which selects `c_k`.

This boundary is sharp.  Fix any cube point `x_*`, choose `y!=x_*`, and let

```text
f_1=indicator_(x_*),
f_0=indicator_y.                                     (57)
```

For `d>=2`, both tables have odd support and ANF degree `d`.  Their Walsh
magnitudes agree exactly:

```text
|What(empty)|=2^d-2,
|What(B)|=2 for B nonempty,                          (58)
```

and every coefficient is nonzero.  Yet `f_1(x_*)=1` and `f_0(x_*)=0`.
Because every odd arc operator is invertible, both tables occur as terminal
profiles of some phase currents.  This is an ambient operator hostile, not a
claim that either synthetic current is a Rule 30 current.

The restart theorem cleanly separates the remaining Prize-2 target:

```text
k>=P_k: c_k=(floor(k/P_k) mod 2)epsilon_k;
k<P_k:  high Pascal face / signed Walsh alignment remains. (59)
```

## 6. Preservation ledger, boundaries, and next target

| map or identity | preserves / proves | destroys or leaves open | required sidecar |
|---|---|---|---|
| one `P_k` restart | every lower bit and the one lift toggle | frequency of future wraps | width `k` and physical time |
| wrapped depth -> lift bit | exact center for `P_k<=k<=2P_k` | every unwrapped center | comparison of `k` with `P_k` |
| profile -> pointed Pascal face | exact marked coefficient | all moments outside the face | standard phase representative `p-s` |
| arc image -> low Hasse constraints | exact repeated-root loss | the independent high observer face | full high Hasse coordinates |
| innovation cube -> ANF/Walsh spectra | maximal degree and full spectral support | signed value at `x_k^*` | calibrated address and coefficient signs |
| Haar phase law -> physical recursion | spatial averages | temporal law of the deterministic prefix | cross-depth current/address coupling |

The next exact target is not another unpointed rank computation.  It is one
of the following equivalent hard-regime tasks:

1. control the sparse high-face parity in (35) for `k<P_k`;
2. control the terminal-current face in (36) beyond its independent ambient
   degrees of freedom; or
3. control the signed Walsh alignment (56) along the self-addressed prefixes
   (48)--(50).

The number of profile-face coordinates can be small -- for example it is two
when `s-1` is a power of two -- but their orders lie near `P_k`.  Sparsity of
the coordinate list is not a bound on the cost or discrepancy of its values.

This theorem proves no eventual absence or abundance of wraps, no center
nonperiodicity, no center density, and no random-access complexity bound.
The companion universe through depth `30` has physical wraps only at depths
`2,3,4`; this is labelled finite evidence and is not extrapolated.

## 7. Exact verification

Reproduce the companion with

```bash
python3 04-computation/rule30_packed_restart_pointed_pascal_thm3489.py
python3 -O 04-computation/rule30_packed_restart_pointed_pascal_thm3489.py
```

The abstract audit checks every pair

```text
p=2^d, 1<=d<=8,
1<=s<p,                                               (60)
```

for `502` total pairs.  It represents both sides of (22) and (27) as exact
bit-packed linear functionals, so equality covers **every** profile and
current, not a sample.  For `p<=64`, it additionally constructs all `120`
direct arc matrices, verifies their ranks, and finds `966` explicit
one-Hasse-coordinate image profiles witnessing independent load-bearing.

The Rule 30 audit checks every restart phase `0<=t<=P_k` for

```text
2<=k<=30,                                             (61)
```

giving `13,445` restart comparisons.  An independent centered local-rule
evolution agrees with the packed rows.  The audit also checks every physical
terminal current, both Pascal formulas, all available wrap cases, the
innovation-address recursion, ANF point faces, Walsh inversion, and the
low-two-adic-layer blindness.  Actual wrap controls in this universe are

```text
(k,p,epsilon,c)=(2,2,0,0),(3,2,1,1),(4,4,1,1).       (62)
```

There is no observed `k=2p` control in this finite universe; that endpoint is
proof-controlled by (16), not advertised as computationally instantiated.

Finally, the Boolean audit exhausts every odd-support table in dimensions
`2` and `3`, and every translated-singleton hostile pair in dimensions
`2,...,6`.  All gates use explicit `check` calls, so optimized mode removes
none of them.  Universal claims are the proofs above; the computation audits
implementations and sharp finite controls only.
