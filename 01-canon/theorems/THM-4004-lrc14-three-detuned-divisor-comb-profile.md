---
id: THM-4004
title: "LRC(14) three-detuned divisor-comb profile and the t<U scale floor"
status: >
  PROVED RELATIVE TO CITED LRCUpTo13 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. In THM-3818/3878's exact rank-eleven `11+2` branch, a divisor of
  t that divides eight body coordinates leaves a ten-speed harmonic pack and
  three labelled detuned combs. Their exact common-branch count gives a new
  necessary prime profile: every prime ell>=5 dividing t divides at most
  seven body coordinates, 3 divides at most eight, and 2 divides at most
  nine; equality at 2 is confined to scale one and a reduced odd exception
  pair of sum >7. In the untouched t<U branch, every survivor also has
  `U>=3,208,300,859`, while the literal component-swapped width criterion can
  never fire. These are reductions, not an LRC(14) proof.
source: bounded_residual + root / LRC14 long session, 2026-08-24
audit: >
  INDEPENDENT MATHEMATICAL AND EXACT AUDIT PASS. The open danger-comb count,
  branch multiplicities, strict union gate, prime thresholds, scale
  restrictions, crossing floor and role-swap no-go were independently
  rederived. The audit caught and removed one discovery-only false adjective:
  `(3,7)` is a valid two-lift hostile but not a minimal one; `(3,5)` is
  already hostile. Normal, optimized and independent stored outputs match.
depends_on:
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
  - THM-3878-lrc14-eleven-plus-two-harmonic-absorption-seam-collapse
  - THM-668-detuned-harmonic-dispatch
  - LRCUpTo13
related:
  - THM-3910-lrc14-auxiliary-center-erosion-and-t-sheet-variance-response
  - THM-4002-lrc14-signed-endpoint-cross-phase-and-fixed-scale-two-family
  - THM-4003-lrc14-scale-two-component-erosion-boundary-strip
script: 04-computation/lrc14_tltu_divisor_comb_profile_thm4004.py
output: 05-knowledge/results/lrc14_tltu_divisor_comb_profile_thm4004.out
independent_audit_script: 04-computation/lrc14_tltu_divisor_comb_profile_independent_audit_thm4004.py
independent_audit_output: 05-knowledge/results/lrc14_tltu_divisor_comb_profile_independent_audit_thm4004.out
script_sha256: 3563112bb3a9f8022635faf38109eb773d425d7fccf7cd9d4ddd490ac58961e8
output_sha256: 1ae46b86df5965f14a63d7dd860c2a23c61a58188cf10560b13f62f04e8adb4e
independent_audit_script_sha256: f3b25ed7184113383c79488cf778d89fee157c06751919b8629ad88202287dd6
independent_audit_output_sha256: 6bf2461d4a784ef4ca85491d6e231161d1c8d814ea67a1580e36469d353c5366
hash_basis: raw LF bytes
---

# THM-4004 -- three-detuned divisor combs and the `t<U` scale floor

**PROVED RELATIVE TO CITED LRCUpTo13 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED.** The theorem narrows both relative-scale lanes but closes neither.
LRC(14) remains open.

## 1. Typed row and the common branch label

Work in THM-3818/3878's exact rank-eleven, two-component `11+2` branch:

```text
n=s(u_1,...,u_11) direct-sum t(p,q),
s in {1,2},  gcd(s,t)=gcd(p,q)=1,
u primitive, p<q, p+q<=356.                           (1)
```

All thirteen speeds are distinct. Fix a divisor `d|t`. Since `gcd(s,t)=1`,
one has `gcd(s,d)=1`. Suppose `d` divides `k` of the eleven body coordinates.
Those `k` actual body speeds and both pair speeds are divisible by `d`; after
division they form a `(k+2)`-speed pack `H`. If `y` is a cited LRC witness for
that pack, all of the lifts

```text
x_a=(y+a)/d,                  0<=a<d,                  (2)
```

preserve every pack clearance. The same branch label `a` is retained for all
exceptional body owners. That shared finite label, rather than an unpaired
marginal count, is the essential carrier.

For one exceptional actual speed `delta=s u_j`, put

```text
g=gcd(d,delta),               m=d/g.                  (3)
```

Its phases on `(2)` form an order-`m` equally spaced orbit, each point repeated
`g` times among the labelled branches. The strict danger set
`||delta x_a||<1/14` is an open arc of length `1/7`. Its exact worst-case
branch count is

```text
b_d(delta)=g*(m/7)             if 7|m,
           g*ceil(m/7)         otherwise.             (4)
```

When `7|m`, openness prevents both endpoints of a block of `m/7` spacings
from being counted; otherwise the usual ceiling is attained. Formula `(4)`
is uniform in the cited pack phase `y`.

## 2. One and two detuned owners

If `d` divides exactly ten body coordinates, the divided pack has twelve
speeds and one owner is detuned. THM-668 applies directly and closes the row.
Equivalently, every unresolved row satisfies all eleven deletion conditions

```text
gcd(t,u_i:i!=j)=1              for every j.            (5)
```

Prime by prime, every prime divisor of `t` must miss at least two body
coordinates. Primitivity already forbids zero misses; harmonic detuning
forbids exactly one.

If `d` divides exactly nine body coordinates, the divided pack has eleven
speeds, hence cited clearance at least `1/12`, and two owners are detuned.
For orbit order `m>=2`, the normalized count from `(4)` is `1/2` only at
`m=2` and is strictly below `1/2` for every `m>=3`. Thus the two bad sets
cannot cover all `d` branches unless both exceptional orders equal two.

In that equality case, `d/2` divides both exceptional coordinates as well as
the nine divisible coordinates. Body primitivity forces `d/2=1`, hence

```text
d=2.                                                     (6)
```

Now `2|t` and `gcd(s,t)=1` force `s=1`; the two exceptions are odd. Divide
them by their gcd to obtain a primitive odd pair `(a,b)`. THM-3878's exact
two-lift geometry supplies a safe branch whenever `a+b<=7`. Therefore a
survivor with nine even body coordinates must obey

```text
s=1,                    a+b>7.                         (7)
```

The familiar pair `(3,7)` is a valid hostile to a phase-uniform selector,
but it is not minimal: `(3,5)` at phase `3/16` already gives clearance
`1/16` on both selected lifts. This correction changes neither `(7)` nor the
validity of the `(3,7)` witness.

## 3. Three detuned owners and the prime profile

Suppose `d` divides exactly eight body coordinates. The eight divided body
owners and the two divided pair owners form a ten-speed pack. Cited
`LRCUpTo13` supplies a phase `y` with clearance at least `1/11`. Three actual
body speeds `delta_1,delta_2,delta_3` remain detuned. Since their bad subsets
all live on the same branch label in `(2)`, the union bound proves the
owner-typed sufficient certificate

```text
b_d(delta_1)+b_d(delta_2)+b_d(delta_3)<d
  ==> the full row is 1/14-lonely.                     (8)
```

The inequality is sufficient, not necessary; equality carries no
conclusion.

Let a prime `ell>=5` divide `t` and exactly eight body coordinates. Each
exception has full order `ell`, because `ell` divides neither its body
coordinate nor `s`. Hence each contributes `ceil(ell/7)` bad branches, and

```text
3 ceil(ell/7)<ell                for every prime ell>=5. (9)
```

Indeed `ceil(ell/7)<=(ell+6)/7`, and `3(ell+6)<7ell` for `ell>=5`. Combining
`(5)--(9)` gives the necessary prime-incidence profile

```text
ell>=5:  ell divides at most seven body coordinates;
ell=3:   3 divides at most eight body coordinates;
ell=2:   2 divides at most nine body coordinates,
         with equality only under (7).                (10)
```

The `ell=3` boundary is sharp for this phase-free strict count. Take

```text
s=1, t=3, (p,q)=(1,4),
u=(1,10,11,6,9,15,18,21,24,27,30).                    (11)
```

The divided pack is `{1,...,10}`. At its valid phase `y=1/11`, actual owners
`1,11,10` spoil the three lifts at distances `1/33,0,1/33`. This is not a
counterexample: the same full row has clearance `1/11` at `x=4/33`. It is a
hostile only to choosing an arbitrary cited pack phase and then using the
strict count.

At the first closed prime, the control

```text
s=1, t=5, (p,q)=(1,4),
u=(1,2,3,10,15,25,30,35,40,45,50)                     (12)
```

again has divided pack `{1,...,10}`. From `y=1/11`, its five full-row
clearances are

```text
(1/55,1/11,1/11,1/11,1/11),                          (13)
```

so four labelled branches win. The small controls `(11)--(13)` test the
selector mechanism; they are not claimed to survive THM-3818's inherited
high-height filter.

## 4. Every `t<U` survivor lies above a large crossing floor

Now impose the previously untouched relative-scale lane

```text
t<U=max_i u_i.                                         (14)
```

Put `Q=91^6=567869252041`. For every body owner `i`, the two actual speeds
`s u_i` and `tp` have the primitive crossing relation with coefficients

```text
(tp/g_i,-s u_i/g_i),          g_i=gcd(tp,su_i).        (15)
```

THM-3818's unresolved `W=V_dec` branch forbids every support-at-most-three
crossing row of height at most `Q`. Therefore

```text
max(tp,su_i)/g_i>Q             for every i.            (16)
```

In particular,

```text
tp>Q,  or  s min_i u_i>Q.                              (17)
```

Using the maximum owner in `(16)` and `t<U` gives

```text
Q<max(tp,sU)/g_i<=max(p,s)U.                           (18)
```

The atlas inequalities `p<q`, `p+q<=356` imply `p<=177`, while `s<=2`.
Thus every row in this exact `t<U` branch satisfies

```text
U>Q/177,
U>=3,208,300,859.                                     (19)
```

This is a scale floor, not a loneliness certificate.

## 5. Why literal component swapping stops

Let `lambda(u)` be the maximum length of a connected component of `G(u)`.
The maximum body runner alone has safe components of length `6/(7U)`, and
`G(u)` is contained in their union. Hence

```text
lambda(u)<=6/(7U).                                     (20)
```

The component-swapped version of THM-3818's one-component coherent-grid
criterion would require `t lambda(u)>=1`. But `(14)` and `(20)` give

```text
t lambda(u)<6/7<1.                                    (21)
```

So the literal scalar-width replay can never fire in `t<U`. This refutes only
that operation. A shifted `t`-grid may still hit a union of many short body
components; the phase-labelled union remains open.

## 6. Scope, incoming signal and replay

THM-4002 retains signed endpoint cross-phase and closes fixed scale-two
bodies in the opposite lane `t>=max E`. It supplies neither the divisor
profile `(10)` nor the `t<U` scale floor. Conversely, `(8)` forgets which
pack phase should be selected and does not reproduce THM-4002's fixed-body
closure. The two signals are orthogonal and compatible.

The remaining live object is the phase-labelled projection

```text
B_(u;s,t)={z:exists a in {0,...,t-1},
              s(z+a)/t in G(u)}.                       (22)
```

A swapped proof asks whether `G(p,q)` meets `(22)`, after imposing `(10)`.
Total safe mass, component count, or `lambda(u)` alone discards the branch
label and cannot answer this. The prime-two and prime-three equality cases,
multi-component arrival, the full `t<U` branch and LRC(14) remain open.

Reproduce from the repository root:

```text
python3 -B 04-computation/lrc14_tltu_divisor_comb_profile_thm4004.py
python3 -B -O 04-computation/lrc14_tltu_divisor_comb_profile_thm4004.py
python3 -B 04-computation/lrc14_tltu_divisor_comb_profile_independent_audit_thm4004.py
python3 -B -O 04-computation/lrc14_tltu_divisor_comb_profile_independent_audit_thm4004.py
```

The normal and optimized primary outputs and the independent audit match the
frozen raw-LF streams. **QED.**
