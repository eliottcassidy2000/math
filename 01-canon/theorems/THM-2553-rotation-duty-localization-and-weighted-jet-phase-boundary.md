---
id: THM-2553
title: "Rotation-duty localization stops at augmentation: the weighted jet-phase boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Prime rotation gives an exact scalar duty formula, but its nonzero residue
  lies in the augmentation component.  Nontrivial characters have zero
  augmentation, and THM-2356's required Gram weights replace the scalar count
  by an uncontrolled line Fourier coefficient.  A graph-preserving 13-orbit
  has identical complete chirp intensities while its target-zero member has
  zero detector and its twelve other members have positive detector.  This
  does not refute HYP-9050/9055's unweighted counts; it scopes their jet
  handoff.  No masked landing, row exclusion, or LRC(14) conclusion follows.
source: codex-2026-07-27-holotopy
depends_on:
  - THM-2356-finite-field-chirp-gram-tomography-and-bockstein-pairing
  - THM-2363-planar-graph-detector-dominates-word-support-energy
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
related:
  - HYP-9050
  - HYP-9055
  - MISTAKE-282
script: 04-computation/lrc14_rotation_localization_weighted_jet_boundary_thm2553.py
output: 05-knowledge/results/lrc14_rotation_localization_weighted_jet_boundary_thm2553.out
script_sha256: 0d9de6e7a39e3cf983f69f949a83744cf48014dcca148de2831c305f91089c44
output_sha256: fb81e4c2d572205e2e29913be718f44f43a6c8456149299364c1e0fbad640168
hash_basis: LF-normalized bytes
---

# THM-2553 -- rotation localizes scalar duty, not charged jet phase

HYP-9050 and HYP-9055 found exact `13`- and `91`-rotation congruences.  This
theorem determines precisely what those congruences can see after the
THM-2356/2363 target--jet weights are restored.

The answer has a sharp augmentation/charge split:

```text
unweighted scalar profile -> exact duty localization;
nontrivial character      -> zero augmentation, hence no mod-p forcing;
actual Gram weight         -> an arbitrary line Fourier coefficient.       (1)
```

## 1. Exact prime-rotation formula

Let `p` be prime, let `f:F_p->R` take values in any commutative ring, and let
`v_1,...,v_n in F_p`.  Put

```text
A=sum_(r in F_p)f(r),
z=#{i:v_i=0},
N_f(v)=sum_i sum_(j in F_p)f(v_i j).                       (2)
```

Every nonzero `v_i` permutes `F_p`, while a zero `v_i` contributes `p f(0)`.
Therefore

```text
N_f(v)=(n-z)A+p z f(0).                                   (3)
```

When `n=p`,

```text
N_f(v)=-zA mod p.                                         (4)
```

HYP-9050's band is `f=1_{0,+-1}`, so `A=3`; at `p=13`, (4) is exactly

```text
N_f(v)=10 z mod 13.                                       (5)
```

This is an orbit-permutation identity, not a general Burnside statement
about every rotation-covariant functional.

For the zero-intercept congruence (5), the sharp hypotheses are equal
coordinate weights, exactly `p` coordinates, and augmentation `A` invertible
modulo `p` if one wants to recover `z`.  For known arbitrary `n`, (3) still
gives an affine localization of `z` whenever `A` is invertible.  With fixed
coordinate weights `lambda_i`, the exact variant is

```text
N_(f,lambda)
 =A sum_(v_i!=0)lambda_i
  +p f(0)sum_(v_i=0)lambda_i.                              (6)
```

It localizes to `z` alone only under an additional total-weight law.  The
atomic gauge in THM-2337 (38)--(39), with exact ratio `-1/9456`, shows that
canonical atomic weights are not generally gauge-constant; it is not claimed
to be a typed unequal-weight witness along this particular rotation orbit.

## 2. Charged profiles are rotation-invisible modulo p

For a nontrivial additive character `f(r)=zeta^r`,

```text
A=sum_r zeta^r=0.                                         (7)
```

Thus (3), with equal weights and `n=p`, gives

```text
N_f(v)=p z                                                (8)
```

exactly in the cyclotomic ring.  It is divisible by `p`.  Any nonzero
residue in (4) comes from the constant/augmentation component, which carries
no target or jet phase.  A nonzero *census of degenerate labels* must not be
confused with a nonzero charged detector value.

## 3. The F169 rotation is an annihilator projector

Write

```text
K=F_169,          phi(q)=q^2/2,
D_h phi(y)=phi(y+h)-phi(y)-phi(h)=hy,                       (9)
```

and let `psi(x)=zeta^Tr(x)` be a nontrivial additive character of `K`.
For `e!=0`, rotation along the line `F_13 e` gives

```text
D_h phi(y+t e)=D_h phi(y)+t h e.                          (10)
```

Hence, for `a,h in K`,

```text
sum_(t in F_13) psi(a D_h phi(y+t e))
 =13 psi(aD_h phi(y))       if Tr(ahe)=0,
 =0                         otherwise.                     (11)
```

For fixed nonzero `h,e`, the annihilator

```text
{a in K:Tr(ahe)=0}                                        (12)
```

is one `F_13`-line: thirteen labels, twelve nonzero.  Its label count is
`-1 mod 13`, but every surviving character sum in (11) has the literal
factor `13`.  One rotation is merely the subgroup-averaging projector onto
that annihilator.  Since `K` is `C_13^2`, not cyclic `C_169`, a second
independent direction is needed even to resolve all characters.

The corresponding unweighted derivative incidence is universal.  For every
`e!=0` and `a in K`,

```text
I_e(a)=#{(t,y) in F_13 x K:D_(t e)phi(y)=a}
      =181  if a=0,
      =12   if a!=0.                                      (13)
```

Indeed `t=0` contributes all `169` points to zero; each of the twelve
nonzero `t` makes `y->tey` a bijection.  Thus every value is `12 mod 13`,
independent of the signal being detected.

## 4. The THM-2356 weight is the missing line Fourier coefficient

The graph detector does not sum (11) with constant weight.  Its pair product
has the form

```text
g_t=Z(y+t e+h) conjugate(Z(y+t e)).                         (14)
```

The rotated contribution is therefore, up to the displayed harmless scalar
phase,

```text
psi(bh+a(D_h phi(y)+phi(h)))
  sum_t g_t zeta^(t Tr(ahe)),                              (15)
```

an arbitrary length-thirteen Fourier coefficient of the unknown word
`(g_t)`.  HYP-9050 controls only `g_t=1`.  At `h=0`, (15) contains exactly
the labelled singleton line mass

```text
sum_t |Z(y+t e)|^2.                                       (16)
```

This is the first failed implication in the proposed rotation-to-jet bridge.
The missing sidecar is a controlled nonzero pair-product line coefficient,
or an absolute singleton location, not another unweighted duty count.

## 5. Sharp graph-preserving rotation hostile

Fix a graph label `c in K`.  For the three THM-2337 word masks choose

```text
e_{ {a} }=(1,0),       e_{ {b} }=(0,1),
e_{ {a,b} }=(1,1)                                      (17)
```

in the canonical `F_13^2` target coordinates.  For `t in F_13`, let the
joint array `A_t^sigma` have its sole nonzero coefficient at

```text
(q,z)=(t e_sigma, (t e_sigma)^2/2+c).                       (18)
```

These thirteen arrays are one orbit of the graph-preserving Heisenberg
translation

```text
T_e(q,z)=(q+e,z+eq+e^2/2).                                (19)
```

Their graph signals are `Z_t=delta_(t e_sigma)`.  For every one of the
`169^2=28561` chirp masks, the amplitude is one root of unity, so every
intensity is exactly `1`, independently of `t`.  Nevertheless THM-2363's
coefficient detector and word energy satisfy

```text
t=0:     D_graph=0,       E_sigma=0;
t!=0:    D_graph=2,       E_sigma=1.                       (20)
```

The unique target-zero member and all twelve live members are therefore
indistinguishable by the complete intensity/incidence table.  Across the
three masks this is an exact hostile bank of `39` signals and
`1,113,879` chirp entries.  Orbit summation produces thirteen identical
intensity contributions, hence zero modulo thirteen.  Marking the exceptional
member separately is precisely the absent singleton-location sidecar.

This hostile is compatible with, and conceptually refines, THM-2356's
existing exhaustive singleton-location boundary.  It does not refute the
unweighted sector formulas in HYP-9050/9055.

## 6. Strongest positive survivor

There is one clean conditional statement.  If the support of a nonzero joint
array is invariant under a nontrivial graph-preserving `T_e`, then any
supported coefficient generates a thirteen-point orbit on one graph.  At
most one point of that orbit has target `q=0`; at least twelve have nonzero
target.  Therefore

```text
D_graph>0.                                                (21)
```

Current canon supplies no such support covariance.  Physical full-factor
translation preserves different data, while the factor-only jet gauge of
THM-2337 is not weight-preserving.  Even (21) gives refined graph survival,
not a coarse THM-2334 landing or a fixed word mask unless those extra labels
are separately retained.

No row is excluded.  LRC(14) remains open.

## 7. Exact companion

The dependency-free companion verifies all `8192*14` Boolean scalar profiles,
the unequal-weight hostile, the full `28,392` derivative-incidence table,
all `28,224` annihilator lines, `28,392` graph-preserving orbits, and all
`1,113,879` singleton chirp intensities together with literal graph-label,
word-mask, `D_graph`, and `E_sigma` calculations in (20).  Run

```bash
python3 04-computation/lrc14_rotation_localization_weighted_jet_boundary_thm2553.py
python3 -O 04-computation/lrc14_rotation_localization_weighted_jet_boundary_thm2553.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_rotation_localization_weighted_jet_boundary_thm2553.out
```

after LF normalization.  Every check raises explicitly under optimized
Python.

## 8. Independent hostile audit

An independent audit rederived the scalar and weighted prime-rotation
formulas, the augmentation/charge split, the `F_169` trace derivative and
`181/12` incidence law, every annihilator line, and the graph-preserving
singleton orbit.  It specifically checked the constant phase
`psi(bh+a(D_h phi(y)+phi(h)))`, the distinction between zero-intercept and
affine localization, all three THM-2337 word masks, and literal calculations
of `D_graph` and `E_sigma`.  Every executable test uses an explicit raising
check; ordinary and optimized runs byte-match the stored output, and the
audit independently reproduced the LF hashes above.  No row or masked
landing consequence is inferred.
