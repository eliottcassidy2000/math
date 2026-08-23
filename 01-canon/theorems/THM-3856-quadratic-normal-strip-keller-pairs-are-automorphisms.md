---
id: THM-3856
title: "Quadratic normal-strip Keller pairs are automorphisms"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over a
  characteristic-zero field, every polynomial
  Keller pair in k[s,z] whose two transverse z-degrees are at most two is a
  polynomial automorphism.  After a constant target change, every genuinely
  quadratic pair is exactly C=b(s)+beta*z and
  A=rho*C^2+d*C-(lambda/beta)*s+a0.  Thus no such pair can realize the
  self-identifying Russell arm required by THM-3843.  This does not cover
  cubic, infinite formal, or rational normal expansions.
source: root / Russell higher-normal algebraization lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (higher_normal_lift and
  rule30_tariff_audit, 2026-08-23).  Both auditors independently derived the
  four Jacobian buckets, checked the zero-top, zero-linear, one-sided, and
  proportional-top cases, reduced the quadratic row by target shears, and
  recovered explicit polynomial inverses.  The independent exact companion
  additionally checks all minimal nodal Bezout rows, the conductor and
  self-collision, six bounded Groebner hostiles, positive linear/quadratic
  automorphisms, and the nonzero debt left by a Catalan truncation.  Normal
  and optimized runs byte-match both frozen LF transcripts.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3843-russell-arm-birational-immersion-and-forced-self-identification
  - THM-3846-formal-arm-darboux-lift-and-algebraization-gate
related:
  - THM-3795-r-independent-quadratic-normal-nodal-carriers-have-critical-points
  - THM-3812-nodal-arm-coefficient-second-normal-profile-nonentry
  - THM-3860-russell-higher-normal-rational-lifts-and-vertical-pole-barrier
script: 04-computation/jc2_quadratic_normal_strip_keller_thm3856.py
output: 05-knowledge/results/jc2_quadratic_normal_strip_keller_thm3856.out
script_sha256: c7fc82175a2daa9a2bc81df24c31c1cddb0f67cd2263ee26ea6257c02d6fb144
output_sha256: ce2cd8f92df66caeb6e58dbdaddb6edf8c4e90084552e0c707768d2c35626d9b
semantic_sha256: 02fce0fd2eea8d8ec75ee05bb8264d4bd82278d3f01db5c79ebf54e37ae3f82a
independent_audit_script: 04-computation/jc2_quadratic_normal_strip_keller_independent_audit_thm3856.py
independent_audit_output: 05-knowledge/results/jc2_quadratic_normal_strip_keller_independent_audit_thm3856.out
independent_audit_script_sha256: a38a5f450c8fe9744ec9953d0f39b388b3845ff33101cecac402f450487bb150
independent_audit_output_sha256: 27edc28aea3c9bfb28c2f0164e9fa6c020cc59ecd85a0ba45f0bad38442dc224
independent_semantic_sha256: e76ba373b48b84a69d24c2650ac5bd92fc697d39f7686928f9ab721a62ace210
hash_basis: raw LF bytes
---

# THM-3856 -- quadratic transverse depth is still triangular

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Let `k` be a
field of characteristic zero and let

```text
A,C in k[s,z],             deg_z A, deg_z C <=2,                 (1)
```

satisfy

```text
J_(z,s)(A,C)=A_z C_s-A_s C_z=lambda in k*.                      (2)
```

Then `(A,C):A2_(s,z)->A2` is a polynomial automorphism.  More precisely,
after a constant target `GL_2(k)` change, exactly one of the following holds.

1. The pair is linear in `z`, and a constant target combination is affine
   linear in `s`.
2. There are `rho,beta in k*`, `d,a0 in k`, and `b in k[s]` such that

   ```text
   C=b(s)+beta z,
   A=rho C^2+d C-(lambda/beta)s+a0.                              (3)
   ```

In the second case the inverse is displayed by

```text
s=(beta/lambda)(rho C^2+d C+a0-A),
z=(C-b(s))/beta.                                                (4)
```

Consequently, in the canonical Russell completed-arm coordinates of
THM-3846, no pair lying in the polynomial subring `k[s,z]` with both normal
degrees at most two can have the self-identifying arm forced by THM-3843.
A hypothetical Russell Darboux pair must either use polynomial normal degree
at least three or have a nonpolynomial (rational or infinite formal) strip
expansion.

## 1. The four exact Jacobian buckets

Write

```text
A=a+alpha z+u z^2,                 C=b+beta z+v z^2,             (5)
```

with all six coefficients in `k[s]`.  Equating the four coefficients of
`(2)` gives

```text
z^0: alpha b'-a' beta=lambda,                                  (6)
z^1: alpha beta'-alpha' beta+2(u b'-a'v)=0,                    (7)
z^2: alpha v'+2u beta'-2alpha'v-u' beta=0,                     (8)
z^3: 2(uv'-u'v)=0.                                             (9)
```

There is no quotient or completion loss here: these are identities in the
ordinary polynomial ring `k[s,z]`.

If `(u,v)!=(0,0)`, equation `(9)` says that the rational function `u/v`
has derivative zero whenever `v!=0`; in characteristic zero the constant
field of `k(s)` is `k`.  The zero-component cases give the same conclusion.
Thus

```text
(u,v)=h(s)(U,V),                 (U,V) in k^2\{0}.               (10)
```

A constant target `SL_2(k)` change sends `(U,V)` to `(1,0)`.  It preserves
`lambda`, so the genuinely quadratic case may be written

```text
A=a+alpha z+u z^2,                C=b+beta z,        u!=0.       (11)
```

Equations `(6)--(9)` become

```text
alpha b'-a'beta=lambda,                                      (12)
alpha beta'-alpha'beta+2u b'=0,                              (13)
2u beta'-u'beta=0.                                           (14)
```

## 2. The quadratic case has a polynomial inverse

If `beta=0`, equation `(12)` makes `alpha` and `b'` polynomial units, while
`(13)` then forces `u=0`.  This contradicts the genuine quadratic assumption
in `(11)`.  Hence `u beta!=0`.

Equation `(14)` is exactly

```text
(beta^2/u)'=0.                                                (15)
```

Therefore `u=rho beta^2` for some `rho in k*`.  Divide `(13)` by `beta^2`:

```text
(alpha/beta)'=2rho b'.                                       (16)
```

It follows that

```text
alpha=beta(2rho b+d),                         d in k.           (17)
```

Substitution in `(12)` yields the polynomial product identity

```text
beta (rho b^2+d b-a)'=lambda.                                (18)
```

Both factors must be units of `k[s]`.  Thus `beta in k*`, and integration
of `(18)` gives

```text
a=rho b^2+d b-(lambda/beta)s+a0.                              (19)
```

Equations `(11),(17),(19)` assemble to `(3)`, and direct substitution gives
the polynomial inverse `(4)`.

## 3. The linear boundary is also an automorphism

Now let `u=v=0`.  Equations `(6),(7)` reduce to

```text
alpha b'-a'beta=lambda,            alpha beta'-alpha'beta=0.  (20)
```

If `alpha=0`, the first identity makes `a'` and `beta` nonzero constants;
`A` recovers `s`, and then `C` recovers `z`.  If `alpha!=0`, the second
identity gives `beta=mu alpha` for some `mu in k`.  The first becomes

```text
alpha(b'-mu a')=lambda.                                     (21)
```

Hence `alpha in k*` and

```text
C-mu A=(lambda/alpha)s+kappa.                               (22)
```

Again the target combination recovers `s`, and `A=a(s)+alpha z` then
recovers `z`.  This proves the theorem in every degree-drop case.

## 4. Russell consequence and exact boundary

THM-3846 identifies the completed arm with

```text
Bhat=k[s][[z]],                           {z,s}=1.              (23)
```

The restriction of a polynomial automorphism of `A2_(s,z)` to the line
`z=0` is injective, because the ambient map itself is injective.  But
THM-3843 proves that the arm restriction of every hypothetical global
Darboux pair on the Russell pseudo-plane is noninjective.  Thus `(1)` cannot
hold for that completed pair.

The qualification "polynomial subring" is load-bearing.  A global element
of the Russell ring generally has an infinite `z`-series in `(23)`, and a
rational higher-normal correction need not lie in `k[s,z]`.  THM-3860 is the
positive hostile: it constructs a rational exact lift crossing the canonical
square-root gate.  The present theorem claims neither cubic-normal
nonexistence nor a result for rational or arbitrary formal pairs.  Positive
controls `(3)` exist for every polynomial `b(s)` and arbitrarily high
`s`-degree; the obstruction is transverse depth, not arm degree.  No global
Russell Darboux pair is constructed, and `JC(2)` remains open.
**QED.**

## Reproduction

```bash
python3 -B 04-computation/jc2_quadratic_normal_strip_keller_thm3856.py
python3 -B -O 04-computation/jc2_quadratic_normal_strip_keller_thm3856.py
python3 -B 04-computation/jc2_quadratic_normal_strip_keller_independent_audit_thm3856.py
python3 -B -O 04-computation/jc2_quadratic_normal_strip_keller_independent_audit_thm3856.py
```

The first two executions must byte-match
`05-knowledge/results/jc2_quadratic_normal_strip_keller_thm3856.out`; the
last two must byte-match the corresponding independent-audit output.
