---
id: THM-4253
title: "W=0 degree-forty-two norm-three profile exclusion"
status: >
  PROVED COROLLARY TO THM-4249 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  The degree-42 residual profile (N(d),K)=(3,13), equivalently
  (q(v),q(ell))=(12,156), has 672 vectors in 28 size-24 source-target
  symmetry orbits. Its visible torsion envelope has the unique gate-admissible
  ratio 1/3, which THM-4249 proves is not in S_42. Hence the whole profile is
  impossible. The degree-42 frontier sharpens from 3,168 vectors/132 orbits
  to 2,496 vectors/104 orbits. Its already-final 648 incidence count is
  unchanged because the 28 profile incidences were among the 132 common
  ratio-1/3 incidences already removed by THM-4249. W=0, M=12, entry, JC(2),
  and DC(2) remain OPEN.
source: root/jc-planar/2026-08-26
depends_on:
  - THM-4249-w0-cyclic-projector-missing-eigenline-attachment-squeeze
related:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4241-w0-hidden-cm-saturation-and-visible-hidden-index-four-gluing
  - THM-4247-w0-involution-degree-twelve-attachment-exclusion
mistake_firewall:
  - MISTAKE-521
scripts:
  - 04-computation/jc23_w0_degree42_norm3_profile_exclusion_thm4253.py
  - 04-computation/jc23_w0_degree42_norm3_profile_exclusion_independent_audit_thm4253.py
outputs:
  - 05-knowledge/results/jc23_w0_degree42_norm3_profile_exclusion_thm4253.out
  - 05-knowledge/results/jc23_w0_degree42_norm3_profile_exclusion_independent_audit_thm4253.out
script_sha256:
  - b9cf49556a8f4229a02d8e4e0d6300172596286ce9dbf185d535d22d952daf5e
  - 6196a7c1a8cb674e32266a213a624da2c901032b6724f97f422f36cf87ed236d
output_sha256:
  - 09714e2d7f846582dcd86cbcbf716ac6a403de0e8e02dc8c5185415f25ee3f35
  - 7d74814f8d1f706448572dd234e9329cafd6ddb7c1ab2a3631d364c6c3dd9154
hash_basis: raw LF bytes
audit: >
  PASS/ACCEPT. The primary reconstructs the THM-4249 residual and symmetry
  profile, proves the ideal/kernel corollary, audits the incidence overlap,
  and independently performs 2,688 hidden-origin tests in each of four good
  fields. The clean-room path imports no maintained computation and directly
  rebuilds all 3,168 residual vectors, 132 size-24 orbits, profile counts,
  Eisenstein ideal arithmetic, torsion-orbit sizes, and incidence ledger.
  Normal, optimized, and fixed-hash-seed streams byte-match frozen outputs.
---

# THM-4253 -- degree-forty-two norm-three profile exclusion

**PROVED COROLLARY TO THM-4249 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement

Use the `W=0` notation and integral cyclic-projector residual of THM-4249.
Every degree-`42` candidate has zero `u` coordinate and a profile

```text
(N(d),K),             q(v)=4N(d),             q(ell)=12K.   (1)
```

> **Theorem.** No attachment-collapse candidate has profile
>
> ```text
> (N(d),K)=(3,13),     equivalently (q(v),q(ell))=(12,156). (2)
> ```
>
> This removes exactly `672` raw vectors, or `28` size-`24` orbits. The
> degree-`42` residual therefore becomes
>
> ```text
> 2496 vectors = 104 source-target symmetry orbits.                 (3)
> ```

The torsion-ratio envelope still has cardinality `34`, and the surviving
degree-`42` map-ratio incidence count remains `648`. The latter is an overlap
statement, not a failure to subtract.

## 2. The norm-three kernel has one admissible ratio

Put

```text
O=Z[omega],                   pi=omega^2-1=-2-omega,
N(pi)=3,                      pi^2=-3omega^2.          (4)
```

The Eisenstein ring is Euclidean and `pi` is the unique prime above `3`.
Every element `d` with `N(d)=3` is a unit multiple of `pi`; hence

```text
(d*pi)=(3),                   ker[d*pi]=E0[3].         (5)
```

THM-4249 proves that the visible point `P=v(Q_0)` of a collapsing map lies in
`ker[d*pi]`, and that the marked ratio is

```text
R=U/Z=-1/(X(P)^3+1).                                  (6)
```

The third division polynomial on `E0:Y^2=X^3+1` is

```text
psi_3(X)=3X(X^3+4).                                   (7)
```

The target-unit action splits `E0[3]` into orbits of sizes `1,2,6`:

- `O` is nonaffine and corresponds to a vanished endpoint;
- the two nonzero `pi`-torsion points have `X=0`, so `R=-1`, the wall
  `U+Z=0`; and
- the six remaining points have `X^3=-4`, so `R=1/3`.

Thus the profile `(2)` has the singleton gate-admissible envelope

```text
{1/3}.                                                 (8)
```

THM-4249 proves `1/3 notin S_42`. Equations `(5)--(8)` therefore exclude the
whole profile.

## 3. Exact vector and orbit ledger

The clean-room audit enumerates the full `a_u=0` degree-`42` lattice and
applies exactly the proved THM-4249 necessary filters. It obtains

```text
(N(d),K)   vectors   size-24 orbits
(3,13)       672          28
(9,11)       576          24
(12,10)      864          36
(21,7)       768          32
(27,5)       288          12
total       3168         132.                         (9)
```

Deleting the first line gives `(3)`. The generated action is the actual
`mu_6 x <T>` action on the full glued basis; no rational eigenspace quotient
or assumed freeness enters the count.

## 4. Incidence overlap

For one map orbit, the raw admissible ratio counts in `E0[d*pi]/mu_6` are

```text
N(d):      3   9   12   21   27
ratios:    1   4    5   10   13.                     (10)
```

Before THM-4249's ratio-`1/3` exclusion these contribute `780` incidences.
The `28` norm-three orbits contribute exactly `28` of them, all at `1/3`.
Hence two equivalent ledgers are

```text
780-132=648,
(780-28)-(132-28)=752-104=648.                       (11)
```

The profile deletion must not be subtracted again from `648`: its incidences
were already in the common-ratio deletion. The refined frontier is therefore

```text
degree 34: 4224 vectors, 176 orbits, 864 incidences;
degree 42: 2496 vectors, 104 orbits, 648 incidences.  (12)
```

## 5. Independent hostile control

At `R=1/3`, the attachment parameter satisfies

```text
t^4-14t^2+1=0.                                        (13)
```

The primary certificate enumerates all `672` compatible hidden projections
and all four roots of `(13)` in each good THM-4247 field

```text
(q,z,p,s)=(313,29,135,21),(349,24,246,28),
          (373,69,297,33),(397,157,161,27).           (14)
```

It performs `672*4=2688` exact hidden-origin tests per field and finds zero
hits, while rechecking the coefficient-field, source, gate, separability, and
simple-root conditions. This is an independent contrapositive control on the
already-proved THM-4249 ratio exclusion, not a modular zero-to-characteristic-
zero inference.

## 6. Scope and reproduction

The remaining `104` degree-`42` classes, all `648` degree-`42` incidences,
the `176/864` degree-`34` class/incidence frontier, emptiness of `S_34,S_42`,
`W=0`, exact `M=12`, entry, `JC(2)`, and `DC(2)` remain open.

```bash
python3 -B 04-computation/jc23_w0_degree42_norm3_profile_exclusion_thm4253.py
python3 -B -O 04-computation/jc23_w0_degree42_norm3_profile_exclusion_thm4253.py
python3 -B 04-computation/jc23_w0_degree42_norm3_profile_exclusion_independent_audit_thm4253.py
python3 -B -O 04-computation/jc23_w0_degree42_norm3_profile_exclusion_independent_audit_thm4253.py
```

Normal, optimized, and fixed-seed streams byte-match. **QED.**
