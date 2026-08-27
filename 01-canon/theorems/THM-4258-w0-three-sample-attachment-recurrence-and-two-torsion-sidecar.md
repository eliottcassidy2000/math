---
id: THM-4258
title: "W=0 three-sample attachment recurrence and two-torsion sidecar"
status: >
  PROVED COROLLARY TO THM-4249/4253 + VERIFIED-EXACT. On the residual
  a_u=0 Hom-lattice lane, three consecutive attachment differences vanish
  iff all twelve vanish. Thus a uniform observer table for the 1,512 live
  map-ratio incidences uses 4,536 group-value rows rather than the naive
  18,144. Observer-map degrees are at most 120 in degree 34 and 150 in
  the post-THM-4253 degree-42 shell. Conditional on both visible and hidden
  projectors already collapsing, the residual differences form at most 64
  F4-valued confluent two-torsion patterns. No incidence is excluded; W=0,
  M=12, seam entry, JC(2), and DC(2) remain OPEN.
source: codex-padic-density-cartier-20260826
depends_on:
  - THM-4249-w0-cyclic-projector-missing-eigenline-attachment-squeeze
  - THM-4253-w0-degree-forty-two-norm-three-profile-exclusion
related:
  - THM-4241-w0-hidden-cm-saturation-and-visible-hidden-index-four-gluing
  - THM-4255-specialization-kernel-and-transverse-hasse-jet-repair
mistake_firewall:
  - MISTAKE-528
script: 04-computation/jc23_w0_three_sample_attachment_observer_thm4258.py
output: 05-knowledge/results/jc23_w0_three_sample_attachment_observer_thm4258.out
script_sha256: be5a1c821e85011b3226313f7c0f0461c10afc1609ff8e867be697e8f1d8d3d3
output_sha256: 9042b3a1017ac2205b9a76ff89b2624d4109cf757a530442bf6efe6c5e1fbfb0
hash_basis: raw LF bytes
audit: >
  PASS. The companion imports the dependency-free THM-4249 lattice auditor,
  reconstructs the 4,224/2,496 post-sieve vectors and 176/104 size-24
  orbits, verifies the cubic annihilator on every vector, and freezes both
  map-orbit and incidence-weighted observer-degree histograms. A separate
  standard-library F4 path exhausts all 64 recurrence patterns, proves
  three-prefix injectivity, exhibits four-to-one two-prefix fibres, and
  checks the period profile. Normal, optimized, and fixed-hash-seed outputs
  byte-match.
---

# THM-4258 -- `W=0` three-sample attachment recurrence

**PROVED COROLLARY TO THM-4249/4253 + VERIFIED-EXACT. NO CANDIDATE IS
EXCLUDED BY THIS THEOREM.**

## 1. Exact recurrence and complete three-sample observer

Use the `W=0` gate, elliptic target `E_0`, Eisenstein ring
`O=Z[omega]`, order-twelve source action `T`, attachments
`Q_j=tau^j Q_0`, and normalized basis `[u,f,g,h]` of THM-4249. Put

```text
M_0={m=a_u u+b f+c g+d h : a_u=0},
D_j=[Q_(j+1)-Q_j] in J(C_0),
delta_j=m(D_j) in E_0,                   j mod 12.       (1)
```

Attachment collapse for `m` is exactly `delta_j=O` for all twelve `j`.
For any base point `Q_*`, set

```text
m_tilde(Q)=m([Q-Q_*]).                                  (1a)
```

Changing `Q_*` translates all values by one target point, so equality of
the twelve `m_tilde(Q_j)` is intrinsic.

THM-4249 proves

```text
P_u=-(T-omega^2)(T^2+omega),       P_u(m)=a_u u.        (2)
```

Hence every `m in M_0` satisfies

```text
P_0(T)m=0,
P_0(X)=(X-omega^2)(X^2+omega)
      =X^3-omega^2 X^2+omega X-1.                       (3)
```

Since `(T^k m)(D_j)=m(D_(j+k))`, evaluation of `(3)` gives the exact
elliptic-group recurrence

```text
delta_(j+3)-[omega^2]delta_(j+2)
             +[omega]delta_(j+1)-delta_j=O.             (4)
```

Define the translation-normalized difference map on the curve by

```text
F_m(Q)=m([tau Q-Q]).                                    (5)
```

Its induced Hom class is `(T-1)m`; this definition removes the arbitrary
translation constant of a curve-map representative. Define
`T^jF_m=F_m composed_with tau^j`. Then

```text
delta_j=(T^j F_m)(Q_0).                                 (6)
```

> **Theorem 1 (three samples are complete).** For every `m in M_0`,
>
> ```text
> m_tilde(Q_0)=m_tilde(Q_1)=...=m_tilde(Q_11)
> iff F_m(Q_0)=(TF_m)(Q_0)=(T^2F_m)(Q_0)=O.             (7)
> ```
>
> Equivalently, the `O`-linear observer
>
> ```text
> Obs_3:M_0 -> E_0^3,
> Obs_3(m)=(delta_0,delta_1,delta_2)                     (8)
> ```
>
> has kernel exactly the twelve-attachment collapse submodule.

### Proof

If all attachments have the same image, all three displayed differences
vanish. Conversely, three consecutive zero values in `(4)` force the next
one to vanish. Iterating from `j=0` gives `delta_3,...,delta_11=O`.
Equation `(6)` identifies the three values with `(8)`. QED.

The same proof works from any three consecutive edges. It uses an intrinsic
binary relation—successive attachment equality—and retains the full target
group value, so there are no artificial tournament ties or quotient losses.

## 2. Exact observer degrees and workload

Let `q` be the map-degree Hermitian form, applied to the induced Hom class
`(T-1)m`. Precomposition by `T` and target unit scaling preserve `q`, so
`q(F_m)` is constant on each source-target
orbit. Exact enumeration after every current filter gives

```text
degree 34, 176 map orbits:
{30:8,36:4,42:10,48:6,54:7,60:8,66:21,72:6,
 78:19,84:13,90:16,96:10,102:20,108:9,114:17,120:2};  (9)

degree 42 after THM-4253, 104 map orbits:
{36:1,42:2,54:8,66:6,72:6,78:3,84:4,90:14,96:4,
 102:10,108:9,114:4,120:4,126:12,132:4,138:10,
 144:2,150:1}.                                        (10)
```

Thus the exact observer-degree ceilings are

```text
q(F_m)<=120 in degree 34,          q(F_m)<=150 in degree 42. (11)
```

There is a conceptual bound before enumeration. On `M_0`, the unitary
operator `T` has eigenvalue `omega^2` and the two roots of `X^2=-omega`.
Their squared chord lengths from `1` are `3`, `2-sqrt(3)`, and
`2+sqrt(3)`. Therefore

```text
q((T-1)m)<=(2+sqrt(3))q(m).                           (11a)
```

Pulling the Hermitian Gram back by `T-1` on `[f,g,h]` gives

```text
[[18,          -20-10omega, -10+4omega],
 [-10+10omega, 18,            4-10omega],
 [-14-4omega,  14+10omega,   12]].                   (11b)
```

The diagonals are divisible by six. For each upper off-diagonal `a`, both
`Tr(a)` and `Tr(a*omega)` are divisible by six, so

```text
q(F_m) in 6Z.                                         (11c)
```

Equations `(11a)--(11c)` already give `q(F_m)<=126,156`; exact enumeration
sharpens them to `120,150`. The cruder quadratic triangle bound alone gives
`136,168`.

Weighting each map orbit by its current marked-ratio incidences gives
`864+648=1,512` live incidences. A uniform direct observer table from `(7)`
has

```text
1,512*12=18,144  to  1,512*3=4,536,                    (12)
```

a reduction of `13,608` group-value rows relative to the naive twelve-edge
table. This is sufficient, not a global minimality theorem: telescoping
already makes one of twelve raw differences dependent, and shared symbolic
work may compress further. The full incidence-weighted degree histograms are
frozen in the companion output. No incidence is excluded.

## 3. Conditional two-torsion confluent sidecar

For `m in M_0`, THM-4249's integral reconstruction reads

```text
2m=dv+H,                 H=2ell in O f+O g.             (13)
```

Suppose, as an additional gate, that both projected maps `dv` and `H`
already collapse the attachment orbit. Then `(13)` gives

```text
[2]delta_j=O for every j,       so delta_j in E_0[2].  (14)
```

This hypothesis is essential. The visible torsion incidence enforces the
`dv` condition, but the hidden `H` condition is not yet part of the current
`1,512`-incidence ledger.

As an `O`-module, `E_0[2]` is free of rank one over `O/(2)=F_4`. Put
`alpha=omega^2`. Modulo two,

```text
P_0(S)=(S-alpha)(S^2+omega)=(S-alpha)^3.               (15)
```

Twist

```text
epsilon_j=alpha^(-j)delta_j=omega^j delta_j.            (16)
```

With `Delta epsilon_j=epsilon_(j+1)-epsilon_j`, equation `(4)` is equivalent
to

```text
Delta^3 epsilon=0.                                     (17)
```

Every such sequence has the unique Newton form

```text
epsilon_j=c_0+binom(j,1)c_1+binom(j,2)c_2,
c_0,c_1,c_2 in E_0[2] ~= F_4.                          (18)
```

Here the binomial coefficients are reduced modulo two.

Consequently the abstract recurrence module has exactly `4^3=64` patterns,
all satisfying `sum_(j=0)^11 delta_j=O`. Geometry realizes an unknown subset,
hence at most 64. The three initial values recover `(c_0,c_1,c_2)` exactly.
The exhaustive period profile is

```text
{period 1:1, period 3:3, period 6:12, period 12:48}.    (19)
```

Three samples are ambient-sharp. For nonzero `P in E_0[2]`,

```text
epsilon_j=binom(j,2)P,
delta_j=alpha^j binom(j,2)P                             (20)
```

satisfies `(4)`, has `delta_0=delta_1=O`, and has
`delta_2!=O`. This proves that two observations cannot classify the abstract
recurrence module. It does **not** prove that `(20)` is realized by one of
the geometric Hom maps, nor that all 64 patterns occur geometrically.

## 4. The exact THM-4255 connection

THM-4255 studies a chosen-section quotient with normal kernel `(u-g)` and
repairs the lost direction by transverse Hasse jets. Here no formal-section
map occurs. Instead, reduction modulo two makes the cyclic operator acquire
the confluent factor `(S-alpha)^3`; one scalar eigenvalue forgets a length-
three generalized eigenline, and three consecutive observations recover it.

The typed correspondence is therefore

```text
normal power (u-g)^3          <-> confluent factor (S-alpha)^3,
Hasse coefficients           <-> forward-difference coefficients,
three transverse jets        <-> three consecutive attachment samples. (21)
```

It preserves the principle “retain the nilpotent sidecar before quotienting”
but not the underlying object: formal arcs and elliptic attachment orbits are
different categories. No p-adic-zeta or Jacobian theorem transfers through
the analogy alone.

## 5. Concrete normalization not supplied here

Equations `(7)--(12)` reduce the task to three evaluations, but this theorem
does not supply them in concrete normalized coordinates. THM-4241 gives the
explicit degree-four map `H_lambda`, then says only that scaling by an
`O`-unit and changing representatives by `D=V direct-sum L` normalizes its
glue class to

```text
2h=v+omega^2 f+g.                                      (22)
```

THM-4241 itself does not freeze the change-of-basis dictionary from
`H_lambda` to the specific `[u,f,g,h]` coefficients used by the residual
census. A separate finite-coordinate certificate can supply it. Until that
certificate is canonized, halving the explicit value of
`v+omega^2f+g` leaves the two-torsion ambiguity measured by Section 3.

The cheapest decisive bridge is now finite and explicit:

1. compute the visible and hidden projections of `H_lambda` and express
   `2H_lambda-v_0=A f+B g` for the correct visible unit multiple `v_0`;
2. identify the target unit and `D`-translate sending that concrete glue
   generator to `(22)`;
3. freeze the resulting basis-change and `T` matrices; and
4. evaluate only `F_m(Q_0),(TF_m)(Q_0),(T^2F_m)(Q_0)` using the
   translation-normalized definition `(5)` for the `1,512`
   incidence rows, with the existing ratio-overlap ledger retained.

## 6. Scope and reproduction

This theorem proves no marked-ratio set empty, removes no map or incidence,
does not close `W=0` or exact `M=12`, and does not prove seam entry,
`JC(2)`, or `DC(2)`.

```bash
python3 -B 04-computation/jc23_w0_three_sample_attachment_observer_thm4258.py
python3 -B -O 04-computation/jc23_w0_three_sample_attachment_observer_thm4258.py
PYTHONHASHSEED=4258 python3 -B \
  04-computation/jc23_w0_three_sample_attachment_observer_thm4258.py
```

All three streams byte-match the frozen output. **QED.**
