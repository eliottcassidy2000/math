---
id: THM-4265
title: "W=0 reduced wall factor and transverse-Jacobian reduction"
status: >
  PROVED RELATIVE TO THM-4260 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  For every final degree-34/42 map class, the unavoidable t=+1 and t=-1
  factors of the reduced hidden coefficient denominator are simple. Hence
  any coherent one-parameter off-W=0 lift has unique formal root sections
  through the wall points, whose first-order separation is exactly measured
  by a two-by-two transverse Jacobian. No lift, transverse derivative,
  off-fibre exclusion, M=12 entry, JC(2), or DC(2) is proved.
source: root/cross-frontier-bridge/2026-08-26
depends_on:
  - THM-4260-w0-canonical-node-reciprocal-denominator-attachment-exclusion
related:
  - THM-4255-specialization-kernel-and-transverse-hasse-jet-repair
  - THM-4263-moving-multigraph-filtered-jet-and-finite-factor-density-transport
  - THM-4264-w0-visible-incidence-two-edge-attachment-observer
script: 04-computation/jc23_w0_reduced_wall_factor_thm4265.py
output: 05-knowledge/results/jc23_w0_reduced_wall_factor_thm4265.out
script_sha256: 41849ca35164391ec08820b63e9d4ae89de6358696501c7e51508ca5c4e7b54b
output_sha256: 9e5225503e81ea71a0ea2d222c60949f809c6910936c209af28c3ac968fa4ef0
independent_script: 04-computation/jc23_w0_reduced_wall_factor_independent_audit_thm4265.py
independent_output: 05-knowledge/results/jc23_w0_reduced_wall_factor_independent_audit_thm4265.out
independent_script_sha256: 75027020442d59e6168ff4f82941f043a417ceae1020871974901a965cb85d51
independent_output_sha256: f35ba8fb733df7d15519f408ebc53fd87e4f18c3a3f5f7ac3e72f30ff317fe4a
hash_basis: raw LF bytes
audit: >
  PASS/ACCEPT. A clean-room SymPy/GF(t) reconstruction recovered all 280
  classes independently of the candidate checker and verified simple wall
  roots, odd derivative parity, the reciprocal gcd, and the exact degree
  census over q=397 and q=577. Normal, optimized, and fixed-hash-seed runs
  byte-match; extra-factor and missing-root hostile controls are detected.
---

# THM-4265 -- `W=0` reduced wall factor and transverse-Jacobian reduction

**PROVED RELATIVE TO THM-4260 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
THE OFF-`W=0` PROBLEM REMAINS OPEN.**

## 1. Statement

Use THM-4260's final `W=0` class universe and notation. It has `176`
degree-`34` and `104` degree-`42` source-target classes. For a class `m`, put

```text
ell=(1-iota)m,                 q(ell)=12K,
d_m(t)=reduced coefficient denominator of X_ell/x,
d_m^*(t)=t^deg(d_m) d_m(-1/t).                         (1)
```

THM-4260 proves in characteristic zero that

```text
d_m(t)=tD_m(t^2),       deg d_m=4K-1,       ord_t d_m=1,
gcd(d_m,d_m^*)=t^2-1.                                 (2)
```

Then, for every one of the `280` classes,

```text
(t-1)||(d_m),             (t+1)||(d_m),                (3)
```

and the same holds for `d_m^*`. Thus the two unavoidable geometric wall
points are reduced in the denominator divisor; the gcd in `(2)` is not the
visible part of a higher-order collision.

## 2. Characteristic-zero proof and exact corroboration

Write `mu_+` and `mu_-` for the multiplicities of `+1` and `-1` in `d_m`.
The form `d_m=tD_m(t^2)` gives `mu_+=mu_-=:mu`. The reciprocal involution
`t -> -1/t` swaps the two points, so `d_m^*` has multiplicity `mu` at each
one as well. Therefore the multiplicity of either root in
`gcd(d_m,d_m^*)` is `mu`. The last equality of `(2)` says that this gcd is
the squarefree polynomial `t^2-1`; hence `mu=1`. This proves `(3)` directly
from THM-4260.

Two exact derivative ledgers independently corroborate the conclusion. The
primary checker reconstructs THM-4260's class universe and records

```text
d_m(1)=d_m(-1)=0,       d_m'(1)=d_m'(-1)!=0            (4)
```

for all `280` classes over each of `q=397` and `q=577`. The derivative
equality is the expected parity control because `d_m` is odd. Both fields
give the degree census

```text
degree(d_m):classes =
11:8, 19:36, 23:24, 27:64, 35:52, 39:72, 43:24.        (5)
```

The complete primary wall-ledger hashes are

```text
q=397: 58809890f97d24628406ad9f8d28b31f61a7f45baae34c6c7c0f61ec67d49571
q=577: e50dc1d8b22eba1d9b12e9c5a984d4740b52c3f99c5765bf9bec91aee81fb9da.
```

A clean-room referee rebuilt the `176+104` classes from the frozen THM-4260
SymPy/GF(t) compiler, independently repeated every elliptic-curve addition,
and recovered `(4)--(5)` together with `gcd(d_m,d_m^*)=t^2-1`. Its class
ledger hashes are

```text
q=397: e8cf8fa7c5fd6f3632015e1455a1d7fb66e70e345d849a482a5ffdbb505d8d2b
q=577: 9985416e7f14c9f23e0baa0efe3464cff6ec0b3fbc64c84707259181f24ca09c.
```

These computations are hostile controls and reproducible corroboration; the
load-bearing reducedness proof is the characteristic-zero argument above.

## 3. Formal local consequence

The following conditional consequence explains why reducedness matters
without assuming that the required geometry exists. Suppose a class and its
reciprocal condition admit coherent normalized lifts over a characteristic-
zero field,

```text
D(W,t), E(W,t) in K[[W,t-t_0]],       t_0 in {1,-1},  (6)
```

whose `W=0` fibres are `d_m,d_m^*` up to units. By `(3)`, both `t`
derivatives at `(0,t_0)` are nonzero. The formal implicit-function theorem
therefore gives unique root sections

```text
t=a(W) for D=0,              t=b(W) for E=0.           (7)
```

Set

```text
J_m(t_0)=det [[partial_t D, partial_W D],
              [partial_t E, partial_W E]]_(0,t_0).     (8)
```

Implicit differentiation gives the exact identity

```text
a'(0)-b'(0)=J_m(t_0)/(D_t(0,t_0) E_t(0,t_0)).          (9)
```

Thus `J_m(t_0)!=0` is precisely first-order transverse separation: the two
divisors have only an isolated closed-point intersection there. If `J_m`
vanishes, the remaining local question is the finite or infinite
`W`-adic order of `a(W)-b(W)`; higher transverse order may still be needed.
No additional `t`-jet is needed merely to resolve a multiple special-fibre
root.

## 4. Scope and reproduction

This theorem does not construct `D,E`, prove that a `W=0` Hom class extends,
calculate a `W` derivative, or exclude either wall off the special fibre.
The other coefficient walls, the hidden-Hom neighbourhood, exact `M=12`
seam entry, `JC(2)`, and `DC(2)` remain open. THM-4255 and THM-4263 are
specialization firewalls and observer templates, not dependencies: any
application must retain the specialization kernel, transverse class/lattice
sidecar, and `W`-order.

From the repository root:

```bash
PYTHONHASHSEED=0 python3 -B \
  04-computation/jc23_w0_reduced_wall_factor_thm4265.py --repo .

PYTHONHASHSEED=0 python3 -B \
  04-computation/jc23_w0_reduced_wall_factor_independent_audit_thm4265.py \
  --repo .
```

Normal, optimized, and fixed-hash-seed runs are byte-identical to their
frozen outputs. **QED.**
