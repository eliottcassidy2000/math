---
id: THM-3678
title: "Russell-cylinder Q_dagger actual-ring J_0 lift"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For the
  degree-eight zero-second-and-fourth-debt fold Q_dagger from THM-3677, the
  zero-stable pair U=c,V=e has a rational actual-target-ring lift through
  J_0=1.  In the raw b^i c^j e^k restriction-degree filtration, cutoff 195
  is sharp: cutoff 192 has exact ranks 208 and 209, while cutoff 195 has
  exact operator and augmented rank 214 and supplies a deterministic RREF
  certificate.  This is only the first stable gate; no J_1,J_2, Keller-pair,
  or JC(2) conclusion is claimed.
source: kps-s194 / Q_dagger actual-lift scout, 2026-08-21
audit: >
  PASS -- kps-s195 independently checked the J_0 pullback equation, actual
  target-ring typing of the raw monomial columns, both exact rank jumps, the
  gcd-three cutoff minimality scope, deterministic free-zero RREF
  certificate, restriction/normal-derivative identities, retained values,
  and every pinned hash.  Normal and optimized companions returned all 238
  gates byte-identically.  The sharpness claim remains explicitly confined
  to the raw weighted packet filtration.  No correction was required.
depends_on:
  - THM-3677-russell-cylinder-degree-eight-fourth-debt-parabola
related:
  - THM-3642-russell-cylinder-zero-debt-actual-lift-and-fourth-stable-closure
script: 04-computation/jc2_russell_cylinder_qdagger_actual_j0_lift_thm3678.py
output: 05-knowledge/results/jc2_russell_cylinder_qdagger_actual_j0_lift_thm3678.out
script_sha256: 63f59e111e298652073bfeb7187eebaeb99c8d5e8def6efbdb5420b099d92d46
output_sha256: f62db26ffde4c167b63cbb984d9bc5a38c46833ce236e0b41ec17c67d6f0fa2b
hash_basis: raw LF bytes
---

# THM-3678 -- Q_dagger actual-ring J_0 lift

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The sharp
fourth-debt candidate from THM-3677 passes the first actual-ring gate.  The
statement is deliberately limited to `J_0`; the coupled stable equations have
not yet been certified.

All rings are over `C`, while the displayed certificate is rational.

## 1. Setup

Use the Russell cylinder and its fixed polynomial compiler

```text
R=C[b,c,e]/(c^2e-b(b+4)),

D=1+x^2q,
b=(D-1)(D+2)^2,       c=xD(D+2),       e=q(D+3).       (1)
```

The principal fold is

```text
Q_dagger=
 (22868x^8-89583x^6+2916x^5+123684x^4-5832x^3
  -63530x^2+2916x-2187)/2916.                          (2)
```

It is the `r=0` point of the zero-fourth-debt parabola in THM-3677.  Write
`gamma:R -> Q[x]` for restriction to `q=Q_dagger(x)` and `delta` for first
normal differentiation with respect to `q`, followed by the same restriction.
The three raw target generators have restriction degrees

```text
(deg gamma(b),deg gamma(c),deg gamma(e))=(30,21,18).    (3)
```

For the zero-stable choice

```text
U=c,                         V=e,                       (4)
```

and target coefficients `F_1,G_1 in R`, the constant stable coefficient is

```text
J_0=c' gamma(G_1)-gamma(F_1)e'.                        (5)
```

Here and below primes mean `x` derivatives after restriction.

## 2. Exact sharp rank transition

Let `M_N` be the ordered raw packet

```text
M_N={b^i c^j e^k: 30i+21j+18k <= N},                  (6)
```

listed in nested `(i,j,k)` order.  Form the rational coefficient matrix with
columns

```text
[c' gamma(m)]_(m in M_N),       [-e' gamma(m)]_(m in M_N), (7)
```

and augment it by the coefficient vector of `1`.  Exact rational elimination
gives

```text
 N    |M_N|   matrix       rank(operator)  rank(augmented)
192    173    213 x 346          208              209
195    179    216 x 358          214              214.    (8)
```

Thus cutoff `192` is impossible and cutoff `195` is feasible.  Moreover every
weight in `(3)` is divisible by three, so `M_194=M_192`; all lower packets are
subsets of `M_192`.  Consequently `195` is the **minimal cutoff in this raw
restriction-degree filtration**, not merely the first cutoff sampled.

The derivative gcd control is also exact:

```text
gcd(c',e')=1 in Q[x].                                  (9)
```

It is necessary but not sufficient for the target-ring lift; the rank jump in
`(8)` is the nontrivial part.

## 3. Frozen actual-ring certificate

At cutoff `195`, set all RREF free variables to zero.  Interpreting the two
coefficient blocks against the actual monomials `(6)` defines elements
`F_1,G_1 in R`, not merely arbitrary source polynomials.  Their exact source
restrictions satisfy

```text
c' gamma(G_1)-gamma(F_1)e'=1 in Q[x].                 (10)
```

The frozen certificate data are

```text
                    support   degree   SHA-256
F_1 target vector      107       --    587dc21d...068cf72
G_1 target vector      104       --    2aed39df...2cdca6c
gamma(F_1)              --      195    1180c242...59c36cb
gamma(G_1)              --      192    942601ed...b5bcc3
delta(F_1)              --      187    fca72977...e78aba
delta(G_1)              --      184    0e3b6b7a...58fa683. (11)
```

The companion prints and gates the full hashes.  It also verifies the retained
values

```text
gamma(F_1)(-1,0,1)=(0,0,0),
gamma(G_1)(-1,0,1)=(1/3,1/3,1/3).                    (12)
```

These values are useful controls for the next coupled solve, but `(10)` is the
actual global polynomial identity.

## 4. Boundary of the result

This theorem proves exactly one thing: `Q_dagger` does not die at the
`J_0=1` target-ring gate.  It does **not** provide `F_2,G_2,F_3,G_3`, does not
establish `J_1=J_2=0`, and does not produce a Keller pair.  The next decisive
test is the coupled `J_1,J_2` rank calculation with the coefficient vectors in
`(11)` frozen, because their normal derivatives enter that system.

## 5. Reproduction

```bash
python -B 04-computation/jc2_russell_cylinder_qdagger_actual_j0_lift_thm3678.py
python -O -B 04-computation/jc2_russell_cylinder_qdagger_actual_j0_lift_thm3678.py
```

Both runs return the stored transcript with `RESULT=PASS;gates=238`.
