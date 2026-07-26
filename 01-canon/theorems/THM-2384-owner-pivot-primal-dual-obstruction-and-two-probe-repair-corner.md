---
id: THM-2384
title: "Owner-pivot primal-dual obstruction and two-probe repair corner"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In every THM-2309
  owner-aligned pivot, the balanced omitted-unit/helper relation is a
  pure nonzero primal class in K/L but is not a target covector in
  L^perp/<w>; an ungrafted star row detects it, and adding a scalar
  gauge multiple cannot repair the defect. The true pure-target dual
  is the target-blocker/helper dipole. Independently, the false
  two-factor shift still gives an exact positive duplicate-probe
  corner: for a deepest-safe, guard/unit-failed weight, some mixed
  coefficient is at least rho/2704 in the ordinary case and
  7rho/24336 in the guard case, with an exact energy floor. This is
  not a lawful owner-target current. No scalar-row exclusion, ledger
  decrement, or LRC(14) consequence follows.
source: codex-2026-07-25-owner-pivot-primal-dual-obstruction
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
related:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2379-anchored-guard-unit-deletion-factor-repair-current
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
script: 04-computation/lrc14_owner_pivot_primal_dual_repair_thm2384.py
output: 05-knowledge/results/lrc14_owner_pivot_primal_dual_repair_thm2384.out
script_sha256: 5d1eab37a5df1218d6fcdd1f36ac19d728ec5d2234133dec10f36a3dda324e76
output_sha256: 5927fa9ab1484a32c0fa801310d42cb3d4fd5d4ba41af44e6cace48ec06615aa
hash_basis: working-tree bytes (LF)
---

# THM-2384 -- a pure relation class is not its target covector

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The most tempting repair of THM-2379 is type-incorrect. A balanced
two-coordinate relation can represent a pure target in the **primal**
quotient, while the same coordinate vector fails to define any action
on that quotient:

```text
pure class in K/L
  does not imply
lawful co-shift in (K/L)^*=L^perp/<w>.                         (1)
```

The failure is universal and is detected by one ungrafted owner-pivot
row. The analytic two-probe table that suggested the false
identification remains useful, but only as a duplicate-probe
coefficient with no canonical target meaning.

## 1. The owner-pivot spaces have opposite variance

Work over `F_p` with `p=13`. Use THM-2309's notation

```text
E=U disjoint_union {j,a,b},             |U|=6,

wbar=(w_e mod p)_(e in E),

w_u!=0                    for u in U,

w_j=w_a=w_b=0.                                      (2)
```

Fix an omitted unit `u_0`, and put

```text
P={j} union (U minus {u_0}).
```

Choose distinct graft helpers `k_a,k_b in U\{u_0}` for the two target
blockers. Let `L` be the row space of the resulting owner-aligned pivot
and let

```text
K=wbar^perp.
```

THM-2309 proves

```text
K=L direct_sum span(e_a,e_b).                       (3)
```

Thus `K/L` is the two-dimensional space of **primal relation
classes**. Its character group is instead

```text
(K/L)^*=L^perp/span(wbar).                          (4)
```

The quotient in (3) and the quotient in (4) are paired dual spaces.
There is no canonical identification between their displayed
coordinate vectors.

## 2. The balanced helper relation is pure primal

Put

```text
v=w_(u_0),             k=w_(k_a),                   (5)

h=k e_(u_0)-v e_(k_a).
```

All scalars in (5) are nonzero. The vector `h` lies in `K`. The
`a`-grafted owner-pivot row is

```text
h+t_a in L,

t_a=-v e_a                                               mod p. (6)
```

Consequently,

```text
[h]=v[e_a]                         in K/L,            (7)

[k^(-1)h]=(v/k)[e_a].
```

This is a genuine pure nonzero primal target class.

It is not a target co-shift. Normalize

```text
g=k^(-1)h=e_(u_0)-(v/k)e_(k_a).                    (8)
```

Choose any unit label

```text
k' in U\{u_0,k_a,k_b}.
```

Its ungrafted star row is

```text
ell=w_(k')e_(u_0)-v e_(k') in L.                   (9)
```

Then

```text
g.ell=w_(k')!=0.                                    (10)
```

Hence

```text
g notin L^perp.                                     (11)
```

The obstruction is gauge-invariant. Since every owner-pivot row lies
in `K`,

```text
wbar.ell=0,
```

so

```text
(g+alpha wbar).ell=g.ell!=0                         (12)
```

for every `alpha in F_p`. No representative modulo `span(wbar)` repairs
the proposed two-coordinate action.

The same witness shows that shifting the omitted failed factor alone is
not lawful:

```text
e_(u_0).ell=w_(k')!=0.                              (13)
```

## 3. The true pure-target dual

The actual `a`-target covector is

```text
eta_a=e_a-e_(k_a).                                  (14)
```

Indeed, `eta_a` pairs to zero with every ungrafted row, with the
`b`-grafted row, and with the `a`-grafted row because that row has equal
coefficient `-v` in the `a` and `k_a` columns. Thus

```text
eta_a in L^perp,
```

and its class modulo `span(wbar)` is the pure `a` character in (4).
Similarly,

```text
eta_b=e_b-e_(k_b).                                  (15)
```

Equations (14)--(15) are the owner-pivot dipoles of THM-2350. They
shift a target blocker and its graft helper. They do **not** shift the
omitted failed factor. Any lawful dual representative with a nonzero
`u_0` coordinate must add the pivot-coordinate corrections forced by
the ungrafted equations (9), so the two-factor counting argument below
cannot be transferred unchanged.

## 4. The surviving analytic two-probe corner

The bad target interpretation does not invalidate the positive
analysis which led to it. On the circle `T=R/Z`, put

```text
d_L(y)=1_(||y||<L/14),        u_L=1-d_L,

L in {1,2},                   p=13.                 (16)
```

All identities are almost everywhere with strict-open endpoints. Let
`c,v,k` be positive integers, let `lambda in F_p^*`, and let
`omega>=0` be integrable with mass

```text
rho=int_T omega(x)dx>0
```

and support

```text
support(omega)
 subset {d_1(cx)=0} intersection {d_L(vx)=1}.       (17)
```

Define the nonnegative duplicate-probe table

```text
J(r,t)
 =int_T omega(x)
    d_1(cx-r/p)
    u_L(vx-t/p)
    u_1(kx+lambda t/p) dx,             r,t in F_p. (18)
```

The anchors (17) give two zero axes:

```text
J(0,t)=0                  for every t,

J(r,0)=0                  for every r.              (19)
```

The first is deepest safety; the second is the failed guard/unit
factor at its original shift.

For every phase `x`, the exact shift laws give

```text
sum_r d_1(cx-r/p)=2-d_1(pcx)>=1,                    (20)

sum_t u_L(vx-t/p)u_1(kx+lambda t/p)

 =11-2L+d_L(pvx)+d_1(pkx)
   +sum_t d_L(vx-t/p)d_1(kx+lambda t/p)

 >=11-2L.                                           (21)
```

The only use of `lambda!=0` is that multiplication by `lambda`
permutes the thirteen helper shifts. Multiplying (20)--(21) pointwise
and integrating gives

```text
sum_(r,t)J(r,t)>=(11-2L)rho.                        (22)
```

## 5. Exact coefficient and energy floors

Use the normalized transform

```text
Jhat(a,b)
 =1/p^2 sum_(r,t)J(r,t)zeta^(ar+bt),

zeta=exp(2*pi*i/p).                                 (23)
```

Finite character orthogonality and the two zero axes in (19) give

```text
sum_(a!=0,b!=0)Jhat(a,b)
 =1/p^2 sum_(r,t)J(r,t)
 >=(11-2L)rho/p^2.                                 (24)
```

There are `(p-1)^2=144` mixed colours. Some `a,b!=0` therefore satisfies

```text
Re Jhat(a,b)
 >=(11-2L)rho/(p^2(p-1)^2).                        (25)
```

At `p=13`,

```text
failed ordinary unit, L=1:
  Re Jhat(a,b)>=9rho/24336=rho/2704;

failed guard, L=2:
  Re Jhat(a,b)>=7rho/24336.                         (26)
```

Cauchy--Schwarz also gives

```text
sum_(a!=0,b!=0)|Jhat(a,b)|^2
 >=(11-2L)^2 rho^2/(p^4(p-1)^2).                   (27)
```

These are exact positive duplicate-probe coefficients. They are
strictly weaker than THM-2379's single-complement floors because the
helper-safe probe can consume two of the thirteen shift slots.

## 6. Why the coefficient is not a target current

If

```text
lambda=v/k                         mod p,
```

then the two shifts in (18) have the coordinate pattern `g` from (8).
Equation (10) proves that their `t`-colour does not descend through the
owner-pivot row space. It can detect a relation address whose canonical
target is zero.

There is a second, physical distinction. The weight `omega` in (17)
already retains the unshifted deepest-safe and failed-danger anchors.
The `c`-danger and `v`-safe probes in (18) are intersections against
fixed `c`-safe and `v`-danger anchor roles, not translations replacing
the unique canonical endpoint factors. If the helper-safe factor is
already retained inside `omega`, it too is duplicated.

Thus (24)--(27) prove a positive factor/helper probe, not:

- a THM-2365 target current;
- a pure `a`-target landing;
- a canonical owner current;
- the endpoint-matched cross-word observable required by THM-2380.

The cheapest honest target-side replacement starts with the true dual
dipole (14), shifts **every** occurrence of its two endpoint
coordinates, and then rebuilds a positive mass floor without keeping
fixed anchor copies. That mass problem is not solved here.

## 7. Scope and stopping rule

This theorem records the first failed implication in a tempting bridge:

```text
balanced relation h
  -> pure class [h] in K/L                  true;

same coordinate vector h
  -> target action on K/L                   false;

two-probe mixed Fourier coefficient         true;

that coefficient has canonical target       false.            (28)
```

The source of the error is variance, not a missing scalar factor:
relations are primal and target shifts are covectors. Whenever a
relation atlas is reused as a Fourier-shift atlas, the mandatory check
is annihilation of the entire chosen row space.

No scalar row is excluded. The ledger remains `165`, and LRC(14)
remains open.

## 8. Exact companion

Run

```bash
python 04-computation/lrc14_owner_pivot_primal_dual_repair_thm2384.py
python -O 04-computation/lrc14_owner_pivot_primal_dual_repair_thm2384.py
```

and compare both transcripts byte-for-byte, after LF normalization,
with
`05-knowledge/results/lrc14_owner_pivot_primal_dual_repair_thm2384.out`.
The dependency-free companion checks:

- `248,832` projectively normalized owner-weight tuples, representing
  `2,985,984` raw tuples, for the primal identity (7), universal
  witness (9)--(13), gauge invariance, and both true dual dipoles
  (14)--(15);
- all `12` nonzero helper slopes;
- `196,560` ordinary-unit and `365,040` guard rational-cell phase
  profiles on the full `52`-cell refinement;
- both zero axes, the exact slot identity (21), corner identity (24),
  and energy floor (27); and
- the constants in (26).

Ordinary and optimized runs are byte-identical to the stored transcript.
The LF SHA-256 hashes are recorded in the frontmatter. An independent
hostile audit reconstructed the owner-pivot quotient and dual,
checked every sign in (5)--(15), replayed both zero-axis and
shift-count arguments, and verified all coefficient and energy
constants against the exact companion. QED.
