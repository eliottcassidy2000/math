---
id: THM-2420
title: "Affine-shell cross-reference composition and complete zero-reference hostile"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. On one
  nonzero affine shell, a charged coefficient
  packet and a residue-zero reference packet compose by oriented
  correlation into a homogeneous charged relation packet. The total,
  support, Fourier product, sharp max, and sharp energy bounds are
  exact. A three-coordinate Boolean projector hostile has a nonzero
  all-unit sideband at X=26 while every coordinate-residue-zero
  projector sector and the self-Gram are constant. Thus a same-shell
  address, complete zero-character bank, and positive local Gram do
  not construct the needed reference amplitude. The theorem is
  abstract: no canonical nine-coordinate reference, original linear
  LRC current, row exclusion, or LRC(14) conclusion is proved.
source: codex-2026-07-26-affine-shell-reference-audit
depends_on:
  - THM-2419-valuation-normalized-homogenization-of-affine-sideband-shells
related:
  - THM-2325-complete-character-address-tomography-and-invisible-full-support-boundary
  - THM-2333-canonical-root-target-abel-fibre-decomposition-and-zero-target-hostile
  - THM-2334-marked-current-residue-inversion-and-inverse-character-boundary
  - THM-2383-polarized-complete-subcube-gram-tomography
  - THM-2410-full-coordinate-projector-local-gram-and-integrated-phase-boundary
  - THM-2416-zero-current-or-bounded-sideband-prony-dichotomy
script: 04-computation/lrc14_same_shell_cross_reference_thm2420.py
output: 05-knowledge/results/lrc14_same_shell_cross_reference_thm2420.out
script_sha256: 18cf938ac7a3ff9036b1b9f412a2687391ee856387fd34c22dee917768df2829
output_sha256: 2ca80104b72ef36b30150f6c0c945f530a7a715bdbf68b999fde3c04a5e1a4a2
hash_basis: working-tree bytes (LF)
---

# THM-2420 -- a same-shell reference is sufficient, but no reference is free

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2416 produces a nonzero affine sideband and THM-2419 identifies its
canonical finite torsor.  Homogenizing it inside the physical lattice
requires a reference with nonzero aggregate amplitude, not merely an
address or a cancelling collection of coefficients, in the residue-zero
sector of the same shell.  The exact operation is a polarized
cross-correlation:

```text
charged affine packet at X
  + residue-zero reference packet at the same X
  -> charged homogeneous cross-current,

but

all residue-zero projector sectors may have zero amplitude at X.       (1)
```

The second line is realized by a three-coordinate Boolean projector.

## 1. Finite same-shell cross-reference composition

Let `w in Z^n` be primitive, let `X!=0`, put `M=|X|`, and let the prime
`p` divide `M`.  Define

```text
K_M=ker(w:(Z/M)^n -> Z/M),
K_p=ker(w:(Z/p)^n -> Z/p).                                  (2)
```

Reduction gives a surjection

```text
pi:K_M -> K_p
```

whose every fibre has cardinality

```text
N_0=(M/p)^(n-1).                                             (3)
```

Indeed a unimodular change of coordinates sends primitive `w` to
`(1,0,...,0)`, where both assertions are immediate.

Fix `q in K_p`.  Let `C,R:K_M -> C` satisfy

```text
supp(C) subset pi^(-1)(q),
supp(R) subset pi^(-1)(0),                                  (4)
```

and put

```text
S_C=sum_u C(u),              S_R=sum_v R(v).
```

Define the oriented correlation

```text
H(t)=sum_(v in K_M) C(t+v) conjugate(R(v)).                  (5)
```

Then

```text
supp(H) subset pi^(-1)(q),                                  (6)

sum_t H(t)=S_C conjugate(S_R).                              (7)
```

For the unnormalized transform

```text
fhat(chi)=sum_x f(x) conjugate(chi(x)),
```

one also has

```text
Hhat(chi)=Chat(chi) conjugate(Rhat(chi)).                    (8)
```

The support statement follows from
`pi(t)=pi(t+v)-pi(v)=q`; (7) follows from the bijection
`(t,v)<->(u=t+v,v)`; and (8) is the same substitution with a
character inserted.

Consequently `S_C S_R!=0` forces some relation fibre
`t mod p=q` with `H(t)!=0`; it is charged when `q!=0`.  More precisely,

```text
max_(pi(t)=q) |H(t)|
 >= |S_C S_R|/N_0,                                          (9)

sum_(pi(t)=q) |H(t)|^2
 >= |S_C S_R|^2/N_0.                                       (10)
```

Both constants are sharp.  Take `R=S_R delta_(v_0)` on the zero fibre
and take `C=S_C/N_0` uniformly on `pi^(-1)(q)`.  After replacing `R`
by `R/S_R`, the total cross-current in (7) is exactly `S_C`: aggregate
charged amplitude, not merely nonvanishing, is preserved.

## 2. Exact-address and Abel formulation

Let

```text
Lambda={r in Z^n:r.w=0},
Lambda_X={a in Z^n:a.w=X}.                                  (11)
```

At a positive Abel radius `0<rho<1`, let
`c_rho,d_rho in ell^1(Lambda_X)` and define

```text
h_rho(r)
 =sum_(b in Lambda_X)
    c_rho(b+r) conjugate(d_rho(b)),          r in Lambda.    (12)
```

Absolute convergence gives

```text
||h_rho||_1 <= ||c_rho||_1 ||d_rho||_1,                     (13)

sum_r h_rho(r)
 =(sum_a c_rho(a)) conjugate(sum_b d_rho(b)).                (14)
```

Push `c_rho,d_rho` to `Lambda_X/(M Lambda)=K_M` as
`C_rho,R_rho`.  Then the mod-`M` pushforward of (12) is exactly

```text
sum_(r=t mod M Lambda) h_rho(r)
 =sum_v C_rho(t+v) conjugate(R_rho(v)).                      (15)
```

Thus coordinatewise boundary limits `C_rho(u)->C(u)` and
`R_rho(v)->R(v)`, taken along the same Abel radii and in one common
gauge, give the finite correlation in Section 1, with the limiting
support conditions (4).  If in addition `S_C S_R!=0`, a nonzero
limiting `H(t)` is a genuine nonzero Abel relation-residue fibre.

Three qualifications are load-bearing:

1. a residue-zero **address**, or even one nonzero coefficient inside a
   cancelling reference packet, is insufficient for the aggregate
   criterion; one needs `S_R!=0` or an independently nonzero
   cross-correlation coefficient;
2. a finite residue fibre need not select one fixed exact relation
   address at `rho=1` without an undamped `ell^1` or no-escape sidecar;
3. (12) is a polarized quadratic cross-current.  It is not
   automatically one of the original linear LRC currents unless a
   physical realization theorem supplies the labelled cross-phase.

Separate intensity banks do not determine (12); the missing datum is
the same polarized Gram coordinate isolated abstractly by THM-2383.

## 3. A complete residue-zero bank can miss the sideband

Put `p=13`, `zeta=exp(2 pi i/13)`, and on the circle let

```text
J=1_[0,1/13),               I=1-J.
```

For an integer speed `v` and residues `k,m in F_13`, define

```text
A_(v,k)(x)
 =13^(-1) sum_(s in F_13) I(vx-s/13) zeta^(ks),

D_(v,m)(x)
 =13^(-1) sum_(s in F_13) J(vx-s/13) zeta^(ms).             (16)
```

Away from endpoints, let

```text
sigma_v(x)=floor(13 {vx}).
```

The shifted `J` cells partition the circle, so

```text
A_(v,0)=12/13,
A_(v,k)=-13^(-1) zeta^(k sigma_v)       for k!=0,
D_(v,m)=13^(-1) zeta^(m sigma_v).                           (17)
```

Take the three present speeds, safe characters, and deepest character

```text
v=(1,2,26),                 k=(1,6,2),       m=12.           (18)
```

Then

```text
k.v=0 mod13,
q=k+m e_3=(1,6,1) in (F_13^*)^3,
m!=-k_3.                                                        (19)
```

The eligible all-unit endpoint/deep packet is

```text
F(x)=D_(26,12) A_(1,1) A_(2,6) A_(26,2).                    (20)
```

It is `1/13`-periodic and

```text
|F(x)|^2=13^(-8).                                           (21)
```

In the quotient coordinate `t={13x}`, the phase word on the twenty-six
equal cells is

```text
0,1,...,12, 6,7,...,12,0,1,...,5.                           (22)
```

Every residue occurs twice, and therefore

```text
Fhat(0)=0.                                                   (23)
```

At the physical sideband `X=26`, direct cell integration gives

```text
Fhat(26)
 =-(1-zeta^(-1))(1+zeta^6)/(4 pi i 13^3) !=0.               (24)
```

The first factor is nonzero because `zeta!=1`; the second is nonzero
because an odd-order root of unity cannot equal `-1`.

Now enumerate the entire coordinate-residue-zero projector sector.
Its characters must have

```text
(k_1,k_2,k_3,m)=(0,0,-t,t),               t in F_13.
```

Every possible reference is consequently

```text
R_t=D_(26,t) A_(1,0) A_(2,0) A_(26,-t),                     (25)
```

and (17) gives the constants

```text
R_0=12^3/13^4,
R_t=-12^2/13^4                for t!=0.                      (26)
```

Thus

```text
Rhat_t(X)=0                    for every X!=0 and every t,   (27)
```

including `X=26`.  The self-Gram `|F|^2` is constant by (21), so it
also has no nonzero sideband.  Hence all of the following can hold
simultaneously:

- a nonzero all-unit affine sideband;
- positive local Gram energy;
- the complete coordinate-residue-zero projector bank; and
- no same-shell residue-zero reference amplitude.

This is a genuine hostile inside the exact Boolean character-projector
operation.  It is not a literal LRC row: it has three rather than nine
present coordinates, no wide guard, no terminal word, and `1/13` model
cells rather than the final interval geometry.  It refutes an inference
from projector symmetry, local Gram energy, and Prony sideband survival;
it does not refute a construction using the omitted LRC geometry.

Within this formal architecture, three present coordinates are minimal:
with only one `13`-unit speed and one `13`-divisible deepest speed, a
present-character annihilator has zero coefficient on the unit speed,
so the affine residue cannot be all-unit.

## 4. Exact frontier

The theorem makes the post-THM-2419 target precise:

```text
Omega(C,R)
 =sum_(X in 13Z, X!=0)
    |Fhat_q(X)|^2 |Rhat_0(X)|^2 >0,                          (28)
```

with a common gauge, Abel regulator, labelled phase, and physical
realization.  For an all-`91`-unit conclusion, the same packet must
retain the septimal residue/word phase.  THM-2418 shows that the latter
can independently collapse to a flat carry class.

No canonical same-shell reference is constructed here, no original
linear current is recovered, no row is excluded, the scalar ledger
remains `165`, and LRC(14) remains open.

## 5. Exact companion

The dependency-free companion:

- verifies the reduction `K_M -> K_13`, fibre sizes, support, total,
  sharp max, and sharp energy bounds on four primitive finite shells;
- verifies the exact group-ring Fourier product on sampled characters;
- checks an exact Gaussian coefficient-`ell^1` majorant and the
  mod-`M` address pushforward identities on a finite shell;
- verifies (18)--(24), including the twenty-six-cell phase census and
  the nonzero cyclotomic numerator; and
- enumerates all thirteen zero-residue sectors in (25)--(27) and the
  constant self-Gram.

Run

```bash
python3 04-computation/lrc14_same_shell_cross_reference_thm2420.py
python3 -O 04-computation/lrc14_same_shell_cross_reference_thm2420.py
```

and compare both transcripts, after LF normalization, with

```text
05-knowledge/results/lrc14_same_shell_cross_reference_thm2420.out
```

## 6. Independent hostile audits

Two independent audits accepted the theorem for promotion.  The proof
audit rederived the shell-fibre size, oriented correlation, Fourier
product, sharp inequalities, Abel pushforward, every zero-reference
sector, the nonzero sideband, and three-coordinate minimality.  The
exact audit independently replayed normal and optimized modes against
the stored transcript, checked all `676` shell fibres and `256`
ambient-character samples, verified strict hostile variants, and
confirmed the repaired twenty-six-cell self-Gram derivation.  The
audited LF hashes are the frontmatter hashes above.
