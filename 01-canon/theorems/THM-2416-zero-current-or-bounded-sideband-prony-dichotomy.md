---
id: THM-2416
title: "Zero current or bounded sideband Prony dichotomy"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING. A
  nonzero periodic finite step function with zero mean and L effective
  jump nodes has a nonzero Fourier coefficient at some integer
  1<=X<=L-1. Periodicity supplies the missing zeroth Prony moment, and
  the bound L-1 is sharp. Applied to each eligible physical THM-2410
  packet, either its zero-frequency full-coordinate current survives,
  or a bounded nonzero physical sideband survives; in the latter branch
  the total sideband energy retains the full exponent-20 Gram floor.
  No uniform single-sideband amplitude follows without endpoint
  separation, and the sideband frequency is not yet a relation current.
source: codex-2026-07-26-full-coordinate-sideband-prony
depends_on:
  - THM-2410-full-coordinate-projector-local-gram-and-integrated-phase-boundary
related:
  - THM-2286-endpoint-prony-lift-bank-and-sharp-owner-multiplier-landing
  - THM-2408-endpoint-prony-resultant-clock-separation-and-shared-node-boundary
script: 04-computation/lrc14_zero_current_or_sideband_prony_thm2416.py
output: 05-knowledge/results/lrc14_zero_current_or_sideband_prony_thm2416.out
script_sha256: dd5de18c529fc74744fab7657f6a7467c9cbe976e08d1745f8e0ad58628338de
output_sha256: c21d72bf86fe9d764eb92012866bff710d94d426e29bd5eff6929fc933381129
hash_basis: working-tree bytes (LF)
---

# THM-2416 -- zero current or bounded sideband

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-2410 proves a positive local Gram for every eligible full-coordinate
packet, but correctly stops because its integrated linear coefficient can
cancel.  Finite endpoint structure makes that cancellation much more
rigid than an arbitrary phase loss:

```text
nonzero finite step packet
  -> nonzero mean
     or a nonzero sideband among 1,...,L-1,                  (1)
```

where `L` is the number of distinct effective jump nodes.  The mean-zero
branch retains the entire Gram floor as nonzero-sideband energy.  This
repairs loss of the scalar integral, but not relation neutrality: the new
frequency is an additional physical exponential that a later sidecar
must absorb.

## 1. Periodic jump-Prony lemma

Let `f:T=R/Z -> C` be a finite step function, with endpoint values
ignored.  Combine coincident jumps and discard zero jumps.  Write the
remaining distinct jump nodes and nonzero oriented jumps as

```text
x_1,...,x_L in T,          gamma_1,...,gamma_L in C*.        (2)
```

Use the Fourier convention

```text
fhat(q)=integral_T f(x)exp(-2pi i qx)dx.                     (3)
```

In distributions,

```text
f'=sum_(j=1)^L gamma_j delta_(x_j).
```

Therefore, for every nonzero integer `q`,

```text
2pi i q fhat(q)
 =A_q:=sum_(j=1)^L gamma_j z_j^q,

z_j=exp(-2pi i x_j).                                        (4)
```

Periodicity supplies the extra identity

```text
A_0=sum_j gamma_j=0.                                        (5)
```

Suppose `f` is nonzero and has zero mean.  It cannot be constant, so
`L>=2`.  If

```text
fhat(1)=...=fhat(L-1)=0,                                    (6)
```

then (4)--(6) give `A_q=0` for the `L` consecutive exponents
`q=0,...,L-1`.  The coefficient matrix is the Vandermonde matrix

```text
[z_j^q]_(0<=q<L,1<=j<=L),
```

whose determinant is

```text
product_(i<j)(z_j-z_i)!=0.                                  (7)
```

Hence every `gamma_j=0`, contradicting (2).  We have proved:

> **Periodic sideband lemma.**  If a nonzero periodic finite step
> function has zero mean and `L` effective jump nodes, then `L>=2` and
> some integer `X` satisfies
>
> ```text
> 1<=X<=L-1,             fhat(X)!=0.                         (8)
> ```

This is the zero-residue counterpart of THM-2286's endpoint-Prony lift:
there the residue window has length at most the endpoint count; here
periodicity supplies `A_0=0` and saves one frequency.

Parseval also gives the exact energy identity

```text
sum_(q!=0)|fhat(q)|^2=||f||_2^2                if fhat(0)=0. (9)
```

## 2. Both sharp boundaries

The index `L-1` in (8) is sharp.  Choose any `L>=2` distinct unit-circle
nodes `z_j` and put

```text
P(Z)=product_j(Z-z_j),       gamma_j=1/P'(z_j).                (10)
```

The standard partial-fraction/Vandermonde identities are

```text
sum_j gamma_j z_j^q=0       for 0<=q<=L-2,

sum_j gamma_j z_j^(L-1)=1.                                  (11)
```

In particular `sum gamma_j=0`, so the `gamma_j` are the jumps of a
periodic step function.  Adjusting its additive constant makes its mean
zero without changing any nonzero Fourier coefficient.  Equations
(4) and (11) then give

```text
fhat(1)=...=fhat(L-2)=0,        fhat(L-1)!=0.                 (12)
```

There is no lower bound on the magnitude of the selected sideband in
terms of `||f||_2` and `L` alone.  For `0<epsilon<1`, put

```text
f_epsilon=1_[0,epsilon)-epsilon.                              (13)
```

It has mean zero, two effective jumps, and

```text
||f_epsilon||_2^2=epsilon(1-epsilon),

|fhat_epsilon(1)|=sin(pi epsilon)/pi<epsilon.                 (14)
```

Thus

```text
|fhat_epsilon(1)|^2/||f_epsilon||_2^2
 <epsilon/(1-epsilon) ->0.                                   (15)
```

Endpoint separation, a phase cone, or a conditioned denominator is
necessary for a quantitative single-mode floor.

## 3. A crude but explicit jump invoice

Let `J(g)` denote the number of effective circular jump nodes of a
finite step function `g`.  Directly from one-sided values,

```text
J(gh)<=J(g)+J(h),          J(g(R .))<=R J(g).                 (16)
```

In THM-2410, for an eligible `(m,k)`, put

```text
F_(m,k)(x)
 =Q(Rx)D_m(x)product_(i=1)^9 A_(i,k_i)(x).                   (17)
```

Each integer-speed Boolean factor of speed `v` has at most `2v`
circular endpoints.  Each `A_i` or `D_m` is a sum of thirteen translates,
so a union bound gives

```text
J(A_(i,k_i))<=26w_i,           J(D_m)<=26c.                  (18)
```

Consequently the effective jump count `L_(m,k)=J(F_(m,k))` obeys

```text
L_(m,k)
 <=R J(Q)+26(sum_(i=1)^9 w_i+c)
 =:L_crude.                                                   (19)
```

Coincident endpoints and product cancellations only improve this bound.

## 4. Application to every eligible THM-2410 packet

THM-2410 proves

```text
B(m,k)=Fhat_(m,k)(0),                                        (20)

||F_(m,k)||_2^2
 >=mu(Q)(sin(pi/13)/13)^20
 >mu(Q)(2/169)^20                                           (21)
```

for every

```text
m!=0,           every k_i!=0,           k.wbar=0.            (22)
```

Therefore exactly one of the following holds.

### Mean branch

```text
B(m,k)!=0.                                                   (23)
```

Then THM-2410's conditional Abel extraction supplies an all-coordinate
unit residue modulo thirteen and an exact deep multiplier coprime to
`91`.

### Sideband branch

```text
B(m,k)=0.                                                    (24)
```

The Gram floor (21) makes `F_(m,k)` nonzero.  The periodic sideband lemma
gives

```text
1<=X_side<=L_(m,k)-1<=L_crude-1,

Fhat_(m,k)(X_side)!=0,                                       (25)
```

and Parseval preserves the full quantitative statement

```text
sum_(X!=0)|Fhat_(m,k)(X)|^2
 =||F_(m,k)||_2^2
 >=mu(Q)(sin(pi/13)/13)^20.                                 (26)
```

Equation (25) is a physical total-frequency defect.  Multiplying (17)
by `exp(-2pi i X_side x)` changes the total frequency balance; it does
not produce THM-2410's zero-frequency relation current, retain its exact
Abel address, or imply an all-`91`-unit relation.  A common endpoint or
pure-exponential sidecar must absorb `X_side`.

## 5. Exact scope

This theorem proves:

- the abstract sharp `L-1` sideband lemma for complex finite step
  functions;
- the exact Parseval energy retention in the zero-mean branch;
- the explicit THM-2410 jump invoice (19); and
- the mean-or-bounded-sideband dichotomy for every eligible physical
  full-coordinate packet.

It does not prove:

- that a THM-2410 mean is nonzero;
- any uniform magnitude for the bounded sideband;
- that the sideband is relation-neutral or shares a preselected owner,
  triangle, or terminal endpoint;
- a mod-seven word-phase repair; or
- a scalar-row exclusion or LRC(14).

The scalar ledger remains `165`.

## 6. Exact companion

The dependency-free companion:

- verifies the Vandermonde determinant and the sharp weights
  `1/P'(z_j)` over exact rational node controls for `2<=L<=9`;
- checks the `L-1` moment boundary and a nonzero hostile control;
- verifies the exact two-jump norm and the rational upper bound
  `epsilon/(1-epsilon)` from (15);
- checks the product/composition jump invoice on finite cyclic step
  tables; and
- records the exact rational THM-2410 exponent-`20` Gram floor.

Run:

```bash
python3 04-computation/lrc14_zero_current_or_sideband_prony_thm2416.py
python3 -O 04-computation/lrc14_zero_current_or_sideband_prony_thm2416.py
```

Both commands must reproduce

```text
05-knowledge/results/lrc14_zero_current_or_sideband_prony_thm2416.out
```

with the LF hashes in the frontmatter.
