---
id: THM-2164
title: "Relative Selberg-packet rank harvesting for LRC(14)"
status: >
  PROVED from THM-2085 and THM-2144. On the codimension-one subtorus
  defined by any primitive relation supported on at least three coordinates,
  the degree-43 signed Selberg tensor has a positive Haar packet. Therefore
  a zero-safe row has a second independent relation of height at most 43.
  For a signed-unit relation the residue-sensitive degree-34 packet suffices.
  Consequently every zero-safe distinct thirteen-speed row has either
  rank(W_43)>=2 or a two-term lock of coefficient height at most 29, and
  unconditionally rank(W_105)>=2. These are relation-rank advances, not a
  proof of LRC(14).
source: codex-2026-07-24-relation-carry-spectrum
depends_on:
  - THM-2085
  - THM-2144
related:
  - THM-2052
  - THM-2054
script: 04-computation/lrc14_relative_packet_rank_harvesting_referee_codex_20260724.py
output: 05-knowledge/results/lrc14_relative_packet_rank_harvesting_referee_codex_20260724.out
script_sha256: edd4a3981053011e6fca3eda57970a323b1ac87b81ed2dfe46de5c32218952dc
output_sha256: e86b9c733a5e52ed9117dda36c1b6274a2465b33ca14427cd7d839b82668eb22
hash_basis: working-tree bytes (LF)
---

# THM-2164 -- relative Selberg-packet rank harvesting

Let

```text
J=[1/14,13/14],             chi=1_J,
S(v)={t:(v_1t,...,v_13t) in J^13}.                    (1)
```

For a positive integer row `v`, write

```text
Lambda(v)={m in Z^13:m.v=0},
W_H(v)=span_Q{m in Lambda(v):||m||_infinity<=H}.       (2)
```

The result turns one known relation into a relative Fourier ambient space.
It is the codimension-one packet specialization of THM-2054's
whole-product philosophy, but the positive packet estimates below are new.

## 1. The whole-subtorus packet lemma

Let `0!=w in Lambda(v)` be primitive and define

```text
T_w={x in T^13:w.x=0 mod 1}.                           (3)
```

The annihilator of `T_w` is exactly `Z w`. Let `P_H` be the equal-degree
signed Selberg tensor from THM-2085:

```text
P_H(x)
 =product_i U(x_i)
  -sum_i D(x_i) product_(j!=i)U(x_j),
D=U-L,                                                 (4)
```

where `L<=chi<=U` has degree `H`. Thus

```text
P_H<=product_i chi,
supp(P_H hat) subset [-H,H]^13.                        (5)
```

**Packet lemma.** If

```text
integral_(T_w) P_H(x) dx>0                             (6)
```

and `mu(S(v))=0`, then there is

```text
r in Lambda(v)\Qw,          ||r||_infinity<=H.         (7)
```

### Proof

Assume no `r` in (7) exists. A Fourier character in the support box
survives integration on the line `t|->vt` exactly when it belongs to
`Lambda(v)`. By the assumption, every such character lies in `Qw`; because
it is integral and `w` is primitive, it lies in `Zw`. Conversely every
multiple of `w` in the support box survives both averages. Hence finite
Fourier expansion gives

```text
integral_T P_H(vt)dt=integral_(T_w)P_H(x)dx>0.         (8)
```

But (5) and the zero-safe hypothesis give

```text
integral_T P_H(vt)dt
 <=integral_T product_i chi(v_i t)dt
 =mu(S(v))=0,                                          (9)
```

a contradiction. QED.

This comparison retains all old multiples of `w`; it does not discard them
termwise. That is why it avoids the invalid heuristic that a constant
Fourier coefficient should dominate every relation channel separately.

## 2. A coefficient-free packet majorant

Put

```text
alpha=6/7,        epsilon=1/(H+1),
u=alpha+epsilon,  d=2epsilon.                         (10)
```

The order relations for the Selberg pair imply

```text
U-chi>=0,       integral(U-chi)=epsilon,
D>=0,           integral D=d.                         (11)
```

Therefore, at every nonzero integer frequency `n`,

```text
|U hat(n)|<=|chi hat(n)|+epsilon,
|D hat(n)|<=d.                                         (12)
```

No explicit formula or sign for a Vaaler coefficient is used.

The zero coefficient of (4) is

```text
B_H=u^12(alpha-25epsilon).                             (13)
```

Suppose `w` has support `s>=3`. For `l>=1`, every nonzero coordinate of
`lw` is an integer of magnitude at least `l`. Since

```text
|chi hat(n)|
 =|sin(pi n/7)|/(pi|n|)
 <5/(16|n|),                                           (14)
```

put

```text
A_l=5/(16l)+epsilon.                                   (15)
```

Taking absolute values only after the complete coefficient of (4) has been
formed gives

```text
|P_H hat(lw)|
 <=A_l^s u^(13-s)
   +s d A_l^(s-1)u^(13-s)
   +(13-s)d A_l^s u^(12-s).                           (16)
```

Frequencies outside the tensor support simply vanish, so summing all
`1<=l<=H` is a valid overcount. Equations (13)--(16) give

```text
integral_(T_w)P_H
 >=B_H-2 sum_(l=1)^H [
      A_l^s u^(13-s)
     +s d A_l^(s-1)u^(13-s)
     +(13-s)d A_l^s u^(12-s)].                        (17)
```

At `H=43`, the exact rational minimum of (17) over `3<=s<=13` occurs at
`s=3` and is greater than

```text
1/634>0.                                               (18)
```

The same support-three ledger is negative at `H=42`; this is only the
boundary of the displayed majorant, not an optimality claim.

For completeness, (14) follows from the integer residue classes modulo
seven. With the elementary bounds `157/50<pi<22/7`,

```text
|chi hat(n)| <=
  0             if n=0 mod 7,
  1/(7|n|)      if n=+/-1 mod 7,
  1/(4|n|)      if n=+/-2 mod 7,
  5/(16|n|)     if n=+/-3 mod 7.                      (19)
```

The first line is exact and the second uses `sin x<x`. For the third,
`sin(2pi/7)=cos(3pi/14)`; the quartic alternating upper bound for cosine,
monotonicity of that polynomial below one, and

```text
3pi/14>471/700,
1-(471/700)^2/2+(471/700)^4/24<157/200<pi/4           (20)
```

give the claim. For the fourth,

```text
sin(3pi/7)=cos(pi/14)
 <1-(3/14)^2/2+(3/14)^4/24
 <49/50<5pi/16.                                       (21)
```

Thus the rational wrapper in the exact referee proves every constant used
in (14)--(19).

Combining (17)--(18) with the packet lemma proves:

```text
mu(S(v))=0 and primitive w in Lambda(v), |supp w|>=3
  implies some r in Lambda(v)\Qw with ||r||_infinity<=43.  (22)
```

## 3. The signed-unit packet

If every nonzero coordinate of `w` is `+1` or `-1`, the sharper
residue-sensitive table (19) can be retained at each multiple `lw`. Define

```text
C_l=
  0          for l=0 mod 7,
  1/(7l)     for l=+/-1 mod 7,
  1/(4l)     for l=+/-2 mod 7,
  5/(16l)    for l=+/-3 mod 7,
A_l=C_l+epsilon.                                       (23)
```

Use these `A_l` in (17). At `H=34`, the exact rational minimum over
`3<=s<=13` again occurs at support three and is greater than

```text
1/212>0.                                               (24)
```

The corresponding support-three certificate is negative at `H=33`.
Therefore

```text
mu(S(v))=0,
w primitive with w_i in {-1,0,1}, |supp w|>=3
  implies some r in Lambda(v)\Qw with ||r||_infinity<=34.  (25)
```

## 4. Two rank consequences

Assume now that `v_1,...,v_13` are pairwise distinct positive integers and
`mu(S(v))=0`.

THM-2144 first supplies a nonzero height-29 relation. Replace it by its
primitive generator `w`. A one-coordinate relation is impossible. If
`|supp w|>=3`, (22) supplies an independent height-43 relation. If
`|supp w|=2`, positivity gives a primitive two-term lock

```text
a v_i=b v_j,              1<=a,b<=29.                 (26)
```

Hence

```text
dim_Q W_43(v)>=2
or a bounded two-term lock (26) occurs.                (27)
```

This fork isolates the only obstruction to immediate low-height rank two.

The obstruction disappears at THM-2144's anisotropic height. Suppose for
contradiction that `dim_Q W_105(v)=1`, and let `w` generate the saturated
rank-one lattice `W_105(v) intersect Z^13`. The thirteen cap profiles
`(1,105,...,105)` used in THM-2144 show that every coordinate of `w` lies in
`{-1,0,1}`. Both signs occur. Pairwise distinctness rules out support two,
since that would say `v_i=v_j`; thus `|supp w|>=3`. Equation (25) now gives
an independent relation of height at most `34<=105`, a contradiction.
Therefore

```text
dim_Q W_105(v)>=2.                                     (28)
```

Equation (28) removes THM-2144's signed-subset-sum alternative for the
distinct rows relevant to LRC(14). It does not make the speed search finite,
force sparse rank two, or prove the lonely-runner inequality.

## 5. Exact referee and boundaries

The companion uses only `Fraction` arithmetic. It verifies (20)--(21), every
support value `3,...,13` in the two packet sums, the positive bounds
(18) and (24), and the adjacent negative ledgers at degrees `42` and `33`.
It contains explicit raising checks, so optimized Python does not remove the
validity gate. Normal and optimized executions are transcript-identical.
