---
id: THM-2336
title: "Prime-target Gordian owner fan and exact translation-bypass split"
status: >
  PROVED CANDIDATE UNDER INDEPENDENT AUDIT. For any labelled knot packet,
  any partition of its labels, and any nontrivial prime target J,
  THM-2330's lift cost is the root block cost plus the minimum of
  d_G(K_B,J)-u(K_B) over the current partition blocks. The minimizing
  blocks form a tied tropical owner fan on prime targets. If one block is
  J, it is automatically a minimizer. For a two-knot packet (K,J), the
  target obstruction at J is exactly THM-2176's directional translation
  term C_J(K), while the drop from the unknot obstruction to the J
  obstruction is exactly the directional bypass term B_J(K). Thus the
  three landmark values at U,K,J recover both directional splits of the
  symmetric connected-sum defect. No new unknotting-number value, knot
  classification, or positive catalyst is asserted.
source: codex-2026-07-25-prime-target-owner-fan
depends_on:
  - THM-2330-partition-lattice-gordian-lift-spectrum
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
related:
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2248-higher-interaction-defect-complex-and-tropical-trace-spectrum
---

# THM-2336 -- prime-target Gordian owner fan

**PROVED CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2330 replaces the single number `u(#_i K_i)` by a target- and
partition-indexed lifting spectrum. At a general target its endpoint fibre
contains every distribution of the target's prime factors among the current
blocks. At a prime target there is only one token to distribute. This makes
the entire lift cost explicit and turns the surviving binary comparison into
a tied tropical owner relation rather than a tournament.

## 1. Setup

Let `I` be a nonempty finite label set, let

```text
x=(K_i)_(i in I)
```

be a packet of oriented knots, and let `pi` be a partition of `I`. For a
block `B in pi`, put

```text
K_B=#_(i in B)K_i.
```

Retain THM-2330's product Gordian graph, target fibre, lift cost, and
obstruction:

```text
F_pi(J)
 ={(L_B)_(B in pi): #_(B in pi)L_B=J},

Lambda_x(pi;J)
 =min_(L in F_pi(J))
    sum_(B in pi)d_G(K_B,L_B),

Omega_x(pi;J)
 =Lambda_x(pi;J)-d_G(K_I,J).                       (1)
```

Write

```text
u(A)=d_G(A,U).
```

Let `J` be a nontrivial prime knot.

## 2. Prime-target formula

> **Prime-target owner formula.**
>
> ```text
> Lambda_x(pi;J)
>  =min_(B in pi)
>      [d_G(K_B,J)+sum_(C in pi, C!=B)u(K_C)]       (2)
> ```
>
> and equivalently
>
> ```text
> Lambda_x(pi;J)
>  =sum_(C in pi)u(K_C)
>    +min_(B in pi)[d_G(K_B,J)-u(K_B)].             (3)
> ```

### Proof

Schubert prime decomposition says that a connected-sum factorization of the
nontrivial prime `J` has exactly one nontrivial factor, equal to `J`; every
other factor is `U`. Consequently every endpoint in `F_pi(J)` is uniquely of
the form

```text
L_B=J,
L_C=U                         for C!=B              (4)
```

for one block `B in pi`. Its product distance from the prescribed start is

```text
d_G(K_B,J)+sum_(C!=B)d_G(K_C,U)
 =d_G(K_B,J)+sum_(C!=B)u(K_C).                     (5)
```

Minimizing (5) over the possible owner block gives (2), and adding and
subtracting `u(K_B)` gives (3). QED.

Define the prime-target score and owner set

```text
s_B(J)=d_G(K_B,J)-u(K_B),

Own_pi(J)=argmin_(B in pi)s_B(J).                   (6)
```

The triangle inequality through the unknot gives the sharp universal range

```text
-u(J)<=s_B(J)<=u(J).                               (7)
```

Thus the functions `s_B` form a tropical lower envelope on the set of prime
knots. The cells

```text
V_B={J prime:B in Own_pi(J)}                        (8)
```

are a tied Gordian owner fan. Formula (3), not an arbitrary orientation of
the ties, is the invariant statement.

## 3. A prime block owns itself

Suppose one block `D in pi` satisfies

```text
K_D=J.                                              (9)
```

Then

```text
s_D(J)=-u(J).                                       (10)
```

Equation (7) shows that `D in Own_pi(J)` and hence

```text
Lambda_x(pi;J)
 =sum_(C in pi, C!=D)u(K_C).                       (11)
```

Other blocks may tie with `D`; uniqueness is neither needed nor asserted.
Equality in the lower triangle bound is exactly the condition for such a
tie.

For the discrete partition of a two-knot packet

```text
x=(K,J),
```

equation (11) becomes

```text
Lambda_(K,J)(0hat;J)=u(K).                          (12)
```

This identity does not require `K` to be prime.

## 4. The target spectrum recovers translation and bypass

THM-2176 defines the symmetric interaction defect and its directional split:

```text
sigma(K,J)
 =u(K)+u(J)-u(K#J),

C_J(K)
 =u(K)-d_G(K#J,J),

B_J(K)
 =d_G(K#J,J)+u(J)-u(K#J),

sigma(K,J)=C_J(K)+B_J(K).                          (13)
```

At the unknot, THM-2330 gives

```text
Omega_(K,J)(0hat;U)=sigma(K,J).                    (14)
```

At the prime target `J`, equations (1) and (12) give

```text
Omega_(K,J)(0hat;J)
 =u(K)-d_G(K#J,J)
 =C_J(K).                                          (15)
```

Subtracting (15) from (14) yields

```text
Omega_(K,J)(0hat;U)
 -Omega_(K,J)(0hat;J)

 =d_G(K#J,J)+u(J)-u(K#J)
 =B_J(K).                                          (16)
```

Both quantities are nonnegative for their original metric reasons:
connected-sum translation is nonexpansive for `C_J(K)`, and the path

```text
K#J -> J -> U
```

proves nonnegativity of `B_J(K)`.

If both `K` and `J` are prime, the three target landmarks recover the two
directional decompositions:

```text
Omega(U)=sigma(K,J),

Omega(J)=C_J(K),
Omega(K)=C_K(J),

B_J(K)=Omega(U)-Omega(J),
B_K(J)=Omega(U)-Omega(K).                          (17)
```

The symmetric root defect therefore does not lose the directional
mechanisms once the two prime-summand target probes are retained.

## 5. The Brittenham--Hermiller pair

Let

```text
K=T(2,7),                 J=mirror(K).
```

Both are prime. THM-2176's signature calibration proves

```text
C_J(K)=C_K(J)=0,                                   (18)
```

while the cited shortcut gives

```text
sigma(K,J)>=1.                                     (19)
```

Equations (14)--(17) re-express this as

```text
Omega(U)>=1,
Omega(K)=Omega(J)=0,

B_J(K)=B_K(J)=Omega(U)>=1.                         (20)
```

Thus the first connected-sum counterexample is a nontrivial root landmark
with two vanishing prime-summand landmarks. Its shortcut is entirely the
cost of bypassing both separated-summand waypoints, not translation
contraction. No exact value of `u(K#J)` is used or inferred.

## 6. Why the owner fan is not a tournament

The intrinsic comparison at a fixed prime target is

```text
s_B(J) <= s_C(J).                                  (21)
```

It fails the tournament type for four separate reasons.

1. The comparator depends on the external target `J`; changing `J` can
   reverse it.
2. Equality is meaningful: all minimizing blocks are genuine endpoint
   factorizations and must be retained.
3. Formula (3) needs the score magnitude, not only the winning orientation.
4. After partition coarsening the vertices become new composite blocks
   `K_B`; a graph on the original labels has discarded them.

The lossless prime-target shadow is therefore the full score vector

```text
(s_B(J))_(B in pi)
```

modulo addition of a common scalar if only the owner cell is required. The
owner set is its argminimum hyperedge. A tournament is lawful only after an
additional asymmetric target-independent observable resolves every tie and
is proved to preserve the desired lift cost.

For a composite target, its prime multiset must be distributed among all
blocks and the one-token formula (2) becomes a genuine min-plus allocation
problem. The prime owner fan is the exact boundary case, not a
classification of the full target spectrum.

## 7. Scope

- The theorem evaluates an existing spectrum; it computes no previously
  unknown Gordian distance or unknotting number.
- It detects translation catalysis at a **prime summand target**. A general
  catalyst in THM-2191 need not be prime or already occur as a summand.
- Owner ties are retained and may be nontrivial.
- The Brittenham--Hermiller application uses only the already audited
  signature calibration and upper-bound shortcut.
- No binary relation on prime knots is claimed to classify knots, Gordian
  distance, or unknotting number.
