---
id: THM-2644
title: "Odd-torsor purity/return gate and nonlinear fixed-branch decoder"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For a
  nonnegative weight mu on a finite group, put M=sum mu, E=sum mu^2,
  delta=M^2-E, and R=(mu*mu)(e).  The squared mass on involutions is at least
  R-delta.  Hence on an odd group R>delta forces positive identity mass.  If
  additionally E=M^2 and M>0, mu is a singleton and the return condition
  identifies the identity exactly.  On C_p the k-fold version works precisely
  when gcd(k,p)=1; k divisible by p is the sharp hostile.  This is a nonlinear
  fixed-branch decoder, but it requires two compositions of one nonnegative
  transition through a common physical middle fibre.  Reverse/Gram data alone,
  unrelated endpoint marginals, and complex signed currents do not supply it.
source: wild-holotopy-mining-2026-07-28-purity-return-gate
depends_on: []
related:
  - THM-2315-marked-target-gain-corolla-and-pairwise-composition-boundary
  - THM-2634-endpoint-pair-two-carry-cospan-and-single-carry-no-go
  - THM-2637-derangement-character-fixed-branch-holotopy-principle
  - THM-2642-cyclic-difference-relation-saturation-and-thick-holotopy-no-go
script: 04-computation/lrc14_odd_torsor_purity_return_gate_thm2644.py
output: 05-knowledge/results/lrc14_odd_torsor_purity_return_gate_thm2644.out
script_sha256: 9444f551aeae69688b2f84bbddc8ecec100a4dbe031572501256d978c84092cf
output_sha256: ac059a1cf020778abbb9f6e411d39b9775089b077ff3b8edb938de71164e74af
hash_basis: LF-normalized bytes
---

# THM-2644 -- two quadratic rungs recover an odd-torsor fixed branch

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Private coefficient rows are not the only algebraic way to locate a branch.
There is a sharp nonlinear alternative: first prove that a nonnegative
transition is pure, then test whether its unique branch returns under an
oriented two-step composition.  The two rungs are different.  The first is a
reverse/Gram loop and forgets the branch; the second keeps orientation and
detects the surviving torsor element.

## 1. The involution-energy inequality

Let `G` be a finite group with identity `e`, and let

```text
mu:G -> R_>=0.                                             (1)
```

Use unnormalized convolution and put

```text
M=sum_g mu(g),                 E=sum_g mu(g)^2,
delta=M^2-E,                   R=(mu*mu)(e)
                                  =sum_g mu(g)mu(g^-1),     (2)

Inv(G)={h in G:h=h^-1},        J=sum_(h in Inv(G))mu(h)^2. (3)
```

Then

```text
J >= R-delta.                                               (4)
```

Indeed nonnegativity gives the exact off-diagonal expansion

```text
delta=sum_(g!=h) mu(g)mu(h).                                (5)
```

Splitting the return sum at the involutions gives

```text
R=J+sum_(g notin Inv(G))mu(g)mu(g^-1).                      (6)
```

Every term in the last sum is one of the ordered off-diagonal terms in (5),
so that sum is at most `delta`.  This proves (4).

The equality boundary is exact.  Equality in (4) holds precisely when the
positive support of `mu` has at most one point, or is one inverse pair
`{a,a^-1}` with `a!=a^-1`.  In particular the strict sign in

```text
R>delta                                                     (7)
```

is sharp: an inverse pair has `J=0` and `R=delta`.

## 2. Odd groups and the robust fixed-branch gate

If `G` has odd order, its only involution is `e`.  Equation (4) becomes

```text
mu(e)^2 >= R-(M^2-E).                                      (8)
```

Consequently (7) forces a positive identity branch, quantitatively

```text
mu(e) >= sqrt(R-(M^2-E)).                                  (9)
```

This conclusion does not require purity.  For example, on `C_13` the weight
`3 delta_0+delta_1` has `(M,E,delta,R)=(4,10,6,9)`, so the robust gate fires
even though two branches remain.

The exact private-branch specialization uses

```text
E=M^2,                   M>0.                              (10)
```

By (5), (10) is equivalent to `mu` having one positive support point, say

```text
mu=M delta_a.                                              (11)
```

Then `R>0` says `a^2=e`; on an odd group this forces `a=e`.  Thus

```text
E=M^2 and R>0
   iff mu=M delta_e                              (M>0).     (12)
```

This is an exact nonlinear decoder of the identity row even when no linear
private-row functional has been supplied.

## 3. The complete cyclic k-rung version

On the additive torsor `C_p`, with `p` prime, define

```text
R_k=(mu^(*k))(0).                                         (13)
```

After the purity rung (10), equation (11) gives

```text
R_k = M^k if k a=0 mod p, and 0 otherwise.                (14)
```

Therefore, whenever `gcd(k,p)=1`, the condition `R_k>0` forces `a=0`.
If `p|k`, every singleton translation has `R_k=M^k`, so the conclusion fails
for every nonzero `a`.  For odd `p`, `k=2` is the cheapest nontrivial return
and is exactly (12).

## 4. Why the two compositions are not interchangeable

Write the weighted transition in the real group algebra as

```text
A=sum_g mu(g)[g].                                         (15)
```

The coefficient of `e` in the reverse composition `A^*A` is `E`; the
coefficient of `e` in the same-orientation composition `A^2` is `R`.
On `C_p`, for the unnormalized Fourier transform,

```text
E=(1/p)sum_k |muhat(k)|^2,
R=(1/p)sum_k muhat(k)^2.                                  (16)
```

Thus the second rung is phase-sensitive.  Every singleton `delta_a` has the
same full Gram spectrum `|muhat(k)|^2=1`, while only `a=0` survives the
two-step return when `p` is odd.  Repeating a reverse/Gram composition cannot
replace the oriented return.

In path language, `A^*A` is a backtracking two-gon and `A^2` is an oriented
two-gon.  Purity collapses a thick relation to one transport branch; the
oriented two-gon then tests whether that branch has trivial holonomy.  This is
the nonlinear thin/thick bridge between THM-2637 and THM-2642.

## 5. Sharp hostile controls

Every hypothesis and inequality direction is load-bearing.

1. **Purity alone.**  On `C_13`, `mu=delta_5` has `E=M^2` but `R=0`.
2. **Return alone.**  `mu=delta_5+delta_8` has no identity mass and satisfies
   `R=delta=2`; hence `R>=delta` is insufficient.
3. **Odd order.**  On `C_2`, `mu=delta_1` has `E=M^2` and `R>0`, but its
   branch is not the identity.
4. **Nonnegativity.**  On `C_13`, the signed weight
   `mu=delta_1-delta_2` has `M=0`, `E=2`, `delta=-2`, and `R=0>delta`, while
   `mu(0)=0`.
5. **Coprime return length.**  For `k=13`, every `delta_a` has `R_k=1`.
6. **A private linear row is genuinely different.**  Sampling `mu(0)` would
   settle the issue immediately; (12) reconstructs the same conclusion only
   from the two quadratic compositions and positivity.

## 6. Physical typing boundary

The algebra does **not** authorize an LRC application by itself.  To use
(8) or (12), one must construct on one physical object

```text
one nonnegative transition A,
one common source and target torsor gauge,
a lawful middle fibre on which both A^*A and A^2 compose,
and the augmentation M of that same A.                    (17)
```

THM-2634 supplies a two-carry complex endpoint cospan, not a single
nonnegative `A` with the common composability in (17).  THM-2642 supplies an
abstract relation-convolution theorem and explicitly does not promote the
current static coefficient rows to physical transitions.  Separate endpoint
marginals, separate positive rows, or a complex signed current can reproduce
one rung while violating the other.

The cheapest decisive LRC test is therefore to find one same-base transition
packet, compute both the reverse and oriented length-two identity
coefficients before marginalization, and compare `R` with `M^2-E`.  A strict
gap gives a genuine fixed branch by (8); exact purity upgrades it to a private
identity branch by (12).  No such packet is constructed here.  No scalar row
is removed and the LRC ledger remains `165`.

## 7. Exact reproduction

Run

```bash
python 04-computation/lrc14_odd_torsor_purity_return_gate_thm2644.py
python -O 04-computation/lrc14_odd_torsor_purity_return_gate_thm2644.py
```

The dependency-free companion uses explicit optimization-safe guards.  It
checks every `{0,1,2}`-valued weight on `C_3,C_5,C_7`, all `8,192` Boolean
weights on `C_13`, every positive integer weight of support at most three and
height at most three on `C_13`, the equality classification, both group-
algebra compositions, all six sharp hostiles, and every cyclic singleton
return for `2<=k<=26`.  Normal and optimized executions must byte-match the
stored transcript and end in `PASS`.

Two independent hostile audits rederived the involution inequality and its
equality cases, the odd-group and cyclic return gates, every sharp hostile,
the group-algebra/Fourier signs, the physical typing boundary, and both exact
transcripts and hashes.

QED.
