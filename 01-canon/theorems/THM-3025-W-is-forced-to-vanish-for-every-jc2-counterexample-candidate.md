---
id: THM-3025
title: "JC(2): the subleading correction W vanishes identically for every counterexample candidate"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT. Answers the
  first open question named in THM-3016 section 5 ("whether W = 0 is forced in
  general"). Setup as in THM-3016: a planar Jacobian pair with leading forms
  P_n = c H^a, Q_m = c' H^b, H homogeneous of degree g, n = ga, m = gb,
  gcd(a,b) = 1, and W = P_{n-1} - kappa H^{a-b} Q_{m-1} of degree n-1.
  (W1) THM-3016 (R) gives J * Jac(W,H) = 0, and J = Jac(P,Q) is a NONZERO
  constant, so Jac(W,H) = 0 outright.
  (W2) For binary forms Jac(W,H) = 0 with W != 0 is equivalent to
  W^g = c H^{n-1}, i.e. W and H are powers of one common form.
  (W3) ARITHMETIC OBSTRUCTION: gcd(g, n-1) = gcd(g, ga-1) = 1 for ALL g,a.
  Writing H = prod l_i^{e_i}, (W2) forces g | e_i(n-1), hence g | e_i, hence
  (since sum e_i = g, e_i >= 1) exactly one factor with e_i = g. So W != 0
  forces H = l^g, a pure power of a SINGLE linear form -- that is K = 1, the
  classically resolved one-place-at-infinity configuration (Abhyankar-Moh).
  (W4) THEREFORE: K >= 2  =>  W = 0, unconditionally. Since a JC(2)
  counterexample requires K >= 2 (HYP-9070 gate), W vanishes identically for
  every counterexample candidate.
  (W5) CONSEQUENCE: THM-3016 section 4b(B)'s cross-term tower -- the degree
  n+m-4 relation with all three terms carrying explicit powers H^{a-1},
  H^{a-b-1}, H^{b-1}, and the Euclidean-reducing iteration -- was derived
  UNDER the hypothesis W = 0. That hypothesis is now discharged: the tower
  holds unconditionally on the counterexample locus.
  The dichotomy is SHARP, not vacuous: at K = 1 the solution space is
  genuinely nonzero, W = c l^{ga-1}, verified for (g,a) = (2,2) and (3,2).
source: death-star-2026-07-31-coinC2
depends_on:
  - THM-3016
related:
  - HYP-9070
  - THM-3003
external:
  - "Abhyankar, Moh: embeddings of the line in the plane (one place at infinity)."
  - "Jung; van der Kulk: the structure of the plane automorphism group."
script: 04-computation/jc2_W_forced_thm3025.py
output: 05-knowledge/results/jc2_W_forced_thm3025.out
---

# THM-3025 -- `W = 0` is forced

## 1. Setup

As in THM-3016. `(P,Q)` is a planar Jacobian pair, `J = Jac(P,Q)` a nonzero
constant, with forced leading forms

```text
P_n = c H^a,   Q_m = c' H^b,   deg H = g,   n = ga,   m = gb,   gcd(a,b) = 1,
```

and the subleading correction

```text
W = P_{n-1} - kappa H^{a-b} Q_{m-1},      deg W = n-1.
```

THM-3016's rigidity identity (R) is `J * Jac(W,H) = 0`. Since `J` is a
**nonzero** constant,

```text
Jac(W, H) = 0.                                                        (W1)
```

## 2. Binary forms: vanishing Jacobian means common power (W2)

For homogeneous `W, H` in two variables, `Jac(W,H) = 0` holds **iff** `W` and
`H` are powers of a common form, equivalently

```text
W^{deg H} = c * H^{deg W},   i.e.   W^g = c H^{n-1}    (for W != 0).    (W2)
```

Verified in both directions on samples: a common linear form and a common
quadratic satisfy both sides; a non-proportional pair fails both.

## 3. The arithmetic obstruction (W3)

Here is the point. Write `H = prod_i l_i^{e_i}` with distinct linear forms
`l_i` and `sum_i e_i = g`. By (W2), `g e'_i = e_i (n-1)` where `W = c' prod
l_i^{e'_i}`, so `g | e_i (n-1)` for every `i`. But

```text
gcd(g, n-1) = gcd(g, ga-1) = gcd(g, -1) = 1        for ALL g, a >= 1,
```

(checked exhaustively for `1 <= g, a <= 39`: no exceptions). Hence `g | e_i`
for every `i`. Since every `e_i >= 1` and `sum e_i = g`, there is **exactly
one** factor, with `e_1 = g`:

```text
W != 0   =>   H = l^g   for a single linear form l.                    (W3)
```

Then `P_n = c l^{ga}` and `Q_m = c' l^{gb}` have a single root direction:
`K = 1`, the **one-place-at-infinity** configuration, which is classically
resolved (Abhyankar-Moh). We cite that; we do not reprove it.

## 4. The conclusion (W4)

```text
K >= 2   =>   W = 0,   unconditionally.
```

A JC(2) counterexample requires `K >= 2` (HYP-9070's gate; `K = 1` is the
resolved case). So **`W` vanishes identically on the entire counterexample
locus.**

Direct confirmation -- for `H` with at least two distinct linear factors, the
only homogeneous solution of `Jac(W,H) = 0` in degree `n-1 = ga-1` is `W = 0`:

```text
H = (x+y)(x-y),        g=2, a=2, deg W=3,   K=2  ->  only W = 0
H = (x+y)(x-y),        g=2, a=3, deg W=5,   K=2  ->  only W = 0
H = (x+y)^2(x-y),      g=3, a=2, deg W=5,   K=2  ->  only W = 0
H = x y (x+y),         g=3, a=3, deg W=8,   K=3  ->  only W = 0
H = (x+y)(x-y)(x+2y),  g=3, a=4, deg W=11,  K=3  ->  only W = 0
```

**The dichotomy is sharp, not vacuous.** At `K = 1` the solution space really
is nonzero:

```text
H = (x+y)^2, deg W = 3:  W = c (x+y)^3
H = (x+y)^3, deg W = 5:  W = c (x+y)^5
```

so (W3) is an exact characterisation of when `W` may survive, not a
one-sided bound.

## 5. Consequence for the tower (W5)

THM-3016 section 4b(B) derived, **assuming `W = 0`**, that
`Jac(P_{n-1},Q_{m-1}) = kappa(a-b) H^{a-b-1} Q_{m-1} J`, hence that in the
degree `n+m-4` relation

```text
c a H^{a-1} Jac(H, Q_{m-2}) + Jac(P_{n-1},Q_{m-1}) + c' b H^{b-1} Jac(P_{n-2},H) = 0
```

all three terms carry explicit powers of `H` (exponents `a-1`, `a-b-1`,
`b-1`), so dividing by `H^{min(a-b-1,b-1)}` produces the next relation with a
Euclidean-reduced exponent pair.

That derivation was conditional. **It is now unconditional on the
counterexample locus**: by (W4) every counterexample candidate has `W = 0`.
The subleading tower therefore runs without a case split, and the Euclidean
reduction on `(a,b)` is available at every step -- which is exactly the
structure HYP-9070's Euclidean-depth search order was built to exploit.

## 6. Scope

(W1)-(W4) are proofs. (W2) is the classical binary-form fact, verified here
rather than cited. (W3) is elementary arithmetic, checked exhaustively in a
finite box and proved in general by `gcd(g, ga-1) = 1`. The appeal to
Abhyankar-Moh for `K = 1` is a **citation**, not a result of this file, and
the `K >= 2` gate is HYP-9070's proposal, not a theorem -- so the honest
reading of (W4) is: *`W = 0` holds except in the one-place-at-infinity
configuration.* **Nothing here decides JC(2)**, and no bridge or reduction to
another repo lane is claimed (cf. MISTAKE-237). What is gained is the removal
of a hypothesis from THM-3016's tower and an exact characterisation of the
only configuration in which `W` can be nonzero.

Referee: `jc2_W_forced_thm3025.py` -- the (W2) equivalence on samples, the
exhaustive `gcd(g, ga-1) = 1` check, the five `K >= 2` solution-space
computations, and the two `K = 1` controls.
