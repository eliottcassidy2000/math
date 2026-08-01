---
id: THM-3022
title: "The two-slot moment dichotomy for weight sequences, and the sharpness of the factorial"
status: >
  PROVED + VERIFIED-EXACT. For a positive weight sequence w and the functional
  L_w(s^j) = w_j, the two-slot support {a,b} carries the invariant
  Q_w(a,b) = w_2a w_b^2 - 2 w_a w_b w_(a+b) + w_a^2 w_2b. DICHOTOMY:
  Q_w(a,b) != 0 forces f = 0 from L_w(f) = L_w(f^2) = 0 alone (two moments
  suffice), while Q_w(a,b) = 0 makes the second moment vacuous. Q_w vanishes
  identically iff w is geometric, i.e. iff L_w is an evaluation. COROLLARY
  (the factorial case, any number of variables, arbitrary distinct
  multi-indices): coordinatewise log-convexity of Gamma gives
  (alpha+beta)!^2 <= (2alpha)!(2beta)! strictly off the diagonal, and AM-GM
  then gives Q > 0, so two moments force f = 0. This is a short, positivity-
  only proof of two-slot SFC, previously carried externally (Edo-van den
  Essen) or by resultant/ODE arguments (THM-3019/3020). SHARPNESS: the bound
  is not soft. For any w obeying the Fibonacci recurrence, f = s(1-s) gives
  L_w(f^m) = (-1)^m w_0 for EVERY m, because Delta w_n = w_(n-1) makes
  Delta^m w_n = w_(n-m); so the FIBONACCI weight (w_0 = 0) is an exact
  two-slot counterexample with ALL moments vanishing, while LUCAS (w_0 = 2)
  gives exactly +-2 at every moment. What the factorial supplies is precisely
  that its differences never return into the sequence: Delta(n!) = n * n!.
source: opus-2026-07-31-amm12592-writeup
depends_on: []
related:
  - THM-2173  # sparse projective factorial moment floor: the matching lower half, by Krull height
  - THM-2824  # arbitrary three-slot first-window detection
  - THM-3010  # ballot-column Newton ratios: Q_w is a Newton/Turan ratio
  - THM-3019
  - THM-3020
  - HYP-9076
external:
  - "E. Edo, A. van den Essen, The strong factorial conjecture (two slots)."
script: 04-computation/fc_two_slot_threshold_proof.py
output: 05-knowledge/results/fc_two_slot_threshold_proof.out
---

# THM-3022 -- the two-slot moment dichotomy

## 1. Setting

Let `w = (w_j)_(j>=0)` be positive reals and `L_w(s^j) = w_j`, extended
linearly. For a two-slot `f = c_1 s^a + c_2 s^b` with `a != b`, define

```text
Q_w(a,b) = w_(2a) w_b^2 - 2 w_a w_b w_(a+b) + w_a^2 w_(2b).          (1)
```

## 2. The dichotomy

**Theorem.** If `L_w(f) = 0` and `f != 0` then `c_1 != 0`,
`c_2 = -c_1 w_a/w_b`, and

```text
L_w(f^2) = (c_1^2 / w_b^2) * Q_w(a,b).                               (2)
```

Hence:

1. `Q_w(a,b) != 0`  =>  `L_w(f) = L_w(f^2) = 0` forces `f = 0`. Two moments
   suffice on that support. (Note the SIGN of `Q_w` is irrelevant: over `C`,
   `c_1^2 Q = 0` with `Q != 0` real already gives `c_1 = 0`.)
2. `Q_w(a,b) = 0`  =>  the second moment is vacuous on that support and the
   threshold exceeds two.
3. `Q_w(a,b) = 0` for ALL `a != b`  <=>  `w_j = w_0 r^j` is geometric. Then
   `L_w(f^m) = f(r)^m`, `L_w` is an evaluation, and every `f` with `f(r) = 0`
   kills all moments at once.

*Proof.* (2) is direct substitution. For (3): geometric `w` gives
`Q = r^(2a+2b)(1 - 2 + 1) = 0`; conversely `Q_w(a,a+1) = 0` for all `a` is the
recursion `w_(2a) w_(a+1)^2 - 2 w_a w_(a+1) w_(2a+1) + w_a^2 w_(2a+2) = 0`,
which forces the ratio `w_(j+1)/w_j` constant. QED

## 3. The factorial corollary

**Corollary.** Let `f = c_1 x^alpha + c_2 x^beta` with `alpha != beta`
multi-indices in ANY number of variables, `L(x^gamma) = gamma!`. Then `L(f)`
and `L(f^2)` cannot both vanish.

*Proof.* Coordinatewise **log-convexity of Gamma** gives
`((alpha_i+beta_i)!)^2 <= (2alpha_i)!(2beta_i)!` with equality iff
`alpha_i = beta_i`; multiplying over `i`, `C^2 <= D E` strictly unless
`alpha = beta`. **AM-GM** gives `D B^2 + A^2 E >= 2 A B sqrt(D E)`. So

```text
Q = D B^2 - 2 A B C + A^2 E >= 2 A B ( sqrt(D E) - C ) > 0,
```

and part 1 applies. QED

Verified on `29112` multi-index pairs in 1, 2, 3 variables: no violation, and
`Q = 0` exactly in the excluded case `alpha = beta`.

`C^2 <= D E` is a **Turan-type inequality** for the factorial sequence; the
repo contains no other Turan inequality of this kind, and `Q_w` is a Newton
ratio in the sense of THM-3010 (`R_k = h_k^2/(h_(k-1)h_(k+1))`), evaluated at
the midpoint of `[2alpha, 2beta]`.

### A free certificate for translated supports

For `alpha = n+a`, `beta = n+b` in one variable, `d = b-a >= 1`,
`N_0 = 2n+2a`:

```text
C^2/(D E) = prod_(i=1)^(d) (N_0+i)/(N_0+d+i)  <  1,
```

each factor `< 1`, so positivity holds for ALL translations at once with no
resultant and no case analysis -- the statement THM-2922 needs a 197-term
Gregory--Newton certificate for at four slots is free at two.

## 4. Sharpness: the factorial is doing real work

The corollary is not soft. Let `w` obey the Fibonacci recurrence
`w_(n+1) = w_n + w_(n-1)` and take `f = s(1-s)`. Then
`Delta w_n = w_(n+1) - w_n = w_(n-1)`: **the difference operator is the
backward shift**, so `Delta^m w_n = w_(n-m)` and

```text
L_w(f^m) = sum_j binom(m,j)(-1)^j w_(m+j) = (-1)^m w_0   for every m >= 1.
```

Therefore:

```text
FIBONACCI (w_0 = 0):  L_w(f^m) = 0 for EVERY m -- an exact two-slot
                      counterexample; the threshold is infinite.
LUCAS     (w_0 = 2):  L_w(f^m) = +-2 for every m -- never zero.
```

Both verified to `m = 10`. So two-slot detection genuinely depends on the
weight, and what the factorial supplies is exactly that its differences never
return into the sequence: `Delta(n!) = (n+1)! - n! = n * n!`.

## 5. Scope

This settles two slots only. THM-2173 supplies the matching lower half for all
slot counts (a `t`-slot envelope always admits a nonzero witness killing
`t-1` moments, by Krull height). The upper half beyond two slots is the open
SFC program: THM-2824 at three slots, THM-2908/2921/2922/2940/2849 at four in
pieces, nothing at five. **No general "N slots need N moments" theorem exists,
and equation counting cannot supply one** -- that inference is MISTAKE-246.
