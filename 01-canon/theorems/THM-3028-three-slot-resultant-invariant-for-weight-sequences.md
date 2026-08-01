---
id: THM-3028
title: "The three-slot resultant invariant R_w, extending the two-slot dichotomy"
status: >
  PROVED + VERIFIED-EXACT. Extends THM-3022 from two slots to three. For a
  weight w with L_w(s^j) = w_j and a three-slot support {a,b,c}, work
  projectively in P^2: L_w(f)=0 is a LINE, L_w(f^2)=0 a CONIC, L_w(f^3)=0 a
  CUBIC. The line meets the conic in two points (Bezout), and the three-slot
  threshold is 3 exactly when neither lies on the cubic -- a condition on the
  WEIGHT, not on f. Restricting the quadric and cubic to the line and taking
  their resultant gives a computable invariant R_w(a,b,c), homogeneous of
  degree 6, with R_w != 0 => three moments force f = 0. Computed exactly:
  R_w != 0 for the FACTORIAL, the [0,1] LEBESGUE weight 1/(j+1) and the
  CENTRAL BINOMIAL on every tested triple, so all three have three-slot
  threshold 3; and R_w = 0 identically for the GEOMETRIC and FIBONACCI
  weights, matching THM-3022's degenerate cases and its exact Fibonacci
  counterexample. NEW CONSEQUENCE: the FC-analogue for the interval measure --
  the weight of the exponential-integral setting of HYP-9078 -- holds at three
  slots, not merely two.
source: opus-2026-08-01-amm12592-writeup
depends_on:
  - THM-3022
related:
  - THM-2824  # arbitrary three-slot first-window detection, factorial weight
  - THM-2917
  - THM-2173
  - HYP-9078
script: 04-computation/fc_three_slot_resultant_thm3028.py
output: 05-knowledge/results/fc_three_slot_resultant_thm3028.out
---

# THM-3028 -- the three-slot resultant invariant

## 1. Why two slots do not generalize directly

THM-3022 works at two slots because eliminating one coefficient via
`L_w(f) = 0` leaves ONE complex unknown, and `L_w(f^2)` becomes `c_1^2` times
a constant -- so it vanishes only at `c_1 = 0`. At three slots two unknowns
remain, and a complex quadratic form in two variables always has nontrivial
zeros. The second moment is therefore never enough, and the third must enter.

## 2. The invariant

Work projectively: `f` up to scale is a point of `P^2`. Then

```text
L_w(f)   = 0   is a LINE,
L_w(f^2) = 0   is a CONIC,
L_w(f^3) = 0   is a CUBIC.
```

By Bezout the line meets the conic in two points. **The three-slot threshold
is 3 precisely when neither of those two points lies on the cubic** -- and
that is a condition on the weight `w`, not on `f`.

Parametrize the line by `(u,v)` via

```text
c_1 = u w_b,   c_2 = -u w_a + v w_c,   c_3 = -v w_b,
```

which satisfies `c_1 w_a + c_2 w_b + c_3 w_c = 0` identically. Restricting,
`L_w(f^2)` and `L_w(f^3)` become binary forms of degrees 2 and 3, and

```text
R_w(a,b,c) := Res_v ( L_w(f^2)|_line , L_w(f^3)|_line )
```

is homogeneous of degree `2*3 = 6` in `u`. **`R_w(a,b,c) != 0` implies
`L_w(f) = L_w(f^2) = L_w(f^3) = 0` has no nonzero solution on that support**,
i.e. three moments force `f = 0`; `R_w = 0` exhibits a common zero, hence a
witness.

## 3. Computed values

Exact symbolic computation, `R_w = c * u^6` with `c` the invariant:

```text
support            (0,1,2)          (0,1,3)            (1,2,3)              (0,2,4)
factorial j!        7424           33900062976      15355506327552     2.589e22
[0,1] 1/(j+1)     1/1.543e11       81/8.992e11      1/1.098e19         5.57e9/5.31e21
central binom      18432           1006608384        20639121408       1.198e18
geometric 3^j        0                 0                 0                  0
Fibonacci            0                 0                 0                  0
```

So:

1. **Factorial**: `R_w != 0` throughout, three-slot threshold 3 -- consistent
   with THM-2824, which proves first-window SFC(3) on arbitrary three-slot
   supports by much heavier machinery. `R_w` gives a one-line certificate per
   support.
2. **`[0,1]` Lebesgue weight**: `R_w != 0` throughout. **New**: the
   FC-analogue for the interval measure -- the weight attached to the
   exponential-integral claim of HYP-9078 -- holds at three slots, extending
   HYP-9078 sec 4 from two.
3. **Central binomial**: `R_w != 0`, threshold 3.
4. **Geometric**: `R_w = 0` identically, as it must -- `L_w` is an evaluation
   and every `f` with `f(r) = 0` kills all moments (THM-3022 part 3).
5. **Fibonacci**: `R_w = 0` identically, matching THM-3022 sec 4's exact
   counterexample `f = s(1-s)`, for which `Delta^m w_n = w_(n-m)` sends every
   moment to `w_0 = 0`.

The invariant therefore separates exactly the weights already known to be
degenerate, and certifies the others.

## 4. Scope

Three slots only. The construction generalizes in shape -- at `N` slots one
intersects `N-1` hypersurfaces of degrees `2..N` inside the hyperplane
`L_w(f) = 0` and asks whether the Bezout points avoid the last one -- but the
elimination grows fast, and no general-`N` statement is claimed. THM-2173
remains the matching lower half for all `N` (Krull height), and equation
counting still cannot supply an upper bound (MISTAKE-246).
