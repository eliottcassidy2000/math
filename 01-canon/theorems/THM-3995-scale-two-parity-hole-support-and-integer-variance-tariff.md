---
id: THM-3995
title: "Scale-two parity holes shrink support and strengthen the integer variance tariff"
status: >
  PROVED IN THE THM-3878/3910 CONDITIONAL SLICE + VERIFIED-EXACT +
  INDEPENDENTLY HOSTILE-AUDITED. For the sole scale-two (2,1,9) survivor with
  odd t>=U, parity forbids the correctly
  oriented body-wall event in four one-sided intervals just inside the exact
  quotient obstruction. Any failure therefore has sheet-count support at
  most s_t=4(t-1)/(63t), rather than 4/63. A sharp integer support-envelope
  lemma gives Var(N_t)>=m^2(1-s_t)/s_t+s_t*theta_t*(1-theta_t), including
  k=0 and failures away from the old tariff equality. Combined with the BV
  upper bound this gives the sufficient gate displayed below. It does not
  close the survivor, the t<U branch, or LRC(14).
source: root + lrc_incidence_filter / Hopf repair-response transfer, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (thm3995_hostile_audit, 2026-08-24). An
  independent clockwise wall sweep checked every orientation for odd t<=51
  and u<=t; the proof audit checked masking, null walls, disjoint hole widths,
  the support cap, the integer envelope including k=0, the BV algebra, hashes,
  dependencies, and scope. Normal and optimized outputs LF-normalized
  byte-match the frozen transcript.
depends_on:
  - THM-3878-lrc14-eleven-plus-two-harmonic-absorption-seam-collapse
  - THM-3910-lrc14-auxiliary-center-erosion-and-t-sheet-variance-response
related:
  - THM-3878-lrc14-eleven-plus-two-harmonic-absorption-seam-collapse
  - THM-3910-lrc14-auxiliary-center-erosion-and-t-sheet-variance-response
  - THM-3990-componentwise-harmonic-obstruction-and-repair-quotient
script: 04-computation/lrc14_scale_two_parity_hole_tariff_thm3995.py
output: 05-knowledge/results/lrc14_scale_two_parity_hole_tariff_thm3995.out
script_sha256: 342ca76c34bb09d2f19f0237d973edb8872d2cf43c8c637a6efbfc99804ff8c6
output_sha256: 2973b970aa75a6234f5aca3af02293476a35d9d64cee7b7dda9a83f3d6cf6a9c
hash_basis: raw LF bytes
---

# THM-3995 -- scale-two parity holes and the integer variance tariff

**PROVED IN THE THM-3878/3910 CONDITIONAL SLICE + VERIFIED-EXACT +
INDEPENDENTLY HOSTILE-AUDITED.** Retain all notation and hypotheses of those
theorems. In particular, work on
`T=R/Z` with

```text
D_u={x: ||ux||<1/14},
G=intersection_i D_(u_i)^c,       U=max_i u_i,
N_t(w)=sum_(a=0)^(t-1) 1_G((w+a)/t),
m=integral_T N_t=t*meas(G).                                  (1)
```

Consider only the conditional `t>=U` scale-two survivor
`(s,p,q)=(2,1,9)`. Coprimality makes `t` odd. THM-3878 computes its exact
open quotient obstruction as

```text
C=(2/21,8/63) union (55/63,19/21),       meas(C)=4/63.  (2)
```

The theorem proves that a failure cannot use four one-sided subintervals of
`C`. This strengthens the support and variance invoices without assuming
equality in THM-3910's old tariff.

## 1. The four correctly oriented event numerators are odd

An individual owner-`u` safe interval is entered at

```text
y=(14a+1)/(14u)
```

and exited at `y=(14a-1)/(14u)`. Multiplication by `t` preserves orientation
on the circle. To turn a zero sheet count just outside `C` into a positive
count inside its left endpoint, an entering event is necessary. To return to
zero across a right endpoint, an exiting event is necessary.

Clear denominators between these correctly oriented owner walls and the four
successive boundaries in `(2)`. For arbitrary integers `a,M`, the four
difference numerators are

```text
3t(14a+1)-4u-42uM,        denominator 42u,              (3a)
9t(14a-1)-16u-126uM,      denominator 126u,             (3b)
9t(14a+1)-110u-126uM,     denominator 126u,             (3c)
3t(14a-1)-38u-42uM,       denominator 42u.              (3d)
```

Every subtracted term is even, while `t` and `14a+/-1` are odd. Hence all
four numerators are odd and in particular nonzero. The nearest possible
correctly oriented events have respective distances at least

```text
1/(42u), 1/(126u), 1/(126u), 1/(42u).                  (4)
```

Since every body owner satisfies `u<=U<=t`, the same lower bounds hold with
`u` replaced by `t`.

## 2. Every failure has four interior parity holes

Let `H_2=T minus C`. By THM-3910, the full safe mass is positive exactly when
`integral_(H_2) N_t>0`. Thus a failure has

```text
N_t=0 almost everywhere on H_2.                         (5)
```

Between owner-wall images, `N_t` is constant. If it becomes positive after a
left endpoint of `C`, its first positive jump contains at least one entering
owner event; simultaneous masked events do not change this necessity. The
first and third bounds in `(4)` therefore force zero intervals inside the two
left endpoints. Reading from the right, the last positive cell would require
an exiting event, so the second and fourth bounds force zero intervals inside
the two right endpoints. Isolated wall values have measure zero.

The four intervals are disjoint in the inherited `t>1` range. Their total
measure is at least

```text
1/(42t)+1/(126t)+1/(126t)+1/(42t)=4/(63t).              (6)
```

Consequently, for the essential positive support `S={N_t>0}` of any failure,

```text
meas(S)<=s_t:=4/63-4/(63t)=4(t-1)/(63t).                (7)
```

This is strictly smaller than the quotient-obstruction measure for every
finite allowed `t`.

## 3. A sharp integer support-envelope lemma

More generally, let `N` be any nonnegative integer-valued `L^2` function on a
probability space, with finite mean `m` and positive support of measure at
most `s in (0,1]`. Equivalently, one may read the inequalities below in the
extended-real sense. Put

```text
k=floor(m/s),                 theta=m/s-k.              (8)
```

For every nonnegative integer `n`, adjacency of `k,k+1` gives the pointwise
inequality

```text
n^2 >= (2k+1)n-k(k+1)*1_(n>0).                          (9)
```

Indeed, for `n>0` its difference is `(n-k)(n-k-1)>=0`; there is no integer
strictly between `k` and `k+1`. Integrating `(9)` and using
`meas(N>0)<=s` yields

```text
integral N^2
 >=(2k+1)m-k(k+1)s
 =m^2/s+s*theta*(1-theta).                              (10)
```

After subtracting `m^2`,

```text
Var(N)>=m^2*(1-s)/s+s*theta*(1-theta).                  (11)
```

On a nonatomic probability space, including `T`, the bound is sharp:
distribute the mass on `s` measure between the adjacent values `k,k+1` in
proportions `1-theta,theta` (when `k=0`, the positive support may be smaller).
In particular `(11)` remains valid at `k=0`; no exposed-endpoint equality
assumption is present.

Apply `(11)` to `(7)` and write

```text
theta_t={m/s_t}.
```

Every scale-two failure in the stated slice satisfies

```text
Var(N_t)>=m^2*(1-s_t)/s_t+s_t*theta_t*(1-theta_t),
s_t=4(t-1)/(63t).                                      (12)
```

## 4. The strengthened BV gate and exact residual

If `G` has `r_G` positive-length components, THM-3910 gives

```text
Var(N_t)<=r_G^2/3.                                     (13)
```

Moreover

```text
(1-s_t)/s_t=(59t+4)/(4(t-1)).                          (14)
```

Dropping only the nonnegative integer correction in `(12)`, a sufficient
condition for positive full safe mass is therefore

```text
t*meas(G)/r_G > sqrt(4(t-1)/(3(59t+4))).                (15)
```

The full formula `(12)` is stronger whenever `m/s_t` is nonintegral.

This theorem uses the exact scale-two obstruction, odd `t`, and `t>=U`. It
does not control `r_G/meas(G)` uniformly, so it does not close `(2,1,9)`.
It says nothing about `t<U`, the sixteen scale-one types, or other component
shapes, and it does not prove `LRC(14)`. Its gain is that the scale-two row now
has a parity-sensitive tariff valid for arbitrary failures, not only the old
abstract equality locus. **QED.**

## Reproduction

```bash
python3 04-computation/lrc14_scale_two_parity_hole_tariff_thm3995.py
python3 -O 04-computation/lrc14_scale_two_parity_hole_tariff_thm3995.py
sha256sum 04-computation/lrc14_scale_two_parity_hole_tariff_thm3995.py \
  05-knowledge/results/lrc14_scale_two_parity_hole_tariff_thm3995.out
python3 agents/check_docs.py
```
