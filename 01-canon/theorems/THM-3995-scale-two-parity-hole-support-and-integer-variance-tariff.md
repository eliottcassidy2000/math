---
id: THM-3995
title: "Scale-two parity holes shrink support and strengthen the integer variance tariff"
status: >
  PROVED IN THE THM-3878/3910 CONDITIONAL SLICE + VERIFIED-EXACT +
  INDEPENDENTLY HOSTILE-AUDITED. For the sole scale-two (2,1,9) survivor with
  odd t>=U, parity forbids the correctly oriented body-wall event in four
  one-sided intervals just inside the exact quotient obstruction. Retaining
  the owner bound u<=U, rather than weakening it to u<=t, shows that any
  failure has sheet-count support at most s_U=4(U-1)/(63U). A sharp integer
  support-envelope lemma gives the resulting variance, discrepancy, and BV
  gates, including the k=0 plateau. It does not close the survivor, the t<U
  branch, or LRC(14).
source: root + lrc_incidence_filter / Hopf repair-response transfer, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS, INCLUDING BODY-MAXIMUM REFINEMENT
  (2026-08-24). The original independent wall audit checked masking, null
  walls, disjoint hole widths, the integer envelope, dependencies, and scope.
  A second exact audit reconstructed all oriented wall formulas, swept 575,664
  owner events for 11<=U<=50 and odd t satisfying U<=t<4U/3, tested even-t and masking
  hostiles, separated eventwise from global sharpness, checked an asymptotic
  U/U-1 equality family, the k=0 plateau, and the discrepancy/BV algebra.
  Normal and optimized outputs LF-normalized byte-match the frozen transcript.
depends_on:
  - THM-3878-lrc14-eleven-plus-two-harmonic-absorption-seam-collapse
  - THM-3910-lrc14-auxiliary-center-erosion-and-t-sheet-variance-response
related:
  - THM-3878-lrc14-eleven-plus-two-harmonic-absorption-seam-collapse
  - THM-3910-lrc14-auxiliary-center-erosion-and-t-sheet-variance-response
  - THM-3990-componentwise-harmonic-obstruction-and-repair-quotient
script: 04-computation/lrc14_scale_two_parity_hole_tariff_thm3995.py
output: 05-knowledge/results/lrc14_scale_two_parity_hole_tariff_thm3995.out
script_sha256: 4249cc9cc0597482d7634969d5e1c6a5f453a9c1f34ad8ad971517a986055b99
output_sha256: 3d7e6bd375231d3d6d179276d9e715c808c2ed7fee9016ecfdf50fdaf90b3e33
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
`C`. Keeping the actual body maximum `U` strengthens the support and variance
invoices without assuming equality in THM-3910's old tariff.

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

Since every body owner satisfies `u<=U`, the same lower bounds hold with `u`
replaced by `U`. This is the load-bearing refinement: the earlier `t`-based
version discarded the available inequality `U<=t`. Oddness is required of
`t`, not of `U`. If even `t` were admitted, the legal arithmetic control
`(U,t,u,a)=(11,12,9,0)` would put an entering wall exactly at `2/21` and an
exiting wall exactly at `19/21`, so the corresponding holes could vanish.

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
the two right endpoints. Simultaneous masking is possible--for example, owner
`1` enters while owner `13` exits at `y=1/14`--but it cannot remove the need
for an entering event at the first positive jump or an exiting event at the
last positive jump. Isolated wall values have measure zero.

The body has eleven distinct positive speeds, so `U>=11`. In particular, the
four intervals are disjoint: in either obstruction component their two widths
sum to `2/(63U)<2/63`. Their total measure is at least

```text
1/(42U)+1/(126U)+1/(126U)+1/(42U)=4/(63U).              (6)
```

More precisely, modulo the null set of owner walls, the positive support lies
in the two **closed** trimmed cores

```text
K_U=[2/21+1/(42U), 8/63-1/(126U)]
    union [55/63+1/(126U), 19/21-1/(42U)].              (7)
```

Closed endpoints, or the phrase *modulo null walls*, are necessary: an event
can occur at equality and the point value of the closed-safe indicator can
survive there. Consequently, for the essential positive support
`S={N_t>0}` of any failure,

```text
meas(S)<=s_U:=4/63-4/(63U)=4(U-1)/(63U).                (8)
```

Thus `s_U<=s_t:=4(t-1)/(63t)`, with equality exactly when `t=U`. If `U` is
even, oddness of `t` forces `t>U`, hence strict improvement; if `U` is odd,
`t=U` is allowed and gives no improvement. In every case `s_U<4/63`.

The four constants in `(4)` are individually attained in the allowed
arithmetic range. For rows `(3a)--(3d)`, examples `(U,t,u,a,M)` are

```text
(44,45,44,1,1), (40,41,40,8,8),
(40,41,40,32,32), (44,45,44,43,43),                    (9)
```

with numerators `1,-1,1,-1`, respectively. They cannot all be attained by the
same owner `u=U`: the outer pair would require `3t-4U=1 mod 42`, while the
inner pair would require `9t+16U=1 mod 126`; three times the first congruence
would force `28U=-2 mod 126`, impossible modulo `14`. Hence no exact global
support equality is claimed.

The coefficient in `(6)` is nevertheless asymptotically sharp at the
**oriented-event arithmetic layer**. For

```text
U=126n+122,       t=U+13,
n>=0,             n not congruent to 11 modulo 13,
```

the outer walls attain numerator magnitude one with owner `U`, and the inner
walls do so with owner `U-1`. Their total spacing is

```text
(1+1/(4(U-1)))*4/(63U),
```

whose ratio to `(6)` tends to one. This family neither realizes equality in
the support cap nor supplies an LRC failure. Indeed,
`3t-4U=1 mod 42` and `9t+16(U-1)=1 mod 126`. The exclusion on `n` gives
`gcd(t,U)=gcd(13,U)=1`, while `gcd(t,U-1)=gcd(14,U-1)=1`, so the associated
linear Diophantine equations in `a,M` are soluble. Also `t` is odd and
`U<=t<4U/3`. Taking the other nine body owners to be `1,...,9` cannot insert
a closer oriented wall, because their individual lower bounds in `(4)` are
larger than those for `U-1` and `U`.

## 3. A sharp integer support-envelope lemma

More generally, let `N` be any nonnegative integer-valued `L^2` function on a
probability space, with finite mean `m` and positive support of measure at
most `s in (0,1]`. Equivalently, one may read the inequalities below in the
extended-real sense. Put

```text
k=floor(m/s),                 theta=m/s-k.             (10)
```

For every nonnegative integer `n`, adjacency of `k,k+1` gives the pointwise
inequality

```text
n^2 >= (2k+1)n-k(k+1)*1_(n>0).                         (11)
```

Indeed, for `n>0` its difference is `(n-k)(n-k-1)>=0`; there is no integer
strictly between `k` and `k+1`. Integrating `(11)` and using
`meas(N>0)<=s` yields

```text
 integral N^2
 >=(2k+1)m-k(k+1)s
 =m^2/s+s*theta*(1-theta).                             (12)
```

After subtracting `m^2`,

```text
Var(N)>=m^2*(1-s)/s+s*theta*(1-theta).                 (13)
```

On a nonatomic probability space, including `T`, the bound is sharp. If
`k>=1`, equality requires support of measure exactly `s` and the adjacent
values `k,k+1` in proportions `1-theta,theta`. If `k=0`, `(12)` reduces to
`integral N^2>=m`; equality is attained by `N=1` on a set of measure `m<=s`,
so the positive support may be strictly smaller than `s`. These are abstract
envelope equalities, not realizability claims for a geometric lift count.

Apply `(13)` to `(8)` and write

```text
theta_(t,U)={m/s_U}.
```

Every scale-two failure in the stated slice satisfies

```text
Var(N_t)>=m^2*(1-s_U)/s_U+s_U*theta_(t,U)*(1-theta_(t,U)),
s_U=4(U-1)/(63U).                                      (14)
```

As a function of the permitted support cap, the exact integer envelope is
constant at `integral N^2=m` while `m<=s`. Therefore, even when `t>U` and
`s_U<s_t`, replacing `s_t` by `s_U` gives no exact envelope gain for
`m<=s_U`, including the boundary `m=s_U`; it gives a strict gain exactly when
`m>s_U`. This `k=0` plateau must not be advertised as a strict improvement.

## 4. The strengthened BV gate and exact residual

If `G` has `r_G` positive-length components, THM-3910 gives

```text
Var(N_t)<=r_G^2/3.                                     (15)
```

Moreover

```text
(1-s_U)/s_U=(59U+4)/(4(U-1)).                          (16)
```

In THM-4002's notation,

```text
disc_t(G)=sum_(k not=0)|Ghat(kt)|^2=Var(N_t)/t^2.
```

Thus every failure obeys the exact necessary discrepancy tariff

```text
disc_t(G)>=meas(G)^2*(59U+4)/(4(U-1))
       +(s_U/t^2)*theta_(t,U)*(1-theta_(t,U)).          (17)
```

The strict reverse inequality is sufficient for positive full safe mass;
equality is not. Combining `(15)` with the strict reverse of `(17)` gives the
full BV-plus-integer sufficient gate

```text
r_G^2/(3t^2)<meas(G)^2*(59U+4)/(4(U-1))
              +(s_U/t^2)*theta_(t,U)*(1-theta_(t,U)).  (18)
```

Dropping only the nonnegative integer correction gives the simpler sufficient
condition

```text
t*meas(G)/r_G > sqrt(4(U-1)/(3(59U+4))).                (19)
```

The integer term in `(17)` is positive exactly when `m/s_U` is nonintegral.
Relative to the real-valued `s_U` tariff this is then a strict improvement;
relative to the old `s_t` integer tariff, the `k=0` plateau described above
still applies.

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
