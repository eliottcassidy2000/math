---
id: THM-4038
title: "AP-cover deficit holonomic sixty-phase law"
status: >
  PROVED + FINITE-EXACT PHASE AUDIT. After the THM-4029 twelve-owner
  reduction, the AP seven-sector deficit is an exact 60-periodic rational
  coefficient law in six affine counters. Clearing those counters gives a
  degree-five quasipolynomial of exact period 60, hence an explicit
  polynomial-coefficient recurrence and a D-finite but nonalgebraic ordinary
  generating function. The phase first appears at second order. This is a
  theorem about the consecutive AP cover and does not prove AP extremality or
  LRC(14).
source: codex sun-minimal-sturmian / Sturmian-triangular-Fibonacci audit, 2026-08-24
audit: >
  PASS. The exact companion reconstructs all sixty rational laws from
  THM-4029, aggregates the six shifts and six owner-denominator clocks, tests
  every divisor of 60 for minimality, clears the common denominator, checks
  the sixth 60-step recurrence through n=499, and verifies every triangular
  track identity for q<=6, all phases, and 1<=K<=29.
depends_on:
  - THM-4029-lrc14-ap-cover-twelve-owner-rational-tail
  - THM-536-lrc-seven-sector-sturmian-partial-sum-reframe
related:
  - THM-637
  - THM-4028
script: 04-computation/lrc14_ap_cover_holonomic_sixty_phase_thm4038.py
output: 05-knowledge/results/lrc14_ap_cover_holonomic_sixty_phase_thm4038.out
---

# THM-4038 -- the AP deficit is holonomic with an exact sixty-phase law

**PROVED + FINITE-EXACT PHASE AUDIT.** This theorem gives the positive scalar
classification left open by THM-4029's stopping results. The AP deficit is not
eventually periodic, quasipolynomial, or C-finite. It is nevertheless a
finite periodic family of rational functions, becomes quasipolynomial after
one fixed denominator clearing, and is P-recursive. The phase and the
unbounded affine counter are both load-bearing.

## 1. Object, normalization, and inherited owner law

Use THM-4029's consecutive-AP cover sequence

```text
a(m)=meas{x in [0,1):
  {floor(7 e x) mod 7:0<=e<m}=Z/7},
D(m)=1-a(m).                                             (1)
```

Put `n=m-1`. THM-4029 proves that for every `n>=12`, the twelve rational
owners and their two sides give finitely many terms `C/(n-c)`, where
`0<=c<=5`, and that the selected terms depend only on `n mod 60`. For
`r in {0,...,59}` and `0<=c<=5`, let `A_c(r)` be the sum of the coefficients
`C` of all selected terms having shift `c` on phase `r`. Then the exact
normalized law is

\[
 \boxed{D(n+1)={1\over7}\sum_{c=0}^{5}{A_c(n\bmod60)\over n-c}}
 \qquad(n\ge12).                                      \tag{2}
\]

The factor `1/7` is the change from theta-length on `[0,7)` to x-measure on
`[0,1)`. It must not be absorbed silently into the phase coefficients. In
every phase,

\[
 A_c(r)\ge0,
 \qquad \sum_{c=0}^{5}A_c(r)={127\over5},              \tag{3}
\]

so `(2)` recovers THM-4029's leading constant `127/35` after the `1/7`
normalization.

## 2. Exact phase minimality and its denominator owners

**FINITE-EXACT:** the companion evaluates every `A_c(r)` in rational
arithmetic and tests every divisor of 60. The six coefficient sequences have
the following exact statistics.

| shift `c` | minimal period | nonzero phases | distinct values |
|---:|---:|---:|---:|
| 0 | 60 | 60 | 12 |
| 1 | 60 | 60 | 12 |
| 2 | 60 | 60 | 12 |
| 3 | 60 | 57 | 8 |
| 4 | 30 | 52 | 4 |
| 5 | 6 | 20 | 2 |

The full vector

```text
A(r)=(A_0(r),...,A_5(r))                                (4)
```

has minimal period exactly `60`, and its sixty phase vectors are pairwise
distinct. Thus “sixty-phase” is not a fit or merely an upper bound.

Resolve `(4)` instead by the denominator `q` of the rational owner. The
aggregate coefficient vector contributed by denominator `q` has exact
minimal period and exact number of distinct states

```text
q                         1  2  3  4  5  6
minimal period            1  2  3  4  5  6
distinct coefficient rows 1  2  3  4  5  6.            (5)
```

Consequently the global clock is precisely the synchronization

\[
 \operatorname{lcm}(1,2,3,4,5,6)=60.                 \tag{6}
\]

### Fibonacci/Ostrowski hostile

The Fibonacci convergent denominators at this scale are `1,2,3,5`. Their
combined coefficient skeleton has minimal period `30`. Adding the `q=6`
owners still leaves period `30`; deleting only the `q=4` owners from the full
packet also leaves period `30`. The non-Fibonacci denominator `4` supplies
the missing 2-adic clock and changes `lcm(30,4)` to `60`. Hence Fibonacci or
Zeckendorf data alone cannot explain the exact tail phase.

There is also a direct golden hostile. For

\[
 \alpha={\sqrt5-1\over2},
\]

the lower mechanical word is genuinely Fibonacci/Sturmian, but
`9 alpha<6<10 alpha`. Since its partial sums `floor(e alpha)` start at zero
and increase by zero or one, the values for `0<=e<=10` include every integer
from zero through six. This slope covers all seven residues already at
`m=11`; it is not a persistent tail owner.

## 3. The phase begins at second order

Define the phase moments

\[
 M_k(r)=\sum_{c=0}^{5}c^k A_c(r).                     \tag{7}
\]

Equation `(3)` says `M_0(r)=127/5`. Expanding the six denominators in `(2)`,
uniformly over the finite phase set, gives

\[
 D(n+1)={127\over35n}+{M_1(r)\over7n^2}+O(n^{-3}),
 \qquad r=n\bmod60.                                   \tag{8}
\]

**FINITE-EXACT:** `M_1` has minimal period `60`, takes exactly `25` distinct
values, and has exact range

\[
 {131\over10}\le M_1(r)\le {39\over2}.               \tag{9}
\]

Returning to `m=n+1`, `(8)` becomes

\[
 D(m)={127\over35m}
 +{127/35+M_1(m-1\bmod60)/7\over m^2}+O(m^{-3}).       \tag{10}
\]

Thus there is no phase-independent second-order limit. On the subsequence
`m-1=r mod 60`,

\[
 \lim m^2\left(D(m)-{127\over35m}\right)
 ={127\over35}+{M_1(r)\over7}.                         \tag{11}
\]

These subsequential limits take `25` distinct values, with exact range

\[
 {11\over2}\le {127\over35}+{M_1(r)\over7}
 \le {449\over70}.                                    \tag{12}
\]

This locates the sixty-clock information precisely: all phases share the
leading `1/m` mass, while the phase becomes visible at order `1/m^2`.

## 4. Clearing the six counters gives a quasipolynomial

Let

\[
 Q_6(n)=\prod_{c=0}^{5}(n-c),
 \qquad Y(n)=Q_6(n)D(n+1).                             \tag{13}
\]

Multiplying `(2)` by `Q_6(n)` gives

\[
 Y(n)={1\over7}\sum_{c=0}^{5}
 A_c(r)\prod_{\substack{0\le d\le5\\d\ne c}}(n-d),
 \qquad r=n\bmod60.                                   \tag{14}
\]

For each fixed phase, `(14)` is a polynomial in `n` of degree exactly five,
with leading coefficient `127/35`. Therefore `Y` is an exact quasipolynomial
for every `n>=12`. Its minimal phase period is exactly `60`, and all sixty
phase polynomials are distinct.

Writing those phase polynomials in ascending powers of `n`, the exact minimal
periods of their degree `0,1,2,3,4,5` coefficient sequences are

```text
60, 60, 60, 60, 60, 1,                                (15)
```

and the corresponding numbers of distinct coefficient values are

```text
12, 59, 60, 60, 25, 1.                                (16)
```

This does not contradict THM-4029: the original bounded sequence `a(m)` is
not eventually quasipolynomial. The fixed polynomial multiplier `Q_6` is
essential.

### Sharpness of the common desingularizer

For each phase `r`, regard the right side of `(2)` as the rational function

\[
 b_r(n)={1\over7}\sum_{c=0}^{5}{A_c(r)\over n-c}.      \tag{17}
\]

Whenever `A_c(r)>0`, this function has a genuine simple pole at `n=c`, with
residue `A_c(r)/7`: every term with a different shift is regular there, so
the pole cannot cancel. The support column in Section 2 proves that for every
`c in {0,...,5}` there is at least one phase with `A_c(r)>0`.

Suppose `P in Q[n]` has the property that `P(n)b_r(n)` is a polynomial for
every phase `r`. More generally, it is enough to assume agreement with a
polynomial at every sufficiently large integer in that phase: equality at
infinitely many points forces equality of the rational continuations. Choose
an active phase for shift `c`. Cancellation of its genuine pole forces
`P(c)=0`. This holds for all six distinct shifts, hence

\[
 \prod_{c=0}^{5}(n-c)\mid P(n).                        \tag{18}
\]

Since `Q_6` itself works by `(14)`, it is the unique monic common
desingularizer of least degree. Moreover,

\[
 Q_6(n)=n(n-1)\cdots(n-5)=6!{n\choose6}.              \tag{19}
\]

This is the precise falling-factorial/Pascal bridge to the Sun
`2-4-6-8` work: the six affine pole locations force exactly the degree-six
binomial-basis atom. The bridge preserves the falling-factorial divisor; it
does not identify the LRC deficit with a Sun representation count.

## 5. Explicit P-recursive law

On a fixed phase, `Y(n)` has degree five. Its sixth forward difference with
step 60 therefore vanishes. Substituting `(13)` gives the exact recurrence

\[
 \boxed{\sum_{j=0}^{6}(-1)^{6-j}{6\choose j}
 Q_6(n+60j)D(n+60j+1)=0}
 \qquad(n\ge12).                                      \tag{20}
\]

This is a homogeneous polynomial-coefficient recurrence of order `360`
with zero coefficients at the unused intermediate shifts. Hence the deficit
sequence is P-recursive over `Q`. The companion independently evaluates the
left side of `(20)` exactly for every `12<=n<500`; its vanishing for all
`n>=12` follows from the degree-five phase-polynomial proof, not from finite
extrapolation.

## 6. Ordinary generating function: D-finite but nonalgebraic

Let

\[
 \mathcal D(z)=\sum_{n\ge12}D(n+1)z^n.                \tag{21}
\]

For a phase `r`, let `n_r` be the least integer at least `12` congruent to
`r mod 60`, and put `t_{r,c}=n_r-c`. Since `t_{r,c}>=7`, define

\[
 L_t(z)=\sum_{k\ge0}{z^{t+60k}\over t+60k}.           \tag{22}
\]

Its derivative is rational:

\[
 L_t'(z)={z^{t-1}\over1-z^{60}}.                      \tag{23}
\]

Reindexing `(2)` phase by phase gives the exact finite expression

\[
 \mathcal D(z)={1\over7}\sum_{r=0}^{59}\sum_{c=0}^{5}
 A_c(r)z^c L_{t_{r,c}}(z).                             \tag{24}
\]

Every `L_t` is an antiderivative of a rational function, hence is D-finite;
finite sums and multiplication by polynomials preserve D-finiteness. Thus
`mathcal D` is D-finite. Over `C`, partial fractions in `(23)` also express
each `L_t` by rational functions and logarithms at sixtieth roots of unity.

The function is not algebraic. Uniformly in the phase, `(8)` gives

\[
 D(n+1)={127\over35n}+O(n^{-2}).                       \tag{25}
\]

Therefore, as real `z` increases to one,

\[
 \mathcal D(z)={127\over35}\log {1\over1-z}+O(1).    \tag{26}
\]

An algebraic function over `C(z)` has a Newton--Puiseux expansion at a finite
singularity; an unbounded branch has power-law, not logarithmic, growth.
Equation `(26)` is impossible for an algebraic function. Hence `mathcal D` is
D-finite but nonalgebraic, and in particular nonrational. This is consistent
with THM-4029's stronger rejection of a constant-coefficient recurrence.

## 7. Exact triangular track bridge

For an owner of denominator `q<=6` and a track `0<=s<q`, THM-4029 uses the
largest time in that residue class,

\[
 E_s(n)=n-((n-s)\bmod q).                              \tag{27}
\]

Let `N=qK+r`, where `K>=1` and `0<=r<q`. Then

\[
 \boxed{\sum_{n=q}^{N}E_s(n)
 =q^2 T_{K-1}+(r+1)(qK+s)-q\min(r+1,s)},               \tag{28}
\]

where `T_{K-1}=K(K-1)/2`.

Indeed, on the complete block `n=qt+u`, `0<=u<q`, the first `s` values have
`E_s(n)=q(t-1)+s` and the remaining `q-s` values have `E_s(n)=qt+s`. Their
sum is exactly `q^2 t`, independent of `s`. Summing the complete blocks
`1<=t<K` gives `q^2 T_{K-1}`, and the final partial block gives

```text
(r+1)(qK+s)-q min(r+1,s).
```

This is a lawful connection to triangular-number clocks, but it is not a
second derivation of the deficit. Block summation preserves total arrival
time and phase. It destroys the reciprocal costs `1/E_s`, the inner and
outer max-min selectors, the owner side, and the missing-sector label. Those
sidecars are exactly what `(2)` retains.

## 8. Strict word terminology and scope

For fixed `theta=j+alpha`, the binary part of

```text
floor((e+1)theta)-floor(e theta)
```

is a lower mechanical word of slope `alpha`.

- If `alpha` is irrational, the infinite word is Sturmian in the strict
  sense: aperiodic with factor complexity `ell+1`.
- If `alpha` is rational, the word is periodic, equivalently a rational
  mechanical word or a repetition of a Christoffel block. It is not
  Sturmian under the strict aperiodic definition.

The twelve persistent centres have `x=p/q`, `q<=6`. Under `theta=7x`, the
fractional mechanical slope has exact denominator `q`, because
`gcd(7,q)=1`; its word has exact period `q`. Irrational slopes occur in every
finite noncover neighborhood, but no irrational is a persistent owner.

Finally, `D(m)` is a scalar measure, not a mechanical word. Its coefficient
controller has sixty phases, but its values are not periodic because the
affine counters `n-c` grow without bound. Forgetting that counter would turn
the true holonomic law into the false finite-state or C-finite scalar model
already excluded by THM-4029.

This theorem concerns only the consecutive AP cover. It does not prove that
the AP maximizes seven-sector cover measure among sparse sets, does not
classify the larger-span residual in THM-536, and does not prove LRC(14).

## 9. Replay

From the repository root:

```text
python -B 04-computation/lrc14_ap_cover_holonomic_sixty_phase_thm4038.py
```

The companion imports the proved THM-4029 owner law and performs only exact
rational and integer checks. **QED.**
