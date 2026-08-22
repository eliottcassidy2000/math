---
id: THM-3032
title: "Sharpened half-tail critical-run extractor and the shell-four Pareto frontier"
status: >
  PROVED + VERIFIED-EXACT. The THM-2160 S5 half-tail extractor for AMM 12592
  admits a stratum-(ii) completion that is oblivious to the last tail bit,
  given in closed form by 2S(u) = sum_{i<m}(-u)^i + (1+u)^(m-2)
  - u(1-u^2)^(t-1), t = m/2. The resulting exactly fair rule has deadline
  T(1)=2 and, for m<=n<=2m-1 with m=2^r>=2: T(n)=3m/2 (n=m), 2m-1
  (m<n<=2m-2), 2m=n+1 (n=2m-1). Hence T(n)<=max(n+1,2n-2) for all n and
  T(n)<=max(n+1,2n-3) for n>=5; the floor T(n)=n+1 is attained at n=1 and at
  n=2^(r+1)-2, 2^(r+1)-1. The sharpening makes the half-tail rule dominate
  the THM-2225 checksum rule POINTWISE, strictly exactly at n=2^r (r>=2),
  where it saves 2^(r-1)-1 flips; without it the two are incomparable.
  One further ignorable coordinate is impossible (S(-1)=m/2!=0, and globally
  c<=1 by THM-2160 S6.3). Exhaustive over all 1036800 shell-balanced rules on
  the shell m=4: the achievable deadline vectors (T(4),T(5),T(6),T(7)) are
  exactly T(4),T(5) in {6,7,8}, T(6) in {7,8}, T(7)=8 with (T(4),T(5))!=(6,6);
  two Pareto-minimal points (6,7,7,8) [this rule] and (7,6,7,8). So T(4)>=6
  for every shell-balanced method, and T(4)=5 would require abandoning shell
  balance. Slope C=2 is unchanged: this is a D-improvement, not a C-improvement.
source: opus-2026-07-31-amm12592-writeup
depends_on:
  - THM-2160
  - THM-2225
related:
  - THM-2253  # dyadic contrast, n+2^nu_2(n); Pareto-incomparable profile
  - THM-2966  # spine normal form; C* = 1 + gamma*
  - HYP-9061  # the minimal linear deadline C*, OPEN
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
script: 04-computation/amm12592_refined_half_tail_referee.py
output: 05-knowledge/results/amm12592_refined_half_tail_referee.out
writeup: 06-writeups/amm12592-solution.tex
---

# THM-3032 -- the sharpened half-tail extractor

Bits are independent with `P(0)=p`, `P(1)=q=1-p`, `0<p<1` unknown; `n` is the
length of the maximal constant initial run. For a nonconstant stream let `m`
be the power of two with `m<=n<=2m-1`, and let

```text
S_m={w in {0,1}^(2m): w_1=...=w_m, w nonconstant}
```

be its shell. Put `t=m/2`, `b=X_1`, and `z_i=X_(m+i) xor b` for `1<=i<=m`;
membership in `S_m` says `z!=0`, and `z_1=1` iff `n=m`. On the branch `b=1`
give `z` the verdict opposite to its branch-`0` verdict. By THM-2225's
composition lemma, exact fairness reduces to the single layer condition

```text
#{z in H: |z|=j} = (1/2) binom(m,j),          1<=j<=m-1,        (L)
```

with `z=1^m` free (the middle pair `0^m 1^m / 1^m 0^m` self-cancels).

## 1. The rule

```text
(i)  z_1=1 (n=m):    heads iff z_2+...+z_t = 0 mod 2;  stop at flip 3m/2.
(ii) z_1=0 (n>m):    heads iff (z_2,...,z_(m-1)) in S;  stop at flip
                     2m-1 when m>=4 and n<=2m-2, else at flip 2m.
```

THM-2160 S5 supplies stratum (i) and the required completion counts

```text
R_j = (1/2)[binom(m-1,j) - e_j],   e_j = [u^j] u(1+u)(1-u^2)^(t-1),
e_(2a+1) = e_(2a+2) = (-1)^a binom(t-1,a),                        (R)
```

integral and in range by Lucas and Vandermonde.

## 2. The new step: the completion can ignore `z_m`

**Theorem.** For `m=2^r>=4` there is a set `S` in `{0,1}^(m-2)` realizing
(R) through `R_j = s_j + s_(j-1)`, `s_j = #{y in S: |y|=j}`.

*Proof.* The condition is `(1+u)S(u) = s_0 + R(u)` with
`R(u)=sum_(j=1)^(m-1) R_j u^j` (no requirement at `j=0`, because the partner
`z=0^m` of `y=0^(m-2)` is not a shell word). From (R) and `e_0=0`,

```text
2R(u) = [(1+u)^(m-1) - 1] - [u(1+u)(1-u^2)^(t-1) - e_m u^m],
e_m = (-1)^(t-1) = -1                     (t=m/2>=2 a power of two, so even).
```

Taking `s_0=1` and dividing by `1+u` (legal since `m` is even, using
`1-u^m=(1+u) sum_(i<m) (-u)^i`):

```text
2 S(u) = sum_(i=0)^(m-1) (-u)^i + (1+u)^(m-2) - u(1-u^2)^(t-1).    (S)
```

Hence the closed form

```text
s_j       = (1/2)[binom(m-2,j) + 1],                       j even,
s_(2a+1)  = (1/2)[binom(m-2,2a+1) - 1 - (-1)^a binom(t-1,a)].
```

`m-2=2^r-2` has binary expansion `1...10`, so by Lucas `binom(m-2,j)` is odd
for even `j<=m-2` and even for odd `j`: both are integers. The range
`0<=s_j<=binom(m-2,j)` reduces, for `j=2a+1`, to
`binom(m-2,2a+1) >= 1 + binom(t-1,a)` when `a` is even (which implies the
`a` odd requirement `>= binom(t-1,a)-1`). Since `t` is even, `a=t-1` is odd
and there `binom(m-2,m-1)=0=binom(t-1,t-1)-1`. For `a<=t-2`, writing
`m-2=(t-1)+(t-1)`, Vandermonde gives
`binom(m-2,2a+1) >= binom(t-1,a) binom(t-1,a+1) >= 2 binom(t-1,a)` unless
`a=t-2`, where `binom(m-2,m-3)=2t-2>=1+(t-1)`. QED

Note `s_(m-1)=0`, as it must be.

## 3. Deadline profile and consequences

```text
T(1)=2;   for m=2^r>=2 and m<=n<=2m-1:
T(n) = 3m/2      (n=m),
       2m-1      (m<n<=2m-2, so m>=4),
       2m = n+1  (n=2m-1).

T(n) <= max(n+1, 2n-2)   for all n>=1,
T(n) <= max(n+1, 2n-3)   for all n>=5,
T(n) = n+1               at n=1, 2^(r+1)-2, 2^(r+1)-1.
```

Enumerated: `2,3,4,6,7,7,8,12,15,15,15,15,15,15,16,24,31,...`

**Pointwise domination.** The exact profiles are `2m` (THM-2225 part (a)),
`2m-1` except `2m` at `n=2m-1` and `2` at `n=1` (THM-2225 checksum), and the
display above. So the checksum beats compression outside `{1} u {2m-1}`, and
this rule beats the checksum **exactly at `n=2^r`, `r>=2`**, saving
`2^(r-1)-1` flips. Without section 2 the half-tail rule costs `2m` on all of
stratum (ii) and is therefore *incomparable* to the checksum (worse at
`n=5,6,9,...`): the sharpening is what converts a trade into a domination.

**One ignorable coordinate is maximal.** THM-2160 S6.3 gives `c<=1` for
globally ignored tail coordinates. Inside this scheme the obstruction is
explicit: dropping `z_(m-1)` too needs `(1+u)^2 | s_0+R(u)`, i.e. `S(-1)=0`,
but (S) gives `S(-1)=m/2`.

## 4. Exhaustive shell-four frontier

Call a rule *shell-balanced* if it decides by flip `2m` and bisects every
composition class of every shell (sufficient for exact fairness; not
necessary, since shells could compensate one another). On `m=4` there are
exactly `1036800` shell-balanced rules. Enumerating all of them:

```text
achievable (T(4),T(5),T(6),T(7)):
    T(4),T(5) in {6,7,8},  T(6) in {7,8},  T(7)=8,
    subject to (T(4),T(5)) != (6,6).
Pareto-minimal:  (6,7,7,8)  [this rule]   and   (7,6,7,8).
THM-2253 dyadic contrast realizes (8,6,8,8): NOT minimal.
```

Consequences: `T(4)>=6`, `T(6)>=7`, `T(7)>=8` for every shell-balanced rule
(the last two are the floor `n+1`); the componentwise minimum `(6,6,7,8)` is
realized by no single rule; and any method with `T(4)=5` must abandon shell
balance. Combined with THM-2253 (floor at every odd `n`) and section 3 (floor
at `n=2^(r+1)-2`), the exact optimum `T(n)=n+1` is now settled for every odd
`n` and every `n=2^(r+1)-2`; the first undetermined value is
`T(4) in {5,6}`.

**Later supersession.** The final sentence records the frontier when this
theorem was proved. THM-3337 subsequently gave `T_opt(4)=5`, and THM-3340
proved `T_opt(n)=n+1` for every positive integer `n`. The shell-balanced
enumeration and its two Pareto points remain valid exactly in their stated
class.

## 5. Scope

This is a `D`-improvement, not a `C`-improvement: the slope stays `2` (worst
case `n=m+1`, `T=2m-1=2n-3`), so HYP-9061's `C*` question is untouched. The
shell-four enumeration is a lower bound only within the shell-balanced class.

## 6. Referee

```bash
python3 04-computation/amm12592_refined_half_tail_referee.py
```

Checks A-H: composition balance and deadlines for all three constructions
through the shell `2m=32`; exact fairness of the sharpened rule as a
polynomial identity at `p=1/3,2/7,9/11`; causality at every claimed stop;
the closed form of section 2 legal for `m=4..4096`; and the full shell-four
enumeration. All assertions pass. QED.
