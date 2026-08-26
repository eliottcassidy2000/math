---
id: THM-4193
title: "Cycle-first transitive-tail crossing and transitive-context positivity"
status: >
  PROVED exact all-order singleton crossing and strict transitive-context
  prefix positivity for C3 followed by a transitive tail of length at least
  five, with the sharp failure boundary at tail lengths zero through four +
  FINITE-EXACT unique source-free singleton survivor through order eight and
  positive arbitrary-context factor-class census through order seven +
  VERIFIED-EXACT + INDEPENDENTLY AUDITED. The universal arbitrary-context
  threshold for this cycle-prefix family is subsequently proved in THM-4212.
  Classification for general left factors, general (OS+), and the
  order-eleven asymmetric bank remain OPEN.
source: root-tournament-os-plus-extension-20260826
depends_on:
  - THM-4184-path-cover-parity-ordinal-cocycle-and-lollipop-positivity
  - THM-4187-universal-source-padding-supermodularity-and-transitive-left-ordinal-positivity
related:
  - THM-4181-ordinal-sum-capacity-transfer-and-parity-component-exchange
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
  - THM-4212-cycle-prefix-uniform-arbitrary-context-threshold-and-tail-five-lower-bound
script: 04-computation/tournament_cycle_first_transitive_tail_crossing_thm4193.py
output: 05-knowledge/results/tournament_cycle_first_transitive_tail_crossing_thm4193.out
independent_audit_script: 04-computation/tournament_cycle_first_transitive_tail_crossing_independent_audit_thm4193.py
independent_audit_output: 05-knowledge/results/tournament_cycle_first_transitive_tail_crossing_independent_audit_thm4193.out
script_sha256: 3264096090d6aab569f96120d064f6fa42b4e2fd150304438f929abdaff5ec5c
output_sha256: 82292a9eb8256286966a3f2ad9bbc9f1af90d69bfa877aa538c046d1373d99cc
independent_audit_script_sha256: bf83bf2bf17f2d2c35a502ae793636ec8b3ff13431e9100189aa0bf74cafcfb8
independent_audit_output_sha256: f43b9bfa1e5e8d53c4913caa5bb87429110204fca3534fcb24ef687103372e38
hash_basis: raw LF bytes
primary_audit: >
  PASS. The inherited exact ordinal-capacity engine checks 576 direct closed-
  formula rows, 167,936 algebra rows, 163,776 tail increments, the complete
  7,412-class singleton census through order eight, and 283,024 arbitrary-
  context presentations for the first survivor. Normal, optimized, and
  fixed-hash streams byte-match.
independent_audit: >
  ACCEPT. A standalone actual-child engine imports no tournament code and
  rebuilds capacities from complementary Hamilton-path tables. It checks
  twelve transitive gates, eight singleton rows, 120 direct context rows,
  and 167,936 arithmetic formula rows. Normal, optimized, and fixed-hash
  streams byte-match.
---

# THM-4193 -- cycle-first transitive-tail crossing and transitive-context positivity

**PROVED all-order sharp crossing in transitive contexts + FINITE-EXACT
minimality sidecars + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4187 proved nonnegative local interaction for every transitive left
factor and exhibited the cycle-first negative family

```text
Theta(C3,P_b,P_c)<0.
```

That hostile does not remain negative after every transitive suffix is
absorbed into the left factor.  There is a sharp and unexpectedly late
crossing.  Put

```text
A_n=C3 triangleright P_n,                n>=0,             (1)
```

where `P_n` is the transitive tournament of order `n` and `P_0` denotes an
empty suffix only in `(1)`.  Every `A_n` is nontransitive and source-free.
For a fixed left prefix `A`, define its literal contextual defect

```text
F_A(B,C)=R_+(A triangleright B,C)-R_+(B,C).               (2)
```

Then `A_0,...,A_4` fail already at `(B,C)=(1,1)`, while `A_n` is strictly
positive against **every pair of transitive contexts** once `n>=5`:

```text
F_(A_n)(P_b,P_c)>0                 for all b,c>=1
                                   if and only if n>=5.   (3)
```

Thus `A_5=C3 triangleright P_5`, of order eight, is the first proved
nontransitive left prefix beyond THM-4187's transitive cone in this exact
context family.  The mechanism is a competition between a positive
`4^(n+b+c)` transitive-gate term and the cycle-first negative `2^(n+b+c)`
curvature.  Five tail vertices are exactly where the former wins uniformly.

Equation `(3)` does **not**, inside this theorem, handle arbitrary `B,C` or
classify all `A` satisfying `(2)`. THM-4212 subsequently proves the sharp
arbitrary-context threshold for this `A_n` family. Classification for general
left factors and general `(OS+)` remain open.

## 1. Conventions and inherited algebra

All quantified tournaments are finite and nonempty.  Use THM-4184/4187's
Hamilton count `H`, exposure capacity `c`, gate `G_+`, ordinal sum
`triangleright`, remainder `R_+`, and local interaction `Theta`.  Write

```text
T_m=G_+(c(P_m)).                                           (4)
```

THM-4184 proves

```text
T_1=0,
T_m=[4*4^m-9(m+2)2^m+24m+32]/18,          m>=2,           (5)
```

and therefore

```text
r_(a,c):=R_+(P_a,P_c)=T_(a+c)-T_a-T_c.                    (6)
```

THM-4187 proves, for `a,c>=1`,

```text
theta_(a,c):=Theta(C3,P_a,P_c)
 =144-(27c+153)2^c,                                      a=1,
 =-2^(a-1)[(27(a+c)+126)2^c-(27a+126)],                  a>=2. (7)
```

Every value in `(7)` is negative.  The point of this theorem is that `(7)`
is only one term in the contextual defect after the transitive suffix has
joined the left prefix.

## 2. Exact singleton crossing

> **Theorem 1 (sharp singleton formula).** For every `n>=0`,
>
> ```text
> F_(A_n)(1,1)
>  =12[2*4^n-(3n+21)2^n+1].                              (8)
> ```
>
> It is negative exactly for `0<=n<=4` and positive exactly for `n>=5`.
> The boundary rows are
>
> ```text
> n:                    0     1     2      3      4       5
> F_(A_n)(1,1):      -216  -468  -900  -1332   -180   10764. (9)
> ```

### Proof

Since `H(C3)=3`, the definition of `Theta` gives

```text
R_+(C3 triangleright P_m,P_c)
 =9r_(m,c)+theta_(m,c).                                  (10)
```

Set `m=n+1` and `c=1` in `(10)`.  Since `R_+(1,1)=0`, equations
`(2)` and `(5)--(7)` simplify exactly to `(8)`, including the `m=1`
endpoint.  Direct substitution gives the five negative rows and the first
positive row in `(9)`.

For `n>=5`,

```text
2^(n+1)>3n+21                                             (11)
```

holds at `n=5` (`64>36`) and persists because the left side doubles while
the right side increases by three.  Thus the bracket in `(8)` is positive.
This proves the sign classification. QED.

The negative value at `n=4` is a canonical hostile to shortening the tail.
The crossing is not a source-padding consequence: every `A_n` is
source-free, since no cycle vertex dominates the other two and every tail
vertex loses to all three cycle vertices.

## 3. Exact transitive-context formula

> **Theorem 2 (cycle-first buffered transitive contexts).** For all
> `n>=0` and `b,c>=1`,
>
> ```text
> F_n(b,c):=F_(A_n)(P_b,P_c)
>  =9r_(n+b,c)-r_(b,c)+theta_(n+b,c).                     (12)
> ```
>
> Moreover,
>
> ```text
> F_n(b,c)>0 for every b,c>=1    iff n>=5.                (13)
> ```

### Proof of the formula

Associativity gives

```text
A_n triangleright P_b=C3 triangleright P_(n+b).
```

Apply `(10)` with `m=n+b` and subtract `R_+(P_b,P_c)=r_(b,c)`.
This is exactly `(12)`.  Notice that the factor `9=H(C3)^2` is
load-bearing; discarding it would compare different normalizations. QED.

### The base tail `n=5`

It remains to prove positivity uniformly in the two unbounded context
orders.  Define

```text
Q_j=[4*4^j-9(j+2)2^j+24j+32]/18,          j>=1,
epsilon_j=1 if j=1, and 0 otherwise.                      (14)
```

Thus `Q_1=1` and `T_j=Q_j-epsilon_j`.  Expanding `(12)` at `n=5` gives

```text
18F_5(b,c)=N(b,c)+144epsilon_c-18epsilon_b,               (15)
```

where

```text
N(b,c)=
 36860*4^(b+c)-36860*4^b-32*4^c
 -10359b*2^(b+c)-10359c*2^(b+c)-93294*2^(b+c)
 +10359b*2^b+93294*2^b+72c*2^c+144*2^c-256.              (16)
```

For `b=c=1`, equation `(15)` gives `F_5(1,1)=10764`.  Now put
`N_0=b+c>=3`.  Dropping the positive terms on the last line of `(16)` and
using `c>=1`, `b>=1` gives

```text
18F_5(b,c)
 >=27637*4^N_0-(10359N_0+93294)2^N_0-274.                (17)
```

Indeed

```text
36860(4^N_0-4^b)>=27645*4^N_0,
32*4^c<=8*4^N_0,
144epsilon_c-18epsilon_b>=-18.                           (18)
```

The right side of `(17)` is `773526` at `N_0=3`.  Its forward difference is

```text
2^N_0[82911*2^N_0-(10359N_0+114012)]>0,     N_0>=3.      (19)
```

The inequality in `(19)` holds at three and then strengthens because its
exponential term doubles while the competing linear term grows by `10359`.
Thus `F_5(b,c)>0` for every `b,c>=1`.

### Every longer tail

Subtract `(12)` at consecutive tail lengths.  Exact cancellation yields,
with `s=n+b`,

```text
F_(n+1)(b,c)-F_n(b,c)
 =6*2^s B(s,c),                                          (20)

B(s,c)=2^s(4^c-1)-3[2^c(s+c+6)-(s+6)].                  (21)
```

For `n>=5`, one has `s>=6`.  At fixed `c`,

```text
B(s+1,c)-B(s,c)
 =2^s(4^c-1)-3(2^c-1)>0.                                (22)
```

Here `(4^c-1)=(2^c-1)(2^c+1)`, so positivity is immediate for
`s>=6,c>=1`.

Hence it suffices to use `s=6`, where

```text
B(6,c)=2^c(64*2^c-3c-36)-28.                            (23)
```

This is `150` at `c=1`; its forward difference is

```text
2^c(192*2^c-3c-42)>0.                                   (24)
```

Therefore `(20)` is strictly positive for every `n>=5,b,c>=1`.  The base
case proves the positive direction of `(13)`.  For `0<=n<=4`, Theorem 1's
choice `b=c=1` is negative, proving the converse. QED.

## 4. Exact minimality census through order eight

**FINITE-EXACT.**  The primary certificate uses one `gentourng`
representative of every tournament class through order eight.  It records
whether `A` has a universal source and the exact singleton defect
`F_A(1,1)`:

| `|A|` | classes | source, positive | source, nonpositive | source-free, negative | source-free, nonnegative |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 0 | 1 | 0 | 0 |
| 2 | 1 | 1 | 0 | 0 | 0 |
| 3 | 2 | 1 | 0 | 1 | 0 |
| 4 | 4 | 2 | 0 | 2 | 0 |
| 5 | 12 | 4 | 0 | 8 | 0 |
| 6 | 56 | 12 | 0 | 44 | 0 |
| 7 | 456 | 56 | 0 | 400 | 0 |
| 8 | 6,880 | 456 | 0 | 6,423 | 1 |

Thus the unique source-free class of order at most eight with nonnegative
singleton defect is

```text
A_5=C3 triangleright P_5,
gentourng pair-bit label 1011111111111111111111111111,
F_(A_5)(1,1)=10764.                                     (25)
```

There are no zero source-free rows.  The ordered singleton-census digest is

```text
2bd16f607e89f4bb8efb09c319497da8a42f7d553b83898f44ccf7dedbbcefed. (26)
```

This proves minimality only for the singleton test and only through order
eight.  It does not classify the universal contextual property `(2)`.

## 5. Arbitrary-context stress sidecar

**FINITE-EXACT, NON-LOAD-BEARING.**  For the fixed survivor `A_5`, take one
tournament-class representative of every `B,C` of orders one through seven.
There are `532` factor classes and

```text
532^2=283024                                              (27)
```

ordered presentations.  Every row satisfies

```text
F_(A_5)(B,C)>0,                                          (28)
```

with exact minimum `10764` at `B=C=1`.  The stronger finite inequality

```text
F_(A_5)(B,C)>=10764 H(B)^2H(C)^2                         (29)
```

also holds throughout this bank, again with equality only at `(1,1)`.

Equations `(28)--(29)` are evidence within this theorem, not its all-order
proof. THM-4212 subsequently proves `(29)` for all tournaments `B,C` and
propagates it to every `n>=5`. Neither result implies general `(OS+)`.

## 6. Audits and type firewall

The primary audit checks the exact ordinal capacity construction, `(8)`,
`(12)`, `(15)`, and `(20)` on independent finite heads before running the
complete singleton and arbitrary-context censuses.  Its class rows are
ordered factor presentations, not child-isomorphism classes.

The independent audit imports no tournament module, no `gentourng` data, and
no ordinal capacity-transfer implementation.  It rebuilds the capacity of
each actual child directly from complementary Hamilton-path tables and checks
the crossing and context formulas on all declared direct rows.  Normal,
optimized, and fixed-hash executions of both scripts byte-match their frozen
outputs.

The connection contract is

```text
source:       actual ordered factors C3,P_n,P_b,P_c,
map:          associate C3 triangleright P_(n+b) before evaluating R_+,
preserved:    H, the actual capacity tensor, G_+, R_+, and factor order,
mechanism:    exact competition of 4-power gate growth with 2-power
              cycle-first curvature,
lost by order-only transitive contexts:
              arbitrary B/C rooted parity states and capacity geometry,
needed sidecar:
              THM-4184's exact T_m and THM-4187's theta_(a,c),
hostile:      (C3 triangleright P_4,1,1), defect -180,
first survivor:
              (C3 triangleright P_5,1,1), defect 10764.
```

Nothing proved inside THM-4193 alone establishes arbitrary-context prefix
positivity for `A_5`; THM-4212 subsequently does so and proves the sharp
threshold for all `A_n`. A classification of general left factors satisfying
`(2)`, arbitrary-left local positivity, general `(OS+)`, THM-4177's
no-sink/no-source gate law, the order-eleven asymmetric bank, exact Johnson
cosets, and actual response maximizers remain open.

## 7. Replay

Primary exact audit:

```bash
python3 -B \
  04-computation/tournament_cycle_first_transitive_tail_crossing_thm4193.py
python3 -O -B \
  04-computation/tournament_cycle_first_transitive_tail_crossing_thm4193.py
PYTHONHASHSEED=4193 python3 -B \
  04-computation/tournament_cycle_first_transitive_tail_crossing_thm4193.py
```

Independent direct-capacity audit:

```bash
python3 -B \
  04-computation/tournament_cycle_first_transitive_tail_crossing_independent_audit_thm4193.py
python3 -O -B \
  04-computation/tournament_cycle_first_transitive_tail_crossing_independent_audit_thm4193.py
PYTHONHASHSEED=4193 python3 -B \
  04-computation/tournament_cycle_first_transitive_tail_crossing_independent_audit_thm4193.py
```

**QED here for the exact singleton crossing, sharp all-order
transitive-context positivity threshold, and the stated finite-exact
sidecars. THM-4212 is the subsequent arbitrary-context closure.**
