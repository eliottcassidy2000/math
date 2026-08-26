---
id: THM-4216
title: "Tail-five non-singleton right-factor positivity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For
  Q_n=C3 triangleright P_n, bare ordinal remainder is positive against every
  non-singleton right factor exactly when n>=5 and against every nonempty
  right factor exactly when n>=6. At tail five the exact four-debt
  decomposition gives the strict polynomial floor and the uniform
  non-singleton bound 2730H(C)^2. Consequently THM-4213's separated
  multi-cycle language needs final tail five, not six. General (OS+),
  adjacent-cycle thresholds, and arbitrary-left bare positivity remain open.
source: root-frontier-synthesis-20260826
depends_on:
  - THM-4187-universal-source-padding-supermodularity-and-transitive-left-ordinal-positivity
  - THM-4208-cycle-prefix-arbitrary-context-recurrence-endpoint-energy-and-eventual-positivity
  - THM-4212-cycle-prefix-uniform-arbitrary-context-threshold-and-tail-five-lower-bound
  - THM-4213-uniform-prefix-ordinal-semigroup-and-tail-five-cycle-language
related:
  - THM-4193-cycle-first-transitive-tail-crossing-and-transitive-context-positivity
script: 04-computation/tournament_tail_five_nonsingleton_right_factor_positivity_thm4216.py
output: 05-knowledge/results/tournament_tail_five_nonsingleton_right_factor_positivity_thm4216.out
independent_audit_script: 04-computation/tournament_tail_five_nonsingleton_right_factor_positivity_independent_audit_thm4216.cpp
independent_audit_output: 05-knowledge/results/tournament_tail_five_nonsingleton_right_factor_positivity_independent_audit_thm4216.out
script_sha256: 4e189e468b580fa5bacc8ff80fceabe4ecb72002196fc989a447e40012a7c5e6
output_sha256: 35221cfa04c7fd9f9e89bd0f1ffa3febc188c53e4dd9f76603e3bbb729640921
independent_audit_script_sha256: b2b85eafe4bc6a823913a0098fe02c7dd68b109de481f110a44d5f3c4b9d4e52
independent_audit_output_sha256: 7a11ecba317d3367cf29acfcca5d1597b06e6f3442bc4228d96be2b932e64ac4
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact symbolic coordinate change reconstructs the tail-five unary
  response, turns the incoming-current term into a nonnegative edge-avoidance
  term, and verifies the four nonnegative debts in the advertised floor.
  The inherited exact transfer engine checks the identity and bound on all
  1,099 labelled right factors through order five, both sharp hostile lists,
  2,198 later-tail rows, and 34 final-tail-five no-sink controls. Normal and
  python -O streams byte-match.
independent_audit: >
  ACCEPT. A standalone C++17 referee reconstructs tournaments from literal
  adjacency, uses its own subset path DP, and imports neither ordinal transfer
  nor response jets. It independently recovers the Q5 jet, unary coefficients,
  fan rewrite, exact four-debt identity, strict tau boundary, all-adjacencies
  marked-Hamilton bound, hostile tables, later-tail propagation, and 34
  final-tail-five no-sink controls. O0, O3, and ASan/UBSan streams byte-match.
---

# THM-4216 -- tail-five non-singleton right-factor positivity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The result sharpens THM-4213's final-tail-six `(OS+)` corollary to
final tail five. The mechanism is an exact positive decomposition hidden in
THM-4208's unary response: after the incoming fan is traded for Hamilton
starts, the signed current becomes an edge-avoidance mass.

## 1. Statement and inheritance pass

All tournaments quantified below are finite and nonempty. Use the Hamilton
count `H`, capacity mass `W`, ordinal remainder `R_+`, and ordinal sum
`triangleright` from THM-4181/4184. Put

```text
Q_n=C3 triangleright P_n,                    n>=0,       (1)
```

where `P_0` is empty notation only in `(1)`.

> **Theorem 1 (exact tail-five positive decomposition).** For
> every right factor `C`, write `h=H(C)` and `w=W(C)`. Then
>
> ```text
> R_+(Q_5,C)
>  >=6[-33h^2+274hw+214w^2].                           (2)
> ```
>
> If `C` is non-singleton, then `w>=h` and consequently
>
> ```text
> R_+(Q_5,C)>=2730h^2>0.                               (3)
> ```
>
> The excluded singleton has the exact value
>
> ```text
> R_+(Q_5,P_1)=-180.                                   (4)
> ```

> **Theorem 2 (two sharp right-factor thresholds).** One has
>
> ```text
> R_+(Q_n,C)>0 for every non-singleton C   iff n>=5,
> R_+(Q_n,C)>0 for every nonempty C        iff n>=6.    (5)
> ```

> **Corollary 3 (final-tail-five multi-cycle `(OS+)`).** Every
> THM-4213 strong-component word
>
> ```text
> P_(a_0) triangleright C3 triangleright P_(b_1)
>  triangleright ... triangleright C3 triangleright P_(b_m),
> a_0>=0, m>=1, b_1,...,b_m>=5,                        (6)
> ```
>
> satisfies `R_+(X,C)>0` for every no-sink right factor `C`.

The closest proved mechanism is THM-4208's four-mode unary response and its
endpoint-energy inequality. The canonical hostile is `(Q_5,P_1)`: a
uniformly positive contextual prefix can still have negative bare remainder
because the empty middle factor is not in THM-4212's quantifiers. The
corrected near miss is therefore to infer bare `(OS+)` directly from a
positive `F_A(B,C)`. The least-used sidecar is the edge-avoidance mass
`W-d_i`, retained only after replacing the incoming fan by endpoint mass
minus Hamilton starts.

The live concept board was

```text
unary response | endpoint energy | incoming fan | Hamilton starts
edge avoidance | contextual floor | ordinal ideal.                    (7)
```

## 2. The exact incoming-fan rewrite

For `C`, let `d_i,r_i` be the unsigned and incoming capacity masses. Retain
THM-4208's right rooted states and put

```text
v_i=V_i^0+V_i^1,
t_i=V_i^1-V_i^0=Start_i(C),
m=sum_i v_i^2,       p=sum_i v_i t_i,
tau=sum_i t_i^2.                                      (8)
```

All `v_i,t_i` are nonnegative, and THM-4208's exact fan/endpoint identities
give

```text
v_i=r_i+t_i,
sum_i v_i=w+h,               sum_i t_i=h.              (9)
```

Write

```text
E=sum_i(1143v_i+9t_i)(w-d_i)>=0.                       (10)
```

The inequality is literal: `w-d_i` is the mass of all capacity edges
avoiding `i`.

THM-4208's unary formula at tail five is

```text
R_+(Q_5,C)
 =353088h^2+472302hw+118656w^2
  +1134L_0+1152L_1
  -341856M_00-695052M_01-353268M_11,                  (11)
```

where

```text
L_e=sum_i V_i^e(w-d_i-4r_i),
M_ef=sum_i V_i^eV_i^f.                                 (12)
```

Changing from `(V^0,V^1)` to `(v,t)` gives exactly

```text
1134V_i^0+1152V_i^1=1143v_i+9t_i,                     (13)

-341856M_00-695052M_01-353268M_11
 =-347544m-5706p-18tau.                                (14)
```

Now use `r_i=v_i-t_i` inside `(12)`. The current term and `(14)` combine,
without an inequality, to give

```text
R_+(Q_5,C)
 =353088h^2+472302hw+118656w^2
  +E-352116m-1170p+18tau.                              (15)
```

This is the load-bearing step. Bounding `r_i` by `w` before making the
replacement in `(9)` destroys the sign.

## 3. Four nonnegative debts

THM-4208's endpoint energy is

```text
Delta_V=(w+h)(w+3h)-3m>=0.                             (16)
```

Nonnegativity and the totals in `(9)` also give

```text
p<= (sum_i v_i)(sum_i t_i)=(w+h)h.                    (17)
```

Subtraction of the floor in `(2)` from the exact identity `(15)`
produces the exact certificate

```text
R_+(Q_5,C)-6[-33h^2+274hw+214w^2]
 =E+117372Delta_V
   +1170[(w+h)h-p]+18tau.                              (18)
```

Every term on the right is nonnegative, proving `(2)`.
In fact `tau>0`, because every nonempty tournament has a Hamilton path and
the Hamilton-start counts sum to `h>0`.  Thus `(2)` is always strict, although
the weaker displayed form is the useful compositional certificate.

If `C` has at least two vertices, mark each of the `|C|-1` adjacencies in
every Hamilton path:

```text
w>=(|C|-1)h>=h.                                        (19)
```

Writing `w=h+u`, `u>=0`, the floor itself satisfies

```text
6[-33h^2+274hw+214w^2]-2730h^2
 =12u(351h+107u)>=0.                                   (20)
```

This proves `(3)`. Direct evaluation of the singleton rooted coordinates in
`(11)` gives `(4)`.

## 4. Sharp thresholds

For the failure direction in the first line of `(5)`, the fixed
non-singleton right factor `P_2` gives

```text
n:                         0     1      2      3      4
R_+(Q_n,P_2):           -288  -684  -1368  -2232  -1512. (21)
```

Theorem 1 proves the `n=5` entry. For `n=5+r`, `r>=1`, write

```text
Q_n=Q_5 triangleright P_r.                             (22)
```

By the definition of THM-4212's contextual defect,

```text
R_+(Q_n,C)=F_5(P_r,C)+R_+(P_r,C).                     (23)
```

THM-4212 gives `F_5(P_r,C)>=10764h^2`, while THM-4187 gives
`R_+(P_r,C)>=0`. Thus

```text
R_+(Q_n,C)>=10764h^2>0,                    n>=6,       (24)
```

even for singleton `C`.

For sharpness of the second line in `(5)`, the singleton rows are

```text
n:                         0     1     2     3      4     5
R_+(Q_n,P_1):            -72  -216  -468  -900  -1332  -180. (25)
```

Equations `(3)`, `(21)`, and `(23)--(25)` prove both classifications in
`(5)`.

## 5. Removing the final extra source

Let `X` be a word `(6)`, and split off its final cycle-tail block:

```text
X=A triangleright Q_(b_m).                              (26)
```

The prefix `A` may be absent. If it is present, THM-4213's semigroup/ideal
theorem gives `A in N`, so

```text
F_A(Q_(b_m),C)>=0.                                     (27)
```

A no-sink tournament is non-singleton. Theorem 1 handles `b_m=5`, and
`(24)` handles `b_m>=6`. Therefore

```text
R_+(X,C)
 =F_A(Q_(b_m),C)+R_+(Q_(b_m),C)>0,                    (28)
```

with the first term omitted when `A` is absent. This proves Corollary 3 and
removes the additional final singleton required by THM-4213 Corollary 4.

## 6. Scope firewall

The theorem proves a bare-remainder sign only for the cycle-first
transitive-tail family and for THM-4213's tail-five-separated language. It
does not prove that every member of the uniformly positive ideal has positive
bare remainder; the singleton hostile `(4)` explicitly refutes that
shortcut. It does not classify the neutral semigroup, settle adjacent-cycle
thresholds, prove the all-order no-sink/no-source gate law, or prove general
`(OS+)`.

The connection contract is

```text
source:       THM-4208's exact Q5 unary response,
map:          (V0,V1,r) -> (v,t,r=v-t),
preserved:    exact R_+, endpoint energy, Hamilton starts, edge incidence,
destroyed:    individual rooted-path addresses after scalar summation,
sidecar:      E=sum_i(1143v_i+9t_i)(W-d_i),
decisive tests:
              P1 at tail five and P2 through tails zero to four.          (29)
```

## 7. Exact replay

```bash
PYTHONPATH=04-computation python3 -B \
  04-computation/tournament_tail_five_nonsingleton_right_factor_positivity_thm4216.py

PYTHONPATH=04-computation python3 -O -B \
  04-computation/tournament_tail_five_nonsingleton_right_factor_positivity_thm4216.py
```

Both streams byte-match the frozen primary output. Independent audit remains
available through the standalone literal engine:

```bash
clang++ -std=c++17 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_tail_five_nonsingleton_right_factor_positivity_independent_audit_thm4216.cpp \
  -o /tmp/thm4216-independent
/tmp/thm4216-independent | diff -u \
  05-knowledge/results/tournament_tail_five_nonsingleton_right_factor_positivity_independent_audit_thm4216.out -
```

The independent O0, O3, and ASan/UBSan streams byte-match the frozen output.
