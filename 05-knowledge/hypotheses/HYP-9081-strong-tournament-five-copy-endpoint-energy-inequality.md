---
id: HYP-9081
title: "Strong-tournament five-copy endpoint-energy inequality"
status: >
  OPEN HYPOTHESIS. Conjectures the five-copy endpoint
  dispersion inequality for every finite strong tournament, with equality
  only at C3. Conditional on it, every finite no-sink tournament satisfies
  Delta_V>=3(W+2H)^2/8, with equality exactly at a transitive prefix followed
  by C3. FINITE-EXACT support through strong order nine and exact sharp
  near-ordinal families are recorded in THM-4219; neither is an
  all-order proof.
source: codex-endpoint-gap-session-20260826
depends_on:
  - THM-4208-cycle-prefix-arbitrary-context-recurrence-endpoint-energy-and-eventual-positivity
related:
  - THM-4219-no-sink-endpoint-energy-floor-and-near-ordinal-sharpness
---

# HYP-9081 -- strong-tournament five-copy endpoint-energy inequality

**OPEN HYPOTHESIS.**

## 1. The open strong inequality

For a finite tournament `X`, retain THM-4208's

```text
H=H(X),                         W=W(X),
v_i=V_i^0+V_i^1,               m=sum_i v_i^2,
b=W+2H,
Delta_V=b^2-H^2-3m.                                  (1)
```

> **Hypothesis 1 (five-copy strong endpoint dispersion).** Every finite
> strong tournament satisfies
>
> ```text
> 5m<=(W+H)(W+3H)=b^2-H^2,                            (2)
> ```
>
> with equality exactly for `C3`.

The coefficient `5` is conjecturally best possible. THM-4219's exact
near-ordinal strong families have the ratio on the right of `(2)` divided by
`m` tending to `5` from above.

## 2. Exact one-defect form

Let `D_i` count permutations with exactly one backward adjacency, located
immediately after `i`, and put `D=sum_iD_i`. THM-4219's split bijection gives

```text
v_i=H+D_i,                     W=(n-1)H+D.             (3)
```

Thus Hypothesis 1 is exactly

```text
5sum_iD_i^2
 <=D^2+2(n-4)HD+n(n-3)H^2.                            (4)
```

This is the unresolved structural core: the unique-defect label must be
globally dispersed. Componentwise canonical-slack bounds do not suffice. In
the unique strong order-four tournament,

```text
H=5,
v=(7,5,6,9),
sigma=(13,7,4,7),                                    (5)
```

so both `sigma_i>=v_i` and `sigma_i>=H` fail pointwise.

An endpoint-pair matrix also does not close `(2)` by degree floors alone.
Indeed, let `N` be an abstract symmetric loopless matrix on
`{0,1,2,3,4}`, take endpoint counts `End_i=1`, put every off-diagonal
`N_ij=1`, and add mass `L` to `N_12`. Then

```text
H=degree(0)=4,
(v_1,v_2,v_3,v_4)=(L+4,L+4,4,4),                     (6)
```

while the genuine endpoint-matrix floors
`N_ij>=max(End_i,End_j)` and `v_i>=H` hold. Nevertheless
`sum_i v_i^2/(sum_i v_i)^2` tends to `1/2`, far above the five-copy scale
`1/5`. Any proof must use additional path-cover switching structure, not
only the symmetric degree model.

## 3. Conditional three-eighths theorem

Hypothesis 1 is a sufficient condition for, but is not equivalent to, the
following consequence.

> **Conditional theorem 2.** If Hypothesis 1 holds, then every finite
> no-sink tournament `X` satisfies
>
> ```text
> Delta_V(X)>=3(W(X)+2H(X))^2/8,                       (7)
> ```
>
> with equality exactly for
>
> ```text
> X=P_a triangleright C3,                    a>=0.     (8)
> ```

### Conditional proof

First let `S` be strong. Splitting every Hamilton path immediately after
each vertex gives `v_i>=H`, so if `n=|S|>=3`,

```text
b=sum_i v_i+H>=(n+1)H>=4H.                            (9)
```

Under Hypothesis 1,

```text
Delta_V=b^2-H^2-3m
 >=2(b^2-H^2)/5
 >=3b^2/8.                                             (10)
```

Equality in the last step forces `b=4H`, hence `n=3`; the only strong
three-vertex tournament is `C3`, and it attains equality.

Now let `X` have no sink. The strongly connected components of a tournament
are totally ordered, so write

```text
X=A triangleright S,                                  (11)
```

where `S` is the terminal strong component. It has at least three vertices,
since a singleton terminal component would be a sink. THM-4208's exact
ordinal law gives

```text
Delta_V(X)=H(S)^2Delta_V(A)+b_A^2Delta_V(S),
b_X=b_Ab_S.                                            (12)
```

The unconditional THM-4208 inequality has `Delta_V(A)>=0`, while `(10)`
applies to `S`; hence `(7)`. Equality in `(12)` forces
`Delta_V(A)=0`, which THM-4208 classifies as transitivity, and forces
`S=C3`. This is exactly `(8)`, including the empty prefix `P_0`. Conversely,
the same ordinal law verifies equality for every tournament in `(8)`.

## 4. The logical gap is nonzero

Put

```text
G=b^2-H^2-5m.                                         (13)
```

Direct algebra gives

```text
Delta_V-3b^2/8=(b^2-16H^2+24G)/40.                    (14)
```

Thus when `b>4H`, the three-eighths conclusion can hold even with `G<0`.
Equivalently, its allowed upper bound on `m` exceeds the bound in `(2)` by

```text
(b^2-16H^2)/120.                                       (15)
```

This is not merely formal. For every `a>=1`,

```text
Delta_V(P_a triangleright C3)=3b^2/8,
G(P_a triangleright C3)=-6(4^a-1)<0.                  (16)
```

Therefore Hypothesis 1 must remain scoped to **strong** tournaments, and
the words "sufficient" and "conditional" in section 3 cannot be replaced
by "equivalent".

## 5. Exact evidence, not proof

THM-4219 records an exact census of every strong tournament isomorphism
class through order nine. The class counts are

```text
1, 1, 6, 35, 353, 6008, 178133                       (17)
```

at orders `3,...,9`; there is no failure of `(2)`, and equality occurs only
at `C3`. The minimum gaps are

```text
0, 44, 510, 3186, 19842, 122778, 753810.              (18)
```

The exact `T(n,1)` and `T(n,2)` formulae in THM-4219 prove positivity on
those all-order families and show sharpness of the coefficient `5`. They do
not control arbitrary strong tournaments, so Hypothesis 1 and Conditional
theorem 2 remain open.
