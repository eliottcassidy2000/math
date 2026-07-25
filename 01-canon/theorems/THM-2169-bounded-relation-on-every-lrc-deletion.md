---
id: THM-2169
title: "A bounded relation on every deletion of a zero-safe LRC(14) row"
status: >
  PROVED from THM-2144 and THM-2164. For every coordinate i of a distinct
  positive thirteen-speed row with zero-measure safe set, the other twelve
  speeds carry a nonzero integer relation of height at most 1247. The proof
  contracts one primitive height-29 relation against deletion-adapted
  anisotropic relations; six exact Kraft-admissible cap profiles cover all
  possible heights. Thus every deletion core relevant to a hypothetical
  LRC(14) counterexample lies on an explicit bounded relation hyperplane.
  After primitive normalization its l1 norm is at most 14963, with at most
  14962
  positive-row carry values; quotient-owner masks remain necessary. This is
  a structural reduction, not a finite speed bound, Q108, or LRC(14).
source: opus-2026-07-24-puzzle-atlas
depends_on:
  - THM-2144-anisotropic-selberg-kraft-relation-box
  - THM-2164-relative-packet-rank-harvesting
related:
  - THM-2162
  - THM-2163
  - THM-2166
  - THM-2167
script: 04-computation/lrc14_deletion_anisotropic_profiles_referee_codex_20260724.py
output: 05-knowledge/results/lrc14_deletion_anisotropic_profiles_referee_codex_20260724.out
script_sha256: e486ab21c1b7596de319998ebd8c353882ea654ece3e9eab73adcc9309c4ec56
output_sha256: 4034186b1a45b5ef2bdc8e63d8cc7b0b6d32fe22e39f465bb228de63cc693d12
hash_basis: working-tree bytes (LF)
---

# THM-2169 -- bounded relations on every deletion

Put

```text
J=[1/14,13/14],
S(v)={t in T:v_i t in J for every i},

Lambda(v)={m in Z^13:m.v=0}.                         (1)
```

Let `v=(v_1,...,v_13)` have pairwise distinct positive integer coordinates,
and suppose

```text
mu(S(v))=0.                                           (2)
```

## 1. The coordinatewise dichotomy

Fix a coordinate `i`. THM-2144, applied with cap `1` at coordinate `i` and
cap `105` elsewhere, supplies a nonzero relation

```text
r in Lambda(v),
|r_i|<=1,              ||r||_infinity<=105.           (3)
```

Divide by the gcd of its coefficients and call the resulting primitive
relation `w`. The bounds in (3) persist.

If `w_i=0`, then `w` itself is already a relation on the deletion
`v\{v_i}`, of height at most `105`.

Suppose `w_i!=0`; then `w_i=+1` or `-1`.

### Support at least three

If `|supp(w)|>=3`, THM-2164's whole-subtorus packet lemma gives

```text
s in Lambda(v)\Qw,              ||s||_infinity<=43.   (4)
```

Eliminate coordinate `i`:

```text
u=s_i w-w_i s.                                        (5)
```

Then

```text
u in Lambda(v),        u_i=0,
```

and `u!=0`, because otherwise (5) would make `s` rationally proportional to
`w`. Coordinatewise,

```text
|u_k|
 <=|s_i||w_k|+|w_i||s_k|
 <=43*105+43
 =4558.                                               (6)
```

Thus the deletion has a nonzero height-`4558` relation.

### Support two

If `|supp(w)|=2`, write its support as `{i,j}`. Positivity of the speeds
forces the two coefficients to have opposite signs. Since `|w_i|=1`,

```text
v_i=q v_j,                 q=|w_j|.                   (7)
```

Distinctness excludes `q=1`, while (3) gives

```text
2<=q<=105.                                            (8)
```

We have proved the sharp coordinatewise fork:

```text
either
  Lambda(v) contains nonzero u with u_i=0,
  ||u||_infinity<=4558,

or
  v_i=q v_j for some j!=i and 2<=q<=105.              (9)
```

The lock in the second branch is directed: it says the selected speed
`v_i`, not merely one of the pair, is the bounded multiple.

## 2. Deletion-adapted anisotropy sharpens every deletion

The coordinatewise fork is useful, but its cap profile is not adapted to the
relation that it is trying to eliminate. THM-2144 supplies one nonzero
height-29 relation. Divide by its coefficient gcd and call it

```text
g in Lambda(v),             ||g||_infinity<=29.       (10)
```

Write

```text
M=||g||_infinity.                                     (11)
```

Choose `p` with `|g_p|=M`. If `M=1`, then `g` is signed-unit. Positivity and
distinctness rule out support two, so THM-2164 supplies an independent
relation `s` of height at most `34`. For a coordinate with `g_i=0`, use `g`
itself. Otherwise

```text
u^(i)=s_i g-g_i s.                                    (12)
```

It is nonzero, deletes `i`, and has height at most `68`.

Now suppose `2<=M<=29`. We repeatedly apply THM-2144 with a cap at `p`
strictly below `M`. After primitive normalization the resulting relation
`t` is independent of `g`: proportional primitive integer vectors differ
only by sign, whereas `|t_p|<|g_p|`.

### Deleting the maximal coordinate

For `i=p`, use one of the following two profiles:

| height of `g` | `H_p` | the other twelve caps | Kraft sum |
|---|---:|---:|---:|
| `2<=M<=7` | `1` | `105` | `12194/12217` |
| `8<=M<=29` | `7` | `36` | `12446/12595` |

Both displayed fractions are strictly less than one. Contracting

```text
u^(p)=t_p g-g_p t                                     (13)
```

gives respectively

```text
||u^(p)||_infinity<=106M<=742,
||u^(p)||_infinity<=43M<=1247.                        (14)
```

### Deleting any other supported coordinate

Fix `i!=p`. If `g_i=0`, then `g` itself is a height-29 relation on that
deletion. If `g_i!=0`, use the appropriate row of this table:

| height of `g` | `H_p` | `H_i` | each other cap | Kraft sum |
|---|---:|---:|---:|---:|
| `2<=M<=7` | `1` | `35` | `126` | `3257870/3258253` |
| `8<=M<=22` | `7` | `9` | `46` | `1061102/1064965` |
| `23<=M<=27` | `22` | `8` | `36` | `2022566/2025505` |
| `28<=M<=29` | `27` | `9` | `34` | `51854/51925` |

Again every Kraft sum is strictly below one and every `p`-cap is below
`M`. The contraction

```text
u^(i)=t_i g-g_i t                                     (15)
```

is nonzero, since otherwise `g_i!=0` would make `t` proportional to `g`.
For every `k`,

```text
|u^(i)_k|<=H_i M+|g_i|H_k.                            (16)
```

Taking the largest non-`i` cap in each row gives, in order,

```text
||u^(i)||_infinity
 <=161M<=1127,
 <= 55M<=1210,
 <= 44M<=1188,
 <= 43M<=1247.                                       (17)
```

The six profiles are not a numerical search assumption: their rational
Kraft sums are displayed exact, and they cover the finite integer range
`2<=M<=29`. Combining them with the signed-unit branch proves:

> **Deletion theorem.** For every `i=1,...,13`, there is a nonzero
> `u^(i) in Lambda(v)` such that
>
> ```text
> u^(i)_i=0,             ||u^(i)||_infinity<=1247.    (18)
> ```

Equivalently, every twelve-speed deletion of `v` lies on an integer relation
hyperplane of coefficient height at most `1247`.

## 3. Exact interfaces

### Q108 / deletion induction

Any zero-safe row relevant to LRC(14) therefore has **all thirteen** of its
twelve-speed deletion cores in the bounded-relation class (18). A
deletion-based proof of LRC(14) does not need a Q108-type stability theorem
for arbitrary twelve-sets first; it is enough to prove the required
measure/stability statement on twelve-sets carrying a height-`1247`
relation, together with the directed locks (7).

This is a genuine restriction but not a finite search. A bounded relation
still leaves eleven rational degrees of freedom and unbounded speeds.

### Radix carry

Regard `u^(i)` as a twelve-coordinate relation after deleting its zero
coordinate, and divide by the coefficient gcd. Primitivity rules out all
twelve absolute coefficients being `1247`, so

```text
||u^(i)||_1<=11*1247+1246=14963.                      (19)
```

THM-2163 converts it, in every base `q>=2`, into an exact terminating carry
path with fewer than

```text
14963-1=14962                                         (20)
```

possible integer carry states, using THM-2163's sign-sharp positive-row
count. The symmetric general-purpose estimate remains valid but is weaker.

These are state-alphabet bounds, not depth bounds. THM-2163's nested
quotient-owner masks remain indispensable for termination; (18) does not
make the speed search finite by itself.

## 4. Boundary audit

Each hypothesis is used:

- zero safe measure invokes the Selberg relations;
- positivity turns a support-two relation into a divisibility lock;
- distinctness changes `1<=q<=105` into `2<=q<=105`;
- primitivity of `w` makes `w_i` exactly `+/-1` when nonzero;
- independence makes every contracted deletion relation nonzero.

The final constants are literal contraction costs:

```text
106*7=742,
43*29=1247,
161*7=1127,
55*22=1210,
44*27=1188.
```

No determinant division, saturation denominator, or unproved sparsity
claim is hidden in them. QED.
