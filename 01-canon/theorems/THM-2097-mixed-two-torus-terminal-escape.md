---
id: THM-2097
title: "Every rank-two seven-terminal plane has a mixed-threshold torus escape"
status: >
  PROVED. For one odd guard character and seven distinct positive terminal
  characters spanning a two-torus, there is a torus point outside the
  radius-1/7 guard and outside all seven radius-1/14 dangers. The proof splits
  off guard-proportional characters using THM-2080, while the transverse case
  uses an exact seven-band mass equality and a pair-rigidity lemma. Primitive
  geodesic density then leaves only finitely many obstructing directions in
  each rank-two coefficient template. Rank-one terminal restrictions freeze
  and THM-2077 bounds their outer rows. Hence every rank-seven depth-four
  template, including THM-2087's bounded guard-ratio branch, is finite
  template-by-template. The finite rows are not discharged, so this is not
  LRC(14).
source: codex-2026-07-22-LRC-mixed-two-torus-terminal
depends_on:
  - THM-2052
  - THM-2053
  - THM-2077
  - THM-2080
related:
  - THM-2092
  - THM-2093
  - THM-2095
---

# THM-2097 -- mixed-threshold escape on every terminal two-plane

Let

```text
c_0,c_1,...,c_7 in Z^2                                (1)
```

be character columns on `T^2`. Assume they span `Q^2`, and assume some
primitive `d in Z^2` gives

```text
h=c_0.d>0 odd,
q_i=c_i.d>0 pairwise distinct,          1<=i<=7.       (2)
```

Then the mixed safe cell

```text
Omega={X in T^2:
       ||c_0.X||>1/7 and ||c_i.X||>1/14 for all i}     (3)
```

is nonempty and open.

## 1. Guard-proportional columns

Call `c_i` vertical if it is rationally proportional to `c_0`, and let `r`
be the number of vertical terminal columns. A unimodular coordinate change
puts

```text
c_0=(0,b),                  b odd.                     (4)
```

After reversing the second coordinate if needed, the vertical columns are

```text
c_i=(0,a_i),                a_i>0 distinct.            (5)
```

Indeed (2) makes the second parameter coordinate odd and of fixed sign, so
positivity and distinctness transfer to the `a_i`.

Let

```text
C_b={y:||by||>1/7}.                                    (6)
```

For a vertical column, THM-2080 gives

```text
measure(D_(a_i) intersect E_b)>=1/42,                 (7)
```

with equality only at `a_i=6b`. Hence its danger mass inside `C_b` is at
most

```text
1/7-1/42=5/42.                                        (8)
```

If `r<=5`, the union of the vertical dangers has mass at most
`5r/42<5/7=measure(C_b)`. If `r=6`, equality in that union bound would force
all six instances of (7) to be equalities, hence all six `a_i` to equal
`6b`, contrary to distinctness. Thus for every `1<=r<=6` there is a `y` with

```text
||by||>1/7,             ||a_i y||>1/14                (9)
```

for all vertical columns. Null endpoints do not affect the strict positive
mass gap, so `y` may be chosen with strict inequalities.

Every remaining column has form `(A_i,B_i)` with `A_i!=0`. At the fixed `y`,
its closed danger set as a function of `x` has measure exactly `1/7`, because
`x |-> A_i x+B_i y` preserves Haar measure. There are `7-r<=6` such sets, so
their union cannot cover the `x`-circle. Choosing `x` outside them and using
(9) gives a point of `Omega`.

The case `r=7` would make all eight columns proportional and contradict the
rank-two hypothesis. It remains only to prove the case `r=0`.

## 2. Pair rigidity in the transverse case

We first isolate the needed character lemma.

> **Pair-rigidity lemma.** Let `f,f',g:T^2->T` be integer characters such
> that `(f,g)` and `(f',g)` are surjective. If
>
> ```text
> ||f(X)||<1/14 and ||f'(X)||<1/14
>       imply ||g(X)||<=1/7,                            (10)
> ```
>
> then
>
> ```text
> g=epsilon f+eta f',           epsilon,eta in {+-1}.  (11)
> ```

### Proof

Put `Phi=(f,f')`. If `Phi` has rank one, its connected kernel contains the
kernel of a nonzero primitive character. Surjectivity of `(f,g)` makes `g`
nonconstant, indeed surjective, on that kernel. This contradicts (10), since
both `f` and `f'` vanish there.

Thus `Phi` has rank two and is a finite surjective torus map. On `ker(Phi)`,
equation (10) confines the finite subgroup `g(ker(Phi))` to the closed arc of
radius `1/7`. A nontrivial finite subgroup of the circle contains a point of
distance at least `1/3` from zero, so this subgroup is trivial. Therefore `g`
factors through `Phi`:

```text
g=a f+b f',                    a,b in Z.               (12)
```

Both coefficients are nonzero, since otherwise one of `(f,g)` and `(f',g)`
would have rank one.

Put `K=|a|+|b|`. If `K>=3`, choose

```text
1/(7K)<t<min(1/14,1/(2K)).                             (13)
```

Surjectivity of `Phi` supplies an `X` with

```text
f(X)=sign(a)t,                 f'(X)=sign(b)t.
```

Both characters are dangerous, while (12)--(13) give

```text
1/7<||g(X)||=Kt<1/2,
```

contradicting (10). Hence `K<=2`; since `a,b` are nonzero, both are `+-1`.
This proves (11). QED.

## 3. Seven transverse bands cannot cover

Assume `r=0`. Put

```text
A_i={X:||c_i.X||<=1/14},
C={X:||c_0.X||>1/7}.                                  (14)
```

Every map `(c_i,c_0):T^2->T^2` has nonzero integer determinant and is
surjective. Consequently

```text
measure(A_i intersect C)=(1/7)(5/7)=5/49.             (15)
```

The seven masses in (15) sum to `5/7=measure(C)`. If the `A_i` covered `C`,
equality in the union bound would force every pair intersection in `C` to
have measure zero. The corresponding intersection of the two **open** danger
bands with the open set `C` would then be an open null set, hence empty.
Thus every pair `(c_i,c_j)` would satisfy the implication (10), with
`g=c_0`.

Fix `c_1`. The pair-rigidity lemma says that every `j=2,...,7` obeys

```text
c_0=epsilon_j c_1+eta_j c_j.
```

Up to changing the irrelevant sign of `c_j`, its danger band is therefore
one of only

```text
c_0+c_1,                 c_0-c_1.                     (16)
```

Among six bands, two coincide. That repeated open band has positive
intersection with `C`, because its character and `c_0` form a surjective pair.
This contradicts the pairwise-null conclusion. Hence the `A_i` do not cover
`C`, and any point in their complement inside `C` lies in `Omega`. This
completes the proof of (3). QED.

## 4. Primitive geodesics leave only a finite disk

Choose `X_0 in Omega` and let `sigma>0` be the minimum slack in its eight
strict inequalities. Put

```text
L=max_(0<=i<=7)||c_i||_2.                              (17)
```

THM-2053's elementary primitive-geodesic estimate says that every point of
`T^2` lies within distance `1/(2||d||_2)` of

```text
C_d={(d_1t,d_2t):t in T}.                              (18)
```

Every character `c_i` is `||c_i||_2`-Lipschitz. Therefore

```text
||d||_2>L/(2sigma)                                     (19)
```

puts some point of `C_d` inside `Omega`. At the corresponding time `t`,
equations (2)--(3) say

```text
||ht||>1/7,                  ||q_i t||>1/14 for all i. (20)
```

Thus a terminal containment `G_Q subset E_h` forces `d` into the finite
lattice disk complementary to (19).

## 5. Consequence for the depth-four LRC branch

Apply THM-2052 to an original thirteen-speed row. Its bounded relation code
has rank twelve, in which case the row is already finite, or rank eleven,
with a two-dimensional rational solution plane. There are only finitely many
such coefficient templates.

Inside a fixed depth-four chart, pass to the finite-index lattice on which the
normalized terminal coordinates

```text
(h,q_1,...,q_7)                                        (21)
```

are integral. If the terminal restriction has rank two, Sections 1--4 leave
only finitely many primitive parameter directions. Primitivity of the
original row makes its direction primitive in this depth lattice.

If the terminal restriction has rank one, its image is a rational line
through the primitive vector `(h,Q)`. Any other admissible normalized terminal
is `lambda(h,Q)`. Integrality and `gcd(Q)=1` give `lambda in Z`; terminal-core
primitivity and positivity force `lambda=1`. Thus the terminal block is
literally frozen. THM-2077 then gives

```text
max(S)<=128 max(Q)/3,                                  (22)
```

so only finitely many outer rows occur.

Hence **every rank-seven depth-four coefficient template is finite**. This
includes THM-2087's bounded guard-ratio branch: a persistent bounded-ratio
relation creates a vertical column handled by Section 1, a nonpersistent
relation cuts the plane to rank one, and a rank-one terminal restriction
freezes as above. THM-2093 gives a stronger uniform box on the no-pair global
star; the present result is a mixed-threshold, template-specific closure that
also covers the separate guard-ratio branch. THM-2095's p-adic program may
still sharpen this to a small uniform commensurability scale.

This is finiteness, not finite discharge. Rank-seven terminal rows inside the
resulting disks and terminal sizes eight through eleven at smaller tower depth
remain to be decided.

## 6. Two exact interval taxes

The same guard-local viewpoint gives cheap prefilters for that finite work.

### Parity tax

Let `o` be the number of odd terminal speeds and put

```text
M_2=max({h} union {q:q odd}).                           (23)
```

On

```text
I=[1/2-5/(14M_2),1/2+5/(14M_2)],
```

the odd guard and every odd speed are safe, so only even speeds can cover.
The exact interval discrepancy

```text
measure(D_q intersect I)<=|I|/7+6/(49q)               (24)
```

therefore forces

```text
M_2 sum_(2|q)1/q>=5o/6.                                (25)
```

### Seven-carrier tax

Assume `7` does not divide `h`, put `k=#{q:7|q}`, and let

```text
M_7=max({h} union {q:7 does not divide q}).             (26)
```

Choose `a` with `ah=+-3 (mod 7)` and use

```text
J=[a/7,a/7+1/(14M_7)].
```

The guard and every noncarrier speed are safe throughout `J`; only the `k`
seven-divisible combs can cover it. Applying (24) gives

```text
M_7 sum_(7|q)1/q>=(7-k)/12.                            (27)
```

In particular their smallest speed obeys

```text
q_min<=2M_7, (24/5)M_7, 9M_7, 16M_7
```

for `k=1,2,3,4`, respectively. Endpoint equalities lie outside the strict
danger sets, and measure identities ignore the finite boundary set.

## 7. Assumption challenge and Tournament Analysis

The challenged assumption is that a rank-two terminal needs either a global
scalar margin or another bounded relation. The mixed torus retains the two
different radii and the proportional/transverse split; its open cell gives a
direct finite-direction terminal.

Candidate vertices were characters, proportionality classes, danger bands,
pair-cover obligations, and primitive directions. An orientation of the seven
terminal characters loses the vertical count and the equality ledger (15).
The faithful carriers are the character matroid plus actual integer columns
and the two radii. In the transverse case the pair-rigidity proof reduces six
bands to two sign classes; this is a pigeonhole obstruction, not a tournament
cycle or score invariant. QED.
