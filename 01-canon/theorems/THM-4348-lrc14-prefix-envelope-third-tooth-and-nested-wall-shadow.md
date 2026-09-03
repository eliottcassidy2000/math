---
id: THM-4348
title: "LRC(14) prefix-envelope third-tooth criterion and nested wall-shadow reduction"
status: >
  PROVED ELEMENTARY PREFIX-ENVELOPE, EXACT CENTERED-RESIDUE SUCCESSOR,
  AGGREGATE SIEVE, NESTED-MULTIPLE, WALL-SHADOW, AND SHARP WALL-CAPACITY
  LEMMAS + FINITE-EXACT
  CONTROLS. Farthest-reach greedy is exactly a prefix-record recurrence. A
  proposed successor has an integer necessary-and-sufficient domination
  test retaining open endpoints and deterministic ties. Located capacity
  factors through renewal obligation, successor survival, and prefix-orbit
  reachability. Odd non-fourteen-divisible multiples of one speed preserve
  every two-tooth chain involving that speed, but are safe on all of its
  1/14 walls. Each residual speed r covers at most
  2gcd(r,u)ceil(u/[7gcd(r,u)]) of the 2u signed u-walls. This yields a
  quantitative residual boundary-cover necessity. One residual nonmultiple
  cannot cover unless the anchor is totally resonant; two residuals have one
  completely classified quotient-three equality pattern. This also yields a
  primitive full twelve-tail family on the 420|h wall with at least h/21 actual uses of
  one coprime pair, despite the remaining denominator necessities and failure
  of both half-turn clocks; every row in the family has an explicit 1/14
  witness. The result refutes unconditional bounded third-tooth reuse, not a
  counterexample-conditional competition theorem. LRC(14) remains open.
source: third_tooth / LRC14 continuation session, 2026-09-02
depends_on:
  - THM-4330-lrc14-affine-two-adic-root-types-and-anchored-pool-entry-sieve
  - THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal
related:
  - THM-4331-lrc14-safe-component-endpoint-denominator-odd-wall-escape
script: 04-computation/lrc14_third_tooth_competition_probe_third_tooth_20260902.py
output: 05-knowledge/results/lrc14_third_tooth_competition_probe_third_tooth_20260902.out
reflection: 07-reflections/lrc14-third-tooth-competition-third_tooth_20260902.md
script_sha256: 28d170e65faee92103baf279d91b4b7705e10be9223a46b8a47eea7a48486540
output_sha256: 80ea67fb849fee0c36ccef5b75283d3c4caf90218ab920c66f11113513b43b39
hash_basis: raw LF bytes
audit: >
  PASS. Two independent exact implementations of greedy selection, literal
  active teeth and the prefix maximum among all begun teeth, agree on 2,170
  small component instances. The centered-residue successor agrees with
  literal rational wall comparison on every audited selected and candidate
  edge. Visibility caps, exact ties, strict open endpoints, and the c=14
  boundary failure are checked. The wall-count identity and capacity inequality
  are checked on 14,336 exact (u,s) pairs through u=63, including a sharp
  control. The quotient-three two-residual classification is checked on 20,384
  scaled rows (128 covers), including its least-above-anchor hostile. Both
  frozen h=420 controls are reconstructed exactly. Normal, optimized, and
  hash-seeded runs reproduce the frozen transcript byte-for-byte.
---

# THM-4348 -- prefix-envelope third-tooth criterion and nested wall shadow

**PROVED ELEMENTARY LEMMAS + FINITE-EXACT CONTROLS. LRC(14) REMAINS OPEN.**

## 1. Objects, endpoint convention, and selection order

Let `h` be a positive integer and let `W` be a finite set of distinct positive
odd integers.  The physical safe components of the anchor `{2h}` are

```text
I_k=[L_k,R_k]
   =[(14k+1)/(28h),(14k+13)/(28h)],       0<=k<2h.    (1)
```

For `w in W` and `n in Z`, write its **open** danger tooth as

```text
T(w,n)=(a(w,n),b(w,n)),
a(w,n)=(14n-1)/(14w),       b(w,n)=(14n+1)/(14w).     (2)
```

Only finitely many of these teeth meet a fixed `I_k`.  Order them by the
lexicographic key

```text
(b(T),-a(T),-w(T),-n(T)).                             (3)
```

Thus the tooth reaching farthest right wins.  If right walls agree, the wider
tooth, equivalently the smaller speed, wins.  This tie rule agrees with the
deterministic exact controls in THM-4335; the theorem below does not turn a
left or right endpoint of an open tooth into an active point.

## 2. Prefix-envelope and visibility theorem

For a real frontier `x in I_k`, define `Psi_k(x)` to be the maximum-key tooth
among all teeth satisfying

```text
a(T)<x.                                               (4)
```

If there is no such tooth, or if `b(Psi_k(x))<=x`, declare a gap.  Otherwise
define recursively

```text
x_0=L_k,       T_j=Psi_k(x_j),       x_(j+1)=b(T_j). (5)
```

> **Prefix-envelope theorem.** Recurrence `(5)` is exactly the farthest-reach
> greedy interval-cover algorithm on `I_k`, including strict activity and the
> tie rule `(3)`.  It terminates with `x_(m)>R_k` if and only if the tail teeth
> cover `I_k`.  Its selected teeth then form the THM-4335 minimum-cardinality
> chain.

For a tooth `T`, define its visibility cap relative to the teeth meeting
`I_k` by

```text
c_k(T)=min(
  {b(T)} union {a(Z):Z outranks T and a(Z)<b(T)}
).                                                    (6)
```

Then, whenever `x in I_k`, the following are equivalent:

```text
T=Psi_k(x) and b(T)>x;
a(T)<x<b(T) and x<=c_k(T).                            (7)
```

The inequality at the visibility cap is nonstrict.  If `x=a(Z)`, the open
competitor `Z` has not yet begun and does not dominate `T` at `x`.

### Proof

Suppose some tooth is active at `x`.  Every tooth satisfying `(4)` but no
longer active has right wall at most `x`, while every active tooth has right
wall greater than `x`.  The maximum begun tooth is therefore exactly the
maximum active tooth.  If no active tooth exists, either no tooth has begun or
the prefix maximum has right wall at most `x`, which is exactly the stated gap
test.  This proves the recurrence and its endpoint semantics.  The standard
frontier exchange from THM-4335 then gives minimum cardinality.

A tooth `T` wins at `x` precisely when it is active and no outranking tooth
has left wall strictly below `x`.  Taking the first left wall of all
outrankers gives `(6)--(7)`. **QED.**

## 3. Exact centered-residue successor theorem

Fix an addressed outgoing tooth `A=T(u,n)` and put

```text
x=b(u,n)=(14n+1)/(14u).                              (8)
```

For every `z in W`, let `e_z` be the unique centered representative

```text
e_z = z(14n+1) mod 14u,             -7u<e_z<=7u.     (9)
```

Then a `z`-tooth is strictly active at `x` if and only if

```text
|e_z|<u.                                             (10)
```

When `(10)` holds, the active address and normalized rightward reach are

```text
ell_z=[z(14n+1)-e_z]/(14u),
R_z=(u-e_z)/z,
b(z,ell_z)-x=R_z/(14u).                              (11)
```

Let a proposed proper-crossing successor be `B=T(v,m)` and put

```text
D=um-vn,                    |u-v|<14D<u+v.           (12)
```

Then

```text
e_v=v-14D,                 R_v=(u-v+14D)/v.          (13)
```

> **Centered-residue successor theorem.** `B` is selected immediately after
> `A` if and only if no distinct `z in W` satisfies `(10)` together with
>
> ```text
> R_z>R_v,                                               (14a)
> ```
>
> or satisfies `R_z=R_v` and `z<v`.                     `(14b)`

### Proof

The active condition

```text
a(z,ell)<x<b(z,ell)
```

becomes

```text
-u<z(14n+1)-14u ell<u,                               (15)
```

which is `(10)` and uniquely determines `ell_z`.  Subtracting `x` from the
right wall gives `(11)`.  Substitution of `D` gives `(13)`.  Maximizing the
right wall is therefore equivalent to maximizing `R_z`.  At equal right wall,
the left wall is `b-1/(7z)`, so the smaller `z` is wider and wins.  This is
exactly `(14)`. **QED.**

## 4. Necessary-and-sufficient edge sieve and aggregate identity

Call a proper addressed crossing

```text
e:T(u,n)->T(v,m)                                     (16)
```

**located** when its closed handoff window satisfies

```text
[a(v,m),b(u,n)] subset I_k                           (17)
```

for some `k`.  Endpoint equality is allowed in `(17)`; activity remains
strict.  Let `E_h(W)` be the finite set of all located candidates.  Define

```text
O_h(W)={e in E_h(W):I_k is covered and not one-tooth-spanned},
S_h(W)={e in E_h(W):the proposed successor passes (14)},
A_h(W)={edges actually selected in the greedy chains}. (18)
```

Then

```text
A_h(W) subset S_h(W) intersect O_h(W) subset O_h(W). (19)
```

More strongly, a located event `e=A->B` belongs to `A_h(W)` if and only if

1. `I_k` is covered and nonspanning;
2. `B=Psi_k(b(A))`, equivalently `(14)` holds; and
3. the recurrence `(5)` reaches `A` at one of its frontiers.

Thus `(18)` is an exact obligation / successor / upstream-reach sieve, not
only a chain of upper bounds.

Group `E_h(W)` by the triple `(k,u,n)` of component and outgoing addressed
tooth, and let `G_h(W)` be the set of nonempty groups.  At most one candidate
successor survives in each group, and hence

```text
|E_h(W)|-|S_h(W)|
 =(|E_h(W)|-|G_h(W)|)+(|G_h(W)|-|S_h(W)|).           (20)
```

The first term counts excess proposed successors at one frontier; the second
counts groups whose actual frontier winner is not a located candidate in that
group.

### Proof

The first two inclusions follow from Sections 2--3 and the definition of a
renewal chain.  Conditions 1--3 say exactly that recurrence `(5)` selects `A`
and then `B`, proving the equivalence.  A deterministic maximum has at most
one winner at each frontier, proving `(20)`. **QED.**

## 5. Nested-multiple preservation and wall shadow

Let `u` be positive and let `c>=2` be an integer with `14` not dividing `c`.
At either boundary of any `u`-tooth,

```text
t=(14n+-1)/(14u),
||(cu)t||=||c/14||>=1/14.                            (21)
```

> **Nested-multiple preservation theorem.** Let `I` be a closed interval and
> let a farthest-reach greedy cover for some finite tooth bank consist of
> exactly two teeth `A,B`, one of which has speed `u`.  Adjoin any finite set
> of new, distinct speeds `c_i u`, where every `c_i>=2` and
> `14` does not divide `c_i`.  Then the greedy chain on `I` remains exactly
> `(A,B)`, with the same addresses and order.

> **Wall-shadow theorem.** Define the finite `u`-wall set
>
> ```text
> X_u={ (14n+sigma)/(14u) mod 1:
>       n in Z/uZ, sigma in {+1,-1} }.                (22)
> ```
>
> Let `C` be any finite set of positive integers not divisible by `14`, and
> let `R` be any finite set of positive integer speeds.  For
>
> ```text
> S={2h,u} union {cu:c in C} union R,                 (23)
> ```
>
> strict failure `G_S=empty` implies the following quantified residual cover:
>
> ```text
> for every x in X_u with ||2hx||>=1/14,
> there exists r in R with ||rx||<1/14.               (24)
> ```

For positive integers `u,s`, put

```text
d_s=gcd(u,s),             q_s=u/d_s,       S_s=s/d_s,
C_q(a)=#{y in Z:-q<y<q and y=a mod 14},
N_u(s)=d_s[C_(q_s)(S_s)+C_(q_s)(-S_s)],              (25a)
B_u(s)=2d_s ceil(q_s/7).                             (25b)
```

> **Wall-capacity theorem.** Among the `2u` distinct points of `X_u`, the
> speed `s` is strictly dangerous at exactly `N_u(s)` points, and
> `N_u(s)<=B_u(s)`. Therefore
> strict failure in `(23)` necessarily implies
>
> ```text
> 2u <= B_u(2h)+sum_(r in R)B_u(r).                   (26)
> ```
>
> In particular, the strict reverse inequality is a sufficient safety
> certificate for `(23)`.

### Proof

Equation `(21)` says that no new `c_i u` danger tooth crosses an `u`-tooth
boundary.  A new tooth active inside an `u`-tooth is contained in that parent
tooth and ends no later; at an equal right wall, the parent speed `u<c_i u`
wins the tie.

If the `u`-tooth is first in the two-tooth chain, every new tooth active at the
initial frontier is contained in it and cannot beat it, while `(21)` makes all
new teeth inactive at its outgoing wall.  If the `u`-tooth is second, a new
tooth active before its left wall expires no later than that wall and cannot
beat the first tooth; after that wall it is contained in the `u`-tooth and
cannot beat its parent.  This proves preservation for either order and any
finite collection.

For `x in X_u`, the speed `u` is equality-safe and every `cu` is safe by
`(21)`.  If the anchor is also safe and no `r in R` is dangerous, then `x` is
a physical witness for the full row `(23)`.  Taking the contrapositive proves
`(24)`. **QED.**

For a fixed sign `sigma`, the danger condition is equivalent to

```text
y=14S_s n+S_s sigma-14q_s m,             |y|<q_s.    (26a)
```

The residue `y=S_s sigma mod 14` uniquely determines `n mod q_s`, and each
such residue is represented `d_s` times modulo `u`. Summing the two signs
proves the exact formula `N_u(s)` in `(25a)`.

Equivalently, as `n` runs modulo `u`, each residue modulo `q_s` occurs exactly
`d_s` times, and

```text
s(14n+sigma)/(14u)
 = S_s n/q_s + S_s sigma/(14q_s).                    (27)
```

Since `gcd(S_s,q_s)=1`, this is a translate of the uniform `q_s`-point grid,
with multiplicity `d_s`. An open circular interval of length `1/7` contains
at most `ceil(q_s/7)` points of that grid. Summing the two signs proves the
bound `B_u(s)`. Under strict failure, the anchor and residual danger sets must
cover all `2u` walls by `(24)`; the union bound gives `(26)`. **QED.**

The capacity bound is attained, for example, by

```text
u=31,                 s=433=14u-1,
#(X_u intersect D_s)=10=B_u(s).                      (28)
```

Because `u` is odd, every proper quotient `q_s>1` is odd and at least three.
The elementary inequality

```text
ceil(q_s/7)<=q_s/3                                   (28a)
```

is strict except at `q_s=3`. Consequently every speed not divisible by `u`
satisfies

```text
N_u(s)<=2u/3.                                        (28b)
```

Equality in `(28b)` forces `q_s=3`. For odd `s`, equality then holds exactly
when

```text
S_s=s/d_s = +-1 mod 14.                              (28c)
```

For the even anchor `s=2h`, equality at quotient three holds exactly when
`2h/d_s=0,+-2 mod 14`.

> **One-residual corollary.** Let `u` and `r` be positive odd integers with
> `u` not dividing `r`, and let `C` be any finite set of positive integers not
> divisible by `14`. If
>
> ```text
> 7u does not divide h,                              (28d)
> ```
>
> then the row
>
> ```text
> {2h,u,r} union {cu:c in C}                         (28e)
> ```
>
> is safe at a point of `X_u`.

Indeed, if `u` does not divide `h`, the anchor kills at most `2u/3` walls and
at least `4u/3` are anchor-safe. The residual `r` kills at most `2u/3` of
them. If `u|h` but `(28d)` holds, the anchor is safe on every `u`-wall, while
the same residual bound applies. The excluded case is exact: when `7u|h`,
the anchor is dangerous on all of `X_u`.

> **Two-residual equality classification.** Let `u,r_1,r_2` be positive odd
> integers, with `u` dividing neither residual, and assume `(28d)`. Then the
> anchor and the two residual danger sets cover all of `X_u` if and only if
> there are a positive odd `d` and positive integers `H,R_1,R_2` such that, after possibly
> exchanging the residuals,
>
> ```text
> u=3d,                 h=dH,              r_i=dR_i,
> H=+-7 mod 21,
> R_1=+-1 mod 42,       R_2=+-13 mod 42.              (28f)
> ```

For necessity, `u|h` would make the anchor safe on all walls under `(28d)`,
so two residuals could cover at most `4u/3` walls. Thus `u` does not divide
`h`. The union bound and `(28b)` now force equality for the anchor and both
residuals, and their three size-`2u/3` masks must be disjoint. Hence all three
gcd quotients equal three. On the six quotient walls

```text
{+-1,+-13,+-15} mod 42,                              (28g)
```

an odd equality residual kills either the `+-1` pair (`R=+-1 mod 42`) or the
`+-13` pair (`R=+-13 mod 42`). It never kills the `+-15` pair. The anchor must
therefore kill `+-15`, which is exactly `2H=+-14 mod 42`, or
`H=+-7 mod 21`. These three masks are disjoint and exhaustive, also proving
sufficiency.

The least positive representatives are `(r_1,r_2)=(1,13)`. If both residuals
are required to exceed the anchor `2h`, the least representative hostile is

```text
(h,u,r_1,r_2)=(7,3,29,41).                           (28h)
```

Its anchor and residual masks are the three disjoint two-point classes in
`(28g)`. It is a hostile to the wall certificate, not an unsafe LRC row.

There is a useful primitive `420|h` consequence. Under `(28d)`, if `420|h`
and the complete row `(23)` with two residual nonmultiples is primitive, then
the equality case `(28f)` is impossible. Indeed all displayed speeds would
share `d`; primitivity forces `d=1`, after which `420|h` contradicts
`H=+-7 mod 21`. Thus such a row is safe on `X_u`. The exact untreated
wall-shadow resonance is `7u|h`.

The divisibility hypothesis is the exact boundary-crossing condition.  In
the general, not necessarily odd-tail setting, it cannot be dropped:

```text
h=26, k=35, base chain (3,2)->(29,20);
adjoining 42=14*3 changes it to (3,2)->(42,29).       (29)
```

## 6. Primitive full twelve-tail unbounded-reuse family

For every integer `M>=1`, put

```text
N=10M,                 h=420M,
u=140M+1,              v=140M+3,
C=(9,11,13,39,41,43,45,47,49,51),
W_M={u,v} union {cu:c in C}.                          (30)
```

Then:

1. `W_M` consists of twelve distinct positive odd speeds;
2. `{2h} union W_M` is primitive because `gcd(u,v)=1`;
3. `420|h`, and the tails `9u,11u,13u` satisfy the remaining THM-4330
   divisibility necessities for `9,11,13`;
4. the multipliers `39,41,43,45,47` obey

   ```text
   12h<cu<16h,                                         (31)
   ```

   so neither THM-4330 half-turn clock is a witness;
5. the same unordered coprime pair `{u,v}` occurs in at least

   ```text
   2N=20M=h/21                                        (32)
   ```

   actual greedy transitions across the `2h` anchor components;
6. nevertheless the full row is safe at

   ```text
   t_M=1/(14u),              min_(s in {2h} union W_M)||s t_M||=1/14. (33)
   ```

### Proof of the transition count

Before adjoining the multiples, consider the two-speed bank `{u,v}`.  Since
`gcd(u,v)=1`, every determinant

```text
d=N+1,...,2N                                           (34)
```

gives one proper crossing in each orientation.  For `u->v`, when `d=2r`, take

```text
n=u-r,              m=v-r,              k=84N-6r.    (35)
```

When `d=2r-1`, take

```text
n=7N+1-r,          m=7N+2-r,           k=42N-6r+3.  (36)
```

Direct substitution gives `um-vn=d`.  In anchor-normalized coordinates, the
incoming left and outgoing right walls are respectively

```text
d=2r:     6(3r-N)/v,          6(r+N)/u;
d=2r-1:   3(6r-2N-3)/v,      3(2N+2r-1)/u.          (37)
```

For the range `(34)` these values lie in `[1/14,13/14]`.  The normalized
`u`-tooth width is `12N/u`, so its left wall is strictly below `1/14`; the
normalized `v`-tooth width is `12N/v`, so its right wall is strictly above
`13/14`.  Hence `(35)--(36)` give literal two-tooth greedy covers.  Reflection
gives the opposite orientation.  The component addresses are distinct:
the four parity/orientation classes are distinct modulo six.  Thus there are
at least `2N` pair uses.

Every `c in C` is odd and hence is not divisible by `14`.  The nested-multiple
preservation theorem keeps all `2N` covers after the ten speeds `cu` are
adjoined, proving `(32)`.

Finally, at `t_M` the speed `u` is tight, every `cu` is safe by `(21)`, and

```text
v t_M=1/14+1/(7u),
||(2h)t_M||=60M/(140M+1),                            (38)
```

where the last quantity lies in `[1/14,1/2)`. This proves `(33)`. The
wall-capacity certificate independently applies with `R={v}`: here
`gcd(u,v)=1` and `gcd(u,2h)` is either `1` or `3`, so

```text
B_u(2h)+B_u(v)<2u.                                  (39)
```

This also proves safety. **QED.**

## 7. Exact controls and failure anatomy

For the two inherited `h=420` rows, the exact all-candidate / renewal-
obligation / successor-survivor / reached-edge sieves are

```text
P=1287:       21,176 -> 50 -> 10 -> 6,
P=9009:       27,404 -> 86 -> 11 -> 6.               (40)
```

Both have `726` uncovered, `110` one-tooth-spanned, and `4` renewal
components.  The pair `{525,1365}` has 105 located candidates in each
orientation, but each orientation splits as

```text
90 uncovered-component candidates,
14 one-tooth-span candidates,
 1 renewal candidate, which survives and is reached. (41)
```

Thus its apparent `210 -> 2` collapse is almost entirely the obligation
filter, not third-tooth domination.  On all pairs, both later stages in `(40)`
remain genuine.

For `(30)`, exact controls at `M=1,2,5` retain `24,50,122` base pair uses and
have `26,50,124` pair uses after all ten multiples are adjoined.  They return
the exact witnesses `(33)` with clearance `1/14`. Their anchor-plus-residual
wall-capacity bounds are respectively

```text
84<282,                 164<562,                 404<1402. (42)
```

## 8. Scope and next reduction

The theorem proves a precise dichotomy target, not the minority branch.

- For many nonnested tails, use `(14)` and the visibility cap to bound
  successors and upstream prefix-orbit reach on genuine renewal components.
- For a large nesting class, use `(24)--(26)` to delete that class on `X_u`
  and bound the capacity of the small residual nonmultiple bank on every
  anchor-safe wall.

The family `(30)` is safe and therefore does not refute a theorem conditioned
on strict failure.  It does refute every unconditional `O(1)` reuse claim
based only on primitivity, pair coprimality, twelve-tail competition, `420|h`,
the `9/11/13` necessities, and failure of the two half-turn clocks.  Neither a
uniform residual wall-cover bound, the full `2+12` branch, arbitrary entry,
nor LRC(14) follows.

## 9. Reproduction

From the repository root:

```bash
python3 -B 04-computation/lrc14_third_tooth_competition_probe_third_tooth_20260902.py \
  | diff -u 05-knowledge/results/lrc14_third_tooth_competition_probe_third_tooth_20260902.out -
python3 -O -B 04-computation/lrc14_third_tooth_competition_probe_third_tooth_20260902.py \
  | diff -u 05-knowledge/results/lrc14_third_tooth_competition_probe_third_tooth_20260902.out -
PYTHONHASHSEED=913 python3 -B \
  04-computation/lrc14_third_tooth_competition_probe_third_tooth_20260902.py \
  | diff -u 05-knowledge/results/lrc14_third_tooth_competition_probe_third_tooth_20260902.out -
```
