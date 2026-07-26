---
id: THM-2385
title: "Two-top septimal blocker-collision closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In
  THM-2367's only-c_3-dominant k=2, (t,b)=(2,0) alternative,
  the guard, three lower q masks, and two lower blocker masks have
  total width exactly one and therefore form a one-fold partition off
  the two top q masks and c_3. A thirteen-root fibre then forces the
  quotient containment D_C1 intersect D_C2 inside closure(D_C3).
  Unequal lower septimal depths give exactly N/49 intersection points
  on every generic N=7^(M+1) fibre and hence an exact 6/343
  anti-shield, so the two lower blockers would have to have equal
  7-adic depth. After scaling that common depth, a uniform two-source
  lemma proves D_a intersect D_b is never contained in closure(D_c)
  when 7 does not divide ab and 49 divides c. Hence the entire (2,0)
  alternative is empty. Independently, every strict thirteen-adic row
  already fails after a second thirteen-root descent and an exact
  6/49 anti-shield. Equal-depth repeated-first hostiles show why the
  uniform joint-comb lemma is genuinely stronger than either marginal
  fibre argument. No entire thirteen-adic row, ledger entry, other
  septimal alternative, target, or LRC(14) conclusion is removed.
source: codex-2026-07-26-two-top-blocker-collision
depends_on:
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
related:
  - THM-2377-septimal-valuation-collision-and-bockstein-carry-gate
  - THM-2378-hard-septimal-root-stalk-closure
  - THM-2381-one-top-one-blocker-septimal-root-stalk-closure
  - THM-2382-saturated-septimal-seven-bin-root-fibre-closure
script: 04-computation/lrc14_two_top_blocker_collision_reduction_thm2385.py
output: 05-knowledge/results/lrc14_two_top_blocker_collision_reduction_thm2385.out
script_sha256: c1bf13bd3d4d56399a4647e71e5d3b42471e39f30ae085ae7e272bc9e9c5dca3
output_sha256: a86473f6f0e61fffa034edb2806b5661e4d2911811ecbe4313bc795cb138645f
hash_basis: working-tree bytes (LF)
---

# THM-2385 -- the two-top septimal blocker collision cannot occur

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2367 leaves the nonsaturated alternative

```text
k=2,                    (t,b)=(2,0).                (1)
```

Here two `q` masks lie at the septimal top layer, both first blockers
lie strictly below it, and `c_3` lies strictly above it. The faithful
object is no longer a single low-blocker word as in THM-2381, but the
intersection word of the two low blockers:

```text
seven-bin exactness
  -> two low blockers cannot be simultaneously dangerous away from c_3
  -> D_C1 intersect D_C2 is shielded by D_C3;

unequal lower 7-depths
  -> transverse intersection of density 1/49 on every fibre
  -> exact 6/343 shield failure;

strict 13-depths
  -> expose a second thirteen-root stalk
  -> one nonempty unit word against a constant danger/safe pair
  -> exact 6/49 shield failure;

equal lower 7-depths
  -> scale to two 7-unit source combs and one 49-divisible target
  -> common zeros + tooth centres + boundary perturbation
  -> the joint source intersection cannot fit in the target.        (2)
```

The two marginal descents isolate a doubly ramified collision in which
the blockers have both repeated-first thirteen-adic depth and equal
lower septimal depth. A new joint-comb lemma then excludes that last
collision. Thus (1) is empty.

## 1. The two-top branch

Retain a canonical first-depth-one scalar cover

```text
C_H subset
  union_(i=1)^5 D_(q_i)
  union D_(c_1) union D_(c_2) union D_(c_3)          (3)
```

almost everywhere, where

```text
D_v={x in T:||vx||<1/14},

C_H={x in T:||Hx||>1/7},

E_H={x in T:||Hx||<1/7}.                             (4)
```

The six labels `H,q_1,...,q_5` are units modulo thirteen and all three
blockers are divisible by thirteen. Put

```text
M=max(nu_7(H),nu_7(q_1),...,nu_7(q_5))>0.           (5)
```

Assume

```text
nu_7(H)<M,

nu_7(c_3)>M,

nu_7(q_a)=nu_7(q_b)=M,

nu_7(q_i)<M                         for i notin {a,b},

nu_7(c_1),nu_7(c_2)<M.                              (6)
```

These are precisely THM-2367's `k=2,(t,b)=(2,0)` hypotheses. Define

```text
mathcal A=D_(q_a) union D_(q_b) union D_(c_3),

L
 =1_(E_H)
  +sum_(i notin {a,b})1_(D_(q_i))
  +1_(D_(c_1))
  +1_(D_(c_2)).                                     (7)
```

Measured in ordinary `1/7`-arc units, the lower load has total width

```text
2+3+1+1=7.                                          (8)
```

## 2. Exact lower load off the two top masks

Put

```text
N=7^(M+1)
```

and consider a generic additive fibre

```text
O_z={z+j/N:j in Z/NZ}.                              (9)
```

The word of `c_3` is constant on this fibre because `N|c_3`. Each top
`q` word occupies one complete bin `j mod 7`; the two words may occupy
the same bin. Every lower ordinary word contributes `N/49` points to
each bin, and the lower guard contributes `2N/49`. Thus (8) gives
exactly

```text
2N/49+3N/49+N/49+N/49=N/7                          (10)
```

lower incidences in every bin.

On a `c_3`-safe fibre, any bin occupied by neither top `q` word must
be covered by the lower masks. Its lower incidence count equals its
cardinality, so the coverage is an exact one-fold partition. There are
at least five such bins, whether or not the two top words collide.
Consequently

```text
L=1                         almost everywhere on mathcal A^c. (11)
```

Genericity here excludes the finite endpoint set and the null
exceptional set in the almost-everywhere scalar cover. Equation (11)
does not assert exactness inside either top word or on `D_(c_3)`.

## 3. The first thirteen-root descent forces blocker intersection

Write

```text
c_j=13C_j,                         j=1,2,3.          (12)
```

Choose a representative `y in [0,1)` outside all blocker endpoints
and the finitely many inverse-root pullbacks of the exceptional set in
(11). Its thirteen inverse roots are

```text
x_h=(y+h)/13,                         h=0,...,12.    (13)
```

Every blocker word is constant on this root fibre:

```text
1_(D_(c_j))(x_h)=1_(D_(C_j))(y).                    (14)
```

Suppose

```text
y in D_(C_1) intersection D_(C_2)
      minus closure(D_(C_3)).                       (15)
```

Then both lower blocker summands of `L(x_h)` equal one on all thirteen
roots, while `c_3` is safe. Each top `q` word contains at most two
roots, because its speed is a thirteen-unit and an arc of length
`1/7` contains at most two points of a thirteen-grid. Their union
therefore contains at most four roots. At any of the at least nine
remaining roots, (11) gives `L(x_h)=1`, contradicting the two blocker
contributions.

After removing the null exceptional base set, this proves the
necessary containment

```text
D_(C_1) intersection D_(C_2)
  subset closure(D_(C_3))                 almost everywhere. (16)
```

The statement is deliberately about the intersection word. Neither
individual low blocker is forced below `D_(C_3)`.

## 4. Unequal septimal depths create an exact `6/343` anti-shield

The intersection in (16) is transverse whenever its two lower speeds
have distinct septimal depths.

> **Unequal-depth intersection lemma.** Let positive integers
> `U,V,B` satisfy
>
> ```text
> nu_7(U),nu_7(V)<M<nu_7(B),
>
> nu_7(U)!=nu_7(V).                                  (17)
> ```
>
> Then
>
> ```text
> measure(
>   D_U intersection D_V minus closure(D_B)
> )=6/343.                                           (18)
> ```

**Proof.** Assume

```text
r=nu_7(U)<s=nu_7(V)<M
```

after interchanging `U,V`, and put `N=7^(M+1)`. On every generic
inverse fibre

```text
y_j=(z+j)/N,                         j in Z/NZ,
```

each lower word has `N/7` points. To count their intersection, use the
`U`-phase as coordinate after removing `7^r`.
The `U` word is a cyclic interval of

```text
7^(M+1-r)/7
```

grid points. The `V` word is periodic in that coordinate with period

```text
7^(M+1-s)
```

and occupies one seventh of every period. The number of complete
periods inside the `U` interval is

```text
7^(s-r-1),
```

an integer because `r<s`. Hence the intersection has exactly

```text
N/49                                               (19)
```

points on every generic fibre.

The `B` word is constant on the fibre and equals the word of `B/N` at
the base `z`. Its quotient-safe base set has measure `6/7`, so
disintegration under multiplication by `N` gives

```text
(6/7)(1/49)=6/343.                                  (20)
```

This proves (18). Endpoints are null. QED.

Apply the lemma with `(U,V,B)=(C_1,C_2,C_3)`. Since division by
thirteen does not change septimal depth, (16) and (18) are
incompatible unless

```text
nu_7(c_1)=nu_7(c_2)<M.                              (21)
```

Thus the two-source collision must lie in one repeated septimal layer.

## 5. Strict thirteen-adic rows are impossible

The canonical blocker-depth profiles are

```text
nu_13(c_1,c_2,c_3)=(1,beta,gamma),                  (22)
```

with either

```text
2<=beta<gamma                         [strict],

beta=1<gamma                           [repeated first]. (23)
```

Suppose the row is strict. After (12),

```text
13 does not divide C_1,

C_2=13A_2,

C_3=13A_3.                                          (24)
```

Assume a base phase

```text
y in D_(A_2) minus closure(D_(A_3)).                (25)
```

which is also outside all thirteen inverse-root pullbacks of the null
exceptional set in (16). Such a choice exists whenever (25) has
positive measure; removing finitely many null pullbacks does not change
that measure.

On its thirteen inverse roots, `C_2` is dangerous everywhere and
`C_3` is safe everywhere. But `C_1` is a thirteen-unit, so its ordinary
root word has one or two elements. At any element of that nonempty
word, the containment (16) fails. Therefore (16) would force

```text
D_(A_2) subset closure(D_(A_3))          almost everywhere. (26)
```

This second descent is where strictness is load-bearing: it makes one
lower quotient word constant while leaving the other a nonempty unit
word.

Division by thirteen again preserves septimal depth, and (6) gives

```text
nu_7(A_2)<M<nu_7(A_3).                               (27)
```

The corresponding one-source anti-shield is exact:

```text
measure(D_(A_2) minus closure(D_(A_3)))
 =(6/7)(1/7)
 =6/49.                                             (28)
```

Indeed `A_3` is constant on every `N`-root fibre and is safe on base
mass `6/7`, while the lower `A_2` word occupies exactly `N/7` roots on
each generic fibre. Equations (26) and (28) contradict one another.
Hence no strict thirteen-adic row realizes (1).

Combining Sections 4 and 5 gives the exact reduction:

> If a scalar cover realizes `k=2,(t,b)=(2,0)`, then
>
> ```text
> nu_13(c_1,c_2,c_3)=(1,1,gamma),       gamma>1,
>
> nu_7(c_1)=nu_7(c_2)<M<nu_7(c_3).                  (29)
> ```

## 6. The equal-depth joint-comb lemma closes the residual

It remains to use the intersection in (16) as a joint comb rather than
two marginal fibre counts. Put

```text
r=nu_7(C_1)=nu_7(C_2),

C_1=7^r a,              C_2=7^r b,              C_3=7^r c. (30)
```

Since `r<M<nu_7(C_3)`,

```text
7 does not divide ab,                   49 divides c.       (31)
```

Multiplication by `7^r` is onto and Haar-measure preserving on `T`,
and each mask in (16) is the pullback of the corresponding mask in
`(a,b,c)`. It first gives the following containment almost
everywhere. But

```text
(D_a intersection D_b) minus closure(D_c)
```

is open. Every nonempty open subset of `T` has positive Haar measure,
so an almost-everywhere containment upgrades here to the literal
containment

```text
D_a intersection D_b subset closure(D_c).          (32)
```

The following uniform lemma rules this out.

> **Two-layer single-target lemma.** If `a,b,c` are positive integers
> satisfying
>
> ```text
> 7 does not divide ab,                    49 divides c,
> ```
>
> then
>
> ```text
> D_a intersection D_b minus closure(D_c)
>   contains a nonempty open interval.               (33)
> ```

**Proof.** Suppose the literal containment (32) held, and put

```text
g=gcd(a,b).
```

For every integer `k`, the point `x=k/g` is a common zero of the two
source phases. Hence every element of the cyclic subgroup

```text
{ck/g mod 1:k in Z}                                  (34)
```

would have to lie in the closed centred arc of radius `1/14`. Any
nontrivial cyclic subgroup has an element at centred distance at least
`1/3`: for order two the distance is `1/2`, and for order at least
three take a residue nearest one half. At the corresponding common
source zero, all inequalities are strict and persist on a
neighbourhood, contradicting even almost-everywhere containment.
Therefore (34) is trivial and

```text
g divides c.                                         (35)
```

Divide `a,b,c` by `g` and use the surjective change of variable
`y=gx`. Since `7` does not divide `g`, the hypotheses (31) persist.
We may therefore assume

```text
gcd(a,b)=1.                                          (36)
```

If `a=b`, then (36) gives `a=b=1`. The common central source interval
would force `c<=1`, contradicting `49|c`. Hence, after interchanging
the sources, assume

```text
0<a<b.                                               (37)
```

The interval

```text
0<x<1/(14b)
```

lies in both source dangers. If `c>b`, choose

```text
1/(14c)<x<min(1/(14b),13/(14c)).
```

Then both sources are dangerous while `c x` lies strictly between two
adjacent `c`-teeth. Thus containment forces `c<=b`. Equality is
impossible because `7|c` and `7` does not divide `b`, so

```text
0<c<b,                         b>=50.                (38)
```

Now inspect the centres `x=n/b` of the `b`-teeth. Since `a` is
invertible modulo `b`, for every residue `j` there is an `n` with

```text
an==j mod b.
```

Let `lambda` be the balanced representative of

```text
lambda==c a^(-1) mod b,              |lambda|<=b/2. (39)
```

Containment at those tooth centres would imply

```text
||j/b||<1/14       ->       ||lambda j/b||<=1/14.   (40)
```

Put `A=|lambda|`. If `A>=2`, take

```text
j=floor(b/(14A))+1.
```

Since `b>=50`,

```text
0<j<b/14,

b/14<Aj<=b/14+A<=4b/7.                              (41)
```

If `Aj<=b/2`, its centred residue is already greater than `b/14`; if
`Aj>b/2`, its centred residue is at least `b-4b/7=3b/7`. Either case
contradicts (40). Hence

```text
lambda in {0,+1,-1}.                                (42)
```

The case `lambda=0` gives `b|c`, contradicting `0<c<b`. The case
`lambda=+1` gives `c=a`, contradicting `7|c` and `7` not dividing
`a`. The remaining case `lambda=-1` gives

```text
b=a+c.                                              (43)
```

Write

```text
b=14m+r_0,                         1<=r_0<=13,
```

and choose `n` with `an==m mod b`. Since `c>=49`, one has `c>r_0`.
Put

```text
delta=(c-r_0)/(14b)
```

and choose

```text
0<eta<
 min(
   1/(7b),
   (1/14-delta)/a,
   delta/c
 ).
```

All three bounds are positive. At

```text
x=n/b-1/(14b)+eta
```

the phases, modulo integers, are

```text
b x=-1/14+b eta,

a x=delta+a eta,

c x=-(1/14+delta)+c eta.                            (44)
```

The first two have norm strictly below `1/14`, whereas the last has
norm strictly between `1/14` and `1/7`. Thus `x` belongs to both
source dangers and is strictly `c`-safe. All three inequalities are
strict, so a whole open neighbourhood of `x` lies in

```text
D_a intersection D_b minus closure(D_c).
```

In particular the containment fails on positive measure, not merely at
one exceptional endpoint. This is the final contradiction. QED.

Applying the lemma to (30)--(32) contradicts the quotient containment.
Therefore:

> **Two-top branch closure.** No canonical scalar cover realizes
> `k=2,(t,b)=(2,0)`.

The strict second descent of Section 5 is an independent shorter
exclusion on strict rows; the joint-comb lemma is what also reaches
the repeated-first equal-depth boundary.

## 7. Equal-depth repeated-first marginal hostiles

Neither marginal exclusion extends across both equalities in (29).
The following exact controls use quotient blockers

```text
(C_1,C_2,C_3)=(1,8,637)=(1,8,13*7^2).              (45)
```

The original blockers

```text
(c_1,c_2,c_3)=(13,104,8281)
```

have

```text
nu_13(c_1,c_2,c_3)=(1,1,2),

nu_7(c_1,c_2,c_3)=(0,0,2).                          (46)
```

Take `M=1`, `N=49`, and the generic fibre based at

```text
z=1/97.
```

Here “based at” means the additive orbit `z+j/N` of (9), not the
inverse-fibre convention `(z+j)/N`.

The two lower quotient words are

```text
D_(C_1): {0,1,2,45,46,47,48},

D_(C_2): {6,12,18,24,30,36,42},                    (47)
```

while the `C_3` word is empty. Thus the equal-depth intersection can
be zero on a high-safe fibre. The uniform `N/49` count in (19), and
hence the `6/343` marginal proof, genuinely fails.

For the second descent take instead

```text
y=2/15.
```

The thirteen-root words are

```text
C_1: {0,12},

C_2: {8},

C_3: empty.                                         (48)
```

The two unit words are nonempty but disjoint. In the repeated-first
case neither becomes a constant all-root danger word, so the inference
from (25) to (26) genuinely fails.

The joint-comb lemma nevertheless catches this same hostile packet
immediately. At

```text
x=1/1274=1/(2C_3),
```

both `C_1x` and `C_2x` have norm below `1/14`, while

```text
C_3x=1/2.
```

Thus its marginal words evade the two fibre arguments, but its actual
source intersection is maximally far from the target comb.

These are local/fibre hostiles, not scalar covers and not witnesses
against the joint-comb lemma. They identify the information missing
from both marginal arguments: the actual geometry of the joint
same-layer comb. Section 6 restores it through common zeros, tooth
centres, and one located boundary perturbation.

## 8. Relation to the collision gate and remaining scope

THM-2377 predicts that a nonzero high-target carrier needs a repeated
lower septimal layer. Equation (29) identifies exactly that apparent
two-top survivor, coupled to a repeated thirteen-adic first depth.
Section 6 shows that the scalar cover imposes more than the collision
gate: the entire intersection comb would have to fit inside one
two-layer-deeper target comb. The intrinsic binary datum is therefore
the located intersection word of `c_1,c_2`, not only its valuation,
support size, or Bockstein carry. Replacing it by the two marginal
support sizes loses precisely the distinction exhibited in (47)--(48).

This theorem:

- closes the entire `(2,0)` branch;
- retains the strict second descent and unequal-depth anti-shield as
  independent explanations of the two easier strata;
- removes no entire thirteen-adic row, because a strict row may realize
  another septimal alternative;
- leaves the remaining `(1,0)` septimal alternative untouched; and
- lands no target and does not prove LRC(14).

The scalar ledger therefore remains `165`.

## 9. Exact companion

The dependency-free exact companion:

- checks the lower per-bin counts and width-seven identity behind
  (10)--(11);
- exhausts the ordinary thirteen-root phase cells and verifies that
  two top words occupy at most four roots;
- checks unequal-depth intersections through three septimal layers;
- independently obtains `6/343` as the exact wall-cell measure of
  `D_1 intersect D_7 minus D_343`;
- checks the second descent's nonempty unit word and independently
  obtains `6/49` as the wall-cell measure of
  `D_1 minus D_637`;
- verifies the strict packet's `(1,2,3)` thirteen-adic profile; and
- reproduces both equal-depth repeated-first marginal hostiles
  (47)--(48);
- verifies the cyclic-subgroup escape for every order through `400`
  and the multiplier inequalities for every `50<=b<=400`; and
- constructs an exact source-intersection/target-safe witness for all
  `72,600` triples

  ```text
  1<=a<=b<=140,

  7 does not divide ab,

  c in {49,98,...,490},
  ```

  with positive controls for each of the subgroup, central-gap,
  multiplier, and boundary-perturbation mechanisms.

Every executable assertion uses an explicit optimization-safe
`require` function.

An independent audit reconstructed the a.e.-to-literal openness
upgrade, the gcd quotient, central-gap reduction, balanced-multiplier
classification, and final boundary perturbation. It also reran every
finite control in normal and optimized modes and found no missing
case or endpoint.

Run

```bash
python3 04-computation/lrc14_two_top_blocker_collision_reduction_thm2385.py
python3 -O 04-computation/lrc14_two_top_blocker_collision_reduction_thm2385.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_two_top_blocker_collision_reduction_thm2385.out
```

after LF normalization.
