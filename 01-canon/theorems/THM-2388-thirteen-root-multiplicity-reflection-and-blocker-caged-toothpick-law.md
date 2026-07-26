---
id: THM-2388
title: "Last septimal excess ledger and blocker-caged toothpick specialization"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. This
  specializes THM-2198/2222's already-proved thirteen-root -1
  eigenline and unique anti-defect root to the remaining
  k=2,(t,b)=(1,0) septimal excess. On
  X=D_(q_*)^c intersection D_(c_3)^c, the lower excess R is
  nonnegative, has exactly N/49 mass in every q_*-safe top bin, and
  has total integral 36/343. If G=(K-1)_+ and Z={K=0}, then
  36/343<=mu(Z intersection X)<=72/343 and an exact nonnegative
  decomposition spends the 36/343 current among unit overlap,
  simultaneous low-blocker ownership of a hole, and low-blocker
  incidence off the holes. Unit holes in X map into the divided-low
  blocker union outside the divided c_3 comb. This is a signed
  two-scale carrier, not a branch exclusion, row decrement, target
  landing, or proof of LRC(14).
source: codex-2026-07-26-thirteen-root-multiplicity-reflection
depends_on:
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
related:
  - THM-2234-first-depth-one-private-two-owner-mass-and-two-step-expansion
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
  - THM-2378-hard-septimal-root-stalk-closure
  - THM-2381-one-top-one-blocker-septimal-root-stalk-closure
  - THM-2382-saturated-septimal-seven-bin-root-fibre-closure
  - THM-2385-two-top-septimal-blocker-collision-reduction
script: 04-computation/lrc14_thirteen_root_multiplicity_reflection_thm2388.py
output: 05-knowledge/results/lrc14_thirteen_root_multiplicity_reflection_thm2388.out
script_sha256: 82f66c71d9cf941fcfa36f80ea1fae94c2171625942c100c8be8fe0ca5ee329e
output_sha256: c5c58af72968bcc7e1e2b65d36b2d68a1568aaa0534c7d8fcebd90958148661d
hash_basis: working-tree bytes (LF)
---

# THM-2388 -- the last septimal lane has an exact excess ledger

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Among THM-2367's only-`c_3`-dominant alternatives, the still-live lane
addressed here is

```text
k=2,                         (t,b)=(1,0).            (1)
```

Its lower load has width eight rather than seven, so replacing it by an
exact partition is false. The correct object is the signed multiplicity
current before owner deletion. It obeys an all-scale transfer law:

```text
unit hole at y outside the quotient blockers
  <-> one and only one double-covered inverse root,

one-fold point at y outside the quotient blockers
  <-> thirteen one-fold inverse roots.               (2)
```

The root reflection and unique anti-defect tooth are not new: they are the
sign-reversed form of THM-2198 Section 6 and the transfer-parity eigenline
of THM-2222 Sections 2--3. The new content is their exact coupling to the
septimal `36/343` waste ledger in (1). This inheritance matters: it prevents
an old automaton from being renamed as a new obstruction.

## 1. Recalling the proved unit-multiplicity reflection

On `T=R/Z`, put

```text
D_v={x:||vx||<1/14},

E_H={x:||Hx||<1/7},

K=1_(E_H)+sum_(i=1)^5 1_(D_(q_i)),                  (3)
```

where

```text
13 does not divide Hq_1...q_5.                       (4)
```

Let

```text
T(x)=13x
```

and use the unnormalised root transfer

```text
(P f)(y)=sum_(x:T(x)=y)f(x).                         (5)
```

For every generic `y`, a thirteen-unit ordinary mask has the exact root
count

```text
#{x:T(x)=y, x in D_q}=2-1_(D_q)(y),                  (6)
```

and the guard has

```text
#{x:T(x)=y, x in E_H}=4-1_(E_H)(y).                  (7)
```

To see (6), multiplication by `q` permutes the thirteen roots. On the
standard root grid, a centered arc of length `1/7` contains two points
unless the base phase itself lies in the centered `1/7` arc, when it
contains one. The identical count for an arc of length `2/7` is four
outside the guard arc and three inside it. Aligned endpoints form a finite
null set.

Summing (6)--(7) gives the pointwise reflection already encoded by
THM-2222:

```text
P K=14-K.                                             (8)
```

Since the six masks in (3) have total Haar mass

```text
2/7+5/7=1,
```

the centered multiplicity

```text
F=K-1                                                 (9)
```

is mean zero, integer valued, and satisfies

```text
P F=-F.                                              (10)
```

With Fourier convention

```text
fhat(n)=integral_T f(x) exp(-2 pi i n x)dx,
```

equation (10) is equivalently

```text
13 Fhat(13n)=-Fhat(n),                               (11)

13^r Fhat(13^r n)=(-1)^r Fhat(n).                   (12)
```

Thus the functional form is not an asymptotic drift or a finite-level
coincidence. It holds at every depth.

## 2. The prior scalar-cover cage in multiplicity notation

Suppose three blockers have the canonical form

```text
c_j=13C_j,                         j=1,2,3,           (13)
```

and the nine masks give a scalar cover:

```text
E_H union union_(i=1)^5 D_(q_i)
    union union_(j=1)^3 D_(c_j)=T                    (14)
```

almost everywhere. Put

```text
B=union_(j=1)^3 D_(C_j).                             (15)
```

The blocker bits are constant on a thirteen-root fibre:

```text
1_(D_(c_j))(x)=1_(D_(C_j))(T x).                    (16)
```

Two sign-support inclusions follow.

First, if `K(x)=0`, the scalar cover forces an original blocker at `x`.
By (16),

```text
{K=0} subset T^(-1)(B).                              (17)
```

Second, if `K(y)>=2`, then (10) gives

```text
sum_(T x=y)(K(x)-1)=1-K(y)<0.
```

Some inverse root has `K(x)=0`; (17) then gives `y=T x in B`. Hence

```text
{K>=2} subset B.                                     (18)
```

Equations (17)--(18) are a two-scale statement. Replacing both blocker
families by the same scalar support would lose the direction of `T`.

## 3. The prior toothpick law and its mass corollary

Let `y notin B` be generic. By (16), all three blockers are absent at
all thirteen inverse roots, so (14) gives

```text
K(x)>=1                         for every T x=y.      (19)
```

Now (8) has only two possibilities.

- If `K(y)=1`, the root sum is `13`; (19) forces

  ```text
  K(x)=1                         on all thirteen roots. (20)
  ```

- If `K(y)=0`, the root sum is `14`; (19) forces one root with
  multiplicity two and twelve roots with multiplicity one:

  ```text
  #{x:T x=y,K(x)=2}=1,

  #{x:T x=y,K(x)=1}=12.                              (21)
  ```

No `K(y)>=2` case occurs off `B`, consistently with (18).

Put

```text
S_0={y notin B:K(y)=0},

S_2={x:T x in S_0, K(x)=2}.                          (22)
```

Equation (21) makes

```text
T:S_2 -> S_0
```

a bijection up to null endpoints. Moreover, (18) puts

```text
S_2 subset B.                                        (23)
```

Haar disintegration and the unique root give the exact mass scaling

```text
measure(S_2)=measure(S_0)/13.                        (24)
```

The pointwise statement is THM-2198's one-anti-defect root after replacing
its residual `1-K` by `K-1`. The explicit mass identity (24) is a useful
corollary. This is not a two-to-one label quotient: one unit of missing
multiplicity becomes one unit of double multiplicity on a uniquely selected
root.

For a base inside `B`, the general layer balance retained by (10) is

```text
#{x:T x=y,K(x)=0}
 -sum_(x:T x=y)(K(x)-1)_+
 =K(y)-1.                                            (25)
```

This signed stalk, rather than the support of either side alone, is the
recursive datum.

## 4. The exact `36/343` current in the last septimal lane

Now assume the remaining alternative (1) in THM-2367's scope:

```text
M>0,                         nu_7(H)<M,

nu_7(q_*)=M,                nu_7(q_i)<M for q_i!=q_*,

nu_7(c_1),nu_7(c_2)<M<nu_7(c_3).                    (26)
```

Let

```text
mathcal A=D_(q_*) union D_(c_3),

L
 =1_(E_H)
  +sum_(q_i!=q_*)1_(D_(q_i))
  +1_(D_(c_1))+1_(D_(c_2)),

R=(L-1)1_(mathcal A^c).                              (27)
```

Write

```text
X=mathcal A^c,

b_i=1_(D_(c_i)),                         i=1,2,

F=K-1,                    G=F_+,         Z={K=0}.    (27a)
```

Because `q_*` is safe on `X`,

```text
R=(F+b_1+b_2)1_X>=0.                                (27b)
```

Put `N=7^(M+1)`. On every generic
`c_3`-safe `N`-orbit, `q_*` occupies one complete top bin. Every factor
in `L` is below depth `M`; its width times `N/49` is therefore its exact
incidence in each bin. The total lower width is

```text
2+4+1+1=8.                                          (28)
```

Consequently, in each of the six `q_*`-safe bins `B_r`,

```text
sum_(x in B_r)L(x)=8N/49,

sum_(x in B_r)R(x)=N/49.                             (29)
```

The first equality and coverage imply the second pointwise layer-cake law:

```text
sum_(j=2)^7
  #{x in B_r:L(x)>=j}
=N/49.                                              (30)
```

The `c_3`-safe orbit bases have mass `6/7`. Summing the six safe bins and
disintegrating gives the global exact current

```text
integral_T R
 =(6/7)(6/49)
 =36/343.                                            (31)
```

The same bin count, now before subtracting one, gives

```text
measure(X)=36/49=252/343,

integral_X K=216/343,

integral_X F=-36/343,                                (31a)

integral_X b_1=integral_X b_2=36/343.               (31b)
```

The last equality uses that either low blocker is subtop: it contributes
`N/49` points in each of the six safe bins on every `c_3`-safe orbit.

Since

```text
F=G-1_Z,
```

equation (31a) becomes the exact hole/overlap balance

```text
measure(Z intersection X)
 =36/343+integral_X G.                               (31c)
```

On a unit hole in `X`, the cover forces at least one low blocker. Splitting
the nonnegative current (27b) on and off `Z` yields the sharper ledger

```text
integral_X G
 +measure(Z intersection X intersection D_(c_1)
                          intersection D_(c_2))
 +integral_(X minus Z)(b_1+b_2)
 =36/343.                                            (31d)
```

Every term is nonnegative. Therefore

```text
0<=integral_X G<=36/343,

36/343<=measure(Z intersection X)<=72/343.          (31e)
```

There is also an exact directed support statement. If `x in Z intersection
X`, then the original cover pays the hole with `c_1` or `c_2`, while `c_3`
is safe. By (16),

```text
Z intersection X
 subset T^(-1)(
   (D_(C_1) union D_(C_2)) minus D_(C_3)
 ).                                                  (31f)
```

This is the genuinely new two-scale cage in the last lane. The old
unique-double rule controls fibres over `B^c`, whereas (31f) sends these
holes into `B`; confusing the two directions would create a false
recursion.

## 5. The thirteen-root form of the remaining excess

Write `c_j=13C_j` and define the rooted current globally by

```text
R_13=P R,

R_13(y)=sum_(T x=y)R(x).                            (32)
```

It vanishes when `y in D_(C_3)`, because then `c_3` is dangerous on
all thirteen inverse roots and `R` vanishes there. Now fix a generic
quotient phase `y` at which `C_3` is safe, and put

```text
x_h=(y+h)/13,                         h in F_13,

alpha=1_(D_(C_1))(y),       beta=1_(D_(C_2))(y),

Q(y)={h:x_h in D_(q_*)},              q=|Q(y)| in {1,2},

U_h
 =1_(E_H)(x_h)
  +sum_(q_i!=q_*)1_(D_(q_i))(x_h).                  (32a)
```

On every root outside `Q(y)`, both absorbers in `mathcal A` are absent,
so

```text
R(x_h)=U_h+alpha+beta-1.
```

Therefore, on `D_(C_3)^c`, the exact rooted excess is

```text
R_13(y)=sum_(h notin Q(y))R(x_h)

 =(13-q)(alpha+beta-1)
   +sum_(h notin Q(y))U_h
 >=0.                                               (33)
```

Set `R_13(y)=0` when `y in D_(C_3)`, where the rooted safe-sheet
description above is not being used. Thus every integral of `R_13`
below is over the extended function on the whole circle, equivalently
over `D_(C_3)^c`.

Its mean is fixed by (31):

```text
integral_T R_13(y)dy
 =13 integral_T R
 =468/343.                                           (34)
```

If both lower blockers are dangerous, every `q_*`-safe root is already
overflow, so the overflow word has size at least

```text
13-q>=11.                                           (35)
```

If exactly one is dangerous, overflow at a safe root is exactly the
presence of one of the five lower unit masks. In particular, the guard
word outside `Q(y)` is paid literally, not through an unsigned capacity
estimate.

Equations (30), (33), and (34) are the live interface: a small, uniformly
distributed septimal excess must realize a rooted thirteen-adic word and
the sign-reflecting toothpick law simultaneously.

There is a second exact coupling on the omitted top word. Since

```text
K=U+1_(D_(q_*)),
```

put

```text
O_(q_*)(y)=sum_(h in Q(y))U_h.                       (35a)
```

On a `q_*`-dangerous root one has `F=U`, so the prior eigenline gives

```text
sum_(h notin Q(y))F(x_h)
 =-F(y)-O_(q_*)(y).                                  (35b)
```

The seven-bin calculation on the top bin is as exact as the safe-bin
calculation: on every `c_3`-safe orbit, the lower unit load `U` has width
six and contributes `6N/49` incidences to the `q_*` bin. Therefore

```text
integral_(D_(q_*) intersection D_(c_3)^c) U
 =36/343,                                            (35c)

integral_(D_(C_3)^c) O_(q_*)(y)dy
 =468/343.                                           (35d)
```

Whenever the integrand in (35c) is positive, `q_*` and a lower unit mask
overlap, so (18) forces the physical point into the quotient-blocker cage:

```text
D_(q_*) intersection D_(c_3)^c intersection {U>0}
 subset B.                                           (35e)
```

For a parent `y notin B`, (20)--(21) make `O_(q_*)(y)` the zero-one
indicator that the unique double root of a unit hole lies in the
`q_*` word. Equations (35b)--(35e) are the explicit bridge between the
old thirteen-transfer automaton and the new septimal excess ledger.

## 6. Sharp local boundary and honest residual

Seven-bin capacity alone cannot close (1). At

```text
N=49,                    x_0=1/686,

H=1,

(q_*,q_2,q_3,q_4,q_5)=(7,148,171,172,243),

(c_1,c_2,c_3)=(195,169,215306),                      (36)
```

the roles are exactly

```text
M=1, k=2, (t,b)=(1,0),

(nu_13(c_1),nu_13(c_2),nu_13(c_3))=(1,2,3).
```

On the displayed `49`-orbit, the guard, `c_1`, and four lower unit masks
already form the exact tiling from THM-2367's hostile chamber. The new
`c_2=169` only adds coverage, `c_3` is absent, and `q_*=7` covers its top
bin. The mask word persists on a positive open chamber. It is not a global
scalar cover; it proves that a finish must use cross-orbit coherence of
`R`, the blocker cage, or the rooted word.

The natural binary relation now has vertices equal to unit-pair collision
obligations, not runners. Equation (18) sends every such collision into
the three-blocker union, but selecting a blocker colour discards the
root chosen in (21) and the sign balance (25). Any tournament/Fano
refinement must retain those sidecars; a cosmetic three-colouring is not an
equivalence.

This theorem does not close (1). THM-2385 independently closes the
two-top collision branch, leaving the present `(1,0)` lane as the sole
septimal alternative. No thirteen-adic row is removed here, the ledger
remains `165`, and LRC(14) remains open.

## 7. Exact companion

The dependency-free companion:

- checks (6)--(8) on every generic phase cell for a finite hostile bank
  of thirteen-unit speeds;
- verifies the Fourier/transfer normalization and the complete finite
  composition table behind (20)--(21);
- checks the unique-root mass factor `1/13`;
- verifies the septimal per-bin counts, layer-cake arithmetic, and
  constants `36/343`, `72/343`, and `468/343`;
- verifies the nonnegative decomposition (31d) symbolically and checks
  the exact `X`-integrals (31a)--(31b);
- checks the cross-overlap identities (35b)--(35d) and their exact
  `36/343`/`468/343` normalization;
- reproduces the rooted formula (33) over every blocker-bit and top-word
  size; and
- reconstructs the exact `49`-orbit hostile (36), including its
  valuation roles and positive endpoint margin.

Run

```bash
python3 04-computation/lrc14_thirteen_root_multiplicity_reflection_thm2388.py
python3 -O 04-computation/lrc14_thirteen_root_multiplicity_reflection_thm2388.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_thirteen_root_multiplicity_reflection_thm2388.out
```

after LF normalization. Every executable check raises explicitly under
optimized Python.

## 8. Independent audit

An independent read-only audit reconstructed the inherited
`P F=-F` sign and its Fourier factor `13`, all identities
(31a)--(31e), every directed cage in (17), (18), (31f), and (35e), and
the `36/343`, `72/343`, and `468/343` normalizations. It also rebuilt
the hostile masks and found an exact global hole with positive endpoint
margins, confirming that the local hostile is not a scalar cover.
Normal, optimized, and stored transcripts agree, and no executable
assertion disappears under optimization. QED.
