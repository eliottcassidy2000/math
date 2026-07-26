---
id: THM-2426
title: "Compositional thirteen-root exclusion of the primitive final septimal lane"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Under
  THM-2391's primitive final-lane scalar-cover hypotheses, the remaining
  adjacent top depth M=1 is impossible. Its thirteen-root section has
  mass 6/637 in the guard-collision phase, while a sharp
  narrow-comb/phase overlap lemma loses at most 1/182, leaving 5/1274.
  THM-2391's physical blocker partition then doubles both guard roots,
  contradicting its W8 unique-double word. Together with THM-2417's
  independently audited M=2 orbit dual, the entire primitive final
  septimal lane is empty. A second compositional proof of M=2 lands
  mass at least 36/4459. With THM-2367/2378/2381/2382/2385 this proves
  the deep-blocker/guard dichotomy nu_7(c_3)>M => nu_7(H)=M. Full
  owner-mask cancellation remains open; no scalar profile row is
  removed and LRC(14) is not proved.
source: codex-2026-07-26-final-lane-compositional-root
depends_on:
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
  - THM-2378-hard-septimal-root-stalk-closure
  - THM-2381-one-top-one-blocker-septimal-root-stalk-closure
  - THM-2382-saturated-septimal-seven-bin-root-fibre-closure
  - THM-2385-two-top-septimal-blocker-collision-reduction
  - THM-2391-blocker-caged-septimal-single-layer-address-reduction
  - THM-2415-last-lane-septimal-depth-two-cap
  - THM-2417-middle-depth-two-relaxed-orbit-dual-obstruction
related:
  - THM-2390-septimal-layer-kraft-peeling-and-heavy-word-reduction
  - THM-2405-two-level-septimal-sheet-independence-and-middle-depth-two-cage-elimination
  - THM-2414-thirteen-skew-septimal-word-transport-and-local-stopping-atlas
script: 04-computation/lrc14_final_lane_compositional_root_exclusion_thm2426.py
output: 05-knowledge/results/lrc14_final_lane_compositional_root_exclusion_thm2426.out
script_sha256: 5b502bcdd0c4908114f9c81f85b8c1192eea5cb392aa535bad1710e19d6078ef
output_sha256: 449cae7e426f65796cb82692f126723d7fc16b4954c37d6233a60da90dfcf255
hash_basis: working-tree bytes (LF)
cite_by_filename: true
---

# THM-2426 -- the primitive final septimal lane is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2415 proves uniformly that the top septimal depth in THM-2391's
primitive final lane satisfies `M<=2`. THM-2417 independently excludes
`M=2` by a complete `343`-orbit relaxation and weighted dual. The
remaining canonical task is therefore `M=1`.

The same compositional mechanism also supplies a shorter second proof
at `M=2`. THM-2415's exact hostile shows that one isolated seven-point
phase orbit can miss the required guard chamber. That is a boundary of
the one-scale orbit argument, not of the scalar packet.

Keeping the factor thirteen in the high blocker supplies a transverse
root cospan:

```text
M=2:
  divide by seven
  -> expose thirteen roots
  -> at least six roots survive both transition and phase
  -> lift the unique good phase root back through seven;

M=1:
  expose thirteen roots directly
  -> one exact phase root per good base
  -> subtract at most half of the narrow central comb.               (1)
```

Both constructions leave positive measure in the phase chamber where
the guard words at speeds `H` and `13H` coincide. THM-2391 then makes
the two physical blockers partition those same two guard roots. Both
roots double, contradicting the unique-double W8 word.

## 1. Inherited hypotheses and the common contradiction

Retain all hypotheses and notation of THM-2391 and THM-2415. Thus
there is an almost-everywhere scalar cover

```text
C_H subset
  union_(i=1)^5 D_(q_i)
  union D_(c_1) union D_(c_2) union D_(c_3),                       (2)

D_v={x in R/Z: ||vx||<1/14},
E_v={x in R/Z: ||vx||<1/7},
```

the blocker labels obey

```text
c_j=13C_j,                                                         (3)
```

and one top ordinary label has

```text
nu_7(q_*)=M>=1,                 nu_7(c_3)>M.                        (4)
```

Primitivity and THM-2391 give

```text
nu_7(H)=nu_7(q_i)=nu_7(c_1)=nu_7(c_2)=0,
                                      q_i!=q_*.                    (5)
```

In particular,

```text
7 does not divide H,             13 does not divide H*q_*.         (6)
```

Put

```text
Q=q_*/7,                         D=c_3/7,                           (7)

S=(D_(13Q) minus D_Q) intersection D_D^c,                          (8)

I_H={y: 6/13<{Hy}<7/13}.                                           (9)
```

The divided-base identities in THM-2415 say that every `y in S`
makes `q_*` and `c_3` safe on all seven roots

```text
x_j=(y+j)/7
```

while `13q_*` is dangerous on all seven.

We first isolate the conclusion already proved in THM-2415 Sections
4--5:

> **Phase-collision consequence.** If `mu(S intersection I_H)>0`,
> then the packet is impossible.

Indeed, delete the finitely many comb endpoints and the seven-root
pullbacks of the two inherited null exceptional sets. Choose a
remaining `y`. In the `H`-gauge the guard roots are

```text
W_H={0,6}.                                                         (10)
```

For `r={Hy}` and `b=floor(13r)`, the affine thirteen-drift word is

```text
W_(13H)={b,b+1}.                                                   (11)
```

On (9), `b=6`, so

```text
W_(13H)=W_H.                                                       (12)
```

The lower layer consists of the width-two guard, four singleton lower
ordinary words, and two singleton physical blocker words. Its total
incidence is eight and it covers seven roots, so THM-2391's W8 theorem
allows exactly one doubled root.

On the other hand, the physical pullback identity holds at every root:

```text
1_(D_(c_1))+1_(D_(c_2))=1_(E_(13H))

on D_(13q_*) intersection D_(c_3)^c.                              (13)
```

Equations (12)--(13) make the blockers partition the same two roots
already occupied by the guard. Both distinct guard roots therefore
have multiplicity at least two, a contradiction.

It remains only to prove positive phase mass at both possible depths.

## 2. A reusable different-depth safe section

Let

```text
7 does not divide u*v,                    g>=1,

E(u,v,g)=D_u intersection D_(7^g v)^c.                              (14)
```

Then

```text
mu(E(u,v,g))=6/49.                                                  (15)
```

To see this, disintegrate over the seven roots

```text
t_j=(s+j)/7.
```

The high mask is constant on each root fibre:

```text
1_(D_(7^g v))(t_j)=1_(D_(7^(g-1)v))(s).                            (16)
```

Away from endpoints, the seven-unit mask `D_u` occupies exactly one
root. Hence

```text
mu(D_u intersection D_(7^g v))
 =(1/7)mu(D_(7^(g-1)v))
 =1/49.                                                            (17)
```

Since `mu(D_u)=1/7`, equation (15) follows.

This is the operation-graph content of the proof: the seven-root
section composes with a thirteen-root section before any information
is collapsed to one scalar measure.

## 3. Excluding `M=2`

Assume

```text
M=2.
```

Write

```text
q_*=49u,                  7 does not divide u, 13 does not divide u,

c_3=13*7^L v,             L>=3, 7 does not divide v.                 (18)
```

Then

```text
Q=7u,                     D=13*7^(L-1)v.
```

Set

```text
t=7y,                     g=L-2>=1.                                (19)
```

The transition set (8) is exactly the pullback under multiplication
by seven of

```text
S_0
 =(D_(13u) minus D_u)
   intersection D_(13*7^g v)^c:

S=T_7^(-1)(S_0).                                                     (20)
```

Define the closed phase band

```text
F={t: ||Ht||<=3/13}.                                                (21)
```

Take a generic

```text
s in E(u,v,g).
```

On the thirteen roots

```text
t_j=(s+j)/13,                      j in F_13,                        (22)
```

the following exact counts hold:

1. all thirteen roots lie in `D_(13u)`, because `s in D_u`;
2. all thirteen lie outside `D_(13*7^g v)`, because
   `s notin D_(7^g v)`;
3. exactly one root lies in `D_u`; and
4. exactly seven roots lie in `F^c`.

For item 3, multiplication by the 13-unit `u` permutes the roots. In
the balanced lift `epsilon=us` with `|epsilon|<1/14`, the central root
is dangerous while every nonzero residue root is at distance at least

```text
(1-|epsilon|)/13>=1/14,                                          (22a)
```

with equality only at a deleted endpoint. For item 4, the 13-unit `H`
permutes the roots, while the generic closed arc `F` has length `6/13`.
These root counts, and hence the overlap floor below, are
almost-everywhere statements: the finitely many endpoint bases have
already been removed.

At least

```text
13-1+7-13=6                                                     (23)
```

roots therefore lie in `S_0 intersection F^c`. Root disintegration
and (15) give

 ```text
mu(S_0 intersection F^c)
 >=(6/13)(6/49)
 =36/637.                                                           (24)
```

Multiplication by seven maps the chamber `I_H` bijectively, in the
`H`-phase gauge, onto `F^c`:

```text
r in (6/13,7/13)
  -> {7r} in (3/13,10/13).                                         (25)
```

Because `H` is a seven-unit, every `t in F^c` has exactly one of its
seven preimages in `I_H`. Combining (20), (24), and (25),

```text
mu(S intersection I_H)
 =mu(S_0 intersection F^c)/7
 >=36/4459
 >0.                                                               (26)
```

The phase-collision consequence excludes `M=2`.

This repairs rather than contradicts THM-2415's hostile
`Q=7,D=49,H=1,y=1/91`: that one seven-point orbit misses `I_H`, but
the transverse thirteen-root section proves that the full transition
set has other orbits with positive landed phase mass.

## 4. A narrow-comb half-overlap lemma

The adjacent depth needs one additional estimate.

> **Half-overlap lemma.** If `13` divides neither `H` nor `u`, then
>
> ```text
> mu(D_u intersection D_(13u) intersection I_H)<=1/182.             (27)
> ```

First observe, up to endpoints, that

```text
D_u intersection D_(13u)
 ={x: ||ux||<1/182}.                                                (28)
```

Indeed, inside `||ux||<1/14`, the neighboring components of the
`13u` comb begin exactly at the excluded endpoints `+/-1/14`, leaving
only its central interval.

Put

```text
g_0=gcd(H,u),             H=g_0 a,             u=g_0 b,

gcd(a,b)=1.                                                         (29)
```

Haar pushforward by `g_0` reduces (27) to

```text
mu({||bz||<1/182} intersection I_a)<=1/182.                         (30)
```

If `b>=2`, disintegrate the narrow `b`-comb using

```text
z=(x+j)/b,                    |x|<1/182.                            (31)
```

For each fixed `x`, coprimality makes the `a`-phases of the `b` roots
a translated `b`-grid. An open interval of length `1/13` contains at
most

```text
ceil(b/13)
```

grid points. Since `13` does not divide `b` and `b>=2`,

```text
ceil(b/13)/b<=1/2.                                                  (32)
```

Integrating (31) gives (30).

If `b=1`, put

```text
P=union_(k in Z)(k+6/13,k+7/13).                                   (33)
```

On the positive half-line, translation by `-6/13` injects

```text
P intersection (0,A)
```

measure-preservingly into its complement in `(0,A)`: each component
maps to `(k,k+1/13)`. Hence `P` occupies at most half of `(0,A)`.
Since `P=-P`, it occupies at most half of every centered interval
`(-A,A)`. Apply this after `t=az` to the centered narrow interval
`|z|<1/182`; equation (30) follows.

The constant is sharp at the grid/Farey boundary `b=2`. For example,
odd `a in {1,3,5,7,9,11}` place one whole narrow component in the
half-bin.

## 5. Excluding `M=1`

Assume now

```text
M=1.
```

Write

```text
q_*=7u,                   7 does not divide u, 13 does not divide u,

c_3=13*7^L v,             L>=2, 7 does not divide v.                 (34)
```

Now the transition (8) is already

```text
S
 =(D_(13u) minus D_u)
   intersection D_(13*7^(L-1)v)^c.                                 (35)
```

Let

```text
E=E(u,v,L-1).
```

On each generic thirteen-root fibre

```text
y_j=(s+j)/13,                    s in E,                            (36)
```

all thirteen roots lie in `D_(13u)`, all are outside
`D_(13*7^(L-1)v)`, and exactly one lies in `I_H`. The last statement
uses only that `H` is a 13-unit: an interval of length `1/13` contains
one generic point of a thirteen-grid.

Consequently the section

```text
J=D_(13u) intersection D_(13*7^(L-1)v)^c intersection I_H          (37)
```

has exact mass

```text
mu(J)=mu(E)/13=6/637.                                               (38)
```

The only roots in `J` not belonging to the transition (35) lie in
`D_u`. The half-overlap lemma therefore gives

```text
mu(S intersection I_H)
 >=6/637-1/182
 =5/1274
 >0.                                                               (39)
```

The phase-collision consequence excludes `M=1`.

The bound is attained by the typed arithmetic control

```text
(H,u,v,L)=(5,2,1,2),

q_*=14,                         c_3=637.                            (40)
```

This is a sharp transition/phase control, not a scalar-cover packet.

## 6. Conclusion and scope

THM-2415 leaves only

```text
M in {1,2}.
```

THM-2417 excludes `M=2`, and Section 5 excludes the remaining `M=1`.
Section 3 independently re-proves the `M=2` exclusion by a transverse
root argument. Therefore:

> **Primitive final-lane exclusion.** No packet satisfies all of
> THM-2391's primitive final-lane scalar-cover hypotheses.

The factor `13 does not divide H` is load-bearing. If `H=13,u=v=g=1`,
there are generic bases in (14) for which the six-useful-root count in
Section 3 collapses to zero. This is an exact hostile to deleting the
inherited 13-unit condition, not a hostile within the theorem.

The constants are one instance of a compositional root family rather
than isolated `7,13` numerology. Let `r>=3` and `p=2r-1` be primes,
let `g>=1`, and take reduced labels satisfying

```text
r does not divide u*v*H,              p does not divide u*H.
```

Replace the danger radius by `1/(2r)`, use the safe section

```text
E=D_u intersection D_(r^g v)^c
```

and its `p`-root high mask `D_(p r^g v)`. Replace (9) by

```text
I_(H,r)={(r-1)/p<{Hx}<r/p},
```

The same two disintegrations give

```text
different-r-depth safe mass                  (r-1)/r^2,

depth-one p-root section                     (r-1)/(p r^2),

reflection-fixed central-comb loss at most   1/(2rp),

depth-one landed floor                       (r-2)/(2p r^2),

depth-two landed floor                       (r-1)^2/(p r^3).        (41)
```

The identity `p=2r-1` is what makes the neighboring `D_(pu)` teeth
begin exactly at the boundary of `D_u`; all overlap asymmetry is
therefore concentrated in the reflection-fixed central tooth. The
reduced denominator `2` is the sharp half-overlap boundary. At depth
two the order of operations is forced: the `r`-adic depth first makes
the high mask constant on `r`-roots, and its mandatory `p` factor then
makes it constant on `p`-roots. Both landed floors in (41) are positive,
and they specialize to `5/1274` and `36/4459` at `(r,p)=(7,13)`. It is
a root-cospan mechanism, not a tournament equivalence.

The theorem closes the primitive `k=2,(t,b)=(1,0)` final septimal lane,
including both positive-clean and no-clean packets. It does **not**
remove one of the `165` scalar profile rows: the lane is a downstream
structural alternative rather than a complete profile classifier. It
does not close other LRC routes or prove LRC(14).

## 7. Proof-graph corollary: deep blocker forces a top guard

There is a stronger consequence only after composing this theorem with
the earlier branch closures. Retain a primitive scalar cover and put

```text
M=max(nu_7(H),nu_7(q_1),...,nu_7(q_5)).
```

Then the earlier closures and this theorem imply the clean dichotomy

```text
nu_7(c_3)>M          implies          nu_7(H)=M.                    (42)
```

Indeed, suppose instead that

```text
nu_7(c_3)>M,                    nu_7(H)<M.
```

The second inequality forces `M>0`. THM-2367's complete seven-bin
classification gives:

```text
k=0                                      impossible in THM-2367,

k=1, (t,b)=(1,0)                         empty by THM-2378,

k=2, (t,b)=(1,1)                         empty by THM-2381,

k=2, (t,b)=(2,0)                         empty by THM-2385,

k=2, (t,b)=(5,2)                         empty by THM-2382,

k=2, (t,b)=(1,0)                         empty by this theorem.    (43)
```

Here `k,t,b` are exactly THM-2367's low-blocker, top-ordinary, and
top-blocker counts. The list is exhaustive, proving (42).

Equivalently, every subtop-guard cover satisfies

```text
nu_7(H)<M          implies          nu_7(c_3)<=M.                   (44)
```

Choose an ordinary label `q_*` with `nu_7(q_*)=M`. Since `c_3` is the
eligible deepest thirteen-adic target, (44) puts its isolated lawful
`(c_3,q_*)` graft on the noncirculant side of THM-2367 Section 1. Thus
the old arithmetic role mismatch is gone: in the subtop-guard regime,
the deepest target itself supplies an escaping septimal graft.

This still does not prove drift of the complete bare-owner tensor.
THM-2367's exact `90/91`-mass mask shows that additional owner factors
can restore circulancy. The next obligation is therefore a
source-conditioned no-cancellation theorem for this now-canonical
`(c_3,q_*)` graft, not another search over septimal valuation
alternatives.

## 8. Exact companion

Run:

```text
python3 04-computation/lrc14_final_lane_compositional_root_exclusion_thm2426.py
python3 -O 04-computation/lrc14_final_lane_compositional_root_exclusion_thm2426.py
```

The standard-library `Fraction` companion:

- replays the different-depth mass `6/49`;
- checks the `M=2` thirteen-root floor and exact seven-root landing;
- verifies the `M=1` section mass, overlap ceiling, and sharp
  `5/1274` phase floor on typed controls;
- independently integrates `3367` coprime narrow-comb/phase cases;
- checks `624` generic thirteen-root bases directly;
- verifies the sharp half-overlap value; and
- retains the `13|H` zero-useful-root hostile.

Normal and optimized modes reproduce:

```text
05-knowledge/results/lrc14_final_lane_compositional_root_exclusion_thm2426.out
```

## 9. Lean kernel

The root-imported module

```text
04-computation/lean/TournamentH7/TournamentH7/
  LRCFinalLaneRootCospan.lean
```

formalizes:

- unit-affine permutations of `ZMod 13`;
- the sharp finite `12 intersection 7 >= 6` cospan count;
- the rational identities `36/637`, `36/4459`, and `5/1274>0`; and
- normalized Haar invariance under nonzero integer dilation, including
  dilation by seven.

It builds with

```text
lake build TournamentH7.LRCFinalLaneRootCospan
```

and its printed axiom audit contains only `propext`, `Classical.choice`,
and `Quot.sound`, with no `sorryAx` or `native_decide`.

The module deliberately does not formalize the comb-to-fibre
identifications, generic endpoint deletion, conditional disintegration,
or the THM-2391 scalar-cover contradiction. It is a compile-checked
kernel of this theorem, not a formal proof of the whole lane exclusion.
