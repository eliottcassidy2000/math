---
id: THM-2586
title: "Depth-five arrival/deep/later-root diagonal on all 84 owner cells"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  On the canonical typed row, every one of the 12*7 nonzero-displacement/
  owner-clock cells in THM-2584's sigma={b}, K=2, r=5 packet has a positive
  theta-zero edge among (v,t)=(0,0),(6,12).  The two exact edge-zero sets are
  {(7,4),(7,5),(7,6)} and {(6,4),(6,5),(6,6)}, hence are disjoint.  In the
  rescaled deepest-root coordinate w=7t, theta=0 gives w=v.  Expanding the
  weighted transfer profiles and freezing a positive preimage branch turns
  each edge into one literal rational Boolean ancestry carrier.  Its route-two
  coordinate satisfies x=T^5 X exactly.  At every sufficiently large common
  delay, THM-2583 then inserts a named k_a boundary whose later physical root
  again equals v.  A positive danger half therefore has the genuine triple
  diagonal current arrival v = rescaled deepest root w = later physical root,
  while retaining the b word, collision displacement, owner clock, and deep
  root on one Boolean fibre product.  The fixed sheet is used only as a base
  carrier for a separately labelled future occurrence; it is not asserted
  neutral under the old whole-packet action.  This is not THM-2545's selected
  old-head diagonal, a THM-2334 current, or a covariant packet comparison.  The
  row is a typed non-cover control; no row is excluded and LRC(14) remains open.
source: lrc-semantic-frontier-2026-07-28
depends_on:
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2565-target-active-self-return-and-future-root-overlap-audit
  - THM-2583-self-similar-digit-needle-internalization-and-carrier-boundary
  - THM-2584-b-word-depth-five-absolute-deep-root-tensor
related:
  - THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile
  - THM-2569-stationary-diagonal-conditioned-paired-corner-and-frozen-future-role-boundary
script: 04-computation/lrc14_depth_five_arrival_future_diagonal_thm2586.py
output: 05-knowledge/results/lrc14_depth_five_arrival_future_diagonal_thm2586.out
script_sha256: 5bd3e3958ea5f5c27a2d61b5d9c05ea5369c6ac73ca73ac2ef62583b9b8bf6a0
output_sha256: 0e01ded47fb003978abbfe40a5fc1654a961b13f9fc7e56c95aff520113009b0
hash_basis: working-tree bytes (LF)
---

# THM-2586 -- every depth-five arrival cell reaches a later-role root

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2584 gives a positive support graph between the current arrival root and
the deepest-speed root, but stops before temporal continuation.  THM-2583
gives temporal digit continuation inside any positive rational carrier, but
its live application does not begin on a collision-root cell.  The two
objects compose without multiplying separate marginals:

```text
one theta-zero fine K cell
  -> one literal Boolean transfer/preimage sheet
  -> one future k_a boundary inside that same sheet
  -> arrival root = rescaled deepest root = genuinely later k_a root. (1)
```

This works for every nonzero displacement and owner clock in the packet.  The
conclusion is a real semantic diagonal at the **current-arrival** interface.
That qualifier is load-bearing: it is not yet the selected old head used by
THM-2545.

## 1. A theta-zero rail edge in all 84 cells

Retain THM-2584's canonical typed row, word, and packet:

```text
w=(1,14,27,40,53,66,13,13^3,2*13^5),

sigma={b},                  K=2,                  d=13^5,

f=1_Q P_(13^2)1_E,         U=P_d f,              V=P_d1_E.
                                                                    (2)
```

On its route-two physical coordinate `x`, the fine positive tensor is

```text
K_ell(s,v,t)
 =13 integral_T U(x)V(x-s/13)
      1_(floor(13x)=v)
      1_(floor(26x)=t mod 13)
      1_(frac(169x) in win_ell) dx.                 (3)
```

Here `s in F_13^*`, `ell in F_7`, and `v,t in F_13`.  Exact reconstruction
of all `13^3*7=15,379` entries gives the two zero sets

```text
{(s,ell):K_ell(s,0,0)=0}
  ={(7,4),(7,5),(7,6)},

{(s,ell):K_ell(s,6,12)=0}
  ={(6,4),(6,5),(6,6)}.                              (4)
```

They are disjoint.  Therefore the priority selector

```text
(v_(s,ell),t_(s,ell))
  =(0,0)  if K_ell(s,0,0)>0,
  =(6,12) otherwise                                  (5)
```

is defined for every one of the `12*7=84` pairs `(s,ell)` and obeys

```text
K_ell(s,v_(s,ell),t_(s,ell))
 >=13793820678/2120125746145771>0.                   (6)
```

It chooses `(0,0)` in `81` cells and `(6,12)` in the remaining three.  On
both edges, exactly in `F_13`,

```text
theta=t-2v=0,                  w=t/2=7t=v.           (7)
```

Thus (4), rather than a marginal sum, is the uniform all-cell diagonal-rail
certificate.

Fix an arbitrary `(s,ell) in F_13^* x F_7`, use (5), and abbreviate the
pointwise integrand in (3) by

```text
R_(s,ell)(x)
 =U(x)V(x-s/13)
  1_(floor(13x)=v)
  1_(floor(26x)=t mod 13)
  1_(frac(169x) in win_ell).                         (8)
```

It is a nonnegative rational step function and

```text
13 integral_T R_(s,ell)(x) dx=K_ell(s,v,t)>0.        (9)
```

## 2. Freeze one literal Boolean ancestry sheet

The positivity in (9) is not left as an averaged transfer statement.  Put

```text
R_0=13^2,                         x_s={x-s/13},

X_a(x)=(x+a)/d,                   0<=a<d,

Z_(a,beta)(x)=(X_a(x)+beta)/R_0,  0<=beta<R_0,

Y_e(x)=(x_s+e)/d,                 0<=e<d.             (10)
```

The braces in `x_s` are essential: all Perron preimages use the canonical
circle representative.  Expanding the two Perron averages gives exactly

```text
U(x)
 =1/(d R_0) sum_(a,beta)
    1_Q(X_a(x))1_E(Z_(a,beta)(x)),

V(x-s/13)
 =1/d sum_e 1_E(Y_e(x)).                            (11)
```

For each triple `(a,beta,e)`, let `S_(a,beta,e)` be the Boolean set on which
the three indicators in (11) and all three root/owner indicators in (8)
equal one.  Then

```text
R_(s,ell)
 =1/(R_0 d^2) sum_(a,beta,e)1_(S_(a,beta,e)).        (12)
```

All sets are rational finite unions of intervals.  Equations (9)--(12)
force at least one triple with

```text
measure(S_(a,beta,e))>0.                              (13)
```

Freeze such a triple and write `S=S_(a,beta,e)`.  It is a literal positive
Boolean fibre-product sheet.  Every `x in S`, away from null walls, retains
simultaneously

```text
X=X_a(x) in Q=E_b,                 the pure {b} arrival word;

Z=Z_(a,beta)(x) in E,              prescribed-clock source ancestry;

Y=Y_e(x) in E,                     source collision branch;

s=v-u !=0,                         first-collision displacement;

ell,                               native owner-clock cell;

v=floor(13x) in {0,6},

t=floor(26x) mod 13 in {0,12},

theta=t-2v=0,                      w=7t=v.             (14)
```

No word, root, or owner fact in (14) was recovered from a separately
integrated marginal.

## 3. The route-two coordinate is original time five

The fixed arrival branch in (10) satisfies

```text
dX=x+a.
```

Since `d=13^5`, this is the exact circle identity

```text
T^5 X=x,                    T(y)=13y mod 1.             (15)
```

Consequently, for every positive integer `N` and every physical speed `k`,

```text
k T^N x=k T^(N+5)X                         modulo 1.    (16)
```

The chronology is therefore:

| event | coordinate | time from `X` | retained root |
|:--|:--|--:|:--|
| fixed depth-five arrival | `x=T^5X` | `5` | `v=floor(13x)` |
| separately labelled future occurrence | `T^Nx` | `N+5` | `floor(13T^Nx)` |

Equation (16) does not identify two abstract `F_13` alphabets.  It identifies
two evaluations of the same physical arrival point `X`.

## 4. Insert a genuinely later target-active occurrence

The positive rational carrier `S` lies inside the open root digit
`I_v=(v/13,(v+1)/13)` up to null walls.  Apply the geometric piercing and
cylinder parts of THM-2583 to `S`, with

```text
k=k_a,                       L=L_a,                    (17)
```

where `k_a` is put first in THM-2461's freely ordered future bank.  For every
sufficiently large `N`, there are a labelled future boundary `x_N` in the
interior of `S` and a physical base-thirteen cylinder `H_N subset S` such
that

```text
floor(13x)=floor(13T^N x)=v

                         for every x in H_N,           (18)
```

and `H_N` meets the complete target boundary atlas of

```text
P_q^(N)(x)=d_(L_a)(k_a T^N x+q/13)                    (19)
```

at the single point `x_N`, with one label `q=q_0`.  Both sets

```text
H_N P_(q_0)^(N),             H_N(1-P_(q_0)^(N))       (20)
```

have positive measure.  Choose the danger half in (20).  Putting `k_a`
first makes this a literal `k_a`-labelled first-failure occurrence, hence a
categorically target-active role by THM-2565 Section 3.  By (15)--(19), its
physical root at original time `N+5` is

```text
floor(13T^(N+5)X)
 =floor(13T^N x)
 =v
 =w.                                                   (21)
```

All facts in (14) remain true because the danger half is contained in
`H_N subset S`.  Thus (21) is one positive common-ancestry triple diagonal,
not a coupling inferred from equal marginals.  Since only 84 carriers were
selected and THM-2583 works at every sufficiently large delay, one maximum
of their finite piercing and chronology thresholds gives a single common
`N` if desired.

There are two different actions here and they must not be conflated:

| action | what is fixed | lawful conclusion |
|:--|:--|:--|
| old whole-packet action | nothing is asserted neutral; `S` may contain old moving factors | only the base facts (14) on the frozen sheet |
| new future-occurrence label `q` in (19) | `S` was chosen before and independently of `q` | one labelled boundary and one positive physical danger occurrence |

Accordingly, no target colours or shifted copies of `S` are combined with
the old packet.  This is exactly the separation needed to avoid
MISTAKE-266: the theorem uses one physical danger side, not a frozen old
factor as a covariant target current.

## 5. What this closes and what it does not

For every `(s,ell) in F_13^* x F_7`, the proved chain is

```text
{b}, r=5 first-collision service
  + a theta-zero arrival/deep rail with w=v
  + one literal transfer/preimage branch
  + a sufficiently later named k_a boundary
  -> positive same-ancestry equality
       current arrival root
         = rescaled deepest root
         = later target-active physical root.                       (22)
```

This closes the temporal-sheet question left by THM-2471 on all 84 cells of
the canonical typed packet.  Its route-two arrival root is not permanently
trapped at the current collision time: it can be continued to a literal
later occurrence without losing the selected word, nonzero displacement,
deepest root, or owner clock.

It does **not** close THM-2545.  That theorem's left vertex is the selected
empty old-head root produced by THM-2537, while the left vertex in (21) is
the current word-restricted arrival root `v`.  No map between them is proved.

Likewise:

- the future boundary label `q` is not THM-2365's lawful full-packet
  co-shift;
- neither `v` nor the later digit is THM-2334's relation residue `eta.u`;
- the digit cylinder does not select a physical Fourier frequency `X` or
  survive THM-2568's full-frequency pushforward;
- the fixed sheet is not promoted to a target current or an old-action
  neutral packet; and
- the construction is on one canonical typed non-cover row, not all `165`
  live scalar rows.

The next map is narrower: attach THM-2537's selected old head, or THM-2569's
old-head/future packet, to the current arrival sheet `S` before applying
(21).  A second root census is unnecessary.

The independent hostile audit rederived all 84 quantifiers and the exact
floor, checked `K=13 integral R`, the canonical circle Perron expansion,
the `v,w,t` typing, finite common-delay maximum, `T^5X=x` chronology,
`k_a`-first categorical role, and the old/new-action boundary.  Normal and
optimized runs byte-match the stored output; hashes and bytecode compilation
pass.

No scalar row is excluded and LRC(14) remains open. **QED.**
