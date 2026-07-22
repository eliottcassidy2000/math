---
id: THM-2137
title: "Deep scalar tails force unit-comb boundary complexity"
status: >
  PROVED. If a residual guard/unit set with boundary budget
  B=H+sum(q_i) is covered by danger combs whose coefficients share a factor
  M, then its normalized M-root-count transfer has variation at most 2B/M
  and support inside the divided target combs. Consequently
  B>=M*measure(L)/measure(T). THM-2091's distinct odd-guard charge spectrum
  gives measure(L)>=191/6930 for the six-plus-one scalar tail and
  measure(L)>=961/6930 for the five-plus-two tail. Thus the respective sharp
  consequences of this method are B>=(191/990)13^a and
  B>=(12493/24750)13^a, where a is the common deep 13-adic depth. In
  particular every surviving six-plus-one row has B>=71634. These are
  all-depth height-versus-depth reductions, not by themselves LRC(14).
source: codex-2026-07-22-LRC-deep-fibre-boundary-complexity
depends_on:
  - THM-1166
  - THM-2091
  - THM-2133
  - THM-2135
related:
  - THM-2080
  - THM-2128
---

# THM-2137 -- deep scalar tails force boundary complexity

Put

```text
C_H={t:||Ht||>1/7},       E_H={t:||Ht||<1/7},
D_r={t:||rt||<1/14}.                                  (1)
```

The scalar tails in THM-2135 were first studied one thirteen-root needle at
a time.  The more faithful all-depth object is the complete root fibre of the
largest common thirteen-power in the deep coefficients.  Its multiplicity
function retains both the support of the quotient and the boundary variation
lost by merely projecting the residual set.

## 1. The deep-root boundary-complexity lemma

Let `q_1,...,q_k` be distinct positive integers and put

```text
L=C_H minus union_(i=1)^k D_(q_i),
B=H+sum_(i=1)^k q_i.                                  (2)
```

Suppose that, outside a null set, one has

```text
L subset union_(s=1)^ell D_(M w_s),                  (3)
```

where `M` and the `w_s` are positive integers.  Define the divided target

```text
T=union_(s=1)^ell D_(w_s),             sigma=measure(T), (4)
```

and assume `sigma<1`.  Then

```text
B >= M measure(L),
B >= M measure(L)/sigma.                              (5)
```

The second inequality is the stronger one; the first has a useful direct
empty-fibre proof.

### Proof of the empty-fibre bound

Let `Z` be a null exceptional set outside which (3) holds.  The boundary of
`L` is contained in the guard and unit-comb boundary set, which has at most

```text
2H+2sum_i q_i=2B                                      (6)
```

points.  Modulo those endpoints, `L` is therefore a disjoint union of at
most `B` circular intervals.

Multiplication `[M]:t |-> Mt` sends `Z` to a null set and the finite boundary
of `L` to a finite set.  Since the closure of `T` has measure `sigma<1`, one
may choose

```text
alpha notin closure(T) union [M](Z) union [M](boundary(L)). (7)
```

If `Mx=alpha`, then every deep value is

```text
(M w_s)x=w_s alpha,                                  (8)
```

so `x` is strictly safe for every deep tooth.  It is also nonexceptional.
Thus no such root lies in `L`.

The `M` roots of `alpha` form a shifted grid.  On one circular interval the
number of grid points differs from `M` times its length by at most one.
Summing over the at most `B` components of `L`, and using the zero root count
just obtained, gives

```text
M measure(L)<=B.                                      (9)
```

This proves the first inequality in (5).  Removing `[M](Z)` in (7) is
essential: a null exceptional set can contain one complete finite root
fibre.

### Proof of the support/variation bound

Choose any endpoint convention for `1_L` and define the normalized transfer

```text
F_M(alpha)=(1/M) sum_(Mx=alpha) 1_L(x).               (10)
```

Haar invariance gives

```text
integral F_M=measure(L).                              (11)
```

Outside `[M](Z)`, equation (8) and the cover imply

```text
essential_support(F_M) subset closure(T).             (12)
```

Every jump of `1_L` contributes a jump of size `1/M` after transfer.
Different boundary events may collide, but the triangle inequality still
gives

```text
Var(F_M)<=Var(1_L)/M<=2B/M.                           (13)
```

The nonnegative function `F_M` is supported on a set of measure at most
`sigma`, so (11) gives

```text
essential_sup F_M>=measure(L)/sigma.                  (14)
```

It vanishes on an open interval outside the finite union `closure(T)`.  A
periodic function that rises from zero to height `A` and returns to zero has
total variation at least `2A`.  Equations (13)--(14) therefore yield

```text
2B/M>=2 measure(L)/sigma,                             (15)
```

which proves the second inequality in (5).  This argument is insensitive to
all endpoint conventions because they change only finite sets.  QED.

## 2. The two exact residual-mass constants

For an odd guard `H`, write the THM-2080 overlap as

```text
I(H,q)=measure(E_H intersection D_q)=2/49+e(H,q).     (16)
```

For fixed `H` and distinct `q`, THM-2091 proves that the seven largest
possible deficits `-e(H,q)` are exactly

```text
5/294, 8/539, 3/245, 3/245, 4/441, 4/441, 9/1078.   (17)
```

Hunter's star centred at `E_H` gives

```text
measure(L)
 >=1-(k+2)/7+sum_(i=1)^k I(H,q_i).                   (18)
```

For five unit combs, the five worst distinct charges in (17) give

```text
measure(L)>=delta_5
 =10/49-(5/294+8/539+3/245+3/245+4/441)
 =961/6930.                                           (19)
```

For six unit combs, one also pays the excess baseline `1/7`, and the first
six charges give

```text
measure(L)>=delta_6
 =12/49-(5/294+8/539+3/245+3/245+4/441+4/441)-1/7
 =191/6930.                                           (20)
```

Distinctness is decisive here.  The sharp single-pair overlap `1/42` cannot
be repeated five or six times against one fixed odd guard.

## 3. Application to the five-plus-two tail

In case (II) of THM-2135, put

```text
a=min(nu_13(13v_1),nu_13(13v_2)),       M=13^a,
13v_s=Mw_s.                                           (21)
```

The divided coefficients `w_1,w_2` remain distinct.  By THM-1166 their two
danger combs have intersection measure at least `1/91`; hence

```text
sigma=measure(D_(w_1) union D_(w_2))
 <=2/7-1/91=25/91.                                   (22)
```

Equations (5), (19), and (22) prove

```text
H+sum_(i=1)^5 q_i
 >=(delta_5/(25/91))13^a
 =(12493/24750)13^a.                                 (23)
```

For example, common original deep-tooth valuations at least two, three, and
four force the integer boundary budgets

```text
a>=2: B>=86,             a>=3: B>=1109,
a>=4: B>=14417.                                      (24)
```

The estimate is uniform in the unit parts of the two deep coefficients and
uses no finite residue scan.

## 4. Application to the six-plus-one tail

In case (I), put

```text
a=nu_13(13v),                 M=13^a,       13v=Mw.   (25)
```

There is one divided target, so `sigma=measure(D_w)=1/7`.  Equations (5) and
(20) give

```text
H+sum_(i=1)^6 q_i
 >=7 delta_6 13^a
 =(191/990)13^a.                                     (26)
```

THM-2135 already proves `a>=5`.  Consequently every surviving
six-plus-one scalar row satisfies the explicit all-depth height floor

```text
H+sum_(i=1)^6 q_i>=ceil((191/990)13^5)=71634.         (27)
```

This does not conflict with THM-2135's local width invoice `v<=q_(4)`; that
invoice can itself force several large unit coefficients.  Equation (26) is
instead a quotient-support theorem: the common deep period cannot outrun the
total number of source boundary components.

## 5. Scope and the faithful quotient

The theorem does not bound the unit labels from above, so (23) and (26) do
not alone eliminate either scalar tail.  They do replace every arbitrarily
deep survivor by a proportionally high-boundary row and give a clean split
for subsequent finite-state, lattice, or component arguments.

The challenged assumption is that the quotient vertices should be the raw
roots of one thirteen-needle.  Here the vertices are **boundary events**, and
the quotient is the transfer operator (10).  Its preserved predicate is the
deep-cover support inclusion (12); a raw quotient destroys within-fibre
location, while the multiplicity and total-variation sidecars recover the
amount of residual mass and the cost of concentrating it.  Pairwise boundary
order can be oriented to form a tournament, but its scores, cycles, and
Hamiltonian paths do not determine variation after collisions.  The faithful
carrier is the signed boundary measure pushed through `[M]`, together with
the nonnegative root-count function `F_M`.  QED.
