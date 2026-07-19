---
id: THM-1137
title: The exact interval transfer closes two multiplicatively spread r=6 cones
status: PROVED — an arbitrary one-period window contains a safe subinterval of length at least 3/(7k), not 6/(7k), and the exact transfer is Phi(x)=min(6/7,(x-1/7)/2).  It closes every r=6 row with k1>=168 and adjacent ratios at least 7/3; once k1>=312, the sharper ratio 19/10 suffices.  No Covering assumption or computation is used.
source: codex-2026-07-18-S67 (audit and repair of the incoming four-comb gap-recursion scout)
depends_on:
  - THM-1132  # final open-danger-tooth component absorption; equality correction recorded there
related: [THM-1128, THM-1129, THM-1135, THM-1136, HYP-7607, MISTAKE-169]
---

# THM-1137 — the exact interval-transfer recursion

Let `P` be any subset of `{1,...,12}` and let

```text
k1<k2<k3<k4<k5<k6
```

be positive integers.  If

```text
k1 >= 168,
3*k_(i+1) >= 7*k_i       for i=1,2,3,4,              (1)
```

then `P union {k1,...,k6}` has a time `t` satisfying

```text
||v*t|| >= 1/14
```

for every speed `v` in the family.  Thus the multiplicatively spread cone of
the unbounded `r=6` branch is closed uniformly and without a covering
hypothesis.

There is a second, larger cone.  The same conclusion holds under

```text
k1 >= 312,
10*k_(i+1) >= 19*k_i      for i=1,2,3,4.              (1')
```

## The exact one-period lemma

For a positive integer `k`, put

```text
D_k={t: ||k*t||<1/14}.
```

This is an open periodic comb.  Its safe part in one period has length
`6/(7k)`.

> **Lemma.** If a closed real interval `J` has length `L>=1/k`, then it
> contains a closed subinterval disjoint from `D_k` of length at least
>
> \[
> \frac1k\Phi(kL),\qquad
> \Phi(x)=\min\left(\frac67,\frac{x-1/7}{2}\right).    \tag{2}
> \]
>
> In particular the length is at least `3/(7k)`.

Scale by `k`, so the window has length `x=kL>=1`.  Modulo one, the safe set is
the single closed circular arc

```text
[1/14,13/14]
```

of length `6/7`.  If `x>=13/7`, every window contains a complete translate of
that arc, giving `6/7`.  If `1<=x<=13/7` and the window contains no complete
safe arc, it meets at most two adjacent safe arcs, separated by one complete
danger gap of length `1/7`.  Their total length is at least `x-1/7`, so the
larger has length at least `(x-1/7)/2`.  This proves (2); at `x=1` it gives
`3/7`.

The factor `3/7` is sharp for an arbitrary window: at scale `k=1`,

```text
J=[1/2,3/2]
```

cuts the safe arc into exactly `[1/2,13/14]` and `[15/14,3/2]`, both of
length `3/7`.  This is the endpoint configuration missed by the earlier
`6/(7k)` recursion (MISTAKE-169).

## A universal core interval

The closed interval

```text
I0=[1/14,13/168]                                      (3)
```

has length `1/168` and is safe for every speed `p` in `{1,...,12}`.  Indeed,
for `t in I0`,

```text
1/14 <= p*t <= 13/14,
```

because the lower extreme occurs at `p=1,t=1/14` and the upper extreme at
`p=12,t=13/168`.  Hence (3) is safe for every possible core `P`.

Condition `k1>=168` gives `|I0|>=1/k1`.  Apply the one-period lemma to obtain
a subinterval `I1` safe for `P,k1` with

```text
|I1| >= 3/(7k1).
```

Suppose recursively that `Ii` is safe for `P,k1,...,ki` and has length at
least `3/(7ki)`.  From (1),

```text
3/(7ki) >= 1/k_(i+1).
```

The lemma therefore supplies `I_(i+1)` safe for the next killer and of length
at least `3/(7k_(i+1))`.  Iterating through `k5` gives a closed residual
interval

```text
|I5| >= 3/(7k5) > 1/(7k6).                            (4)
```

Every connected component of the open danger comb `D_(k6)` has length exactly
`1/(7k6)`.  Equation (4) prevents `I5` from lying inside one such component,
so `I5` contains a point safe for `k6`.  That point is safe for all thirteen
speeds, proving the theorem.

## The sharper `19/10` cone

Formula (2) retains the normalized interval width instead of discarding it at
`3/7`.  Under (1'), the initial normalized length is

```text
k1*|I0| >= 312/168 = 13/7,
```

so after removing `k1` we may start with `c1=6/7`, where `|I1|>=c1/k1`.
Since `Phi` is monotone, four successive ratios at least `r=19/10` give the
exact lower ledger

```text
c1 =       6/7,
c2 =      26/35,
c3 =     111/175,
c4 =    1859/3500,
c5 =  212247/490000 > 1/7.                            (5)
```

Here `c_(i+1)=Phi(r*c_i)`.  The last application is legal because

```text
r*c4 = 35321/35000 > 1,
```

and all previous normalized inputs are larger.  Equation (5) gives

```text
|I5| >= c5/k5 > 1/(7k5) > 1/(7k6),
```

so the same final-tooth argument closes the family.  This proves (1').  The
numbers in (4) are exact rational consequences of one one-dimensional transfer
function, not fitted constants.

## Frontier left by the theorem

In the standard seven-core/six-killer stratum, every remaining tuple obeys

```text
k1 <= 167; or
168 <= k1 <= 311 and 3*k_(i+1) < 7*k_i for some i; or
k1 >= 312 and 10*k_(i+1) < 19*k_i for some i.          (6)
```

Thus failure of this proof is itself structural: the first five killers must
start in a short bottom strip or contain a multiplicatively clustered edge.
THM-1136 cuts (6) in a transverse coordinate, by the position of the first
modulo-13 carrier and the number of distinct nonzero residue classes.
THM-1135 cuts it again by harmonic load and much larger proper-prefix ratios.
The honest residual is their intersection, not the original unstructured
six-dimensional tail.

## Tournament and assumption challenge

Ordering the killers gives the usual transitive tournament, with score
histogram `0,1,2,3,4,5`, no directed cycle, singleton SCCs, and one Hamiltonian
path.  That tournament forgets both the absolute base gate `k1>=168` and the
metric edge labels `3k_(i+1)-7k_i`, so it does not preserve (1).

The proof-bearing vertices are instead the successive **interval obligations**
`I0 -> I1 -> ... -> I5`.  An edge records a period-fitting inequality and the
amount of interval width that survives.  This challenges the assumption that
runners, danger arcs, or residues must be the vertices: here the faithful
object is a labelled proof-state path.  Its naked orientation is transitive
telemetry; its metric labels prove the LRC predicate.

∎
