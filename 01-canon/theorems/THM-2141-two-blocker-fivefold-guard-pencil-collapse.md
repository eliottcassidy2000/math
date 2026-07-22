---
id: THM-2141
title: "Two-target character rigidity empties the two-blocker fivefold pencil"
status: >
  PROVED. The rank-eight guard-pencil profile with two true terminal
  13-blockers, five guard-parallel nonblockers, and one transverse terminal
  is impossible. A clean 13^2-root column count forces each asymmetric
  guard/aligned phase rectangle into the union of the two divided-blocker
  bands. A two-target character lemma proves that, for every aligned label,
  either one blocker equals that label up to sign or one blocker lies on the
  guard line. Five distinct positive labels therefore put one blocker and at
  least four aligned terminals on the guard line. If the second divided
  blocker is off that rational line, THM-2080 leaves at most three fibre
  bands; if it is collinear, THM-2138's all-depth 5+2 scalar closure applies.
  Both alternatives give a positive-measure escape. The profiles (2,6,0)
  and (3,5,0) remain open.
source: codex-2026-07-22-LRC-two-target-fivefold-pencil
depends_on:
  - THM-1166
  - THM-2080
  - THM-2125
  - THM-2138
  - THM-2139
related:
  - THM-2120
  - THM-2140
---

# THM-2141 -- the two-blocker fivefold guard pencil is empty

Let `Gamma` be a rank-two character lattice and

```text
K=Hom(Gamma,R/Z).                                     (1)
```

Take a nonzero guard, two true blockers, five aligned nonblockers, and one
transverse terminal

```text
g,       c_*1=13u_1, c_*2=13u_2,
c_1,...,c_5,                         d in Gamma.       (2)
```

Assume `g mod 13` is nonzero,

```text
c_i mod 13 in F_13^* g mod 13,          i=1,...,5,
d mod 13 notin F_13 g mod 13,                         (3)
```

and retain the LRC specialization: an integral cocharacter `h` makes `g.h`
positive odd and all eight terminal values positive and pairwise distinct.
Suppose, toward a contradiction, that outside a null set the eight danger
bands cover

```text
{X:||g.X||>1/7}.                                     (4)
```

This is the exact `(b,r,t)=(2,5,1)` profile left by THM-2125.

## 1. Every aligned phase rectangle is absorbed by two blockers

For each `i=1,...,5`, the root-column ledger gives

```text
{Y:||g.Y||<1/7, ||c_i.Y||<1/14}
 subset {Y:||u_1.Y||<=1/14} union {Y:||u_2.Y||<=1/14}. (5)
```

Indeed, suppose the strict reverse of (5) were nonempty.  On a small inverse
domain of `[13]:K->K`, remove the inverse images of the original null
exceptional set under all 169 root sheets and choose a clean phase `Y`.
Both blockers are constant and safe on its full root plane.  The guard
forbids exactly three outer columns, the selected aligned label is a
singleton, and the other four aligned labels forbid at most two columns
each.  Their total capacity is

```text
3+1+4*2=12<13.                                       (6)
```

A surviving column is guard-safe and safe for all five aligned labels.  The
transverse terminal hits at most two of its thirteen inner points, leaving a
clean uncovered root.  This proves (5), including the almost-everywhere
null-set sidecar.

## 2. A two-target asymmetric character lemma

We need the following strengthening of THM-2139's one-target lemma.

> **Lemma.** Let `a,b,w_1,w_2` be nonzero characters of a connected
> rank-two torus and suppose
>
> ```text
> {||a.X||<1/7, ||b.X||<1/14}
>  subset {||w_1.X||<=1/14} union {||w_2.X||<=1/14}. (7)
> ```
>
> Then:
>
> 1. if `a,b` are rationally independent, one of `w_1,w_2` equals `+b` or
>    `-b`;
> 2. if `a,b` are rationally dependent, one of `w_1,w_2` lies on their
>    common saturated rational line.

### The independent case

Let `L=ker(a) intersection ker(b)`.  It is finite.  The danger proportion of
a nontrivial character on a finite group is at most one half: if its image
has order `m>=2`, the closed radius-`1/14` arc contains

```text
2 floor(m/14)+1 <=m/2                               (8)
```

image points, with the small orders checked by the same formula.  Both
danger sets contain the identity.  If both `w_j|L` were nontrivial, their
union would therefore have relative size strictly below one, contradicting
(7) on `L`.  Thus at least one target kills `L`.

Suppose exactly one does, say `w_1`.  Exact annihilator duality gives

```text
w_1=ma+nb,                         m,n in Z.           (9)
```

On every fibre of `(a,b)`, the other restriction remains nontrivial and its
translated danger set cannot contain the whole finite fibre.  If `w_1` were
safe at any source phase, (7) would require precisely that impossible full
containment.  Hence the one band of `w_1` contains the entire source
rectangle.  THM-2139's asymmetric lemma then gives

```text
w_1=+b or w_1=-b.                                    (10)
```

It remains to treat the case in which both targets kill `L`.  Write

```text
w_j=m_j a+n_j b.                                     (11)
```

Use `(x,y)=(a.X,b.X)` as Haar coordinates on the source rectangle

```text
R={|x|<1/7, |y|<1/14},                measure(R)=2/49. (12)
```

If `m_j!=0`, then on every fixed `y`-slice the full `x`-circle danger set has
measure `1/7`; its intersection with `|x|<1/7` has length at most `1/7`.
After integrating over the `y`-interval,

```text
measure(R intersection D_(w_j))<=1/49.               (13)
```

If `m_j=0` and `w_j` is not `+/-b`, then `|n_j|>=2`.  THM-1166's exact pair
formula gives

```text
measure(D_1 intersection D_(n_j))<=1/14              (14)
```

on the `y`-circle.  For completeness, (14) is immediate for
`|n_j|=2,3,4`; for `|n_j|>=5`, the correction term in that formula has
absolute value at most `49/(196|n_j|)`, which gives the same bound.  Multiplying
by the `x`-length `2/7` again proves (13).

If neither target is `+/-b`, both obey (13).  Their intersections with `R`
overlap on a positive neighborhood of `(0,0)`, so their union has measure
strictly less than `2/49`.  It cannot cover `R`.  Therefore one target is
`+/-b`, proving the independent clause.

### The dependent case

Write `a,b` on the saturated line generated by a primitive `alpha`.  The
connected circle `ker(alpha)` lies in the source of (7).  If both target
restrictions were nontrivial, each danger band would have Haar measure `1/7`
on that circle, and their union would have measure at most `2/7<1`.
Therefore one target kills `ker(alpha)` and belongs to `Z alpha`.  This proves
the lemma. QED.

## 3. Five labels force a fourfold rational line

Apply the lemma to (5) with

```text
(a,b,w_1,w_2)=(g,c_i,u_1,u_2).                       (15)
```

If all five `c_i` were independent of `g`, each would equal one of
`+/-u_1,+/-u_2`.  Positivity under `h` excludes the minus signs, and
distinctness leaves at most two labels, a contradiction.  Hence some `c_i`
is proportional to `g`, and the dependent clause puts one blocker on that
line.  Relabel so

```text
u_1 in Qg.                                            (16)
```

For any `c_k` independent of `g`, the independent clause cannot identify it
with `+/-u_1`, so it must equal `+/-u_2`.  Positivity again fixes the plus
sign, and distinctness allows at most one such label.  Consequently

```text
g,u_1 and at least four of c_1,...,c_5
                         lie on one rational line.    (17)
```

## 4. Both positions of the second blocker give an escape

Let `alpha` generate the line in (17).

### The second blocker is off the line

If `u_2 notin Q alpha`, the line contains the guard, the blocker `13u_1`,
and four or five aligned terminal characters.  Their positive coefficients
are distinct and the guard coefficient is odd.  THM-2080 therefore supplies
a positive-measure scalar phase on which the guard and all five or six line
terminal bands are strictly safe.

On an `alpha`-fibre over such a phase, at most three off-line bands remain:

```text
13u_2,
possibly the unique independent c_k=u_2,
and d.                                                (18)
```

Each restricts nontrivially and has fibre measure `1/7`.  Their union has
measure at most `3/7<1`, so Fubini supplies a positive-measure uncovered set.

### The second blocker is on the line

If `u_2 in Q alpha`, no `c_i` can be independent of `g`: the lemma would
identify it with one of the two line characters.  Hence all five `c_i` lie on
the line.  The scalar line data are exactly a unit guard, five unit terminal
coefficients, and the two deeper coefficients belonging to `13u_1,13u_2`.
THM-2138 proves that this full `5+2` scalar containment is impossible at every
thirteen-adic depth.  Choose a positive-measure scalar escape and then avoid
the sole transverse band `D_d` on its circle fibre.  Again Fubini gives a
positive-measure contradiction to (4).

Both cases are impossible.  Therefore the `(2,5,1)` fivefold guard pencil is
empty.

## 5. Scope and Tournament Analysis

Together THM-2139 and this theorem remove the two exact minimal-fivefold
profiles

```text
(1,5,2),                    (2,5,1).                 (19)
```

The fivefold branch now consists of

```text
(1,6,1),(1,7,0),(2,6,0),(3,5,0).                    (20)
```

The challenged assumption is that two blocker bands merely weaken a
one-target cross lemma.  Their restrictions to the finite common kernel form
a two-colour cover obligation; overlap at the identity prevents two
nontrivial colours from covering that kernel.  Candidate tournament vertices
are aligned labels, the two blocker colours, common-kernel sheets, and source
rectangles.  Orienting a label toward one qualifying blocker gives a
bipartite choice tournament with a tie Hamiltonian path, but scores, cycles,
SCCs, edge flips, and path counts forget the continuous rectangle measure in
(13).  The faithful carrier is the two-coloured character-kernel relation
with its finite-fibre and Haar-area sidecars.  QED.
