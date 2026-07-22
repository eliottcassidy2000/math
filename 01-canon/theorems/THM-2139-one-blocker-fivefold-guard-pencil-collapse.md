---
id: THM-2139
title: "One blocker and five guard-parallel terminals collapse to one line"
status: >
  PROVED. In the rank-eight fivefold guard-pencil branch with exactly one
  terminal 13-blocker, exactly five nonblocker terminals on the guard's
  projective mod-13 line, and two transverse terminals, an almost-everywhere
  cover is impossible. On a clean 13^2-root fibre, one aligned singleton
  phase and the three guard-unsafe columns have total column capacity twelve.
  This forces E_g intersect D_c inside the divided blocker band for each of
  the five aligned characters. An integral character-kernel lemma then makes
  the guard, divided blocker, and all five aligned terminals rationally
  collinear. THM-2080 supplies a scalar escape from those six terminal combs,
  and the last two transverse bands cannot cover its circle fibre. This closes
  the minimal profile (blockers,aligned,transverse)=(1,5,2), not the other
  fivefold-pencil profiles or LRC(14).
source: codex-2026-07-22-LRC-fivefold-guard-pencil
depends_on:
  - THM-2080
  - THM-2125
related:
  - THM-2120
  - THM-2123
  - THM-2138
external: Settled Lonely Runner Conjecture for at most six terminal speeds, through THM-2080.
---

# THM-2139 -- the minimal fivefold guard pencil is empty

Let `Gamma` be a rank-two character lattice and

```text
K=Hom(Gamma,R/Z).                                     (1)
```

Take a nonzero guard, one true terminal blocker, five aligned nonblockers,
and two transverse terminals

```text
g,       c_*=13u,       c_1,...,c_5,       d_1,d_2 in Gamma. (2)
```

Assume that these characters span `Gamma tensor Q`, that `g mod 13` is
nonzero, and that

```text
c_i mod 13 in F_13^* g mod 13,          i=1,...,5,
d_j mod 13 notin F_13 g mod 13,         j=1,2.        (3)
```

Retain the LRC specialization: there is an integral cocharacter `h` for
which `g.h` is positive and odd and all eight terminal values in (2) are
positive and pairwise distinct.  Suppose, toward a contradiction, that
outside a null set the terminal danger bands cover the guard-safe region:

```text
{X:||g.X||>1/7}
 subset D_(c_*) union union_(i=1)^5 D_(c_i)
                    union D_(d_1) union D_(d_2),      (4)
D_a={X:||a.X||<1/14}.                                (5)
```

This is exactly the `(b,r,t)=(1,5,2)` boundary case left by THM-2125.

## 1. A singleton aligned phase leaves a whole column

For every aligned label `i`, we claim the phase containment

```text
{Y:||g.Y||<1/7, ||c_i.Y||<1/14}
                  subset {Y:||u.Y||<=1/14}.           (6)
```

Suppose instead that the open set

```text
O_i={Y:||g.Y||<1/7, ||c_i.Y||<1/14,
                         ||u.Y||>1/14}                (7)
```

is nonempty.  Use the finite covering `[13]:K->K` and choose a local inverse
sheet over a small open subset of `O_i`.  Choose a basis `z,v` of `K[13]`
such that

```text
g.z=1/13,                         g.v=0.              (8)
```

The 169 local roots of `Y` are indexed by

```text
X_jk(Y)=sigma(Y)+jz+kv,              (j,k) in F_13^2. (9)
```

Let `E` be the exceptional null set in (4).  Each local sheet in (9) is a
local diffeomorphism, so the set of `Y` for which any one of the 169 roots
lies in `E` is null.  Choose `Y in O_i` outside their finite union.  Every
root in (9) is then clean.

The blocker is constant and safe on all roots:

```text
c_*.X_jk(Y)=u.Y.                                      (10)
```

Every aligned `c_l` kills `v` modulo thirteen, hence is constant on the
inner `k`-column.  Across the thirteen columns, its values form a translated
thirteen-grid.  Such a terminal band occupies at most two columns.  For the
chosen label `i`, the strict phase in (7) makes it occupy exactly one column:
the thirteen values are the roots of `c_i.Y`, and a direct endpoint count
gives one root when `||c_i.Y||<1/14`.

The guard is also constant on each column.  Because `||g.Y||<1/7`, exactly
three of its thirteen root values have norm at most `1/7`.  Thus the columns
forbidden by the guard and the five aligned terminals have total capacity at
most

```text
3+1+4*2=12<13.                                       (11)
```

Choose a column on which the guard and all five aligned terminals are safe.
On that column, transversality in (3) makes each `d_j` a complete translated
thirteen-grid as `k` varies.  Each of their danger bands hits at most two
points, so together they hit at most four.  Some root is safe for both.  It
is also blocker-safe by (10), contradicting the clean cover (4).  Hence (6)
holds for every `i`.

The use of an open counterexample set in (7), followed by all 169 local
sheets, is the null-set sidecar: a single arbitrarily chosen root need not be
clean, but a nonempty phase violation always has a clean full root fibre.

## 2. The asymmetric character-kernel lemma

We use the following rigidity statement.

> **Lemma.** Let `a,b,w` be nonzero characters of a connected rank-two
> torus.  If
>
> ```text
> {||a.X||<1/7, ||b.X||<1/14}
>                       subset {||w.X||<=1/14},       (12)
> ```
>
> then:
>
> 1. if `a,b` are rationally independent, `w=+b` or `w=-b`;
> 2. if `a,b` are rationally dependent, all three characters lie on their
>    common saturated rational line.

### Proof

In the independent case the torus map `(a,b):K->T^2` is surjective with a
finite kernel `L`.  On `L`, the left side of (12) holds with both phases
zero.  The compact subgroup `w(L)` is contained in the closed radius-`1/14`
arc, so it is trivial.  Exact character duality therefore gives

```text
w=ma+nb,                         m,n in Z.            (13)
```

The image of a small real rectangle under the lifted linear form in (13)
must remain inside the target arc.  That lifted image interval is connected
and contains zero, so it cannot wrap to a different lift of the target arc.
Letting its two coordinates approach the source radii with signs matching
`m,n` gives

```text
|m|/7+|n|/14<=1/14,
2|m|+|n|<=1.                                          (14)
```

Since `w` is nonzero, (14) says `m=0,n=+/-1`.

In the dependent case, write `a,b` on the saturated line generated by a
primitive character `alpha`.  The connected circle `ker(alpha)` lies in the
left side of (12).  The image of its restriction under `w` is again a compact
subgroup contained in a radius-`1/14` arc, hence is trivial.  Therefore `w`
also lies in `Z alpha`.  This proves the lemma. QED.

## 3. Five containments force one rational line

Apply the lemma to (6) with

```text
(a,b,w)=(g,c_i,u).                                    (15)
```

For each aligned terminal, either

```text
c_i is rationally proportional to g,
or
u=+c_i or u=-c_i.                                    (16)
```

The five `c_i` are distinct.  Only two characters can equal `+u` or `-u`,
so at least one `c_i` is proportional to `g`.  The dependent clause of the
lemma then puts `u` on the same rational line as `g`.  If some other `c_j`
were independent of `g`, the second alternative in (16) would put it on the
`u`-line and hence on the `g`-line, a contradiction.  Consequently

```text
g,u,c_1,...,c_5 all lie on one saturated rational line. (17)
```

This is stronger than their initial projective agreement modulo thirteen;
one singleton column at a time has lifted the entire fivefold pencil to
characteristic zero.

## 4. The line has a six-comb escape

Let `alpha` be the primitive generator of the line in (17), oriented using
the cocharacter `h`.  Write

```text
g=H alpha,        u=w alpha,        c_i=q_i alpha.    (18)
```

Then `H,w,q_i` are positive, `H` is odd, and the six terminal coefficients

```text
13w,q_1,...,q_5                                      (19)
```

are distinct.  THM-2080 supplies a positive-measure set of scalar phases
`t` on which

```text
||Ht||>1/7,
||13wt||>1/14,             ||q_i t||>1/14 for all i. (20)
```

The two transverse characters `d_1,d_2` restrict nontrivially to every
circle fibre `alpha.X=t`; otherwise their reductions would lie on the guard
line in (3).  On each such fibre either danger band has Haar measure `1/7`,
so their union has measure at most `2/7`.  Fubini applied over the
positive-measure phase set in (20) gives a positive-measure family of points
that avoid both transverse bands while retaining every inequality in (20).
This family lies in the guard-safe set and outside every terminal danger
band, contradicting (4).

Therefore the `(1,5,2)` fivefold-pencil profile is empty.

## 5. Scope and Tournament Analysis

After this theorem, the profiles still allowed by the THM-2125 count are

```text
(b,r,t)=(1,6,1),(1,7,0),(2,5,1),(2,6,0),(3,5,0).    (21)
```

They require higher singleton multiplicity or several divided-blocker target
bands; the one-target containment (6) cannot simply be reused.

The challenged assumption is that the 169 roots should be flattened into
one incidence graph.  The faithful vertices are the thirteen guard columns,
each with an inner thirteen-point kernel fibre.  The pair observable is
whether an aligned terminal spends one or two whole columns; the switch is
its quotient phase `c_i.Y`, and cyclic column order supplies a tie
Hamiltonian path.  Scores, cycles, SCCs, edge flips, and path counts forget
that one singleton changes total outer capacity from thirteen to twelve.
The preserved carrier is the two-level root incidence together with the
integral character-kernel containment (12) and the clean-sheet null-set
sidecar.  QED.
