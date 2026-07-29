# Ordered continuation is a poset, not a tournament

**Status:** theorem-grade abstract lemma and exact hostile controls; no theorem
ID is claimed.  The knot specialization uses only the proved/cited data in
THM-2176.  This note is intended for audit and possible later routing.

## 1. Universal ordered continuation

Let `(M,+,0)` be a commutative monoid, let `(Y,<=)` be a poset, and let
`f:M->Y`.  Define

```text
x <=_f^ctx y
    iff f(x+z)<=f(y+z) for every z in M.                 (1)
```

Then `<=_f^ctx` is the largest translation-compatible preorder contained in
the pullback of `<=` along `f`.

Indeed, reflexivity and transitivity hold pointwise.  If (1) holds, then for
every `a,z`,

```text
f((x+a)+z)=f(x+(a+z))<=f(y+(a+z))=f((y+a)+z),           (2)
```

so translation preserves the preorder.  Taking `z=0` shows that (1) implies
`f(x)<=f(y)`.  Conversely, if `R` is any translation-compatible preorder
such that

```text
x R y  =>  f(x)<=f(y),
```

then `x R y` implies `x+z R y+z` for every `z`, and hence (1).  This proves
maximality.

Its symmetric kernel is exactly THM-2176's continuation congruence:

```text
x <=_f^ctx y and y <=_f^ctx x
  iff f(x+z)=f(y+z) for every z
  iff Pi_f(x)=Pi_f(y).                                  (3)
```

Consequently the continuation quotient

```text
P_f=M/equiv_f
```

is a poset, and

```text
[x] -> Pi_f(x) in Y^M                                  (4)
```

is an order embedding for the pointwise product order.  This is the ordered
Myhill--Nerode refinement of THM-2176's equality theorem.

## 2. The intrinsic crossing relation

Now let `Y` be totally ordered.  For a pair define the complete sign set

```text
S_f(x,y)
 ={sgn(f(x+z)-f(y+z)):z in M}.                         (5)
```

After quotienting continuation ties, there are three intrinsic pair types:

```text
[x]<_f^ctx[y],     [y]<_f^ctx[x],     or
[x] cross_f [y],                                         (6)
```

where crossing means that both signs `+` and `-` occur in (5).  Equivalently,
`cross_f` is incomparability in `P_f`.  It is symmetric, not an orientation.

Every finite induced crossing graph is therefore an incomparability graph of
a finite poset.  Its cliques are contextual antichains and its proper color
classes are chains.  Dilworth's theorem gives

```text
chi(cross_f|A)=omega(cross_f|A)                         (7)
```

for every finite `A subset P_f`; in particular these graphs are perfect.
This is a genuine graph-class restriction, unlike a cosmetic tournament
orientation.

Translation retains the preorder but need not retain crossing.  More
precisely,

```text
S_f(x+a,y+a)
 ={sgn(f(x+t)-f(y+t)):t in a+M}
 subset S_f(x,y).                                      (8)
```

Thus a future cone can resolve an old crossing into dominance or a tie.  If
`M` is a group, `a+M=M` and equality holds in (8); group translations do
preserve every pair type.

## 3. Why the root projection is load-bearing

The same pointwise-order construction is trivial on THM-2242's full
two-ended metric kernels

```text
P_x(a,b)=d(x+a,b).
```

Indeed,

```text
P_x<=P_y pointwise
  => d(x,y)=P_x(0,y)<=P_y(0,y)=0
  => x=y.                                                (9)
```

Thus the two-ended kernel is faithful but has only equality as pointwise
dominance.  The nontrivial contextual poset appears only after projecting to
the root column `f(x)=d(x,0)` and then retaining all its translates.  This
projection is controlled forgetting: it loses enough endpoint information
to create dominance and incomparability, while the full continuation family
repairs scalar ties.  It is not licensed by the bare kernel embedding alone.

## 4. Knot specialization and the first exact crossing

For oriented knots, take `M` under connected sum and `f=u`.  Put

```text
K=T(2,7),              Kbar=mirror(K).
```

THM-2176 records

```text
u(K)=u(Kbar)=3,
u(K#K)=u(Kbar#Kbar)=6,
u(K#Kbar)<=5.                                          (10)
```

The two diagonal equalities use the signature lower bound and separate
unknotting; the off-diagonal inequality is the
Brittenham--Hermiller certificate.  Hence the two contexts `K` and `Kbar`
give opposite strict comparisons:

```text
u(K#K)>u(Kbar#K),
u(K#Kbar)<u(Kbar#Kbar).                                (11)
```

Therefore

```text
[K] cross_u [Kbar].                                    (12)
```

This sharpens “the profiles differ”: the scalar tie splits into an intrinsic
contextual antichain.  Mirroring exchanges the two vertices and the two
witness contexts, so the symmetric crossing edge is mirror-gauge invariant.
No orientation of this pair is mirror invariant.  A tournament shadow would
either break the gauge or identify two continuation-distinct knots.

The statement is orthogonal to THM-2281.  That theorem puts any *finite
metric packet* into one common optimal translated slice.  Ordered
continuation quantifies over the entire future cone.  Its failure is
finitely witnessed (one context for failed dominance, two for crossing),
but dominance itself is not supplied by a finite common-context theorem.
Likewise, in the language of THM-2348, one favorable context is analogous to
one Cartesian optimum: it does not imply the robust all-context condition.

## 5. Exact `V4` control, four-context tax, and the `S4/V4` quotient

Take

```text
M=F_2^2=V4,
f(0)=0,             f(x)=1 for x!=0.
```

This is the root length of the discrete translation-invariant metric.
For distinct `x,y`, context `z=x` makes

```text
f(x+z)=0<f(y+z),
```

while `z=y` reverses the inequality.  Thus every two distinct continuation
classes cross and the crossing graph is

```text
K4.                                                     (13)
```

This gives three exact boundaries at once.

1. `K4` contains triangles, so ordered continuation need not be a partial
   cube.
2. The unweighted crossing graph has automorphism group `S4`.  The affine
   translation subgroup is a normal regular `V4`, and its quotient is
   `S4/V4 isomorphic to S3`.
3. Since every transposition reverses one edge, no tournament orientation of
   `K4` is invariant under the full affine gauge.

There is also an exact finite-query tax.  For a finite packet `A`, call a
context dictionary `D` **order-complete** when

```text
x<=_f^ctx y
 iff f(x+z)<=f(y+z) for every z in D,               x,y in A.  (14)
```

In this `V4` model, the minimum order-complete dictionary for `A=V4` has
size exactly four.  To refute the ordered comparison `x<=y`, the only
possible context is `z=y`: there the right response is zero and the left is
one, whereas at every `z!=y` the right response is one.  Hence every
`y in V4` must belong to `D`, and `D=V4` works.

The abstract order dimension of a four-element antichain is only two.
Thus abstract linear extensions or the three `S3` direction channels do not
realize the physical context coordinates.  The missing fourth datum is the
affine origin.  This is an exact form of the user's “forbidden degree
`|V4|=4`” observation.

This is the rigorous place where the user's `V4`/`S4`/`S3` frame occurs.
Composing `PSL_2(Z)=C2*C3` with its finite `S3` quotient only acts on the
three nonzero affine directions.  It does not recover free-word order, braid
height, or any knot continuation value.

## 6. Sharp non-group hostile

Crossing is not itself an operation congruence.  Let `M=N` and use the
translation-invariant integer metric

```text
d(m,n)=0                 if m=n,
       =1                 if |m-n|=2,
       =2                 otherwise.                (15)
```

Every nontrivial distance is `1` or `2`, so the triangle inequality is
immediate.  Its root length is

```text
f(0)=0, f(2)=1, and f(n)=2 for n>0, n!=2.           (16)
```

The classes `1,2` cross:

```text
f(1)>f(2),               f(2)<f(3).                 (17)
```

After translating both by `1`, however,

```text
2 <_f^ctx 3,                                            (18)
```

because `f(2)<f(3)` at the zero context and
`f(2+z)=f(3+z)=2` for every `z>=1`.  This realizes strict inclusion in
(8), even for a jointly nonexpansive integer metric monoid.

## 7. Information ledger

```text
source:
  the continuation profile Pi_f;

target:
  its pointwise order, continuation quotient poset, and symmetric
  incomparability graph;

preserved:
  every future scalar comparison and every continuation tie;

lost by the crossing graph:
  magnitudes, witness contexts, dominance direction, and monoid operation;

gauge:
  continuation equality; in the mirror example, simultaneous mirror of
  object and context;

ties:
  literal equality of full profiles, not equality of f at the empty context;

cheapest decisive tests:
  one context refutes dominance, two opposite contexts prove crossing;

stopping boundary:
  crossing can disappear under a noninvertible translation, and neither a
  finite common-context slice nor an S3 tournament shadow recovers the full
  continuation order.                                          (19)
```
