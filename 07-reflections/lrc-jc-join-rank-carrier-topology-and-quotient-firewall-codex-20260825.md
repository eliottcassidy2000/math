# LRC--JC join rank, carrier topology, and the quotient firewall

**Status.** The join-rank, quotient-descent, and signed-cover statements below
are **PROVED ELEMENTARY**. Their LRC use is present in THM-4151, with
THM-4148/4136 as the prior carrier mechanism. Their planar-Jacobian inputs are
the proved THM-4138 and THM-4146 formulas. THM-4147 now applies the merger-rank
mechanism to exclude the dense generic eta-only exact-`M=9` chamber. The
zeta-only and mixed chambers, their cubic carrier, and specialization walls
remain **OPEN**.

## 1. One capacity lemma, two different lattices

Let `E` be finite.  In either the Boolean lattice of subsets of `E`, with
rank `r(A)=|A|`, or the partition lattice of `E`, with

```text
r(P)=|E|-number_of_blocks(P),
```

one has

```text
r(join_i x_i) <= sum_i r(x_i).                         (1)
```

For subsets this is the union bound.  For partitions, choose a spanning tree
inside every block of every `x_i`.  The components of the union graph are the
blocks of `join_i x_i`; a graph with that component partition needs at least
`r(join_i x_i)` forest edges.  This proves `(1)`.

For the LRC two-sheet carrier, put `E={0,1}` and, for odd `r`,

```text
K_r(z)={i in E: ||r(z+i/2)||<1/14}.                    (2)
```

The two arguments differ by `1/2`, so `|K_r(z)|<=1`.  Total failure for tails
`p,q` says `K_p(z) union K_q(z)=E`.  Equality is therefore forced throughout

```text
2=|E|<=|K_p(z)|+|K_q(z)|<=2:                           (3)
```

each tail kills one sheet and the killers are complementary.  These are
exactly the two cross-comb orientations.  THM-4151 uses this Boolean equality
before exploiting the stronger anchored odd-wall slack; THM-4148 instead uses
the component-width consequence after quotienting.

For permutations `g_j in S_n`, let `P(g_j)` be the partition into cycles and
let `c(g_j)` count all cycles, including fixed points.  The join of these
partitions is the orbit partition of the generated group.  Hence transitivity
forces the exact merger-capacity inequality

```text
sum_j (n-c(g_j)) >= n-1.                               (4)
```

If `g` has support `d` split into `k` nontrivial cycles, then
`n-c(g)=d-k<=d-1`.  Thus `(4)` recovers THM-4138's support-minus-one lemma and
shows what its coarse version forgets: a 4-cycle has rank `3`, while
`(12)(34)` has the same support but rank `2`.

## 2. Quotient descent is preservation of the fibre partition

For a finite map `f:E->Y`, a permutation `sigma` descends through `f` only if

```text
f(x)=f(y)  implies  f(sigma x)=f(sigma y);             (5)
```

with the corresponding condition for `sigma^{-1}`.  Conversely `(5)` makes
the induced map on fibre classes well-defined.  Sheet count alone is not a
descent certificate.

The LRC total-failure predicate passes this test: the deck involution
`z->z+1/2` swaps the two physical sheets, while “both sheets are bad” is
invariant.  The oriented statement “tail `p` kills sheet zero” is
anti-invariant and lives only upstairs.

THM-4146 supplies the exact hostile.  At `a=-48`, let

```text
R={-6,6,0},                    sigma=(-6 6 0).          (6)
```

The horizontal response `H` has fibre partition

```text
{{-6,6},{0}},
H(-6)=H(6)=(12,0,-8640),
H(0)=(-24,0,-55296).                                  (7)
```

But `H(sigma(-6))=H(6)` and `H(sigma(6))=H(0)` are unequal.  Thus `sigma`
does not preserve `(7)` and cannot descend to the horizontal image.  The
source root divisor and its order-three permutation are real; a horizontal
surface action or BC-cover monodromy is not.

## 3. Why compact-to-open escape does not cross to the JC carrier

On every THM-4148 danger component, the doubling cover restricts over an open
interval.  Its sheet local system is trivial, and a compact image trapped in
that open component must have strictly smaller length than the component.
This is the load-bearing ordered, bounded-component topology behind the
compact-to-open argument.

THM-4138's horizontal normalization instead has

```text
q=q_0+kappa z^2,
q_0=a^3/2,                       kappa=9a^2/16.         (8)
```

On the punctured complex `q`-plane, the compact loop

```text
q(theta)=q_0+kappa R^2 exp(i theta),    0<=theta<=2pi  (9)
```

lifts from `z=R` to `z=-R`.  It remains inside one open connected set and can
have arbitrarily large length, but it creates nontrivial sheet monodromy
rather than escaping.  At `a=1,q_*=1/4`, its two normalization sheets are
`z=+-2i/3`, giving exactly THM-4138's labelled points

```text
Q_+-=(1/18,+-11i/54).                                 (10)
```

Their separate meridians and transpositions are therefore load-bearing.
Degree two is the common shadow; interval order/component width versus the
nontrivial signed local system is the decisive non-bridge.

## 4. Connection contract

```text
source:       LRC physical-sheet kill sets and open carrier components
target:       JC cycle partitions, normalization sheets, labelled meridians
map:          local obstruction -> Boolean/partition-lattice element;
              double cover -> its Z/2 sheet local system
preserved:    top-join obligation, rank capacity, invariant predicates
destroyed:    interval metric, group word, branch divisor, oriented labels,
              q-fibre partition when only the source divisor is retained
sidecar:      LRC tooth endpoints and actual larger tail;
              JC full cycle types, branch set, fibre partition, meridians
decisive test: equations (3)--(5), then sheet monodromy on a base generator.
```

## 5. Exact `M=9` consumer and next parity test

THM-4147's eta-only chamber has two labelled transposition meridians and a
strict merger deficit, so the coarse support-minus-one consequence of `(4)`
already excludes its finite-carrier degree. For the remaining zeta-only and
mixed chambers, the next lawful audit is not another support census. Compute
every labelled generator's complete cycle partition and test `(4)`. Also use

```text
sign(g)=(-1)^(n-c(g)).                                  (11)
```

For a punctured-surface relation

```text
product_i [A_i,B_i] product_l M_l=1,
```

commutators are even, so

```text
sum_l (n-c(M_l)) = 0 mod 2.                             (12)
```

In the THM-4138 model, two simple transposition meridians force the remaining
puncture monodromy to have even cycle rank. In the open cubic-carrier chambers,
`(12)` can shave one merger unit when a support upper bound has the wrong
parity. Any further promotion still requires cubic polynomial-section and
loop-avoidance theorems, labelled-generator completeness, all specialization
walls, frozen outputs, and an independent audit. Neither `(4)` nor `(12)`
supplies entry into the reduced seam or proves `JC(2)`.
