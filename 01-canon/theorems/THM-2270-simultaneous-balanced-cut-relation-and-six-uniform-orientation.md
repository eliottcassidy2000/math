---
id: THM-2270
title: "Simultaneous balanced-cut relation and six-uniform orientation"
status: >
  PROVED + INDEPENDENTLY AUDITED from THM-2145. Every zero-safe row of
  thirteen distinct positive speeds has one integer relation of coefficient
  height at most 504576 whose weighted partial sum is nonzero on every one of
  the 1716 labelled 6+7 cuts simultaneously. In particular, the relation has
  support at least eight. Its cut signs form a nonconstant threshold
  orientation of the labelled six-subsets; after sorting the thirteen
  weighted coordinates, the positive sets form a Gale up-set. The chamber
  graph of these orientations is a partial cube. This is a simultaneous
  high-support certificate, not a tournament, a height improvement, a new
  rank increment, or a proof of LRC(14).
source: codex-2026-07-25-simultaneous-balanced-cut
depends_on:
  - THM-2145-two-block-spectral-crossing-and-6-plus-7-carry
related:
  - THM-2164-relative-packet-rank-harvesting
  - THM-2166-hybrid-core-smoothing-low-carry-crossing
  - THM-2266-depth-one-deep-pair-centered-signed-dual-and-relation-atlas
  - THM-2267-static-owner-coverage-is-flag-and-transition-holonomy-is-a-cut-kernel
---

# THM-2270 -- simultaneous balanced-cut relation

**PROVED + INDEPENDENTLY AUDITED.**

THM-2145 gives a bounded relation separately for every labelled `6+7`
cut. A separate choice for every cut forgets whether the relations can be
made compatible. The coefficient-grid argument below retains all cuts at
once:

```text
1716 cutwise height-584 relations
  -> one height-504576 relation crossing every cut
  -> support at least eight
  -> one global threshold sign chamber.                         (1)
```

The larger height is the price of simultaneity. THM-2164 remains much
stronger for low-height rank two.

## 1. Statement

Let

```text
v=(v_1,...,v_13)
```

be pairwise distinct positive integers, and put

```text
G={x in R/Z: ||x||>=1/14},
S(v)={t: v_i t in G for every i},

Lambda(v)={m in Z^13: sum_i m_i v_i=0}.               (2)
```

For a labelled six-subset `A subset [13]`, define the cut functional

```text
L_A(m)=sum_(i in A)m_i v_i.                            (3)
```

There are

```text
q=binomial(13,6)=1716                                  (4)
```

such cuts. If

```text
mu(S(v))=0,                                            (5)
```

then there is an `m in Lambda(v)` such that

```text
||m||_infinity<=504576,

L_A(m)!=0          for every A in binomial([13],6).    (6)
```

Consequently,

```text
|supp(m)|>=8.                                          (7)
```

The conclusion is simultaneous: the same `m` works in all `1716`
inequalities in (6).

## 2. The cutwise input

Fix a six-subset `A`, and use `A` and its seven-element complement as the
labelled blocks in THM-2145. Equations (5) and THM-2145 give

```text
m_A in Lambda(v),
||m_A||_infinity<=584,
L_A(m_A)!=0.                                           (8)
```

The last inequality is precisely the genuine-crossing clause: since
`m_A.v=0`, the partial sum on the complementary block is `-L_A(m_A)` and
is nonzero as well.

Let

```text
W=span_Q{m_A: |A|=6}.                                  (9)
```

Choose a basis

```text
b_1,...,b_r
```

from the family in (8). Basis extraction therefore preserves

```text
||b_j||_infinity<=584.                                (10)
```

Moreover,

```text
r<=dim_Q(v^perp)=12.                                  (11)
```

For every `A`, the restriction `L_A|W` is nonzero because it is already
nonzero on `m_A in W`.

## 3. Centered coefficient-grid lemma

We use the following elementary finite nonvanishing fact.

> **Coefficient-grid lemma.** Let `W` have rational basis
> `b_1,...,b_r`, and let `ell_1,...,ell_q` be nonzero rational linear
> functionals on `W`. If `q` is even, there are integers `t_1,...,t_r`
> such that
>
> ```text
> ell_a(sum_j t_j b_j)!=0             for every a,
>
> sum_j |t_j|<=q/2+floor(r/2).                         (12)
> ```

### Proof

Form the polynomial

```text
P(T_1,...,T_r)
 =product_(a=1)^q ell_a(sum_j T_j b_j).                (13)
```

Each factor is a nonzero linear polynomial. Since a polynomial ring over
`Q` is an integral domain, `P` is nonzero and homogeneous of total degree
`q`.

Choose a monomial

```text
T_1^(d_1)...T_r^(d_r),        sum_j d_j=q,             (14)
```

with nonzero coefficient in `P`. For each `j`, use the centered consecutive
integer set

```text
S_j={-floor(d_j/2),...,+ceil(d_j/2)}.                 (15)
```

It has `d_j+1` elements. The standard coefficient-grid identity, obtained
by applying one-variable Lagrange interpolation successively in the
variables, says that the coefficient in (14) is a weighted sum of the
values of `P` on

```text
S_1 x ... x S_r.                                      (16)
```

If `P` vanished on the entire grid, that coefficient would be zero.
Therefore there is a grid point `t=(t_1,...,t_r)` with `P(t)!=0`.

Let `o` be the number of odd exponents among the `d_j`. Since their sum
`q` is even, `o` is even. Hence

```text
sum_j |t_j|
 <=sum_j ceil(d_j/2)
 =(q+o)/2
 <=q/2+floor(r/2).                                    (17)
```

Every factor in (13) is nonzero at this point, proving (12). QED.

This is the coefficient form of the Combinatorial Nullstellensatz, but the
interpolation paragraph is a complete proof of the instance used here.

## 4. Simultaneous extraction and the exact height

Apply the lemma to the space `W`, the basis from (10), and the `1716`
restrictions `L_A|W`. Equations (4), (11), and (12) give

```text
sum_j |t_j|<=1716/2+12/2=864.                         (18)
```

Set

```text
m=sum_(j=1)^r t_j b_j.                                (19)
```

Because each basis vector is an integral relation, `m in Lambda(v)`.
Every factor of (13) is nonzero at `t`, so (6) holds. Finally,

```text
||m||_infinity
 <=sum_j |t_j| ||b_j||_infinity
 <=864*584
 =504576.                                             (20)
```

If six coordinates of `m` were zero, choose those coordinates as `A`.
Then `L_A(m)=0`, contradicting (6). Thus at most five coordinates vanish,
which proves (7).

The certificate has a dimension-sensitive refinement. If the cutwise span
in (9) has dimension `r`, the same proof gives

```text
||m||_infinity
 <=584*(858+floor(r/2)).                               (21)
```

No optimality is claimed for (20) or (21).

## 5. The six-uniform orientation

Put

```text
x_i=m_i v_i.                                          (22)
```

Then

```text
sum_i x_i=0,
sum_(i in A)x_i!=0             for every |A|=6.       (23)
```

Define

```text
sigma_m(A)=sign(sum_(i in A)x_i)
           in {-,+}.                                  (24)
```

This is what **six-uniform orientation** means here: a sign on each
labelled six-subset induced by a single additive weight vector. It is not
a tournament orientation. The vertices are not runners with an oriented
edge between each pair. For the associated `6+7` cut, the seven-side
partial sum has the opposite sign because the total sum is zero.

Both signs occur. Indeed,

```text
sum_(|A|=6) sum_(i in A)x_i
 =binomial(12,5) sum_i x_i
 =0,                                                   (25)
```

and none of the `1716` summands is zero.

The orientation is highly constrained. Relabel the coordinates so that

```text
x_1<=x_2<=...<=x_13.                                  (26)
```

For

```text
A={a_1<...<a_6},       B={b_1<...<b_6},
```

write `A<=_G B` in Gale order when `a_j<=b_j` for every `j`. Then

```text
A<=_G B
  implies
sum_(i in A)x_i<=sum_(i in B)x_i.                     (27)
```

Consequently, the negative six-subsets form a Gale down-set and the
positive six-subsets form a Gale up-set. The entire `1716`-sign word is
therefore encoded by a boundary antichain rather than by arbitrary binary
data.

Replacing `m` by `-m` reverses every sign at once. This is the natural
global gauge.

## 6. Chamber graph and partial-cube metric

The orientation also has a path-sensitive realization. In

```text
E={x in R^13: sum_i x_i=0},                            (28)
```

consider the central hyperplane arrangement

```text
H_A={x in E: sum_(i in A)x_i=0},       |A|=6.         (29)
```

Its chambers are exactly the realizable sign words (24). Join two chambers
when their closures share a facet. The resulting chamber graph is a
partial cube under the sign map

```text
C |-> (sign(sum_(i in A)x_i))_(|A|=6).                (30)
```

For completeness, any chamber path between `C` and `D` must cross every
hyperplane on which their signs differ, so its length is at least their
sign Hamming distance. Conversely, choose interior points of `C` and `D`.
Their segment remains in the common open halfspace of every hyperplane on
which the signs agree and crosses every separating hyperplane. A generic
arbitrarily small perturbation separates simultaneous crossings, producing
a chamber path that crosses each separating hyperplane exactly once.
Thus graph distance equals Hamming distance.

This is the full real arrangement chamber graph in `E`. No claim is made
that the chambers reached by bounded integral certificates for one fixed
row induce an isometric subgraph.

This partial-cube carrier is different from THM-2267's static owner flag
complex. It records a global cut-sign chamber and an exact metric for
changing chambers. A future application must still construct a legitimate
dynamic path of LRC states before this metric can price owner transport.

## 7. Relation to the current LRC frontier

The theorem adds simultaneity and support, not low height.

1. THM-2164 already proves `dim_Q W_105(v)>=2`; (20) does not improve that
   height or add another rank.
2. THM-2166 retains a much smaller scalar carry and a sparse core encoding
   on a fixed defect-six split. Combining its cutwise low-carry outputs by
   (13) would destroy those individual carry bounds and the `7`-unit
   coefficient property.
3. If another argument, such as a branch of THM-2266, supplies a two-term
   relation, the relation in (6) is automatically independent of it because
   its support is at least eight. This is a nondegenerate rank-two
   presentation, but THM-2266's relation alternative is not itself a row
   exclusion.
4. The sign word in (24) may serve as a transition sidecar for the owner
   cut kernel of THM-2267. Static sign data alone does not prove that an LRC
   orbit pays positive chamber distance.

The connection and loss ledger is

```text
source:
  one genuine height-584 relation for each labelled 6+7 cut;

target:
  one integral relation crossing every cut simultaneously;

map:
  span a bounded cutwise basis, multiply all restricted cut
  functionals, and evaluate a nonzero top monomial on centered grids;

preserved:
  all 1716 labelled cut nonvanishing predicates at once;

destroyed:
  the cut-specific relation identity, low scalar carry, 7-unit
  coefficient support, and the original height 584;

new invariant:
  a nonconstant shifted six-subset sign word and its partial-cube
  chamber distance;

needed sidecar:
  a legal dynamic lift from owner/root transitions to chamber changes,
  or a use of support at least eight that cannot be replaced by the
  already known low-height rank-two theorem.                           (31)
```

Thus THM-2270 rules out the possibility that all balanced cuts can only be
crossed by mutually incompatible sparse relations. It does not rule out a
zero-safe row and does not prove LRC(14). QED.
