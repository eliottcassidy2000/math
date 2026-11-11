# Every H9 source: one-birth divisibility and its exact lifetime

**PROVED in this note, with FINITE-EXACT independent graph controls.** This
is an all-order classification for H=9. It refines a fixed-count stratum
of the metagraph; no general Chan–Pak coverage follows.

## Inheritance

The complete SCC alphabet at H<=9 is inherited from
[FAMILIES.md](../../uniform_return_metagraph_20261111/metagraph/FAMILIES.md).
The incoming extension-divisibility filtration is
[P13 / EXTENSION_CONTENT.md](../../modular_centers_collatz_pell_20261111/EXTENSION_CONTENT.md),
published at `f322658908`. It defines

```
C_k(T)=gcd{H(U): U is any k-vertex induced extension of T}.
```

Its general theorem proves C_(n-1)=1, with C_(k+1)|C_k. We inherit that
definition and credit the general theorem; the exact H9 classification
below uses an elementary two-path calculation and explicit two births.
The incoming finite census's single H9 obstruction at order7 is the
first member of the infinite family here.

## The classification

Write TT_a for the a-vertex transitive tournament, allowing a=0, and
`A -> B` for ordinal join. Every ordinary H9 tournament is uniquely either

```
TT_a -> S_i -> TT_c,        i in {1,2}, a,c>=0,             (1)
TT_a -> C3 -> TT_b -> C3 -> TT_c,  a,b,c>=0,               (2)
```

where S_1,S_2 are the two strong five-vertex H9 classes, represented by
the one-apex words OIII and OIOI. This is the inherited complete fixed-H
SCC classification, not a claim that an arbitrary tournament is one-apex.

For every source (1), C_1=1. For (2),

```
C_1(T)=1 if b=0;       C_1(T)=3 if b>=1.                  (3)
```

Every source in the latter case satisfies C_2=1. Thus its divisibility
lifetime is exactly2. All other H9 sources have lifetime1.

Consequently the ordinary generating functions by order are

```
all H9:       2z^5/(1-z)^2 + z^6/(1-z)^3,
lifetime2:    z^7/(1-z)^3,
lifetime1:    all H9 - lifetime2.                        (4)
```

For n>=7 there are exactly binom(n-5,2) H9 sources with C_1=3. This
classifies their complete growth with size, including arbitrary exterior
padding. It refines the H9 contribution to A000568 rather than computing
all A000568 classes anew.

## The exact extension polynomial

Let x_i=1 mean the new vertex v beats old vertex i. Let s_i,e_i be the
numbers of Hamilton paths starting/ending at i in T. Let K_ij count
ordered pairs of nonempty disjoint directed paths covering T, with the
first path ending at i and the second starting at j. Partition a
Hamilton path at v. Its exact count is

```
H(T+v_x)=H(T)+sum_j L_j x_j - sum_(i<j)(K_ij+K_ji)x_i x_j,
L_j=s_j-e_j+sum_i K_ij.                                 (5)
```

The endpoint cases give sum_i e_i(1-x_i)+sum_j s_j x_j;
the internal cases give sum_ij K_ij(1-x_i)x_j. These are disjoint and
exhaustive, proving (5). Integer finite differences show that C_1 is the
gcd of H(T), all L_j and all K_ij+K_ji. This is the k1 specialization of
the incoming coefficient-gcd mechanism, derived directly here.

## Proof of (3)

For (1), append I to its displayed apex word, preserving the original
five-vertex source. The respective H values become17 and19. Retain the
exterior padding by uniform joins. Since gcd(9,17)=gcd(9,19)=1, C_1=1.

Now let T have form (2). Name its first triangle A and second triangle B.
For i in A and j in B, the ordered-cover counts are

```
K_ij=1[b=0],          K_ji=9*2^b.                        (6)
```

For the forward order i then j, all vertices of A must belong to the
first path and all of B to the second. A middle singleton could belong
to neither, so the count vanishes when b>0. When b=0, each designated
triangle endpoint has exactly one path, giving1.

For the reverse order, both paths may occupy each intermediate singleton,
giving2^b assignments. The triangle entry and exit split sums are each3:
sum H(A-S)*paths_starting_i(S)=3, and the analogous end sum is3. Exterior
vertices are forced into the correct path and contribute1. This proves
(6), including independence from a,c.

If b=0, the coefficient in (5) is -(1+9)=-10, coprime to9. Hence C_1=1.

Suppose b>=1. Every s_j and e_j is divisible by3: any Hamilton path
traverses both triangle components, and prescribing an endpoint can
remove the factor3 from at most one of them. Every K_ij is divisible by3
as well. Here is the complete placement argument, so no exterior case
is being silently omitted:

- If both endpoints belong to the same triangle, the other triangle is
  wholly in one path and contributes H(C3)=3.
- If they belong to different triangles, (6) applies (or its reversal).
- If at most one endpoint belongs to a triangle, either a whole triangle
  is forced into one path (factor3), or a triangle is split between two
  paths with neither endpoint prescribed inside it (factor Q_C3(2)=12).
  A path cannot end before a skipped component and the other start after
  it; such a placement gives zero. These alternatives exhaust the
  relative SCC positions, including endpoints in exterior or middle TT.

Thus every coefficient and constant in (5) is divisible by3. For a
directed adjacent pair inside A, `K_ij+K_ji=2*H(B)=6`: its triangle
ordered-cover sum is2, the exterior pieces are forced, and B contributes3.
Therefore C_1 divides gcd(9,6)=3. Together these give C_1=3.

## A literal two-birth escape, uniform in every gap

Extend each old triangle C3=(2,4,1) to the strong four-vertex tournament
S4=(2,4,9,3), where the original triangle occupies its first three
vertices. Orient all cross-component incidences by ordinal join. This is
an actual induced extension retaining every old vertex and old edge:

```
TT_a -> S4 -> TT_b -> S4 -> TT_c.
```

Its path count is5*5=25. Appending two global sinks to the original source
instead gives9. Hence C_2 divides gcd(9,25)=1. The same two births work
for arbitrary a,b,c; no graph representative is silently substituted.
This proves the exact lifetime, not only escape from one selected prime.

The escape realizes25, not an arbitrary prescribed odd target such as
29 or31. This distinction is essential when importing the result into
the target-preserving return program.

## Reproduction and independent checks

Run `python -B extension_content.py` here. The checker builds all156
ordinary H9 family representatives of orders5 through12. It derives
the entire extension polynomial from independent subset endpoint DP and
checks all family gcds, cross-cover formulas and the coefficient6.
It makes584 separate literal extension/escape counts, including every
one-vertex cut for the two-triangle sources at orders6 and7. Each escape
also deletes its two new vertices and verifies equality with the old
matrix. `extension_content.json` retains the complete census and small
source matrices. All predicates remain enabled under Python optimization.

Connection: five small H strata -> complete SCC source families -> exact
attachment polynomial -> arithmetic lifetime. Preserved: actual source,
all named old incidences, count and common-divisor predicate. Lost by H
alone: the number and spacing of cyclic components. Sidecar: ordered SCC
word and its middle gap. Cheapest hostile: C3->C3 versus C3->K1->C3,
both H9 but C_1 respectively1 and3.
