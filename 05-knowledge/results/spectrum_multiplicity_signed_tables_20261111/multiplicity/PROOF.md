# Exact ordinary and labelled tournament multiplicities

Status: **PROVED** weighted SCC extension and cumulative formulas;
**FINITE-EXACT** complete orders 1–8 from the inherited class bank;
**CITED / INHERITED** global completeness of the strong inventory at path totals
at most 375. No new tournament census was run.

## Inheritance and counting convention

The fixed-path-total polynomial mechanism is already proved in
[`actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md`](../../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md),
not a new result of this packet. Its load-bearing data are
`extended_census.json` and `n8_classes.csv` in that directory. The current
canon route is **THM-4884**, slug
`definability-lenses-korselt-congruence-of-spectrum-polynomial-strong-multiplicity-census-minimal-rank-sublanguages`;
the complete order-8 atlas is **THM-4810**, slug
`h-spectrum-atlas-strong-fibres-exact-low-spectra-complement-parity-invariant-score-monotonicity-paley-budgets`.

Write $H(T)$ for the actual number of directed Hamilton paths. Define

- $N_n(h)$: ordinary isomorphism classes on $n$ vertices with $H=h$;
- $L_n(h)$: tournaments on the fixed labelled vertex set $[n]$ with $H=h$;
- $S_n(h)$, $R_n(h)$: the corresponding counts restricted to strongly
  connected tournaments, including the one-vertex tournament.

Converse classes are **not** identified. In particular, neither the number of
converse-merged representatives nor $H/|\operatorname{Aut}(T)|$ is one of these
counts. The least-used useful sidecar here is the automorphism order of every
strong class. It turns the existing unlabelled compiler into an exact labelled
compiler without another enumeration. The hostile example is $h=9$: two
asymmetric strong order-5 classes and the ordered pair of triangles carry
different automorphism weights. The corrected near miss is the old suggestion
that a finite census or growing fixed-$h$ multiplicity proves unrestricted
support; the compiler preserves the actual strong-atom existence requirement.

## 1. Exact first-component recurrence

**PROVED.** With $N_0(1)=L_0(1)=1$ and other order-zero values zero,

\[
N_n(h)=\sum_{k=1}^n\sum_{q\mid h} S_k(q)N_{n-k}(h/q),
\]
\[
L_n(h)=\sum_{k=1}^n\sum_{q\mid h}
 \binom nk R_k(q)L_{n-k}(h/q).
\]

Every tournament has a unique ordered list of strongly connected components:
all arcs from an earlier component to a later one point forward. Every Hamilton
path visits these components in that order, so $H$ is their product. Removing
the first component proves the first recurrence. For labels, choose its
$k$-element vertex set, then choose the two labelled tournaments. This gives
the second recurrence with no component permutation factor. Conversely every
term reconstructs a unique tournament. These are finite exact recurrences for
every $n,h$, once the strong inventories through order $n$ are supplied.

For each strong isomorphism class $C$, its labelled contribution is
$|C|!/a(C)$, where $a(C)=|\operatorname{Aut}(C)|$. Hence

\[
R_k(q)=\sum_{[C]:\ |C|=k,\ H(C)=q}\frac{k!}{a(C)}.
\]

## 2. Automorphism-weighted core words

A core word $W=(C_1,\ldots,C_r)$ consists of non-singleton strong
isomorphism classes. Put

\[
m(W)=\sum_i|C_i|,\quad H(W)=\prod_i H(C_i),\quad
a(W)=\prod_i a(C_i).
\]

For $h>1$, let $\mathcal W_h$ contain precisely the core words with
$H(W)=h$. The inherited strong-order bound makes this a finite set. More
basically, the bound $H(C)\ge |C|\ge3$ already gives finiteness for each fixed
$h$. Set

\[
B_{n}(W)=
\begin{cases}\binom{n-m(W)+r}{r},&n\ge m(W),\\0,&n<m(W).
\end{cases}
\]

**PROVED.** For every $n\ge1$,

\[
\boxed{N_n(h)=\sum_{W\in\mathcal W_h}B_n(W)},\qquad
\boxed{\frac{L_n(h)}{n!}=\sum_{W\in\mathcal W_h}\frac{B_n(W)}{a(W)}}.
\tag{1}
\]

Indeed, distribute the $n-m(W)$ singleton components among the $r+1$
ordered gaps of $W$. Distinct distributions have distinct ordered SCC words.
An automorphism fixes the component order, every singleton, and each strong
component setwise; its order is exactly $a(W)$. Orbit-stabilizer proves the
second identity. Repeated isomorphic strong components do not introduce a
factorial, since their positions are fixed by the condensation order.

Thus the map from ordinary core words to labelled mass preserves size,
path total, support, and polynomial degree; it destroys the integer number of
classes and needs precisely the automorphism sidecar. The $h=9$ example below
is a decisive test that replacing every class weight by one would fail.

Let $r_h=\max_{W\in\mathcal W_h}|W|$. Both $N_n(h)$ and $L_n(h)/n!$
are exact polynomials for

\[
n\ge n_0(h):=\max_{W\in\mathcal W_h}\{m(W)-|W|\},
\]

with degree $r_h$. Their leading coefficients respectively are

\[
\frac{\#\{W:|W|=r_h\}}{r_h!},\qquad
\frac1{r_h!}\sum_{|W|=r_h}\frac1{a(W)}.
\tag{2}
\]

For $m-r\le n<m$, the polynomial binomial has an integer upper argument
between 0 and $r-1$, so vanishes as required. This proves the stated onset;
it is a sufficient onset, not a claim that it is minimal. If no word exists,
both sequences vanish and no degree is assigned. For $h=1$, the empty word
gives $N_n(1)=1$, $L_n(1)=n!$.

The action of $\operatorname{Aut}(T)$ on its Hamilton paths is free: fixing
one ordered spanning path fixes all vertices. Therefore $a(T)\mid H(T)$, and

\[
\frac{N_n(h)}h\le\frac{L_n(h)}{n!}\le N_n(h).
\tag{3}
\]

This explains why both counting conventions have the same polynomial degree,
while their modal path total can differ.

## 3. Counting every order up to a vertex budget

**PROVED.** For $h>1$, the exact number of isomorphism classes on at most
$n$ vertices is

\[
\sum_{j=1}^nN_j(h)=\sum_{W\in\mathcal W_h}
 \mathbf1_{n\ge m(W)}\binom{n-m(W)+|W|+1}{|W|+1}.
\tag{4}
\]

This is the hockey-stick identity applied to (1). Its eventual degree is
$r_h+1$, with leading coefficient
$\#\{W:|W|=r_h\}/(r_h+1)!$. For $h=1$, the count is exactly $n$.

For labelled tournaments, “up to $n$” means the disjoint union of tournaments
on $[1],[2],\ldots,[n]$, not all embeddings into $[n]$. Its exact count is
$\sum_{j=1}^nL_j(h)$. For any attained fixed $h$, if
$P_h(n)=L_n(h)/n!$, then

\[
\sum_{j=1}^n L_j(h)
=n!P_h(n)\left(1+\frac1n+O_h(n^{-2})\right).
\tag{5}
\]

The term at $n-1$ has relative weight
$P_h(n-1)/(nP_h(n))=1/n+O_h(n^{-2})$; all earlier terms together have relative
weight $O_h(n^{-2})$. The latter follows by polynomial boundedness for
$j\le n$ and the factorial denominator, with finitely many initial terms
absorbed. In contrast to the ordinary count, summing labelled orders does not
increase the polynomial exponent multiplying $n!$.

## 4. Complete finite table and all-order low totals

**FINITE-EXACT.** [`spectrum_n1_n8.csv`](spectrum_n1_n8.csv) gives every odd
$h\le661$, including zeros, at each order 1–8: 2,648 rows. It includes ordinary,
labelled, strong-only, and both cumulative counts. All even totals and all
totals above 661 are absent at these orders. The table is reconstructed from
the inherited complete class inventory, with no isomorphism testing or new
tournament generation.

| Order $n$ | Ordinary classes | Labelled tournaments | Strong classes | Distinct $H$ | Maximum $H$ |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 1 | 1 |
| 2 | 1 | 2 | 0 | 1 | 1 |
| 3 | 2 | 8 | 1 | 2 | 3 |
| 4 | 4 | 64 | 1 | 3 | 5 |
| 5 | 12 | 1,024 | 6 | 7 | 15 |
| 6 | 56 | 32,768 | 35 | 19 | 45 |
| 7 | 456 | 2,097,152 | 353 | 77 | 189 |
| 8 | 6,880 | 268,435,456 | 6,008 | 320 | 661 |

At order 8, precisely the following odd values at most 661 are absent:

\[
7,21,611,615,617,619,623,625,635,647,655.
\]

Thus every odd total through 609 except 7 and 21 is present at order 8. Adding
a singleton source is an injection on isomorphism classes and preserves $H$,
so support at order $n$ equals support at orders at most $n$. Counts at those
two budgets differ; the CSV supplies both.

The mode distinction is concrete:

| $H$ at order 8 | Ordinary classes | Labelled tournaments |
|---:|---:|---:|
| 281 | 48 | 1,935,360 |
| 345 | 49 | 1,911,168 |

The ordinary mode is 345; the labelled mode is 281.

**INHERITED GLOBAL SCOPE.** The prior packet applies the cited Busch bound
$H(C)\ge5^{(|C|-1)/3}$ for a nontrivial strong tournament. Since
$5^8>73^3$, every strong tournament with $H\le73$ has at most eight
vertices. Consequently the finite bank determines (1) for every vertex count
and every $H\le73$, not merely for orders up to eight.
[`all_order_H_le73.json`](all_order_H_le73.json) records the exact core terms,
both polynomial coefficient lists, degrees, and sufficient onsets.

The broader inherited inventory recovered in Section 5 supplies the same
automorphism data for all strong classes with $H\le375$. It includes the four
order-12 classes at $H=375$; the cited next strong minimum is 625 at order 13.
Consequently [`all_order_H_le375.json`](all_order_H_le375.json) extends both
all-order multiplicity channels to every odd total at most 375. The query
program uses this larger bank. The order-8 table remains a separate complete
census, including totals between 375 and 661.

| $h$ | $N_n(h)$ | $L_n(h)/n!$ | Valid from $n$ |
|---:|---:|---:|---:|
| 3 | $n-2$ | $(n-2)/3$ | 2 |
| 5 | $n-3$ | $n-3$ | 3 |
| 9 | $(n-4)(n-1)/2$ | $(n-4)(n+31)/18$ | 4 |
| 15 | $(n-4)(n-3)$ | $(5n^2-17n-32)/15$ | 5 |
| 33 | $n^2-5n-1$ | $(n^2+7n-67)/3$ | 6 |
| 39 | $(n-6)(n-5)$ | $(n-6)(n-1)/3$ | 6 |
| 63 | $2(n-7)$ | $2(n-7)$ | 7 |

Outside the displayed onset use the truncated binomial formula (1), not an
unrestricted polynomial evaluation. In particular, the formula for $h=63$
is zero through order 7, then gives two classes at order 8. Both global strong
$H=63$ classes have order 8 and trivial automorphism group. Hence the
cumulative ordinary count through $n\ge7$ is $(n-7)(n-6)$.

For the inherited ray $h=7\cdot3^a$, $a\ge2$, the maximal-length words
contain one of these two $H=63$ atoms and $a-2$ triangles. Their length is
$a-1$, and their automorphism order is $3^{a-2}$. Thus

\[
N_n(7\cdot3^a)\sim\frac{2}{(a-2)!}n^{a-1},\qquad
L_n(7\cdot3^a)\sim
\frac{2}{3^{a-2}(a-2)!}n!n^{a-1}.
\tag{6}
\]

The factor $3^{a-2}$ is exactly the automorphism cost of the triangle
components. For comparison,
$N_n(3^a)\sim n^a/a!$ and $L_n(3^a)\sim n!n^a/(3^a a!)$.
The two absent initial members $7,21$ and the first present member 63
therefore induce a precise degree deficit along this multiplicative ray.

## 5. The first 25 primes give 23 nonzero columns

**FINITE-EXACT, global in vertex order.** The inherited
[`session173_modular_clocks_20261103/counts527/inventory_375.json`](../../session173_modular_clocks_20261103/counts527/inventory_375.json)
records the complete strong inventory through path total 375. Its completeness
route is **THM-4744**, slug
`interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve`.
For the present prime extraction only orders through nine are needed:
the inherited minimum at order 10 is 125, and subsequent minima are larger.

Of the first 25 primes, all below 100, precisely 23 are path totals: 2 is
excluded by parity and 7 is absent. All 23 already occur at order seven.
Their exact first orders are:

| First order | Prime path totals |
|---:|---|
| 3 | 3 |
| 4 | 5 |
| 5 | 11, 13 |
| 6 | 17, 19, 23, 29, 31, 37, 41, 43 |
| 7 | 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 |

[`prime_columns_lt100.csv`](prime_columns_lt100.csv) contains all 25 primes,
their minimum orders, global strong-class counts by order, both eventual
linear formulas, and their exact onsets. The corresponding
[`prime_columns_lt100.json`](prime_columns_lt100.json) records all truncated
linear terms and the inherited source hash.

A prime total permits exactly one non-singleton strong component. Thus, if
$s_{k,p}$ is its strong-class count at order $k$ and
$a_{k,p}=\sum_{C}1/|\operatorname{Aut}(C)|$, the all-order formulas are

\[
N_n(p)=\sum_k s_{k,p}(n-k+1)_+,\qquad
\frac{L_n(p)}{n!}=\sum_k a_{k,p}(n-k+1)_+,
\tag{7}
\]

where $x_+=\max(0,x)$. All extracted classes at primes $p\ge5$ are rigid.
At $p=3$ the unique triangle has automorphism order 3. This finite observation
also agrees with the parent packet's general prime-rigidity theorem.

The order-nine sidecar matters at the two largest primes:

| $p$ | Strong counts at orders 7, 8, 9 | Global strong total | $N_n(p)=L_n(p)/n!$, valid $n\ge8$ |
|---:|---:|---:|---:|
| 79 | 5, 6, 0 | 11 | $11n-72$ |
| 83 | 7, 4, 0 | 11 | $11n-70$ |
| 89 | 14, 7, 2 | 23 | $23n-149$ |
| 97 | 6, 10, 7 | 23 | $23n-162$ |

Summing all 23 attained prime columns gives, for every $n\ge8$,

\[
\boxed{\sum_{\substack{p<100\\p\ {\rm prime}}}N_n(p)=172n-1073},\qquad
\boxed{\sum_{\substack{p<100\\p\ {\rm prime}}}L_n(p)
 =n!\frac{514n-3215}{3}}.
\tag{8}
\]

The discrepancy comes entirely from the triangle column. Thus the finite
cardinality 23 comes with an exact all-size multiplicity statement. Its
agreement with the constant 23 in a separate affine construction is only a
numerical coincidence until a map preserving the relevant constructor data
is supplied.

## 6. Reproduction and exact boundary

Run from the repository root:

```bash
python3 -B 05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/multiplicity/build.py
python3 -B 05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/multiplicity/build.py --n 100 --h 39
```

The second command returns exact ordinary and labelled counts at order $n$ and
through order $n$. It deliberately
restricts the all-order query to $1\le h\le375$, the globally complete strong
inventory range. [`receipt.json`](receipt.json) records source SHA256 hashes,
scope, sample queries, and the finite checks. The script checks both SCC
recurrences against every inherited bin and checks

\[
\sum_hL_n(h)=2^{\binom n2},\qquad
\sum_hhL_n(h)=n!\,2^{(n-1)(n-2)/2}.
\]

The second identity counts each Hamilton ordering together with arbitrary
orientations of its nonconsecutive pairs. It is a global control independent
of the SCC recurrence. The script also compares every all-order formula with
all its available order-1–8 bins. It never treats a converse-merged record as
one ordinary class: inherited ordinary multiplicities are explicitly used.

The theorem converts a finite strong inventory into exact multiplicities at
all larger sizes. Proving that all odd $h\notin\{7,21\}$ occur still requires
a construction producing at least one core word for every such $h$. The
finite table, positive coefficients, and growing multiplicities alone do not
supply that missing existence statement.

The companion [RIGIDITY.md](RIGIDITY.md) supplies a precise constructor
interface: 1-ready classes are rigid, their marks are unique, and their
fixed-count fibres eventually stabilize. It proves an injective two-birth
map from these ready fibres to the ordinary fibre at $3h+2$, and separates
the two actual sources at $99\cdot2^k+23$.
