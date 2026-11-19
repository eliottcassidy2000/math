# Star-seam inverse: defect coverage through 4096 and the one-arm core bound

**Proved:** every admissible odd target of canonical defect at most
4096 occurs. The [star product](STAR_SEAM_PRODUCT.md) is unchanged;
one extra vertex on either arm, and a few cases with two, supply the
profiles missing from the previous certificate. Also proved below:
for the balanced one-arm selector, a fixed underlying Hamilton
count cannot produce arbitrarily large odd slopes even with
arbitrary singleton padding.

## The next defect layer really closes

For an arm use the four-tuple $(h,k,e,o)$ defined in the star
product note. The target intercept and slope are
$$
I=h_A(k_B-2e_B)+h_B(k_A-h_A),\qquad
c=(h_A+e_A+o_A)(e_B+o_B).
$$
The first odd part the old bank did not supply was 3697, at
defect 2056. The following actual tagged arms now give it:

| Arm | Outgoing masks | Tag mask | $(h,k,e,o)$ |
|---|---|---:|---|
| A | 132,13,368,133,43,75,283,118,187 | 90 | (1785,8833,6058,6945) |
| B | 4,5,72,67,15,31,179,63 | 131 | (97,575,479,545) |

The intercept is $-683655+683656=1$. The two factors
are 14788 and 1024, yielding $c=3697\cdot2^{12}$.
Any earlier terms on this ray are smaller than $c+2$.

The [extension bank](star_seam_bank_4096.tsv) adds another 1,421
literal tagged pairs, bringing the total to 2,842. The replay
checks the modules by exact subset dynamic programming and reruns
the complete canonical prefix classification. At defect $D$
the final zero run ends with marked count at most $2(D+1)$,
while the marked count after $k$ letters is at least $2^k+1$.
Thus at $D\le4096$ there are at most 13 letters before the
final all-one suffix. Each checked odd part is realized by one
uncapped star ray. The only terms before these rays begin are at most
62,660,609, covered by the
[corrected smaller-height certificate](https://github.com/eliottcassidy2000/math/blob/main/session203/coverage-recertification/README.md),
apart from 7 and 21. The count 1 and all-one words are immediate.
These facts extend the **all-height defect certificate**, not an
unconditional cutoff in the scalar count. Any new missing odd
count has defect at least 4097; Chan–Pak remains open.

## The one-arm reduction and what it cannot do

Allow $A$ to be empty, or a singleton with no tags. Put
$P=e_B+o_B$ and let $J$ count nonempty even-order paths on
subsets of $B$ with **both** endpoints tagged, each weighted by
the Hamilton count of the complementary subtournament.
An odd cycle through the attached apex consists of an even-order
outside path starting tagged and ending untagged. By odd-cycle
packing at weight 2, this gives $k_B=h_B+2(e_B-J)$.
The exact star output is therefore
$$
H= h_B-2J+P2^m.
$$
In this slice intercept one is **equivalent** to
$J=(h_B-1)/2$. As literal controls, masks (24,1,3,6,14) with
tags 6, 12, 24 give $h_B=11,\ J=5$ and slopes $29,37,47$
respectively. The script counts $J$ independently of the apex.

Here is an all-order bound on this slice. Fix $h=H(B)$.
The strongly connected components of a tournament are ordinally
ordered. Call each singleton component padding. There cannot be
two tagged padding vertices: their two-vertex path alone would
contribute $h$ to $J$, exceeding $(h-1)/2$.
If $h>1$, balance implies $J>0$ and there are at least two tags.
Let $k$ be the number of padding vertices strictly between the
first and last tagged components or vertices in the component order.
When $k>0$ the two extremal tags lie in different components.
Use those two tags and any even subset of the $k$ padding vertices
to form an even-order directed path. Its complement has at least
one path. Hence $J\ge2^{k-1}$; the total number of internal
padding vertices, not merely each gap, is at most
$\lfloor\log_2(h-1)\rfloor$.

Strip any untagged padding before every tag or after every tag.
Stripping before the first tag leaves $P,J$ unchanged: a path
starting tagged cannot go back to that vertex and the complementary
path count does not change. Each stripped vertex after all tags
divides $P$ by two, since for every path starting tagged that
vertex can independently be included or omitted in its unique
ordinal position; it leaves $J$ unchanged because double-tagged
paths end before it. This argument applies even when such
padding is between nontrivial components. The normalized tournament
therefore has at most $1+\lfloor\log_2(h-1)\rfloor$ singleton
vertices and exactly the original nontrivial components.

The number of nontrivial components is at most $\log_3 h$.
Their total order is at most $\log_3 h+3\log_5 h$ by
[Busch's strong-tournament bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
and multiplicativity over components. For $h=1$ the source is
transitive, with at most one tag and dyadic nonzero $P$.
Thus for each $h$ one can enumerate a finite set of normalized
sources accounting for **all** possible nonzero odd parts of $P$
in this balanced one-arm slice. In particular, since on $n$
vertices $P\le n\cdot n!$ (split a permutation after a nonempty
prefix), any large odd part $q$ forces
$\log q=O(\log h\log\log h)$. A one-arm selector needs genuinely
larger cores, not just longer singleton padding.

This is a constraint on one slice, not on the two-arm star
product. In fact the exact pair for 3697 puts its prime factor
on **A**, while B supplies a power of two, avoiding the one-arm
equation. The remaining problem is a target-dependent module
selection rule for unbounded defects, not more padding.

Run python3 star_4096_replay.py in the result directory.
