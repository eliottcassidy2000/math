# Prime induced jumps of twelve and the 13 gluing row

**Theorem.** Among prime sources $q$ at least 11, an induced increase $q\to q+12$ occurs exactly for $q\in\{11,13,17,19\}$. For every such increase, after deleting the host's singleton strong components, the host is the old source's unique nontrivial strong core with one vertex added. A host of count $q+12$ can contain induced witnesses of counts $q$ and $13$ if and only if $q=13$.

In the positive self-pair the order-six tournament of orientation code 1105 gives the exact identity $13+13-1=25$, with two distinct deletions giving 13. Here code 1105 uses bits $i\to j$, then its reverse if zero, ordered by increasing $j$ and within that by increasing $i<j$. The theorem only concerns **induced** gluing and does not exclude a rewired construction of $q+12$.

## Localization

Use the [positive odd-cycle-packing expansion](https://arxiv.org/abs/2307.05569) of the Hamilton-path count and the repo's [triangle budget and small-gap classifications](https://github.com/eliottcassidy2000/math/blob/c3cf550db9ba849e8d1dfb9cf2c6eefc61b550e5/05-knowledge/results/triangular_backbone_chanpak_20261116/CORE.md). There is no positive induced increase of 2, 4 or 6 from a count at least 11. Insert the additional vertices one at a time, starting with the old prime witness. There is exactly one strict step, since two would cost at least 16. Neutral additions create no directed triangle; they only add singleton components. The old prime witness has a unique nontrivial strong component $C$, of count $q$. At the strict step with vertex $x$, the triangle budget bounds the number of triangles in $C$ by six. The strong order bound of $|C|-2$ triangles gives $|C|\le8$.

As in the [11-row theorem](../eleven_row/README.md), uniform contact from $x$ to $C$ preserves the module. A strict step then creates a triangle in the quotient, whose paths each lift to $q$ paths, costing at least $2q>12$. Hence contact is mixed, and there is a triangle through $x$ in $C+x$. Let $d=(H(C+x)-q)/2$. In the positive expansion the contributions to the strict half-increase from cycles confined to $C+x$ sum to $d$. Any vertices outside this subgraph were singleton pads just before $x$ arrived.

## An exterior triangle bound

Let $E_-$ be the sum of $H(C-z)$ over old vertices $z$ that beat $x$, and $E_+$ the sum over old vertices $z$ that $x$ beats. In the ordinal order of the pads and $C$, write $b_-$ for the number of pads before $C$ that $x$ beats, and $b_+$ for the number after $C$ that beat $x$. Write $b_0$ for the number of pairs of pads $y_i,y_j$ in that order with $x\to y_i$ and $y_j\to x$. At the strict step, writing $r$ for its full half-increase, we have
$$
r-d\ \ge\ b_- E_- + b_+ E_+ + b_0 q.
$$

Indeed any counted earlier pad makes a triangle with $x$ and **every** $z$ of the incoming side in $C$; the later pads use the outgoing side instead. Their positive packing weights are $H(C-z)$, since the other old parts are ordinal singleton pads. A counted pair of pads yields a distinct triangle with weight $q$. These are disjoint *families of rooted triangles*, so their contributions add even when individual triangles share vertices. The right-hand side is exactly the rooted triangle contribution from outside $C+x$, before any longer cycles are charged.

If any positive exterior increment occurs, at least one such triangle must occur: otherwise all triangles in the extension are in $C+x$, so the other vertices remain singleton components, and no directed cycle can visit them. Thus $r-d$ is either zero, or at least $\min(E_-,E_+,q)$, independent of the number of available outside pads. This filter is reusable for other bounded-gap rows.

## The finite certificate

Relabel every core along a Hamilton path. Keeping the edges $i\to i+1$ fixed leaves just 21 free bits at order eight, without losing any isomorphism class. For each core with at most six triangles, compute path counts by direct subset DP, then test every mixed contact using the same DP on the augmented graph. The following table lists all prime core counts at least 11 passing the triangle bound:

| $q$ | 11 | 13 | 17 | 19 | 23 | 29 | 31 | 37 | 41 | 53 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Minimum local half-gap $d$ | 4 | 5 | 5 | 6 | 7 | 8 | 10 | 16 | 14 | 14 |

The only eligible contacts with $d<6$ are $(q,d)=(11,4),(13,5),(17,5)$. On these three groups the minima of the exterior threshold are 5, 5, and 9, respectively. The slack $6-d$ is at most two, so none works. Thus every actual +12 increase comes from a local contact of half-gap six, using no exterior increment. All later additions are neutral, so the host after trimming singleton components is exactly $C+x$.

There are 66, 52, 68, and 38 exact contacts in this labeling convention for $q=11,13,17,19$ respectively. These are numbers of encoded core/contact pairs, not isomorphism classes. Every induced subset of these cores was counted. The rows for $q=11,17,19$ contain no induced count 13. For existence of the +12 steps, first exact-contact witnesses $(n,\text{core-code},\text{contact-mask})$ are $(5,5,6),(5,2,4),(6,0,2),(6,5,4)$ in the same order. The core-code convention omits the fixed consecutive arcs; the mask records the vertices beaten by $x$. The separate code-1105 witness demonstrates the two distinct same-count deletions for $q=13$.

Run the exact certificate from this folder:
```sh
c++ -O3 -std=c++17 check_thirteen.cpp -o /tmp/check_thirteen
/tmp/check_thirteen
```
