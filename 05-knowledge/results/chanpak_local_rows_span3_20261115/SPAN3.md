# Full-alphabet span-three floor bands factor exactly

Write a span-3 ribbon in its given order: backward arcs go back by at most three positions. All eight letters are allowed. At the start, instructions referring to nonexistent vertices are ignored. Write $H$ for Hamiltonian path count and $f(n)$ for the Busch minimum among strong tournaments of order $n$.

**Proved.** Any contiguous interval of length at least 20 whose induced tournament is strong and at the Busch floor is a rank-one reset in the full minimal span-3 transfer, regardless of its internal letters or the two ambient sides. If its order is $m$, it forces divisibility of the whole count by $5^{\max(0,\lceil(m-27)/3\rceil)}$. In particular, prime words cannot contain such an interval. Independently, all five length-two trace-3 and thirteen length-three trace-5 cyclic reset classes are tabulated with their exact phase exponents below. The all-alphabet growth bound $1.7699^n$ is not claimed; §4 gives a scoped exponential premium and a counterexample to a tempting covering shortcut.

## 1. A floor interval has an intrinsic safe junction

The [Moon–Yang equality characterization](https://citeseerx.ist.psu.edu/document?doi=78c085212810b466efaf3f903baeb284a695f719&repid=rep1&type=pdf) says that a strong Busch-floor tournament is a chain of nearly transitive links of orders 3, 4, or 5. In the intrinsic chain order all arcs go forward except for one backward bottom-to-top arc in each link. Successive links share a junction. The equality floor budget permits at most two links of order 3 or one of order 5, not both. Thus each vertex participates in at most one exceptional arc coming from a later vertex and one to an earlier vertex.

Label the physical positions inside the interval $U$ from 0 to $m-1$. For a vertex $x$, an intrinsically later vertex appearing before it in the ribbon either supplies the unique exceptional arc toward $x$ or must lie within three physical positions of $x$ (otherwise their intrinsically forward arc would be a forbidden long backward arc). Thus at most four intrinsically later vertices appear before $x$. Dually at most four intrinsically earlier vertices appear after it. The physical position minus the intrinsic rank is the difference of these numbers and has absolute value at most four.

Choose an intrinsic junction $t$ of rank in $[8,m-9]$. For $m\ge20$ the interval contains at least four consecutive ranks, and junction gaps are at most four. For arbitrary ambient ribbon prefix $P$ and suffix $S$, split $PUS$ into induced factors $A,B$ attached at $t$: put $P$ and the intrinsic floor prefix through $t$ in $A$, and the floor suffix from $t$ and $S$ in $B$.

All cross arcs $A-t\to B-t$ point forward in that direction. For the intrinsic pieces this is the chain property. The earliest possible physical position in $U$ of a vertex in the right piece is $t+1-4\ge5$, so every $P$ vertex is too distant for a backward arc; the left piece versus $S$ is dual. The outside words are separated by $U$. Also $t$'s arcs to the outside go forward in the expected directions, since its physical position is at least $t-4\ge4$ and at most $t+4\le m-5$.

The vertex of intrinsic rank $t-1$ is an interior vertex of the preceding link. Every other vertex of the intrinsic $A-t$ beats it. It is beaten by all of $P$ as well, since its physical position inside $U$ is at least $t-5\ge3$, at distance at least four from the last prefix vertex. Hence it is a sink of $A-t$; dually rank $t+1$ is a source of $B-t$. Only these nearest interiors need to be insulated from outside effects.

## 2. Exact product and rank-one consequences

The [one-way seam formula](../chanpak_oneway_seams_20261115/README.md) is
$$
 H(PUS)=H(A)H(B)+a(H(B)-S_B)+b(H(A)-E_A)+a_+b_+.
$$
Here $a,a_+$ count bad nonempty $A$-prefixes ending in an outneighbor of $t$, for which the other path of the two-path cover is empty or a valid suffix; $b,b_+$ are the dual bad-suffix cases. $S_B,E_A$ count paths starting and ending at $t$. This formula depends on the one-way cross arcs, not on strength of the factors.

Within $A$, $t$ has exactly one outneighbor, the preceding link's top. It cannot end a spanning path of $A-t$ since rank $t-1$ is a different vertex and a sink. It also cannot both end a bad prefix and start a valid nonempty suffix in a vertex-disjoint cover. Thus $a=a_+=0$, and dually $b=b_+=0$ for $B$. This gives the universal separated product
$$
 H(PUS)=H(A(P))H(B(S)).
$$
Each factor contains an intact adjacent link, hence a directed triangle. Its count is at least three by deterministic Rédei insertion of the other vertices. A prime count is therefore impossible, even if the floor band touches one or both outer endpoints.

For the matrix consequence, use the full minimal integer representation of [THM-4880](../../01-canon/theorems/THM-4880-ribbon-tournaments-bounded-span-back-arcs-give-a-finite-state-transfer-for-hamiltonian-paths-strong-connectivity-is-regular-the-span-2-series-is-an-irreducible-14-dimensional-integer-representation-47-needs-span-3.md). Minimality means that prefix vectors $M_P\mathrm{init}$ span the state space and suffix observables $\mathrm{fin}M_S$ span its dual. Choose bases among them. Their evaluations $\mathrm{fin}M_S M_U M_P\mathrm{init}$ form a rank-one matrix by the separated product, so $\operatorname{rank}M_U=1$. Initial missing bits pose no restriction: the truncation to the induced interval is independent of the ambient sides. Positive counts exclude rank zero.

Apply the exact rule at the first and last safe junctions to split off the intact intervening floor chain. Its count is the product of its link counts. The first safe junction has rank at most 11 and the last at least $m-12$, so their intrinsic separation is at least $m-23$. At most four order-cost units can come from non-order-4 links, by the floor budget. Hence the intervening chain has at least $\max(0,\lceil(m-27)/3\rceil)$ links of order 4, each contributing a factor 5. If only one safe junction exists, the exponent in the claimed bound is zero.

## 3. Exact short reset periods in the full representation

This calculation is separate from the floor-band proof. For each word of length two up to cyclic rotation, test eventual rank one with sole nonzero eigenvalue $\lambda=3$; for length three use $\lambda=5$. A necessary condition is $\operatorname{tr}(M^j)=\lambda^j$ for every positive integer $j$. The replay rejects all omitted classes by a trace-power mismatch and checks the listed classes by the exact pivot identity for rank one together with $M^{e+1}=\lambda M^e$. Entries in each row are the minimal exponents $e$ for successive left cyclic rotations. Letter 0 is empty, and digits inside one letter are the backward arc spans.

| period-2 class | exponents |
|---|---|
| 0.2 | 3,4 |
| 0.12 | 4,4 |
| 0.123 | 5,5 |
| 1.2 | 4,4 |
| 2.13 | 5,5 |

| period-3 class | exponents |
|---|---|
| 0.0.3 | 3,2,3 |
| 0.0.13 | 2,3,2 |
| 0.0.123 | 3,3,3 |
| 0.1.3 | 3,4,3 |
| 0.1.23 | 2,3,2 |
| 0.1.123 | 3,4,3 |
| 0.12.3 | 3,3,2 |
| 0.3.1 | 3,2,3 |
| 0.3.12 | 3,4,4 |
| 0.13.2 | 3,4,4 |
| 0.13.12 | 2,3,3 |
| 1.2.3 | 4,3,3 |
| 1.12.3 | 4,4,3 |

In each phase $U$, a factor $U^{e+1}$ forces divisibility of the whole count by $\lambda$, even at an endpoint, by the integer matrix identity. Counts coprime to 15 avoid all such factors. THM-5202 already supplied a subset of the period-3 rows and the transfer factorization; the extension here is the full phase table and the independent floor-band criterion. The calculation does not require a short-period class to be bulk strong.

## 4. What this gives toward full-alphabet growth

A scoped all-length growth statement follows. Suppose a prime-count word admits a partition into consecutive strong intervals of orders between 20 and fixed $R$, meeting successively in exactly one endpoint, with only boundedly many end vertices omitted. By §2 every such interval is nonfloor, so contributes at least $f(m)+2$ (positive excess is even). Busch's two-path merge injects products of counts across a shared single vertex, even with arbitrary cross orientations. Hence
$$
 H(T)\ge C_R\,\Big(5^{1/3}(1+2/f(R))^{1/(R-1)}\Big)^n ,
$$
where $C_R$ may also depend on the bounded endpoint budget. Indeed $f(m)\ge5^{(m-1)/3}$; the sum of the costs $m-1$ is the covered order minus one, and there are at least this sum divided by $R-1$ intervals. This is an exponential premium over the Busch rate in the indicated language, not the constant 1.7699.

The covering premise is **not** automatic. At each even order $n\ge4$, with vertices numbered from zero, use just the backward arcs $j\to j-3$ for odd $j\ge3$. Its consecutive induced windows $[0,3],[2,5],\ldots,[n-4,n-1]$ are nearly transitive and strong; they overlap, making the full ribbon strong. Every backward-arc tail is odd and every head even. No vertex can be the shared boundary of two nontrivial strong consecutive intervals: strength on the left needs an arc from that vertex back to the left, while strength on the right requires a backward arc into it from the right. Thus proving global growth still requires larger shared interfaces or a transfer-cone argument; the one-vertex cover shortcut fails even at fixed span.

## Replay

Run python verify_span3.py PATH_TO_SPAN3_INTREP, using the existing repository file 04-computation/experiments/ribbon_spectrum_20261104/span3_intrep.txt. The phase table is checked exactly in the existing rank-62 representation on 212 cyclic classes (every omitted survivor of the first trace check fails at the second trace). There are 25 exact reordered floor-band controls with arbitrary endpoint perturbations, five checking nontrivial forced powers of 5. The structural argument and the obstruction in §4 are all-length proofs independent of the experiments. No density-one or unrestricted Chan–Pak claim follows.
