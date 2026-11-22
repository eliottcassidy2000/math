# Literal bridges over 7 and 21, and a homogeneous-port obstruction

A bridge here is a *marked tournament*: an unordered vertex pair is specified in the lower-count tournament, and reversing its arc gives the higher-count tournament. Equivalence is a directed isomorphism of the lower tournament preserving the marked pair, not separate isomorphisms of the two endpoints.

## Complete literal bridges at the two holes

| Counts | Order of reduced core | Marked isomorphism classes | Converse orbits |
|---|---:|---:|---:|
| 5 and 9 | 5 | 4 | 2 |
| 19 and 23 | 6 | 6 | 3 |

Every pair in the table can be enlarged by any number of common transitive singleton components before and after its core. This gives **all** pairs, and those numbers are unique up to the indicated marked isomorphism. In particular, the generating functions by full order are $4z^5/(1-z)^2$ and $6z^6/(1-z)^2$. Reversal of all arrows acts by pairing classes of each reduced order, not by fixing a class.

**Reduction.** The strongly connected components of any tournament are linearly ordered, and the Hamilton count factors as the product of their Hamilton counts. An arc reversal within a component leaves its outside components unchanged. If the arc instead crosses from component $i$ to component $j>i$, either two adjacent singletons swap (which leaves the count unchanged) or their entire interval becomes a single strong component: the reversed arc closes the forward chain, and internal strong connectivity of the original components gives access to each endpoint. Components outside the interval are unchanged. As each pair of displayed counts is coprime, every unchanged outside component has count 1 and hence is a singleton (any nontrivial strong component has a Hamilton cycle and at least three Hamilton paths).

The exact small-order lower bounds for strong tournaments, from [Busch](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3), are $3,5,9,15,25$ at orders 3 through 7, and increasing thereafter. If the endpoints are 5 and 9, the reduced core is at most order 5. In a cross-component reversal the merged strong side cannot have count 5: it would have order four, whereas a nonstrong tournament with count 9 requires either a strong component of order at least five or two strong components of order at least three. Therefore the split side has count 5, with its sole strong component on four vertices, and the merged component has exactly five vertices. Any same-component case can be checked through order 5. If the endpoints are 19 and 23, whichever is the split side must have a single strong component on six vertices, by primality of both counts, the bounds, and the exact count scan through five vertices. A cross-component merger would have at least seven vertices, where the minimum is 25. Therefore the reduced pair is strong on both sides and has order six. Outside singleton padding cannot introduce identifications because the core is the unique interval containing the marked edge and all nonsingleton components.

**Exact finite remainder.** In lexicographic order of pairs $(i,j)$ with $i<j$, a mask bit one specifies $i\to j$, and the marked edge has endpoints 0 and 1. Minimize the mask over all relabelings sending the unordered edge to $\{0,1\}$. The [replay](farey-small-gap.py) uses subset path DP, checks strong connectivity, reverses the marked bit, and enumerates every labeled tournament up to order six. The canonical masks and their converse pairs are
- 5 to 9: $(9,111)$, $(259,319)$; all four lower cores are nonstrong.
- 19 to 23: $(57,248)$, $(2085,2294)$, $(2099,2255)$.

No assumption that independently chosen representatives on either side are adjacent is needed.

## Why a homogeneous plug cannot give the universal signed affine step

Let a fixed outside tournament receive a nonempty substituent tournament $A$ in a module: each outside vertex has the same arc direction to every vertex of $A$. Then the resulting Hamilton count is
$$\sum_{j=1}^d a_j q_j(A)$$
where $q_j$ counts coverings by $j$ ordered, nonempty directed paths. Indeed, remove the vertices of $A$ from a full path and record its ordered outside vertices and the nonempty gaps for $A$. The outside edges between adjacent outside vertices and the contact edges to the used gaps are checked entirely on the fixed outside shell. Each accepted skeleton contributes exactly $q_j(A)$ with $j$ occupied gaps; these coefficients are nonnegative integers. Equivalently one can use allowed-empty slots, writing $q_j(A)=\sum_{k=1}^j(-1)^{j-k}\binom jk Q_k(A)$, where $Q_k$ is the allowed-empty count, retaining no $Q_0$ term for nonempty $A$.

Suppose the shell gives an answer depending only on $H(A)$ on a class closed under ordinally prepending singletons. For any one $A$ in the class, appending $m$ such singleton positions at the start holds $Q_1=H(A)$ constant and replaces $Q_j$ by $j^m Q_j(A)$, since each singleton can independently go in front of any allowed-empty slot. As $Q_j(A)>0$, the linear independence of $1^m,2^m,\ldots,d^m$ forces all coefficients of $Q_j$ with $j\ge2$ to vanish. The answer is thus a fixed multiple of $H(A)$. This is true for each orientation of any marked edge confined to the fixed outside shell; it remains true if the two fixed homogeneous contact patterns differ.

Consequently fixed homogeneous shells cannot implement even one of $2H(A)-1$ and $2H(A)+1$ on arbitrary tournaments, or on the class of all prime-Hamilton-count tournaments (which is closed under singleton prepending and has distinct counts 3 and 5). Nor can they yield an unchanging gap of two while the input count varies. A universal signed step must therefore expose internal endpoint structure, constrain the input order or choose representatives, or change its contact pattern as the input varies. The finite hole bridges are compatible with this obstruction: their roles are marked fixed vertices, not arbitrary homogeneous plugs.
