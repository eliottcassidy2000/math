# A physical carry and the literal second-step boundary

The [splice](README.md) gives an unbounded family with profile $P=(H,E,O,g)$, and two strong children $A,B$ with profiles $3P$ and $3P+(2,0,0,0)$. The children differ by one arc. In the formulas below, a row mask for a new vertex specifies exactly which old vertices it beats. If there are two new vertices the earlier one beats the later.

There are two different ways to carry the +2 choice into the next place without losing the rooted data:

* From $B$, append a new vertex with mask 124. The new profile is $3P(B)=9P+(6,0,0,0)$.
* From $A$, append two vertices with masks $(60,60)$. The new profile is $3P(A)+(6,0,0,0)=9P+(6,0,0,0)$.

These are different tournaments, and even have different orders, but the absolute counts and all the specified rooted coordinates coincide, for **all** sizes of both inflated sites. For the base source both reach $(159,81,108,27)$. The surplus 6 is obtained either at the earlier stage as 2 and then tripled, or deposited directly at the next stage. It is a carried digit, not a new independent low digit.

There is also a pure strong zero step from $A$: append masks $(189,189)$ to reach $9P$. So the original $+2$ separation can be propagated to a $+6$ separation with the same continuation, using strong hosts on both branches. Finally from the mask-124 successor of $B$, append masks $(509,60)$ to reach $27P+(18,0,0,0)$, giving $(477,243,324,81)$ at the base size. This is a third strong step in that branch. More equivalent choices for each of the first two steps are in the verifier.

All these templates leave vertices 2 and 3 with the same unique exterior exit, respectively vertex 3 and vertex 4. Thus the same run caps used in the first proof hold; the exact coefficient tables again have degree at most two in each module's run count. The parity-sensitive two-by-two interpolation is therefore an all-size proof. As in the first splice it even extends to arbitrary modules, replacing each occurrence of the scalar 2 by $2H(M)H(N)$. Independent literal subset DP checks several actual inflations. Run **python3 check_carry.py**.

## What cannot be renewed at the literal seeds

For the six-vertex $C$ and the two eight-vertex children with the two modules left as singletons, the following classifications hold for **all orders** of the strong host. Counts in the columns are numbers of extensions labeled relative to the fixed old vertices, with the one or two new vertices in increasing orientation order. They are not counts of isomorphism classes.

| Old count | Strong hosts at $3H$ | At $3H+2$ | At $3H+4$ |
|--:|--:|--:|--:|
| 17 | 10 | 14 | 0 |
| 51 | 11 | 2 | 1 |
| 53 | 10 | 4 | 0 |

[Busch's lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) gives at least $5^{(n-1)/3}$ paths for a strong tournament on $n$ vertices. At orders 9 and 11 this already exceeds 55 and 163 respectively. Hence the displayed searches only need up to two new vertices. Their tournament is automatically transitive, so labeling them by its order loses no extension. The checker tries all orientations to the old vertices, with no scalar pruning, and computes Hamilton counts directly by subset DP.

For $3H+2$ and $3H+4$ no disconnected strong-component decomposition can supply a different host: the SCC containing the old strong source has count at least $H$, and any other nonsingleton SCC would force a multiplier of 3 (since $3H+4<5H$), contrary to these targets modulo 3. Singleton pads do not change the count. In particular the literal source cannot acquire its missing +4 sibling at any order, and the literal count-53 source cannot give 163.

The literal count-51 source does give 155 and 157, and the literal count-53 source gives 161. But **none** of these nonzero second digits has the required tripled $(E,O,g)=(81,108,27)$, even checking every root in the trimmed strong host. Within either parent, no two distinct target values in the displayed row have any matching rooted continuation among their strong hosts, even permitting arbitrary choice of root. Five count-153 placements and six count-159 placements do triple the rooted data when retaining root 4 for their respective parents 51 and 53. Adding singleton pads and taking a new root there is outside the strong-host continuation classification.

To replay:
```sh
c++ -std=c++17 -O3 check_literal_boundary.cpp -o /tmp/literal_boundary
/tmp/literal_boundary
```

These results show why the physical carry helps, but does not yet give an inductive digit system. The first two branches can meet exactly after changing their extensions, and the zero branch can continue. A fresh nonzero digit with the checked continuation does not reappear at the second step of the literal construction. Nothing here excludes changing the old arcs or finding another source family.
