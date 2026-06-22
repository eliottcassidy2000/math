# u(21) = 57 is proven; the anatomy of the certificate (S711)

The question was whether the maximum number of unit distances among 21 points is 57, and whether we
can prove it. The answer is yes on both counts, but the second yes needs an honest qualifier.

**It is proven, not conjectured.** Alexeev, Mixon and Parshall (and the closely related Alexeev-
Tikhonov note, arXiv 2412.11914) determine `u(n)` exactly through `n = 21`. The verified sequence is

```
u(n) = 0,1,3,5,7,9,12,14,18,20,23,27,30,33,37,41,43,46,50,54,57, ...   (n = 1..21)
```

so `u(20) = 54`, `u(21) = 57`, and `n = 22` is the first open case, `60 <= u(22) <= 61`. The cluster
already held this (S614 reflection) and an explicit 21-point certificate ("AMP extremal core 1"); this
session pins down what that certificate *is*.

**The lower bound is a certificate we can check by hand, and I did.** The 21 stored coordinates realize
exactly 57 pairwise distances equal to 1, to machine precision (max deviation `2.2e-16`), with no
coincident points and the next-closest distance a full `0.24` away from 1. So `u(21) >= 57` is not a
delicate near-miss; it is a robust, faithful unit-distance graph.

**The proof of the matching upper bound is genuinely computer-assisted, and that is the qualifier on
"can we prove it."** The hard half is `u(21) <= 57`: show no 21-point set has 58 unit distances. AMP do
this in two stages. Combinatorially, they generate all graphs avoiding a set `F` of forbidden
subgraphs — structures no unit-distance graph can contain — efficiently enough to enumerate every
candidate with too many edges. Algebraically, they decide embeddability with a custom solver faster in
practice than cylindrical algebraic decomposition, and kill every over-dense candidate that the graph
bound left standing. For `n <= 21` the combinatorial ceiling already meets the construction, so the
enumeration closes the gap. This is not a one-page argument; it is a verified computation. We can
*reproduce and check the certificate*, and we can *understand the mechanism*, but "prove it from
scratch by hand" is not on the table for this `n` — the result lives in the same computer-assisted
regime as the four-color theorem.

**The anatomy of the extremal graph (new here).** Building the abstract graph on the 57 unit edges and
measuring it:

- degree sequence: **eighteen 5's and three 8's**, sum `114 = 2 * 57`, average `5.43`;
- clique number `omega = 3` — forced, since `K_4` is not a unit-distance graph (no four points are
  mutually at distance 1 in the plane);
- **`K_{2,3}`-free** (zero copies) — exactly the forbidden bipartite structure the faithful-embedding
  bound leans on; a faithful UDG cannot contain it, and the extremal graph saturates right up against
  that wall;
- 25 triangles, 54 four-cycles.

The three degree-8 vertices are the signature. A degree-8 vertex has eight neighbors all at distance 1
from it — eight points on a unit circle about it. Three such hubs means three overlapping unit-circle
pencils, the Moser-ring motif that makes these small extremal configurations dense: you reuse the same
points as circle-points for several centers at once. The eighteen degree-5 vertices are the shared rim
doing double and triple duty. The whole object is the densest way to make centers and rims overlap
without ever completing a `K_4` or a `K_{2,3}` — the extremal graph is the tightest `{K_4, K_{2,3},
...}`-free graph on 21 vertices that still embeds.

**Connection to the cluster's frontier.** This is the same lesson S614 drew for `u(22)` and for LRC
`n=14`: the visible graph is not the master object; the carrier is the graph *plus* the embeddability
probe. The `u(21)` case is the last `n` where the cheap graph-only ceiling and the geometry agree, so
the enumeration suffices. At `u(22)` they diverge — the `F`-free bound gives 62, but the two 62-edge
graphs carry totally-unfaithful subgraphs, dropping the truth to `<= 61`, and the one open bit is
whether an embeddable 61-edge graph exists. The degree structure found here sharpens the S614 deletion
argument concretely: a 61-edge UDG on 22 vertices, minus a min-degree vertex, leaves `>= 56` edges on
21 vertices, i.e. a near-extremal 21-core that must look like a one-vertex perturbation of *this*
eighteen-5's-three-8's graph (or its few siblings). The hunt for `u(22)` is a hunt for embeddable
degree-4-or-5 extensions of the `u(21)` extremal deck.

Sources: Alexeev-Mixon-Parshall and Alexeev-Tikhonov, "The Erdős unit distance problem for small point
sets," arXiv:2412.11914; Engel-Hammond-Lee-Su-Varga-Zsamboki, arXiv:2406.15317.
