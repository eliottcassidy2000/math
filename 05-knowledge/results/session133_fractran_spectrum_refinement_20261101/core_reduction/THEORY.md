# Exact singleton-SCC core floor and bounded-H kernel

**Status:** PROVED from the cited Moon–Busch strong minimum; independent audit
requested. The finite checker is corroboration, not the all-order argument.
No external novelty claim. The permanent omissions 7 and 21 and effective
finite decidability are inherited results, not new resolutions here.

## 1. Object, inheritance, and primary source

`H(T)` counts directed Hamiltonian vertex words of a finite tournament, not
path orbits under automorphisms. Set `H(empty)=1`. Write `C(T)` for the induced
subtournament on all vertices in nonsingleton strongly connected components
(SCCs), and `c(T)=|C(T)|`. Such a tournament is called a **core** here; it may
be reducible, and this is not graph-homomorphism terminology.

Inherited sources at repository commit `335577ae6d58`:

* `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`: all-order
  omissions, and its explicit retraction of a novelty claim.
* `05-knowledge/results/rotation_owners_spectrum_20261020e/spectrum.md`,
  sections 3–4: exact strong floor, singleton-SCC deletion,
  `c(T)<=floor(3 log_3 H(T))`, equality for triangle chains, and
  `2^{O((log h)^2)}` finite spectrum membership algorithm.
* `05-knowledge/results/session124_tournament_cut_response_20261031/README.md`,
  “Equality and strong components”: precisely which one-vertex extensions
  preserve H. Our canonical induced-kernel statement below is a convenient
  equivalent consequence, not an assertion that strictness was unknown.
* `01-canon/theorems/THM-002-ocf.md`: positive odd-cycle collection formula.

Primary source independently opened: Arthur H. Busch, *A Note on the Number
of Hamiltonian Paths in Strong Tournaments*, Electronic Journal of
Combinatorics 13 (2006), N3,
[Theorem 1 and Corollary 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf).
For strong tournaments of order n>=3 the exact minimum is

```
f(3k)   = 3*5^(k-1),
f(3k+1) = 5^k,
f(3k+2) = 9*5^(k-1),            k>=1.
```

Busch is the sole author. The paper also states SCC multiplicativity. The
strong minimum does not assert completeness of any H spectrum.

The odd-cycle formula is independently confirmed in Irving–Omar,
[*Revisiting the Rédei–Berge Symmetric Functions via Matrix Algebra*,
Corollary 20](https://arxiv.org/html/2412.10572v3#S3.SS4), attributed there to
Grinberg–Stanley. The displayed complement has the same Hamiltonian count
as the tournament by path reversal. Collections are of distinct directed
odd cycles of length at least three, pairwise vertex-disjoint, with weight
`2^(number of cycles)`; the empty collection has weight one.

## 2. The canonical induced deletion kernel

**Lemma.** For every vertex subset S,

```
H(T[S]) = H(T)  iff  V(C(T)) is contained in S.
```

**Proof.** SCCs are totally ordered by dominance; a Hamiltonian path
concatenates paths of the SCCs in that order, giving
`H(T)=product H(component)`. Singleton factors are one. Thus deletion of
any collection of singleton SCCs preserves H.

Conversely, every vertex v in a nonsingleton SCC lies in a directed
triangle. To see this, choose a shortest directed cycle through v. If its
length were at least four, the chord joining v to the second vertex after
v, in either orientation, would give a shorter cycle through v. Thus a
shortest such cycle has length three. If S omits v, the positive odd-cycle
collection sum loses this triangle's one-cycle collection (weight two),
while every collection on S is still one on V(T). Hence H(T[S])<H(T).
This also proves `H(T)-H(T-v)>=2` for every nonsingleton-SCC vertex. QED.

Consequently C(T) is the unique inclusion-minimal **induced**
H-preserving subtournament, with the empty convention for H=1. This is
not a smallest-order witness among unrelated tournaments: two ordered
triangles give H=9 on six vertices, while a strong order-five witness
also has H=9.

**Hostile to twin/module contraction.** Start with a directed triangle and
replace one vertex by a transitive module of size t>=1. Its count is
`2^t+1`: this is the nearly transitive tournament on t+2 vertices obtained
by reversing the extreme arc of a total order, as in Busch's proof. A
path avoiding the reversed arc is unique; paths using it correspond to
subsets of the t internal vertices assigned before/after that arc.
Contracting the module back to a vertex gives H=3. Thus even transitive
directed twins cannot generally be contracted while preserving H.
The canonical kernel avoids exactly that loss; scalar H alone does not
remember module multiplicity.

There is no corresponding bound on the total tournament order: an
arbitrarily long transitive prefix joined before C3 still has H=3. The
bounded object must discard singleton SCCs or retain their multiplicities
separately.

## 3. Exact maximum core size under a path-count budget

For every integer M>=1 define

```
K(M) = max { c(T) : H(T)<=M }.
```

**Theorem A.**

```
K(M) = max { 3a+4b : a,b>=0, 3^a*5^b<=M }
     = max_{b in {0,1,2}, 5^b<=M}
           (4b + 3*floor(log_3(M/5^b))).                    (A)
```

Compute the logarithmic floors with integer powers; no floating-point
rounding is needed.

**Proof.** Delete singleton SCCs, preserving H. For each remaining strong
component of order n, replace it by an order join of the following small
strong components, using Busch's lower bound:

| Original size | Replacement | New size | Replacement H <= original H |
|---|---|---:|---|
| 3k | k directed triangles | 3k | `3^k <= 3*5^(k-1)` |
| 3k+1 | k-1 triangles and one strong order-four block | 3k+1 | `5*3^(k-1) <= 5^k` |
| 3k+2 | k+1 triangles | 3k+3 | `3^(k+1) <= 9*5^(k-1)` |

All rows have k>=1. A strong order-four block with H=5 is the nearly
transitive order-four tournament. Multiplicativity shows that the
resulting core has at least as many vertices and no larger H. Hence its
3/4-block budget bounds the size of the original core from above.
Conversely, every a,b is attained by the indicated ordered join, proving
the first equality. Replacing three order-four blocks by four triangles
preserves twelve vertices and lowers their product from 125 to 81.
Repeatedly make this exchange to obtain b<=2. For each such b, maximize a
by the displayed integer logarithmic floor. QED.

**Important quantifier.** K(M) is the maximum reduced-core size over
`H<=M`, not necessarily over `H=M`. Replacements need not preserve the
original H. For example K(21)=7 is attained with H=15; H=21 itself is
impossible.

Selected values:

| M | 1 | 3 | 5 | 7 | 9 | 15 | 21 | 25 | 27 | 45 | 75 | 81 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| K(M) | 0 | 3 | 4 | 4 | 6 | 7 | 7 | 8 | 9 | 10 | 11 | 12 |

Thus the inherited log bound is tightened from 5 to 4 at M=7 and from 8
to 7 at M=21. These improvements concern deletion-only witness budgets;
the old strong bound already reduces H=21 to order six for one SCC.

## 4. Exact minimum H for a specified core order

Let g(N) be the minimum H over cores with exactly N vertices. Orders one
and two are infeasible. Define g(0)=1.

**Theorem B.**

```
g(3)=3, g(4)=5, g(5)=9,

g(3k)   = 3^k,
g(3k+1) = 5*3^(k-1),
g(3k+2) = 25*3^(k-2),                                  (B)
```

where the last three formulas apply to N>=6 (so k>=2 in the last row).
In particular `g(N+3)=3g(N)` for every N>=6. At N=3 onward the sequence is

```
3,5,9,9,15,25,27,45,75,81,135,225,243,405,675,... .
```

**Proof.** For N=3,4,5 there is only one nonsingleton SCC, so Busch gives
the values. For N>=6, every N can be written 3a+4b. The minimum H among
3/4-block joins has b<=2 by the 125→81 exchange. Since b=N mod3, this
minimum is exactly the corresponding value in (B).

Call that block minimum q(N). It is strictly increasing on N>=6: the
three successive ratios, beginning at N=3k, are 5/3, 5/3, and 27/25.
Apply the replacement in Theorem A to any core of order N>=6. The output
has order N'>=N and count H'<=H. Hence
`H >= H' >= q(N') >= q(N)`.
The block joins attain q(N), proving (B). QED.

**Equality classification for N>=6.** Every minimizing core is an ordered
join of a directed triangles and b strong order-four tournaments, where

```
b = N mod3 in {0,1,2},   a=(N-4b)/3.
```

Indeed every component of order at least six incurs a strict decrease in
the replacement table, preventing equality. An order-five component is
replaced by two triangles at no larger cost and increases N by one; in
this case N>=8, where q is strictly increasing, again preventing equality.
Finally three order-four blocks cannot occur at equality because 125>81.
Conversely the indicated joins attain the floor.

There is one strong order-four isomorphism class: its score sequence is
(1,1,2,2); naming the arc between the two score-one vertices fixes the
remaining orientations up to swapping labels. Thus the number of
unlabelled minimizing cores is

```
e(N)=binomial(a+b,b),       N>=6,
e(3k)=1, e(3k+1)=k, e(3k+2)=binomial(k,2).
```

Different orders of SCC types are nonisomorphic; SCC order is invariant.
For labelled vertex set [N], the number is
`e(N)*N!/3^a`: each triangle has three automorphisms, each strong
order-four block has one, and the linear SCC order cannot be permuted.
No classification of all strong minimum cases from Busch is needed here.
The exceptional order-five multiplicity is not asserted by this formula.

## 5. Consequences and scope boundaries

* **Inherited, sharpened:** H=h>1 is attainable iff some core on at most
  K(h) vertices has H=h. This remains an effective finite algorithm,
  not a proof every allowed odd h is accepted. h=1 uses a singleton or
  the empty convention.
* **Inherited:** 7 and 21 are absent at every order. Their proof needs
  the strong floors and small strong spectra; it does not follow from
  the budget K alone. See THM-1370 and the separate Session133 atlas.
* **Fresh refinement:** the exact core-order floor and its inverse
  three-candidate budget eliminate integer slack in the older sharp
  continuous logarithmic inequality. Equality models give explicit
  witnesses and the extremizer-count sequence.
* **Open:** whether every positive odd value other than 7 and 21 occurs.
  None of these bounds or recurrences supplies the unbounded existence
  construction needed for that claim.

## 6. Threshold recognition for a given tournament

**Corollary C (algorithmic; no novelty claim).** Given an n-vertex tournament
T by its adjacency matrix and an integer B>=1, one can decide whether
H(T)<=B and, in the positive case, return H(T) exactly. Put
`alpha=3 log_5(2)` (approximately 1.29203). The cost is

```
O(n^2) graph operations
    + O(B^alpha * log^3(B+2)) bit operations.                 (C)
```

Here a graph operation accesses an adjacency entry or a vertex index.
With conventional bit-cost accounting for O(log(n+2))-bit indices, the
first term may conservatively be written `O(n^2 log(n+2))` bit operations.
This is polynomial in the numeric threshold B, not in its binary length.

**Algorithm and proof.**

1. Find the SCCs and their sizes by a standard graph traversal in O(n^2)
   graph operations. Ignore singleton SCCs, whose H factors are one. Let
   their remaining sizes be s_1,...,s_r, with sum c(T).
2. Compute K(B) using (A) and integer arithmetic. If c(T)>K(B), reject:
   Theorem A proves H(T)>B. Otherwise the total number of remaining
   vertices is O(log(B+2)); in particular r=O(log(B+2)).
3. For each remaining SCC, reject if
   `5^(s_i-1)>B^3`. Busch gives `H(SCC)^3>=5^(s_i-1)`, so this test is
   sound. Powers can be evaluated with early saturation at B^3+1.
   Every surviving component has
   `s_i<=L(B):=1+floor(3 log_5 B)`.
4. For each surviving SCC use the standard subset DP, storing
   `D[S,v]=min(B+1, number of directed paths with vertex set S ending v)`.
   The singleton entries are one. For |S|>=2, sum `D[S-{v},u]` over
   `u in S-{v}` with `u -> v`, saturating every addition at B+1. Sum the
   full-set endpoints, again saturated, to obtain the component count.
   Multiply component counts with saturation at B+1. Return that product
   when at most B; otherwise report H(T)>B. If no nonsingleton SCC
   remains, return one.

Saturation is exact for this purpose: for nonnegative integers,
`min(B+1,x+y)` and `min(B+1,xy)` are unchanged when each argument is
first capped at B+1. Induction over the DP states therefore gives the
true capped component counts. SCC multiplicativity gives the true capped
whole-tournament count, and an answer at most B is consequently exact.

For the time bound, component i uses O(s_i^2*2^s_i) additions on
O(log(B+2))-bit integers. After the two early tests,

```
2^s_i <= 2*B^alpha,
sum_i s_i = c(T) <= K(B) = O(log(B+2)),
sum_i s_i^2 <= (sum_i s_i)^2 = O(log^2(B+2)).
```

Summing the DP costs proves (C). Computing K(B), capped strong-bound
powers, and the final capped product takes at most an additional
polynomial in log(B+2), absorbed by (C). The cases B=1,2 can simply test
whether the core is empty. The DP working memory, apart from the input,
can be bounded by `O(B^alpha log^2(B+2))` bits by processing components
one at a time.

**Scope firewall.** This is **recognition/counting for a supplied T**.
It does not search for a tournament with prescribed H=h. The latter
existence problem still has only the inherited finite-search bound here,
and unrestricted positive-odd spectrum completeness remains open. The
n^2 input work cannot be removed merely because the cyclic core is small:
arbitrarily many singleton SCCs can pad a fixed-H tournament.

## 7. Verification and connection contract

Run from this directory:

```
python3 check_core.py --out normal.json
python3 -O check_core.py --out optimized.json
cmp normal.json optimized.json
```

The checker compares g with a separate component-partition DP using f
through core order 400, compares restricted/unrestricted resource budgets
and the inverse floor for 1<=M<=100000, enumerates every labelled
tournament through order six, checks every induced subset's equality
criterion, and compares Hamiltonian DP with literal permutations through
order five. It checks 3/4-block word extremizer counts through order 25, explicit extremizing
budgets, and transitive-module
hostiles. No complete order-seven census or large-target spectrum
completeness is asserted. Large witness H is certified by SCC factors
where direct DP is explicitly marked false.

Connection: tournament → ordered nonsingleton SCC core → 3/4-block
resource replacement. First arrow preserves H exactly; second preserves
only inequalities (size nondecrease, H nonincrease). The lost exact H is
the crucial sidecar. The exact optimization is valid for a budget, not
for the fibre H=M. The cheapest hostile is M=21 with optimizer H=15.
