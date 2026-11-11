# A lawful domination-to-tournament adapter for chain graphs

> **Continuation:** the general inequality and exact equality criterion left open below are now proved in [BIPARTITE_DEFECT.md](BIPARTITE_DEFECT.md). The earlier exploration is retained as a discovery record.

**Status: PROVED adapter, formula and digit operations. PROVED obstruction with a tiny finite base. The broader bipartite inequality at the end is only FINITE-EXACT evidence and OPEN.**

## 1. Inheritance and the precise target

THM-5040 (`sharp-even-graph-chan-pak`), THM-5055 (`even-graph-chan-pak-theorem-dominating-set-counts-of-even-graphs-are-all-odd-numbers-except-3-7-13-with-sharp-orders`) and THM-5056 (`odd-graph-chan-pak-dominating-set-counts-of-odd-graphs-are-all-odd-numbers-except-1-5-37-on-ceil-log2-vertices`) solve the stated graph-domination spectra. Their parity notion is von Broemssen's automorphism edge-sign, not vertex-degree parity. None supplies a d-to-H bijection. THM-4996 already warns that evaluating independence polynomials of even graphs does not repair this.

The inherited one-apex tournament formulas of THM-4748 and the OCF, THM-002, are the positive mechanism. The new object here is an explicit graph subclass on which domination count really is the same statistic. The conflict graph is not the input graph: this distinction avoids the forbidden universal I(G,2) claim (K3 would give7).

Concept board: nested bipartite neighbourhoods; an apex over a transitive order; odd cycles through one marked vertex; domination boundary conditions; lawful two-register digits; the first graph outside the adapter.

## 2. Actual graph and tournament constructors

Take a word w=w_1...w_n over {A,B}.

- G_w is the undirected bipartite graph on positions1,...,n, with edge ij exactly when i<j, w_i=A and w_j=B.
- T_w has the transitive tournament 1->2->...->n and one new vertex z. Set z->i for letter A and i->z for letter B.

Every bipartite chain graph (equivalently, a bipartite graph with nested neighbourhoods on one shore) has such a word representation. Arrange the A vertices by decreasing neighbourhood inclusion and interleave each B at its neighbourhood threshold. Leading B letters and trailing A letters represent isolated vertices. Conversely every word clearly gives nested A-neighbourhoods.

**Theorem.** For every binary word,

    d(G_w) = H(T_w)
           = 1 + sum_{i<j, w_i=A,w_j=B} c(j-i),             (1)

where c(1)=2 and c(t)=2^(t-1) for t>=2. Thus every chain graph has a literal count-preserving tournament realizer using one additional vertex.

### Tournament side: exact OCF cycle description

Every directed cycle of T_w contains z, because deletion of z is transitive. A directed cycle is z followed by an increasing subsequence whose first vertex is A and last is B. The cycle is odd exactly when that subsequence has even cardinality. All such cycles intersect at z, so at most one belongs to any odd-cycle packing.

Fix endpoints i<j. If j=i+1, there is exactly one such cycle. Otherwise there are 2^(j-i-2) even-cardinality choices from the j-i-1 interior positions. OCF multiplies the cycle count by2 and adds the empty packing, giving (1).

This makes Omega(T_w) a complete graph on the allowed odd cycles, but only particular clique sizes arise. It does not assert that arbitrary complete graphs are realizable conflict graphs.

### Domination side: a boundary recurrence proves the same formula

Leading B and trailing A vertices are isolated in G_w, so deleting them leaves d unchanged. They occur in no term of (1). Reduce to a core word beginning A and ending B. Empty and one-letter words have d=1, and the two-letter core AB has d=3.

For a core of length n>=3, let w' delete its first and last letters, and put

    P=sum_{j<n,w_j=B} 2^(j-2),
    Q=sum_{i>1,w_i=A} 2^(n-i-1).

Classify dominating sets by membership of the first A and last B.

- Both selected: 2^(n-2) sets, because these two vertices dominate every vertex.
- First selected, last omitted: 1+P sets. With no selected B, all A vertices must be selected. Otherwise let j be the last selected B; every A after j is forced, and the earlier choices give 2^(j-2).
- Last selected, first omitted: symmetrically 1+Q sets.
- Both omitted: d(G_w') minus the dominating sets of G_w' contained wholly in one shore. There is exactly one A-only such set iff w_2=A, and one B-only such set iff w_(n-1)=B. Since w' is nonempty, these sets are distinct.

Therefore

    d(G_w)-d(G_w')
      =2^(n-2)+P+Q+[w_2=B]+[w_(n-1)=A].                   (2)

The difference of the right side of (1) between w and w' is exactly (2): count endpoint pairs involving the first A or last B, subtract their double-counted common pair, and retain the exceptional adjacent-pair coefficient c(1)=2. Induction proves (1).

## 3. Concrete lawful subclasses and composition

For w=A^a B^b, G_w=K_(a,b) and T_w=C3(TT_a,TT_b,1), giving

    d(K_(a,b))=H(T_w)=(2^a-1)(2^b-1)+2.

Thus the arithmetic lane's two-run family is an actual domination adapter, not just an equal-looking formula. For w=A^(a-1) B A B^(b-1), the graph is K_(a,b) with one edge removed, and the tournament is its boundary crossing flip from CROSSING_FLIP.md. The equality of their decrements is explained by the same adapter.

Disjoint unions of chain graphs lift to ordinal sums of their tournaments:

    d(G_1 disjoint-union G_2)=d(G_1)d(G_2),
    H(T_1 -> T_2)=H(T_1)H(T_2).

This is a lawful multiplicative closure with explicit graphs. It is not an addition operation or a closure under the even-graph hub moves.

## 4. A genuine two-register digit lift

Let x be the integer represented by the word with A=1 and B=0, preserving leading zeroes only as harmless isolates. Start at (h,x)=(1,0). The actual word append operations give

    append A: (h,x) -> (h, 2x+1),
    append B: (h,x) -> (h+x+(x mod2), 2x).                 (3)

To derive this, let E and O count respectively even and odd nonempty subsequences beginning with A. Appending a letter sends

    E' = E+O,
    O' = E+O+[new letter=A].

Inductively E=floor(x/2), O=ceil(x/2). Appending A creates no new subsequence ending B; appending B creates O even subsequences beginning A and ending B. OCF therefore increases H by2O=x+(x mod2). Equation (1) transfers these very same operations to graph domination.

The sidecar x is essential: h alone does not determine the next increment. Even h=1 has states x=0 (empty word), x=1 (A), x=3 (AA), whose next B creates h=1,3,5. Appending B is a hub over the A shore, not a universal graph cone and not the PH move of THM-5055.

## 5. Why the solved even-graph theorem does not transfer wholesale

The adapter itself exposes a sharp family boundary. For a core of length n>=3, its first A/last B pair alone contributes2^(n-2), so

    H(T_w)=d(G_w)>=1+2^(n-2).                              (4)

Consequently a chain graph with d=29 would have a core of length at most6. Enumerating the 31 core words of lengths2,...,6 directly by (1) gives exactly

    {3,5,9,11,15,17,19,23,25,27,31,33,37,39,45,47,51}.

Thus29 is absent from chain graphs at every order. Since29 is prime, it is absent from disjoint unions of chain graphs too. The same argument excludes13 (already using length at most5). These are **PROVED using the explicit bounded finite base**, not an extrapolation from a length cutoff.

But29 has the even graph realizer with edges

    01,02,03,34,45.

This six-vertex tree has29 dominating sets. Its only nonidentity automorphism exchanges leaves1 and2, which reverses no edge orientations under this labelling, so it is even. Therefore the even-graph spectrum already leaves the adapter at29, even though29 is a perfectly valid tournament count by other constructions.

A structural operator obstruction is visible earlier: PH(P3, pendant,1), one of THM-5055's lawful operations, contains a triangle (the new hub and the old pendant's neighbour and opposite leaf). It is not bipartite, hence not a chain graph. It may have a different chain graph with the same d, but that would require a new constructor or a target-dependent rebuild. The hub proof does not supply one.

The strongest positive survivor is consequently a full explicit adapter and its multiplicative closure, together with the two-register grammar (3). It does not lift all graph digit moves, nor assert Chan–Pak completeness.

## 6. A further positive signal, explicitly OPEN

For an arbitrary ordered bipartite graph G with shores A,B, form a tournament on A,B,z as follows: each shore is transitive; orient a cross pair A_i->B_j iff it is an edge of G, and reverse it otherwise; put z->A and B->z. The chain-graph adapter is the case where the tournament after removing z is transitive.

A small scout checked every graph with shore sizes(2,2),(2,3),(3,3), respectively16,64,512 graphs. In all cases

    H(T_G)>=d(G),

and equality occurred in6,10,20 cases. These numbers equal the respective counts of compatible total-order interleavings. This is **FINITE-EXACT evidence only**: neither the general inequality nor its equality characterization is proved here. The promising question is whether a domination-set injection produces a nonnegative defect indexed by cycles of T_G-z. Such an identity would retain a real graph interface while recording the extra cycles OCF must add. It would still need a defect-elimination constructor to transfer the solved graph spectrum.

## 7. Reproduction and connection ledger

`domination_check.py` checks the adapter by direct dominating-subset enumeration and independent tournament subset DP for all511 words of length at most8, as well as the register formula, the finite core values, the tree29 witness, and the stated exploratory bipartite universe. Results: `domination_receipt.json`.

Source: chain-graph dominating sets. Target: apex-over-order tournament Hamilton counts. Map: an interleaving of the two nested shores, then one apex. Preserved predicate: the exact count. Lost information: the domination sets themselves are not identified bijectively with Hamilton paths by the proof; graph automorphism sign need not be retained. Required sidecar: the word/interleaving, or x for (3). Hostile: the even29 tree and the nonbipartite PH move. Next test: prove or refute the proposed nonnegative defect formula for arbitrary bipartite graphs.
